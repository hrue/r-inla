#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#if defined(__linux__)
#       if !defined(_GNU_SOURCE)
#              define _GNU_SOURCE			       // Required for CPU affinity macros
#       endif
#       include <sched.h>
#       include <ftw.h>
#       include <unistd.h>
#endif

int inla_ncpu(void)
{
#if defined(_SC_NPROCESSORS_ONLN)			       /* Linux, Solaris, AIX */
	return (int) sysconf(_SC_NPROCESSORS_ONLN);
#elif defined(__APPLE__)				       /* MacOSX */
	int count = -1;
	size_t size = sizeof(count);
	sysctlbyname("hw.ncpu", &count, &size, NULL, 0);
	return count;
#elif defined(_WIN32)
	SYSTEM_INFO SystemInfo;
	GetSystemInfo(&SystemInfo);
	return SystemInfo.dwNumberOfProcessors;
#else
	return -1;
#endif
}


#if defined(__linux__)
int inla_remove_dir_callback(const char *dirname, const struct stat *UNUSED(sb), int typeflag, struct FTW *UNUSED(ftwbuf))
{
	if (typeflag == FTW_F || typeflag == FTW_SL) {
		unlink(dirname);
	} else if (typeflag == FTW_DP) {
		rmdir(dirname);
	}
	return 0;
}
void inla_remove_dir(char *dirname)
{
	nftw(dirname, inla_remove_dir_callback, 10, FTW_DEPTH | FTW_PHYS);
}
#else
void inla_remove_dir(char *UNUSED(dirname))
{
}
#endif

int inla_mkdir(const char *dirname)
{
#if defined(_WIN32)
	return mkdir(dirname);
#else
	return mkdir(dirname, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
#endif
}

// this is from https://stackoverflow.com/questions/2513505/how-to-get-available-memory-c-g
// return RAM in Mb
#if defined(_WIN32)
unsigned long long getTotalSystemMemory()
{
	MEMORYSTATUSEX status;
	status.dwLength = sizeof(status);
	GlobalMemoryStatusEx(&status);
	return ((status.ullTotalPhys / 1024L / 1024L));
}
#else
unsigned long long getTotalSystemMemory()
{
	long pages = sysconf(_SC_PHYS_PAGES);
	long page_size = sysconf(_SC_PAGE_SIZE);
	return ((pages * page_size) / 1024L / 1024L);
}
#endif

#if defined(_WIN32)
void inla_signal(int UNUSED(sig))
{
	/*
	 * ... 
	 */
}
#else
void inla_signal(int sig)
{
	switch (sig) {
	case SIGUSR1:
		GMRFLib_write_state = 1;
		fprintf(stdout, "\n\n*** set GMRFLib_write_state = 1\n\n");
		break;
	case SIGUSR2:
		break;
	default:
		_exit(sig);
		break;
	}
	return;
}
#endif

int inla_endian(void)
{
	int x = 1;
	return ((*(char *) &x) ? INLA_LITTLE_ENDIAN : INLA_BIG_ENDIAN);
}

int inla_parse_libR(inla_tp *mb, dictionary *ini, int sec)
{
	/*
	 * parse section = libR
	 */
	char *secname = NULL, *env = NULL;

	if (mb->verbose) {
		printf("\tinla_parse_libR...\n");
	}
	secname = Strdup(iniparser_getsecname(ini, sec));
	if (mb->verbose) {
		printf("\t\tsection[%s]\n", secname);
	}

	mb->libR_R_HOME = Strdup(iniparser_getstring(ini, inla_string_join(secname, "R_HOME"), NULL));
	inla_set_R_home(mb->libR_R_HOME);
	if (mb->verbose) {
		printf("\t\t\tR_HOME=[%s]\n", mb->libR_R_HOME);
	}

	if (mb->libR_R_HOME) {
		// set the R_HOME variable
		GMRFLib_sprintf(&env, "R_HOME=%s", mb->libR_R_HOME);
		my_setenv(env, 0);
		Free(env);
	}

	return INLA_OK;
}

int inla_tolower(char *string)
{
	if (string) {
		int i;
		for (i = 0; i < (int) strlen(string); i++) {
			string[i] = (char) tolower((int) string[i]);
		}
	}
	return GMRFLib_SUCCESS;
}

#define NUM_P_CORES_DEFAULT() IMAX(omp_get_max_threads(), omp_get_num_procs())

#if defined(__linux__)
// Automatically parses standard Linux core range strings (e.g., "0-11" or "0-7,16-23")
static void parse_and_add_cpus(const char *str, cpu_set_t *cpuset)
{
	char *dup = strdup(str);
	char *token = strtok(dup, ",\n");
	while (token != NULL) {
		int start, end;
		// Check if token is a range (e.g., "0-7") or a single core (e.g., "0")
		if (sscanf(token, "%d-%d", &start, &end) == 2) {
			for (int i = start; i <= end; i++) {
				CPU_SET(i, cpuset);
			}
		} else if (sscanf(token, "%d", &start) == 1) {
			CPU_SET(start, cpuset);
		}
		token = strtok(NULL, ",\n");
	}
	free(dup);
}

int inla_lock_to_p_cores(void)
{
	cpu_set_t cpuset;
	CPU_ZERO(&cpuset);

	// 1. Open the Linux kernel file that stores the P-core mappings
	FILE *f = fopen("/sys/devices/cpu_core/cpus", "r");
	if (!f) {
		// Fallback: If the file isn't there, it might be an older CPU with no E-cores
		perror("Could not read P-core layout (non-hybrid CPU?)");
		return 1;
	}

	char buf[256];
	if (fgets(buf, sizeof(buf), f) != NULL) {
		printf("Detected P-core range string: %s", buf);
		parse_and_add_cpus(buf, &cpuset);
	}
	fclose(f);

	// 2. Apply the parsed P-core mask to our running thread
	if (sched_setaffinity(0, sizeof(cpu_set_t), &cpuset) != 0) {
		perror("sched_setaffinity failed");
		return 1;
	}
	return 0;
}

static int parse_max(const char *str)
{
	char *dup = strdup(str);
	char *token = strtok(dup, ",\n");
	while (token != NULL) {
		int start, end;
		// Check if token is a range (e.g., "0-7") or a single core (e.g., "0")
		if (sscanf(token, "%d-%d", &start, &end) == 2) {
			free(dup);
			return end;
		} else if (sscanf(token, "%d", &start) == 1) {
			free(dup);
			return start;
		}
		token = strtok(NULL, ",\n");
	}
	free(dup);
	return 0;
}
int inla_num_p_cores(void)
{
	cpu_set_t cpuset;
	CPU_ZERO(&cpuset);

	// 1. Open the Linux kernel file that stores the P-core mappings
	FILE *f = fopen("/sys/devices/cpu_core/cpus", "r");
	if (!f) {
		// Fallback: If the file isn't there, it might be an older CPU with no E-cores
		perror("Could not read P-core layout (non-hybrid CPU?)");
		return NUM_P_CORES_DEFAULT();
	}

	int num_p = 0;
	char buf[256];
	if (fgets(buf, sizeof(buf), f) != NULL) {
		// printf("Detected P-core range string: %s", buf);
		num_p = parse_max(buf);
	}
	fclose(f);
	return (num_p > 0 ? num_p : NUM_P_CORES_DEFAULT());
}

#       if 0
int main(void)
{
	printf("Starting program...\n");

	if (lock_to_p_cores() == 0) {
		printf("Success! Program successfully locked to P-cores.\n");
	} else {
		printf("Running on default OS cores due to fallback.\n");
	}
	return 0;
}
#       endif
#endif

#if defined(__APPLE__)
#       include <stdio.h>
#       include <sys/types.h>
#       include <sys/sysctl.h>
#       include <pthread.h>

int inla_num_p_cores(void)
{
	int p_cores = 0;
	size_t size = sizeof(p_cores);
	sysctlbyname("hw.perflevel0.physicalcpu", &p_cores, &size, NULL, 0);
	if (p_cores > 0) {
		return p_cores;
	} else {
		return NUM_P_CORES_DEFAULT();
	}
}

int inla_lock_to_p_cores(void)
{
	// lock == 'bind' here
	// Elevate the current thread to the absolute highest performance tier.
	// This strictly forces macOS to run your code on the P-Cores.
	int result = pthread_set_qos_class_self_np(QOS_CLASS_USER_INTERACTIVE, 0);
	if (result != 0) {
		perror("Failed to set Mac QoS class");
		return 1;
	}
	return 0;
}

#       if 0
int main(void)
{
	// 1. Get the P-Core count using sysctl
	int p_cores = 0;
	size_t size = sizeof(p_cores);
	sysctlbyname("hw.perflevel0.physicalcpu", &p_cores, &size, NULL, 0);
	printf("Detected P-Cores available: %d\n", p_cores);

	// 2. Force this thread onto the P-cores
	inla_lock_to_p_cores();
	// --- Run your max double vector benchmarks here ---
	// macOS will execute this loop on the ultra-fast P-cores natively.
	return 0;
}
#       endif
#endif

#if defined(_WIN32)
#       include <windows.h>
int inla_lock_to_p_cores(void)
{
	// not yet implemented. very different on Windows, not sure its worth while
	return 0;
}
int inla_num_p_cores(void)
{
	int p_cores_default = NUM_P_CORES_DEFAULT();
	DWORD bufferSize = 0;

	// First call to determine the required buffer size
	if (!GetLogicalProcessorInformationEx(RelationProcessorCore, NULL, &bufferSize)) {
		if (GetLastError() != ERROR_INSUFFICIENT_BUFFER) {
			fprintf(stderr, "Error determining buffer size. Code: %lu\n", GetLastError());
			return p_cores_default;
		}
	}

	PSYSTEM_LOGICAL_PROCESSOR_INFORMATION_EX buffer = (PSYSTEM_LOGICAL_PROCESSOR_INFORMATION_EX) malloc(bufferSize);
	if (!buffer) {
		fprintf(stderr, "Memory allocation failed.\n");
		return p_cores_default;
	}

	// Second call to actually populate the buffer
	if (!GetLogicalProcessorInformationEx(RelationProcessorCore, buffer, &bufferSize)) {
		fprintf(stderr, "Error retrieving processor information. Code: %lu\n", GetLastError());
		free(buffer);
		return p_cores_default;
	}

	int totalPhysicalCores = 0;
	BYTE maxEfficiency = 0;

	// Step 1: Find the maximum EfficiencyClass value across all cores
	unsigned char *ptr = (unsigned char *) buffer;
	unsigned char *end = ptr + bufferSize;

	while (ptr < end) {
		PSYSTEM_LOGICAL_PROCESSOR_INFORMATION_EX info = (PSYSTEM_LOGICAL_PROCESSOR_INFORMATION_EX) ptr;
		if (info->Relationship == RelationProcessorCore) {
			totalPhysicalCores++;
			if (info->Processor.EfficiencyClass > maxEfficiency) {
				maxEfficiency = info->Processor.EfficiencyClass;
			}
		}
		ptr += info->Size;
	}

	// Step 2: Count how many cores belong to that maximum efficiency class
	int pCoreCount = 0;
	int POSSIBLY_UNUSED(eCoreCount) = 0;
	ptr = (unsigned char *) buffer;
	while (ptr < end) {
		PSYSTEM_LOGICAL_PROCESSOR_INFORMATION_EX info = (PSYSTEM_LOGICAL_PROCESSOR_INFORMATION_EX) ptr;
		if (info->Relationship == RelationProcessorCore) {
			// If the max efficiency is 0, the system is symmetric (all cores are the same)
			if (maxEfficiency == 0) {
				pCoreCount = totalPhysicalCores;
				break;
			} else {
				if (info->Processor.EfficiencyClass == maxEfficiency) {
					pCoreCount++;
				} else {
					eCoreCount++;
				}
			}
		}
		ptr += info->Size;
	}
	free(buffer);

	return pCoreCount;

#       if 0
	if (maxEfficiency > 0) {
		printf("Performance Cores (P-Cores): %d\n", pCoreCount);
		printf("Efficiency Cores (E-Cores):  %d\n", eCoreCount);
	} else {
		printf("Performance Cores (P-Cores): %d (Symmetric CPU Architecture)\n", pCoreCount);
		printf("Efficiency Cores (E-Cores):  0\n");
	}
	return 0;
#       endif
}
#endif

#if !defined(__linux__) && !defined(__APPLE__) && !defined(_WIN32)
int inla_lock_to_p_cores(void)
{
	return 0;
}
int inla_num_p_cores(void)
{
	return 0;
}
#endif

