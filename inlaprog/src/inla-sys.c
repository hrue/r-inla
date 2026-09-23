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
static int p_cores_verbose = 0;

#if defined(__linux__)
#       include <stdio.h>
#       include <stdlib.h>
#       include <unistd.h>
#       include <sys/stat.h>
#       include <sched.h>
#       include <hwloc.h>
#       include <hwloc/cpukinds.h>			       // Handles modern Intel & AMD hybrid architectures

#       define MAX_CORES (2*1024)
static int *g_p_core_pu_ids = NULL;
static int g_p_core_count = 0;

int inla_num_p_cores(void)
{
        g_p_core_count = 0;
        Free(g_p_core_pu_ids);
        g_p_core_pu_ids = Calloc(MAX_CORES, int);

        hwloc_topology_t topology;

        if (hwloc_topology_init(&topology) < 0)
                return NUM_P_CORES_DEFAULT();
        hwloc_topology_load(topology);

        // Query different core types. hwloc sorts kinds automatically by efficiency.
        // The highest-performance group (P-Cores) is always the LAST index.
        int num_kinds = hwloc_cpukinds_get_nr(topology, 0);
        int total_cores = hwloc_get_nbobjs_by_type(topology, HWLOC_OBJ_CORE);

        if (num_kinds > 1) {
                // --- HYBRID ARCHITECTURE FOUND VIA OS (Intel P/E or AMD Zen/Zen-c) ---
                hwloc_bitmap_t p_core_cpuset = hwloc_bitmap_alloc();
                int best_kind_index = num_kinds - 1;

                hwloc_cpukinds_get_info(topology, best_kind_index, p_core_cpuset, NULL, NULL, NULL, 0);

                for (int i = 0; i < total_cores && g_p_core_count < MAX_CORES; i++) {
                        hwloc_obj_t core_obj = hwloc_get_obj_by_type(topology, HWLOC_OBJ_CORE, i);

                        if (core_obj && hwloc_bitmap_intersects(core_obj->cpuset, p_core_cpuset)) {
                                // Harvest the lowest logical PU ID belonging to this physical core (skips hyperthreads)
                                g_p_core_pu_ids[g_p_core_count] = hwloc_bitmap_first(core_obj->cpuset);
                                g_p_core_count++;
                        }
                }
                hwloc_bitmap_free(p_core_cpuset);
                if (p_cores_verbose) {
                        printf("[Linux Detected] Hybrid layout found via OS kinds. Using %d performance cores.\n", g_p_core_count);
                }
        } else {
                // --- FALLBACK INTERMEDIATE STEP: DYNAMIC MAX L2 CACHE SCAN ---
                size_t max_l2_size = 0;

                // Step A: Find the absolute largest L2 cache size active on this chip
                for (int i = 0; i < total_cores; i++) {
                        hwloc_obj_t core_obj = hwloc_get_obj_by_type(topology, HWLOC_OBJ_CORE, i);
                        if (core_obj) {
                                hwloc_obj_t parent = core_obj->parent;
                                while (parent && parent->type != HWLOC_OBJ_L2CACHE) {
                                        parent = parent->parent;
                                }
                                if (parent && parent->attr && parent->attr->cache.size > max_l2_size) {
                                        max_l2_size = parent->attr->cache.size;
                                }
                        }
                }

                // Step B: Determine if there is an unequal/heterogeneous cache layout
                int is_heterogeneous_cache = 0;
                if (max_l2_size > 0) {
                        for (int i = 0; i < total_cores; i++) {
                                hwloc_obj_t core_obj = hwloc_get_obj_by_type(topology, HWLOC_OBJ_CORE, i);
                                if (core_obj) {
                                        hwloc_obj_t parent = core_obj->parent;
                                        while (parent && parent->type != HWLOC_OBJ_L2CACHE) {
                                                parent = parent->parent;
                                        }
                                        if (parent && parent->attr && parent->attr->cache.size < max_l2_size) {
                                                is_heterogeneous_cache = 1; // Found an efficiency core with a smaller cache!
                                                break;
                                        }
                                }
                        }
                }

                if (is_heterogeneous_cache) {
                        // Gather ONLY the cores attached to the maximum discovered L2 Cache tier (P-Cores)
                        for (int i = 0; i < total_cores && g_p_core_count < MAX_CORES; i++) {
                                hwloc_obj_t core_obj = hwloc_get_obj_by_type(topology, HWLOC_OBJ_CORE, i);
                                if (core_obj) {
                                        hwloc_obj_t parent = core_obj->parent;
                                        while (parent && parent->type != HWLOC_OBJ_L2CACHE) {
                                                parent = parent->parent;
                                        }
                                        if (parent && parent->attr && parent->attr->cache.size == max_l2_size) {
                                                g_p_core_pu_ids[g_p_core_count] = hwloc_bitmap_first(core_obj->cpuset);
                                                g_p_core_count++;
                                        }
                                }
                        }
                        if (p_cores_verbose) {
                                printf("[Linux Detected] Hybrid layout found via Dynamic Cache inspection (Max L2: %zu KB). Using %d performance cores.\n", 
                                       max_l2_size / 1024, g_p_core_count);
                        }
                } else {
                        // --- STANDARD SYMMETRIC ARCHITECTURE FOUND (All caches are identical) ---
                        for (int i = 0; i < total_cores && g_p_core_count < MAX_CORES; i++) {
                                hwloc_obj_t core_obj = hwloc_get_obj_by_type(topology, HWLOC_OBJ_CORE, i);

                                if (core_obj) {
                                        g_p_core_pu_ids[g_p_core_count] = hwloc_bitmap_first(core_obj->cpuset);
                                        g_p_core_count++;
                                }
                        }
                        if (p_cores_verbose) {
                                printf("[Linux Detected] Symmetric architecture found. Using all %d physical cores.\n", g_p_core_count);
                        }
                }
        }
        hwloc_topology_destroy(topology);

        if (p_cores_verbose) {
                printf("Target Logical OS PU IDs for pinning: ");
                for (int i = 0; i < g_p_core_count; i++)
                        printf("%d ", g_p_core_pu_ids[i]);
                printf("\n\n");
        }

        return g_p_core_count;
}

int inla_lock_to_p_cores(void)
{
	return 0;
}

int inla_lock_thread_to_p_core(int tid)
{
	if (g_p_core_count == 0)
		return 0;
	if (LEGAL(tid, g_p_core_count)) {
		int target_pu = g_p_core_pu_ids[tid];
		cpu_set_t mask;

		CPU_ZERO(&mask);
		CPU_SET(target_pu, &mask);
		if (pthread_setaffinity_np(pthread_self(), sizeof(cpu_set_t), &mask) == 0) {
			if (p_cores_verbose) {
				int actual_pu = sched_getcpu();

				printf("[Thread %d] Fixed to P-Core (PU %d). OS Verification: Active on PU %d\n", tid, target_pu, actual_pu);
			}
			return 0;
		} else {
			perror("Linux pthread_setaffinity_np failed");
			return 1;
		}
	}
	return 0;
}
#endif

#if defined(__APPLE__)
#       include <stdio.h>
#       include <sys/types.h>
#       include <sys/sysctl.h>
#       include <pthread.h>

int inla_num_p_cores(void)
{
	// query macOS kernel attributes directly via sysctl
	int64_t p_cores = 0;
	size_t size = sizeof(p_cores);

	if (sysctlbyname("hw.perflevel0.physicalcpu", &p_cores, &size, NULL, 0) != 0) {
		// Fallback to total physical cores if the metric is missing
		int64_t total_cores = 0;
		size_t total_size = sizeof(total_cores);

		sysctlbyname("hw.physicalcpu", &total_cores, &total_size, NULL, 0);
		p_cores = total_cores;
	}
	return (p_cores > 0 ? p_cores : NUM_P_CORES_DEFAULT());
}

int inla_lock_to_p_cores(void)
{
	// lock == 'bind' here
	// Elevate the current thread to the absolute highest performance tier.
	// This strictly forces macOS to run your code on the P-Cores.
	int result = pthread_set_qos_class_self_np(QOS_CLASS_USER_INITIATED, 0);

	if (result != 0) {
		perror("Failed to set Mac QoS class");
		return 1;
	}
	return 0;
}

int inla_lock_thread_to_p_core(int tid)
{
	return 0;
}
#endif

#if defined(_WIN32)
#       include <windows.h>
int inla_lock_to_p_cores(void)
{
	// not yet implemented. very different on Windows, not sure its worth while
	return 0;
}

int inla_lock_thread_to_p_core(int tid)
{
	return 0;
}

int inla_num_p_cores(void)
{
	int p_cores_default = NUM_P_CORES_DEFAULT();
	DWORD bufferSize = 0;

	// First call to determine the required buffer size
	if (!GetLogicalProcessorInformationEx(RelationProcessorCore, NULL, &bufferSize)) {
		if (GetLastError() != ERROR_INSUFFICIENT_BUFFER) {
			if (p_cores_verbose) {
				fprintf(stderr, "Error determining buffer size. Code: %lu\n", GetLastError());
			}
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

	if (p_cores_verbose) {
		if (maxEfficiency > 0) {
			printf("Performance Cores (P-Cores): %d\n", pCoreCount);
			printf("Efficiency Cores (E-Cores):  %d\n", eCoreCount);
		} else {
			printf("Performance Cores (P-Cores): %d (Symmetric CPU Architecture)\n", pCoreCount);
			printf("Efficiency Cores (E-Cores):  0\n");
		}
	}
	return (pCoreCount > 0 ? pCoreCount : NUM_P_CORES_DEFAULT());
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
