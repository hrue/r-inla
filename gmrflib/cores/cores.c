
/*
 * Copyright (c) Przemyslaw Skibinski, Yann Collet, Facebook, Inc.
 * All rights reserved.
 *
 * This source code is licensed under both the BSD-style license (found in the
 * LICENSE file in the root directory of this source tree) and the GPLv2 (found
 * in the COPYING file in the root directory of this source tree).
 * You may select, at your option, one of the above-listed licenses.
 */

#include "GMRFLib/cores.h"
#include <stdlib.h>					       /* malloc, realloc, free */
#include <stdio.h>					       /* fprintf */
#include <time.h>					       /* clock_t, clock, CLOCKS_PER_SEC, nanosleep */
#include <errno.h>
#include <assert.h>

#if defined(_WIN32)
#include <sys/utime.h>					       /* utime */
#include <io.h>						       /* _chmod */
#else
#include <unistd.h>					       /* chown, stat */
#define PLATFORM_POSIX_VERSION _POSIX_VERSION
#if PLATFORM_POSIX_VERSION < 200809L || !defined(st_mtime)
#include <utime.h>					       /* utime */
#else
#include <fcntl.h>					       /* AT_FDCWD */
#include <sys/stat.h>					       /* utimensat */
#endif
#endif

#if defined(_MSC_VER) || defined(__MINGW32__) || defined (__MSVCRT__)
#include <direct.h>					       /* needed for _mkdir in windows */
#endif

#if defined(__linux__) || (PLATFORM_POSIX_VERSION >= 200112L)  /* opendir, readdir require POSIX.1-2001 */
#include <dirent.h>					       /* opendir, readdir */
#include <string.h>					       /* strerror, memcpy */
#endif							       /* #ifdef _WIN32 */

#if defined(_MSC_VER)
#define chmod _chmod
#endif




/*-****************************************
*  count the number of cores
******************************************/

#if defined(_WIN32) || defined(WIN32)

#include <windows.h>

typedef BOOL(WINAPI * LPFN_GLPI) (PSYSTEM_LOGICAL_PROCESSOR_INFORMATION, PDWORD);

DWORD CountSetBits(ULONG_PTR bitMask)
{
	DWORD LSHIFT = sizeof(ULONG_PTR) * 8 - 1;
	DWORD bitSetCount = 0;
	ULONG_PTR bitTest = (ULONG_PTR) 1 << LSHIFT;
	DWORD i;

	for (i = 0; i <= LSHIFT; ++i) {
		bitSetCount += ((bitMask & bitTest) ? 1 : 0);
		bitTest /= 2;
	}
	return bitSetCount;
}

int UTIL_countCores(int logical)
{
	static int numCores = 0;
	if (numCores != 0)
		return numCores;

	{
		LPFN_GLPI glpi;
		BOOL done = FALSE;
		PSYSTEM_LOGICAL_PROCESSOR_INFORMATION buffer = NULL;
		PSYSTEM_LOGICAL_PROCESSOR_INFORMATION ptr = NULL;
		DWORD returnLength = 0;
		size_t byteOffset = 0;

#if defined(_MSC_VER)

/* Visual Studio does not like the following cast */
#pragma warning( disable : 4054 )			       /* conversion from function ptr to data ptr */
#pragma warning( disable : 4055 )			       /* conversion from data ptr to function ptr */
#endif
		glpi = (LPFN_GLPI) (void *) GetProcAddress(GetModuleHandle(TEXT("kernel32")), "GetLogicalProcessorInformation");

		if (glpi == NULL) {
			goto failed;
		}

		while (!done) {
			DWORD rc = glpi(buffer, &returnLength);
			if (FALSE == rc) {
				if (GetLastError() == ERROR_INSUFFICIENT_BUFFER) {
					if (buffer)
						free(buffer);
					buffer = (PSYSTEM_LOGICAL_PROCESSOR_INFORMATION) malloc(returnLength);

					if (buffer == NULL) {
						perror("zstd");
						exit(1);
					}
				} else {
					/*
					 * some other error 
					 */
					goto failed;
				}
			} else {
				done = TRUE;
			}
		}

		ptr = buffer;

		while (byteOffset + sizeof(SYSTEM_LOGICAL_PROCESSOR_INFORMATION) <= returnLength) {

			if (ptr->Relationship == RelationProcessorCore) {
				if (logical)
					numCores += CountSetBits(ptr->ProcessorMask);
				else
					numCores++;
			}

			ptr++;
			byteOffset += sizeof(SYSTEM_LOGICAL_PROCESSOR_INFORMATION);
		}

		free(buffer);

		return numCores;
	}

      failed:
	/*
	 * try to fall back on GetSystemInfo 
	 */
	{
		SYSTEM_INFO sysinfo;
		GetSystemInfo(&sysinfo);
		numCores = sysinfo.dwNumberOfProcessors;
		if (numCores == 0)
			numCores = 1;			       /* just in case */
	}
	return numCores;
}

#elif defined(__APPLE__)

#include <sys/sysctl.h>

/* Use apple-provided syscall
 * see: man 3 sysctl */
int UTIL_countCores(int logical)
{
	static S32 numCores = 0;			       /* apple specifies int32_t */
	if (numCores != 0)
		return numCores;

	{
		size_t size = sizeof(S32);
		int const ret = sysctlbyname(logical ? "hw.logicalcpu" : "hw.physicalcpu", &numCores, &size, NULL, 0);
		if (ret != 0) {
			if (errno == ENOENT) {
				/*
				 * entry not present, fall back on 1 
				 */
				numCores = 1;
			} else {
				perror("zstd: can't get number of cpus");
				exit(1);
			}
		}

		return numCores;
	}
}

#elif defined(__linux__)

/* parse /proc/cpuinfo
 * siblings / cpu cores should give hyperthreading ratio
 * otherwise fall back on sysconf */
int UTIL_countCores(int logical)
{
	static int numCores = 0;

	if (numCores != 0)
		return numCores;

	numCores = (int) sysconf(_SC_NPROCESSORS_ONLN);
	if (numCores == -1) {
		/*
		 * value not queryable, fall back on 1 
		 */
		return numCores = 1;
	}
	/*
	 * try to determine if there's hyperthreading 
	 */  {
		FILE *const cpuinfo = fopen("/proc/cpuinfo", "r");
#define BUF_SIZE 80
		char buff[BUF_SIZE];

		int siblings = 0;
		int cpu_cores = 0;
		int ratio = 1;

		if (cpuinfo == NULL) {
			/*
			 * fall back on the sysconf value 
			 */
			return numCores;
		}

		/*
		 * assume the cpu cores/siblings values will be constant across all * present processors 
		 */
		while (!feof(cpuinfo)) {
			if (fgets(buff, BUF_SIZE, cpuinfo) != NULL) {
				if (strncmp(buff, "siblings", 8) == 0) {
					const char *const sep = strchr(buff, ':');
					if (sep == NULL || *sep == '\0') {
						/*
						 * formatting was broken? 
						 */
						goto failed;
					}

					siblings = atoi(sep + 1);
				}
				if (strncmp(buff, "cpu cores", 9) == 0) {
					const char *const sep = strchr(buff, ':');
					if (sep == NULL || *sep == '\0') {
						/*
						 * formatting was broken? 
						 */
						goto failed;
					}

					cpu_cores = atoi(sep + 1);
				}
			} else if (ferror(cpuinfo)) {
				/*
				 * fall back on the sysconf value 
				 */
				goto failed;
			}
		}
		if (siblings && cpu_cores && siblings > cpu_cores) {
			ratio = siblings / cpu_cores;
		}

		if (ratio && numCores > ratio && !logical) {
			numCores = numCores / ratio;
		}

	      failed:
		fclose(cpuinfo);
		return numCores;
	}
}

#elif defined(__FreeBSD__)

#include <sys/param.h>
#include <sys/sysctl.h>

/* Use physical core sysctl when available
 * see: man 4 smp, man 3 sysctl */
int UTIL_countCores(int logical)
{
	static int numCores = 0;			       /* freebsd sysctl is native int sized */
#if __FreeBSD_version >= 1300008
	static int perCore = 1;
#endif
	if (numCores != 0)
		return numCores;

#if __FreeBSD_version >= 1300008
	{
		size_t size = sizeof(numCores);
		int ret = sysctlbyname("kern.smp.cores", &numCores, &size, NULL, 0);
		if (ret == 0) {
			if (logical) {
				ret = sysctlbyname("kern.smp.threads_per_core", &perCore, &size, NULL, 0);
				/*
				 * default to physical cores if logical cannot be read 
				 */
				if (ret == 0)
					numCores *= perCore;
			}
			return numCores;
		}
		if (errno != ENOENT) {
			perror("zstd: can't get number of cpus");
			exit(1);
		}
		/*
		 * sysctl not present, fall through to older sysconf method 
		 */
	}
#else
	/*
	 * suppress unused parameter warning 
	 */
	(void) logical;
#endif

	numCores = (int) sysconf(_SC_NPROCESSORS_ONLN);
	if (numCores == -1) {
		/*
		 * value not queryable, fall back on 1 
		 */
		numCores = 1;
	}
	return numCores;
}

#elif defined(__NetBSD__) || defined(__OpenBSD__) || defined(__DragonFly__) || defined(__CYGWIN__)

/* Use POSIX sysconf
 * see: man 3 sysconf */
int UTIL_countCores(int logical)
{
	static int numCores = 0;

	/*
	 * suppress unused parameter warning 
	 */
	(void) logical;

	if (numCores != 0)
		return numCores;

	numCores = (int) sysconf(_SC_NPROCESSORS_ONLN);
	if (numCores == -1) {
		/*
		 * value not queryable, fall back on 1 
		 */
		return numCores = 1;
	}
	return numCores;
}

#else

int UTIL_countCores(int logical)
{
	/*
	 * suppress unused parameter warning 
	 */
	(void) logical;

	/*
	 * assume 1 
	 */
	return 1;
}
#endif

int UTIL_countPhysicalCores(void)
{
	return UTIL_countCores(0);
}

int UTIL_countLogicalCores(void)
{
	return UTIL_countCores(1);
}
