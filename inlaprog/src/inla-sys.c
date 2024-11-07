
/* inla-sys.c
 * 
 * Copyright (C) 2007-2024 Havard Rue
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or (at
 * your option) any later version.
 *  * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 *
 * The author's contact information:
 *
 *        Haavard Rue
 *        CEMSE Division
 *        King Abdullah University of Science and Technology
 *        Thuwal 23955-6900, Saudi Arabia
 *        Email: haavard.rue@kaust.edu.sa
 *        Office: +966 (0)12 808 0640
 *
 */

int inla_ncpu(void)
{
#if defined(_SC_NPROCESSORS_ONLN)			       /* Linux, Solaris, AIX */
	return (int) sysconf(_SC_NPROCESSORS_ONLN);
#elif defined(__APPLE__)				       /* MacOSX */
	int count = -1;
	size_t size = sizeof(count);
	sysctlbyname("hw.ncpu", &count, &size, NULL, 0);
	return count;
#elif (defined(WIN32) || defined(WINDOWS))		       /* Windows */
	SYSTEM_INFO SystemInfo;
	GetSystemInfo(&SystemInfo);
	return SystemInfo.dwNumberOfProcessors;
#else
	return -1;
#endif
}

int inla_mkdir(const char *dirname)
{
#if defined(WINDOWS)
	return mkdir(dirname);
#else
	return mkdir(dirname, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
#endif
}

// this is from https://stackoverflow.com/questions/2513505/how-to-get-available-memory-c-g
// return RAM in Mb
#if defined(WIN32) || defined(WINDOWS)
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

#if defined(WINDOWS)
void inla_signal(int UNUSED(sig))
{
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


int inla_check_pardiso(void)
{
	// check if PARDISO-lib is installed and working
	if (GMRFLib_pardiso_check_install(1, 1) == GMRFLib_SUCCESS) {
		printf("SUCCESS: PARDISO IS INSTALLED AND WORKING\n");
		fflush(stdout);
	} else {
		printf("FAILURE: PARDISO IS NOT INSTALLED OR NOT WORKING\n");
		fflush(stdout);
		GMRFLib_pardiso_check_install(0, 0);
	}
	return GMRFLib_SUCCESS;
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
