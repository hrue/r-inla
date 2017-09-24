
/* my.c
 * 
 * Copyright (C) 2016-2016 Havard Rue
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
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
#ifndef HGVERSION
#define HGVERSION
#endif

static const char RCSId[] = "file: " __FILE__ "  " HGVERSION;

#include <string.h>
#include <stdio.h>
#if !defined(__FreeBSD__)
#include <malloc.h>
#endif
#include <stdlib.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

#include "inla.h"
#include "my.h"
#include "GMRFLib/GMRFLibP.h"
#include "GMRFLib/GMRFLib.h"

int my_file_exists(const char *filename)
{
	struct stat sb;
	return ((stat(filename, &sb) == 0 && S_ISREG(sb.st_mode)) ? INLA_OK : !INLA_OK);
}
int my_dir_exists(const char *dirname)
{
	struct stat sb;
	return ((stat(dirname, &sb) == 0 && S_ISDIR(sb.st_mode)) ? INLA_OK : !INLA_OK);
}
int my_setenv(char *str, int prefix)
{
	/*
	 * set a variable in the enviroment; if PREFIX prepend with inla_, so that a=b yields inla_a=b. 
	 */
	char *p = NULL, *var = NULL;
	int debug = 0;

	if (debug)
		printf("enter my_setenv with [%s]\n", str);

	p = strchr(str, '=');
	if (!p) {
		fprintf(stderr, "*** Error: Argument is void [%s]\n", str);
		exit(EXIT_FAILURE);
	}
	*p = '\0';
#if defined(WINDOWS)
	if (prefix) {
		GMRFLib_sprintf(&var, "inla_%s=%s", str, p + 1);
	} else {
		GMRFLib_sprintf(&var, "%s=%s", str, p + 1);
	}
	putenv(var);
	if (debug)
		printf("putenv \t%s\n", var);
#else
	if (prefix) {
		GMRFLib_sprintf(&var, "inla_%s", str);
	} else {
		GMRFLib_sprintf(&var, "%s", str);
	}
	setenv(var, p + 1, 1);
	if (debug)
		printf("\tsetenv %s=%s\n", var, p + 1);
#endif
	Free(var);
	return INLA_OK;
}
