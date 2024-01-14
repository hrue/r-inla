
/* timer.c
 *
 * Copyright (C) 2001-2024 Havard Rue
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

#include <stddef.h>
#include <math.h>
#include <assert.h>
#include <string.h>
#include <stdio.h>
#include <time.h>
#include <stdlib.h>

#include "GMRFLib/GMRFLib.h"
#include "GMRFLib/GMRFLibP.h"
#include "GMRFLib/hashP.h"


#if defined(WINDOWS)

// define the function in the high-prec-timer.cpp

#else							       // if defined(WINDOWS)

#if defined(_OPENMP)

#include <sys/time.h>
#include <omp.h>
double GMRFLib_cpu_default(void)
{
	return (omp_get_wtime());
}

#else							       // if defined(_OPENMP)

#include <sys/time.h>
#include <sys/resource.h>
#include <unistd.h>
double GMRFLib_cpu_default(void)
{
	struct timeval time1;
	double time;
	gettimeofday(&time1, NULL);
	time = time1.tv_sec + time1.tv_usec * 1.0e-6;
	return (time);
}

#endif							       // if defined(_OPENMP)
#endif							       // if defined(WINDOWS)
