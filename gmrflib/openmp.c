
/* openmp.c
 * 
 * Copyright (C) 2007 Havard Rue
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
 *       H{\aa}vard Rue
 *       Department of Mathematical Sciences
 *       The Norwegian University of Science and Technology
 *       N-7491 Trondheim, Norway
 *       Voice: +47-7359-3533    URL  : http://www.math.ntnu.no/~hrue  
 *       Fax  : +47-7359-3524    Email: havard.rue@math.ntnu.no
 *
 */

#ifndef HGVERSION
#define HGVERSION
#endif
static const char RCSId[] =  "file: " __FILE__ "  " HGVERSION; 
/* Pre-hg-Id: $Id: openmp.c,v 1.7 2007/07/17 05:55:02 hrue Exp $ */ 

#if !defined(__FreeBSD__)
#include <malloc.h>
#endif
#include <stdlib.h>

#include "GMRFLib/GMRFLib.h"

/**
 *  \file openmp.c
 *  \brief This file contains some OpenMP support files for GMRFLib.
 *
 *  The #GMRFLib_rng object containing the RNG-structure, is safe to use in any threaded environment. Its contents is
 *  threadprivate (one object for each thread), and is automatically initialised in each thread. If a spesific seed is required
 *  (\c GMRFLib_rng_init()), then this must be done for each thread. The used seed, is still available by \c GMRFLib_rng_seed,
 *  which is threadprivate.
 *
 */

/* 
   just implement dummy functions for a serial version without openmp-support.
 */

#ifndef _OPENMP

void omp_set_num_threads(int num)
{
	num = 1;
	return;
}
int omp_get_num_threads(void)
{
	return 1;
}
int omp_get_max_threads(void)
{
	return 1;
}
int omp_get_thread_num(void)
{
	return 0;
}
int omp_get_thread_num_(void)
{							       /* used from fortran */
	return 0;
}
int omp_get_num_procs(void)
{
	return 1;
}
int omp_in_parallel(void)
{
	return 0;
}
void omp_set_dynamic(int val)
{
	val = 0;
	return;
}
int omp_get_dynamic(void)
{
	return 0;
}
void omp_set_nested(int val)
{
	val = 0;
	return;
}
int omp_get_nested(void)
{
	return 0;
}
double omp_get_wtime(void)
{
	return GMRFLib_cpu();
}
double omp_get_wtick(void)
{
	return 0.0;
}

#endif
