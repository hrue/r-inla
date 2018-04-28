
/* openmp.c
 * 
 * Copyright (C) 2007-2018 Havard Rue
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

int GMRFLib_openmp_implement_strategy(GMRFLib_openmp_place_tp place, void *arg, GMRFLib_smtp_tp * smtp)
{
	int nt;
	int ntmax = GMRFLib_MAX_THREADS;
	int ntmax_inner, ntmax_outer;
	int strategy = (GMRFLib_openmp ? GMRFLib_openmp->strategy : GMRFLib_OPENMP_STRATEGY_MEDIUM);
	int nested = 0;
	int *nhyper = (int *) arg;

	// once only
	static int pardiso_ok = -1;
	if (pardiso_ok < 0) {
		pardiso_ok = (GMRFLib_pardiso_check_install(1, 1) == GMRFLib_SUCCESS ? 1 : 0);
	}

	static GMRFLib_smtp_tp smtp_store = GMRFLib_SMTP_INVALID;
	if (smtp) {
		smtp_store = *smtp;
	}

	P(ntmax);
	
	if (pardiso_ok && smtp_store == GMRFLib_SMTP_PARDISO) {
		FIXME1("PARDISO is installed and working, set openmp->strategy = 'PARDISO'");
		if (strategy == GMRFLib_OPENMP_STRATEGY_DEFAULT) {
			strategy = GMRFLib_OPENMP_STRATEGY_PARDISO_SERIAL;
		} else {
			assert(strategy == GMRFLib_OPENMP_STRATEGY_PARDISO_SERIAL ||
			       strategy == GMRFLib_OPENMP_STRATEGY_PARDISO_PARALLEL);
		}
	}

	if (smtp_store == GMRFLib_SMTP_PARDISO) {
		switch(strategy) {
		case GMRFLib_OPENMP_STRATEGY_PARDISO_SERIAL:
			ntmax_outer = ntmax;
			ntmax_inner = 1;
			break;
		case GMRFLib_OPENMP_STRATEGY_PARDISO_PARALLEL:
			ntmax_outer = 1;
			ntmax_inner = ntmax;
			break;
		default:
			assert(0 == 1);
		}
	} else {
		ntmax_outer = ntmax;
		ntmax_inner = 1;
	}

	// default values
	nt = ntmax;
	nested = 0;
	GMRFLib_openmp->max_threads_inner = ntmax_inner;
	GMRFLib_openmp->max_threads_outer = ntmax_outer;

	switch (place) {
	case GMRFLib_OPENMP_PLACES_BUILD_MODEL:
		switch (strategy) {
		case GMRFLib_OPENMP_STRATEGY_SMALL:
			nt = 1;
			break;
		case GMRFLib_OPENMP_STRATEGY_MEDIUM:
			nt = 1;
			break;
		case GMRFLib_OPENMP_STRATEGY_LARGE:
		case GMRFLib_OPENMP_STRATEGY_DEFAULT:
			break;
		case GMRFLib_OPENMP_STRATEGY_HUGE:
		case GMRFLib_OPENMP_STRATEGY_PARDISO_SERIAL: 
		case GMRFLib_OPENMP_STRATEGY_PARDISO_PARALLEL: 
			GMRFLib_openmp->max_threads_outer = ntmax;
			break;
		case GMRFLib_OPENMP_STRATEGY_NONE:
		default:
			assert(0 == 1);
		}
		break;

	case GMRFLib_OPENMP_PLACES_OPTIMIZE:
		switch (strategy) {
		case GMRFLib_OPENMP_STRATEGY_SMALL:
			nt = 1;
			break;
		case GMRFLib_OPENMP_STRATEGY_MEDIUM:
			nt = 1;
			break;
		case GMRFLib_OPENMP_STRATEGY_LARGE:
		case GMRFLib_OPENMP_STRATEGY_DEFAULT:
			nt = IMIN(*nhyper + 1, ntmax);
			break;
		case GMRFLib_OPENMP_STRATEGY_HUGE:
			nt = IMIN(*nhyper + 1, ntmax);
			break;
		case GMRFLib_OPENMP_STRATEGY_PARDISO_SERIAL: 
		case GMRFLib_OPENMP_STRATEGY_PARDISO_PARALLEL: 
			nt = ntmax;
			break;
		case GMRFLib_OPENMP_STRATEGY_NONE:
		default:
			assert(0 == 1);
		}
		break;

	case GMRFLib_OPENMP_PLACES_HESSIAN:
		switch (strategy) {
		case GMRFLib_OPENMP_STRATEGY_SMALL:
			nt = 1;
			break;
		case GMRFLib_OPENMP_STRATEGY_MEDIUM:
			nt = ntmax;
			break;
		case GMRFLib_OPENMP_STRATEGY_LARGE:
		case GMRFLib_OPENMP_STRATEGY_DEFAULT:
			nt = ntmax;
			break;
		case GMRFLib_OPENMP_STRATEGY_HUGE:
			nt = ntmax;
			break;
		case GMRFLib_OPENMP_STRATEGY_PARDISO_SERIAL: 
		case GMRFLib_OPENMP_STRATEGY_PARDISO_PARALLEL: 
			break;
		case GMRFLib_OPENMP_STRATEGY_NONE:
		default:
			assert(0 == 1);
		}
		break;

	case GMRFLib_OPENMP_PLACES_HESSIAN_SCALE:
		switch (strategy) {
		case GMRFLib_OPENMP_STRATEGY_SMALL:
			nt = 1;
			break;
		case GMRFLib_OPENMP_STRATEGY_MEDIUM:
			break;
		case GMRFLib_OPENMP_STRATEGY_LARGE:
		case GMRFLib_OPENMP_STRATEGY_DEFAULT:
			break;
		case GMRFLib_OPENMP_STRATEGY_HUGE:
			break;
		case GMRFLib_OPENMP_STRATEGY_PARDISO_SERIAL: 
		case GMRFLib_OPENMP_STRATEGY_PARDISO_PARALLEL: 
			break;
		case GMRFLib_OPENMP_STRATEGY_NONE:
		default:
			assert(0 == 1);
		}
		break;

	case GMRFLib_OPENMP_PLACES_INTEGRATE:
		switch (strategy) {
		case GMRFLib_OPENMP_STRATEGY_SMALL:
			break;
		case GMRFLib_OPENMP_STRATEGY_MEDIUM:
			break;
		case GMRFLib_OPENMP_STRATEGY_LARGE:
		case GMRFLib_OPENMP_STRATEGY_DEFAULT:
			break;
		case GMRFLib_OPENMP_STRATEGY_HUGE:
			break;
		case GMRFLib_OPENMP_STRATEGY_PARDISO_SERIAL: 
		case GMRFLib_OPENMP_STRATEGY_PARDISO_PARALLEL: 
			break;
		case GMRFLib_OPENMP_STRATEGY_NONE:
		default:
			assert(0 == 1);
		}
		break;

	case GMRFLib_OPENMP_PLACES_INTEGRATE_HYPERPAR:
		switch (strategy) {
		case GMRFLib_OPENMP_STRATEGY_SMALL:
			nt = (*nhyper > 2 ? ntmax : 1);
			break;
		case GMRFLib_OPENMP_STRATEGY_MEDIUM:
			nt = *nhyper;
			break;
		case GMRFLib_OPENMP_STRATEGY_LARGE:
		case GMRFLib_OPENMP_STRATEGY_DEFAULT:
			nt = *nhyper;
			break;
		case GMRFLib_OPENMP_STRATEGY_HUGE:
		case GMRFLib_OPENMP_STRATEGY_PARDISO_SERIAL: 
		case GMRFLib_OPENMP_STRATEGY_PARDISO_PARALLEL: 
			GMRFLib_openmp->max_threads_outer = ntmax;
			break;
		case GMRFLib_OPENMP_STRATEGY_NONE:
		default:
			assert(0 == 1);
		}
		break;

	case GMRFLib_OPENMP_PLACES_COMBINE:
		switch (strategy) {
		case GMRFLib_OPENMP_STRATEGY_SMALL:
		case GMRFLib_OPENMP_STRATEGY_MEDIUM:
		case GMRFLib_OPENMP_STRATEGY_LARGE:
		case GMRFLib_OPENMP_STRATEGY_DEFAULT:
		case GMRFLib_OPENMP_STRATEGY_HUGE:
		case GMRFLib_OPENMP_STRATEGY_PARDISO_SERIAL: 
		case GMRFLib_OPENMP_STRATEGY_PARDISO_PARALLEL: 
			GMRFLib_openmp->max_threads_outer = ntmax;
			break;
		case GMRFLib_OPENMP_STRATEGY_NONE:
		default:
			assert(0 == 1);
		}
		break;

	case GMRFLib_OPENMP_PLACES_DEFAULT:
		switch (strategy) {
		case GMRFLib_OPENMP_STRATEGY_SMALL:
			nt = 1;
			break;
		case GMRFLib_OPENMP_STRATEGY_MEDIUM:
			nt = 1;
			break;
		case GMRFLib_OPENMP_STRATEGY_LARGE:
		case GMRFLib_OPENMP_STRATEGY_DEFAULT:
			break;
		case GMRFLib_OPENMP_STRATEGY_HUGE:
			break;
		case GMRFLib_OPENMP_STRATEGY_PARDISO_SERIAL: 
		case GMRFLib_OPENMP_STRATEGY_PARDISO_PARALLEL: 
			break;
		case GMRFLib_OPENMP_STRATEGY_NONE:
		default:
			assert(0 == 1);
		}
		break;

	case GMRFLib_OPENMP_PLACES_NONE:
		switch (strategy) {
		case GMRFLib_OPENMP_STRATEGY_NONE:	       /* only one option allowed */
			break;
		case GMRFLib_OPENMP_STRATEGY_SMALL:
		case GMRFLib_OPENMP_STRATEGY_MEDIUM:
		case GMRFLib_OPENMP_STRATEGY_LARGE:
		case GMRFLib_OPENMP_STRATEGY_DEFAULT:
		case GMRFLib_OPENMP_STRATEGY_HUGE:
		case GMRFLib_OPENMP_STRATEGY_PARDISO_SERIAL: 
		case GMRFLib_OPENMP_STRATEGY_PARDISO_PARALLEL: 
		default:
			assert(0 == 1);
		}
		break;
	}

	FIXME("SET nt = inla_ncpu()");
	nt = IMAX(1, IMIN(ntmax, nt));
	nt = inla_ncpu();
	
	omp_set_num_threads(nt);
	omp_set_nested(nested);

	P(GMRFLib_openmp->max_threads_inner);
	P(GMRFLib_openmp->max_threads_outer);

	return GMRFLib_SUCCESS;
}
