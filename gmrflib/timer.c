
/* timer.c
 *
 * Copyright (C) 2001-2006 Havard Rue
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

/*!
  \file timer.c
  \brief Functions to measure (CPU) time and display statistics on the time spent in various
  routines and number times they are called.
 
  The normal usage is to use \c GMRFLib_cpu() to get the CPU time used since a fixed reference,
  and \c GMRFLib_timer_full_report() to display statistics of the most computational demanding
  routines in \c GMRFLib (if  \c GMRFLib_collect_timer_statistics = \c TRUE).
 
  \sa GMRFLib_collect_timer_statistics
*/

#include <stddef.h>
#include <math.h>
#include <assert.h>
#include <string.h>
#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#if !defined(__FreeBSD__)
#include <malloc.h>
#endif

#include "GMRFLib/GMRFLib.h"
#include "GMRFLib/GMRFLibP.h"
#include "GMRFLib/hashP.h"

#ifndef HGVERSION
#define HGVERSION
#endif
static const char RCSId[] = "file: " __FILE__ "  " HGVERSION;

/* Pre-hg-Id: $Id: timer.c,v 1.58 2010/03/16 22:39:25 hrue Exp $ */

static map_strvp *GMRFLib_timer_hashtable;

/* 
   if we have openmp, then use this function
 */
#if defined(_OPENMP)
#include <sys/time.h>
#include <omp.h>
double GMRFLib_cpu_default(void)
{
	static double ref = 0.0;
	if (!ref) {
		ref = omp_get_wtime();
	}
	return (omp_get_wtime() - ref);
}
#else

/* 
   else, choose default timer according to arch
*/
#if defined(__linux__) || defined(__linux) || defined(__sun__) || defined(__sun) || defined(__FreeBSD__) || defined(__FreeBSD) || defined(__APPLE__)
#include <sys/time.h>
#include <sys/resource.h>
#include <unistd.h>
double GMRFLib_cpu_default(void)
{
	if (1) {
		struct timeval time1;
		double time;
		static double ref = 0.0;

		gettimeofday(&time1, NULL);
		time = time1.tv_sec + time1.tv_usec * 1.0e-6;
		if (!ref) {
			ref = time;
		}
		return (time - ref);
	} else {
		struct rusage a;
		double time;
		static double ref = 0.0;
		
		getrusage(RUSAGE_SELF, &a);
		time = (double) ((a.ru_utime).tv_sec + (a.ru_stime).tv_sec) + (double) ((a.ru_utime).tv_usec + (a.ru_stime).tv_usec) * 1.0e-6;
		if (!ref) {
			ref = time;
		}
		return (time - ref);
	}
}
#else
#include <time.h>
double GMRFLib_cpu_default(void)
{
	static clock_t ref = 0;
	if (!ref) {
		ref = clock();
	}
	return (double) (clock() - ref) / (double) CLOCKS_PER_SEC;
}
#endif							       /* if defined(__linux__)... */
#endif							       /* if defined(_OPENMP)... */

/*
  here are functions for storing and reporting timing information
*/
int GMRFLib_timer_compare(const void *a, const void *b)
{
	const GMRFLib_timer_hashval_tp *aa, *bb;

	aa = (const GMRFLib_timer_hashval_tp *) (((const map_strvp_element *) a)->value);
	bb = (const GMRFLib_timer_hashval_tp *) (((const map_strvp_element *) b)->value);

	/*
	 * sort by name if both are empty 
	 */
	if (!aa->ntimes && !bb->ntimes)
		return strcmp(aa->name, bb->name);

	/*
	 * othewise, by average time spent 
	 */
	if (aa->ntimes && !bb->ntimes)
		return -1;
	if (bb->ntimes && !aa->ntimes)
		return 1;
	return (aa->ctime_acc / aa->ntimes > bb->ctime_acc / bb->ntimes ? -1 : 1);
}
int GMRFLib_timer_enter(const char *name)
{
//#pragma omp critical
	{
		GMRFLib_timer_hashval_tp *p;
		void *vpp;
		char *cname;

		cname = GMRFLib_strdup(name);

		if (!GMRFLib_timer_hashtable) {
			{
				if (!GMRFLib_timer_hashtable) {
					int i;
					map_strvp *tmp;

					tmp = Calloc(GMRFLib_MAX_THREADS, map_strvp);
					for (i = 0; i < GMRFLib_MAX_THREADS; i++) {
						map_strvp_init_hint(&tmp[i], 30);	/* about the number of elmements in the hash-table */
					}
					GMRFLib_timer_hashtable = tmp;
				}
			}
		}
		if ((vpp = map_strvp_ptr(&GMRFLib_timer_hashtable[omp_get_thread_num()], cname))) {
			p = *((GMRFLib_timer_hashval_tp **) vpp);
			Free(cname);
		} else {
			p = Calloc(1, GMRFLib_timer_hashval_tp);
			p->name = GMRFLib_strdup(name);
			map_strvp_set(&GMRFLib_timer_hashtable[omp_get_thread_num()], cname, (void *) p);
		}

		/*
		 * if ctime_ref > 0.0, then this routine is already initialized. in this case, we keep the first. 
		 */
		if (p->ctime_ref <= 0.0) {
			{
				if (p->ctime_ref <= 0.0) {
					p->ctime_ref = GMRFLib_cpu();
				}
			}
		}
	}
	return GMRFLib_SUCCESS;
}
int GMRFLib_timer_leave(const char *name)
{
//#pragma omp critical
	{
		GMRFLib_timer_hashval_tp *p;
		void *vpp;
		double used;
		char *cname;
		int ret = 0;

		cname = GMRFLib_strdup(name);

		if ((vpp = map_strvp_ptr(&GMRFLib_timer_hashtable[omp_get_thread_num()], cname))) {
			p = *((GMRFLib_timer_hashval_tp **) vpp);
			Free(cname);
		} else {
			/*
			 * this happen if _leave is executed but with no matching _enter, for example in the _EWRAP() macro. 
			 */
			Free(cname);
			ret = 1;
		}

		if (!ret) {
			if (p->ctime_ref < 0.0) {
				/*
				 * this is an ``illegal instruction''. _timer_leave is called without a corresponding call to _timer_enter.
				 * this happens, with purpose, with some of the `__intern' routines. 
				 */
			} else {
				used = GMRFLib_cpu() - p->ctime_ref;
				used = DMAX(0.0, used);	       /* yes */
				p->ctime_acc += used;
				p->ctime_acc2 += SQR(used);
				if (p->ntimes) {
					p->ctime_min = DMIN(p->ctime_min, used);
					p->ctime_max = DMAX(p->ctime_max, used);
				} else {
					p->ctime_min = p->ctime_max = used;
				}

				p->ctime_ref = -1.0;	       /* flag it specially */
				p->ntimes++;
			}
		}
	}
	return GMRFLib_SUCCESS;
}
const char *GMRFLib_timer_strip_store(const char *name)
{
	/*
	 * if name ends with `_store', then remove it from name.
	 * 
	 * this routine returns a ptr to a static storage containing the modified name, use 'strdup' alternatively. 
	 */
	size_t idx, len_name, len_suffix;
	const char *suffix = "_store";

#define STRLEN 2048
	static char nm[STRLEN + 1];

#pragma omp threadprivate(nm)

	if (!name) {
		return name;
	}

	len_name = strlen(name);
	GMRFLib_ASSERT_RETVAL(len_name <= STRLEN, GMRFLib_ESNH, name);
	len_suffix = strlen(suffix);

	strcpy(nm, name);
	idx = len_name - len_suffix;

	if (idx > 0 && (strcmp(&nm[idx], suffix) == 0)) {      /* yes, i want '> 0' */
		nm[idx] = '\0';
		return (const char *) nm;
	} else {
		return name;
	}
#undef STRLEN
}
const char *GMRFLib_timer_strip__intern(const char *name)
{
	/*
	 * if name ends with `__intern', then remove it from name.
	 * 
	 * this routine returns a ptr to a static storage containing the modified name, use 'strdup' alternatively. 
	 */
	size_t idx, len_name, len_suffix;
	const char *suffix = "__intern";

#define STRLEN 2048
	static char nm[STRLEN + 1];

#pragma omp threadprivate(nm)

	if (!name) {
		return name;
	}

	len_name = strlen(name);
	GMRFLib_ASSERT_RETVAL(len_name <= STRLEN, GMRFLib_ESNH, name);
	len_suffix = strlen(suffix);

	strcpy(nm, name);
	idx = len_name - len_suffix;

	if (idx > 0 && (strcmp(&nm[idx], suffix) == 0)) {      /* yes, i want '> 0' */
		nm[idx] = '\0';
		return (const char *) nm;
	} else {
		return name;
	}
#undef STRLEN
}
const char *GMRFLib_timer_strip(const char *name)
{
	const char *nm;

	nm = GMRFLib_timer_strip__intern(name);
	nm = GMRFLib_timer_strip_store(nm);
	nm = GMRFLib_timer_strip__intern(nm);

	return nm;
}

int GMRFLib_timer_print_entry(FILE * ffp, GMRFLib_timer_hashval_tp * p)
{
	fprintf(ffp, "%-41s %10.3f %6d %8.3f %8.3f %8.3f %8.3f\n",
		p->name, (p->ntimes ? p->ctime_acc / p->ntimes : 0.0),
		(int) p->ntimes, p->ctime_acc,
		(p->ntimes ? sqrt(DMAX(0.0, p->ctime_acc2 / p->ntimes - SQR(p->ctime_acc / p->ntimes))) : 0.0), p->ctime_min, p->ctime_max);

	return GMRFLib_SUCCESS;
}

/*!
  \brief Write statistics collected for the function named \c name to \c fp.
                                                                                                                  
  \param[in] fp Pointer to (an already open) file the report is written to.
  \param[in] name The name of the function for which statistics is to be reported. If \c name =
  \c NULL, then the statistics for all functions are displayed
                                                                                                                  
  \sa GMRFLib_timer_full_report
*/
int GMRFLib_timer_report(FILE * fp, const char *name)
{
	GMRFLib_timer_hashval_tp *p;
	void *vpp;
	FILE *ffp;
	const char *sep = "-----------------------------------------------------------------------------------------------";
	char *cname;
	int k;

	if (!GMRFLib_timer_hashtable) {
		return GMRFLib_SUCCESS;
	}
	ffp = (fp ? fp : stdout);

	for (k = 0; k < GMRFLib_MAX_THREADS; k++) {
		fprintf(ffp, "\n\nGMRFLib report on time usage for thread %1d\n%-41s %10s %6s %8s %8s %8s %8s\n%s\n",
			k, "Function", "Mean", "N", "Total", "Stdev", "Min", "Max", sep);
		if (name) {
			cname = GMRFLib_strdup(name);
			if ((vpp = map_strvp_ptr(&GMRFLib_timer_hashtable[k], cname))) {
				p = *((GMRFLib_timer_hashval_tp **) vpp);
				GMRFLib_timer_print_entry(ffp, p);
			}
			Free(cname);
		} else {
			map_strvp_element *all;
			mapkit_size_t count, i;

			map_strvp_getall(&GMRFLib_timer_hashtable[k], &all, &count);
			qsort(all, (size_t) count, sizeof(map_strvp_element), GMRFLib_timer_compare);
			for (i = 0; i < count; i++) {
				GMRFLib_timer_print_entry(ffp, (GMRFLib_timer_hashval_tp *) (all[i].value));
			}
			free(all);
			all = NULL;
		}
		fprintf(ffp, "%s\n", sep);
	}
	return GMRFLib_SUCCESS;
}

/*!
  \brief Write  statistics collected for all functions to \c fp.
                                                                                                                  
  This is the routine to use to print internal statistics of number of times each routine is called
  and its CPU usage.
                                                                                                                  
  \remark This function is simply a wrapper for \c GMRFLib_timer_report(fp,NULL).
                                                                                                                  
  \sa  GMRFLib_collect_timer_statistics
*/
int GMRFLib_timer_full_report(FILE * fp)
{
	if (!GMRFLib_timer_hashtable) {
		return GMRFLib_SUCCESS;
	}
	return GMRFLib_timer_report(fp, NULL);
}
void GMRFLib_timer_report__signal(int sig)
{
	/*
	 * a version to be installed using signal(...) 
	 */
	GMRFLib_timer_report(stdout, NULL);
	return;
}
