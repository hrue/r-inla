
/* timer.c
 *
 * Copyright (C) 2001-2022 Havard Rue
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
#if !defined(__FreeBSD__)
#include <malloc.h>
#endif

#include "GMRFLib/GMRFLib.h"
#include "GMRFLib/GMRFLibP.h"
#include "GMRFLib/hashP.h"

#ifndef GITCOMMIT
#define GITCOMMIT
#endif
static const char GitID[] = "file: " __FILE__ "  " GITCOMMIT;

static map_strvp *GMRFLib_timer_hashtable = NULL;

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
	return (aa->ctime_acc > bb->ctime_acc ? -1 : 1);
}

int GMRFLib_timer_init(void)
{
	return GMRFLib_timer_enter(NULL);
}

int GMRFLib_timer_enter(const char *name)
{
	GMRFLib_timer_hashval_tp *p;
	void *vpp;
	char *cname;

	cname = GMRFLib_strdup(name);
	if (!GMRFLib_timer_hashtable) {
		if (!GMRFLib_timer_hashtable) {
			int i;
			map_strvp *tmp;

			tmp = Calloc(GMRFLib_MAX_THREADS(), map_strvp);
			for (i = 0; i < GMRFLib_MAX_THREADS(); i++) {
				map_strvp_init_hint(&tmp[i], 30);	/* about the number of elmements in the hash-table */
			}
			GMRFLib_timer_hashtable = tmp;
		}
	}
	if (!name) {
		return GMRFLib_SUCCESS;
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
	return GMRFLib_SUCCESS;
}

int GMRFLib_timer_leave(const char *name)
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
			used = DMAX(0.0, used);		       /* yes */
			p->ctime_acc += used;
			p->ctime_acc2 += SQR(used);
			if (p->ntimes) {
				p->ctime_min = DMIN(p->ctime_min, used);
				p->ctime_max = DMAX(p->ctime_max, used);
			} else {
				p->ctime_min = p->ctime_max = used;
			}

			p->ctime_ref = -1.0;		       /* flag it specially */
			p->ntimes++;
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

int GMRFLib_timer_print_entry(FILE * ffp, GMRFLib_timer_hashval_tp * p, double total_time)
{
	char *nm = (strncmp(p->name, "GMRFLib_", (size_t) 8L) == 0 ? p->name + 8L : p->name);

	fprintf(ffp, "%-39s %10.1f %4.1f%% %10d %6.2f %6.2f %6.2f %6.2f\n",
		nm, p->ctime_acc, (total_time > 0.0 ? p->ctime_acc / total_time * 100.0 : 0.0),
		(int) p->ntimes, (p->ntimes ? p->ctime_acc / p->ntimes : 0.0),
		(p->ntimes ? sqrt(DMAX(0.0, p->ctime_acc2 / p->ntimes - SQR(p->ctime_acc / p->ntimes))) : 0.0), p->ctime_min, p->ctime_max);

	return GMRFLib_SUCCESS;
}

int GMRFLib_timer_report_OLD(FILE * fp)
{
	FILE *ffp;
	const char *sep = "--------------------------------------------------------------------------------------------------";
	int k;

	if (!GMRFLib_timer_hashtable) {
		return GMRFLib_SUCCESS;
	}
	ffp = (fp ? fp : stdout);

	for (k = 0; k < GMRFLib_MAX_THREADS(); k++) {
		fprintf(ffp, "\n\nGMRFLib report on time usage for thread %1d\n%-42s   %11s %10s %6s %6s %6s %6s\n%s\n",
			k, "Function", "Total(s) %", "N", "Mean", "Stdev", "Min", "Max", sep);
		if (1) {
			map_strvp_element *all;
			mapkit_size_t count, i;

			map_strvp_getall(&GMRFLib_timer_hashtable[k], &all, &count);
			qsort(all, (size_t) count, sizeof(map_strvp_element), GMRFLib_timer_compare);
			for (i = 0; i < count; i++) {
				GMRFLib_timer_print_entry(ffp, (GMRFLib_timer_hashval_tp *) (all[i].value), -1.0);
			}
			free(all);
			all = NULL;
		}
		fprintf(ffp, "%s\n", sep);
	}
	return GMRFLib_SUCCESS;
}

int GMRFLib_timer_report(FILE * fp)
{
	// this function is used just before exit(), hence its ok that this function leak
	fp = (fp ? fp : stdout);

	char **namelist = NULL;
	int namelist_len = 0;
	int namelist_len_max = 1024;
	int k, j, found;
	GMRFLib_timer_hashval_tp *val;
	map_strvp_element *all;
	mapkit_size_t count, i;

	namelist = Calloc(namelist_len_max, char *);
	for (k = 0; k < GMRFLib_MAX_THREADS(); k++) {

		map_strvp_getall(&GMRFLib_timer_hashtable[k], &all, &count);
		if (all) {
			for (i = 0; i < count; i++) {
				val = (GMRFLib_timer_hashval_tp *) (all[i].value);
				for (j = found = 0; j < namelist_len; j++) {
					if (!strcmp(namelist[j], val->name)) {
						found = 1;
						break;
					}
				}
				if (!found) {
					assert(namelist_len < namelist_len_max);
					namelist[namelist_len] = strdup(val->name);
					namelist_len++;
				}
			}
			Free(all);
		}
	}

	map_strvp *total = Calloc(1, map_strvp);
	map_strvp_init_hint(total, namelist_len);

	double total_acc_time = 0.0;
	for (int ilist = 0; ilist < namelist_len; ilist++) {
		char *nm = namelist[ilist];
		GMRFLib_timer_hashval_tp *value = Calloc(1, GMRFLib_timer_hashval_tp);

		value->name = strdup(nm);
		for (k = 0; k < GMRFLib_MAX_THREADS(); k++) {
			map_strvp_element *all;
			mapkit_size_t count, i;

			map_strvp_getall(&GMRFLib_timer_hashtable[k], &all, &count);
			if (all) {
				for (i = 0, found = 0; i < count && !found; i++) {
					val = (GMRFLib_timer_hashval_tp *) (all[i].value);
					if (!strcmp(nm, val->name)) {
						total_acc_time += val->ctime_acc;
						value->ctime_acc += val->ctime_acc;
						value->ctime_acc2 += val->ctime_acc2;
						value->ctime_min = (value->ctime_min > 0.0 ?
								    DMIN(value->ctime_min, val->ctime_min) : val->ctime_min);
						value->ctime_max = DMAX(value->ctime_max, val->ctime_max);
						value->ntimes += val->ntimes;
						found = 1;
					}
				}
			}
		}
		map_strvp_set(total, value->name, value);
	}

	const char *sep = "-----------------------------------------------------------------------------------------------";

	fprintf(fp, "\n%-42s   %11s %10s %6s %6s %6s %6s\n%s\n", "Function", "Total(s) %", "N", "Mean", "Stdev", "Min", "Max", sep);

	map_strvp_getall(total, &all, &count);
	qsort(all, (size_t) count, sizeof(map_strvp_element), GMRFLib_timer_compare);
	for (i = 0; i < count; i++) {
		GMRFLib_timer_print_entry(fp, (GMRFLib_timer_hashval_tp *) (all[i].value), total_acc_time);
	}
	fprintf(fp, "%s\n", sep);

	Free(namelist);
	Free(total);

	return GMRFLib_SUCCESS;
}

int GMRFLib_timer_full_report(FILE * fp)
{
	if (!GMRFLib_timer_hashtable) {
		return GMRFLib_SUCCESS;
	}
	return GMRFLib_timer_report(fp);
}

void GMRFLib_timer_report__signal(int UNUSED(sig))
{
	/*
	 * a version to be installed using signal(...) 
	 */
	GMRFLib_timer_report(stdout);
	return;
}

int GMRFLib_timer_exit(void)
{
	if (GMRFLib_timer_hashtable) {
		for (int i = 0; i < GMRFLib_MAX_THREADS(); i++) {
			map_strvp_free(&(GMRFLib_timer_hashtable[i]));
		}
		Free(GMRFLib_timer_hashtable);
		GMRFLib_timer_hashtable = NULL;
	}
	return GMRFLib_SUCCESS;
}
