
/* timer.h
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
 *
 */

/*! 
  \file timer.h
  \brief Typedefs and defines for the timer-utility in GMRFLib.
*/

#ifndef __GMRFLib_TIMER_H__
#define __GMRFLib_TIMER_H__

#include <time.h>
#include <stdlib.h>
#if !defined(__FreeBSD__)
#include <malloc.h>
#endif

#undef __BEGIN_DECLS
#undef __END_DECLS
#ifdef __cplusplus
#define __BEGIN_DECLS extern "C" {
#define __END_DECLS }
#else
#define __BEGIN_DECLS					       /* empty */
#define __END_DECLS					       /* empty */
#endif

__BEGIN_DECLS

/* 
   dummy comment
 */
    typedef struct {
	char *name;					       /* function name */
	double ntimes;					       /* number of times called */
	double ctime_acc;				       /* accumulated ctime */
	double ctime_ref;				       /* last reference time */

	double ctime_min;				       /* minimum time used */
	double ctime_max;				       /* maximum times used */
	double ctime_acc2;				       /* accumulated ctime^2 */
} GMRFLib_timer_hashval_tp;

/*!
  \brief Macro to be placed at the entry point of each routine which CPU time is be monitored.
*/
#define GMRFLib_ENTER_ROUTINE if (1) { if(0) printf("==>Enter [%s]\n", __GMRFLib_FuncName); if (GMRFLib_collect_timer_statistics) GMRFLib_timer_enter(GMRFLib_timer_strip(__GMRFLib_FuncName));}

/*!
  \brief Macro to be placed at \em each exit point of each routine which CPU time is be monitored.
*/
#define GMRFLib_LEAVE_ROUTINE if (1) { if (0) printf("==>    Leave [%s]\n", __GMRFLib_FuncName); if (GMRFLib_collect_timer_statistics) GMRFLib_timer_leave(GMRFLib_timer_strip(__GMRFLib_FuncName));}

double GMRFLib_cpu_default(void);
int GMRFLib_timer_table_expand(void);
int GMRFLib_timer_find_entry(const char *name);
int GMRFLib_timer_compare(const void *a, const void *b);
int GMRFLib_timer_table_sort(void);
int GMRFLib_timer_enter(const char *name);
int GMRFLib_timer_leave(const char *name);
int GMRFLib_timer_print_entry(FILE * ffp, GMRFLib_timer_hashval_tp * p);
int GMRFLib_timer_report(FILE * fp, const char *name);
int GMRFLib_timer_full_report(FILE * fp);
const char *GMRFLib_timer_strip_store(const char *name);
const char *GMRFLib_timer_strip__intern(const char *name);
const char *GMRFLib_timer_strip(const char *name);
void GMRFLib_timer_report__signal(int sig);

__END_DECLS
#endif
