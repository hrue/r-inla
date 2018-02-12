
/* smtp-pardiso.h
 * 
 * Copyright (C) 2018 Havard Rue & Alexander Litvinenko
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
  \file smtp-band.h
  \brief GMRFLib interface to the band-solver in the LAPACK library
*/

#ifndef __GMRFLib_SMTP_PARDISO_H__
#define __GMRFLib_SMTP_PARDISO_H__

#if !defined(__FreeBSD__)
#include <malloc.h>
#endif
#include <stdlib.h>

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
   
 */

typedef struct {
	int base;					       /* 0 or 1 */
	int n;
	int na;
	int *ia;
	int *ja;
	double *a;
} GMRFLib_csr_tp;

typedef enum {
	GMRFLib_PARDISO_FLAG_REORDER=1,
	GMRFLib_PARDISO_FLAG_SYMFACT,
	GMRFLib_PARDISO_FLAG_CHOL,
	GMRFLib_PARDISO_FLAG_INV,
	GMRFLib_PARDISO_FLAG_SOLVE
} GMRFLib_pardiso_flag_tp;


typedef struct {
	void *pt[64];
	int iparm[64];
	double dparm[64];

	int init_done;
	
	char *var;					       /* Auxiliary variables. */
	double dummy;
	int err_code;
	int idummy;
	int maxfct;
	int mnum;
	int msglvl;
	int mtype;
	int nrhs;
	int num_procs;					       /* Number of processors. */
	int phase;
	int solver;
	int *my_perm;
	int *perm;
	int L_nnz;
	double log_det_Q;
	GMRFLib_csr_tp *L;
	GMRFLib_csr_tp *Q;
	GMRFLib_csr_tp *Qinv;
} GMRFLib_pardiso_store_tp;

#define STDOUT_TO_DEV_NULL_START(_silent)			\
	int XX_stdout_dupfd;					\
	int XX_silent = _silent;				\
	FILE *XX_temp_out;					\
	if (XX_silent) {					\
		XX_stdout_dupfd = dup(1);			\
		XX_temp_out = fopen("/dev/null", "w");		\
		dup2(fileno(XX_temp_out), 1);			\
	}
	

#define STDOUT_TO_DEV_NULL_END				\
	if (XX_silent) {				\
		fflush(stdout);				\
		fclose(XX_temp_out);			\
		dup2(XX_stdout_dupfd, 1);		\
		close(XX_stdout_dupfd);			\
	}	


int GMRFLib_Q2csr(GMRFLib_csr_tp ** csr, GMRFLib_graph_tp * graph, GMRFLib_Qfunc_tp * Qfunc, void *Qfunc_arg);
int GMRFLib_print_csr(FILE * fp, GMRFLib_csr_tp * csr);
int GMRFLib_csr2Q(GMRFLib_tabulate_Qfunc_tp ** Qtab, GMRFLib_graph_tp ** graph, GMRFLib_csr_tp * csr);
int GMRFLib_duplicate_csr(GMRFLib_csr_tp ** csr_to, GMRFLib_csr_tp * csr_from);
int GMRFLib_free_csr(GMRFLib_csr_tp ** csr);
int GMRFLib_pardiso_setparam(GMRFLib_pardiso_flag_tp flag, GMRFLib_pardiso_store_tp * store);
int GMRFLib_pardiso_check_install(int quiet);
int GMRFLib_csr_base(int base, GMRFLib_csr_tp *M);
int GMRFLib_csr_convert(GMRFLib_csr_tp *M);
int GMRFLib_Q2csr_check(GMRFLib_csr_tp *M);
int GMRFLib_pardiso_init(GMRFLib_pardiso_store_tp *store);
int GMRFLib_pardiso_symfact(GMRFLib_pardiso_store_tp *store);
int GMRFLib_pardiso_chol(GMRFLib_pardiso_store_tp *store, GMRFLib_graph_tp *graph,
			 GMRFLib_Qfunc_tp *Qfunc, void *Qfunc_arg);

double GMRFLib_pardiso_Qfunc_default(int i, int j, void *arg);
int GMRFLib_pardiso_reorder(GMRFLib_pardiso_store_tp *store, GMRFLib_graph_tp *graph, int *reordering);


void pardisoinit(void *, int *, int *, int *, double *, int *);
void pardiso(void *, int *, int *, int *, int *, int *, double *, int *, int *, int *, int *, int *, int *, double *, double *, int *, double *);
void pardiso_chkmatrix(int *, int *, double *, int *, int *, int *);
void pardiso_chkvec(int *, int *, double *, int *);
void pardiso_printstats(int *, int *, double *, int *, int *, int *, double *, int *);
void pardiso_get_factor_csc(void **, double *, int *, int *, double *, int *, int *, int *, int);
void pardiso_get_inverse_factor_csc(void **, double *, int *, int *, int *, int);


__END_DECLS
#endif
