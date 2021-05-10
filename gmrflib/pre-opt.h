
/* pre-opt.h
 * 
 * Copyright (C) 2021-2021 Havard Rue
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

#ifndef __GMRFLib_PREOPT_H__
#define __GMRFLib_PREOPT_H__

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
 * 
 */

typedef enum {
	GMRFLib_PREOPTM_TP_F = 1,
	GMRFLib_PREOPTM_TP_BETA,
	GMRFLib_PREOPTM_TP___VOID = -1
} GMRFLib_preoptm_type_types_tp;

typedef struct {
	GMRFLib_preoptm_type_types_tp tp;
	int idx;
	int tp_idx;
} GMRFLib_preoptm_type_tp;

typedef struct {
	int n;						       /* length of the linear predictor */

	int *idx_map_f;
	int *idx_map_beta;

	GMRFLib_tabulate_Qfunc_tp *latent_Qtab;
	GMRFLib_graph_tp *latent_graph;
	
	int nf;
	GMRFLib_graph_tp **f_graph;
	GMRFLib_Qfunc_tp **f_Qfunc;
	GMRFLib_Qfunc_tp ***ff_Qfunc;			       /* interaction */
	void **f_Qfunc_arg;
	void ***ff_Qfunc_arg;

	int nbeta;
	double **covariate;
	double *prior_precision;

	GMRFLib_preoptm_type_tp *what_type;
} GMRFLib_preoptm_arg_tp;

typedef struct {

	GMRFLib_graph_tp *graph;
	GMRFLib_Qfunc_tp *Qfunc;
	void *Qfunc_arg;
	GMRFLib_constr_tp *constr;
} GMRFLib_preoptm_tp;

int GMRFLib_preopt_init(GMRFLib_preoptm_tp ** preoptm, int n, 
			double *logprec_unstruct, double **logprec_unstruct_omp,
			int nf, int **c, double **w,
			GMRFLib_graph_tp ** f_graph, GMRFLib_Qfunc_tp ** f_Qfunc,
			void **f_Qfunc_arg, char *f_sumzero, GMRFLib_constr_tp ** f_constr,
			GMRFLib_Qfunc_tp *** ff_Qfunc, void ***ff_Qfunc_arg,
			int nbeta, double **covariate, double *prior_precision, 
			GMRFLib_ai_param_tp * UNUSED(ai_par));

GMRFLib_preoptm_type_tp GMRFLib_preoptm_what_type(int node, GMRFLib_preoptm_arg_tp * a);
double GMRFLib_preoptm_Qfunc(int node, int nnode, double *values, void *arg);
int GMRFLib_free_preoptm(GMRFLib_preoptm_tp * preoptm);

__END_DECLS
#endif


