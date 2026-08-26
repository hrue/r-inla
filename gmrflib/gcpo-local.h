/* gcpo-local.h: the local gcpo group build and the store-based pair covariances (gcpo-local.c) */
#ifndef __GMRFLib_GCPO_LOCAL_H__
#define __GMRFLib_GCPO_LOCAL_H__

#include "GMRFLib/GMRFLib.h"

typedef struct {
	GMRFLib_graph_tp *lg;			       /* factor graph of the build store */
	int nlatent;
	GMRFLib_idx_tp **touch;			       /* lift: latent -> data-nodes touching it */
	char *hub;				       /* in H: reachable, never expanded (== cond) */
	char *cond;				       /* H: certificate conditioning set */
	char *lift;				       /* contributes its data points as candidates */
	char *ctouch;				       /* latents in a constraint support: a straddling constraint
						        * re-couples the separated sides, no separator equality */
	int *dist;				       /* walk marks: -1 untouched, else the BFS ring */
	int *stack;				       /* visited latents in discovery order */
	int *mark;				       /* per-node stamps for the support union (demands) */
	int *blist;				       /* the support union of the candidates */
	GMRFLib_idxval_tp **Aidx;		       /* drop: data-node -> A-support */
	GMRFLib_idx_tp **kp;			       /* demanded Qinv pairs (smaller latent -> partners) */
} GMRFLib_gcpo_local_ctx_tp;

typedef struct {
	int radius;				       /* build radius (0 = radius build off) */
	int ab;					       /* A/B mode: run both paths, compare */
	int verbose;
	int cert_ok;				       /* certificate usable (hubs within budget, SHH spd, tab present, pilot not failed) */
	int nhub;
	GMRFLib_idx_tp *hubs;			       /* H as a list */
	GMRFLib_idx_tp **cand;			       /* per-node proposal: the candidate data-nodes */
	GMRFLib_idx_tp **ballI;			       /* per-node K: the certificate interior */
	char *sep;				       /* per-node: walk exhausted its component (H seals the core) */
	char *over;				       /* per-node: candidates overflowed, solve path */
	char *certflag;				       /* A/B only: 1=certified, 2=refused (shadow kept) */
	GMRFLib_idx_tp *fb_sel;			       /* the nodes left for the solve loop */
	double *rhoH;				       /* per-data-node hub share rho_i */
	double rhoH_max;
	double *r1;				       /* |H| == 1: signed hub correlation r_i = cor(eta_i, x_h), rhoH = |r1| */
	double *zh;				       /* |H| == 1: the hub column Sigma e_h over the latents */
	double varh;				       /* |H| == 1: Var(x_h) */
	int *ord;				       /* |H| == 1: data nodes by rhoH descending, index ascending */
	int nord;
	GMRFLib_gcpo_local_ctx_tp ctx;				       /* the walk state */
} GMRFLib_gcpo_local_tp;

int GMRFLib_gcpo_form_levels(int node, int nc, int *cidx, int full, double *cor, double *cor_abs, size_t *largest, GMRFLib_gcpo_param_tp *gcpo_param,
			     GMRFLib_idxval_tp **group, double *v_last, size_t *ntrunc, int *ltrunc);
int GMRFLib_gcpo_local_prepare(GMRFLib_gcpo_local_tp *rb, GMRFLib_ai_store_tp *build_ai_store, GMRFLib_idxval_tp **Aidx,
			       GMRFLib_gcpo_param_tp *gcpo_param, GMRFLib_idx_tp *d_idx, int Npred);
void GMRFLib_gcpo_local_hubs(GMRFLib_gcpo_local_tp *rb, GMRFLib_problem_tp *pb, GMRFLib_idxval_tp **Aidx, GMRFLib_idx_tp *d_idx, double *isd, int mnpred);
GMRFLib_idx_tp *GMRFLib_gcpo_local_lookup(int thread_id, GMRFLib_gcpo_local_tp *rb, GMRFLib_ai_store_tp *build_ai_store, GMRFLib_idxval_tp **Aidx,
					  GMRFLib_gcpo_param_tp *gcpo_param, GMRFLib_idx_tp *d_idx, GMRFLib_idx_tp *selection, double *isd, double min_sd,
					  GMRFLib_idxval_tp **groups, int nt_outer, int Npred, GMRFLib_idxval_tp ***groups_rb);
void GMRFLib_gcpo_local_finish(GMRFLib_gcpo_local_tp *rb, GMRFLib_idxval_tp **groups_rb, GMRFLib_idxval_tp **groups, GMRFLib_idx_tp *selection, int Npred);
void GMRFLib_gcpo_local_demands(GMRFLib_ai_store_tp *ai_store, GMRFLib_idxval_tp **Aidx, GMRFLib_idx2_tp **missing, int Npred, int verbose);
GMRFLib_idx2_tp *GMRFLib_gcpo_cov_lookup(GMRFLib_problem_tp *pb, GMRFLib_idxval_tp **Aidx, GMRFLib_gcpo_groups_tp *groups, GMRFLib_idx_tp *node_idx,
					 GMRFLib_gcpo_elm_tp **gcpo, double *sd, double *lpred_variance, int timing);
void GMRFLib_gcpo_cov_gram(GMRFLib_problem_tp *pb, GMRFLib_idxval_tp **Aidx, GMRFLib_gcpo_groups_tp *groups, GMRFLib_idx2_tp *pair_fb,
			   GMRFLib_gcpo_elm_tp **gcpo, double *sd, double *lpred_variance, int Npred, int nt, int timing);
void GMRFLib_gcpo_cov_solve(GMRFLib_ai_store_tp *ai_store_id, GMRFLib_preopt_tp *preopt, GMRFLib_gcpo_groups_tp *groups, GMRFLib_idx_tp *node_idx,
			    GMRFLib_gcpo_elm_tp **gcpo, double *sd, double *lpred_variance, unsigned char *skip, int n_skip, int nt_outer, int nt_inner,
			    int serial, GMRFLib_gcpo_param_tp *gcpo_param, int detailed_output, int gcpo_timing);
void GMRFLib_gcpo_cov_pairs(GMRFLib_ai_store_tp *ai_store_id, GMRFLib_preopt_tp *preopt, GMRFLib_gcpo_groups_tp *groups, GMRFLib_idx_tp *node_idx,
			    GMRFLib_gcpo_elm_tp **gcpo, double *sd, double *lpred_variance, unsigned char *skip, int n_skip, int nt_outer, int nt_inner,
			    int serial, GMRFLib_gcpo_param_tp *gcpo_param, int detailed_output, int timing);

#endif
