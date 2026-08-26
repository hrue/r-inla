/* gcpo-local.c
 *
 * The local build of the gcpo (leave-group-out CV) groups and the pair
 * covariances read from the Qinv store. Called from GMRFLib_gcpo_build() and
 * GMRFLib_gcpo() in approx-inference.c.
 *
 * A group is the top level sets of |cor(eta_i, eta_j)| over all data nodes j.
 * Instead of one sparse solve per data node against every other node:
 *
 *   walk      from supp(a_i) on the latent graph ('build.radius' hops; hub
 *             latents -- intercept-like, shared by many data points -- settle
 *             but never expand), giving the interior K; the proposal is every
 *             data node touching K
 *   lookup    the candidate correlations come from the Qinv store, whose demand
 *             set (GMRFLib_qinv_keep_pairs, smtp-taucs.c) is the candidates'
 *             support products
 *   certify   for every non-candidate j, with H the hub set and
 *             sig2 = Var(eta_i | x outside K) / Var(eta_i) from one small
 *             Cholesky of Q_KK:
 *               |cor(eta_i,eta_j)| <= min( sqrt(1 - sig2),
 *                                          rho_i rho_max + sqrt(1 - sig2 - rho_i^2) )
 *             strictly below the deepest accepted level, outside the tie band,
 *             proves the group complete. a single hub sealing the core gives
 *             the far field exactly (rank one, cor = r_i r_j) and needs no bound
 *   fail fast a node that cannot be concluded goes to the solve loop; a pilot
 *             switches the certificate off when it concludes too few nodes
 *
 * Per configuration the in-group covariances are read from the store as well;
 * what it misses is computed by Gram half-solves or plain solves.
 */

#include <assert.h>
#include <float.h>
#include <math.h>
#include <omp.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include <time.h>

#include "GMRFLib/GMRFLib.h"
#include "GMRFLib/hashP.h"
#include "GMRFLib/gcpo-local.h"

// tiny dense in-place Cholesky (row-major lower factor) and forward solve, for
// the radius-build separator-certificate (matrix sizes = candidate-ball sizes)
static int rb_chol_(int n, double *A)
{
	for (int j = 0; j < n; j++) {
		double d = A[j * n + j];
		for (int k = 0; k < j; k++) {
			d -= A[j * n + k] * A[j * n + k];
		}
		if (d <= 0.0) {
			return 1;
		}
		d = sqrt(d);
		A[j * n + j] = d;
		for (int i = j + 1; i < n; i++) {
			double s = A[i * n + j];
			for (int k = 0; k < j; k++) {
				s -= A[i * n + k] * A[j * n + k];
			}
			A[i * n + j] = s / d;
		}
	}
	return 0;
}

static void rb_fsolve_(int n, double *L, double *b)
{
	for (int i = 0; i < n; i++) {
		double s = b[i];
		for (int k = 0; k < i; k++) {
			s -= L[i * n + k] * b[k];
		}
		b[i] = s / L[i * n + i];
	}
}

// canonical order for truncating a tie level-set at size_max: exact ties carry
// no information to prefer one member over another, so the fp sort order must
// not decide membership -- the smallest data indices win, making both build
// paths and repeated runs select the same subset
typedef struct {
	int idx;
	double val;
} gcpo_iv_tp_;

static int gcpo_iv_cmp_(const void *a, const void *b)
{
	return (((const gcpo_iv_tp_ *) a)->idx - ((const gcpo_iv_tp_ *) b)->idx);
}

// level-set formation, shared by the radius-lookup and the solve paths. the
// candidates are cor[0..nc-1] (signed) / cor_abs (absolute) with data-node ids
// cidx[]; largest is scratch of length nc. walking down the sorted |cor|, a
// new level opens when equal_cor says the value differs; weights accumulate
// until num_level_sets is reached or size_max is hit. the size_max cap never
// splits an equal-cor tie by fp sort order: the level that overflows the cap
// is truncated canonically (all tied members collected, the smallest data
// indices win) and deeper levels are dropped; ntrunc/ltrunc count those.
// full=1 says the nc candidates are ALL data nodes (solve path, separated
// node): an exhausted scan then simply concludes. full=0 (lookup path): an
// exhausted scan concludes only if every requested level set was found (the
// certificate then supplies the band-completion witness), else returns 0 =
// inconclusive. *v_last = the deepest accepted |cor| (the certificate target).
// returns 2 instead of 1 when the conclusion rested on 'full' with fewer level
// sets than requested (the list ended)
int GMRFLib_gcpo_form_levels(int node, int nc, int *cidx, int full, double *cor, double *cor_abs, size_t *largest,
			     GMRFLib_gcpo_param_tp *gcpo_param, GMRFLib_idxval_tp **group, double *v_last, size_t *ntrunc, int *ltrunc)
{
#define W(node_) (gcpo_param->weights[node_])
#define LEGAL_TO_ADD(node_) (!(gcpo_param->group_selection) ? 1 :	\
			     GMRFLib_iwhich_sorted(node_, gcpo_param->group_selection->idx, (unsigned int) gcpo_param->group_selection->n) >= 0)

	int nls = IABS(gcpo_param->num_level_sets);
	int levels_ok = 0, exhausted = 0;
	double levels_magnify = 1.0;
	double cor_abs_prev = 1.0;

	while (!levels_ok && !exhausted) {
		(*group)->n = 0;
		int siz_g = IMIN(nc, (int) (levels_magnify * (nls + 4L)));
		int capped = (siz_g == nc);
		levels_magnify *= 4.0;
		gsl_sort_largest_index(largest, (size_t) siz_g, cor_abs, (size_t) 1, (size_t) nc);

		double sumw = W(node);
		cor_abs_prev = 1.0;
		int i_prev = cidx[(int) largest[0]];
		GMRFLib_idxval_add(group, i_prev, cor_abs_prev);
		int lvl_start = 0, lvl = 1, cut = 0;
		for (int i = 1; i < siz_g && !levels_ok; i++) {
			int i_new_l = (int) largest[i];
			int i_new = cidx[i_new_l];
			double cor_abs_new = cor_abs[i_new_l];
			if (LEGAL_TO_ADD(i_new)) {
				if (!GMRFLib_equal_cor(cor_abs_new, cor_abs_prev, gcpo_param)) {
					// a new level: we had to look one further to collect all equal ones first
					lvl_start = (*group)->n;
					if ((sumw >= nls) || (gcpo_param->size_max > 0 && (*group)->n >= gcpo_param->size_max)) {
						levels_ok = 1;
					} else {
						lvl++;
						sumw += W(i_new);
						i_prev = i_new;
						cor_abs_prev = cor_abs_new;
						GMRFLib_idxval_add(group, i_new, cor[i_new_l]);
					}
				} else {
					cor_abs[i_new_l] = cor_abs_prev;
					cor[i_new_l] = DSIGN(cor[i_new_l]) * cor_abs_prev;
					GMRFLib_idxval_add(group, i_new, cor[i_new_l]);
					if (W(i_new) > W(i_prev)) {
						// use the maximum weight when they are equal
						sumw += W(i_new) - W(i_prev);
						i_prev = i_new;
					}
					if (gcpo_param->size_max > 0 && (*group)->n > gcpo_param->size_max) {
						cut = 1;	       /* level 'lvl' overflows the cap */
						levels_ok = 1;
					}
				}
			}
		}
		if (cut) {
			// rebuild the included part of the overflowing level from ALL tied
			// candidates in index order (one equal_cor pass, independent of the fp sort)
			int keep = gcpo_param->size_max - lvl_start;
			int nband = 0;
			gcpo_iv_tp_ *band = Malloc(nc, gcpo_iv_tp_);
			for (int c = 0; c < nc; c++) {
				int j = cidx[c];
				if (j != node && LEGAL_TO_ADD(j) && GMRFLib_equal_cor(cor_abs[c], cor_abs_prev, gcpo_param)) {
					band[nband].idx = j;
					band[nband].val = DSIGN(cor[c]) * cor_abs_prev;
					nband++;
				}
			}
			QSORT_FUN(band, (size_t) nband, sizeof(gcpo_iv_tp_), gcpo_iv_cmp_);
			(*group)->n = lvl_start;
			if (lvl_start == 0) {
				// the |cor|=1 band itself overflows: the node keeps its seat
				GMRFLib_idxval_add(group, node, 1.0);
				keep--;
			}
			for (int c = 0; c < IMIN(keep, nband); c++) {
				GMRFLib_idxval_add(group, band[c].idx, band[c].val);
			}
			Free(band);
			(*ntrunc)++;
			*ltrunc = IMAX(*ltrunc, lvl);
		}
		if (!levels_ok) {
			if (sumw > nls) {
				levels_ok = 1;
			} else if (capped) {
				// the whole candidate list was scanned
				if (full || sumw >= nls) {
					levels_ok = (full && sumw < nls ? 2 : 1);	/* 2: concluded only because the list was complete */
				} else {
					exhausted = 1;
				}
			}
		}
	}
	if (levels_ok && gcpo_param->verbose) {
		printf("%s[%1d]: for node=%1d : num.nodes %1d, deepest |cor| %g\n", __GMRFLib_FuncName, omp_get_thread_num(), node,
		       (*group)->n, cor_abs_prev);
	}
	*v_last = cor_abs_prev;
	return levels_ok;
#undef W
#undef LEGAL_TO_ADD
}

// core-first radius build: ONE walk on the factor graph is read out twice --
// the certificate interior K = visited \ cond ('ballI'), and the proposal =
// lift(visited \ hub) = the data-nodes touching the walk ('cand'). growing
// the walk grows both in lockstep, so the certificate hypothesis "any data
// node touching K is a candidate" holds by construction:
// lift(K) subset-of proposal, since K subset-of visited\cond subset-of visited\hub.
// RB_KCAP is the dense-Cholesky budget of the certificate: the interior is
// truncated to its RB_KCAP closest members (stack order), which is always
// legal (any K works, only tightness changes). RB_NHUB_MAX caps the
// conditioning set H (|H| exact column solves and an |H|x|H| Cholesky).
// RB_FP_MARGIN is the rounding-error allowance inside the certificate's square
// roots. RB_PILOT nodes are tried first: if the certificate concludes fewer
// than RB_PILOT_MIN of them, the rest skips straight to the solves (the build
// must never cost more than the solves it replaces)
#define RB_KCAP 256
#define RB_NHUB_MAX 256
#define RB_FP_MARGIN 1.0e-6
#define RB_PILOT 256
#define RB_PILOT_MIN 0.10
#define RB_COND_TOUCH 64			       /* a latent touched by more data points than this is a hub (in H) */
#define RB_LIFT_MAX 256				       /* a hub contributes candidates only up to this many data points */
#define RB_CAND_MAX 4096			       /* candidate cap per node: beyond it the node goes to the solve path */


// grow the walk of 'node' by hop-BFS from its A-support, hubs settle but do not
// expand. two-sided growth targets (0 = off): keep growing past grow_radius
// until the interior holds min_ball latents AND the lifted proposal is
// estimated to hold min_cand data nodes; the retry rounds pass the previous
// sizes doubled, so the climb is geometric and self-paced, the certificate
// arbitrates. returns 1 iff the walk exhausted its component: visited\cond is
// closed in the graph minus cond, i.e. H seals the core off ('separator
// equality': every non-candidate correlation is then EXACTLY the hub part),
// 0 otherwise, and -1 when the candidates overflowed RB_CAND_MAX (nothing kept,
// the node goes to the solve path)
static int rb_grow_node_(GMRFLib_gcpo_local_ctx_tp *c, int node, int grow_radius, int min_ball, int min_cand, GMRFLib_idx_tp **cand, GMRFLib_idx_tp **ballI)
{
	GMRFLib_graph_tp *lg = c->lg;
	GMRFLib_idxval_tp *va = c->Aidx[node];
	int closed = 1;

	// a regrow (retry round) rebuilds the readouts from scratch; stale
	// demands in c->kp are merely extra kept pairs
	if (cand[node]) {
		GMRFLib_idx_free(cand[node]);
		cand[node] = NULL;
	}
	if (ballI[node]) {
		GMRFLib_idx_free(ballI[node]);
		ballI[node] = NULL;
	}

	int ns = 0, nball = 0, ncest = 0;
	for (int ka = 0; ka < va->n; ka++) {
		int a = va->idx[ka];
		if (c->dist[a] < 0) {
			c->dist[a] = 0;
			c->stack[ns++] = a;
			nball += !c->cond[a];
			if (c->lift[a] && c->touch[a]) {
				ncest += c->touch[a]->n;
			}
		}
	}
	int lo = 0, hi = ns;
	for (int r = 1; r <= c->nlatent; r++) {
		if (r > grow_radius && nball >= min_ball && ncest >= min_cand) {
			break;
		}
		if (nball >= RB_KCAP) {
			break;				       /* interior budget reached */
		}
		int added = 0;
		for (int t = lo; t < hi; t++) {
			int a = c->stack[t];
			if (c->hub[a]) {
				continue;		       /* reachable, but does not expand */
			}
			for (int kk = 0; kk < lg->nnbs[a]; kk++) {
				int b = lg->nbs[a][kk];
				if (c->dist[b] < 0) {
					c->dist[b] = r;
					c->stack[ns++] = b;
					nball += !c->cond[b];
					if (c->lift[b] && c->touch[b]) {
						ncest += c->touch[b]->n;
					}
					added++;
				}
			}
		}
		if (!added) {
			break;				       /* component exhausted */
		}
		lo = hi;
		hi = ns;
	}
	// the two readouts. proposal: every data node touching a non-hub visited
	// latent. interior K: the non-cond visited latents (H excluded: Var(eta |
	// rest) is computed from the Q-submatrix over K), truncated at RB_KCAP in
	// stack order (closest first)
	int ncand_raw = 0, overflow = 0;
	for (int t = 0; t < ns; t++) {
		int a = c->stack[t];
		GMRFLib_idx_tp *tc = c->touch[a];
		if (tc && c->lift[a]) {
			ncand_raw += tc->n;
			if (ncand_raw > RB_CAND_MAX) {
				// too many candidates to store and demand: this node is
				// better served by one exact solve. nothing is kept
				overflow = 1;
				break;
			}
			GMRFLib_idx_nadd(&(cand[node]), tc->n, tc->idx);
		}
		if (!c->cond[a] && (!ballI[node] || ballI[node]->n < RB_KCAP)) {
			GMRFLib_idx_add(&(ballI[node]), a);
		}
	}
	// separator-equality detection: visited\cond is closed under the graph
	// minus cond iff no non-cond neighbour was left unvisited (uses the dist
	// marks, so it must run before they are cleared)
	for (int t = 0; t < ns && closed; t++) {
		int a = c->stack[t];
		if (c->cond[a]) {
			continue;
		}
		if (c->ctouch && c->ctouch[a]) {
			closed = 0;
			break;
		}
		for (int kk = 0; kk < lg->nnbs[a]; kk++) {
			int b = lg->nbs[a][kk];
			if (!c->cond[b] && c->dist[b] < 0) {
				closed = 0;
				break;
			}
		}
	}
	for (int t = 0; t < ns; t++) {
		c->dist[c->stack[t]] = -1;
	}
	if (overflow) {
		GMRFLib_idx_free(cand[node]);
		cand[node] = NULL;
		GMRFLib_idx_free(ballI[node]);
		ballI[node] = NULL;
		return -1;
	}
	if (cand[node]) {
		GMRFLib_idx_sort(cand[node]);
		GMRFLib_idx_uniq(cand[node]);
		// every candidate pair is demanded (and closure-computed if outside
		// the fill): the group forms from exact values only, no miss-as-zero.
		// pairs on the Q-pattern are always kept by the store and are not
		// demanded (hubs are neighbours of everything on the posterior graph:
		// with wide A-supports those pairs would dominate the demand list)
		// the candidates' supports overlap heavily (they share the hubs): take
		// their union once, then demand supp(a_i) x union
		int stamp = node + 1, nb = 0;
		for (int cc = 0; cc < cand[node]->n; cc++) {
			int nnode = cand[node]->idx[cc];
			if (nnode == node) {
				continue;
			}
			GMRFLib_idxval_tp *vb = c->Aidx[nnode];
			for (int kb = 0; kb < vb->n; kb++) {
				int b = vb->idx[kb];
				if (c->mark[b] != stamp) {
					c->mark[b] = stamp;
					c->blist[nb++] = b;
				}
			}
		}
		for (int ka = 0; ka < va->n; ka++) {
			int a = va->idx[ka];
			for (int kb = 0; kb < nb; kb++) {
				int b = c->blist[kb];
				if (a != b && !GMRFLib_graph_is_nb(a, b, lg)) {
					GMRFLib_idx_add(&(c->kp[IMIN(a, b)]), IMAX(a, b));
				}
			}
		}
	}
	if (ballI[node]) {
		GMRFLib_idx_sort(ballI[node]);
	}
	return closed;
}

// the separator certificate for one node (level-1 theorem of the notes): for
// every data node j outside the proposal,
//   |cor(eta_i,eta_j)| <= rho_i rho_max + sqrt(1 - sig2_i - rho_i^2),
// rho_i = hub share of eta_i (exact hub columns), sig2_i = Var(eta_i | x_F) /
// Var(eta_i) = a_K' Q_KK^{-1} a_K isd_i^2 from one small Cholesky of the Q-block
// over the interior K (Markov: conditioning on the ball complement separates).
// the level-0 cap sqrt(1 - sig2_i) is valid as well, so the min of the two is
// used. constrained problems ('constrained locality'): given x_F the active
// constraints reduce to A_K x_K, their explained part is conditioned away on
// the same local factor, sig2 <- sig2 - w'(A_K Q_KK^{-1} A_K')^{-1} w with
// w = A_K Q_KK^{-1} a_K; inactive rows are constants given x_F and drop out,
// dependent active rows refuse the node; soft constraints over-correct, hence
// conservative.
// returns 1 iff bnd < v_last strictly outside the tie band (the group is
// complete); *perm = 1 marks a refusal no regrown interior can lift
static int rb_certify_(int thread_id, GMRFLib_problem_tp *pb, GMRFLib_idx_tp *bi, GMRFLib_idxval_tp *va, double isd_i,
		       double rh, double rhoH_max, double v_last, GMRFLib_gcpo_param_tp *gcpo_param, int *perm)
{
	*perm = 0;
	// fp-margin guard: near the |cor|=1 pile-up the equal_cor band (abs
	// halfwidth ~ eps(1-v^2)/2) gets narrower than the fp agreement of two
	// exact-in-theory computation paths; group identity is then not
	// fp-well-defined in ANY implementation, so refuse and let the solve
	// path decide. v_last (the nls-th distinct level value) can only grow
	// as candidates are added, so this refusal is permanent: never retry
	if (0.5 * gcpo_param->epsilon * (1.0 - v_last * v_last) <= 1.0e-9) {
		*perm = 1;
		return 0;
	}
	if (!bi || bi->n == 0 || bi->n > RB_KCAP || !pb->tab) {
		return 0;
	}
	GMRFLib_tabulate_Qfunc_tp *tab = pb->tab;
	GMRFLib_graph_tp *lg = pb->sub_graph;
	int nb = bi->n;
	int ok = 0;
	double *QII = Calloc((size_t) nb * nb, double);
	double *aI = Calloc(nb, double);
	for (int r = 0; r < nb; r++) {
		int a = bi->idx[r];
		QII[(size_t) r * nb + r] = tab->Qfunc(thread_id, a, a, NULL, tab->Qfunc_arg);
		for (int kk = 0; kk < lg->nnbs[a]; kk++) {
			int b = lg->nbs[a][kk];
			int c2 = GMRFLib_iwhich_sorted(b, bi->idx, (unsigned int) nb);
			if (c2 >= 0) {
				QII[(size_t) r * nb + c2] = tab->Qfunc(thread_id, a, b, NULL, tab->Qfunc_arg);
			}
		}
	}
	for (int ka = 0; ka < va->n; ka++) {
		int c2 = GMRFLib_iwhich_sorted(va->idx[ka], bi->idx, (unsigned int) nb);
		if (c2 >= 0) {
			aI[c2] = va->val[ka];
		}
	}
	if (!rb_chol_(nb, QII)) {
		rb_fsolve_(nb, QII, aI);
		double vc = 0.0;
		for (int r = 0; r < nb; r++) {
			vc += aI[r] * aI[r];
		}
		GMRFLib_constr_tp *rcn = pb->sub_constr;
		int cfail = 0;
		if (rcn && rcn->nc > 0) {
			int ncr = rcn->nc, nact = 0;
			double *G = Malloc((size_t) ncr * nb + ncr + (size_t) ncr * ncr, double);
			double *w = G + (size_t) ncr * nb;
			double *M = w + ncr;
			for (int cc = 0; cc < ncr; cc++) {
				double *g = G + (size_t) nact * nb;
				int nz = 0;
				for (int r = 0; r < nb; r++) {
					g[r] = rcn->a_matrix[cc + (size_t) bi->idx[r] * ncr];
					nz += (g[r] != 0.0);
				}
				if (nz) {
					rb_fsolve_(nb, QII, g);
					nact++;
				}
			}
			for (int c1 = 0; c1 < nact; c1++) {
				double sw = 0.0;
				for (int r = 0; r < nb; r++) {
					sw += G[(size_t) c1 * nb + r] * aI[r];
				}
				w[c1] = sw;
				for (int c2 = 0; c2 <= c1; c2++) {
					double m = 0.0;
					for (int r = 0; r < nb; r++) {
						m += G[(size_t) c1 * nb + r] * G[(size_t) c2 * nb + r];
					}
					M[c1 * nact + c2] = M[c2 * nact + c1] = m;
				}
			}
			if (nact > 0) {
				if (rb_chol_(nact, M)) {
					cfail = 1;
				} else {
					rb_fsolve_(nact, M, w);
					double vcorr = 0.0;
					for (int c1 = 0; c1 < nact; c1++) {
						vcorr += w[c1] * w[c1];
					}
					vc = DMAX(0.0, vc - vcorr);
				}
			}
			Free(G);
		}
		if (!cfail) {
			// the leak variance 1 - sig2 - rho^2 is a difference of O(1) quantities
			// that come out of sparse solves on Q; when the core is (nearly)
			// sealed it is ~0 and the square root amplifies any rounding error
			// in them (an ill-conditioned Q, e.g. 'copy' with precision 1e8,
			// leaves ~1e-6 relative error in Sigma, which the root turned into a
			// 2e-4 bound violation). RB_FP_MARGIN inside the roots absorbs that
			// at a tightness cost of at most sqrt(RB_FP_MARGIN)
			double sig2 = TRUNCATE(vc * isd_i * isd_i, 0.0, 1.0);
			double bnd0 = sqrt(1.0 - sig2 + RB_FP_MARGIN);	/* level 0: condition on x_F only */
			double bnd1 = rh * rhoH_max + sqrt(DMAX(0.0, 1.0 - sig2 - rh * rh) + RB_FP_MARGIN);	/* level 1: hubs split off */
			double bnd = DMIN(bnd0, bnd1);
			ok = (bnd < v_last) && !GMRFLib_equal_cor(bnd, v_last, gcpo_param);
		}
	}
	Free(QII);
	Free(aI);
	return ok;
}

// the radius-build state, shared by the stages below and freed by rb_free_

// stage 1, the walk: lift map (touch), the two hub roles, round-0 walks for
// every data node, and the demanded Qinv pairs installed as the global
// keep-pairs so that the store computed next keeps (and closure-computes)
// exactly the pairs the lookup will ask for
static void rb_setup_(GMRFLib_gcpo_local_tp *rb, GMRFLib_problem_tp *pb, GMRFLib_idxval_tp **Aidx, GMRFLib_gcpo_param_tp *gcpo_param,
		      GMRFLib_idx_tp *d_idx, int Npred)
{
	double tref = GMRFLib_timer();
	GMRFLib_graph_tp *lg = pb->sub_graph;
	int nlatent = lg->n;
	GMRFLib_gcpo_local_ctx_tp *c = &(rb->ctx);

	c->lg = lg;
	c->nlatent = nlatent;
	c->Aidx = Aidx;
	c->touch = Calloc(nlatent, GMRFLib_idx_tp *);
	for (int k = 0; k < d_idx->n; k++) {
		int node = d_idx->idx[k];
		GMRFLib_idxval_tp *va = Aidx[node];
		for (int ka = 0; ka < va->n; ka++) {
			GMRFLib_idx_add(&(c->touch[va->idx[ka]]), node);
		}
	}
	// the conditioning set H: every latent touched by more than RB_COND_TOUCH
	// data points (intercept, fixed effects, shared effects) plus all latents of
	// any model component with dim <= 64 (the vb_nodes philosophy: few global
	// effects). the walk lives on the POSTERIOR graph, where two latents are
	// neighbours as soon as they share a data point, so a latent shared by
	// thousands of data points makes everything a 2-hop neighbour of everything:
	// the walk never expands through H (conditioning on H is exactly what removes
	// those paths). a latent of H still contributes its data points as
	// candidates when it is small enough (touch <= RB_LIFT_MAX): the members of
	// a small shared effect are then found exactly; larger shared effects are
	// left to the hub channel of the certificate and the separator equality
	c->hub = Calloc(nlatent, char);
	c->cond = Calloc(nlatent, char);
	c->lift = Calloc(nlatent, char);
	for (int a = 0; a < nlatent; a++) {
		c->cond[a] = (c->touch[a] && c->touch[a]->n > RB_COND_TOUCH);
	}
	if (gcpo_param->idx_tot > 0 && gcpo_param->idx_tag && gcpo_param->idx_start && gcpo_param->idx_n) {
		int n_offset = 0, jfirst = 0;
		if (!strcmp(gcpo_param->idx_tag[0], "APredictor")) {
			n_offset = gcpo_param->idx_n[0] + gcpo_param->idx_n[1];
			jfirst = 2;
		} else if (!strcmp(gcpo_param->idx_tag[0], "Predictor")) {
			n_offset = gcpo_param->idx_n[0];
			jfirst = 1;
		}
		for (int j = jfirst; j < gcpo_param->idx_tot; j++) {
			if (gcpo_param->idx_n[j] > 64) {
				continue;		       /* large component: stays local */
			}
			for (int i = 0; i < gcpo_param->idx_n[j]; i++) {
				int k = gcpo_param->idx_start[j] - n_offset + i;
				if (k >= 0 && k < nlatent) {
					c->cond[k] = 1;
				}
			}
		}
	}
	for (int a = 0; a < nlatent; a++) {
		c->hub[a] = c->cond[a];
		c->lift[a] = (!c->cond[a] || !c->touch[a] || c->touch[a]->n <= RB_LIFT_MAX);
		if (c->cond[a]) {
			GMRFLib_idx_add(&(rb->hubs), a);
		}
	}
	// a straddling constraint re-couples the separated sides: components
	// touched by one lose the separator equality
	if (pb->sub_constr && pb->sub_constr->nc > 0) {
		GMRFLib_constr_tp *rcn = pb->sub_constr;
		c->ctouch = Calloc(nlatent, char);
		for (int cc = 0; cc < rcn->nc; cc++) {
			int jlo = (rcn->jfirst ? rcn->jfirst[cc] : 0);
			int jhi = (rcn->jlen ? jlo + rcn->jlen[cc] : nlatent);
			for (int j = jlo; j < jhi; j++) {
				if (rcn->a_matrix[cc + (size_t) j * rcn->nc] != 0.0) {
					c->ctouch[j] = 1;
				}
			}
		}
	}
	c->kp = Calloc(nlatent, GMRFLib_idx_tp *);
	c->dist = Malloc(nlatent, int);
	c->stack = Malloc(nlatent, int);
	c->mark = Calloc(nlatent, int);
	c->blist = Malloc(nlatent, int);
	for (int i = 0; i < nlatent; i++) {
		c->dist[i] = -1;
	}
	// one walk per node: reach 'radius' at least, and keep going until the
	// proposal can hold the requested level sets (no retries later)
	int min_cand = 2 * IABS(gcpo_param->num_level_sets) + 2;
	rb->cand = Calloc(Npred, GMRFLib_idx_tp *);
	rb->ballI = Calloc(Npred, GMRFLib_idx_tp *);
	rb->sep = Calloc(Npred, char);
	rb->over = Calloc(Npred, char);
	size_t n_over = 0;
	for (int k = 0; k < d_idx->n; k++) {
		int node = d_idx->idx[k];
		int rc = rb_grow_node_(c, node, rb->radius, 0, min_cand, rb->cand, rb->ballI);
		rb->sep[node] = (char) (rc > 0);
		rb->over[node] = (char) (rc < 0);
		n_over += (rc < 0);
	}
	size_t kp_n = 0;
	for (int i = 0; i < nlatent; i++) {
		if (c->kp[i]) {
			GMRFLib_idx_sort(c->kp[i]);
			GMRFLib_idx_uniq(c->kp[i]);
			kp_n += (size_t) c->kp[i]->n;
		}
	}
	// install as the global keep-pairs (owned by the global from here on:
	// rb_free_ does not free c->kp) and force the store to be recomputed
	if (GMRFLib_qinv_keep_pairs) {
		for (int i = 0; i < GMRFLib_qinv_keep_pairs_n; i++) {
			GMRFLib_idx_free(GMRFLib_qinv_keep_pairs[i]);
		}
		Free(GMRFLib_qinv_keep_pairs);
	}
	GMRFLib_qinv_keep_pairs = c->kp;
	GMRFLib_qinv_keep_pairs_n = nlatent;
	GMRFLib_free_Qinv(pb);
	if (rb->verbose) {
		printf("[gcpo-timing] build: radius-%1d candidates+demands %.4f s (%zu latent pairs beyond the Q-graph, %1d hubs, %zu nodes overflowed)\n",
		       rb->radius, GMRFLib_timer() - tref, kp_n, (rb->hubs ? rb->hubs->n : 0), n_over);
	}
}

// stage 2, the hub channel: one batched solve for the |H| hub columns Sigma
// e_h (constraint-corrected by GMRFLib_Qsolves), the Cholesky of Sigma_HH, and
// per data node the hub share rho_i = |L^{-1} c_H(eta_i)| isd_i with c_H =
// (Sigma a_i)_H; rho_max over all data nodes is the j-free cap of the hub
// channel. with a single hub the far field of a sealed core is rank one,
// cor(eta_i, eta_j) = r_i r_j: the signed r_i, the hub column and the data
// nodes in hub order are kept for that (rank-1) path
static double *rb_ord_rho_ = NULL;
static int rb_ord_cmp_(const void *a, const void *b)
{
	int i = *((const int *) a), j = *((const int *) b);
	if (rb_ord_rho_[i] > rb_ord_rho_[j]) {
		return -1;
	}
	if (rb_ord_rho_[i] < rb_ord_rho_[j]) {
		return 1;
	}
	return (i - j);
}

static void rb_hub_columns_(GMRFLib_gcpo_local_tp *rb, GMRFLib_problem_tp *pb, GMRFLib_idxval_tp **Aidx, GMRFLib_idx_tp *d_idx, double *isd, int mnpred)
{
	double tref = GMRFLib_timer();
	int nlatent = pb->sub_graph->n;
	int nhub = (rb->hubs ? rb->hubs->n : 0);

	rb->cert_ok = (pb->tab != NULL);
	if (nhub > RB_NHUB_MAX) {
		// conditioning set too large for the exact-column budget: refuse
		// the certificate (everything falls back, nothing wrong)
		rb->cert_ok = 0;
		nhub = 0;
	}
	rb->nhub = nhub;
	if (nhub == 0) {
		return;
	}
	double *Z = Calloc((size_t) nlatent * nhub, double);
	for (int h = 0; h < nhub; h++) {
		Z[(size_t) h * nlatent + rb->hubs->idx[h]] = 1.0;
	}
	GMRFLib_stiles_idx_tp sidx = { 0, -1, nhub };
	GMRFLib_Qsolves(Z, nhub, pb, &sidx);
	double *SHH = Calloc(nhub * nhub, double);
	for (int h = 0; h < nhub; h++) {
		for (int h2 = 0; h2 < nhub; h2++) {
			SHH[h * nhub + h2] = Z[(size_t) h2 * nlatent + rb->hubs->idx[h]];
		}
	}
	if (rb_chol_(nhub, SHH)) {
		rb->cert_ok = 0;			       /* refuse: all nodes will fall back */
	} else {
		double *ub = Calloc(nhub, double);
		rb->rhoH = Calloc(mnpred, double);
		if (nhub == 1) {
			rb->r1 = Calloc(mnpred, double);
			rb->zh = Malloc(nlatent, double);
			Memcpy(rb->zh, Z, nlatent * sizeof(double));
			rb->varh = SHH[0] * SHH[0];
		}
		for (int k = 0; k < d_idx->n; k++) {
			int node = d_idx->idx[k];
			GMRFLib_idxval_tp *va = Aidx[node];
			for (int h = 0; h < nhub; h++) {
				double sum = 0.0;
				for (int ka = 0; ka < va->n; ka++) {
					sum += va->val[ka] * Z[(size_t) h * nlatent + va->idx[ka]];
				}
				ub[h] = sum;
			}
			rb_fsolve_(nhub, SHH, ub);
			double s2 = 0.0;
			for (int h = 0; h < nhub; h++) {
				s2 += ub[h] * ub[h];
			}
			rb->rhoH[node] = sqrt(TRUNCATE(s2 * isd[node] * isd[node], 0.0, 1.0));
			rb->rhoH_max = DMAX(rb->rhoH_max, rb->rhoH[node]);
			if (nhub == 1) {
				rb->r1[node] = TRUNCATE(ub[0] * isd[node], -1.0, 1.0);
			}
		}
		Free(ub);
		if (nhub == 1) {
			rb->nord = d_idx->n;
			rb->ord = Malloc(rb->nord, int);
			Memcpy(rb->ord, d_idx->idx, rb->nord * sizeof(int));
			rb_ord_rho_ = rb->rhoH;
			qsort(rb->ord, (size_t) rb->nord, sizeof(int), rb_ord_cmp_);
			rb_ord_rho_ = NULL;
		}
	}
	Free(SHH);
	Free(Z);
	if (rb->verbose) {
		printf("[gcpo-timing] build: hub columns %.4f s (%1d hubs, max hub-share %.4f, cert %s%s)\n",
		       GMRFLib_timer() - tref, nhub, rb->rhoH_max, (rb->cert_ok ? "ok" : "REFUSED"), (nhub == 1 ? ", rank-1 far field" : ""));
	}
}

// stage 3, one node: the candidate correlations are read from the Qinv store,
// the level sets are formed, and the node is concluded in one of two ways:
// rank-1 (|H| = 1 and H seals the core): the far field is EXACTLY r_i r_j, so
// the candidates are extended by the data nodes in hub order -- enough to hold
// the level sets and the truncation -- and the group is complete by equality;
// otherwise the separator certificate bounds the far field. anything else
// (candidates overflowed, pairs missing, levels inconclusive, bound too weak)
// goes to the solve path. returns 1 certified, 2 certified by the rank-1
// equality, 0 solve path
static int rb_lookup_node_(int thread_id, GMRFLib_gcpo_local_tp *rb, GMRFLib_problem_tp *pb, GMRFLib_idxval_tp **Aidx, GMRFLib_gcpo_param_tp *gcpo_param,
			   GMRFLib_idx_tp *d_idx, int node, double *isd, double min_sd, GMRFLib_idxval_tp **groups, size_t *ntrunc, int *ltrunc,
			   size_t *nmiss)
{
	GMRFLib_idx_tp *cd = rb->cand[node];
	int ncand = (cd ? cd->n : 0);
	int nls = IABS(gcpo_param->num_level_sets);
	int size_max = gcpo_param->size_max;
	int r1path = (rb->nhub == 1 && rb->sep[node] && rb->ord && rb->r1 && rb->zh);

	if (rb->over[node] || ncand == 0 || (!rb->cert_ok && !r1path)) {
		return 0;
	}
	// the far-field prefix taken in hub order: enough for the level sets and
	// for the canonical truncation of the last tie run
	int take = (r1path ? 4 * (nls + 4) + IMAX(0, size_max) + 1 : 0);
	int cap = ncand + (r1path ? take + (size_max > 0 ? size_max + 1 : RB_CAND_MAX) : 0);
	int *cidx = Malloc(cap, int);
	double *cor = Malloc(2 * cap, double);
	double *cor_abs = cor + cap;
	size_t *largest = Malloc(cap, size_t);
	GMRFLib_idxval_tp *va = Aidx[node];
	double zs_eps = 1.0E-3 * min_sd / isd[node];
	int nc = 0;
	size_t node_miss = 0;

	for (int c = 0; c < ncand; c++) {
		int nnode = cd->idx[c];
		cidx[nc] = nnode;
		if (nnode == node) {
			cor[nc] = cor_abs[nc] = 1.0;
			nc++;
			continue;
		}
		GMRFLib_idxval_tp *vb = Aidx[nnode];
		double sum = 0.0;
		int hit = 1;
		for (int kb = 0; kb < vb->n && hit; kb++) {
			// assemble the latent covariance entry (Sigma a_node)_b and pass
			// it through the same zero_small gate as the solve path, so
			// band-edge ties resolve identically
			double sv = 0.0;
			for (int ka = 0; ka < va->n; ka++) {
				double *q = GMRFLib_Qinv_get(pb, va->idx[ka], vb->idx[kb]);
				if (!q) {
					hit = 0;
					break;
				}
				sv += va->val[ka] * (*q);
			}
			if (ABS(sv) > zs_eps) {
				sum += vb->val[kb] * sv;
			}
		}
		if (hit) {
			sum *= isd[node] * isd[nnode];
			cor[nc] = TRUNCATE(sum, -1.0, 1.0);
			cor_abs[nc] = ABS(cor[nc]);
		} else {
			cor[nc] = cor_abs[nc] = 0.0;	       /* should not happen: every candidate pair is demanded */
			node_miss++;
		}
		nc++;
	}
	*nmiss += node_miss;
	int full = 0, fail = 0, t_end = -1;
	double ci = 0.0, rho_i = 0.0;
	size_t dummy_nt = 0;
	int dummy_lt = 0;

	// value of a far-field data node j by equality: (Sigma a_i)_k = c_i
	// Sigma_{h,k} / Var(x_h) for every latent k outside the core, same
	// zero_small gate as the solve path
#define R1_ADD(j_)							\
	{								\
		GMRFLib_idxval_tp *vb_ = Aidx[j_];				\
		double sum_ = 0.0;						\
		for (int kb_ = 0; kb_ < vb_->n; kb_++) {			\
			double sv_ = ci * rb->zh[vb_->idx[kb_]];		\
			if (ABS(sv_) > zs_eps) {				\
				sum_ += vb_->val[kb_] * sv_;			\
			}							\
		}							\
		sum_ *= isd[node] * isd[j_];				\
		cidx[nc] = j_;						\
		cor[nc] = TRUNCATE(sum_, -1.0, 1.0);			\
		cor_abs[nc] = ABS(cor[nc]);				\
		nc++;							\
	}
#define R1_SKIP(j_) ((j_) == node || GMRFLib_iwhich_sorted(j_, cd->idx, (unsigned int) cd->n) >= 0)

	if (r1path && node_miss == 0) {
		// the far field in hub order: the prefix holds the largest far
		// correlations, enough to find the level sets
		ci = rb->r1[node] / (isd[node] * sqrt(rb->varh));
		rho_i = rb->rhoH[node];
		int taken = 0;
		for (int t = 0; t < rb->nord && taken < take; t++) {
			int j = rb->ord[t];
			if (R1_SKIP(j)) {
				continue;
			}
			R1_ADD(j);
			taken++;
			t_end = t;
		}
		full = 1;
	} else if (nc >= d_idx->n && node_miss == 0) {
		full = 1;
	}
	int ok = 0;
	double v_last = 1.0;
	int levels_ok = 0;
	if (!fail) {
		levels_ok = GMRFLib_gcpo_form_levels(node, nc, cidx, full, cor, cor_abs, largest, gcpo_param, &(groups[node]), &v_last,
					      (r1path ? &dummy_nt : ntrunc), (r1path ? &dummy_lt : ltrunc));
	}
	if (!fail && r1path && node_miss == 0 && levels_ok) {
		// the band of the deepest level among the far field: all j with
		// equal_cor(rho_i rho_j, v_last), a contiguous range of 'ord'. it must
		// be complete in cidx (the canonical truncation keeps its smallest
		// indices, the level-set logic needs all of it), or the node is not
		// concluded here
		int t_in = -1;
		for (int t = t_end; t >= 0 && t_in < 0; t--) {
			int j = rb->ord[t];
			if (!R1_SKIP(j) && GMRFLib_equal_cor(rho_i * rb->rhoH[j], v_last, gcpo_param)) {
				t_in = t;
			}
		}
		if (t_in >= 0) {
			// binary search for the ends of the band around t_in (the predicate
			// is monotone in the hub order on either side of t_in)
			int lo = t_in, hi = t_in + 1, a, b;
			a = 0;
			b = t_in;
			while (a < b) {
				int m = (a + b) / 2;
				if (GMRFLib_equal_cor(rho_i * rb->rhoH[rb->ord[m]], v_last, gcpo_param)) {
					b = m;
				} else {
					a = m + 1;
				}
			}
			lo = a;
			a = t_in + 1;
			b = rb->nord;
			while (a < b) {
				int m = (a + b) / 2;
				if (GMRFLib_equal_cor(rho_i * rb->rhoH[rb->ord[m]], v_last, gcpo_param)) {
					a = m + 1;
				} else {
					b = m;
				}
			}
			hi = a;
			if (hi - 1 > t_end) {
				// the band continues past the prefix
				int nband = hi - lo;
				int gn = groups[node]->n;
				int capped = (size_max > 0 && gn >= size_max);
				if (!capped && nband + gn <= (size_max > 0 ? size_max : RB_CAND_MAX) && nc + nband <= cap) {
					// the whole band fits in the group: include it entirely
					for (int t = t_end + 1; t < hi; t++) {
						int j = rb->ord[t];
						if (!R1_SKIP(j)) {
							R1_ADD(j);
						}
					}
				} else if (size_max > 0 && nc + size_max + 1 <= cap) {
					// the cap cuts this band: the canonical truncation keeps its
					// smallest data indices -- collect the size_max+1 smallest, by
					// enumerating the band when it is small, by scanning the data
					// indices upward when it is large (bounded scan)
					double rho_hi = rb->rhoH[rb->ord[lo]], rho_lo = rb->rhoH[rb->ord[hi - 1]];
					int got = 0;
					if (nband <= 4 * size_max) {
						for (int t = lo; t < hi && got <= size_max; t++) {
							int j = rb->ord[t];
							if (!R1_SKIP(j) && (t > t_end)) {
								R1_ADD(j);
								got++;
							} else if (!R1_SKIP(j)) {
								got++;	       /* already in the prefix */
							}
						}
					} else {
						int scanned = 0;
						for (int k = 0; k < d_idx->n && got <= size_max && scanned < RB_CAND_MAX; k++, scanned++) {
							int j = d_idx->idx[k];
							if (R1_SKIP(j)) {
								continue;
							}
							double rj = rb->rhoH[j];
							if (rj >= rho_lo && rj <= rho_hi) {
								// in the band; in the prefix already iff its ord position <= t_end
								int t = -1;
								for (int tt = t_end; tt >= 0 && tt >= t_end - take; tt--) {
									if (rb->ord[tt] == j) {
										t = tt;
										break;
									}
								}
								if (t < 0) {
									R1_ADD(j);
								}
								got++;
							}
						}
						if (got <= size_max) {
							fail = 1;	       /* scan budget exhausted: solve path */
						}
					}
				} else {
					fail = 1;
				}
				if (!fail) {
					levels_ok = GMRFLib_gcpo_form_levels(node, nc, cidx, 1, cor, cor_abs, largest, gcpo_param, &(groups[node]), &v_last, &dummy_nt,
								      &dummy_lt);
				}
			}
		}
		if (!fail && levels_ok) {
			// count the truncation once, on the final formation
			size_t nt2 = 0;
			int lt2 = 0;
			levels_ok = GMRFLib_gcpo_form_levels(node, nc, cidx, 1, cor, cor_abs, largest, gcpo_param, &(groups[node]), &v_last, &nt2, &lt2);
			*ntrunc += nt2;
			*ltrunc = IMAX(*ltrunc, lt2);
		}
	}
#undef R1_ADD
#undef R1_SKIP
	if (!fail && levels_ok) {
		if (levels_ok == 2 && r1path && nc < d_idx->n) {
			levels_ok = 0;			       /* the prefix ended before the level sets did: solve path */
		}
	}
	if (!fail && levels_ok) {
		GMRFLib_idxval_nsort_x(&(groups[node]), 1, 1, 0, 0);
		if (GMRFLib_iwhich_sorted(node, groups[node]->idx, (unsigned int) groups[node]->n) < 0) {
			GMRFLib_idxval_add(&(groups[node]), node, 1.0);
			GMRFLib_idxval_nsort_x(&(groups[node]), 1, 1, 0, 0);
		}
		if (full) {
			ok = (r1path ? 2 : 1);		       /* complete information: nothing outside to certify against */
		} else if (rb->cert_ok) {
			int perm = 0;
			ok = rb_certify_(thread_id, pb, rb->ballI[node], va, isd[node], (rb->rhoH ? rb->rhoH[node] : 0.0), rb->rhoH_max, v_last, gcpo_param, &perm);
		}
		if (rb->certflag) {
			rb->certflag[node] = (ok ? 1 : 2);	       /* 2: refused, the shadow group is kept for scoring */
		}
	}
	if (!ok && groups[node] && !(rb->certflag && rb->certflag[node] == 2)) {
		groups[node]->n = 0;
	}
	Free(cidx);
	Free(cor);
	Free(largest);
	return ok;
}

// stage 3, all nodes: a pilot over RB_PILOT evenly spaced nodes first -- if the
// certificate concludes fewer than RB_PILOT_MIN of those it is switched off
// and the rest goes straight to the solves (rank-1 nodes are unaffected: they
// need no bound). returns the nodes left for the solve path (caller frees)
static GMRFLib_idx_tp *rb_lookup_(int thread_id, GMRFLib_gcpo_local_tp *rb, GMRFLib_problem_tp *pb, GMRFLib_idxval_tp **Aidx, GMRFLib_gcpo_param_tp *gcpo_param,
				  GMRFLib_idx_tp *d_idx, GMRFLib_idx_tp *selection, double *isd, double min_sd, GMRFLib_idxval_tp **groups, int nt_outer)
{
	double tref = GMRFLib_timer();
	int nsel = selection->n;
	char *done = Calloc(nsel, char);
	GMRFLib_idx_tp **fbl = Calloc(nt_outer, GMRFLib_idx_tp *);
	size_t n_ok = 0, n_r1 = 0, n_miss = 0, ntrunc = 0;
	int ltrunc = 0;

	if (rb->cert_ok && nsel > 4 * RB_PILOT) {
		int step = nsel / RB_PILOT;
		size_t p_ok = 0, p_r1 = 0, p_tried = 0;
#pragma omp parallel for num_threads(nt_outer) schedule(dynamic, 8) reduction(+: p_ok, p_r1, p_tried, n_miss, ntrunc) reduction(max: ltrunc)
		for (int k = 0; k < RB_PILOT; k++) {
			int is = k * step;
			int node = selection->idx[is];
			int tnum = omp_get_thread_num();
			int rc = rb_lookup_node_(thread_id, rb, pb, Aidx, gcpo_param, d_idx, node, isd, min_sd, groups, &ntrunc, &ltrunc, &n_miss);
			done[is] = 1;
			if (rc == 2) {
				p_r1++;
			} else if (rc == 1) {
				p_ok++;
				p_tried++;
			} else {
				p_tried += (!rb->over[node] && rb->cand[node] && rb->cand[node]->n > 0);
				GMRFLib_idx_add(&(fbl[tnum]), node);
			}
		}
		n_ok += p_ok + p_r1;
		n_r1 += p_r1;
		if (p_tried >= 32 && (double) p_ok < RB_PILOT_MIN * (double) p_tried) {
			rb->cert_ok = 0;
			if (rb->verbose) {
				printf("[gcpo-timing] build: pilot certified %zu of %zu tried (%zu rank-1): certificate off, the rest solves\n", p_ok, p_tried, p_r1);
			}
		}
	}
#pragma omp parallel for num_threads(nt_outer) schedule(dynamic, 64) reduction(+: n_ok, n_r1, n_miss, ntrunc) reduction(max: ltrunc)
	for (int is = 0; is < nsel; is++) {
		if (done[is]) {
			continue;
		}
		int node = selection->idx[is];
		int tnum = omp_get_thread_num();
		int rc = rb_lookup_node_(thread_id, rb, pb, Aidx, gcpo_param, d_idx, node, isd, min_sd, groups, &ntrunc, &ltrunc, &n_miss);
		if (rc) {
			n_ok++;
			n_r1 += (rc == 2);
		} else {
			GMRFLib_idx_add(&(fbl[tnum]), node);
		}
	}
	GMRFLib_idx_tp *fb_sel = NULL;
	for (int t = 0; t < nt_outer; t++) {
		if (fbl[t]) {
			GMRFLib_idx_nadd(&fb_sel, fbl[t]->n, fbl[t]->idx);
			GMRFLib_idx_free(fbl[t]);
		}
	}
	Free(fbl);
	Free(done);
	if (!fb_sel) {
		GMRFLib_idx_create_x(&fb_sel, 1);	       /* empty */
	}
	if (rb->verbose) {
		printf("[gcpo-timing] build: radius-lookup groups %.4f s (%zu of %1d certified (%zu rank-1), %1d to solve-fallback; %zu fill-miss pairs)\n",
		       GMRFLib_timer() - tref, n_ok, nsel, n_r1, fb_sel->n, n_miss);
	}
	if (ntrunc) {
		printf("[gcpo] WARNING: size.max=%1d truncated a tie level-set at %zu of %1d nodes (deepest at level %1d); "
		       "tied members kept by smallest index, deeper levels dropped\n", gcpo_param->size_max, ntrunc, nsel, ltrunc);
	}
	return fb_sel;
}

// A/B: compare the radius groups (groups_rb) with the solve-path groups; refused
// nodes kept their shadow group so refusals can be scored as true/false positives
static void rb_ab_compare_(GMRFLib_gcpo_local_tp *rb, GMRFLib_idxval_tp **groups_rb, GMRFLib_idxval_tp **groups, GMRFLib_idx_tp *selection)
{
	size_t n_cmp = 0, n_mismatch = 0, n_ref_cmp = 0, n_ref_wrong = 0;
	for (int is = 0; is < selection->n; is++) {
		int node = selection->idx[is];
		GMRFLib_idxval_tp *g_rb = groups_rb[node];
		if (!(g_rb && g_rb->n > 0)) {
			continue;			       /* radius path fell back with no shadow: nothing to compare */
		}
		GMRFLib_idxval_tp *g_ref = groups[node];
		int eq = (g_rb->n == g_ref->n);
		for (int i = 0; eq && i < g_rb->n; i++) {
			eq = (g_rb->idx[i] == g_ref->idx[i]);
		}
		if (rb->certflag && rb->certflag[node] == 2) {
			n_ref_cmp++;
			n_ref_wrong += !eq;
			continue;
		}
		n_cmp++;
		if (!eq) {
			n_mismatch++;
			if (n_mismatch <= 10) {
				printf("[gcpo-timing] build: A/B MISMATCH node %1d: radius(n=%1d) vs solve(n=%1d)\n", node, g_rb->n, g_ref->n);
			}
		}
	}
	printf("[gcpo-timing] build: A/B compare %zu certified nodes: %zu group-mismatches\n", n_cmp, n_mismatch);
	if (n_ref_cmp) {
		printf("[gcpo-timing] build: A/B refused-shadow: %zu of %zu refusals would have given a WRONG group\n", n_ref_wrong, n_ref_cmp);
	}
}

static void rb_free_(GMRFLib_gcpo_local_tp *rb, int Npred)
{
	GMRFLib_gcpo_local_ctx_tp *c = &(rb->ctx);
	for (int i = 0; i < Npred; i++) {
		GMRFLib_idx_free(rb->cand[i]);
		GMRFLib_idx_free(rb->ballI[i]);
	}
	Free(rb->cand);
	Free(rb->ballI);
	GMRFLib_idx_free(rb->hubs);
	Free(rb->rhoH);
	Free(rb->sep);
	Free(rb->r1);
	Free(rb->zh);
	Free(rb->ord);
	Free(rb->certflag);
	// c->kp is NOT freed: it is owned by the global GMRFLib_qinv_keep_pairs
	if (c->touch) {
		for (int i = 0; i < c->nlatent; i++) {
			GMRFLib_idx_free(c->touch[i]);
		}
		Free(c->touch);
	}
	Free(c->hub);
	Free(c->cond);
	Free(c->lift);
	Free(rb->over);
	Free(c->ctouch);
	Free(c->dist);
	Free(c->stack);
	Free(c->mark);
	Free(c->blist);
}

// the pairs the per-configuration gcpo lookup will ask for: the support products
// of every in-group pair (node, nnode). installed as the global Qinv demand set,
// the store then keeps exactly those fill entries (memory bounded by the group
// structure, not by the factorization fill) and closure-computes the rest
#define GCPO_LOOKUP_SUPP_MAX 1024L		       /* wider support products go to the solves */

static GMRFLib_idx_tp **gcpo_demand_set_(GMRFLib_idxval_tp **Aidx, GMRFLib_idx2_tp **missing, int Npred, GMRFLib_graph_tp *lg, int verbose, size_t *kp_n_out)
{
	int nlatent = lg->n;
	double tref = GMRFLib_timer();
	GMRFLib_idx_tp **kp = Calloc(nlatent, GMRFLib_idx_tp *);
	for (int node = 0; node < Npred; node++) {
		if (missing[node]->n == 0) {
			continue;
		}
		GMRFLib_idxval_tp *va = Aidx[node];
		for (int k = 0; k < missing[node]->n; k++) {
			int nnode = missing[node]->idx[0][k];
			if (nnode == node) {
				continue;
			}
			GMRFLib_idxval_tp *vb = Aidx[nnode];
			if ((long) va->n * (long) vb->n > GCPO_LOOKUP_SUPP_MAX) {
				continue;
			}
			for (int ka = 0; ka < va->n; ka++) {
				int a = va->idx[ka];
				for (int kb = 0; kb < vb->n; kb++) {
					int b = vb->idx[kb];
					if (a != b && !GMRFLib_graph_is_nb(a, b, lg)) {	/* Q-pattern pairs are kept anyway */
						GMRFLib_idx_add(&kp[IMIN(a, b)], IMAX(a, b));
					}
				}
			}
		}
	}
	size_t kp_n = 0;
	for (int i = 0; i < nlatent; i++) {
		if (kp[i]) {
			GMRFLib_idx_sort(kp[i]);
			GMRFLib_idx_uniq(kp[i]);
			kp_n += (size_t) kp[i]->n;
		}
	}
	if (verbose) {
		printf("[gcpo-timing] build: demand-set %.4f s (%zu latent pairs beyond the Q-graph)\n", GMRFLib_timer() - tref, kp_n);
	}
	*kp_n_out = kp_n;
	return kp;
}

// install 'kp' as the global Qinv demand set, freeing the previous one
static void gcpo_install_demand_(GMRFLib_idx_tp **kp, int nlatent)
{
	if (GMRFLib_qinv_keep_pairs && GMRFLib_qinv_keep_pairs != kp) {
		for (int i = 0; i < GMRFLib_qinv_keep_pairs_n; i++) {
			GMRFLib_idx_free(GMRFLib_qinv_keep_pairs[i]);
		}
		Free(GMRFLib_qinv_keep_pairs);
	}
	GMRFLib_qinv_keep_pairs = kp;
	GMRFLib_qinv_keep_pairs_n = nlatent;
}


// one in-group pair into the cov-matrix of group cm_idx: both diagonals from the
// predictor variances, the off-diagonal symmetric and clipped to |cor| <= 1.
// cov is ignored for node == nnode
static void gcpo_cov_set_(GMRFLib_gcpo_elm_tp **gcpo, int cm_idx, int node, int nnode, double cov, double *sd, double *lpred_variance)
{
	gsl_matrix *mat = gcpo[cm_idx]->cov_mat;
	int ii = GMRFLib_iwhich_sorted(node, (int *) gcpo[cm_idx]->idxs->idx, (unsigned int) gcpo[cm_idx]->idxs->n);
	int jj = GMRFLib_iwhich_sorted(nnode, (int *) gcpo[cm_idx]->idxs->idx, (unsigned int) gcpo[cm_idx]->idxs->n);
	assert(ii >= 0 && jj >= 0);
	gsl_matrix_set(mat, ii, ii, lpred_variance[node]);
	if (jj != ii) {
		double f = sd[node] * sd[nnode];
		double c = TRUNCATE(cov / f, -1.0, 1.0) * f;
		gsl_matrix_set(mat, jj, jj, lpred_variance[nnode]);
		gsl_matrix_set(mat, ii, jj, c);
		gsl_matrix_set(mat, jj, ii, c);
	}
}

// lookup path: the partial inverse on the Q-pattern (plus the demanded pairs) is
// in the store and constraint-corrected, so cov(eta_i, eta_j) = sum_kl a_ik a_jl
// Qinv[k,l] is a plain lookup whenever every (k,l) is present. returns the
// missed pairs as (node, k-of-missing[node]), NULL if none
GMRFLib_idx2_tp *GMRFLib_gcpo_cov_lookup(GMRFLib_problem_tp *pb, GMRFLib_idxval_tp **Aidx, GMRFLib_gcpo_groups_tp *groups,
					 GMRFLib_idx_tp *node_idx, GMRFLib_gcpo_elm_tp **gcpo, double *sd, double *lpred_variance, int timing)
{
	double tref = GMRFLib_timer();
	GMRFLib_idx2_tp *pair_fb = NULL;
	size_t nhit = 0, nmiss = 0;

	for (int i = 0; i < node_idx->n; i++) {
		int node = node_idx->idx[i];
		GMRFLib_idxval_tp *va = Aidx[node];
		for (int k = 0; k < groups->missing[node]->n; k++) {
			int nnode = groups->missing[node]->idx[0][k];
			int cm_idx = groups->missing[node]->idx[1][k];
			if (nnode == node) {
				gcpo_cov_set_(gcpo, cm_idx, node, nnode, 0.0, sd, lpred_variance);
				continue;
			}
			GMRFLib_idxval_tp *vb = Aidx[nnode];
			int ok = ((long) va->n * (long) vb->n <= GCPO_LOOKUP_SUPP_MAX);
			double sum = 0.0;
			for (int ka = 0; ka < va->n && ok; ka++) {
				for (int kb = 0; kb < vb->n; kb++) {
					double *q = GMRFLib_Qinv_get(pb, va->idx[ka], vb->idx[kb]);
					if (!q) {
						ok = 0;
						break;
					}
					sum += va->val[ka] * vb->val[kb] * (*q);
				}
			}
			if (ok) {
				nhit++;
				gcpo_cov_set_(gcpo, cm_idx, node, nnode, sum, sd, lpred_variance);
			} else {
				nmiss++;
				GMRFLib_idx2_add(&pair_fb, node, k);
			}
		}
	}
	if (timing) {
		printf("[gcpo-timing] gcpo: LOOKUP %.4f s (%zu pair-hits, %zu pair-misses to fallback)\n", GMRFLib_timer() - tref, nhit, nmiss);
	}
	return pair_fb;
}

// Gram path for the pairs the lookup missed: with Q = LL^T, cov(eta_i, eta_j) =
// (L^-1 a_i) . (L^-1 a_j), so forward solves are enough. w = L^-1 a for every
// node of a missed pair (blocked, in the mapped ordering where the dots are
// invariant, stored sparse), then one sparse dot per pair. constraints: with
// x_c = Q^-1 a - constr_m (C Q^-1 a) the corrected covariance is
// w_i.w_j - (constr_m' a_j).(C Q^-1 a_i), and C Q^-1 a_i = (L^-1 C')' w_i --
// the nc constraint rows are forward-solved once, the correction is a dot of
// two nc-vectors per pair
void GMRFLib_gcpo_cov_gram(GMRFLib_problem_tp *pb, GMRFLib_idxval_tp **Aidx, GMRFLib_gcpo_groups_tp *groups, GMRFLib_idx2_tp *pair_fb,
			   GMRFLib_gcpo_elm_tp **gcpo, double *sd, double *lpred_variance, int Npred, int nt, int timing)
{
	double tref = GMRFLib_timer();
	int nn = pb->sub_graph->n;

	unsigned char *needw = Calloc(Npred, unsigned char);
	for (int i = 0; i < pair_fb->n; i++) {
		int node = pair_fb->idx[0][i];
		int k = pair_fb->idx[1][i];
		needw[node] = 1;
		needw[groups->missing[node]->idx[0][k]] = 1;
	}
	GMRFLib_idx_tp *wlist = NULL;
	for (int i = 0; i < Npred; i++) {
		if (needw[i]) {
			GMRFLib_idx_add(&wlist, i);
		}
	}

	int nc = (pb->sub_constr && pb->sub_constr->nc > 0 ? pb->sub_constr->nc : 0);
	double *Vg = NULL, *zc = NULL, *mc = NULL;
	if (nc > 0) {
		assert(pb->constr_m);
		Vg = Malloc((size_t) nn * nc, double);
		double *tmp = Malloc(nn, double);
		double *vwork = Malloc((size_t) nn * nc, double);
		for (int cc = 0; cc < nc; cc++) {
			for (int j = 0; j < nn; j++) {
				tmp[j] = pb->sub_constr->a_matrix[cc + (size_t) j * nc];
			}
			GMRFLib_convert_to_mapped(Vg + (size_t) cc * nn, tmp, pb->sub_graph, pb->sub_sm_fact.remap);
		}
		GMRFLib_taucs_Lsolve_blocked(pb->sub_sm_fact.TAUCS_L, Vg, nc, vwork);
		Free(tmp);
		Free(vwork);
		zc = Calloc((size_t) Npred * nc, double);
		mc = Calloc((size_t) Npred * nc, double);
	}

	int **widx = Calloc(Npred, int *);
	double **wval = Calloc(Npred, double *);
	int *wnnz = Calloc(Npred, int);
	int gnt = IMAX(1, nt);
	int GB = GMRFLib_taucs_get_block_size();
	int nblocks = (wlist->n + GB - 1) / GB;
	double **gbuf = Calloc(gnt, double *);
	int **gis = Calloc(gnt, int *);
	for (int i = 0; i < gnt; i++) {
		gbuf[i] = Malloc(2 * (size_t) GB * nn, double);
		gis[i] = Malloc(nn, int);
	}

#pragma omp parallel for num_threads(gnt) if(gnt > 1) schedule(dynamic, 1)
	for (int blk = 0; blk < nblocks; blk++) {
		int tn = (gnt > 1 ? omp_get_thread_num() : 0);
		double *bin = gbuf[tn];
		double *bwork = gbuf[tn] + (size_t) GB * nn;
		int c0 = blk * GB;
		int nb = IMIN(GB, wlist->n - c0);

		// scatter the sparse A-rows into the mapped positions; the forward
		// solve only fills indices >= the first nonzero, so remember the
		// per-column start to limit the compression scan below
		int jm[nb];
		GMRFLib_dfill(nb * nn, 0.0, bin);
		int *rmap = pb->sub_sm_fact.remap;
		for (int c = 0; c < nb; c++) {
			GMRFLib_idxval_tp *v = Aidx[wlist->idx[c0 + c]];
			double *bc = bin + (size_t) c * nn;
			int jm_c = nn;
			for (int k = 0; k < v->n; k++) {
				int jj = rmap[v->idx[k]];
				bc[jj] = v->val[k];
				jm_c = IMIN(jm_c, jj);
			}
			jm[c] = jm_c;
		}

		GMRFLib_taucs_Lsolve_blocked(pb->sub_sm_fact.TAUCS_L, bin, nb, bwork);

		for (int c = 0; c < nb; c++) {
			int node = wlist->idx[c0 + c];
			double *b1 = bin + (size_t) c * nn;
			int *is = gis[tn];
			double *vs = bwork;		       /* free after the solve: reuse as scratch */
			int cnt = 0;
			for (int j = jm[c]; j < nn; j++) {
				if (b1[j] != 0.0) {
					is[cnt] = j;
					vs[cnt] = b1[j];
					cnt++;
				}
			}
			widx[node] = Malloc(IMAX(1, cnt), int);
			wval[node] = Malloc(IMAX(1, cnt), double);
			Memcpy(widx[node], is, IMAX(1, cnt) * sizeof(int));
			Memcpy(wval[node], vs, IMAX(1, cnt) * sizeof(double));
			wnnz[node] = cnt;

			if (nc > 0) {
				GMRFLib_idxval_tp *v = Aidx[node];
				for (int cc = 0; cc < nc; cc++) {
					double pp = 0.0;
					double *Vc = Vg + (size_t) cc * nn;
					for (int k = 0; k < cnt; k++) {
						pp += vs[k] * Vc[is[k]];
					}
					zc[(size_t) node * nc + cc] = pp;

					double mm = 0.0;
					double *Mc = pb->constr_m + (size_t) cc * nn;
					for (int k = 0; k < v->n; k++) {
						mm += v->val[k] * Mc[v->idx[k]];
					}
					mc[(size_t) node * nc + cc] = mm;
				}
			}
		}
	}

#pragma omp parallel for num_threads(gnt) if(gnt > 1) schedule(dynamic, 8)
	for (int i = 0; i < pair_fb->n; i++) {
		int node = pair_fb->idx[0][i];
		int k = pair_fb->idx[1][i];
		int nnode = groups->missing[node]->idx[0][k];
		double sum = 0.0;
		int a = 0, b = 0, na = wnnz[node], nb = wnnz[nnode];
		int *ia = widx[node], *ib = widx[nnode];
		double *va = wval[node], *vb = wval[nnode];
		while (a < na && b < nb) {
			if (ia[a] == ib[b]) {
				sum += va[a] * vb[b];
				a++;
				b++;
			} else if (ia[a] < ib[b]) {
				a++;
			} else {
				b++;
			}
		}
		for (int cc = 0; cc < nc; cc++) {
			sum -= mc[(size_t) nnode * nc + cc] * zc[(size_t) node * nc + cc];
		}
		gcpo_cov_set_(gcpo, groups->missing[node]->idx[1][k], node, nnode, sum, sd, lpred_variance);
	}

	if (timing) {
		size_t wtot = 0;
		for (int i = 0; i < wlist->n; i++) {
			wtot += (size_t) wnnz[wlist->idx[i]];
		}
		printf("[gcpo-timing] gcpo: GRAM half-solve %.4f s (%1d w-columns, %1d fallback-pairs, avg w-nnz %.0f)\n",
		       GMRFLib_timer() - tref, wlist->n, pair_fb->n, (double) wtot / IMAX(1, wlist->n));
	}

	for (int i = 0; i < wlist->n; i++) {
		Free(widx[wlist->idx[i]]);
		Free(wval[wlist->idx[i]]);
	}
	Free(widx);
	Free(wval);
	Free(wnnz);
	for (int i = 0; i < gnt; i++) {
		Free(gbuf[i]);
		Free(gis[i]);
	}
	Free(gbuf);
	Free(gis);
	Free(needw);
	Free(Vg);
	Free(zc);
	Free(mc);
	GMRFLib_idx_free(wlist);
}

// the full-solve path (non-TAUCS solvers): one column Sigma a_node per node,
// dotted against the group partners
void GMRFLib_gcpo_cov_solve(GMRFLib_ai_store_tp *ai_store_id, GMRFLib_preopt_tp *preopt, GMRFLib_gcpo_groups_tp *groups, GMRFLib_idx_tp *node_idx,
			    GMRFLib_gcpo_elm_tp **gcpo, double *sd, double *lpred_variance, unsigned char *skip, int n_skip,
			    int nt_outer, int nt_inner, int serial, GMRFLib_gcpo_param_tp *gcpo_param, int detailed_output, int gcpo_timing)
{
#define A_idx(node_) (preopt->pAA_idxval ? preopt->pAA_idxval[node_] : preopt->A_idxval[node_])
	int nn = preopt->n;
	int use_stiles = (GMRFLib_smtp == GMRFLib_SMTP_STILES);
	int nrhs = 1;
	if (use_stiles) {
		nrhs = IMIN(GMRFLib_stiles_get_block_size(), 1 + node_idx->n / nt_outer);
	} else {
		nrhs = IMIN(GMRFLib_taucs_get_block_size(), 1 + node_idx->n / nt_outer);
	}

	int Swork_len = nt_inner;
	double **Swork = NULL;
	Swork = Malloc(Swork_len, double *);
	for (int i = 0; i < Swork_len; i++) {
		Swork[i] = Malloc(nn * nrhs, double);
	}

	GMRFLib_ptr_tp *split = GMRFLib_idx_split(node_idx, nrhs);

	int use_group = 0;
	int tnum_parent = omp_get_thread_num();
	if (serial && use_stiles) {
		GMRFLib_stiles_rescale_start(1);
		use_group = GMRFLib_stiles_rescale_group();
	}

	double gcpo_tref = (gcpo_timing ? GMRFLib_timer() : 0.0);

	int run_parallel = !use_stiles || (use_stiles && serial);
#pragma omp parallel for num_threads(nt_inner) if(run_parallel) schedule(static)
	for (int kk = 0; kk < split->n; kk++) {

		int tnum = omp_get_thread_num();
		GMRFLib_idx_tp *lnode_idx = (GMRFLib_idx_tp *) split->ptr[kk];
		GMRFLib_stiles_idx_tp stiles_idx = { use_group, (serial ? tnum : tnum_parent), lnode_idx->n };
		if (GMRFLib_smtp == GMRFLib_SMTP_STILES) {
			// bind is set only if its not already set
			GMRFLib_stiles_set_idx(&stiles_idx, lnode_idx->n);
			GMRFLib_stiles_bind(&stiles_idx);
		}

		double *Saa = Swork[tnum];
		GMRFLib_dfill(nn * nrhs, 0.0, Saa);

		for (int inode = 0; inode < lnode_idx->n; inode++) {
			int node = lnode_idx->idx[inode];
			GMRFLib_idxval_tp *v = A_idx(node);
			GMRFLib_unpack(v->n, v->val, Saa + inode * nn, v->idx);
		}

		GMRFLib_Qsolves(Saa, lnode_idx->n, ai_store_id->problem, &stiles_idx);

		for (int inode = 0; inode < lnode_idx->n; inode++) {
			int node = lnode_idx->idx[inode];
			if (gcpo_param->verbose || detailed_output) {
				if (skip[node]) {
					printf("%s[%1d]: Skip solve for node %d\n", __GMRFLib_FuncName, omp_get_thread_num(), node);
				} else {
					printf("%s[%1d]: Solve for node %d\n", __GMRFLib_FuncName, omp_get_thread_num(), node);
				}
			}
			for (int k = 0; k < groups->missing[node]->n; k++) {
				int nnode = groups->missing[node]->idx[0][k];
				double cov = (nnode != node ? GMRFLib_sparse_ddot_(A_idx(nnode), Saa + inode * nn) : 0.0);
				gcpo_cov_set_(gcpo, groups->missing[node]->idx[1][k], node, nnode, cov, sd, lpred_variance);
			}
		}
	}

	if (gcpo_timing) {
		printf("[gcpo-timing] gcpo: cov solve loop %.4f s for %1d columns (nrhs %1d, nt_inner %1d, skip %1d)\n",
		       GMRFLib_timer() - gcpo_tref, node_idx->n, nrhs, nt_inner, n_skip);
	}
	if (serial && use_stiles) {
		// this will also do unbind()
		GMRFLib_stiles_rescale_end();
	}

	for (int i = 0; i < Swork_len; i++) {
		Free(Swork[i]);
	}
	Free(Swork);
	GMRFLib_idx_split_free(split);
#undef A_idx
}


// ---------------------------------------------------------------------------
// the stages as called from GMRFLib_gcpo_build() and GMRFLib_gcpo()
// ---------------------------------------------------------------------------

// gate + stage 1. returns 1 when the local build is on
int GMRFLib_gcpo_local_prepare(GMRFLib_gcpo_local_tp *rb, GMRFLib_ai_store_tp *build_ai_store, GMRFLib_idxval_tp **Aidx,
			       GMRFLib_gcpo_param_tp *gcpo_param, GMRFLib_idx_tp *d_idx, int Npred)
{
	memset(rb, 0, sizeof(GMRFLib_gcpo_local_tp));
	// radius build (gcpo_param->build_radius = r > 0; INLA_GCPO_BUILD_RADIUS
	// overrides it for experiments): walk the factor graph
	// r hops out of each data-node's A-support (hubs do not expand), take the
	// data-nodes touching the walk as candidates and demand their support-
	// product pairs into the Qinv-store computed just below (which is computed
	// anyway for the sd's: no extra Takahashi pass). group formation then reads
	// the correlations from the store, and the separator certificate proves
	// that no data-node outside the candidates can enter the group (or, with
	// a single hub sealing the core, the far field is known exactly); a node
	// that cannot be concluded goes to the solve-loop -- no retries, the
	// build must stay cheaper than the solves it replaces.
	// INLA_GCPO_BUILD_AB=1 runs both paths and compares.
	// constrained problems are supported: GMRFLib_Qsolves and
	// GMRFLib_compute_Qinv both apply the constraint correction, so the hub
	// columns, the candidate correlations and the sds are already
	// constrained; the certificate conditions on the active constraints
	// locally, and the separator equality is disabled for components
	// touched by a constraint
		rb->verbose = (gcpo_param->verbose || getenv("INLA_GCPO_TIMING") != NULL);
	rb->radius = IMAX(0, gcpo_param->build_radius);
	if (getenv("INLA_GCPO_BUILD_RADIUS")) {
		rb->radius = IMAX(0, atoi(getenv("INLA_GCPO_BUILD_RADIUS")));
	}
	if (!(GMRFLib_smtp == GMRFLib_SMTP_TAUCS && !(gcpo_param->friends) && gcpo_param->num_level_sets != -1 && d_idx)) {
		rb->radius = 0;		       /* the radius build needs the TAUCS Qinv store and plain level-set groups */
	}
	rb->ab = (rb->radius > 0 && getenv("INLA_GCPO_BUILD_AB") != NULL);
	if (rb->radius == 0 && rb->verbose) {
		printf("[gcpo-timing] build: radius build off (build.radius %1d, smtp %s, friends %s, num.level.sets %1d)\n",
		       gcpo_param->build_radius, (GMRFLib_smtp == GMRFLib_SMTP_TAUCS ? "taucs" : "other"),
		       (gcpo_param->friends ? "yes" : "no"), gcpo_param->num_level_sets);
	}
	if (rb->radius > 0 && gcpo_param->any_rankdef) {
		// unconstrained intrinsic components: the posterior is only weakly
		// identified, the two build paths' cor values then differ beyond
		// the equal_cor band and group identity is not fp-well-defined in
		// ANY implementation -- refuse the fast path entirely
		if (rb->verbose) {
		printf("[gcpo-timing] build: radius-lookup disabled (rank-deficient component without constraint)\n");
		}
		rb->radius = 0;
	}
	if (rb->radius > 0) {
		rb_setup_(rb, build_ai_store->problem, Aidx, gcpo_param, d_idx, Npred);
	}
	return (rb->radius > 0);
}

// stage 2 (after the Qinv store and the predictor sd's exist)
void GMRFLib_gcpo_local_hubs(GMRFLib_gcpo_local_tp *rb, GMRFLib_problem_tp *pb, GMRFLib_idxval_tp **Aidx, GMRFLib_idx_tp *d_idx, double *isd, int mnpred)
{
	if (rb->radius > 0) {
		rb_hub_columns_(rb, pb, Aidx, d_idx, isd, mnpred);
	}
}

// stage 3. returns the nodes the solve loop has to do: the unconcluded ones, or
// the whole selection (local build off, or A/B mode where the solve loop redoes
// everything and *groups_rb keeps the local groups for the comparison)
GMRFLib_idx_tp *GMRFLib_gcpo_local_lookup(int thread_id, GMRFLib_gcpo_local_tp *rb, GMRFLib_ai_store_tp *build_ai_store, GMRFLib_idxval_tp **Aidx,
					  GMRFLib_gcpo_param_tp *gcpo_param, GMRFLib_idx_tp *d_idx, GMRFLib_idx_tp *selection, double *isd, double min_sd,
					  GMRFLib_idxval_tp **groups, int nt_outer, int Npred, GMRFLib_idxval_tp ***groups_rb)
{
	*groups_rb = NULL;
	if (rb->radius == 0) {
		return selection;
	}
	rb->certflag = (rb->ab ? Calloc(Npred, char) : NULL);
	rb->fb_sel = rb_lookup_(thread_id, rb, build_ai_store->problem, Aidx, gcpo_param, d_idx, selection, isd, min_sd, groups, nt_outer);
	if (!rb->ab) {
		return rb->fb_sel;
	}
	// A/B: save the local groups and let the solve loop redo ALL nodes
	*groups_rb = Calloc(Npred, GMRFLib_idxval_tp *);
	for (int is = 0; is < selection->n; is++) {
		int node = selection->idx[is];
		if (groups[node] && groups[node]->n > 0) {
			GMRFLib_idxval_create_x(&((*groups_rb)[node]), groups[node]->n);
			for (int i = 0; i < groups[node]->n; i++) {
				GMRFLib_idxval_add(&((*groups_rb)[node]), groups[node]->idx[i], groups[node]->val[i]);
			}
		}
		groups[node]->n = 0;
	}
	return selection;
}

// after the solve loop: the A/B comparison, then everything is freed
void GMRFLib_gcpo_local_finish(GMRFLib_gcpo_local_tp *rb, GMRFLib_idxval_tp **groups_rb, GMRFLib_idxval_tp **groups, GMRFLib_idx_tp *selection, int Npred)
{
	if (rb->radius == 0) {
		return;
	}
	if (rb->ab && groups_rb) {
		rb_ab_compare_(rb, groups_rb, groups, selection);
		for (int i = 0; i < Npred; i++) {
			GMRFLib_idxval_free(groups_rb[i]);
		}
		Free(groups_rb);
	}
	GMRFLib_idx_free(rb->fb_sel);
	rb_free_(rb, Npred);
}

// the end of the build: the in-group pairs become the Qinv demand set for the
// stores the configurations compute next (one Takahashi per configuration, for
// the predictor variances; the pairs ride along). the store computed during the
// build is kept for the mode: it holds every pair of a certified group (group
// subset-of candidates) but not necessarily those of a fallback node's group,
// which sit far apart in the elimination order -- if any is missing the store
// is recomputed once with the union (the closure budget in smtp-taucs.c bounds
// that pass; pairs it cannot close are solved for in GMRFLib_gcpo_cov_pairs)
void GMRFLib_gcpo_local_demands(GMRFLib_ai_store_tp *ai_store, GMRFLib_idxval_tp **Aidx, GMRFLib_idx2_tp **missing, int Npred, int verbose)
{
	if (GMRFLib_smtp != GMRFLib_SMTP_TAUCS) {
		return;
	}
	int nlatent = ai_store->problem->sub_graph->n;
	size_t kp_n = 0, miss = 0;
	GMRFLib_idx_tp **kp = gcpo_demand_set_(Aidx, missing, Npred, ai_store->problem->sub_graph, verbose, &kp_n);
	if (ai_store->problem->sub_inverse) {
		for (int i = 0; i < nlatent; i++) {
			for (int k = 0; kp[i] && k < kp[i]->n; k++) {
				miss += (GMRFLib_Qinv_get(ai_store->problem, i, kp[i]->idx[k]) == NULL);
			}
		}
	} else {
		miss = kp_n + 1;
	}
	if (miss > 0) {
		double tref = GMRFLib_timer();
		GMRFLib_idx_tp **un = Calloc(nlatent, GMRFLib_idx_tp *);
		for (int i = 0; i < nlatent; i++) {
			if (GMRFLib_qinv_keep_pairs && i < GMRFLib_qinv_keep_pairs_n && GMRFLib_qinv_keep_pairs[i]) {
				GMRFLib_idx_nadd(&(un[i]), GMRFLib_qinv_keep_pairs[i]->n, GMRFLib_qinv_keep_pairs[i]->idx);
			}
			if (kp[i]) {
				GMRFLib_idx_nadd(&(un[i]), kp[i]->n, kp[i]->idx);
			}
			if (un[i]) {
				GMRFLib_idx_sort(un[i]);
				GMRFLib_idx_uniq(un[i]);
			}
		}
		gcpo_install_demand_(un, nlatent);
		GMRFLib_free_Qinv(ai_store->problem);
		GMRFLib_ai_add_Qinv_to_ai_store(ai_store);
		if (verbose) {
			printf("[gcpo-timing] build: store recomputed with the union, %zu of %zu group pairs were missing (%.4f s)\n", miss, kp_n,
			       GMRFLib_timer() - tref);
		}
	}
	gcpo_install_demand_(kp, nlatent);	       /* the configurations to come need the in-group pairs */
}

// per configuration, TAUCS: lookup first; the pairs the store misses go through
// the Gram half-solves, or through the plain solves when the half-solved
// columns would not fit in memory (GCPO_GRAM_MEM_MAX bytes)
#define GCPO_GRAM_MEM_MAX ((size_t) 2 << 30)

void GMRFLib_gcpo_cov_pairs(GMRFLib_ai_store_tp *ai_store_id, GMRFLib_preopt_tp *preopt, GMRFLib_gcpo_groups_tp *groups, GMRFLib_idx_tp *node_idx,
			    GMRFLib_gcpo_elm_tp **gcpo, double *sd, double *lpred_variance, unsigned char *skip, int n_skip, int nt_outer, int nt_inner,
			    int serial, GMRFLib_gcpo_param_tp *gcpo_param, int detailed_output, int timing)
{
	GMRFLib_idxval_tp **Aidx = (preopt->pAA_idxval ? preopt->pAA_idxval : preopt->A_idxval);
	GMRFLib_problem_tp *pb = ai_store_id->problem;
	int Npred = preopt->Npred, nn = preopt->n;

	GMRFLib_ai_add_Qinv_to_ai_store(ai_store_id);	       /* no-op: the configuration computed it for the variances */
	GMRFLib_idx2_tp *pair_fb = GMRFLib_gcpo_cov_lookup(pb, Aidx, groups, node_idx, gcpo, sd, lpred_variance, timing);
	if (pair_fb) {
		// the nodes of the missed pairs, and the memory their half-solved columns would take
		unsigned char *needw = Calloc(Npred, unsigned char);
		size_t nw = 0;
		for (int i = 0; i < pair_fb->n; i++) {
			int node = pair_fb->idx[0][i];
			int nnode = groups->missing[node]->idx[0][pair_fb->idx[1][i]];
			nw += !needw[node];
			needw[node] = 1;
			nw += !needw[nnode];
			needw[nnode] = 1;
		}
		taucs_ccs_matrix *L = (taucs_ccs_matrix *) pb->sub_sm_fact.TAUCS_L;
		size_t nnz_per_col = (L && nn > 0 ? (size_t) L->colptr[nn] / (size_t) nn : 0);
		size_t mem = nw * nnz_per_col * (sizeof(int) + sizeof(double));
		if (mem <= GCPO_GRAM_MEM_MAX) {
			GMRFLib_gcpo_cov_gram(pb, Aidx, groups, pair_fb, gcpo, sd, lpred_variance, Npred, nt_inner, timing);
		} else {
			GMRFLib_idx_tp *wlist = NULL;
			for (int i = 0; i < Npred; i++) {
				if (needw[i]) {
					GMRFLib_idx_add(&wlist, i);
				}
			}
			if (timing) {
				printf("[gcpo-timing] gcpo: %zu half-solved columns would take %.1f GB: plain solves for %1d nodes instead\n", nw,
				       (double) mem / 1073741824.0, wlist->n);
			}
			GMRFLib_gcpo_cov_solve(ai_store_id, preopt, groups, wlist, gcpo, sd, lpred_variance, skip, n_skip, nt_outer, nt_inner, serial, gcpo_param,
					       detailed_output, timing);
			GMRFLib_idx_free(wlist);
		}
		Free(needw);
		GMRFLib_idx2_free(pair_fb);
	}
	if (getenv("INLA_GCPO_COV_AB")) {
		// debug: the stored pair covariances against direct solves, on a sample of nodes
		int nn_ = nn, ns = IMIN(256, node_idx->n), step = IMAX(1, node_idx->n / IMAX(1, ns));
		double *Sa = Malloc(nn_, double), maxabs = 0.0, maxrel = 0.0;
		size_t npairs = 0;
		GMRFLib_stiles_idx_tp sx = { 0, -1, 1 };
		for (int k = 0; k < ns; k++) {
			int node = node_idx->idx[k * step];
			GMRFLib_idxval_tp *v = Aidx[node];
			GMRFLib_dfill(nn_, 0.0, Sa);
			GMRFLib_unpack(v->n, v->val, Sa, v->idx);
			GMRFLib_Qsolves(Sa, 1, ai_store_id->problem, &sx);
			for (int kk = 0; kk < groups->missing[node]->n; kk++) {
				int nnode = groups->missing[node]->idx[0][kk];
				int cm = groups->missing[node]->idx[1][kk];
				if (nnode == node) {
					continue;
				}
				double ref = GMRFLib_sparse_ddot_(Aidx[nnode], Sa);
				int ii = GMRFLib_iwhich_sorted(node, (int *) gcpo[cm]->idxs->idx, (unsigned int) gcpo[cm]->idxs->n);
				int jj = GMRFLib_iwhich_sorted(nnode, (int *) gcpo[cm]->idxs->idx, (unsigned int) gcpo[cm]->idxs->n);
				double got = gsl_matrix_get(gcpo[cm]->cov_mat, ii, jj);
				double dd = ABS(got - ref);
				maxabs = DMAX(maxabs, dd);
				maxrel = DMAX(maxrel, dd / DMAX(ABS(ref), 1.0e-300));
				npairs++;
			}
		}
		printf("[gcpo-timing] gcpo: COV-AB %zu pairs on %1d nodes: max |store - solve| %.3e, max relative %.3e\n", npairs, ns, maxabs, maxrel);
		Free(Sa);
	}
}
