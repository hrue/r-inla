#include <math.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>

#define NDEBUG
#include <assert.h>

#include "taucs.h"

#define FALSE 0
#define TRUE  1

typedef struct {
	int sn_size;
	int n;
	int *rowind;

	int up_size;
	int *sn_vertices;
	int *up_vertices;

	double *f1;
	double *f2;
	double *u;

} supernodal_frontal_matrix;

#define SFM_F1 f1
#define SFM_F2 f2
#define SFM_U   u

typedef struct {
	int flags;

	char uplo;					       /* 'u' for upper, 'l' for lower, ' ' don't know; prefer lower. */
	int n;						       /* size of matrix */
	int n_sn;					       /* number of supernodes */

	int *parent;					       /* supernodal elimination tree */
	int *first_child;
	int *next_child;

	int *sn_size;					       /* size of supernodes (diagonal block) */
	int *sn_up_size;				       /* size of subdiagonal update blocks */
	int **sn_struct;				       /* row structure of supernodes */

	int *sn_blocks_ld;				       /* lda of supernode blocks */
	double **sn_blocks;				       /* supernode blocks */

	int *up_blocks_ld;				       /* lda of update blocks */
	double **up_blocks;				       /* update blocks */
} supernodal_factor_matrix;

#ifdef TAUCS_CORE_GENERAL

static int *compare_indirect_map;
static int compare_indirect_ints(const void *vx, const void *vy)
{
	int *ix = (int *) vx;
	int *iy = (int *) vy;

	if (compare_indirect_map[*ix] < compare_indirect_map[*iy])
		return -1;
	if (compare_indirect_map[*ix] > compare_indirect_map[*iy])
		return 1;
	return 0;
}

#endif							       /* TAUCS_CORE_GENERAL */


#ifndef TAUCS_CORE_GENERAL

static supernodal_factor_matrix *multifrontal_supernodal_create()
{
	supernodal_factor_matrix *L;

	L = (supernodal_factor_matrix *) taucs_malloc(sizeof(supernodal_factor_matrix));
	if (!L)
		return NULL;

	L->flags = TAUCS_DOUBLE;
	L->uplo = 'l';
	L->n = -1;					       /* unused */

	L->sn_struct = NULL;
	L->sn_size = NULL;
	L->sn_up_size = NULL;
	L->parent = NULL;
	L->first_child = NULL;
	L->next_child = NULL;
	L->sn_blocks_ld = NULL;
	L->sn_blocks = NULL;
	L->up_blocks_ld = NULL;
	L->up_blocks = NULL;

	return L;
}

void taucs_dtl(supernodal_factor_free) (void *vL) {
	supernodal_factor_matrix * L = (supernodal_factor_matrix *) vL;
	int sn;

	if (!L)
		return;

	taucs_free(L->parent);
	taucs_free(L->first_child);
	taucs_free(L->next_child);

	taucs_free(L->sn_size);
	taucs_free(L->sn_up_size);
	taucs_free(L->sn_blocks_ld);
	taucs_free(L->up_blocks_ld);

	if (L->sn_struct)
		for (sn = 0; sn < L->n_sn; sn++)
			taucs_free(L->sn_struct[sn]);

	if (L->sn_blocks)
		for (sn = 0; sn < L->n_sn; sn++)
			taucs_free(L->sn_blocks[sn]);

	if (L->up_blocks)
		for (sn = 0; sn < L->n_sn; sn++)
			taucs_free(L->up_blocks[sn]);

	taucs_free(L->sn_struct);
	taucs_free(L->sn_blocks);
	taucs_free(L->up_blocks);

	taucs_free(L);
}

void taucs_dtl(supernodal_factor_free_numeric) (void *vL) {
	supernodal_factor_matrix * L = (supernodal_factor_matrix *) vL;
	int sn;

	for (sn = 0; sn < L->n_sn; sn++) {
		taucs_free(L->sn_blocks[sn]);
		L->sn_blocks[sn] = NULL;
		taucs_free(L->up_blocks[sn]);
		L->up_blocks[sn] = NULL;
	}
}

taucs_ccs_matrix *taucs_dtl(supernodal_factor_to_ccs) (void *vL) {
	supernodal_factor_matrix * L = (supernodal_factor_matrix *) vL;
	taucs_ccs_matrix *C;
	int n, nnz;
	int i, j, ip, jp, sn, next;
	double v;
	int *len;

	n = L->n;

	len = (int *) taucs_malloc(n * sizeof(int));
	if (!len)
		return NULL;

	nnz = 0;
	/*
	 * for (sn=0; sn<L->n_sn; sn++) { for (jp=0; jp<(L->sn_size)[sn]; jp++) { j = (L->sn_struct)[sn][jp]; len[j] =
	 * (L->sn_up_size)[sn] - jp; nnz += len[j]; } } 
	 */

	for (sn = 0; sn < L->n_sn; sn++) {
		for (jp = 0; jp < (L->sn_size)[sn]; jp++) {
			j = (L->sn_struct)[sn][jp];
			len[j] = 0;

			for (ip = jp; ip < (L->sn_size)[sn]; ip++) {
				i = (L->sn_struct)[sn][ip];
				v = (L->sn_blocks)[sn][jp * (L->sn_blocks_ld)[sn] + ip];

				if (taucs_re(v) || taucs_im(v)) {
					len[j]++;
					nnz++;
				}
			}
			for (ip = (L->sn_size)[sn]; ip < (L->sn_up_size)[sn]; ip++) {
				i = (L->sn_struct)[sn][ip];
				v = (L->up_blocks)[sn][jp * (L->up_blocks_ld)[sn] + (ip - (L->sn_size)[sn])];

				if (taucs_re(v) || taucs_im(v)) {
					len[j]++;
					nnz++;
				}
			}
		}
	}

	C = taucs_dtl(ccs_create) (n, n, nnz);
	if (!C) {
		taucs_free(len);
		return NULL;
	}
	C->flags = TAUCS_DOUBLE;
	C->flags |= TAUCS_TRIANGULAR | TAUCS_LOWER;

	(C->colptr)[0] = 0;
	for (j = 1; j <= n; j++)
		(C->colptr)[j] = (C->colptr)[j - 1] + len[j - 1];

	taucs_free(len);

	for (sn = 0; sn < L->n_sn; sn++) {
		for (jp = 0; jp < (L->sn_size)[sn]; jp++) {
			j = (L->sn_struct)[sn][jp];

			next = (C->colptr)[j];

			/*
			 * memcpy((C->rowind) + next, ((L->sn_struct)[sn]) + jp, ((L->sn_up_size)[sn] - jp) * sizeof(int));
			 * memcpy(C->values + next, ((L->sn_blocks)[sn]) + (jp*(L->sn_blocks_ld)[sn] + jp),
			 * ((L->sn_size)[sn] - jp) * sizeof(double)); next += ((L->sn_size)[sn] - jp);
			 * memcpy(C->values + next, ((L->up_blocks)[sn]) + jp*(L->up_blocks_ld)[sn], ((L->sn_up_size)[sn] - 
			 * (L->sn_size)[sn]) * sizeof(double)); 
			 */

			for (ip = jp; ip < (L->sn_size)[sn]; ip++) {
				i = (L->sn_struct)[sn][ip];
				v = (L->sn_blocks)[sn][jp * (L->sn_blocks_ld)[sn] + ip];

				if (!taucs_re(v) && !taucs_im(v))
					continue;
				/*
				 * if (v == 0.0) continue;
				 */

				(C->rowind)[next] = i;
				C->values[next] = v;
				next++;
			}
			for (ip = (L->sn_size)[sn]; ip < (L->sn_up_size)[sn]; ip++) {
				i = (L->sn_struct)[sn][ip];
				v = (L->up_blocks)[sn][jp * (L->up_blocks_ld)[sn] + (ip - (L->sn_size)[sn])];

				if (!taucs_re(v) && !taucs_im(v))
					continue;
				/*
				 * if (v == 0.0) continue;
				 */

				(C->rowind)[next] = i;
				C->values[next] = v;
				next++;
			}
		}
	}

	return C;
}

double *taucs_dtl(supernodal_factor_get_diag) (void *vL) {
	supernodal_factor_matrix * L = (supernodal_factor_matrix *) vL;
	int j, ip, jp, sn;				       /* i,next omer */
	double v;
	double *diag;

	diag = (double *) taucs_malloc((L->n) * sizeof(double));
	if (!diag)
		return NULL;

	for (sn = 0; sn < L->n_sn; sn++) {
		for (jp = 0; jp < (L->sn_size)[sn]; jp++) {
			j = (L->sn_struct)[sn][jp];
			ip = jp;			       /* we just want the diagonal */
			v = (L->sn_blocks)[sn][jp * (L->sn_blocks_ld)[sn] + ip];
			diag[j] = v;
		}
	}

	return diag;
}

static supernodal_frontal_matrix *supernodal_frontal_create(int *UNUSED(firstcol_in_supernode), int sn_size, int n, int *rowind)
{
	supernodal_frontal_matrix *tmp;

	tmp = (supernodal_frontal_matrix *) taucs_malloc(sizeof(supernodal_frontal_matrix));
	if (tmp == NULL)
		return NULL;

	tmp->sn_size = sn_size;
	tmp->n = n;

	tmp->rowind = rowind;

	tmp->n = n;
	tmp->sn_size = sn_size;
	tmp->up_size = n - sn_size;

	tmp->sn_vertices = rowind;
	tmp->up_vertices = rowind + sn_size;

	/*
	 * on some platforms, malloc(0) fails, so we avoid such calls 
	 */

	tmp->SFM_F1 = tmp->SFM_F2 = tmp->SFM_U = NULL;

	if (tmp->sn_size) {
		tmp->SFM_F1 = (double *) taucs_calloc((tmp->sn_size) * (tmp->sn_size), sizeof(double));
	}

	if (tmp->sn_size && tmp->up_size) {
		tmp->SFM_F2 = (double *) taucs_calloc((tmp->up_size) * (tmp->sn_size), sizeof(double));
	}

	if (tmp->up_size) {
		tmp->SFM_U = (double *) taucs_calloc((tmp->up_size) * (tmp->up_size), sizeof(double));
	}

	if ((tmp->SFM_F1 == NULL && tmp->sn_size)
	    || (tmp->SFM_F2 == NULL && tmp->sn_size && tmp->up_size)
	    || (tmp->SFM_U == NULL && tmp->up_size)) {
		taucs_free(tmp->SFM_U);
		taucs_free(tmp->SFM_F1);
		taucs_free(tmp->SFM_F2);
		taucs_free(tmp);
		return NULL;
	}

	assert(tmp);
	return tmp;
}

static void supernodal_frontal_free(supernodal_frontal_matrix *to_del)
{
	if (to_del) {
		taucs_free(to_del->SFM_F1);
		taucs_free(to_del->SFM_F2);
		taucs_free(to_del->SFM_U);
		taucs_free(to_del);
	}
}

static int
multifrontal_supernodal_front_factor(int sn,
				     int *firstcol_in_supernode,
				     int sn_size, taucs_ccs_matrix *A, supernodal_frontal_matrix *mtr, int *bitmap, supernodal_factor_matrix *snL)
{
	int i, j;
	int *ind;
	double *re;
	int INFO;

	for (i = 0; i < mtr->sn_size; i++)
		bitmap[mtr->sn_vertices[i]] = i;
	for (i = 0; i < mtr->up_size; i++)
		bitmap[mtr->up_vertices[i]] = mtr->sn_size + i;

	for (j = 0; j < (mtr->sn_size); j++) {
		ind = &(A->rowind[A->colptr[*(firstcol_in_supernode + j)]]);
		re = &(A->values[A->colptr[*(firstcol_in_supernode + j)]]);
		for (i = 0; i < A->colptr[*(firstcol_in_supernode + j) + 1]
		     - A->colptr[*(firstcol_in_supernode + j)]; i++) {
			if (bitmap[ind[i]] < mtr->sn_size)
				mtr->SFM_F1[(mtr->sn_size) * j + bitmap[ind[i]]] =
				    taucs_add(mtr->SFM_F1[(mtr->sn_size) * j + bitmap[ind[i]]], re[i]);
			else
				mtr->SFM_F2[(mtr->up_size) * j + bitmap[ind[i]] - mtr->sn_size] =
				    taucs_add(mtr->SFM_F2[(mtr->up_size) * j + bitmap[ind[i]] - mtr->sn_size], re[i]);
		}
	}

	if (mtr->sn_size) {
		taucs_potrf("L", &(mtr->sn_size), mtr->SFM_F1, &(mtr->sn_size), &INFO, F_ONE);
	}

	if (mtr->up_size && mtr->sn_size) {
		taucs_trsm("R",
			   "L",
			   "C",
			   "N",
			   &(mtr->up_size), &(mtr->sn_size),
			   &taucs_one_const, mtr->SFM_F1, &(mtr->sn_size), mtr->SFM_F2, &(mtr->up_size), F_ONE, F_ONE, F_ONE, F_ONE);
	}

	(snL->sn_blocks)[sn] = mtr->SFM_F1;
	(snL->sn_blocks_ld)[sn] = mtr->sn_size;

	(snL->up_blocks)[sn] = mtr->SFM_F2;
	(snL->up_blocks_ld)[sn] = mtr->up_size;

	if (mtr->up_size && mtr->sn_size) {
		taucs_herk("L",
			   "N",
			   &(mtr->up_size), &(mtr->sn_size),
			   &taucs_minusone_real_const,
			   mtr->SFM_F2, &(mtr->up_size), &taucs_one_real_const, mtr->SFM_U, &(mtr->up_size), F_ONE, F_ONE);
	}

	mtr->SFM_F1 = NULL;				       /* so we don't free twice */
	mtr->SFM_F2 = NULL;				       /* so we don't free twice */

	return 0;
}

static void multifrontal_supernodal_front_extend_add(supernodal_frontal_matrix *parent_mtr, supernodal_frontal_matrix *my_mtr, int *bitmap)
{
	int j, i, parent_i;
	double v;

	int sn_size = parent_mtr->sn_size;
	int up_size = parent_mtr->up_size;

	double *F1 = parent_mtr->SFM_F1;
	double *F2 = parent_mtr->SFM_F2;
	double *U = parent_mtr->SFM_U;

	for (i = 0; i < parent_mtr->sn_size; i++)
		bitmap[parent_mtr->sn_vertices[i]] = i;

	for (i = 0; i < parent_mtr->up_size; i++)
		bitmap[parent_mtr->up_vertices[i]] = sn_size + i;

	if (1) {
		for (j = 0; j < my_mtr->up_size; j++) {
			double *vv = &((my_mtr->SFM_U)[my_mtr->up_size * j]);
			int parent_jj = bitmap[my_mtr->up_vertices[j]];
			int smaller = (parent_jj < sn_size);
			int n1 = sn_size * parent_jj;
			int n2 = up_size * parent_jj;
			int off_ii = parent_jj - sn_size;
			double *F11 = F1 + n1;
			double *FF11 = F1 + parent_jj;
			double *F22 = F2 + n2;
			double *FF22 = F2 + off_ii;
			double *U1 = U + up_size * (parent_jj - sn_size);
			double *U2 = U + off_ii;

			if (smaller) {
				for (i = j; i < my_mtr->up_size; i++) {
					v = vv[i];
					parent_i = bitmap[my_mtr->up_vertices[i]];

					if (parent_jj > parent_i) {
						if (parent_i < sn_size) {
							FF11[sn_size * parent_i] += v;
						} else {
							int off_j = parent_i - sn_size;
							U2[up_size * off_j] += v;
						}

					} else {
						if (parent_i < sn_size) {
							F11[parent_i] += v;
						} else {
							int off_i = parent_i - sn_size;
							F22[off_i] += v;
						}
					}
				}
			} else {
				for (i = j; i < my_mtr->up_size; i++) {
					v = vv[i];
					parent_i = bitmap[my_mtr->up_vertices[i]];

					if (parent_jj > parent_i) {
						if (parent_i < sn_size) {
							FF22[up_size * parent_i] += v;
						} else {
							int off_j = parent_i - sn_size;
							U2[up_size * off_j] += v;
						}

					} else {
						int off_i = parent_i - sn_size;
						U1[off_i] += v;
					}
				}
			}
		}
	} else {
		for (j = 0; j < my_mtr->up_size; j++) {
			double *vv = &((my_mtr->SFM_U)[my_mtr->up_size * j]);
			int parent_jj = bitmap[my_mtr->up_vertices[j]];
			int smaller = (parent_jj < sn_size);
			int n1 = sn_size * parent_jj;
			int n2 = up_size * parent_jj;
			int off_ii = parent_jj - sn_size;
			double *F11 = F1 + n1;
			double *FF11 = F1 + parent_jj;
			double *F22 = F2 + n2;
			double *FF22 = F2 + off_ii;
			double *U1 = U + up_size * (parent_jj - sn_size);
			double *U2 = U + off_ii;

			for (i = j; i < my_mtr->up_size; i++) {
				v = vv[i];
				parent_i = bitmap[my_mtr->up_vertices[i]];

				if (parent_jj > parent_i) {
					if (parent_i < sn_size) {
						if (smaller) {
							FF11[sn_size * parent_i] += v;
						} else {
							FF22[up_size * parent_i] += v;
						}
					} else {
						int off_j = parent_i - sn_size;
						U2[up_size * off_j] += v;
					}

				} else {
					if (smaller) {
						if (parent_i < sn_size) {
							F11[parent_i] += v;
						} else {
							int off_i = parent_i - sn_size;
							F22[off_i] += v;
						}
					} else {
						int off_i = parent_i - sn_size;
						U1[off_i] += v;
					}
				}

			}
		}
	}
}


#endif							       /* #ifndef TAUCS_CORE_GENERAL */

#ifdef TAUCS_CORE_GENERAL

/* UNION FIND ROUTINES */

static int uf_makeset(int *uf, int i)
{
	uf[i] = i;
	return i;
}
static int uf_find(int *uf, int i)
{
	if (uf[i] != i)
		uf[i] = uf_find(uf, uf[i]);
	return uf[i];
}
static int uf_union(int *uf, int s, int t)
{
	if (uf_find(uf, s) < uf_find(uf, t)) {
		uf[uf_find(uf, s)] = uf_find(uf, t);
		return (uf_find(uf, t));
	} else {
		uf[uf_find(uf, s)] = uf_find(uf, t);
		return (uf_find(uf, t));
	}
}

static
void recursive_postorder(int j, int first_child[], int next_child[], int postorder[], int ipostorder[], int *next)
{
	int c;
	for (c = first_child[j]; c != -1; c = next_child[c]) {
		recursive_postorder(c, first_child, next_child, postorder, ipostorder, next);
	}
	if (postorder)
		postorder[*next] = j;
	if (ipostorder)
		ipostorder[j] = *next;
	(*next)++;
}

static int ordered_uf_makeset(int *uf, int i)
{
	uf[i] = i;
	return i;
}
static int ordered_uf_find(int *uf, int i)
{
	if (uf[i] != i)
		uf[i] = uf_find(uf, uf[i]);
	return uf[i];
}
static int ordered_uf_union(int *uf, int s, int t)
{
	assert(uf[t] == t);
	assert(uf[s] == s);
	assert(t > s);
	if (t > s) {
		uf[s] = t;
		return t;
	} else
		uf[t] = s;
	return s;
}

static void tree_level(int j, int isroot, int first_child[], int next_child[], int level[], int level_j)
{
	int c;

	if (!isroot)
		level[j] = level_j;
	for (c = first_child[j]; c != -1; c = next_child[c]) {
		tree_level(c, FALSE, first_child, next_child, level, level_j + 1);
	}
}

void tree_first_descendant(int j, int isroot, int first_child[], int next_child[], int ipostorder[], int first_descendant[])
{
	int c;
	int fd = ipostorder[j];

	for (c = first_child[j]; c != -1; c = next_child[c]) {
		tree_first_descendant(c, FALSE, first_child, next_child, ipostorder, first_descendant);
		if (first_descendant[c] < fd)
			fd = first_descendant[c];
	}
	if (!isroot)
		first_descendant[j] = fd;
}

int taucs_ccs_etree(taucs_ccs_matrix * A, int *parent, int *l_colcount, int *l_rowcount, int *l_nnz);

static int
recursive_symbolic_elimination(int j,
			       taucs_ccs_matrix *A,
			       int first_child[],
			       int next_child[],
			       int *n_sn,
			       int sn_size[],
			       int sn_up_size[],
			       int *sn_rowind[],
			       int sn_first_child[],
			       int sn_next_child[], int rowind[], int column_to_sn_map[], int map[], int do_order, int ipostorder[]
    )
{
	int i, ip, c, c_sn;
	int in_previous_sn;
	int nnz = 0;					       /* just to suppress the warning */

	for (c = first_child[j]; c != -1; c = next_child[c]) {
		if (recursive_symbolic_elimination(c, A, first_child, next_child, n_sn, sn_size, sn_up_size, sn_rowind, sn_first_child, sn_next_child, rowind,	/* temporary 
																				 */
						   column_to_sn_map, map, do_order, ipostorder)
		    == -1)
			return -1;
	}

	in_previous_sn = 1;
	if (j == A->n)
		in_previous_sn = 0;			       /* this is not a real column */
	else if (first_child[j] == -1)
		in_previous_sn = 0;			       /* this is a leaf */
	else if (next_child[first_child[j]] != -1)
		in_previous_sn = 0;			       /* more than 1 child */
	else {
		c = first_child[j];
		for (ip = (A->colptr)[j]; ip < (A->colptr)[j + 1]; ip++) {
			i = (A->rowind)[ip];
			in_previous_sn = in_previous_sn && (map[i] == c);
		}
	}

	if (in_previous_sn) {
		c = first_child[j];
		c_sn = column_to_sn_map[c];
		column_to_sn_map[j] = c_sn;
		for (ip = sn_size[c_sn]; ip < sn_up_size[c_sn]; ip++)
			if (sn_rowind[c_sn][ip] == j)
				break;
		assert(ip < sn_up_size[c_sn]);
		sn_rowind[c_sn][ip] = sn_rowind[c_sn][sn_size[c_sn]];
		sn_rowind[c_sn][sn_size[c_sn]] = j;
		for (ip = sn_size[c_sn]; ip < sn_up_size[c_sn]; ip++)
			map[sn_rowind[c_sn][ip]] = j;

		sn_size[c_sn]++;
		return 0;
	}

	if (j < A->n) {
		nnz = 1;
		rowind[0] = j;
		map[j] = j;

		for (c = first_child[j]; c != -1; c = next_child[c]) {
			c_sn = column_to_sn_map[c];
			for (ip = sn_size[c_sn]; ip < sn_up_size[c_sn]; ip++) {
				i = sn_rowind[c_sn][ip];
				if (i > j && map[i] != j) {    /* new row index */
					map[i] = j;
					rowind[nnz] = i;
					nnz++;
				}
			}
		}

		for (ip = (A->colptr)[j]; ip < (A->colptr)[j + 1]; ip++) {
			i = (A->rowind)[ip];
			if (map[i] != j) {		       /* new row index */
				map[i] = j;
				rowind[nnz] = i;
				nnz++;
			}
		}
	}

	for (c = first_child[j]; c != -1; c = next_child[c]) {
		c_sn = column_to_sn_map[c];
		/*
		 * printf("%d ",c_sn);
		 */
		if (c == first_child[j])
			sn_first_child[*n_sn] = c_sn;
		else {
			sn_next_child[c_sn] = sn_first_child[*n_sn];
			sn_first_child[*n_sn] = c_sn;
		}
	}

	if (j < A->n) {
		column_to_sn_map[j] = *n_sn;
		sn_size[*n_sn] = 1;
		sn_up_size[*n_sn] = nnz;
		sn_rowind[*n_sn] = (int *) taucs_malloc(nnz * sizeof(int));
		if (!(sn_rowind[*n_sn]))
			return -1;
		for (ip = 0; ip < nnz; ip++)
			sn_rowind[*n_sn][ip] = rowind[ip];
		if (do_order) {
			/*
			 * Sivan and Vladimir: we think that we can sort in 
			 */
			/*
			 * column order, not only in etree postorder.  
			 */
			// radix_sort(sn_rowind [*n_sn],nnz);
			// qsort(sn_rowind [*n_sn],nnz,sizeof(int),compare_ints); 
			compare_indirect_map = ipostorder;
			qsort(sn_rowind[*n_sn], nnz, sizeof(int), compare_indirect_ints);
		}
		assert(sn_rowind[*n_sn][0] == j);
		(*n_sn)++;
	}

	return 0;
}

#endif							       /* #ifdef TAUCS_CORE_GENERAL */

#ifndef TAUCS_CORE_GENERAL

static void extend_add_wrapper(supernodal_frontal_matrix *child_matrix,
			       supernodal_frontal_matrix **my_matrix_ptr,
			       int is_root, int *v, int sn_size, int sn_up_size, int *rowind, int *bitmap, int *fail)
{

	if (*fail) {
		if (*my_matrix_ptr)
			supernodal_frontal_free(*my_matrix_ptr);
		return;
	}

	if (!is_root) {
		if (!(*my_matrix_ptr)) {
			*my_matrix_ptr = supernodal_frontal_create(v, sn_size, sn_up_size, rowind);
			if (!(*my_matrix_ptr)) {
				*fail = TRUE;
				supernodal_frontal_free(child_matrix);
				return;
			}
		}
		multifrontal_supernodal_front_extend_add(*my_matrix_ptr, child_matrix, bitmap);
	}

	/*
	 * moved outside "if !is_root"; Sivan 27 Feb 2002 
	 */
	supernodal_frontal_free(child_matrix);
}

static supernodal_frontal_matrix *recursive_multifrontal_supernodal_factor_llt(int sn,	/* this supernode */
									       int is_root,	/* is v the root? */
									       int **bitmaps,
									       taucs_ccs_matrix *A, supernodal_factor_matrix *snL, int *fail)
{
	supernodal_frontal_matrix *my_matrix = NULL;
	supernodal_frontal_matrix *ret_matrix = NULL;

	int child;
	int *v;
	int sn_size;
	int *first_child = snL->first_child;
	int *next_child = snL->next_child;

	if (!is_root) {
		sn_size = snL->sn_size[sn];
		v = &(snL->sn_struct[sn][0]);
	} else {
		sn_size = -1;
		v = NULL;				       /* not used */
	}

	for (child = first_child[sn]; child != -1; child = next_child[child]) {
		ret_matrix = recursive_multifrontal_supernodal_factor_llt(child, FALSE, bitmaps, A, snL, fail);
		if (!is_root)
			extend_add_wrapper(ret_matrix, &my_matrix, is_root, v, sn_size, snL->sn_up_size[sn], snL->sn_struct[sn], bitmaps[0], fail);
		else
			extend_add_wrapper(ret_matrix, &my_matrix, is_root, NULL, 0, 0, 0, bitmaps[0], fail);

		if (*fail) {
			return NULL;
		}
	}

	if (!is_root && !my_matrix) {
		my_matrix = supernodal_frontal_create(v, sn_size, snL->sn_up_size[sn], snL->sn_struct[sn]);
		if (!my_matrix) {
			*fail = TRUE;
			return NULL;
		}
	}

	if (!is_root) {
		int rc;
		rc = multifrontal_supernodal_front_factor(sn, v, sn_size, A, my_matrix, bitmaps[0], snL);
		if (rc) {
			*fail = TRUE;
			supernodal_frontal_free(my_matrix);
			return NULL;
		}
	}
	return my_matrix;
}

void *taucs_dtl(ccs_factor_llt_mf) (taucs_ccs_matrix * A) {
	void *p;
	p = taucs_dtl(ccs_factor_llt_mf_maxdepth) (A, 0);
	return p;
}

static void recursive_multifrontal_supernodal_factor_llt_caller(int n_sn,	/* this supernode */
								int UNUSED(is_root),	/* is v the root? */
								taucs_ccs_matrix *A, supernodal_factor_matrix *snL, int *fail)
{
	int **maps;
	int i, j;

	maps = (int **) taucs_malloc(sizeof(int *));
	if (!maps) {
		taucs_supernodal_factor_free(snL);
		assert(0);
		return;
	}

	for (i = 0; i < 1; i++) {
		maps[i] = (int *) taucs_malloc((A->n + 1) * sizeof(int));
		if (!maps[i]) {
			for (j = 0; j < i; j++)
				taucs_free(maps[j]);
			taucs_free(maps);
			taucs_supernodal_factor_free(snL);
			assert(0);
			return;
		}
	}

	recursive_multifrontal_supernodal_factor_llt(n_sn, TRUE, maps, A, snL, fail);
	for (i = 0; i < 1; i++)
		taucs_free(maps[i]);
	taucs_free(maps);

}

void *taucs_dtl(ccs_factor_llt_mf_maxdepth) (taucs_ccs_matrix * A, int max_depth) {
	supernodal_factor_matrix * L;

	int fail;
	L = multifrontal_supernodal_create();
	if (!L)
		return NULL;

	fail = taucs_ccs_symbolic_elimination(A, L, FALSE /* don't sort row indices */ ,
					      max_depth);
	if (fail == -1) {
		taucs_supernodal_factor_free(L);
		return NULL;
	}

	fail = FALSE;
	recursive_multifrontal_supernodal_factor_llt_caller((L->n_sn), TRUE, A, L, &fail);
	if (fail) {
		taucs_supernodal_factor_free(L);
		return NULL;
	}

	return (void *) L;
}


void *taucs_dtl(ccs_factor_llt_symbolic) (taucs_ccs_matrix * A) {
	return taucs_dtl(ccs_factor_llt_symbolic_maxdepth) (A, 0);
}

void *taucs_dtl(ccs_factor_llt_symbolic_maxdepth) (taucs_ccs_matrix * A, int max_depth) {
	supernodal_factor_matrix * L;
	int fail;

	L = multifrontal_supernodal_create();
	if (!L)
		return NULL;

	fail = taucs_ccs_symbolic_elimination(A, L, FALSE /* don't sort row indices */ ,
					      max_depth);
	if (fail == -1) {
		taucs_supernodal_factor_free(L);
		return NULL;
	}
	return L;
}

int taucs_dtl(ccs_factor_llt_numeric) (taucs_ccs_matrix * A, void *vL) {
	supernodal_factor_matrix * L = (supernodal_factor_matrix *) vL;
	int fail;
	fail = FALSE;
	recursive_multifrontal_supernodal_factor_llt_caller((L->n_sn), TRUE, A, L, &fail);
	if (fail) {
		taucs_supernodal_factor_free_numeric(L);
		return -1;
	}
	return 0;
}


static void
recursive_leftlooking_supernodal_update(int J, int K, int bitmap[], double *dense_update_matrix, taucs_ccs_matrix *A, supernodal_factor_matrix *L)
{
	int i, j, ir;
	int child;
	int *first_child = L->first_child;
	int *next_child = L->next_child;
	int sn_size_father = (L->sn_size)[J];
	int sn_up_size_father = (L->sn_up_size)[J];
	int sn_size_child = (L->sn_size)[K];
	int sn_up_size_child = (L->sn_up_size)[K];
	int exist_upd = 0;
	int first_row = 0;
	int row_count = 0;
	int PK, M, N, LDA, LDB, LDC;

	for (i = sn_size_child; i < sn_up_size_child; i++)
		if (bitmap[L->sn_struct[K][i]]
		    && L->sn_struct[K][i] <= L->sn_struct[J][sn_size_father - 1]) {
			if (!exist_upd)
				first_row = i;
			row_count++;
			exist_upd = 1;
		}

	if (exist_upd) {
		LDA = LDB = (L->up_blocks_ld)[K];
		M = sn_up_size_child - first_row;	       /* +-1 ? */
		LDC = sn_up_size_father;
		N = row_count;
		PK = L->sn_size[K];

		taucs_herk("L",
			   "N",
			   &N, &PK,
			   &taucs_one_real_const,
			   &(L->up_blocks[K][first_row - sn_size_child]), &LDA, &taucs_zero_real_const, dense_update_matrix, &LDC, F_ONE, F_ONE);

		if (M - N > 0) {
			int newM = M - N;

			taucs_gemm("N",
				   "C",
				   &newM, &N, &PK,
				   &taucs_one_const,
				   &(L->up_blocks[K][first_row - sn_size_child + N]), &LDA,
				   &(L->up_blocks[K][first_row - sn_size_child]), &LDB,
				   &taucs_zero_const, dense_update_matrix + N, &LDC, F_ONE, F_ONE);
		}

		for (j = 0; j < row_count; j++)
			for (ir = j; ir < row_count; ir++) {

				L->sn_blocks[J][(bitmap[L->sn_struct[K][first_row + j]] - 1) * sn_size_father +
						(bitmap[L->sn_struct[K][first_row + ir]] - 1)] =
				    taucs_sub(L->sn_blocks[J][(bitmap[L->sn_struct[K][first_row + j]] - 1) * sn_size_father +
							      (bitmap[L->sn_struct[K][first_row + ir]] - 1)], dense_update_matrix[j * LDC + ir]);
			}

		for (j = 0; j < row_count; j++)
			for (ir = row_count; ir < M; ir++) {
				L->up_blocks[J][(bitmap[L->sn_struct[K][first_row + j]] - 1) * (L->up_blocks_ld)[J] +
						(bitmap[L->sn_struct[K][ir + first_row]] - 1)] =
				    taucs_sub(L->up_blocks[J][(bitmap[L->sn_struct[K][first_row + j]] - 1) * (L->up_blocks_ld)[J] +
							      (bitmap[L->sn_struct[K][ir + first_row]] - 1)], dense_update_matrix[j * LDC + ir]);

			}
		for (child = first_child[K]; child != -1; child = next_child[child]) {
			recursive_leftlooking_supernodal_update(J, child, bitmap, dense_update_matrix, A, L);
		}
	}
}

static int leftlooking_supernodal_front_factor(int sn, int *bitmap, taucs_ccs_matrix *A, supernodal_factor_matrix *L)
{
	int ip, jp;
	int *ind;
	double *re;
	int INFO;

	int sn_size = (L->sn_size)[sn];
	int up_size = (L->sn_up_size)[sn] - (L->sn_size)[sn];

	for (ip = 0; ip < (L->sn_up_size)[sn]; ip++)
		bitmap[(L->sn_struct)[sn][ip]] = ip;

	for (jp = 0; jp < sn_size; jp++) {
		ind = &(A->rowind[A->colptr[(L->sn_struct)[sn][jp]]]);
		re = &(A->values[A->colptr[(L->sn_struct)[sn][jp]]]);
		for (ip = 0; ip < A->colptr[(L->sn_struct)[sn][jp] + 1]
		     - A->colptr[(L->sn_struct)[sn][jp]]; ip++) {
			if (bitmap[ind[ip]] < sn_size)
				(L->sn_blocks)[sn][(L->sn_blocks_ld)[sn] * jp + bitmap[ind[ip]]] =
				    taucs_add((L->sn_blocks)[sn][(L->sn_blocks_ld)[sn] * jp + bitmap[ind[ip]]], re[ip]);
			else
				(L->up_blocks)[sn][(L->up_blocks_ld)[sn] * jp + bitmap[ind[ip]] - sn_size] =
				    taucs_add((L->up_blocks)[sn][(L->up_blocks_ld)[sn] * jp + bitmap[ind[ip]] - sn_size], re[ip]);
		}
	}
	if (sn_size) {
		taucs_potrf("L", &sn_size, (L->sn_blocks)[sn], &((L->sn_blocks_ld)[sn]), &INFO, F_ONE);
	}

	if (INFO) {
		taucs_printf("\t\tLL^T Factorization: Matrix is not positive definite.\n");
		taucs_printf("\t\t                    nonpositive pivot in column %d\n", (L->sn_struct)[INFO - 1]);
		return -1;
	}

	if (up_size && sn_size)
		taucs_trsm("R",
			   "L",
			   "C",
			   "N",
			   &up_size, &sn_size,
			   &taucs_one_const,
			   (L->sn_blocks)[sn], &((L->sn_blocks_ld)[sn]), (L->up_blocks)[sn], &((L->up_blocks_ld)[sn]), F_ONE, F_ONE, F_ONE, F_ONE);

	return 0;
}

static int recursive_leftlooking_supernodal_factor_llt(int sn, /* this supernode */
						       int is_root,	/* is v the root? */
						       int *bitmap, int *indmap, taucs_ccs_matrix *A, supernodal_factor_matrix *L)
{
	int child;
	int *first_child = L->first_child;
	int *next_child = L->next_child;
	double *dense_update_matrix = NULL;

	if (!is_root) {
		(L->sn_blocks)[sn] = (L->up_blocks)[sn] = NULL;
		if (L->sn_size[sn]) {
			(L->sn_blocks)[sn] = (double *) taucs_calloc(((L->sn_size)[sn]) * ((L->sn_size)[sn]), sizeof(double));
			if (!((L->sn_blocks)[sn]))
				return -1;		       /* the caller will free L */
		}
		(L->sn_blocks_ld)[sn] = (L->sn_size)[sn];

		if (((L->sn_up_size)[sn] - (L->sn_size)[sn]) && (L->sn_size)[sn]) {
			(L->up_blocks)[sn] = (double *) taucs_calloc(((L->sn_up_size)[sn] - (L->sn_size)[sn])
								     * ((L->sn_size)[sn]), sizeof(double));
			if (!((L->up_blocks)[sn]))
				return -1;		       /* the caller will free L */
		}
		(L->up_blocks_ld)[sn] = (L->sn_up_size)[sn] - (L->sn_size)[sn];
	}

	for (child = first_child[sn]; child != -1; child = next_child[child]) {
		if (recursive_leftlooking_supernodal_factor_llt(child, FALSE, bitmap, indmap, A, L)
		    == -1) {
			taucs_free(dense_update_matrix);
			return -1;
		}

		if (!is_root) {
			if (!dense_update_matrix) {
				dense_update_matrix = (double *) taucs_calloc((L->sn_up_size)[sn] * (L->sn_size)[sn], sizeof(double));
				if (!dense_update_matrix)
					return -1;	       /* caller will free L */
			}

			/*
			 * prepare the bitmap. Moved out of the recusive update procedure 20/1/2003. Sivan and Elad 
			 */

			{
				int i;
				int J = sn;
				int sn_size_father = (L->sn_size)[J];
				int sn_up_size_father = (L->sn_up_size)[J];

				for (i = 0; i < sn_size_father; i++)
					bitmap[L->sn_struct[J][i]] = i + 1;
				for (i = sn_size_father; i < sn_up_size_father; i++)
					bitmap[L->sn_struct[J][i]] = i - sn_size_father + 1;
			}

			recursive_leftlooking_supernodal_update(sn, child, bitmap, dense_update_matrix, A, L);

			{
				int i;
				int J = sn;
				int sn_size_father = (L->sn_size)[J];
				int sn_up_size_father = (L->sn_up_size)[J];

				for (i = 0; i < sn_size_father; i++)
					bitmap[L->sn_struct[J][i]] = 0;
				for (i = 0; i < sn_up_size_father; i++)
					bitmap[L->sn_struct[J][i]] = 0;
			}

		}
	}
	taucs_free(dense_update_matrix);

	if (!is_root) {
		if (leftlooking_supernodal_front_factor(sn, indmap, A, L)) {
			return -1;			       /* nonpositive pivot */
		}
	}

	return 0;
}

void *taucs_dtl(ccs_factor_llt_ll) (taucs_ccs_matrix * A) {
	return taucs_dtl(ccs_factor_llt_ll_maxdepth) (A, 0);
}

void *taucs_dtl(ccs_factor_llt_ll_maxdepth) (taucs_ccs_matrix * A, int max_depth) {
	supernodal_factor_matrix * L;
	int *map;
	int *map2;
	int fail;

	L = multifrontal_supernodal_create();
	if (!L)
		return NULL;

	fail = taucs_ccs_symbolic_elimination(A, L, TRUE /* sort row indices */ ,
					      max_depth);

	map = (int *) taucs_malloc((A->n + 1) * sizeof(int));
	map2 = (int *) taucs_calloc((A->n + 1), sizeof(int));

	if (fail == -1 || !map || !map2) {
		taucs_supernodal_factor_free(L);
		taucs_free(map2);
		taucs_free(map);
		return NULL;
	}

	if (recursive_leftlooking_supernodal_factor_llt((L->n_sn), TRUE, map2, map, A, L)
	    == -1) {
		taucs_supernodal_factor_free(L);
		taucs_free(map);
		taucs_free(map2);
		return NULL;
	}

	taucs_free(map);
	taucs_free(map2);

	return (void *) L;
}


static void recursive_supernodal_solve_l(int sn,	       /* this supernode */
					 int is_root,	       /* is v the root? */
					 int *first_child, int *next_child,
					 int **sn_struct, int *sn_sizes, int *sn_up_sizes,
					 int *sn_blocks_ld, double *sn_blocks[],
					 int *up_blocks_ld, double *up_blocks[], double x[], double b[], double t[])
{
	int child;
	int sn_size;					       /* number of rows/columns in the supernode */
	int up_size;					       /* number of rows that this supernode updates */
	int ione = 1;

	double *xdense;
	double *bdense;
	int i;						       /* ip,j,jp omer */

	for (child = first_child[sn]; child != -1; child = next_child[child]) {
		recursive_supernodal_solve_l(child,
					     FALSE,
					     first_child, next_child,
					     sn_struct, sn_sizes, sn_up_sizes, sn_blocks_ld, sn_blocks, up_blocks_ld, up_blocks, x, b, t);
	}

	if (!is_root) {

		sn_size = sn_sizes[sn];
		up_size = sn_up_sizes[sn] - sn_sizes[sn];

		xdense = t;
		bdense = t + sn_size;

		for (i = 0; i < sn_size; i++)
			xdense[i] = b[sn_struct[sn][i]];
		for (i = 0; i < up_size; i++)
			bdense[i] = taucs_zero;

		taucs_trsm("L",
			   "L",
			   "N",
			   "N",
			   &sn_size, &ione, &taucs_one_const, sn_blocks[sn], &(sn_blocks_ld[sn]), xdense, &sn_size, F_ONE, F_ONE, F_ONE, F_ONE);

		if (up_size > 0 && sn_size > 0) {
			taucs_gemm("N", "N",
				   &up_size, &ione, &sn_size,
				   &taucs_one_const,
				   up_blocks[sn], &(up_blocks_ld[sn]), xdense, &sn_size, &taucs_zero_const, bdense, &up_size, F_ONE, F_ONE);
		}

		for (i = 0; i < sn_size; i++)
			x[sn_struct[sn][i]] = xdense[i];
		for (i = 0; i < up_size; i++)
			/*
			 * b[ sn_struct[ sn ][ sn_size + i ] ] -= bdense[i];
			 */
			b[sn_struct[sn][sn_size + i]] = taucs_sub(b[sn_struct[sn][sn_size + i]], bdense[i]);

	}
}

static void recursive_supernodal_solve_lt(int sn,	       /* this supernode */
					  int is_root,	       /* is v the root? */
					  int *first_child, int *next_child,
					  int **sn_struct, int *sn_sizes, int *sn_up_sizes,
					  int *sn_blocks_ld, double *sn_blocks[],
					  int *up_blocks_ld, double *up_blocks[], double x[], double b[], double t[])
{
	int child;
	int sn_size;					       /* number of rows/columns in the supernode */
	int up_size;					       /* number of rows that this supernode updates */
	int ione = 1;

	double *xdense;
	double *bdense;
	int i;						       /* ip,j,jp omer */

	if (!is_root) {

		sn_size = sn_sizes[sn];
		up_size = sn_up_sizes[sn] - sn_sizes[sn];

		bdense = t;
		xdense = t + sn_size;

		for (i = 0; i < sn_size; i++)
			bdense[i] = b[sn_struct[sn][i]];
		for (i = 0; i < up_size; i++)
			xdense[i] = x[sn_struct[sn][sn_size + i]];

		if (up_size > 0 && sn_size > 0)
			taucs_gemm("C", "N",
				   &sn_size, &ione, &up_size,
				   &taucs_minusone_const,
				   up_blocks[sn], &(up_blocks_ld[sn]), xdense, &up_size, &taucs_one_const, bdense, &sn_size, F_ONE, F_ONE);

		taucs_trsm("L",
			   "L",
			   "C",
			   "N",
			   &sn_size, &ione, &taucs_one_const, sn_blocks[sn], &(sn_blocks_ld[sn]), bdense, &sn_size, F_ONE, F_ONE, F_ONE, F_ONE);

		for (i = 0; i < sn_size; i++)
			x[sn_struct[sn][i]] = bdense[i];
	}

	for (child = first_child[sn]; child != -1; child = next_child[child]) {
		recursive_supernodal_solve_lt(child,
					      FALSE,
					      first_child, next_child,
					      sn_struct, sn_sizes, sn_up_sizes, sn_blocks_ld, sn_blocks, up_blocks_ld, up_blocks, x, b, t);
	}
}

int taucs_dtl(supernodal_solve_llt) (void *vL, void *vx, void *vb) {
	supernodal_factor_matrix * L = (supernodal_factor_matrix *) vL;
	double *x = (double *) vx;
	double *b = (double *) vb;
	double *y;
	double *t;					       /* temporary vector */
	int i;

	y = (double *) taucs_malloc((L->n) * sizeof(double));
	t = (double *) taucs_malloc((L->n) * sizeof(double));
	if (!y || !t) {
		taucs_free(y);
		taucs_free(t);
		taucs_printf("multifrontal_supernodal_solve_llt: out of memory\n");
		return -1;
	}

	for (i = 0; i < L->n; i++)
		x[i] = b[i];

	recursive_supernodal_solve_l(L->n_sn, TRUE,	       /* this is the root */
				     L->first_child, L->next_child,
				     L->sn_struct, L->sn_size, L->sn_up_size,
				     L->sn_blocks_ld, L->sn_blocks, L->up_blocks_ld, L->up_blocks, y, x, t);

	recursive_supernodal_solve_lt(L->n_sn, TRUE,	       /* this is the root */
				      L->first_child, L->next_child,
				      L->sn_struct, L->sn_size, L->sn_up_size,
				      L->sn_blocks_ld, L->sn_blocks, L->up_blocks_ld, L->up_blocks, x, y, t);

	taucs_free(y);
	taucs_free(t);

	return 0;
}
#endif							       /* #ifndef TAUCS_CORE_GENERAL */

#ifdef TAUCS_CORE_GENERAL

void *taucs_ccs_factor_llt_mf_maxdepth(taucs_ccs_matrix *A, int max_depth)
{
	void *p = NULL;

	if (A->flags & TAUCS_DOUBLE)
		p = taucs_dccs_factor_llt_mf_maxdepth(A, max_depth);
	return p;
}

void *taucs_ccs_factor_llt_ll_maxdepth(taucs_ccs_matrix *A, int max_depth)
{

	if (A->flags & TAUCS_DOUBLE)
		return taucs_dccs_factor_llt_ll_maxdepth(A, max_depth);
	assert(0);
	return NULL;
}

void *taucs_ccs_factor_llt_symbolic_maxdepth(taucs_ccs_matrix *A, int max_depth)
{
	if (A->flags & TAUCS_DOUBLE)
		return taucs_dccs_factor_llt_symbolic_maxdepth(A, max_depth);
	assert(0);
	return NULL;
}

void *taucs_ccs_factor_llt_mf(taucs_ccs_matrix *A)
{
	void *p = NULL;
	if (A->flags & TAUCS_DOUBLE)
		p = taucs_dccs_factor_llt_mf(A);
	return p;
}

void *taucs_ccs_factor_llt_ll(taucs_ccs_matrix *A)
{
	if (A->flags & TAUCS_DOUBLE)
		return taucs_dccs_factor_llt_ll(A);
	assert(0);
	return NULL;
}

void *taucs_ccs_factor_llt_symbolic(taucs_ccs_matrix *A)
{
	if (A->flags & TAUCS_DOUBLE)
		return taucs_dccs_factor_llt_symbolic(A);
	assert(0);
	return NULL;
}

int taucs_ccs_factor_llt_numeric(taucs_ccs_matrix *A, void *L)
{
	int rc = TAUCS_ERROR;
	if (A->flags & TAUCS_DOUBLE)
		rc = taucs_dccs_factor_llt_numeric(A, L);
	return rc;
}

int taucs_supernodal_solve_llt(void *L, void *x, void *b)
{
	if (((supernodal_factor_matrix *) L)->flags & TAUCS_DOUBLE)
		return taucs_dsupernodal_solve_llt(L, x, b);
	assert(0);
	return -1;
}

void taucs_supernodal_factor_free(void *L)
{
	if (((supernodal_factor_matrix *) L)->flags & TAUCS_DOUBLE) {
		taucs_dsupernodal_factor_free(L);
		return;
	}
}

void taucs_supernodal_factor_free_numeric(void *L)
{
	if (((supernodal_factor_matrix *) L)->flags & TAUCS_DOUBLE) {
		taucs_dsupernodal_factor_free_numeric(L);
	}
}

taucs_ccs_matrix *taucs_supernodal_factor_to_ccs(void *L)
{
	if (((supernodal_factor_matrix *) L)->flags & TAUCS_DOUBLE)
		return taucs_dsupernodal_factor_to_ccs(L);
	return NULL;
}

void *taucs_supernodal_factor_get_diag(void *L)
{
	if (((supernodal_factor_matrix *) L)->flags & TAUCS_DOUBLE)
		return taucs_dsupernodal_factor_get_diag(L);
	return NULL;
}

int taucs_ccs_etree(taucs_ccs_matrix *A, int *parent, int *l_colcount, int *l_rowcount, int *l_nnz)
{
	int *prev_p;

	/*
	 * int* prev_nbr;omer
	 */
	int *level;

	/*
	 * int* first_descendant;omer
	 */
	int *l_cc;
	int *l_rc;
	int *wt;

	int n = A->n;
	int pprime;					       /* p,q,u omer */
	int ju;
	int *postorder;
	int *ipostorder;
	int *first_child, *next_child;

	int i, j, k, ip, jp, kp;
	int nnz, jnnz;
	int *uf;
	int *rowptr;
	int *colind;
	int *rowcount;
	int *realroot;

	/*
	 * we need the row structures for the lower triangle 
	 */

	nnz = (A->colptr)[n];

	uf = (int *) taucs_malloc(n * sizeof(int));
	rowcount = (int *) taucs_malloc((n + 1) * sizeof(int));
	rowptr = (int *) taucs_malloc((n + 1) * sizeof(int));
	colind = (int *) taucs_malloc(nnz * sizeof(int));

	if (!uf || !rowcount || !rowptr || !colind) {
		taucs_free(uf);
		taucs_free(rowcount);
		taucs_free(rowptr);
		taucs_free(colind);
		return -1;
	}

	for (i = 0; i <= n; i++)
		rowcount[i] = 0;
	for (j = 0; j < n; j++) {
		jnnz = (A->colptr)[j + 1] - (A->colptr)[j];
		for (ip = 0; ip < jnnz; ip++) {
			i = (A->rowind)[(A->colptr)[j] + ip];
			if (j < i)
				rowcount[i]++;
		}
	}

	ip = 0;
	for (i = 0; i <= n; i++) {
		int next_ip = ip + rowcount[i];

		rowcount[i] = ip;
		rowptr[i] = ip;
		ip = next_ip;
	}

	for (j = 0; j < n; j++) {
		jnnz = (A->colptr)[j + 1] - (A->colptr)[j];
		for (ip = 0; ip < jnnz; ip++) {
			i = (A->rowind)[(A->colptr)[j] + ip];
			if (i == j)
				continue;
			assert(rowcount[i] < rowptr[i + 1]);
			colind[rowcount[i]] = j;
			rowcount[i]++;
		}
	}

	/*
	 * now compute the etree 
	 */

	{
		int u, t, vroot;

		realroot = rowcount;			       /* reuse space */

		for (i = 0; i < n; i++) {
			uf_makeset(uf, i);
			realroot[i] = i;
			parent[i] = n;
			vroot = i;
			for (kp = rowptr[i]; kp < rowptr[i + 1]; kp++) {
				k = colind[kp];
				u = uf_find(uf, k);
				t = realroot[u];
				if (parent[t] == n && t != i) {
					parent[t] = i;
					vroot = uf_union(uf, vroot, u);
					realroot[vroot] = i;
				}
			}
		}
	}

	taucs_free(colind);
	taucs_free(rowptr);
	taucs_free(rowcount);

	/*
	 * now only uf remains allocated 
	 */

	/*
	 * compute column counts 
	 */

	if (l_colcount || l_rowcount || l_nnz) {
		int *l_nz;
		int tmp;
		int u, p, q;

		first_child = (int *) taucs_malloc((n + 1) * sizeof(int));
		next_child = (int *) taucs_malloc((n + 1) * sizeof(int));
		postorder = (int *) taucs_malloc(n * sizeof(int));
		ipostorder = (int *) taucs_malloc(n * sizeof(int));
		wt = (int *) taucs_malloc(n * sizeof(int));
		level = (int *) taucs_malloc(n * sizeof(int));
		prev_p = (int *) taucs_malloc(n * sizeof(int));

		/*
		 * we allocate scratch vectors to avoid conditionals 
		 */
		/*
		 * in the inner loop.  
		 */

		if (l_colcount)
			l_cc = l_colcount;
		else
			l_cc = (int *) (int *) taucs_malloc(n * sizeof(int));
		if (l_rowcount)
			l_rc = l_rowcount;
		else
			l_rc = (int *) (int *) taucs_malloc(n * sizeof(int));
		if (l_nnz)
			l_nz = l_nnz;
		else
			l_nz = &tmp;

		if (!first_child || !next_child || !postorder || !ipostorder || !wt || !level || !prev_p || (!l_colcount && !l_cc)
		    || (!l_rowcount && !l_rc)
		    ) {
			taucs_free(uf);

			if (!l_colcount)
				taucs_free(l_cc);
			if (!l_rowcount)
				taucs_free(l_rc);

			taucs_free(postorder);
			taucs_free(ipostorder);
			taucs_free(wt);
			taucs_free(level);
			taucs_free(prev_p);
			return -1;
		}

		/*
		 * for (j=0; j<n; j++) printf("parent[%d] = %d\n",j,parent[j]);
		 */

		/*
		 * compute the postorder 
		 */

		for (j = 0; j <= n; j++)
			first_child[j] = -1;
		for (j = n - 1; j >= 0; j--) {
			next_child[j] = first_child[parent[j]];
			first_child[parent[j]] = j;
		}

		{
			int next = 0;

			recursive_postorder(n, first_child, next_child, postorder, ipostorder, &next);
		}

		/*
		 * sort by postorder of etree 
		 */
		/*
		 * compute level, fst_desc 
		 */

		tree_level(n, TRUE, first_child, next_child, level, -1);

		for (u = 0; u < n; u++)
			prev_p[u] = -1;
		for (u = 0; u < n; u++)
			l_rc[u] = 1;
		for (u = 0; u < n; u++)
			ordered_uf_makeset(uf, u);
		for (u = 0; u < n; u++) {
			if (first_child[u] == -1)
				wt[u] = 1;		       /* leaves */
			else
				wt[u] = 0;		       /* nonleaves */
		}
		taucs_free(first_child);
		taucs_free(next_child);

		for (p = 0; p < n; p++) {
			jp = postorder[p];
			if (parent[jp] != n)
				wt[parent[jp]]--;
			for (ip = (A->colptr)[jp]; ip < (A->colptr)[jp + 1]; ip++) {
				ju = (A->rowind)[ip];
				u = ipostorder[ju];
				if (ju == jp)
					continue;	       /* we only want proper neighbors */
				if (1) {
					wt[jp]++;
					pprime = prev_p[ju];
					if (pprime == -1)
						l_rc[ju] += level[jp] - level[ju];
					else {
						q = ordered_uf_find(uf, pprime);
						l_rc[ju] += level[jp] - level[q];
						wt[q]--;
					}
					prev_p[ju] = jp;
				}
			}
			if (parent[jp] != n) {
				if (!(ipostorder[parent[jp]] > ipostorder[jp])) {
					printf("jp %d parent %d (ipo_j %d ipo_parent %d\n", jp, parent[jp], ipostorder[jp], ipostorder[parent[jp]]);
				}
				assert(ipostorder[parent[jp]] > ipostorder[jp]);
				ordered_uf_union(uf, jp, parent[jp]);
			}
		}

		*l_nz = 0;
		for (u = 0; u < n; u++) {
			l_cc[u] = wt[u];
			*l_nz += wt[u];
		}
		for (u = 0; u < n; u++) {
			if (parent[u] != n) {
				l_cc[parent[u]] += l_cc[u];
				*l_nz += l_cc[u];
			}
		}

		/*
		 * free scrtach vectors 
		 */

		if (!l_colcount)
			taucs_free(l_cc);
		if (!l_rowcount)
			taucs_free(l_rc);

		/*
		 * free other data structures 
		 */

		taucs_free(postorder);
		taucs_free(ipostorder);
		taucs_free(wt);
		taucs_free(level);
		taucs_free(prev_p);
	}
	taucs_free(uf);
	return 0;
}

int taucs_ccs_symbolic_elimination(taucs_ccs_matrix *A, void *vL, int do_order, int max_depth)
{
	supernodal_factor_matrix *L = (supernodal_factor_matrix *) vL;
	int *first_child;
	int *next_child;
	int j;
	int *column_to_sn_map;
	int *map;
	int *rowind;
	int *parent;
	int *ipostorder;

	int depth;

	L->n = A->n;
	/*
	 * use calloc so we can deallocate unallocated entries 
	 */
	L->sn_struct = (int **) taucs_calloc((A->n), sizeof(int *));
	L->sn_size = (int *) taucs_calloc((A->n + 1), sizeof(int));
	L->sn_up_size = (int *) taucs_calloc((A->n + 1), sizeof(int));
	L->first_child = (int *) taucs_calloc((A->n + 1), sizeof(int));
	L->next_child = (int *) taucs_calloc((A->n + 1), sizeof(int));

	column_to_sn_map = (int *) taucs_calloc((A->n + 1), sizeof(int));
	map = (int *) taucs_calloc((A->n + 1), sizeof(int));
	first_child = (int *) taucs_calloc((A->n + 1), sizeof(int));
	next_child = (int *) taucs_calloc((A->n + 1), sizeof(int));
	parent = (int *) taucs_calloc((A->n + 1), sizeof(int));
	rowind = (int *) taucs_calloc((A->n), sizeof(int));

	if (!(L->sn_struct) || !(L->sn_size) || !(L->sn_up_size) ||
	    !(L->first_child) || !(L->next_child) || !column_to_sn_map || !map || !first_child || !next_child || !rowind || !parent) {
		taucs_free(parent);
		taucs_free(rowind);
		taucs_free(next_child);
		taucs_free(first_child);
		taucs_free(map);
		taucs_free(column_to_sn_map);
		taucs_free(L->next_child);
		taucs_free(L->first_child);
		taucs_free(L->sn_up_size);
		taucs_free(L->sn_size);
		taucs_free(L->sn_struct);
		L->sn_struct = NULL;
		L->sn_size = L->sn_up_size = L->first_child = L->next_child = NULL;

		return -1;
	}

	if (taucs_ccs_etree(A, parent, NULL, NULL, NULL) == -1) {
		taucs_free(parent);
		taucs_free(rowind);
		taucs_free(next_child);
		taucs_free(first_child);
		taucs_free(map);
		taucs_free(column_to_sn_map);
		taucs_free(L->next_child);
		taucs_free(L->first_child);
		taucs_free(L->sn_up_size);
		taucs_free(L->sn_size);
		taucs_free(L->sn_struct);
		L->sn_struct = NULL;
		L->sn_size = L->sn_up_size = L->first_child = L->next_child = NULL;
		return -1;
	}

	for (j = 0; j <= (A->n); j++)
		first_child[j] = -1;
	for (j = (A->n) - 1; j >= 0; j--) {
		int p = parent[j];
		next_child[j] = first_child[p];
		first_child[p] = j;
	}

	/*
	 * let's compute the depth of the etree, to bail out if it is too deep 
	 */
	/*
	 * the whole thing will work better if we compute supernodal etrees 
	 */

	{
		int next_depth_count;
		int this_depth_count;
		int child, i;

		int *this_depth = rowind;		       /* we alias rowind */
		int *next_depth = map;			       /* and map */
		int *tmp;

		this_depth[0] = A->n;
		this_depth_count = 1;
		next_depth_count = 0;
		depth = -1;

		while (this_depth_count) {
			for (i = 0; i < this_depth_count; i++) {
				child = first_child[this_depth[i]];
				while (child != -1) {
					next_depth[next_depth_count] = child;
					next_depth_count++;
					child = next_child[child];
				}
			}

			tmp = this_depth;
			this_depth = next_depth;
			next_depth = tmp;

			this_depth_count = next_depth_count;
			next_depth_count = 0;
			depth++;
		}
	}

	taucs_printf("\t\tElimination tree depth is %d\n", depth);

	if (max_depth && depth > max_depth) {
		taucs_printf("taucs_ccs_symbolic_elimination: etree depth %d, maximum allowed is %d\n", depth, max_depth);
		taucs_free(parent);
		taucs_free(rowind);
		taucs_free(next_child);
		taucs_free(first_child);
		taucs_free(map);
		taucs_free(column_to_sn_map);
		taucs_free(L->next_child);
		taucs_free(L->first_child);
		taucs_free(L->sn_up_size);
		taucs_free(L->sn_size);
		taucs_free(L->sn_struct);
		L->sn_struct = NULL;
		L->sn_size = L->sn_up_size = L->first_child = L->next_child = NULL;
		return -1;
	}

	/*
	 * taucs_free(parent); ipostorder = (int*)taucs_calloc((A->n+1), sizeof(int)); 
	 */

	ipostorder = parent;
	{
		int next = 0;

		/*
		 * int* postorder = (int*)taucs_calloc((A->n+1), sizeof(int));
		 */
		recursive_postorder(A->n, first_child, next_child, NULL, ipostorder, &next);
		/*
		 * printf("ipostorder "); for (j=0; j <= (A->n); j++) printf("%d ",ipostorder[j]); printf("\n"); printf(" postorder 
		 * "); for (j=0; j <= (A->n); j++) printf("%d ",postorder[j]); printf("\n"); 
		 */
	}

	L->n_sn = 0;
	for (j = 0; j < (A->n); j++)
		map[j] = -1;
	for (j = 0; j <= (A->n); j++)
		(L->first_child)[j] = (L->next_child)[j] = -1;

	if (recursive_symbolic_elimination(A->n,
					   A,
					   first_child, next_child,
					   &(L->n_sn),
					   L->sn_size, L->sn_up_size, L->sn_struct,
					   L->first_child, L->next_child, rowind, column_to_sn_map, map, do_order, ipostorder)
	    == -1) {
		for (j = 0; j < (A->n); j++)
			taucs_free((L->sn_struct)[j]);

		taucs_free(parent);
		taucs_free(rowind);
		taucs_free(next_child);
		taucs_free(first_child);
		taucs_free(map);
		taucs_free(column_to_sn_map);
		taucs_free(L->next_child);
		taucs_free(L->first_child);
		taucs_free(L->sn_up_size);
		taucs_free(L->sn_size);
		taucs_free(L->sn_struct);
		L->sn_struct = NULL;
		L->sn_size = L->sn_up_size = L->first_child = L->next_child = NULL;
		return -1;
	}

	for (j = 0; j < (A->n); j++)
		map[j] = -1;

	taucs_free(parent);
	taucs_free(rowind);
	taucs_free(map);
	taucs_free(column_to_sn_map);
	taucs_free(next_child);
	taucs_free(first_child);

	L->sn_blocks_ld = (int *) taucs_calloc((L->n_sn), sizeof(int));
	L->sn_blocks = (double **) taucs_calloc((L->n_sn), sizeof(double *));	/* so we can free before allocation */

	L->up_blocks_ld = (int *) taucs_calloc((L->n_sn), sizeof(int));
	L->up_blocks = (double **) taucs_calloc((L->n_sn), sizeof(double *));

	if (!(L->sn_blocks_ld)
	    || !(L->sn_blocks_ld)
	    || !(L->sn_blocks)
	    || !(L->up_blocks_ld)
	    || !(L->up_blocks))
		return -1;				       /* the caller will free L */

	return 0;
}
#endif
