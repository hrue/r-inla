#include "GMRFLib/GMRFLib.h"
#include "GMRFLib/GMRFLibP.h"

int GMRFLib_compute_reordering(GMRFLib_sm_fact_tp *sm_fact, GMRFLib_graph_tp *graph, GMRFLib_global_node_tp *gn)
{
	GMRFLib_ENTER_FUNCTION;
	GMRFLib_global_node_tp lgn, *gn_ptr = NULL;

	if (gn) {
		/*
		 * then this defines the global node definition 
		 */
		gn_ptr = (GMRFLib_global_node_tp *) gn;
	} else {
		lgn = GMRFLib_global_node;
		gn_ptr = &lgn;
	}

	GMRFLib_reorder_tp r = GMRFLib_reorder;
	if (sm_fact->smtp == GMRFLib_SMTP_PARDISO || sm_fact->smtp == GMRFLib_SMTP_STILES) {
		r = GMRFLib_REORDER_DEFAULT;
	} else if ((sm_fact->smtp == GMRFLib_SMTP_TAUCS || sm_fact->smtp == GMRFLib_SMTP_BAND) &&
		   (r == GMRFLib_REORDER_STILES || r == GMRFLib_REORDER_PARDISO)) {
		r = GMRFLib_REORDER_DEFAULT;
	}

	switch (r) {
	case GMRFLib_REORDER_DEFAULT:
	{
		/*
		 * this choice depends on the sparse-solver 
		 */
		switch (sm_fact->smtp) {
		case GMRFLib_SMTP_BAND:
		{
			GMRFLib_EWRAP1(GMRFLib_compute_reordering_BAND(&(sm_fact->remap), graph));
		}
			break;

		case GMRFLib_SMTP_TAUCS:
		{
			GMRFLib_EWRAP1(GMRFLib_compute_reordering_TAUCS(&(sm_fact->remap), graph, r, gn_ptr));
		}
			break;

		case GMRFLib_SMTP_PARDISO:
		{
			if (sm_fact->PARDISO_fact == NULL) {
				GMRFLib_pardiso_init(&(sm_fact->PARDISO_fact));
			}
			GMRFLib_pardiso_reorder(sm_fact->PARDISO_fact, graph);
			sm_fact->remap = Malloc(graph->n, int);
			Memcpy((void *) sm_fact->remap, (void *) sm_fact->PARDISO_fact->pstore[GMRFLib_PSTORE_TNUM_REF]->perm,
			       graph->n * sizeof(int));
		}
			break;

		case GMRFLib_SMTP_STILES:
		{
			GMRFLib_stiles_idx_tp stiles_idx = { 0, 0, 0 };
			sm_fact->remap = Malloc(graph->n, int);
			Memcpy((void *) sm_fact->remap, (void *) GMRFLib_stiles_get_perm(&stiles_idx), graph->n * sizeof(int));
		}
			break;

		default:
			GMRFLib_ASSERT(1 == 0, GMRFLib_ESNH);
			break;
		}
	}
		break;

	case GMRFLib_REORDER_PARDISO:
	{
		switch (sm_fact->smtp) {
		case GMRFLib_SMTP_PARDISO:		       /* same code as above */
		{
			if (sm_fact->PARDISO_fact == NULL) {
				GMRFLib_pardiso_init(&(sm_fact->PARDISO_fact));
			}
			GMRFLib_pardiso_reorder(sm_fact->PARDISO_fact, graph);
			sm_fact->remap = Malloc(graph->n, int);
			Memcpy((void *) sm_fact->remap, (void *) sm_fact->PARDISO_fact->pstore[GMRFLib_PSTORE_TNUM_REF]->perm,
			       graph->n * sizeof(int));
		}
			break;

		default:
			GMRFLib_ASSERT(1 == 0, GMRFLib_ESNH);
			break;
		}
	}
		break;

	case GMRFLib_REORDER_BAND:
		GMRFLib_EWRAP1(GMRFLib_compute_reordering_BAND(&(sm_fact->remap), graph));
		break;

		/*
		 * all the remaining ones are treated by the _TAUCS routine 
		 */
	case GMRFLib_REORDER_IDENTITY:
	case GMRFLib_REORDER_REVERSE_IDENTITY:
	case GMRFLib_REORDER_METIS:
	case GMRFLib_REORDER_GENMMD:
	case GMRFLib_REORDER_AMD:
	case GMRFLib_REORDER_AMDC:
	case GMRFLib_REORDER_AMDBAR:
	case GMRFLib_REORDER_AMDBARC:
	case GMRFLib_REORDER_MD:
	case GMRFLib_REORDER_MMD:
	{
		GMRFLib_EWRAP1(GMRFLib_compute_reordering_TAUCS(&(sm_fact->remap), graph, r, gn_ptr));
	}
		break;

	case GMRFLib_REORDER_STILES:
	{
		GMRFLib_ASSERT(0 == 1, GMRFLib_ESNH);
	}
		break;

	default:
	{
		GMRFLib_ASSERT(0 == 1, GMRFLib_ESNH);
	}
		break;
	}

	if (sm_fact->remap) {				       /* need this still for the wa-routines. FIXME */
		GMRFLib_EWRAP1(GMRFLib_graph_comp_bw(&(sm_fact->bandwidth), graph, sm_fact->remap));
	} else {
		sm_fact->bandwidth = -1;
	}

	GMRFLib_LEAVE_FUNCTION;

	return GMRFLib_SUCCESS;
}

int GMRFLib_free_reordering(GMRFLib_sm_fact_tp *sm_fact)
{
	if (sm_fact) {
		Free(sm_fact->remap);
		sm_fact->bandwidth = 0;
	}
	return GMRFLib_SUCCESS;
}

int GMRFLib_build_sparse_matrix(int thread_id, GMRFLib_sm_fact_tp *sm_fact, GMRFLib_Qfunc_tp *Qfunc, void *Qfunc_arg,
				GMRFLib_graph_tp *graph, GMRFLib_problem_tp *problem)
{
	GMRFLib_ENTER_FUNCTION;
	int ret = 0;

	switch (sm_fact->smtp) {
	case GMRFLib_SMTP_BAND:
	{
		ret = GMRFLib_build_sparse_matrix_BAND(thread_id, &(sm_fact->bchol), Qfunc, Qfunc_arg, graph, sm_fact->remap, sm_fact->bandwidth);
		if (ret != GMRFLib_SUCCESS) {
			return ret;
		}
	}
		break;

	case GMRFLib_SMTP_TAUCS:
	{
		ret = GMRFLib_build_sparse_matrix_TAUCS(thread_id, &(sm_fact->TAUCS_L), Qfunc, Qfunc_arg, graph, sm_fact->remap);
		if (ret != GMRFLib_SUCCESS) {
			return ret;
		}
	}
		break;

	case GMRFLib_SMTP_PARDISO:
	{
		if (sm_fact->PARDISO_fact == NULL) {
			GMRFLib_pardiso_init(&(sm_fact->PARDISO_fact));
			GMRFLib_pardiso_reorder(sm_fact->PARDISO_fact, graph);
		}
		ret = GMRFLib_pardiso_build(thread_id, sm_fact->PARDISO_fact, graph, Qfunc, Qfunc_arg);
		if (ret != GMRFLib_SUCCESS) {
			return ret;
		}
	}
		break;

	case GMRFLib_SMTP_STILES:
	{
		ret = GMRFLib_stiles_build(problem->stiles_idx, thread_id, Qfunc, Qfunc_arg);
		if (ret != GMRFLib_SUCCESS) {
			return ret;
		}
	}
		break;

	default:
		GMRFLib_ASSERT(1 == 0, GMRFLib_ESNH);
		break;
	}

	GMRFLib_LEAVE_FUNCTION;
	return GMRFLib_SUCCESS;
}

int GMRFLib_factorise_sparse_matrix(GMRFLib_sm_fact_tp *sm_fact, GMRFLib_graph_tp *graph, GMRFLib_problem_tp *problem)
{
	int ret;
	GMRFLib_ENTER_FUNCTION;

	switch (sm_fact->smtp) {
	case GMRFLib_SMTP_BAND:
	{
		ret = GMRFLib_factorise_sparse_matrix_BAND(sm_fact->bchol, &(sm_fact->finfo), graph, sm_fact->bandwidth);
		if (ret != GMRFLib_SUCCESS) {
			return ret;
		}
	}
		break;

	case GMRFLib_SMTP_TAUCS:
	{
		ret =
		    GMRFLib_factorise_sparse_matrix_TAUCS(&(sm_fact->TAUCS_L), &(sm_fact->TAUCS_symb_fact), &(sm_fact->TAUCS_cache),
							  &(sm_fact->finfo));
		if (ret != GMRFLib_SUCCESS) {
			return ret;
		}
	}
		break;

	case GMRFLib_SMTP_PARDISO:
	{
		ret = GMRFLib_pardiso_chol(sm_fact->PARDISO_fact);
		if (ret != GMRFLib_SUCCESS) {
			return ret;
		}
	}
		break;

	case GMRFLib_SMTP_STILES:
	{
		ret = GMRFLib_stiles_chol(problem->stiles_idx);
		if (ret != GMRFLib_SUCCESS) {
			return ret;
		}
	}
		break;

	default:
		GMRFLib_ASSERT(1 == 0, GMRFLib_ESNH);
		break;
	}

	GMRFLib_LEAVE_FUNCTION;

	return GMRFLib_SUCCESS;
}

int GMRFLib_free_fact_sparse_matrix(GMRFLib_sm_fact_tp *sm_fact)
{
	if (sm_fact) {
		switch (sm_fact->smtp) {
		case GMRFLib_SMTP_BAND:
		{
			GMRFLib_EWRAP1(GMRFLib_free_fact_sparse_matrix_BAND(sm_fact->bchol));
			sm_fact->bchol = NULL;
		}
			break;

		case GMRFLib_SMTP_TAUCS:
		{
			GMRFLib_free_fact_sparse_matrix_TAUCS(sm_fact->TAUCS_L, sm_fact->TAUCS_LL, sm_fact->TAUCS_symb_fact);
			GMRFLib_taucs_cache_free(sm_fact->TAUCS_cache);
			sm_fact->TAUCS_L = NULL;
			sm_fact->TAUCS_symb_fact = NULL;
			sm_fact->TAUCS_symb_fact = NULL;
			sm_fact->TAUCS_cache = NULL;
		}
			break;

		case GMRFLib_SMTP_PARDISO:
		{
			if (sm_fact->PARDISO_fact) {
				GMRFLib_pardiso_free(&(sm_fact->PARDISO_fact));
				sm_fact->PARDISO_fact = NULL;
			}
		}
			break;

		case GMRFLib_SMTP_STILES:
			break;

		default:
			GMRFLib_ASSERT(1 == 0, GMRFLib_ESNH);
			break;
		}
	}

	return GMRFLib_SUCCESS;
}

int GMRFLib_solve_l_sparse_matrix(double *rhs, int nrhs, GMRFLib_sm_fact_tp *sm_fact, GMRFLib_graph_tp *graph, GMRFLib_problem_tp *problem)
{
	/*
	 * rhs in real world. solve L x=rhs, rhs is overwritten by the solution 
	 */
	GMRFLib_ENTER_FUNCTION;

	switch (sm_fact->smtp) {
	case GMRFLib_SMTP_BAND:
	{
		omp_set_num_threads(GMRFLib_openmp->max_threads_inner);
#pragma omp parallel for num_threads(GMRFLib_openmp->max_threads_inner) schedule(static)
		for (int i = 0; i < nrhs; i++) {
			GMRFLib_solve_l_sparse_matrix_BAND(&rhs[i * graph->n], sm_fact->bchol, graph, sm_fact->remap, sm_fact->bandwidth);
		}
	}
		break;

	case GMRFLib_SMTP_TAUCS:
	{
		omp_set_num_threads(GMRFLib_openmp->max_threads_inner);
#pragma omp parallel for num_threads(GMRFLib_openmp->max_threads_inner) schedule(static)
		for (int i = 0; i < nrhs; i++) {
			GMRFLib_solve_l_sparse_matrix_TAUCS(&rhs[i * graph->n], sm_fact->TAUCS_L, graph, sm_fact->remap);
		}
	}
		break;

	case GMRFLib_SMTP_PARDISO:
	{
		GMRFLib_pardiso_solve_L(sm_fact->PARDISO_fact, rhs, rhs, nrhs);
	}
		break;

	case GMRFLib_SMTP_STILES:
	{
		GMRFLib_stiles_idx_tp stiles_idx = { problem->stiles_idx->in_group, 0, 0 };
		int err = GMRFLib_stiles_set_idx(&stiles_idx, nrhs);
		if (err == GMRFLib_SUCCESS) {
			GMRFLib_stiles_solve_L(&stiles_idx, rhs);
		} else {
			GMRFLib_stiles_set_idx(&stiles_idx, 1);
			for (int k = 0; k < nrhs; k++) {
				GMRFLib_stiles_solve_L(&stiles_idx, rhs + k * graph->n);
			}
		}
	}
		break;

	default:
		GMRFLib_ERROR(GMRFLib_ESNH);
		break;
	}

	GMRFLib_LEAVE_FUNCTION;

	return GMRFLib_SUCCESS;
}

int GMRFLib_solve_lt_sparse_matrix(double *rhs, int nrhs, GMRFLib_sm_fact_tp *sm_fact, GMRFLib_graph_tp *graph, GMRFLib_problem_tp *problem)
{
	/*
	 * rhs in real world. solve L^Tx=rhs, rhs is overwritten by the solution 
	 */
	GMRFLib_ENTER_FUNCTION;

	switch (sm_fact->smtp) {
	case GMRFLib_SMTP_BAND:
	{
		omp_set_num_threads(GMRFLib_openmp->max_threads_inner);
#pragma omp parallel for num_threads(GMRFLib_openmp->max_threads_inner) schedule(static)
		for (int i = 0; i < nrhs; i++) {
			GMRFLib_solve_lt_sparse_matrix_BAND(&rhs[i * graph->n], sm_fact->bchol, graph, sm_fact->remap, sm_fact->bandwidth);
		}
	}
		break;

	case GMRFLib_SMTP_TAUCS:
	{
		omp_set_num_threads(GMRFLib_openmp->max_threads_inner);
#pragma omp parallel for num_threads(GMRFLib_openmp->max_threads_inner) schedule(static)
		for (int i = 0; i < nrhs; i++) {
			GMRFLib_solve_lt_sparse_matrix_TAUCS(&rhs[i * graph->n], sm_fact->TAUCS_L, graph, sm_fact->remap);
		}
	}
		break;

	case GMRFLib_SMTP_PARDISO:
	{
		GMRFLib_pardiso_solve_LT(sm_fact->PARDISO_fact, rhs, rhs, nrhs);
	}
		break;

	case GMRFLib_SMTP_STILES:
	{
		GMRFLib_stiles_idx_tp stiles_idx = { problem->stiles_idx->in_group, 0, 0 };
		int err = GMRFLib_stiles_set_idx(&stiles_idx, nrhs);
		if (err == GMRFLib_SUCCESS) {
			GMRFLib_stiles_solve_LT(&stiles_idx, rhs);
		} else {
			GMRFLib_stiles_set_idx(&stiles_idx, 1);
			for (int k = 0; k < nrhs; k++) {
				GMRFLib_stiles_solve_LT(&stiles_idx, rhs + k * graph->n);
			}
		}
	}
		break;

	default:
		GMRFLib_ERROR(GMRFLib_ESNH);
		break;
	}

	GMRFLib_LEAVE_FUNCTION;

	return GMRFLib_SUCCESS;
}

int GMRFLib_solve_llt_sparse_matrix(double *rhs, int nrhs, GMRFLib_sm_fact_tp *sm_fact, GMRFLib_graph_tp *graph,
				    GMRFLib_problem_tp *problem, GMRFLib_stiles_idx_tp *stiles_idx)
{
	/*
	 * rhs in real world. solve Q x=rhs, where Q=L L^T 
	 */
	GMRFLib_ENTER_FUNCTION;

	static double **wwork = NULL;
	static int *wwork_len = NULL;
	if (!wwork) {
#pragma omp critical (Name_c02cfe7c85f984ba167d3d158f5219787998c27f)
		if (!wwork) {
			wwork_len = Calloc(GMRFLib_CACHE_LEN(), int);
			double **tmp = Calloc(GMRFLib_CACHE_LEN(), double *);
			wwork = tmp;
		}
	}

	int cache_idx = 0;
	GMRFLib_CACHE_SET_IDX(cache_idx);

	int nw = graph->n * nrhs;
	if (nw > wwork_len[cache_idx]) {
		int numa_node = GMRFLib_numa_get_node();
		if (wwork[cache_idx] && wwork_len[cache_idx] > 0) {
			GMRFLib_numa_free(wwork[cache_idx], wwork_len[cache_idx]);
		}
		wwork_len[cache_idx] = nw;
		wwork[cache_idx] = GMRFLib_numa_alloc_onnode(wwork_len[cache_idx] * sizeof(double), numa_node);
	}
	double *work = wwork[cache_idx];

	if (sm_fact->smtp == GMRFLib_SMTP_BAND) {
		omp_set_num_threads(GMRFLib_openmp->max_threads_inner);
#pragma omp parallel for num_threads(GMRFLib_openmp->max_threads_inner) schedule(static)
		for (int i = 0; i < nrhs; i++) {
			int offset = i * graph->n;
			GMRFLib_solve_llt_sparse_matrix_BAND(rhs + offset, sm_fact->bchol, graph, sm_fact->remap,
							     sm_fact->bandwidth, work + offset);
		}
	} else if (sm_fact->smtp == GMRFLib_SMTP_TAUCS) {
		if (nrhs == 1) {
			// simple entry
			GMRFLib_solve_llt_sparse_matrix_TAUCS(rhs, sm_fact->TAUCS_L, sm_fact->TAUCS_LL, graph, sm_fact->remap, work);
		} else {

			int ntt = -1;
			if (omp_get_level() == 0) {
				ntt = GMRFLib_PARDISO_MAX_NUM_THREADS();
			} else {
				ntt = GMRFLib_openmp->max_threads_inner;
			}

			// default
			int block_size = GMRFLib_taucs_get_ctl_ptr()->block_size;
			if (block_size <= 0)
				block_size = 8;

			int target = IMAX(1, block_size);
			ntt = IMIN(ntt, IMAX(1, nrhs / target));
			int chunk_size = nrhs / ntt;	       /* integer division */
			if (chunk_size == 0 || ntt == 1) {
				GMRFLib_solve_llt_sparse_matrix2_TAUCS(rhs, sm_fact->TAUCS_L, graph, sm_fact->remap, nrhs, work);
			} else {
				// split work in ntt threads each doing a chunk
				int *iwork = Malloc(2 * ntt, int);
				int *csize = iwork;
				int *offset = iwork + ntt;

				GMRFLib_ifill(ntt, chunk_size, csize);
				for (int i = 0; i < ntt; i++) {
					if (GMRFLib_isum(ntt, csize) < nrhs) {
						csize[i]++;
					}
				}
				assert(GMRFLib_isum(ntt, csize) == nrhs);
				offset[0] = 0;
				for (int i = 1; i < ntt; i++) {
					offset[i] = offset[i - 1] + csize[i - 1] * graph->n;
				}

#pragma omp parallel for num_threads(ntt) schedule(static)
				for (int i = 0; i < ntt; i++) {
					GMRFLib_solve_llt_sparse_matrix2_TAUCS(rhs + offset[i], sm_fact->TAUCS_L, graph,
									       sm_fact->remap, csize[i], work + offset[i]);
				}
				Free(iwork);
			}
		}
	} else if (sm_fact->smtp == GMRFLib_SMTP_PARDISO) {
		GMRFLib_EWRAP1(GMRFLib_pardiso_solve_LLT(sm_fact->PARDISO_fact, rhs, rhs, nrhs));
	} else if (sm_fact->smtp == GMRFLib_SMTP_STILES) {
		GMRFLib_stiles_idx_tp s_idx = { 0, 0, 0 };
		if (stiles_idx) {
			s_idx.in_group = stiles_idx->in_group;
		} else {
			s_idx.in_group = problem->stiles_idx->in_group;
		}

		int err = GMRFLib_stiles_set_idx(&s_idx, nrhs);
		if (err == GMRFLib_SUCCESS) {
			GMRFLib_stiles_solve_LLT(&s_idx, rhs);
		} else {
			GMRFLib_stiles_set_idx(&s_idx, 1);
			for (int k = 0; k < nrhs; k++) {
				GMRFLib_stiles_solve_LLT(&s_idx, rhs + k * graph->n);
			}
		}
	} else {
		GMRFLib_ERROR(GMRFLib_ESNH);
	}
	GMRFLib_LEAVE_FUNCTION;

	return GMRFLib_SUCCESS;
}

int GMRFLib_solve_llt_sparse_matrix_special(double *rhs, GMRFLib_sm_fact_tp *sm_fact, GMRFLib_graph_tp *graph, int idx, GMRFLib_problem_tp *problem)
{
	/*
	 * rhs in real world. solve Q x=rhs, where Q=L L^T. BUT, here we know that rhs is 0 execpt for a 1 at index idx.
	 */
	GMRFLib_ENTER_FUNCTION;

	switch (sm_fact->smtp) {
	case GMRFLib_SMTP_BAND:
	{
		GMRFLib_solve_llt_sparse_matrix_special_BAND(rhs, sm_fact->bchol, graph, sm_fact->remap, sm_fact->bandwidth, idx);
	}
		break;

	case GMRFLib_SMTP_TAUCS:
	{
		GMRFLib_solve_llt_sparse_matrix_special_TAUCS(rhs, sm_fact->TAUCS_L, graph, sm_fact->remap, idx);
	}
		break;

	case GMRFLib_SMTP_PARDISO:
	{
		GMRFLib_pardiso_solve_LLT(sm_fact->PARDISO_fact, rhs, rhs, 1);
	}
		break;

	case GMRFLib_SMTP_STILES:
	{
		GMRFLib_stiles_idx_tp stiles_idx = { problem->stiles_idx->in_group, 0, 0 };
		GMRFLib_stiles_set_idx(&stiles_idx, 1);
		GMRFLib_stiles_solve_LLT(&stiles_idx, rhs);
	}
		break;

	default:
		GMRFLib_ERROR(GMRFLib_ESNH);
		break;
	}

	GMRFLib_LEAVE_FUNCTION;

	return GMRFLib_SUCCESS;
}

int GMRFLib_solve_lt_sparse_matrix_special(double *rhs, GMRFLib_sm_fact_tp *sm_fact, GMRFLib_graph_tp *graph, int findx, int toindx,
					   int remapped, GMRFLib_problem_tp *problem)
{
	/*
	 * rhs in real world, bchol in mapped world. solve L^Tx=b backward only from rhs[findx] up to rhs[toindx]. note that
	 * findx and toindx is in mapped world. if remapped, do not remap/remap-back the rhs before solving.
	 */
	GMRFLib_ENTER_FUNCTION;
	switch (sm_fact->smtp) {
	case GMRFLib_SMTP_BAND:
	{
		GMRFLib_EWRAP0(GMRFLib_solve_lt_sparse_matrix_special_BAND
			       (rhs, sm_fact->bchol, graph, sm_fact->remap, sm_fact->bandwidth, findx, toindx, remapped));
	}
		break;

	case GMRFLib_SMTP_TAUCS:
	{
		GMRFLib_EWRAP0(GMRFLib_solve_lt_sparse_matrix_special_TAUCS(rhs, sm_fact->TAUCS_L, graph, sm_fact->remap, findx, toindx, remapped));
	}
		break;

	case GMRFLib_SMTP_PARDISO:
	{
		GMRFLib_pardiso_solve_LT(sm_fact->PARDISO_fact, rhs, rhs, 1);
	}
		break;

	case GMRFLib_SMTP_STILES:
	{
		GMRFLib_stiles_idx_tp stiles_idx = { problem->stiles_idx->in_group, 0, 0 };
		GMRFLib_stiles_set_idx(&stiles_idx, 1);
		GMRFLib_stiles_solve_LT(&stiles_idx, rhs);
	}
		break;

	default:
	{
		GMRFLib_ERROR(GMRFLib_ESNH);
	}
		break;
	}
	GMRFLib_LEAVE_FUNCTION;

	return GMRFLib_SUCCESS;
}

int GMRFLib_solve_l_sparse_matrix_special(double *rhs, GMRFLib_sm_fact_tp *sm_fact, GMRFLib_graph_tp *graph, int findx, int toindx,
					  int remapped, GMRFLib_problem_tp *problem)
{
	/*
	 * rhs in real world, bchol in mapped world. solve Lx=b backward only from rhs[findx] up to rhs[toindx]. note that
	 * findx and toindx is in mapped world. if remapped, do not remap/remap-back the rhs before solving.
	 */

	GMRFLib_ENTER_FUNCTION;
	switch (sm_fact->smtp) {
	case GMRFLib_SMTP_BAND:
	{
		GMRFLib_EWRAP0(GMRFLib_solve_l_sparse_matrix_special_BAND
			       (rhs, sm_fact->bchol, graph, sm_fact->remap, sm_fact->bandwidth, findx, toindx, remapped));
	}
		break;

	case GMRFLib_SMTP_TAUCS:
	{
		GMRFLib_EWRAP0(GMRFLib_solve_l_sparse_matrix_special_TAUCS(rhs, sm_fact->TAUCS_L, graph, sm_fact->remap, findx, toindx, remapped));
	}
		break;

	case GMRFLib_SMTP_PARDISO:
	{
		if (remapped) {
			GMRFLib_pardiso_perm(rhs, 1, sm_fact->PARDISO_fact);
		}
		GMRFLib_pardiso_solve_L(sm_fact->PARDISO_fact, rhs, rhs, 1);
	}
		break;

	case GMRFLib_SMTP_STILES:
	{
		if (remapped) {
			FIXME("UNSURE ABOUT THIS THING");
			int *perm = GMRFLib_stiles_get_perm(problem->stiles_idx);
			double *y = Malloc(graph->n, double);
			Memcpy(y, rhs, graph->n * sizeof(double));
			GMRFLib_pack(graph->n, y, perm, rhs);
			Free(y);
		}
		GMRFLib_stiles_idx_tp stiles_idx = { problem->stiles_idx->in_group, 0, 0 };
		GMRFLib_stiles_set_idx(&stiles_idx, 1);
		GMRFLib_stiles_solve_L(&stiles_idx, rhs);
	}
		break;

	default:
		GMRFLib_ERROR(GMRFLib_ESNH);
		break;
	}

	GMRFLib_LEAVE_FUNCTION;
	return GMRFLib_SUCCESS;
}

int GMRFLib_log_determinant(double *logdet, GMRFLib_sm_fact_tp *sm_fact, GMRFLib_graph_tp *graph, GMRFLib_problem_tp *problem)
{
	switch (sm_fact->smtp) {
	case GMRFLib_SMTP_BAND:
	{
		GMRFLib_EWRAP0(GMRFLib_log_determinant_BAND(logdet, sm_fact->bchol, graph, sm_fact->bandwidth));
	}
		break;

	case GMRFLib_SMTP_TAUCS:
	{
		GMRFLib_EWRAP0(GMRFLib_log_determinant_TAUCS(logdet, sm_fact->TAUCS_L));
	}
		break;

	case GMRFLib_SMTP_PARDISO:
	{
		*logdet = GMRFLib_pardiso_logdet(sm_fact->PARDISO_fact);
	}
		break;

	case GMRFLib_SMTP_STILES:
	{
		*logdet = GMRFLib_stiles_logdet(problem->stiles_idx);
	}
		break;

	default:
		GMRFLib_ERROR(GMRFLib_ESNH);
		break;
	}

	return GMRFLib_SUCCESS;
}

int GMRFLib_comp_cond_meansd(double *cmean, double *csd, int indx, double *x, int remapped, GMRFLib_sm_fact_tp *sm_fact, GMRFLib_graph_tp *graph)
{
	GMRFLib_ENTER_FUNCTION;
	switch (sm_fact->smtp) {
	case GMRFLib_SMTP_BAND:
	{
		GMRFLib_EWRAP1(GMRFLib_comp_cond_meansd_BAND
			       (cmean, csd, indx, x, remapped, sm_fact->bchol, graph, sm_fact->remap, sm_fact->bandwidth));
	}
		break;

	case GMRFLib_SMTP_TAUCS:
	{
		GMRFLib_EWRAP1(GMRFLib_comp_cond_meansd_TAUCS(cmean, csd, indx, x, remapped, sm_fact->TAUCS_L, graph, sm_fact->remap));
	}
		break;

	case GMRFLib_SMTP_PARDISO:
	case GMRFLib_SMTP_STILES:
	{
		assert(0 == 1);
	}
		break;

	default:
		GMRFLib_ERROR(GMRFLib_ESNH);
		break;
	}
	GMRFLib_LEAVE_FUNCTION;

	return GMRFLib_SUCCESS;
}

int GMRFLib_bitmap_factorisation(const char *filename_body, GMRFLib_sm_fact_tp *sm_fact, GMRFLib_graph_tp *graph)
{
	switch (sm_fact->smtp) {
	case GMRFLib_SMTP_BAND:
	{
		GMRFLib_EWRAP1(GMRFLib_bitmap_factorisation_BAND(filename_body, sm_fact->bchol, graph, sm_fact->remap, sm_fact->bandwidth));
	}
		break;

	case GMRFLib_SMTP_TAUCS:
	{
		GMRFLib_EWRAP1(GMRFLib_bitmap_factorisation_TAUCS(filename_body, sm_fact->TAUCS_L));
	}
		break;

	case GMRFLib_SMTP_PARDISO:
	{
		GMRFLib_EWRAP1(GMRFLib_pardiso_bitmap());
	}
		break;

	case GMRFLib_SMTP_STILES:
		break;

	default:
		GMRFLib_ERROR(GMRFLib_ESNH);
		break;
	}
	return GMRFLib_SUCCESS;
}

int GMRFLib_compute_Qinv(void *problem)
{
	GMRFLib_ENTER_FUNCTION;
	GMRFLib_problem_tp *p = (GMRFLib_problem_tp *) problem;

	switch (p->sub_sm_fact.smtp) {
	case GMRFLib_SMTP_BAND:
	{
		GMRFLib_EWRAP0(GMRFLib_compute_Qinv_BAND(p));
	}
		break;

	case GMRFLib_SMTP_TAUCS:
	{
		GMRFLib_EWRAP0(GMRFLib_compute_Qinv_TAUCS(p));
	}
		break;

	case GMRFLib_SMTP_PARDISO:
	{
		GMRFLib_pardiso_Qinv_INLA(p);
	}
		break;

	case GMRFLib_SMTP_STILES:
	{
		GMRFLib_stiles_Qinv_INLA(p);
	}
		break;

	default:
		GMRFLib_ERROR(GMRFLib_ESNH);
		break;
	}
	GMRFLib_LEAVE_FUNCTION;
	return GMRFLib_SUCCESS;
}

int GMRFLib_valid_smtp(int smtp)
{
	return ((smtp == GMRFLib_SMTP_BAND || smtp == GMRFLib_SMTP_TAUCS || smtp == GMRFLib_SMTP_PARDISO
		 || smtp == GMRFLib_SMTP_STILES) ? GMRFLib_TRUE : GMRFLib_FALSE);
}

const char *GMRFLib_reorder_name(GMRFLib_reorder_tp r)
{
	switch (r) {
	case GMRFLib_REORDER_DEFAULT:
		return "default";
	case GMRFLib_REORDER_IDENTITY:
		return "identity";
	case GMRFLib_REORDER_REVERSE_IDENTITY:
		return "reverseidentity";
	case GMRFLib_REORDER_BAND:
		return "band";
	case GMRFLib_REORDER_METIS:
		return "metis";
	case GMRFLib_REORDER_GENMMD:
		return "genmmd";
	case GMRFLib_REORDER_AMD:
		return "amd";
	case GMRFLib_REORDER_MD:
		return "md";
	case GMRFLib_REORDER_MMD:
		return "mmd";
	case GMRFLib_REORDER_AMDBAR:
		return "amdbar";
	case GMRFLib_REORDER_AMDC:
		return "amdc";
	case GMRFLib_REORDER_AMDBARC:
		return "amdbarc";
	case GMRFLib_REORDER_PARDISO:
		return "pardiso";
	case GMRFLib_REORDER_STILES:
		return "stiles";
	default:
		fprintf(stderr, "\n\t*** ERROR *** Reordering [%d] not defined.\n", r);
		GMRFLib_ASSERT_RETVAL(0 == 1, GMRFLib_EPARAMETER, "(unknown reording)");
	}

	return "(unknown reording)";
}

int GMRFLib_reorder_id(const char *name)
{
	if (!strcasecmp(name, "default") || !strcasecmp(name, "auto"))
		return GMRFLib_REORDER_DEFAULT;
	else if (!strcasecmp(name, "identity"))
		return GMRFLib_REORDER_IDENTITY;
	else if (!strcasecmp(name, "reverseidentity"))
		return GMRFLib_REORDER_REVERSE_IDENTITY;
	else if (!strcasecmp(name, "band"))
		return GMRFLib_REORDER_BAND;
	else if (!strcasecmp(name, "metis"))
		return GMRFLib_REORDER_METIS;
	else if (!strcasecmp(name, "genmmd"))
		return GMRFLib_REORDER_GENMMD;
	else if (!strcasecmp(name, "amd"))
		return GMRFLib_REORDER_AMD;
	else if (!strcasecmp(name, "amdbar"))
		return GMRFLib_REORDER_AMDBAR;
	else if (!strcasecmp(name, "md"))
		return GMRFLib_REORDER_MD;
	else if (!strcasecmp(name, "mmd"))
		return GMRFLib_REORDER_MMD;
	else if (!strcasecmp(name, "amdc"))
		return GMRFLib_REORDER_AMDC;
	else if (!strcasecmp(name, "amdbarc"))
		return GMRFLib_REORDER_AMDBARC;
	else if (!strcasecmp(name, "pardiso"))
		return GMRFLib_REORDER_PARDISO;
	else if (!strcasecmp(name, "stiles"))
		return GMRFLib_REORDER_STILES;
	else {
		fprintf(stderr, "\n\t*** ERROR *** Reordering [%s] not defined.\n", name);
		GMRFLib_ASSERT_RETVAL(0 == 1, GMRFLib_EPARAMETER, -1);
	}

	return -1;
}
