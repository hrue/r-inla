int GMRFLib_duplicate_pardiso_store_ORIG(GMRFLib_pardiso_store_tp **nnew, GMRFLib_pardiso_store_tp *old, int UNUSED(copy_ptr), int copy_pardiso_ptr)
{
	int tnum = omp_get_thread_num();
	// if copy_pardiso_ptr, then copy the ptr to read-only objects. 'copy_ptr' is NOT USED
	int debug = S.debug, failsafe_mode = 0;
	if (old == NULL) {
		*nnew = NULL;
		return GMRFLib_SUCCESS;
	}

	if (failsafe_mode) {
		// FIXME("-->duplicate by creating a new one each time");
		GMRFLib_pardiso_init(nnew);
		GMRFLib_pardiso_reorder(*nnew, old->graph);
		return GMRFLib_SUCCESS;
	}
	GMRFLib_ENTER_ROUTINE;

	if (copy_pardiso_ptr) {
#define CP(_what) dup->_what = old->_what
#define CP2(_what) dup->pstore[tnum]->_what = old->pstore[tnum]->_what
#define CP2_ref(_what) dup->pstore[GMRFLib_PSTORE_TNUM_REF]->_what = old->pstore[GMRFLib_PSTORE_TNUM_REF]->_what
#define CPv_ref(_what, type, len)					\
		if (old->pstore[GMRFLib_PSTORE_TNUM_REF]->_what) {			\
			dup->pstore[GMRFLib_PSTORE_TNUM_REF]->_what = Calloc(len, type); \
			Memcpy((void *) (dup->pstore[GMRFLib_PSTORE_TNUM_REF]->_what), (void *) (old->pstore[GMRFLib_PSTORE_TNUM_REF]->_what), (len) * sizeof(type)); \
		} else {						\
			dup->pstore[GMRFLib_PSTORE_TNUM_REF]->_what = NULL;		\
		}							\

		GMRFLib_pardiso_store_tp *dup = Calloc(1, GMRFLib_pardiso_store_tp);
		dup->copy_pardiso_ptr = 1;		       /* YES! */
		for (int i = 0; i < GMRFLib_PARDISO_PLEN; i++) {
			CP(pt[i]);
		}
		CP(iparm_default);
		CP(dparm_default);
		CP(maxfct);
		CP(done_with_init);
		CP(done_with_reorder);
		CP(msglvl);
		CP(mtype);
		CP(solver);
		CP(graph);

		if (!dup->pstore)
			dup->pstore = Calloc(GMRFLib_MAX_THREADS(), GMRFLib_pardiso_store_pr_thread_tp *);
		if (!dup->pstore[tnum])
			dup->pstore[tnum] = Calloc(1, GMRFLib_pardiso_store_pr_thread_tp);
		if (!dup->pstore[GMRFLib_PSTORE_TNUM_REF])
			dup->pstore[GMRFLib_PSTORE_TNUM_REF] = Calloc(1, GMRFLib_pardiso_store_pr_thread_tp);
		for (int i = 0; i < GMRFLib_PARDISO_PLEN; i++) {
			CP2_ref(iparm[i]);
			CP2_ref(dparm[i]);
		}
		CP2_ref(done_with_build);
		CP2_ref(done_with_chol);
		CP2(dummy);
		CP2(err_code);
		CP2(idummy);
		CP2(nrhs);
		CP2(phase);
		CP2(L_nnz);
		CP2(perm_identity);

		CPv_ref(perm, int, old->graph->n);
		CPv_ref(iperm, int, old->graph->n);

		CP2(log_det_Q);
		CP2_ref(Q);
		CP2_ref(Qinv);
#undef CP
#undef CP2
#undef CP2_ref
#undef CPv
#undef CPv_ref
		*nnew = dup;
		return GMRFLib_SUCCESS;
	}

	if (S.static_pstores == NULL) {
#pragma omp critical (Name_046c40f5fd2e202479d5c486dfdf986558c6e681)
		{
			if (S.static_pstores == NULL) {
				if (S.s_verbose) {
					printf("==> init static_pstores\n");
				}
				S.busy = Calloc(PSTORES_NUM, int);
				S.static_pstores = Calloc(PSTORES_NUM, GMRFLib_pardiso_store_tp *);
			}
		}
	}

	int found = 0, idx = -1, ok = 0;
#pragma omp critical (Name_29873421a20baf5374230cb83067fb5810469371)
	{
		for (int i = 0; i < PSTORES_NUM && !found; i++) {
			if (!S.busy[i]) {
				S.busy[i] = 1;
				idx = i;
				found = 1;
			}
		}
	}

	assert(found == 1);
	if (S.static_pstores[idx]) {
		if (debug) {
			printf("%s:%1d: static_pstores...iparm[2] = %1d\n", __FILE__, __LINE__,
			       S.static_pstores[idx]->pstore[GMRFLib_PSTORE_TNUM_REF]->iparm[2]);
			printf("%s:%1d: level %d max_threads_nested = %1d\n", __FILE__, __LINE__, omp_get_level(),
			       GMRFLib_openmp->max_threads_nested[1]);
		}
		ok = (S.static_pstores[idx]->pstore[GMRFLib_PSTORE_TNUM_REF]->iparm[2] >= GMRFLib_openmp->max_threads_nested[1]);
	} else {
		ok = 1;
	}
	if (!ok) {
		P(S.static_pstores[idx]->pstore[GMRFLib_PSTORE_TNUM_REF]->iparm[2]);
		P(GMRFLib_openmp->max_threads_nested[1]);
		FIXME("THIS IS NOT TRUE: iparm[2] >= threads_nested[1]");
	}

	if (S.static_pstores[idx] && ok) {
		*nnew = S.static_pstores[idx];
		if (S.s_verbose) {
			printf("==> reuse store[%1d]\n", idx);
		}
	} else {
		GMRFLib_pardiso_init(&(S.static_pstores[idx]));
		GMRFLib_pardiso_reorder(S.static_pstores[idx], old->graph);
		*nnew = S.static_pstores[idx];
		if (S.s_verbose) {
			printf("==> new store[%1d]\n", idx);
		}
	}

	if (S.s_verbose) {
		printf("duplicate: new=%p old=%p i=%1d\n", *((void **) nnew), ((void *) old), idx);
	}

	GMRFLib_LEAVE_ROUTINE;
	return GMRFLib_SUCCESS;
}
int GMRFLib_pardiso_free_ORIG(GMRFLib_pardiso_store_tp **store)
{
	int tnum = omp_get_thread_num();
	if (store == NULL || *store == NULL) {
		return GMRFLib_SUCCESS;
	}

	if ((*store)->copy_pardiso_ptr) {
		// this is special
		if (S.s_verbose) {
			FIXME("Free pardiso store with copy_pardiso_ptr = 1");
		}
		if ((*store)->pstore) {
			for (int i = 0; i < GMRFLib_MAX_THREADS(); i++) {
				if ((*store)->pstore[i]) {
					Free((*store)->pstore[i]->perm);
					Free((*store)->pstore[i]->iperm);
					Free((*store)->pstore[i]);
				}
			}
		}
		Free((*store)->pstore);
		Free((*store));
		*store = NULL;

		return GMRFLib_SUCCESS;
	}

	if (S.s_verbose) {
		PP("free: old=", *store);
	}

	int found = 0;
	if (S.static_pstores != NULL) {
		if (S.s_verbose) {
			for (int i = 0; i < PSTORES_NUM; i++) {
				if (S.busy[i]) {
					printf("in store: i=%1d s=%p\n", i, (void *) S.static_pstores[i]);
				}
			}
		}
		for (int i = 0; i < PSTORES_NUM && !found; i++) {

			if (S.static_pstores[i] == *store) {
				found = 1;
				if (S.busy[i]) {
					S.busy[i] = 0;
					if (S.s_verbose) {
						printf("==> S.busy[%1d] = 1\n", i);
						PP("S.static_pstores[i]", S.static_pstores[i]);
						PP("*store", *store);
						printf("==> free store[%1d]\n", i);
					}
				} else {
					if (S.s_verbose) {
						printf("==> this one is already free [%1d]. ignore\n", i);
					}
				}
			}
		}
	}
	if (!found) {
		if (S.s_verbose) {
			printf("==> free manually as not found\n");
		}

		if ((*store)->pstore[tnum]) {
			(*store)->pstore[tnum]->phase = -1;
			int mnum1 = 1;
			pardiso((*store)->pt, &((*store)->maxfct), &mnum1, &((*store)->mtype),
				&((*store)->pstore[tnum]->phase),
				&((*store)->pstore[tnum]->idummy),
				&((*store)->pstore[tnum]->dummy), &((*store)->pstore[tnum]->idummy),
				&((*store)->pstore[tnum]->idummy),
				&((*store)->pstore[tnum]->idummy),
				&((*store)->pstore[tnum]->nrhs), (*store)->pstore[tnum]->iparm, &((*store)->msglvl),
				NULL, NULL, &((*store)->pstore[tnum]->err_code), (*store)->pstore[tnum]->dparm);

			GMRFLib_csr_free(&((*store)->pstore[GMRFLib_PSTORE_TNUM_REF]->Q));
			GMRFLib_csr_free(&((*store)->pstore[GMRFLib_PSTORE_TNUM_REF]->Qinv));
			Free((*store)->pstore[tnum]);
		}

		GMRFLib_graph_free((*store)->graph);
		Free((*store)->iparm_default);
		Free((*store)->dparm_default);
		Free(*store);
	}

	return GMRFLib_SUCCESS;
}
