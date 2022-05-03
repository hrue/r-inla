//
// *** test programs ***
//

double my_pardiso_test_Q(int i, int j, double *UNUSED(values), void *arg)
{
	GMRFLib_graph_tp *graph = (GMRFLib_graph_tp *) arg;
	return (i == j ? graph->n + i : -1.0);
}

int my_pardiso_test1(void)
{
	int err = 0;
	int id = GMRFLib_thread_id;

	if (1) {
		err = GMRFLib_pardiso_check_install(1, 0);
		if (err == GMRFLib_SUCCESS) {
			printf("PARDISO OK\n");
		} else {
			printf("PARDISO FAIL\n");
		}
	}

	GMRFLib_tabulate_Qfunc_tp *Qtab;
	GMRFLib_graph_tp *g;

	GMRFLib_openmp->strategy = GMRFLib_OPENMP_STRATEGY_PARDISO;
	GMRFLib_openmp_implement_strategy(GMRFLib_OPENMP_PLACES_OPTIMIZE, NULL, NULL);

	// GMRFLib_tabulate_Qfunc_from_file(&Qtab, &g, "Qdense.txt", -1, NULL);
	// GMRFLib_tabulate_Qfunc_from_file(&Qtab, &g, "I5.txt", -1, NULL);
	GMRFLib_tabulate_Qfunc_from_file(&Qtab, &g, "Q.txt", -1, NULL);

	GMRFLib_csr_tp *csr, *csr2;
	GMRFLib_Q2csr(&csr, g, Qtab->Qfunc, Qtab->Qfunc_arg);
	GMRFLib_csr_print(stdout, csr);

	GMRFLib_csr_duplicate(&csr2, csr);
	// GMRFLib_csr_print(stdout, csr2);

	GMRFLib_csr2Q(&Qtab, &g, csr2);
	// GMRFLib_printf_Qfunc(stdout, g, Qtab->Qfunc, Qtab->Qfunc_arg);

	int *perm = NULL;
	int i, k, nrhs;

	// GMRFLib_printf_graph(stdout, g);
	GMRFLib_pardiso_store_tp *store = NULL;
	GMRFLib_pardiso_init(&store);
	GMRFLib_pardiso_reorder(store, g);
	GMRFLib_pardiso_build(store, g, Qtab->Qfunc, Qtab->Qfunc_arg);
	GMRFLib_pardiso_chol(store);
	GMRFLib_pardiso_Qinv(store);

	P(GMRFLib_openmp->max_threads_outer);
	P(GMRFLib_openmp->max_threads_inner);

	nrhs = 1;
	// S.s_verbose=1;
#pragma omp parallel for private(k) num_threads(GMRFLib_openmp->max_threads_outer)
	for (k = 0; k < 100; k++) {
		GMRFLib_thread_id = id;

		printf("this is k= %d from thread %d\n", k, omp_get_thread_num());
		GMRFLib_pardiso_store_tp *store2 = NULL;
		if (0) {
			GMRFLib_pardiso_init(&store2);
			GMRFLib_pardiso_reorder(store2, g);
			GMRFLib_pardiso_build(store2, g, Qtab->Qfunc, Qtab->Qfunc_arg);
			GMRFLib_pardiso_chol(store2);
		} else {
			GMRFLib_duplicate_pardiso_store(&store2, store, -1, 1);
		}

		int view = 1;
		double *x = Calloc(g->n * nrhs, double);
		double *b = Calloc(g->n * nrhs, double);

		for (i = 0; i < g->n * nrhs; i++) {
			b[i] = ISQR(i + 1);
		}
		GMRFLib_pardiso_solve_LLT(store2, x, b, nrhs);
		if (view) {
			printf("solve LLT\n");
			for (int ii = 0; ii < nrhs; ii++) {
				for (i = 0; i < g->n; i++) {
					printf("ii %d i %d  x= %f  b=%f\n", ii, i, x[i + ii * g->n], b[i + ii * g->n]);
				}
			}
		}

		for (i = 0; i < g->n * nrhs; i++) {
			b[i] = ISQR(i + 1);
		}
		GMRFLib_pardiso_solve_LT(store2, x, b, 1);
		if (view) {
			printf("solve LT\n");
			for (int ii = 0; ii < nrhs; ii++) {
				for (i = 0; i < g->n; i++) {
					printf("ii %d i %d  x= %f  b=%f\n", ii, i, x[i + ii * g->n], b[i + ii * g->n]);
				}
			}
		}

		for (i = 0; i < g->n * nrhs; i++) {
			b[i] = ISQR(i + 1);
		}
		GMRFLib_pardiso_solve_L(store2, x, b, 1);
		if (view) {
			printf("solve L\n");
			for (int ii = 0; ii < nrhs; ii++) {
				for (i = 0; i < g->n; i++) {
					printf("ii %d i %d  x= %f  b=%f\n", ii, i, x[i + ii * g->n], b[i + ii * g->n]);
				}
			}
		}

		GMRFLib_pardiso_free(&store2);
		Free(x);
		Free(b);
	}

	printf("call free\n");
	GMRFLib_free_tabulate_Qfunc(Qtab);
	GMRFLib_pardiso_free(&store);
	GMRFLib_csr_free(&csr);
	GMRFLib_csr_free(&csr2);
	GMRFLib_graph_free(g);
	Free(perm);

	exit(0);
}

int my_pardiso_test2(void)
{
	int n = 5, m = 1, nc = 1, i;
	double *var;

	GMRFLib_graph_tp *graph = NULL;
	GMRFLib_graph_mk_linear(&graph, n, m, 0);
	GMRFLib_problem_tp *problem = NULL;
	GMRFLib_constr_tp *constr = NULL;
	GMRFLib_make_empty_constr(&constr);

	constr->nc = nc;
	constr->a_matrix = Calloc(n * nc, double);
	for (i = 0; i < n * nc; i++)
		constr->a_matrix[i] = GMRFLib_uniform();
	constr->e_vector = Calloc(nc, double);
	for (i = 0; i < nc; i++)
		constr->e_vector[i] = GMRFLib_uniform();
	GMRFLib_prepare_constr(constr, graph, 1);

	// GMRFLib_openmp->strategy = GMRFLib_OPENMP_STRATEGY_PARDISO;
	// GMRFLib_smtp = GMRFLib_SMTP_PARDISO;
	GMRFLib_openmp_implement_strategy(GMRFLib_OPENMP_PLACES_BUILD_MODEL, NULL, NULL);

	double *x = Calloc(n, double);
	double *b = Calloc(n, double);
	double *c = Calloc(n, double);
	double *mean = Calloc(n, double);
	for (i = 0; i < n; i++) {
		x[i] = GMRFLib_uniform();
		b[i] = GMRFLib_uniform();
		c[i] = exp(GMRFLib_uniform());
		mean[i] = GMRFLib_uniform();
	}

	GMRFLib_init_problem(&problem, x, b, c, mean, graph, my_pardiso_test_Q, (void *) graph, constr);
	GMRFLib_evaluate(problem);
	GMRFLib_Qinv(problem);

	for (i = 0; i < n; i++) {
		var = GMRFLib_Qinv_get(problem, i, i);
		printf("Qinv[%1d,%1d] = %g\n", i, i, *var);
	}
	GMRFLib_free_problem(problem);

	return 0;
}

int my_pardiso_test3(void)
{
	int id = GMRFLib_thread_id;
	int err = 0;

	FIXME("this is test3");
	if (1) {
		err = GMRFLib_pardiso_check_install(1, 0);
		if (err == GMRFLib_SUCCESS) {
			printf("PARDISO OK\n");
		} else {
			printf("PARDISO FAIL\n");
		}
	}

	GMRFLib_tabulate_Qfunc_tp *Qtab;
	GMRFLib_graph_tp *g;

	GMRFLib_openmp->strategy = GMRFLib_OPENMP_STRATEGY_PARDISO;
	GMRFLib_openmp_implement_strategy(GMRFLib_OPENMP_PLACES_OPTIMIZE, NULL, NULL);

	// GMRFLib_tabulate_Qfunc_from_file(&Qtab, &g, "Qdense.txt", -1, NULL);
	// GMRFLib_tabulate_Qfunc_from_file(&Qtab, &g, "I5.txt", -1, NULL);
	GMRFLib_tabulate_Qfunc_from_file(&Qtab, &g, "Q.txt", -1, NULL);

	GMRFLib_csr_tp *csr;
	GMRFLib_Q2csr(&csr, g, Qtab->Qfunc, Qtab->Qfunc_arg);
	GMRFLib_csr_print(stdout, csr);
	P(csr->n);

	// GMRFLib_printf_graph(stdout, g);
	GMRFLib_pardiso_store_tp *store = NULL;
	GMRFLib_pardiso_init(&store);
	GMRFLib_pardiso_reorder(store, g);
	GMRFLib_pardiso_build(store, g, Qtab->Qfunc, Qtab->Qfunc_arg);
	GMRFLib_pardiso_chol(store);

	GMRFLib_pardiso_Qinv(store);

	P(GMRFLib_openmp->max_threads_outer);
	P(GMRFLib_openmp->max_threads_inner);

	int nrhs = S.nrhs_max, k;
#pragma omp parallel for private(k) num_threads(GMRFLib_openmp->max_threads_outer)
	for (k = 0; k < 1000; k++) {
		GMRFLib_thread_id = id;

		int i;
		GMRFLib_pardiso_store_tp *local_store = NULL;

		printf("this is k= %d from thread %d\n", k, omp_get_thread_num());
		GMRFLib_duplicate_pardiso_store(&local_store, store, -1, 1);

		int view = 1;
		double *x = Calloc(g->n * nrhs, double);
		double *b = Calloc(g->n * nrhs, double);

		for (i = 0; i < g->n * nrhs; i++) {
			b[i] = ISQR(i + 1);
		}
		GMRFLib_pardiso_solve_LLT(local_store, x, b, nrhs);
		if (view) {
			printf("solve LLT\n");
			for (int ii = 0; ii < nrhs; ii++) {
				for (i = 0; i < g->n; i++) {
					printf("ii %d i %d  x= %f  b=%f\n", ii, i, x[i + ii * g->n], b[i + ii * g->n]);
				}
			}
		}

		if (1) {
			for (i = 0; i < g->n; i++) {
				b[i] = ISQR(i + 1);
				x[i] = 0.0;
			}

			double *x2 = Calloc(g->n, double);
			GMRFLib_pardiso_solve_L(local_store, x2, b, 1);
			GMRFLib_pardiso_solve_LT(local_store, x, x2, 1);
			if (view) {
				printf("solve L D LT\n");
				for (int ii = 0; ii < nrhs; ii++) {
					for (i = 0; i < g->n; i++) {
						printf("ii %d i %d  x= %f  b=%f\n", ii, i, x[i + ii * g->n], b[i + ii * g->n]);
					}
				}
			}
		}

		for (i = 0; i < g->n * nrhs; i++) {
			b[i] = ISQR(i + 1);
		}
		GMRFLib_pardiso_solve_LT(local_store, x, b, 1);
		if (view) {
			printf("solve LT\n");
			for (int ii = 0; ii < nrhs; ii++) {
				for (i = 0; i < g->n; i++) {
					printf("ii %d i %d  x= %f  b=%f\n", ii, i, x[i + ii * g->n], b[i + ii * g->n]);
				}
			}
		}

		for (i = 0; i < g->n * nrhs; i++) {
			b[i] = ISQR(i + 1);
		}
		GMRFLib_pardiso_solve_L(local_store, x, b, 1);
		if (view) {
			printf("solve L\n");
			for (int ii = 0; ii < nrhs; ii++) {
				for (i = 0; i < g->n; i++) {
					printf("ii %d i %d  x= %f  b=%f\n", ii, i, x[i + ii * g->n], b[i + ii * g->n]);
				}
			}
		}

		GMRFLib_pardiso_free(&local_store);
		Free(x);
		Free(b);
	}

	GMRFLib_free_tabulate_Qfunc(Qtab);
	GMRFLib_pardiso_free(&store);
	GMRFLib_csr_free(&csr);
	GMRFLib_graph_free(g);

	return GMRFLib_SUCCESS;
}

int my_pardiso_test4(void)
{
	int k;
	GMRFLib_tabulate_Qfunc_tp *Qtab;
	GMRFLib_graph_tp *g;
	GMRFLib_csr_tp *csr;

	GMRFLib_smtp = GMRFLib_SMTP_PARDISO;
	GMRFLib_openmp->strategy = GMRFLib_OPENMP_STRATEGY_PARDISO;
	GMRFLib_openmp_implement_strategy(GMRFLib_OPENMP_PLACES_OPTIMIZE, NULL, NULL);
	P(GMRFLib_openmp->max_threads_inner);

	GMRFLib_tabulate_Qfunc_from_file(&Qtab, &g, "Q-problem-ijformat.txt", -1, NULL);
	GMRFLib_Q2csr(&csr, g, Qtab->Qfunc, Qtab->Qfunc_arg);

	// GMRFLib_csr_write("Q-problem-csr.dat", csr);
	// GMRFLib_csr_read("Q-problem-csr.dat", &csr);

	GMRFLib_pardiso_store_tp *store = NULL;
	GMRFLib_pardiso_init(&store);
	GMRFLib_pardiso_reorder(store, g);
	GMRFLib_pardiso_build(store, g, Qtab->Qfunc, Qtab->Qfunc_arg);

	for (k = 0; k < 1000; k++) {
		GMRFLib_pardiso_chol(store);
		printf("k %d logdet %.12f\n", k, GMRFLib_pardiso_logdet(store));
	}

	exit(0);
}

int my_pardiso_test5(void)
{
	int id = GMRFLib_thread_id;

	S.msglvl = 1;
	S.csr_check = 1;
	GMRFLib_openmp->strategy = GMRFLib_OPENMP_STRATEGY_PARDISO;
	GMRFLib_openmp_implement_strategy(GMRFLib_OPENMP_PLACES_DEFAULT, NULL, NULL);

	P(GMRFLib_openmp->max_threads_outer);
	P(GMRFLib_openmp->max_threads_inner);
	P(omp_get_nested());

	int k;
#pragma omp parallel for private(k) num_threads(GMRFLib_openmp->max_threads_outer)
	for (k = 0; k < 1000; k++) {
		GMRFLib_thread_id = id;

		// I do not free anything here...

		P(k);
		GMRFLib_tabulate_Qfunc_tp *Qtab = NULL;
		GMRFLib_graph_tp *g = NULL;

		GMRFLib_tabulate_Qfunc_from_file(&Qtab, &g, "Qsparse2.dat", -1, NULL);
		GMRFLib_csr_tp *csr = NULL;
		GMRFLib_Q2csr(&csr, g, Qtab->Qfunc, Qtab->Qfunc_arg);

		GMRFLib_pardiso_store_tp *store = NULL;
		GMRFLib_pardiso_init(&store);
		GMRFLib_pardiso_reorder(store, g);
		GMRFLib_pardiso_build(store, g, Qtab->Qfunc, Qtab->Qfunc_arg);
		GMRFLib_pardiso_chol(store);
		GMRFLib_pardiso_Qinv(store);
	}
	exit(0);
}

int my_pardiso_test6(GMRFLib_ai_store_tp * ai_store, GMRFLib_Qfunc_tp * Qfunc, void *Qfunc_arg, double *c)
{
	int id = GMRFLib_thread_id;
	int n = ai_store->problem->sub_graph->n;
	int i;

	assert(omp_get_num_threads() == 1);

	GMRFLib_thread_id = 0;
	GMRFLib_tabulate_Qfunc_tp *tab = NULL;
	GMRFLib_problem_tp *problem = ai_store->problem;
	GMRFLib_pardiso_store_tp *pardiso_store = problem->sub_sm_fact.PARDISO_fact;

	GMRFLib_tabulate_Qfunc(&tab, ai_store->problem->sub_graph, Qfunc, Qfunc_arg, NULL);
	P(GMRFLib_openmp->max_threads_outer);

#pragma omp parallel for private(i) num_threads(GMRFLib_openmp->max_threads_outer)
	for (i = 0; i < n; i++) {
		GMRFLib_thread_id = id;

		int *iparm = Calloc(GMRFLib_PARDISO_PLEN, int);
		double *dparm = Calloc(GMRFLib_PARDISO_PLEN, double);

		Memcpy(iparm, pardiso_store->iparm_default, GMRFLib_PARDISO_PLEN * sizeof(int));
		Memcpy(dparm, pardiso_store->dparm_default, GMRFLib_PARDISO_PLEN * sizeof(double));

		iparm[7] = 0;
		iparm[25] = 0;

		int j;
		int error = 0;
		int maxfct = 1;
		int mnum = 1;
		int msglvl = 0;
		int mtype = -2;
		int phase = 33;
		int idum = 0;
		int one = 1;

		double *work = Calloc(3 * n, double);
		double *b = work;
		double *x = work + n;
		double *res = work + 2 * n;
		double err;
		double fake_a = 0.0;
		int fake_ia1 = 0;
		int fake_ja1 = 0;

		// if I set b[10]=1, then it works in parallel, but not if they are different
		b[i] = 1.0;

		assert(GMRFLib_openmp->max_threads_inner == iparm[2]);
		omp_set_num_threads(iparm[2]);

		pardiso(pardiso_store->pt, &maxfct, &mnum, &mtype, &phase, &n,
			// not in use
			&fake_a, &fake_ia1, &fake_ja1,
			// 
			&idum, &one, iparm, &msglvl, b, x, &error, dparm);

		// res = Q x
		GMRFLib_Qx2(res, x, problem->sub_graph, tab->Qfunc, tab->Qfunc_arg, c);

		for (err = 0.0, j = 0; j < n; j++) {
			err += SQR(res[j] - b[j]);
		}
		err /= n;
		if (err > 1E-4)
			printf("i %d err %g\n", i, err);

		Free(work);
	}
	exit(0);
}

int my_pardiso_test7(void)
{
	int id = GMRFLib_thread_id;

	S.msglvl = 0;
	S.csr_check = 1;
	GMRFLib_openmp->strategy = GMRFLib_OPENMP_STRATEGY_PARDISO;
	GMRFLib_openmp_implement_strategy(GMRFLib_OPENMP_PLACES_DEFAULT, NULL, NULL);

	P(GMRFLib_openmp->max_threads_outer);
	P(GMRFLib_openmp->max_threads_inner);
	P(omp_get_nested());

	int k;
#pragma omp parallel for private(k) num_threads(GMRFLib_openmp->max_threads_outer)
	for (k = 0; k < 100; k++) {
		GMRFLib_thread_id = id;

		// I do not free anything here...

		omp_set_num_threads(GMRFLib_openmp->max_threads_inner);
		P(k);
		GMRFLib_tabulate_Qfunc_tp *Qtab = NULL;
		GMRFLib_graph_tp *g = NULL;

		if (k < 500)
			GMRFLib_tabulate_Qfunc_from_file(&Qtab, &g, "Q.txt", -1, NULL);
		else
			GMRFLib_tabulate_Qfunc_from_file(&Qtab, &g, "I5.txt", -1, NULL);
		GMRFLib_csr_tp *csr = NULL;
		GMRFLib_Q2csr(&csr, g, Qtab->Qfunc, Qtab->Qfunc_arg);

		GMRFLib_pardiso_store_tp *store = NULL;
		GMRFLib_pardiso_init(&store);
		GMRFLib_pardiso_reorder(store, g);
		GMRFLib_pardiso_build(store, g, Qtab->Qfunc, Qtab->Qfunc_arg);
		GMRFLib_pardiso_chol(store);
		GMRFLib_pardiso_Qinv(store);
	}
	exit(0);
}
