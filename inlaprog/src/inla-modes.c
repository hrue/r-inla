
/* inla-modes.c
 * 
 * Copyright (C) 2007-2024 Havard Rue
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or (at
 * your option) any later version.
 *  * This program is distributed in the hope that it will be useful, but
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

int inla_qinv(const char *filename, const char *constrfile, const char *outfile)
{
	/*
	 * Compute the marginal variances for Cij file in FILENAME and output on stdout, the marginal variances 
	 */
	int i, j, jj, k;

	GMRFLib_tabulate_Qfunc_tp *tab;
	GMRFLib_graph_tp *graph;
	GMRFLib_problem_tp *problem = NULL;
	GMRFLib_constr_tp *constr = NULL;
	GMRFLib_matrix_tp *constr_x = NULL;
	FILE *fp;

	GMRFLib_tabulate_Qfunc_from_file(&tab, &graph, filename, -1, NULL);
	fp = fopen(constrfile, "r");
	if (fp) {
		fclose(fp);
		constr_x = GMRFLib_read_fmesher_file(constrfile, 0L, SEEK_CUR);
		if (constr_x->A[0] > 0) {
			constr = Calloc(1, GMRFLib_constr_tp);
			constr->nc = (int) constr_x->A[0];
			constr->a_matrix = &constr_x->A[1];
			constr->e_vector = &constr_x->A[constr->nc * graph->n + 1];
			GMRFLib_prepare_constr(constr, graph, 1);
		}
	}

	if (GMRFLib_smtp == GMRFLib_SMTP_PARDISO) {
		GMRFLib_reorder = GMRFLib_REORDER_PARDISO;
		GMRFLib_openmp->strategy = GMRFLib_OPENMP_STRATEGY_PARDISO;
		GMRFLib_openmp_implement_strategy(GMRFLib_OPENMP_PLACES_DEFAULT, NULL, NULL);
	} else if (GMRFLib_smtp == GMRFLib_SMTP_BAND) {
		GMRFLib_reorder = GMRFLib_REORDER_BAND;
	} else {
		GMRFLib_reorder = GMRFLib_REORDER_DEFAULT;
		GMRFLib_optimize_reorder(graph, NULL, NULL, NULL);
	}
	int thread_id = 0;
	assert(omp_get_thread_num() == 0);
	GMRFLib_init_problem(thread_id, &problem, NULL, NULL, NULL, NULL, graph, tab->Qfunc, tab->Qfunc_arg, constr);
	GMRFLib_Qinv(problem);

	/*
	 * write a fmesher file and just pass the filename 
	 */
	GMRFLib_matrix_tp *M = Calloc(1, GMRFLib_matrix_tp);

	M->nrow = graph->n;
	M->ncol = graph->n;
	M->elems = 0;
	for (i = 0; i < graph->n; i++) {
		M->elems += 1 + graph->nnbs[i];
	}

	M->i = Calloc(M->elems, int);
	M->j = Calloc(M->elems, int);
	M->values = Calloc(M->elems, double);

	k = 0;
	for (i = 0; i < graph->n; i++) {
		M->i[k] = i;
		M->j[k] = i;
		M->values[k] = *GMRFLib_Qinv_get(problem, i, i);
		k++;

		for (jj = 0; jj < graph->nnbs[i]; jj++) {
			j = graph->nbs[i][jj];
			M->i[k] = i;
			M->j[k] = j;
			M->values[k] = *GMRFLib_Qinv_get(problem, i, j);
			k++;
		}
	}
	assert(k == M->elems);

	GMRFLib_write_fmesher_file(M, outfile, (long int) 0, -1);

	return 0;
}

int inla_qsolve(const char *Qfilename, const char *Afilename, const char *Bfilename, const char *method)
{
	/*
	 * Solve Q X = B, L^T X = B, or L X = B
	 */

	GMRFLib_tabulate_Qfunc_tp *tab;
	GMRFLib_graph_tp *graph;
	GMRFLib_problem_tp *problem = NULL;

	/*
	 * I need B to be dense 
	 */
	GMRFLib_matrix_tp *B = GMRFLib_read_fmesher_file(Bfilename, (long int) 0, -1);
	assert(B->i == NULL);				       /* I want B as dense matrix */

	GMRFLib_tabulate_Qfunc_from_file(&tab, &graph, Qfilename, -1, NULL);
	if (GMRFLib_smtp == GMRFLib_SMTP_PARDISO) {
		GMRFLib_reorder = GMRFLib_REORDER_PARDISO;
		GMRFLib_pardiso_set_nrhs(IMIN(GMRFLib_MAX_THREADS(), B->ncol));
		GMRFLib_openmp_implement_strategy(GMRFLib_OPENMP_PLACES_DEFAULT, NULL, NULL);
	} else if (GMRFLib_smtp == GMRFLib_SMTP_BAND) {
		GMRFLib_reorder = GMRFLib_REORDER_BAND;
	} else if (GMRFLib_smtp == GMRFLib_SMTP_TAUCS) {
		if (GMRFLib_reorder == GMRFLib_REORDER_DEFAULT) {
			GMRFLib_optimize_reorder(graph, NULL, NULL, NULL);
		}
	} else {
		assert(0 == 1);
	}

	int thread_id = 0;
	assert(omp_get_thread_num() == 0);
	GMRFLib_init_problem(thread_id, &problem, NULL, NULL, NULL, NULL, graph, tab->Qfunc, tab->Qfunc_arg, NULL);
	assert(problem->n == B->nrow);

	if (!strcasecmp(method, "solve")) {
		GMRFLib_solve_llt_sparse_matrix(B->A, B->ncol, &(problem->sub_sm_fact), problem->sub_graph);
	} else if (!strcasecmp(method, "forward")) {
		assert(GMRFLib_smtp != GMRFLib_SMTP_PARDISO);
		GMRFLib_solve_l_sparse_matrix(B->A, B->ncol, &(problem->sub_sm_fact), problem->sub_graph);
	} else if (!strcasecmp(method, "backward")) {
		assert(GMRFLib_smtp != GMRFLib_SMTP_PARDISO);
		GMRFLib_solve_lt_sparse_matrix(B->A, B->ncol, &(problem->sub_sm_fact), problem->sub_graph);
	} else {
		assert(0 == 1);
	}

	B->iA = NULL;
	GMRFLib_write_fmesher_file(B, Afilename, (long int) 0, -1);

	return 0;
}

int inla_qsample(const char *filename, const char *outfile, const char *nsamples, const char *rngfile,
		 const char *samplefile, const char *bfile, const char *mufile, const char *constrfile,
		 const char *meanfile, const char *selectionfile, int verbose)
{
	int output_every = 100;
	double t_ref = GMRFLib_timer(), t_reff = GMRFLib_timer();
	size_t siz, ret;
	char *state;
	FILE *fp;

	if (verbose) {
		fprintf(stderr, "inla_qsample: start pre...\n");
	}
	fp = fopen(rngfile, "rb");
	if (fp) {
		fseek(fp, 0L, SEEK_END);
		siz = ftell(fp) + 1;
		rewind(fp);
		state = Calloc(siz, char);
		ret = fread((void *) state, (size_t) 1, siz, fp);
		if (ret > 0) {
			GMRFLib_uniform_setstate((void *) state);
		}
		fclose(fp);
		Free(state);
	}

	int i, ns = 0;
	GMRFLib_tabulate_Qfunc_tp *tab;
	GMRFLib_graph_tp *graph;
	GMRFLib_problem_tp *problem = NULL;

	GMRFLib_matrix_tp *M = Calloc(1, GMRFLib_matrix_tp), *S = NULL, *b = NULL, *mu = NULL, *constr_x = NULL, *selection = NULL;
	GMRFLib_constr_tp *constr = NULL;

	inla_sread_ints(&ns, 1, nsamples);
	GMRFLib_tabulate_Qfunc_from_file(&tab, &graph, filename, -1, NULL);

	fp = fopen(samplefile, "r");
	if (fp) {
		fclose(fp);				       /* file exists */
		S = GMRFLib_read_fmesher_file(samplefile, 0L, SEEK_CUR);
	}

	fp = fopen(bfile, "r");
	if (fp) {
		fclose(fp);				       /* file exists */
		b = GMRFLib_read_fmesher_file(bfile, 0L, SEEK_CUR);
	}

	fp = fopen(mufile, "r");
	if (fp) {
		fclose(fp);				       /* file exists */
		mu = GMRFLib_read_fmesher_file(mufile, 0L, SEEK_CUR);
	}

	fp = fopen(selectionfile, "r");
	if (fp) {
		fclose(fp);				       /* file exists */
		selection = GMRFLib_read_fmesher_file(selectionfile, 0L, SEEK_CUR);
	}

	fp = fopen(constrfile, "r");
	if (fp) {
		fclose(fp);
		constr_x = GMRFLib_read_fmesher_file(constrfile, 0L, SEEK_CUR);
		if (constr_x->A[0] > 0) {
			constr = Calloc(1, GMRFLib_constr_tp);
			constr->nc = (int) constr_x->A[0];
			constr->a_matrix = &constr_x->A[1];
			constr->e_vector = &constr_x->A[constr->nc * graph->n + 1];
			GMRFLib_prepare_constr(constr, graph, 1);
		}
	}

	if (GMRFLib_smtp == GMRFLib_SMTP_PARDISO) {
		GMRFLib_reorder = GMRFLib_REORDER_PARDISO;
		GMRFLib_openmp->strategy = GMRFLib_OPENMP_STRATEGY_PARDISO;
		GMRFLib_openmp_implement_strategy(GMRFLib_OPENMP_PLACES_DEFAULT, NULL, NULL);
	} else if (GMRFLib_smtp == GMRFLib_SMTP_BAND) {
		GMRFLib_reorder = GMRFLib_REORDER_BAND;
	} else {
		GMRFLib_reorder = GMRFLib_REORDER_DEFAULT;
		GMRFLib_optimize_reorder(graph, NULL, NULL, NULL);
	}

	if (verbose) {
		fprintf(stderr, "inla_qsample: end pre %.2fs\n", GMRFLib_timer() - t_ref);
	}
	t_ref = GMRFLib_timer();
	if (verbose) {
		fprintf(stderr, "inla_qsample: start prepare the model...\n");
	}

	int thread_id = 0;
	assert(omp_get_thread_num() == 0);
	GMRFLib_init_problem(thread_id, &problem, NULL, (b ? b->A : NULL), NULL, (mu ? mu->A : NULL), graph, tab->Qfunc, tab->Qfunc_arg, constr);

	if (verbose) {
		fprintf(stderr, "inla_qsample: end prepare the model %.2fs\n", GMRFLib_timer() - t_ref);
	}
	t_ref = GMRFLib_timer();

	if (selection) {
		M->nrow = selection->nrow + 1;
	} else {
		M->nrow = graph->n + 1;
	}
	M->ncol = ns;
	M->elems = M->ncol * M->nrow;
	M->A = Calloc(M->nrow * M->ncol, double);

	if (verbose) {
		fprintf(stderr, "inla_qsample: start to sample %1d samples...\n", ns);
	}

	if ((GMRFLib_smtp == GMRFLib_SMTP_PARDISO)) {
		for (i = 0; i < ns; i++) {
			if (!S) {
				GMRFLib_sample(problem);
			} else {
				Memcpy(problem->sample, &(S->A[i * S->nrow]), S->nrow * sizeof(double));
			}
			GMRFLib_evaluate(problem);

			if (!selection) {
				Memcpy(&(M->A[i * M->nrow]), problem->sample, graph->n * sizeof(double));
			} else {
				for (int ii = 0; ii < selection->nrow; ii++) {
					M->A[i * M->nrow + ii] = problem->sample[(int) selection->A[ii]];
				}
			}
			M->A[(i + 1) * M->nrow - 1] = problem->sub_logdens;

			if (verbose && (!((i + 1) % output_every) || i == (ns - 1))) {
				fprintf(stderr, "inla_qsample: done with %1d samples, with %.2f samples/s and %.2fs in total\n",
					i + 1, (i + 1.0) / (GMRFLib_timer() - t_ref), (GMRFLib_timer() - t_ref));
			}
		}
	} else {
		GMRFLib_problem_tp **problems = Calloc(GMRFLib_openmp->max_threads_outer, GMRFLib_problem_tp *);
#pragma omp parallel for private(i) num_threads(GMRFLib_openmp->max_threads_outer)
		for (i = 0; i < ns; i++) {
			int thread = omp_get_thread_num();
			if (problems[thread] == NULL) {
				problems[thread] = GMRFLib_duplicate_problem(problem, 0, 1, 1);
			}
			if (!S) {
				GMRFLib_sample(problems[thread]);
			} else {
				Memcpy(problems[thread]->sample, &(S->A[i * S->nrow]), S->nrow * sizeof(double));
			}
			GMRFLib_evaluate(problems[thread]);

			if (!selection) {
				Memcpy(&(M->A[i * M->nrow]), problems[thread]->sample, graph->n * sizeof(double));
			} else {
				for (int ii = 0; ii < selection->nrow; ii++) {
					M->A[i * M->nrow + ii] = problems[thread]->sample[(int) selection->A[ii]];
				}
			}
			M->A[(i + 1) * M->nrow - 1] = problems[thread]->sub_logdens;
		}
	}

	if (verbose) {
		fprintf(stderr, "inla_qsample: end in %.2fs with %.2f samples/s\n", GMRFLib_timer() - t_ref, (double) ns / (GMRFLib_timer() - t_ref));
	}
	t_ref = GMRFLib_timer();
	if (verbose) {
		fprintf(stderr, "inla_qsample: start post...\n");
	}

	GMRFLib_write_fmesher_file(M, outfile, (long int) 0, -1);

	GMRFLib_matrix_tp *CM = Calloc(1, GMRFLib_matrix_tp);
	CM->nrow = M->nrow - 1;
	CM->ncol = 1;
	CM->elems = CM->ncol * CM->nrow;
	CM->A = Calloc(CM->nrow * CM->ncol, double);
	if (!selection) {
		Memcpy(CM->A, problem->mean_constr, graph->n * sizeof(double));
	} else {
		for (int ii = 0; ii < selection->nrow; ii++) {
			CM->A[ii] = problem->mean_constr[(int) selection->A[ii]];
		}
	}
	GMRFLib_write_fmesher_file(CM, meanfile, (long int) 0, -1);

	state = (char *) GMRFLib_rng_getstate(&siz);
	fp = fopen(rngfile, "wb");
	fwrite((void *) state, (size_t) 1, siz, fp);
	fclose(fp);

	if (verbose) {
		fprintf(stderr, "inla_qsample: end post %.2fs\n", GMRFLib_timer() - t_ref);
		fprintf(stderr, "inla_qsample: total time %.2fs\n", GMRFLib_timer() - t_reff);
	}

	return 0;
}

int inla_finn(const char *UNUSED(filename))
{
	return 0;
}

int inla_qreordering(const char *filename)
{
	/*
	 * return the rordering either given or computed
	 */
	int i;
	GMRFLib_graph_tp *graph;

	if (GMRFLib_is_fmesher_file(filename, (long int) 0, -1) == GMRFLib_SUCCESS) {
		GMRFLib_tabulate_Qfunc_tp *qtab = NULL;
		GMRFLib_tabulate_Qfunc_from_file(&qtab, &graph, filename, -1, NULL);
		GMRFLib_free_tabulate_Qfunc(qtab);
	} else {
		GMRFLib_graph_read(&graph, filename);
	}

	if (G.reorder < 0) {
		GMRFLib_optimize_reorder(graph, NULL, NULL, NULL);
	}
	GMRFLib_sm_fact_tp sm_fact;
	sm_fact.smtp = GMRFLib_SMTP_TAUCS;
	GMRFLib_compute_reordering(&sm_fact, graph, NULL);

	printf("QREORDERING\n");			       /* code used when the output is parsed */
	printf("%s\n", GMRFLib_reorder_name(GMRFLib_reorder));
	printf("%1d\n", GMRFLib_reorder);
	for (i = 0; i < graph->n; i++) {
		printf("%1d\n", sm_fact.remap[i]);
	}

	return 0;
}

int inla_fgn(char *infile, char *outfile)
{
	double H, H_intern, *res;
	int i, k, len, K, nH;

	GMRFLib_matrix_tp *Hm = GMRFLib_read_fmesher_file(infile, 0, -1);
	assert(Hm->ncol == 1);
	nH = Hm->nrow - 1;
	assert(nH >= 1);
	K = (int) GMRFLib_matrix_get(0, 0, Hm);		       // first element is K, then H's.
	len = 2 * K + 1;
	res = Calloc(nH * len, double);
	for (i = k = 0; i < nH; i++, k += len) {
		H = res[k] = GMRFLib_matrix_get(i + 1, 0, Hm);
		H_intern = map_H(H, MAP_BACKWARD, NULL);
		inla_fgn_get(&res[k + 1], &res[k + 1 + K], H_intern, K);
	}
	GMRFLib_matrix_tp *M = Calloc(1, GMRFLib_matrix_tp), *M_t;
	M->ncol = nH;
	M->nrow = len;
	M->elems = M->nrow * M->ncol;
	M->A = res;
	M_t = GMRFLib_matrix_transpose(M);
	GMRFLib_write_fmesher_file(M_t, outfile, 0L, -1);
	GMRFLib_matrix_free(M);
	GMRFLib_matrix_free(M_t);

	return GMRFLib_SUCCESS;
}
