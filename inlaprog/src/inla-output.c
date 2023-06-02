
/* inla-output.c
 * 
 * Copyright (C) 2007-2023 Havard Rue
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

int inla_output_Q(inla_tp *mb, const char *dir, GMRFLib_graph_tp *graph)
{
	GMRFLib_problem_tp *p = NULL;
	char *fnm = NULL, *newdir = NULL;
	FILE *fp = NULL;
	int thread_id = 0;
	assert(omp_get_thread_num() == 0);

	GMRFLib_init_problem(thread_id, &p, NULL, NULL, NULL, NULL, graph, GMRFLib_Qfunc_generic, (void *) graph, NULL);
	GMRFLib_sprintf(&newdir, "%s/Q", dir);
	if (mb->verbose) {
		printf("\t\tstore factorisation results in[%s]\n", newdir);
	}
	if (inla_mkdir(newdir) == INLA_OK) {
		if (mb->output->q) {
			if (mb->verbose) {
				printf("\t\tstore info precision and related matrices in[%s]\n", newdir);
			}
			GMRFLib_sprintf(&fnm, "%s/%s", newdir, "precision-matrix");
			GMRFLib_bitmap_problem((const char *) fnm, p);
			Free(fnm);
		}
		GMRFLib_sprintf(&fnm, "%s/%s", newdir, "factorisation-information.txt");
		fp = fopen(fnm, "w");
		if (fp) {
			GMRFLib_fact_info_report(fp, &(p->sub_sm_fact));
			fclose(fp);
		}
		GMRFLib_free_problem(p);
		Free(fnm);
	}
	Free(newdir);

	return INLA_OK;
}

int inla_output_graph(inla_tp *mb, const char *dir, GMRFLib_graph_tp *graph)
{
	char *fnm = NULL;
	GMRFLib_sprintf(&fnm, "%s/graph.dat", dir);
	if (mb->verbose) {
		printf("\t\tstore graph in[%s]\n", fnm);
	}
	GMRFLib_graph_write_b(fnm, graph);
	Free(fnm);

	return INLA_OK;
}

int inla_output_matrix(const char *dir, const char *sdir, const char *filename, int n, double *matrix, int *order)
{
	char *fnm, *ndir;

	if (sdir) {
		GMRFLib_sprintf(&ndir, "%s/%s", dir, sdir);
	} else {
		GMRFLib_sprintf(&ndir, "%s", dir);
	}

	GMRFLib_sprintf(&fnm, "%s/%s", ndir, filename);

	GMRFLib_matrix_tp *M = Calloc(1, GMRFLib_matrix_tp);

	M->nrow = M->ncol = n;
	M->elems = ISQR(n);
	M->A = Calloc(ISQR(n), double);
	if (order == NULL) {
		Memcpy(M->A, matrix, ISQR(n) * sizeof(double));
	} else {
		int i, j;
		for (i = 0; i < n; i++) {
			for (j = 0; j < n; j++) {
				M->A[order[i] + order[j] * n] = matrix[i + j * n];
			}
		}
	}

	M->offset = 0L;
	M->whence = SEEK_SET;
	M->tell = -1L;

	GMRFLib_write_fmesher_file(M, fnm, M->offset, M->whence);
	GMRFLib_matrix_free(M);

	Free(fnm);
	Free(ndir);

	return INLA_OK;
}

int inla_output_names(const char *dir, const char *sdir, int n, const char **names, const char *suffix)
{
	FILE *fp;
	char *fnm, *ndir;

	GMRFLib_sprintf(&ndir, "%s/%s", dir, sdir);
	GMRFLib_sprintf(&fnm, "%s/NAMES", ndir);

	int i;
	fp = fopen(fnm, "w");
	for (i = 0; i < n; i++) {
		fprintf(fp, "%s%s\n", names[i], (suffix ? suffix : ""));
	}
	fclose(fp);

	Free(fnm);
	Free(ndir);

	return INLA_OK;
}

int inla_output_size(const char *dir, const char *sdir, int n, int N, int Ntotal, int ngroup, int nrep)
{
	char *fnm, *ndir;
	GMRFLib_sprintf(&ndir, "%s/%s", dir, sdir);
	GMRFLib_sprintf(&fnm, "%s/size.dat", ndir);

	Dinit_s(fnm);
	D5W(n, (N > 0 ? N : n), (Ntotal > 0 ? Ntotal : n), (ngroup > 0 ? ngroup : 1), (nrep > 0 ? nrep : 1));
	Dclose();
	Free(fnm);
	Free(ndir);

	return INLA_OK;
}

char *inla_create_hyperid(int id, const char *secname)
{
	char *hyperid = NULL;
	GMRFLib_sprintf(&hyperid, "%1d|%s", id, secname);

	return (hyperid);
}

int inla_output_hyperid(const char *dir, const char *sdir, char *hyperid)
{
	FILE *fp;
	char *fnm, *ndir;

	// fprintf(stderr, "output hyperid %s / %s [%s]\n", dir, sdir, hyperid);

	GMRFLib_sprintf(&ndir, "%s/%s", dir, sdir);
	GMRFLib_sprintf(&fnm, "%s/hyperid.dat", ndir);

	fp = fopen(fnm, "w");
	if (!fp) {
		inla_error_open_file(fnm);
	}
	fprintf(fp, "%s\n", (hyperid ? hyperid : ""));
	fclose(fp);

	Free(fnm);
	Free(ndir);

	return INLA_OK;
}

int inla_output_id_names(const char *dir, const char *sdir, inla_file_contents_tp *fc)
{
	if (!fc) {
		return INLA_OK;
	}

	char *fnm, *ndir;

	GMRFLib_sprintf(&ndir, "%s/%s", dir, sdir);
	GMRFLib_sprintf(&fnm, "%s/id-names.dat", ndir);

	inla_write_file_contents(fnm, fc);

	Free(fnm);
	Free(ndir);

	return INLA_OK;
}

int inla_output(inla_tp *mb)
{
	int n = 0, j, *offsets = NULL, len_offsets, local_verbose = 0;
	assert(mb);

	/*
	 * compute the offset for each pack of the results 
	 */
	len_offsets = 1 + mb->nf + mb->nlinear;
	offsets = Calloc(len_offsets, int);

	j = n = 0;
	offsets[j++] = n;
	n += mb->predictor_n + mb->predictor_m;
	for (int i = 0; i < mb->nf; i++) {
		offsets[j++] = n;
		n += mb->f_graph[i]->n;
	}
	for (int i = 0; i < mb->nlinear; i++) {
		offsets[j++] = n;
		n++;
	}
	if (GMRFLib_inla_mode == GMRFLib_MODE_TWOSTAGE_PART1 || GMRFLib_inla_mode == GMRFLib_MODE_COMPACT) {
		assert(mb->preopt->mnpred == mb->predictor_m + mb->predictor_n);
	} else {
		assert(mb->hgmrfm->graph->n == n);
	}

	assert(j == len_offsets);

	if (mb->verbose) {
		printf("Store results in directory[%s]\n", mb->dir);
	}
	/*
	 * turn off all the info about output-files... 
	 */
	local_verbose = 0;

	// do this here so they can be parallized in '_output_detail'
	if (1) {
		int offset = offsets[0];
		inla_output_detail(mb->dir, &(mb->density[offset]),
				   NULL, mb->predictor_n + mb->predictor_m, 1,
				   mb->predictor_output, mb->predictor_dir, mb->output->return_marginals_predictor,
				   NULL, NULL, NULL, mb->predictor_tag, NULL, local_verbose);
		inla_output_size(mb->dir, mb->predictor_dir, mb->predictor_n, mb->predictor_n,
				 mb->predictor_n + mb->predictor_m, -1, (mb->predictor_m == 0 ? 1 : 2));
	}

	if (mb->predictor_invlinkfunc && mb->predictor_user_scale) {
		char *sdir, *newtag;
		int offset = offsets[0];
		GMRFLib_sprintf(&newtag, "%s in user scale", mb->predictor_tag);
		GMRFLib_sprintf(&sdir, "%s-user-scale", mb->predictor_dir);
		inla_output_detail(mb->dir, &(mb->density[offset]),
				   NULL, mb->predictor_n + mb->predictor_m, 1,
				   mb->predictor_output, sdir, mb->output->return_marginals_predictor,
				   NULL, NULL, mb->transform_funcs, newtag, NULL, local_verbose);
		inla_output_size(mb->dir, sdir, mb->predictor_n + mb->predictor_m, -1, -1, -1, (mb->predictor_m == 0 ? 1 : 2));
	}

#pragma omp parallel for num_threads(IMIN(IMAX(1, mb->nf), GMRFLib_openmp->max_threads_outer))
	for (int ii = 0; ii < mb->nf; ii++) {
		int offset = offsets[ii + 1];
		inla_output_detail(mb->dir, &(mb->density[offset]),
				   mb->f_locations[ii],
				   mb->f_graph[ii]->n, mb->f_nrep[ii] * mb->f_ngroup[ii], mb->f_output[ii],
				   mb->f_dir[ii], mb->output->return_marginals,
				   NULL, NULL, NULL, mb->f_tag[ii], mb->f_modelname[ii], local_verbose);
		inla_output_size(mb->dir, mb->f_dir[ii], mb->f_n[ii], mb->f_N[ii], mb->f_Ntotal[ii], mb->f_ngroup[ii], mb->f_nrep[ii]);
		inla_output_id_names(mb->dir, mb->f_dir[ii], mb->f_id_names[ii]);
	}


#pragma omp parallel for num_threads(GMRFLib_MAX_THREADS())
	for (int k = 3; k < 9; k++) {
		int ii;
		if (k == 3) {
			char *fnm;
			GMRFLib_sprintf(&fnm, "%s/totaloffset", mb->dir);
			inla_mkdir(fnm);
			Free(fnm);
			GMRFLib_sprintf(&fnm, "%s/totaloffset/totaloffset.dat", mb->dir);
			Dinit(fnm);
			Free(fnm);
			for (ii = 0; ii < mb->predictor_n + mb->predictor_m; ii++) {
				D1W(OFFSET3(ii));
			}
			Dclose();
		} else if (k == 4) {
			for (ii = 0; ii < mb->nlinear; ii++) {
				int offset = offsets[mb->nf + 1 + ii];
				inla_output_detail(mb->dir, &(mb->density[offset]),
						   NULL, 1, 1, mb->linear_output[ii], mb->linear_dir[ii], mb->output->return_marginals,
						   NULL, NULL, NULL, mb->linear_tag[ii], NULL, local_verbose);
				inla_output_size(mb->dir, mb->linear_dir[ii], 1, -1, -1, -1, -1);
			}
		} else if (k == 5) {
			if (mb->density_lin) {
				char *newtag2, *newdir2;
				ii = 0;
				GMRFLib_sprintf(&newtag2, "lincombs.derived.all");
				GMRFLib_sprintf(&newdir2, "lincombs.derived.all");
				inla_output_detail(mb->dir, &(mb->density_lin[ii]), mb->lc_order, mb->nlc,
						   1, mb->lc_output[ii], newdir2, mb->output->return_marginals,
						   NULL, NULL, NULL, newtag2, NULL, local_verbose);
				inla_output_size(mb->dir, newdir2, mb->nlc, -1, -1, -1, -1);
				inla_output_names(mb->dir, newdir2, mb->nlc, (const char **) ((void *) mb->lc_tag), NULL);

				Free(newtag2);
				Free(newdir2);
			}
			if (mb->density_hyper) {
				for (ii = 0; ii < mb->ntheta; ii++) {
					char *sdir;
					GMRFLib_sprintf(&sdir, "hyperparameter-1-%.6d-%s", ii, mb->theta_dir[ii]);
					inla_output_detail(mb->dir, &(mb->density_hyper[ii]), NULL, 1, 1, mb->output, sdir,
							   mb->output->return_marginals, NULL, NULL, NULL, mb->theta_tag[ii], NULL, local_verbose);
					inla_output_hyperid(mb->dir, sdir, mb->theta_hyperid[ii]);
					inla_output_size(mb->dir, sdir, 1, -1, -1, -1, -1);

					GMRFLib_sprintf(&sdir, "hyperparameter-2-%.6d-%s-user-scale", ii, mb->theta_dir[ii]);
					inla_output_detail(mb->dir, &(mb->density_hyper[ii]), NULL, 1, 1, mb->output, sdir,
							   mb->output->return_marginals,
							   mb->theta_map[ii], mb->theta_map_arg[ii], NULL,
							   mb->theta_tag_userscale[ii], NULL, local_verbose);
					inla_output_hyperid(mb->dir, sdir, mb->theta_hyperid[ii]);
					inla_output_size(mb->dir, sdir, 1, -1, -1, -1, -1);
				}
			}
			if (GMRFLib_ai_INLA_userfunc0_density && GMRFLib_ai_INLA_userfunc0_dim > 0) {
				/*
				 * we need to create the corresponding normal as well 
				 */
				char *sdir;
				sdir = GMRFLib_strdup("random.effect.UserFunction0");
				inla_output_detail(mb->dir, GMRFLib_ai_INLA_userfunc0_density, NULL,
						   GMRFLib_ai_INLA_userfunc0_dim, 1, mb->output, sdir, mb->output->return_marginals,
						   NULL, NULL, NULL, "UserFunction0", NULL, local_verbose);
				inla_output_size(mb->dir, sdir, GMRFLib_ai_INLA_userfunc0_dim, -1, -1, -1, -1);
				Free(sdir);
			}
			if (GMRFLib_ai_INLA_userfunc1_density && GMRFLib_ai_INLA_userfunc1_dim > 0) {
				/*
				 * we need to create the corresponding normal as well 
				 */
				char *sdir;
				sdir = GMRFLib_strdup("random.effect.UserFunction1");
				inla_output_detail(mb->dir, GMRFLib_ai_INLA_userfunc1_density, NULL,
						   GMRFLib_ai_INLA_userfunc1_dim, 1, mb->output, sdir, mb->output->return_marginals,
						   NULL, NULL, NULL, "UserFunction1", NULL, local_verbose);
				inla_output_size(mb->dir, sdir, GMRFLib_ai_INLA_userfunc1_dim, -1, -1, -1, -1);

				Free(sdir);
			}
			if (GMRFLib_ai_INLA_userfunc2_density && GMRFLib_ai_INLA_userfunc2_n > 0) {
				for (ii = 0; ii < GMRFLib_ai_INLA_userfunc2_n; ii++) {
					/*
					 * we need to create the corresponding normal as well 
					 */
					char *sdir, *local_tag;

					int dim = GMRFLib_ai_INLA_userfunc2_len[ii];
					GMRFLib_sprintf(&sdir, "spde2.blc.%6.6d", ii + 1);
					GMRFLib_sprintf(&local_tag, "%s", GMRFLib_ai_INLA_userfunc2_tag[ii]);
					inla_output_detail(mb->dir, GMRFLib_ai_INLA_userfunc2_density[ii], NULL, dim, 1,
							   mb->output, sdir, mb->output->return_marginals,
							   NULL, NULL, NULL, local_tag, NULL, local_verbose);
					inla_output_size(mb->dir, sdir, dim, -1, -1, -1, -1);

					Free(sdir);
					Free(local_tag);
				}
			}
			if (GMRFLib_ai_INLA_userfunc3_density && GMRFLib_ai_INLA_userfunc3_n > 0) {
				for (ii = 0; ii < GMRFLib_ai_INLA_userfunc3_n; ii++) {
					/*
					 * we need to create the corresponding normal as well 
					 */
					char *sdir, *local_tag;

					int dim = GMRFLib_ai_INLA_userfunc3_len[ii];
					GMRFLib_sprintf(&sdir, "spde3.blc.%6.6d", ii + 1);
					GMRFLib_sprintf(&local_tag, "%s", GMRFLib_ai_INLA_userfunc3_tag[ii]);
					inla_output_detail(mb->dir, GMRFLib_ai_INLA_userfunc3_density[ii], NULL, dim, 1,
							   mb->output, sdir, mb->output->return_marginals,
							   NULL, NULL, NULL, local_tag, NULL, local_verbose);
					inla_output_size(mb->dir, sdir, dim, -1, -1, -1, -1);

					Free(sdir);
					Free(local_tag);
				}
			}
		} else if (k == 6) {
			if (mb->misc_output) {
				inla_output_misc(mb->dir, mb->misc_output, mb->ntheta, mb->theta_tag, mb->theta_from, mb->theta_to,
						 mb->lc_order, local_verbose, mb);
			}
		} else if (k == 7) {
			if (mb->gcpo) {
				inla_output_detail_gcpo(mb->dir, mb->gcpo, local_verbose);
			}
			if (mb->cpo) {
				inla_output_detail_cpo(mb->dir, mb->cpo, mb->predictor_ndata, local_verbose);
			}
			if (mb->po) {
				inla_output_detail_po(mb->dir, mb->po, mb->predictor_ndata, local_verbose);
			}
			if (mb->dic) {
				inla_output_detail_dic(mb->dir, mb->dic, mb->family_idx, mb->len_family_idx, local_verbose);
			}
			if (mb->output->mlik) {
				inla_output_detail_mlik(mb->dir, &(mb->mlik), local_verbose);
			}
		} else if (k == 8) {
			inla_output_detail_x(mb->dir, mb->x_file, mb->nx_file);
			inla_output_detail_theta(mb->dir, mb->theta, mb->ntheta);
			inla_output_gitid(mb->dir);
			inla_output_linkfunctions(mb->dir, mb);
			if (mb->output->q) {
				if (local_verbose == 0) {
					int save = mb->verbose;
					mb->verbose = 0;
					inla_output_Q(mb, mb->dir, mb->hgmrfm->graph);
					mb->verbose = save;
				}
			}
			if (mb->output->graph) {
				inla_output_graph(mb, mb->dir, mb->hgmrfm->graph);
			}
		}
	}

	int N = -1;
	if (GMRFLib_inla_mode == GMRFLib_MODE_TWOSTAGE_PART1 || GMRFLib_inla_mode == GMRFLib_MODE_COMPACT) {
		N = mb->preopt->n + mb->preopt->mnpred;
	} else {
		N = ((GMRFLib_hgmrfm_arg_tp *) mb->hgmrfm->Qfunc_arg)->N;
	}

#pragma omp parallel for num_threads(GMRFLib_openmp->max_threads_outer)
	for (int i = 0; i < 2; i++) {
		if (i == 0 && mb->density) {
			for (int ii = 0; ii < N; ii++) {
				GMRFLib_free_density(mb->density[ii]);
			}
			Free(mb->density);
		}
	}

	return INLA_OK;
}

int inla_output_detail_gcpo(const char *dir, GMRFLib_gcpo_tp *gcpo, int verbose)
{
	/*
	 * output whatever is requested.... 
	 */
	char *ndir = NULL, *msg = NULL, *nndir = NULL;
	int i, j, n;

	if (!gcpo) {
		return INLA_OK;
	}
	n = gcpo->n;

	GMRFLib_sprintf(&ndir, "%s/%s", dir, "gcpo");
	if (inla_mkdir(ndir) != 0) {
		GMRFLib_sprintf(&msg, "fail to create directory [%s]: %s", ndir, strerror(errno));
		inla_error_general(msg);
	}
	GMRFLib_sprintf(&nndir, "%s/%s", ndir, "gcpo.dat");
	{
		Dinit(nndir);
		if (verbose) {
#pragma omp critical (Name_2f68597610da063cce30cec1d19ac794527bc57e)
			{
				printf("\t\tstore gcpo-results in[%s]\n", nndir);
			}
		}
		D1W(n);
		for (i = 0; i < n; i++) {
			D1W(gcpo->value[i]);
		}
		for (i = 0; i < n; i++) {
			D1W(gcpo->kld[i]);
		}
		for (i = 0; i < n; i++) {
			D1W(gcpo->mean[i]);
		}
		for (i = 0; i < n; i++) {
			D1W(gcpo->sd[i]);
		}
		for (i = 0; i < n; i++) {
			D1W(gcpo->groups[i]->n);
			for (j = 0; j < gcpo->groups[i]->n; j++) {
				D1W(gcpo->groups[i]->idx[j] + 1);	/* back to R-style indexing */
			}
			for (j = 0; j < gcpo->groups[i]->n; j++) {
				D1W(gcpo->groups[i]->val[j]);
			}
		}
		Dclose();
	}

	Free(ndir);
	Free(nndir);

	return INLA_OK;
}

int inla_output_detail_cpo(const char *dir, GMRFLib_ai_cpo_tp *cpo, int predictor_n, int verbose)
{
	/*
	 * output whatever is requested.... 
	 */
	char *ndir = NULL, *msg = NULL, *nndir = NULL;
	int i, n, add_empty = 1;

	if (!cpo) {
		return INLA_OK;
	}
	// n = cpo->n;
	n = predictor_n;				       /* the CPO and PIT are at the first predictor_n */

	GMRFLib_sprintf(&ndir, "%s/%s", dir, "cpo");
	if (inla_mkdir(ndir) != 0) {
		GMRFLib_sprintf(&msg, "fail to create directory [%s]: %s", ndir, strerror(errno));
		inla_error_general(msg);
	}
	GMRFLib_sprintf(&nndir, "%s/%s", ndir, "cpo.dat");
	{
		Dinit(nndir);
		if (verbose) {
#pragma omp critical (Name_6931944003ec5b67e4101e7995e2a5baab470602)
			{
				printf("\t\tstore cpo-results in[%s]\n", nndir);
			}
		}
		D1W(predictor_n);
		for (i = 0; i < n; i++) {
			if (cpo->value[i]) {
				D2W(i, cpo->value[i][0]);
			} else {
				if (add_empty) {
					D2W(i, NAN);
				}
			}
		}
		Dclose();
	}

	if (cpo->pit_value) {
		GMRFLib_sprintf(&nndir, "%s/%s", ndir, "pit.dat");
		Dinit(nndir);
		if (verbose) {
#pragma omp critical (Name_dbf18c17124ff85ed3fce2d13d88ad663accb843)
			{
				printf("\t\tstore pit-results in[%s]\n", nndir);
			}
		}
		D1W(predictor_n);
		for (i = 0; i < n; i++) {
			if (cpo->pit_value[i]) {
				D2W(i, cpo->pit_value[i][0]);
			} else {
				if (add_empty) {
					D2W(i, NAN);
				}
			}
		}
		Dclose();
	}
	if (cpo->failure) {
		GMRFLib_sprintf(&nndir, "%s/%s", ndir, "failure.dat");
		Dinit(nndir);
		if (verbose) {
#pragma omp critical (Name_f742aaba395a3aa285064964ce668d0be3e5e7ef)
			{
				printf("\t\tstore failure-results in[%s]\n", nndir);
			}
		}
		D1W(predictor_n);
		for (i = 0; i < n; i++) {
			if (cpo->failure[i]) {
				D2W(i, cpo->failure[i][0]);
			} else {
				if (add_empty) {
					D2W(i, NAN);
				}
			}
		}
		Dclose();
	}
	GMRFLib_sprintf(&nndir, "%s/%s", ndir, "summary.dat");
	{
		Dinit(nndir);
		if (verbose) {
#pragma omp critical (Name_dfd5be6da1007dda6dbb7abfb479ffb6fdd5c005)
			{
				printf("\t\tstore summary of cpo-results in[%s]\n", nndir);
			}
		}
		D2W(cpo->mean_value, cpo->gmean_value);
		Dclose();
	}

	Free(ndir);
	Free(nndir);

	return INLA_OK;
}

int inla_output_detail_po(const char *dir, GMRFLib_ai_po_tp *po, int predictor_n, int verbose)
{
	/*
	 * output whatever is requested.... 
	 */
	char *ndir = NULL, *msg = NULL, *nndir = NULL;
	int i, n, add_empty = 1;

	if (!po) {
		return INLA_OK;
	}
	n = predictor_n;				       /* the PO are at the first predictor_n */

	GMRFLib_sprintf(&ndir, "%s/%s", dir, "po");
	if (inla_mkdir(ndir) != 0) {
		GMRFLib_sprintf(&msg, "fail to create directory [%s]: %s", ndir, strerror(errno));
		inla_error_general(msg);
	}
	GMRFLib_sprintf(&nndir, "%s/%s", ndir, "po.dat");
	Dinit(nndir);
	if (verbose) {
#pragma omp critical (Name_879176a6e21c1dedb4becf90d4f21dd37210e46f)
		{
			printf("\t\tstore po-results in[%s]\n", nndir);
		}
	}
	D1W(predictor_n);
	for (i = 0; i < n; i++) {
		if (po->value[i]) {
			D3W(i, po->value[i][0], po->value[i][1]);
		} else {
			if (add_empty) {
				D3W(i, NAN, NAN);
			}
		}
	}
	Dclose();
	Free(ndir);
	Free(nndir);

	return INLA_OK;
}

int inla_output_detail_dic(const char *dir, GMRFLib_ai_dic_tp *dic, double *family_idx, int len_family_idx, int verbose)
{
	/*
	 * output whatever is requested.... 
	 */
	char *ndir = NULL, *msg = NULL, *nndir = NULL;
	double *tmp = NULL;

#define _PAD_WITH_NA(xx)						\
	if (1) {							\
		if (!tmp) tmp = Calloc(IMAX(dic->n_deviance, len_family_idx), double); \
		Memcpy(tmp, xx, dic->n_deviance*sizeof(double));	\
		int i;							\
		for(i = dic->n_deviance; i < len_family_idx; i++){	\
			tmp[i] = NAN;					\
		}							\
	}

	if (!dic) {
		return INLA_OK;
	}
	GMRFLib_sprintf(&ndir, "%s/%s", dir, "dic");
	if (inla_mkdir(ndir) != 0) {
		GMRFLib_sprintf(&msg, "fail to create directory [%s]: %s", ndir, strerror(errno));
		inla_error_general(msg);
	}
	GMRFLib_sprintf(&nndir, "%s/%s", ndir, "dic.dat");
	Dinit(nndir);
	if (verbose) {
#pragma omp critical (Name_19769e5d041b18c2ff6bee313efb503ad5ad5aef)
		{
			printf("\t\tstore dic-results in[%s]\n", nndir);
		}
	}
	D4W(dic->mean_of_deviance, dic->deviance_of_mean, dic->p, dic->dic);
	D4W(dic->mean_of_deviance_sat, dic->deviance_of_mean_sat, dic->p, dic->dic_sat);
	Dclose();

	if (dic->n_deviance > 0) {
		GMRFLib_matrix_tp *M = NULL;

		M = Calloc(1, GMRFLib_matrix_tp);
		M->nrow = len_family_idx;
		M->ncol = 1;
		M->elems = M->nrow * M->ncol;

		_PAD_WITH_NA(dic->sign);
		M->A = tmp;
		GMRFLib_sprintf(&nndir, "%s/%s", ndir, "sign.dat");
		GMRFLib_write_fmesher_file(M, nndir, (long int) 0, -1);

		_PAD_WITH_NA(dic->e_deviance);
		M->A = tmp;
		GMRFLib_sprintf(&nndir, "%s/%s", ndir, "e_deviance.dat");
		GMRFLib_write_fmesher_file(M, nndir, (long int) 0, -1);

		_PAD_WITH_NA(dic->e_deviance_sat);
		M->A = tmp;
		GMRFLib_sprintf(&nndir, "%s/%s", ndir, "e_deviance_sat.dat");
		GMRFLib_write_fmesher_file(M, nndir, (long int) 0, -1);

		_PAD_WITH_NA(dic->deviance_e);
		M->A = tmp;
		GMRFLib_sprintf(&nndir, "%s/%s", ndir, "deviance_e.dat");
		GMRFLib_write_fmesher_file(M, nndir, (long int) 0, -1);

		_PAD_WITH_NA(dic->deviance_e_sat);
		M->A = tmp;
		GMRFLib_sprintf(&nndir, "%s/%s", ndir, "deviance_e_sat.dat");
		GMRFLib_write_fmesher_file(M, nndir, (long int) 0, -1);

		M->A = family_idx;
		GMRFLib_sprintf(&nndir, "%s/%s", ndir, "family_idx.dat");
		GMRFLib_write_fmesher_file(M, nndir, (long int) 0, -1);

		Free(tmp);
		M->A = NULL;
		GMRFLib_matrix_free(M);
	}

	Free(tmp);
	Free(ndir);
	Free(nndir);
#undef _PAD_WITH_NA

	return INLA_OK;
}

int inla_output_misc(const char *dir, GMRFLib_ai_misc_output_tp *mo, int ntheta, char **theta_tag, char **theta_from,
		     char **theta_to, double *lc_order, int verbose, inla_tp *mb)
{
	/*
	 * output whatever is requested.... 
	 */
	char *ndir = NULL, *msg = NULL, *nndir = NULL, *nnndir = NULL;
	int i, any;

	if (!mo) {
		return INLA_OK;
	}
	GMRFLib_sprintf(&ndir, "%s/%s", dir, "misc");
	if (inla_mkdir(ndir) != 0) {
		GMRFLib_sprintf(&msg, "fail to create directory [%s]: %s", ndir, strerror(errno));
		inla_error_general(msg);
	}

	if (verbose) {
#pragma omp critical (Name_f2b4e6c2e966340fa5f0b42b80cc8d45f1910623)
		{
			printf("\t\tstore misc-output in[%s]\n", ndir);
		}
	}

	GMRFLib_sprintf(&nndir, "%s/theta-tags", ndir);
	{
		FILE *fp = fopen(nndir, "w");
		for (i = 0; i < ntheta; i++) {
			fprintf(fp, "%s\n", theta_tag[i]);
		}
		fclose(fp);
	}
	Free(nndir);

	for (i = any = 0; i < ntheta; i++) {
		any = (any || theta_from[i]);
	}
	if (any) {
		GMRFLib_sprintf(&nndir, "%s/theta-from", ndir);
		FILE *fp = fopen(nndir, "w");
		for (i = 0; i < ntheta; i++) {
			fprintf(fp, "%s\n", theta_from[i]);
		}
		fclose(fp);
	}
	Free(nndir);

	for (i = any = 0; i < ntheta; i++) {
		any = (any || theta_to[i]);
	}
	if (any) {
		GMRFLib_sprintf(&nndir, "%s/theta-to", ndir);
		FILE *fp = fopen(nndir, "w");
		for (i = 0; i < ntheta; i++) {
			fprintf(fp, "%s\n", theta_to[i]);
		}
		fclose(fp);
		Free(nndir);
	}

	GMRFLib_sprintf(&nndir, "%s/%s", ndir, "covmat-hyper-internal.dat");
	{
		Dinit(nndir);
		D1W(mo->nhyper);
		for (i = 0; i < ISQR(mo->nhyper); i++) {
			D1W(mo->cov_m[i]);
		}
		Dclose();
	}
	Free(nndir);


	GMRFLib_sprintf(&nndir, "%s/%s", ndir, "covmat-eigenvectors.dat");
	{
		Dinit(nndir);
		D1W(mo->nhyper);
		for (i = 0; i < ISQR(mo->nhyper); i++) {
			D1W(mo->eigenvectors[i]);
		}
		Dclose();
	}
	Free(nndir);

	GMRFLib_sprintf(&nndir, "%s/%s", ndir, "covmat-eigenvalues.dat");
	{
		Dinit(nndir);
		D1W(mo->nhyper);
		for (i = 0; i < mo->nhyper; i++) {
			D1W(mo->eigenvalues[i]);
		}
		Dclose();
	}

	Free(nndir);

	GMRFLib_sprintf(&nndir, "%s/%s", ndir, "reordering.dat");
	{
		Dinit(nndir);
		for (i = 0; i < mo->len_reordering; i++) {
			D1W(mo->reordering[i] + 1);	       /* yes, use 1-based indexing. */
		}
		Dclose();
	}
	Free(nndir);

	GMRFLib_sprintf(&nndir, "%s/%s", ndir, "mode-status.dat");
	{
		FILE *fp = fopen(nndir, "w");
		if (!fp) {
			inla_error_open_file(nndir);
		}
		fprintf(fp, "%1d\n", mo->mode_status);
		fclose(fp);
	}
	Free(nndir);

	GMRFLib_sprintf(&nndir, "%s/%s", ndir, "nfunc.dat");
	{
		FILE *fp = fopen(nndir, "w");
		if (!fp) {
			inla_error_open_file(nndir);
		}
		fprintf(fp, "%1d\n", mo->nfunc);
		fclose(fp);
	}
	Free(nndir);

	GMRFLib_sprintf(&nndir, "%s/%s", ndir, "log-posterior-mode.dat");
	{
		FILE *fp = fopen(nndir, "w");
		if (!fp) {
			inla_error_open_file(nndir);
		}
		fprintf(fp, "%.16f\n", mo->log_posterior_mode);
		fclose(fp);
	}
	Free(nndir);

	GMRFLib_sprintf(&nndir, "%s/%s", ndir, "stdev_corr_pos.dat");
	if (mo->stdev_corr_pos) {
		GMRFLib_matrix_tp *M = Calloc(1, GMRFLib_matrix_tp);
		M->nrow = mo->nhyper;
		M->ncol = 1;
		M->elems = M->nrow * M->ncol;
		M->A = Calloc(mo->nhyper, double);
		Memcpy(M->A, mo->stdev_corr_pos, mo->nhyper * sizeof(double));
		GMRFLib_write_fmesher_file(M, nndir, 0L, -1);
		GMRFLib_matrix_free(M);
	}

	GMRFLib_sprintf(&nndir, "%s/%s", ndir, "stdev_corr_neg.dat");
	if (mo->stdev_corr_pos) {
		GMRFLib_matrix_tp *M = Calloc(1, GMRFLib_matrix_tp);
		M->nrow = mo->nhyper;
		M->ncol = 1;
		M->elems = M->nrow * M->ncol;
		M->A = Calloc(mo->nhyper, double);
		Memcpy(M->A, mo->stdev_corr_neg, mo->nhyper * sizeof(double));
		GMRFLib_write_fmesher_file(M, nndir, 0L, -1);
		GMRFLib_matrix_free(M);
	}

	GMRFLib_sprintf(&nndir, "%s/%s", ndir, "opt_directions.dat");
	if (mo->opt_directions) {
		GMRFLib_write_fmesher_file(mo->opt_directions, nndir, 0L, -1);
	}

	if (mo->compute_corr_lin && mo->corr_lin) {
		/*
		 * OOPS: this matrix is in its own internal ordering, so we need to fix it here.
		 */
		assert(lc_order);
		int n = mo->compute_corr_lin;
		int *order = Calloc(n, int);

		for (i = 0; i < n; i++) {
			order[i] = (int) lc_order[i] - 1;
		}
		assert(GMRFLib_imin_value(order, n, NULL) == 0);
		assert(GMRFLib_imax_value(order, n, NULL) == n - 1);
		inla_output_matrix(ndir, NULL, "lincomb_derived_correlation_matrix.dat", n, mo->corr_lin, order);
		inla_output_matrix(ndir, NULL, "lincomb_derived_covariance_matrix.dat", n, mo->cov_lin, order);
		Free(order);
	}

	if (mo->configs) {
		FILE *fp;

		GMRFLib_sprintf(&nndir, "%s/%s", ndir, "config");
		if (inla_mkdir(nndir) != 0) {
			GMRFLib_sprintf(&msg, "fail to create directory [%s]: %s", nndir, strerror(errno));
			inla_error_general(msg);
		}

		GMRFLib_sprintf(&nnndir, "%s/%s", nndir, "theta-tag.dat");
		fp = fopen(nnndir, "w");
		for (i = 0; i < mb->ntheta; i++) {
			fprintf(fp, "%s\n", mb->theta_tag[i]);
		}
		fclose(fp);

		GMRFLib_sprintf(&nnndir, "%s/%s", nndir, "tag.dat");
		fp = fopen(nnndir, "w");
		for (i = 0; i < mb->idx_tot; i++) {
			fprintf(fp, "%s\n", mb->idx_tag[i]);
		}
		fclose(fp);

		GMRFLib_sprintf(&nnndir, "%s/%s", nndir, "start.dat");
		fp = fopen(nnndir, "w");
		for (i = 0; i < mb->idx_tot; i++) {
			fprintf(fp, "%d\n", mb->idx_start[i]);
		}
		fclose(fp);

		GMRFLib_sprintf(&nnndir, "%s/%s", nndir, "n.dat");
		fp = fopen(nnndir, "w");
		for (i = 0; i < mb->idx_tot; i++) {
			fprintf(fp, "%d\n", mb->idx_n[i]);
		}
		fclose(fp);

		GMRFLib_sprintf(&nnndir, "%s/%s", nndir, "configs.dat");
		fp = fopen(nnndir, "wb");
		int id = 0, header = 0, nconfig = 0;
		for (id = 0; id < GMRFLib_MAX_THREADS(); id++) {
			if (mo->configs[id]) {
				nconfig += mo->configs[id]->nconfig;	/* need the accumulated one! */
			}
		}

		for (id = 0; id < GMRFLib_MAX_THREADS(); id++) {
			if (mo->configs[id]) {

				if (!header) {
					header = 1;	       /* do this only once */
					fwrite((void *) &(mo->configs[id]->n), sizeof(int), (size_t) 1, fp);
					fwrite((void *) &(mo->configs[id]->nz), sizeof(int), (size_t) 1, fp);
					fwrite((void *) &(mo->configs[id]->ntheta), sizeof(int), (size_t) 1, fp);
					fwrite((void *) mo->configs[id]->i, sizeof(int), (size_t) mo->configs[id]->nz, fp);	/* 0-based! */
					fwrite((void *) mo->configs[id]->j, sizeof(int), (size_t) mo->configs[id]->nz, fp);	/* 0-based! */
					fwrite((void *) &nconfig, sizeof(int), (size_t) 1, fp);	/* yes!!! */

					if (mo->configs[id]->constr) {
						fwrite((void *) &(mo->configs[id]->constr->nc), sizeof(int), (size_t) 1, fp);
						fwrite((void *) mo->configs[id]->constr->a_matrix, sizeof(double),
						       (size_t) (mo->configs[id]->n * mo->configs[id]->constr->nc), fp);
						fwrite((void *) mo->configs[id]->constr->e_vector, sizeof(double),
						       (size_t) mo->configs[id]->constr->nc, fp);
					} else {
						int zero = 0;
						fwrite((void *) &zero, sizeof(int), (size_t) 1, fp);
					}
				}

				double *off = Calloc(mo->configs[id]->n, double);
				Memcpy(off, &(OFFSET3(0)), (mb->predictor_n + mb->predictor_m) * sizeof(double));

				for (i = 0; i < mo->configs[id]->nconfig; i++) {
					fwrite((void *) &(mo->configs[id]->config[i]->log_posterior), sizeof(double), (size_t) 1, fp);
					fwrite((void *) &(mo->configs[id]->config[i]->log_posterior_orig), sizeof(double), (size_t) 1, fp);
					fwrite((void *) mo->configs[id]->config[i]->theta, sizeof(double), (size_t) mo->configs[id]->ntheta, fp);
					fwrite((void *) mo->configs[id]->config[i]->mean, sizeof(double), (size_t) mo->configs[id]->n, fp);
					fwrite((void *) mo->configs[id]->config[i]->improved_mean, sizeof(double), (size_t) mo->configs[id]->n, fp);
					fwrite((void *) mo->configs[id]->config[i]->skewness, sizeof(double), (size_t) mo->configs[id]->n, fp);
					fwrite((void *) off, sizeof(double), (size_t) mo->configs[id]->n, fp);
					fwrite((void *) mo->configs[id]->config[i]->Q, sizeof(double), (size_t) mo->configs[id]->nz, fp);
					fwrite((void *) mo->configs[id]->config[i]->Qinv, sizeof(double), (size_t) mo->configs[id]->nz, fp);
					fwrite((void *) mo->configs[id]->config[i]->Qprior, sizeof(double), (size_t) mo->configs[id]->n, fp);
				}

				Free(off);
			}
		}
		fclose(fp);
	}

	if (mo->configs_preopt) {

		FILE *fp;
		GMRFLib_sprintf(&nndir, "%s/%s", ndir, "config_preopt");
		if (inla_mkdir(nndir) != 0) {
			GMRFLib_sprintf(&msg, "fail to create directory [%s]: %s", nndir, strerror(errno));
			inla_error_general(msg);
		}

		GMRFLib_sprintf(&nnndir, "%s/%s", nndir, "theta-tag.dat");
		fp = fopen(nnndir, "w");
		for (i = 0; i < mb->ntheta; i++) {
			fprintf(fp, "%s\n", mb->theta_tag[i]);
		}
		fclose(fp);

		GMRFLib_sprintf(&nnndir, "%s/%s", nndir, "tag.dat");
		fp = fopen(nnndir, "w");
		for (i = 0; i < mb->idx_tot; i++) {
			fprintf(fp, "%s\n", mb->idx_tag[i]);
		}
		fclose(fp);

		GMRFLib_sprintf(&nnndir, "%s/%s", nndir, "start.dat");
		fp = fopen(nnndir, "w");
		for (i = 0; i < mb->idx_tot; i++) {
			fprintf(fp, "%d\n", mb->idx_start[i]);
		}
		fclose(fp);

		GMRFLib_sprintf(&nnndir, "%s/%s", nndir, "n.dat");
		fp = fopen(nnndir, "w");
		for (i = 0; i < mb->idx_tot; i++) {
			fprintf(fp, "%d\n", mb->idx_n[i]);
		}
		fclose(fp);

		GMRFLib_sprintf(&nnndir, "%s/%s", nndir, "configs.dat");
		fp = fopen(nnndir, "wb");
		int id = 0, header = 0, nconfig = 0;
		for (id = 0; id < GMRFLib_MAX_THREADS(); id++) {
			if (mo->configs_preopt[id]) {
				nconfig += mo->configs_preopt[id]->nconfig;	/* need the accumulated one! */
			}
		}

		for (id = 0; id < GMRFLib_MAX_THREADS(); id++) {
			if (mo->configs_preopt[id]) {
				if (!header) {
					header = 1;	       /* do this only once */
					fwrite((void *) &(mo->configs_preopt[id]->mpred), sizeof(int), (size_t) 1, fp);
					fwrite((void *) &(mo->configs_preopt[id]->npred), sizeof(int), (size_t) 1, fp);
					fwrite((void *) &(mo->configs_preopt[id]->mnpred), sizeof(int), (size_t) 1, fp);
					fwrite((void *) &(mo->configs_preopt[id]->Npred), sizeof(int), (size_t) 1, fp);
					fwrite((void *) &(mo->configs_preopt[id]->n), sizeof(int), (size_t) 1, fp);
					fwrite((void *) &(mo->configs_preopt[id]->nz), sizeof(int), (size_t) 1, fp);
					fwrite((void *) &(mo->configs_preopt[id]->prior_nz), sizeof(int), (size_t) 1, fp);
					fwrite((void *) &(mo->configs_preopt[id]->ntheta), sizeof(int), (size_t) 1, fp);
					fwrite((void *) mo->configs_preopt[id]->i, sizeof(int), (size_t) mo->configs_preopt[id]->nz, fp);	/* 0-based! 
																		 */
					fwrite((void *) mo->configs_preopt[id]->j, sizeof(int), (size_t) mo->configs_preopt[id]->nz, fp);	/* 0-based! 
																		 */
					fwrite((void *) mo->configs_preopt[id]->iprior, sizeof(int), (size_t) mo->configs_preopt[id]->prior_nz, fp);	/* 0-based! 
																			 */
					fwrite((void *) mo->configs_preopt[id]->jprior, sizeof(int), (size_t) mo->configs_preopt[id]->prior_nz, fp);	/* 0-based! 
																			 */
					fwrite((void *) &nconfig, sizeof(int), (size_t) 1, fp);	/* yes!!! */

					if (mo->configs_preopt[id]->constr) {
						fwrite((void *) &(mo->configs_preopt[id]->constr->nc), sizeof(int), (size_t) 1, fp);
						fwrite((void *) mo->configs_preopt[id]->constr->a_matrix, sizeof(double),
						       (size_t) (mo->configs_preopt[id]->n * mo->configs_preopt[id]->constr->nc), fp);
						fwrite((void *) mo->configs_preopt[id]->constr->e_vector, sizeof(double),
						       (size_t) mo->configs_preopt[id]->constr->nc, fp);
					} else {
						int zero = 0;
						fwrite((void *) &zero, sizeof(int), (size_t) 1, fp);
					}

					double *off = Calloc(mo->configs_preopt[id]->mnpred, double);
					Memcpy(off, &(OFFSET3(0)), mo->configs_preopt[id]->mnpred * sizeof(double));
					fwrite((void *) off, sizeof(double), (size_t) mo->configs_preopt[id]->mnpred, fp);
					Free(off);

					char *A, *pA;
					GMRFLib_sprintf(&A, "%s/%s", nndir, "A.dat");
					GMRFLib_write_fmesher_file(mo->configs_preopt[id]->A, A, (long int) 0, -1);
					GMRFLib_sprintf(&pA, "%s/%s", nndir, "pA.dat");
					GMRFLib_write_fmesher_file(mo->configs_preopt[id]->pA, pA, (long int) 0, -1);
				}

				for (i = 0; i < mo->configs_preopt[id]->nconfig; i++) {
					fwrite((void *) &(mo->configs_preopt[id]->config[i]->log_posterior), sizeof(double), (size_t) 1, fp);
					fwrite((void *) &(mo->configs_preopt[id]->config[i]->log_posterior_orig), sizeof(double), (size_t) 1, fp);
					fwrite((void *) mo->configs_preopt[id]->config[i]->theta, sizeof(double),
					       (size_t) mo->configs_preopt[id]->ntheta, fp);
					fwrite((void *) mo->configs_preopt[id]->config[i]->mean, sizeof(double), (size_t) mo->configs_preopt[id]->n,
					       fp);
					fwrite((void *) mo->configs_preopt[id]->config[i]->improved_mean, sizeof(double),
					       (size_t) mo->configs_preopt[id]->n, fp);
					fwrite((void *) mo->configs_preopt[id]->config[i]->Q, sizeof(double), (size_t) mo->configs_preopt[id]->nz,
					       fp);
					fwrite((void *) mo->configs_preopt[id]->config[i]->Qinv, sizeof(double),
					       (size_t) mo->configs_preopt[id]->nz, fp);
					fwrite((void *) mo->configs_preopt[id]->config[i]->Qprior, sizeof(double),
					       (size_t) mo->configs_preopt[id]->prior_nz, fp);

					double output[2] = {
						(mo->configs_preopt[id]->config[i]->cpodens_moments ? 1.0 : 0.0),
						(mo->configs_preopt[id]->config[i]->gcpodens_moments ? 1.0 : 0.0)
					};

					fwrite((void *) output, sizeof(double), (size_t) 2L, fp);
					if (output[0]) {
						fwrite((void *) mo->configs_preopt[id]->config[i]->cpodens_moments, sizeof(double),
						       (size_t) mo->configs_preopt[id]->Npred * 3, fp);
					}
					if (output[1]) {
						fwrite((void *) mo->configs_preopt[id]->config[i]->gcpodens_moments, sizeof(double),
						       (size_t) mo->configs_preopt[id]->Npred * 3, fp);
					}
					if (mo->configs_preopt[id]->config[i]->arg_str) {
						double one = 1.0;
						fwrite((void *) &one, sizeof(double), 1, fp);

						char **s = mo->configs_preopt[id]->config[i]->arg_str;
						char *null_str = GMRFLib_strdup("");
						for (int ii = 0; ii < mo->configs_preopt[id]->Npred; ii++) {
							char *str = (s[ii] ? s[ii] : null_str);
							fwrite((void *) str, 1, (size_t) strlen(str) + 1, fp);
						}
					} else {
						double zero = 0.0;
						fwrite((void *) &zero, sizeof(double), 1, fp);
					}

					if (mo->configs_preopt[id]->config[i]->ll_info) {
						double one = 1.0;
						fwrite((void *) &one, sizeof(double), 1, fp);
						fwrite((void *) (mo->configs_preopt[id]->config[i]->ll_info),
						       sizeof(double), (size_t) (3 * mo->configs_preopt[id]->Npred), fp);
					} else {
						double zero = 0.0;
						fwrite((void *) &zero, sizeof(double), 1, fp);
					}

					if (mo->configs_preopt[id]->config[i]->lpred_mean) {
						double one = 1.0;
						fwrite((void *) &one, sizeof(double), 1, fp);
						fwrite((void *) (mo->configs_preopt[id]->config[i]->lpred_mean),
						       sizeof(double), (size_t) mo->configs_preopt[id]->mnpred, fp);
						fwrite((void *) (mo->configs_preopt[id]->config[i]->lpred_variance),
						       sizeof(double), (size_t) mo->configs_preopt[id]->mnpred, fp);
					} else {
						double zero = 0.0;
						fwrite((void *) &zero, sizeof(double), 1, fp);
					}
				}
			}
		}
		fclose(fp);
	}

	GMRFLib_sprintf(&nnndir, "%s/%s", ndir, "warnings.txt");
	{
		if (mo->warnings) {
			FILE *fp = fopen(nnndir, "w");
			for (int k = 0;; k++) {
				if (mo->warnings[k]) {
					fprintf(fp, "%s\n", mo->warnings[k]);
				} else {
					break;
				}
			}
			fclose(fp);
		}
	}
	Free(nnndir);

	return INLA_OK;
}

int inla_output_detail_mlik(const char *dir, GMRFLib_ai_marginal_likelihood_tp *mlik, int verbose)
{
	/*
	 * output whatever is requested.... 
	 */
	char *ndir = NULL, *msg = NULL, *nndir = NULL;
	if (!mlik) {
		return INLA_OK;
	}
	GMRFLib_sprintf(&ndir, "%s/%s", dir, "marginal-likelihood");
	if (inla_mkdir(ndir) != 0) {
		GMRFLib_sprintf(&msg, "fail to create directory [%s]: %s", ndir, strerror(errno));
		inla_error_general(msg);
	}
	GMRFLib_sprintf(&nndir, "%s/%s", ndir, "marginal-likelihood.dat");
	Dinit_s(nndir);
	if (verbose) {
#pragma omp critical (Name_9dd76c0b6445188affe98f4c81bec8e5c5e3a95b)
		{
			printf("\t\tstore marginal-likelihood results in[%s]\n", nndir);
		}
	}
	D2W(mlik->marginal_likelihood_integration, mlik->marginal_likelihood_gaussian_approx);
	Dclose();

	Free(ndir);
	Free(nndir);

	return INLA_OK;
}

int inla_output_gitid(const char *dir)
{
	char *nndir = NULL;
	FILE *fp = NULL;

	GMRFLib_sprintf(&nndir, "%s/%s", dir, ".gitid");
	fp = fopen(nndir, "w");
	if (!fp) {
		inla_error_open_file(nndir);
	}
	fprintf(fp, "GITCOMMIT [%s]\n", GITCOMMIT);
	fclose(fp);
	Free(nndir);

	return INLA_OK;
}

int inla_output_linkfunctions(const char *dir, inla_tp *mb)
{
	int i, j;
	char *nndir = NULL;
	FILE *fp = NULL;

	GMRFLib_sprintf(&nndir, "%s/%s", dir, "linkfunctions.names");
	fp = fopen(nndir, "w");
	if (!fp) {
		inla_error_open_file(nndir);
	}

	for (j = 0; j < mb->nds; j++) {
		link_func_tp *lf = mb->data_sections[j].predictor_invlinkfunc;

		if (lf == link_probit) {
			fprintf(fp, "probit\n");
		} else if (lf == link_tan) {
			fprintf(fp, "tan\n");
		} else if (lf == link_cloglog) {
			fprintf(fp, "cloglog\n");
		} else if (lf == link_ccloglog) {
			fprintf(fp, "ccloglog\n");
		} else if (lf == link_log) {
			fprintf(fp, "log\n");
		} else if (lf == link_logit) {
			fprintf(fp, "logit\n");
		} else if (lf == link_identity) {
			fprintf(fp, "identity\n");
		} else if (lf == link_inverse) {
			fprintf(fp, "inverse\n");
		} else if (lf == link_sslogit) {
			fprintf(fp, "sslogit\n");
		} else if (lf == link_robit) {
			fprintf(fp, "robit\n");
		} else if (lf == link_sn) {
			fprintf(fp, "sn\n");
		} else if (lf == link_logoffset) {
			fprintf(fp, "logoffset\n");
		} else if (lf == link_logitoffset) {
			fprintf(fp, "logitoffset\n");
		} else if (lf == link_qpoisson) {
			fprintf(fp, "quantile\n");
		} else if (lf == link_qbinomial) {
			fprintf(fp, "quantile\n");
		} else if (lf == link_pqbinomial) {
			fprintf(fp, "pquantile\n");
		} else if (lf == link_qweibull) {
			fprintf(fp, "quantile\n");
		} else if (lf == link_qweibull) {
			fprintf(fp, "quantile\n");
		} else if (lf == link_qgamma) {
			fprintf(fp, "quantile\n");
		} else if (lf == link_loga) {
			fprintf(fp, "loga\n");
		} else if (lf == link_power_logit) {
			fprintf(fp, "powerlogit\n");
		} else {
			fprintf(fp, "invalid-linkfunction\n");
		}
	}
	fclose(fp);
	Free(nndir);

	GMRFLib_sprintf(&nndir, "%s/%s", dir, "linkfunctions.link");
	fp = fopen(nndir, "wb");
	if (!fp) {
		inla_error_open_file(nndir);
	}

	/*
	 * need to use double as we need NAN
	 */
	double *idx = Calloc(mb->predictor_ndata, double);
	for (i = 0; i < mb->predictor_ndata; i++) {
		int found;

		if (!ISNAN(mb->predictor_family[i])) {
			idx[i] = (int) mb->predictor_family[i];
		} else {
			idx[i] = NAN;
			for (j = found = 0; j < mb->nds && !found; j++) {
				if (mb->data_sections[j].predictor_invlinkfunc == mb->predictor_invlinkfunc[i]) {
					found = 1;
					idx[i] = j;
				}
			}
		}
	}

	fwrite((void *) &(mb->predictor_ndata), sizeof(int), 1, fp);
	fwrite((void *) idx, sizeof(double), (size_t) mb->predictor_ndata, fp);
	Free(idx);

	fclose(fp);
	Free(nndir);

	return INLA_OK;
}

int inla_output_ok(const char *dir)
{
	char *nndir = NULL;
	FILE *fp = NULL;

	GMRFLib_sprintf(&nndir, "%s/%s", dir, ".ok");
	fp = fopen(nndir, "w");
	if (!fp) {
		inla_error_open_file(nndir);
	}
	fprintf(fp, "1");
	fclose(fp);
	Free(nndir);

	return INLA_OK;
}

int inla_output_detail_theta(const char *dir, double ***theta, int n_theta)
{
	/*
	 * write the mode of theta to the file DIR/.theta_mode. This is used only internally... 
	 */
	int i;
	char *nndir = NULL;

	GMRFLib_sprintf(&nndir, "%s/%s", dir, ".theta_mode");
	Dinit_s(nndir);
	D1W(n_theta);
	for (i = 0; i < n_theta; i++) {
		D1W(theta[i][0][0]);
	}
	Dclose();
	Free(nndir);

	return INLA_OK;
}

int inla_output_detail_x(const char *dir, double *x, int n_x)
{
	/*
	 * write the mode of x to the file DIR/.x_mode. This is used only internally... 
	 */
	int i;
	char *nndir = NULL;

	GMRFLib_sprintf(&nndir, "%s/%s", dir, ".x_mode");
	Dinit(nndir);
	D1W(n_x);
	for (i = 0; i < n_x; i++) {
		D1W(x[i]);
	}
	Dclose();
	Free(nndir);

	return INLA_OK;
}


int inla_output_detail(const char *dir, GMRFLib_density_tp **density, double *locations, int n, int nrep, Output_tp *output, const char *sdir,
		       int return_marginals,
		       // Either this
		       map_func_tp *func, void *func_arg,
		       // .. or this
		       GMRFLib_transform_array_func_tp **tfunc,
		       // 
		       const char *tag, const char *modelname, int UNUSED(verbose))
{
	int thread_id = 0;

#define _FUNC (func ? func : NULL)
#define _FUNC_ARG (func ? func_arg : NULL)
#define _TFUNC(_idx) (tfunc ? tfunc[_idx] : NULL)
#define _MAP_DENS(_dens, _x_user, _idx) (func ? (_dens)/(ABS(func(_x_user, MAP_DFORWARD, func_arg))) : \
					 (tfunc ? (_dens)/(ABS(tfunc[_idx]->func(thread_id, _x_user, GMRFLib_TRANSFORM_DFORWARD, tfunc[_idx]->arg, tfunc[_idx]->cov))) : \
					  (_dens)))
#define _MAP_X(_x_user, _idx) (func ? func(_x_user, MAP_FORWARD, func_arg) : \
			       (tfunc ? tfunc[_idx]->func(thread_id, _x_user, GMRFLib_TRANSFORM_FORWARD, tfunc[_idx]->arg, tfunc[_idx]->cov) : \
				(_x_user)))
#define _MAP_INCREASING(_idx) (func ? func(0.0, MAP_INCREASING, func_arg) : \
			       (tfunc ? tfunc[_idx]->func(thread_id, 0.0, GMRFLib_TRANSFORM_INCREASING, tfunc[_idx]->arg, tfunc[_idx]->cov) : 1))
#define _MAP_DECREASING(_idx) (!_MAP_INCREASING(_idx))
#define GMRFLib_MAX_THREADS_LOCAL() (n > 1024 ? GMRFLib_MAX_THREADS() : 2)

	GMRFLib_ENTER_ROUTINE;

	char *ndir = NULL, *ssdir = NULL, *msg = NULL;
	int ndiv;
	int add_empty = 1;
	int plain = ((func || tfunc) ? 0 : 1);

	assert(nrep > 0);
	ndiv = n / nrep;

	double *d_mode = Calloc(n, double);
	for (int i = 0; i < n; i++) {
		d_mode[i] = NAN;
	}

	ssdir = GMRFLib_strdup(sdir);
	GMRFLib_sprintf(&ndir, "%s/%s", dir, ssdir);
	if (inla_mkdir(ndir) != 0) {
		GMRFLib_sprintf(&msg, "fail to create directory [%s]: %s", ndir, strerror(errno));
		inla_error_general(msg);
	}
	Free(ssdir);

	if (1) {
		char *nndir = NULL;
		GMRFLib_sprintf(&nndir, "%s/%s", ndir, "N");
		FILE *fp = fopen(nndir, "w");
		if (!fp) {
			inla_error_open_file(nndir);
		}
		fprintf(fp, "%d\n", n);
		fclose(fp);
		Free(nndir);
	}

	if (tag) {
		char *nndir = NULL;
		GMRFLib_sprintf(&nndir, "%s/%s", ndir, "TAG");
		FILE *fp = fopen(nndir, "w");
		if (!fp) {
			inla_error_open_file(nndir);
		}
		fprintf(fp, "%s\n", tag);
		fclose(fp);
		Free(nndir);
	}

	if (modelname) {
		char *nndir = NULL;
		GMRFLib_sprintf(&nndir, "%s/%s", ndir, "MODEL");
		FILE *fp = fopen(nndir, "w");
		if (!fp) {
			inla_error_open_file(nndir);
		}
		fprintf(fp, "%s\n", modelname);
		fclose(fp);
		Free(nndir);
	}

	if (output->summary) {
		if (inla_computed(density, n)) {
			char *nndir = NULL;
			GMRFLib_sprintf(&nndir, "%s/%s", ndir, "summary.dat");
			Dinit_r(n, 3, nndir);

#define CODE_BLOCK							\
			for (int i = 0; i < n; i++) {			\
				double dm = 0.0, ds = 0.0;		\
				if (density[i]) {			\
					inla_integrate_func(&dm, &ds, &d_mode[i], density[i], _FUNC, _FUNC_ARG, _TFUNC(i)); \
					if (locations) {		\
						D3W_r(i, 0, locations[i % ndiv], dm, ds); \
					} else {			\
						D3W_r(i, 0, i, dm, ds);	\
					}				\
				} else {				\
					if (locations) {		\
						D3W_r(i, 0, locations[i % ndiv], NAN, NAN); \
					} else {			\
						D3W_r(i, 0, i, NAN, NAN); \
					}				\
				}					\
			}

			RUN_CODE_BLOCK(GMRFLib_MAX_THREADS_LOCAL(), 0, 0);
#undef CODE_BLOCK

			Dclose_r();
			Free(nndir);
		}
	}

	if (return_marginals || strncmp("hyperparameter", sdir, 13) == 0) {
		if (inla_computed(density, n)) {
			char *nndir = NULL;
			GMRFLib_sprintf(&nndir, "%s/%s", ndir, "marginal-densities.dat");
			int mm = 0;
			GMRFLib_density_layout_x(NULL, &mm, NULL);
			Dinit_r(n, 2 + mm * 2, nndir);

#define CODE_BLOCK							\
			for (int i = 0; i < n; i++) {			\
				int off = 0;				\
				int nn = mm;				\
				int nn_new;				\
				double *x_user = CODE_BLOCK_WORK_PTR(0); \
				double *dens = CODE_BLOCK_WORK_PTR(1);	\
				double *xx = CODE_BLOCK_WORK_PTR(2);	\
				if (density[i]) {			\
					if (locations) {		\
						D1W_r(i, off, locations[i % ndiv]); \
					} else {			\
						D1W_r(i, off, i);	\
					}				\
					off++;				\
									\
					GMRFLib_density_layout_x(xx, &nn_new, density[i]); assert(nn_new == nn); \
					GMRFLib_density_std2user_n(x_user, xx, nn, density[i]); \
					GMRFLib_evaluate_ndensity(dens, xx, nn, density[i]); \
					D1W_r(i, off, nn);		\
					off++;				\
									\
					if (plain) {			\
						for (int ii = 0; ii < nn; ii++) { \
							double dens_user = dens[ii] / density[i]->std_stdev; \
							D2W_r(i, off, x_user[ii], dens_user); \
							off += 2;	\
						}			\
					} else {			\
						for (int ii = 0; ii < nn; ii++) { \
							double dens_user = dens[ii] / density[i]->std_stdev; \
							D2W_r(i, off, _MAP_X(x_user[ii], i), _MAP_DENS(dens_user, x_user[ii], i)); \
							off += 2;	\
						}			\
					}				\
				} else {				\
					if (locations) {		\
						D1W_r(i, off, locations[i % ndiv]); \
					} else {			\
						D1W_r(i, off, i);	\
					}				\
					off++;				\
									\
					D1W_r(i, off, nn);		\
					off++;				\
									\
					for (int ii = 0; ii < nn; ii++) { \
						D2W_r(i, off, NAN, NAN); \
						off += 2;		\
					}				\
				}					\
			}

			// RUN_CODE_BLOCK(2, 3, mm);
			RUN_CODE_BLOCK(GMRFLib_MAX_THREADS_LOCAL(), 3, mm);
#undef CODE_BLOCK
			Dclose_r();
			Free(nndir);
		}
	}

	if (output->kld) {
		/*
		 * this is ok for _FUNC as well, since the the KL is invariant for parameter transformations. 
		 */
		if (inla_computed(density, n)) {
			char *nndir = NULL;
			GMRFLib_sprintf(&nndir, "%s/%s", ndir, "symmetric-kld.dat");
			Dinit(nndir);
			for (int i = 0; i < n; i++) {
				GMRFLib_density_tp *gd = NULL;
				if (density[i]) {
					double kld;
					GMRFLib_density_create_normal(&gd, 0.0, 1.0, density[i]->std_mean, density[i]->std_stdev, GMRFLib_FALSE);

					if (G.fast_mode) {
						GMRFLib_mkld_sym(&kld, gd, density[i]);
					} else {
						GMRFLib_kld_sym(&kld, gd, density[i]);
					}
					if (locations) {
						D1W(locations[i % ndiv]);
					} else {
						D1W(i);
					}
					D1W(kld);
				} else {
					if (add_empty) {
						if (locations) {
							D1W(locations[i % ndiv]);
						} else {
							D1W(i);
						}
						D1W(NAN);
					}
				}
				GMRFLib_free_density(gd);
			}
			Dclose();
			Free(nndir);
		}
	}

	if (output->nquantiles) {
		if (inla_computed(density, n)) {
			double x_user;
			char *nndir = NULL;
			GMRFLib_sprintf(&nndir, "%s/%s", ndir, "quantiles.dat");
			Dinit(nndir);
			for (int i = 0; i < n; i++) {
				double xp, p;
				if (density[i]) {
					if (locations) {
						D1W(locations[i % ndiv]);
					} else {
						D1W(i);
					}
					D1W(output->nquantiles);
					for (int j = 0; j < output->nquantiles; j++) {
						p = output->quantiles[j];
						if (_MAP_INCREASING(i)) {
							GMRFLib_density_Pinv(&xp, p, density[i]);
						} else {
							GMRFLib_density_Pinv(&xp, 1.0 - p, density[i]);
						}
						x_user = GMRFLib_density_std2user(xp, density[i]);
						D2W(p, _MAP_X(x_user, i));
					}
				} else {
					if (add_empty) {
						if (locations) {
							D1W(locations[i % ndiv]);
						} else {
							D1W(i);
						}
						D1W(output->nquantiles);
						for (int j = 0; j < output->nquantiles; j++) {
							D2W(NAN, NAN);
						}
					}
				}
			}
			Dclose();
			Free(nndir);
		}
	}

	if (output->mode) {
		if (inla_computed(density, n)) {
			char *nndir = NULL;
			GMRFLib_sprintf(&nndir, "%s/%s", ndir, "mode.dat");
			Dinit(nndir);
			for (int i = 0; i < n; i++) {
				if (density[i]) {
					if (locations) {
						D1W(locations[i % ndiv]);
					} else {
						D1W(i);
					}
					D3W(1.0, NAN, d_mode[i]);
				} else {
					if (add_empty) {
						if (locations) {
							D1W(locations[i % ndiv]);
						} else {
							D1W(i);
						}
						D3W(1.0, NAN, NAN);
					}
				}
			}
			Dclose();
			Free(nndir);
		}
	}

	if (output->ncdf) {
		if (inla_computed(density, n)) {
			char *nndir = NULL;
			GMRFLib_sprintf(&nndir, "%s/%s", ndir, "cdf.dat");
			Dinit(nndir);
			for (int i = 0; i < n; i++) {
				double xp, x, p;
				if (density[i]) {
					if (locations) {
						D1W(locations[i % ndiv]);
					} else {
						D1W(i);
					}
					D1W(output->ncdf);
					for (int j = 0; j < output->ncdf; j++) {
						xp = output->cdf[j];
						x = GMRFLib_density_user2std(xp, density[i]);
						GMRFLib_density_P(&p, x, density[i]);
						if (_MAP_DECREASING(i)) {
							p = 1.0 - p;
						}
						D2W(_MAP_X(xp, i), p);
					}
				} else {
					if (add_empty) {
						if (locations) {
							D1W(locations[i % ndiv]);
						} else {
							D1W(i);
						}
						D1W(output->ncdf);
						for (int j = 0; j < output->ncdf; j++) {
							D2W(NAN, NAN);
						}
					}
				}
			}
			Dclose();
			Free(nndir);
		}
	}

	Free(d_mode);

#undef _MAP_DENS
#undef _MAP_X
#undef _MAP_INCREASING
#undef _MAP_DECREASING
#undef _FUNC
#undef _FUNC_ARG
#undef _TFUNC
#undef GMRFLib_MAX_THREADS_LOCAL

	GMRFLib_LEAVE_ROUTINE;

	return INLA_OK;
}

int inla_parse_output(inla_tp *mb, dictionary *ini, int sec, Output_tp **out)
{
	/*
	 * parse the output-options. defaults are given in the type=output-section, which are initialised with program defaults if
	 * the program defaults are NULL. 
	 */
	int i, j, use_defaults = 1, ret, ngroups_eff = 0;
	char *secname = NULL, *tmp = NULL, *gfile = NULL, *sfile = NULL, *ffile = NULL;

	secname = GMRFLib_strdup(iniparser_getsecname(ini, sec));
	if (!mb->output) {
		/*
		 * set default options 
		 */
		assert(mb->output == *out);
		use_defaults = 1;			       /* to flag that we're reading mb->output */
		(*out) = Calloc(1, Output_tp);
		(*out)->gcpo = 0;
		(*out)->cpo = 0;
		(*out)->po = 0;
		(*out)->dic = 0;
		(*out)->summary = 1;
		(*out)->return_marginals = 1;
		(*out)->return_marginals_predictor = 0;
		(*out)->kld = 1;
		(*out)->mlik = 0;
		(*out)->q = 0;
		(*out)->graph = 0;
		(*out)->config = 0;
		(*out)->likelihood_info = 0;
		(*out)->internal_opt = 1;
		(*out)->save_memory = 0;
		(*out)->hyperparameters = (G.mode == INLA_MODE_HYPER ? 1 : 1);
		(*out)->mode = 1;
		(*out)->nquantiles = 0;
		(*out)->ncdf = 0;
		(*out)->quantiles = (*out)->cdf = NULL;
	} else {
		use_defaults = 0;
		*out = Calloc(1, Output_tp);
		(*out)->gcpo = mb->output->gcpo;
		(*out)->cpo = mb->output->cpo;
		(*out)->po = mb->output->po;
		(*out)->dic = mb->output->dic;
		(*out)->summary = mb->output->summary;
		(*out)->kld = mb->output->kld;
		(*out)->mlik = mb->output->mlik;
		(*out)->q = mb->output->q;
		(*out)->graph = mb->output->graph;
		(*out)->config = mb->output->config;
		(*out)->likelihood_info = mb->output->likelihood_info;
		(*out)->internal_opt = mb->output->internal_opt;
		(*out)->save_memory = mb->output->save_memory;
		(*out)->hyperparameters = mb->output->hyperparameters;
		(*out)->mode = mb->output->mode;
		(*out)->return_marginals = mb->output->return_marginals;
		(*out)->return_marginals_predictor = mb->output->return_marginals_predictor;
		(*out)->nquantiles = mb->output->nquantiles;
		if (mb->output->nquantiles) {
			(*out)->quantiles = Calloc(mb->output->nquantiles, double);
			Memcpy((*out)->quantiles, mb->output->quantiles, (size_t) mb->output->nquantiles * sizeof(double));
		}
		(*out)->ncdf = mb->output->ncdf;
		if (mb->output->ncdf) {
			(*out)->cdf = Calloc(mb->output->ncdf, double);
			Memcpy((*out)->cdf, mb->output->cdf, (size_t) mb->output->ncdf * sizeof(double));
		}
	}

	if (!(mb->gcpo_param)) {
		mb->gcpo_param = Calloc(1, GMRFLib_gcpo_param_tp);
		mb->gcpo_param->num_level_sets = iniparser_getint(ini, inla_string_join(secname, "GCPO.NUM.LEVEL.SETS"), -1);
		mb->gcpo_param->size_max = iniparser_getint(ini, inla_string_join(secname, "GCPO.SIZE.MAX"), -1);
		mb->gcpo_param->correct_hyperpar = iniparser_getboolean(ini, inla_string_join(secname, "GCPO.CORRECT.HYPERPAR"), 1);
		mb->gcpo_param->epsilon = iniparser_getdouble(ini, inla_string_join(secname, "GCPO.EPSILON"), GSL_ROOT3_DBL_EPSILON);
		mb->gcpo_param->prior_diagonal = iniparser_getdouble(ini, inla_string_join(secname, "GCPO.PRIOR.DIAGONAL"), 1.0);
		mb->gcpo_param->remove_fixed = iniparser_getboolean(ini, inla_string_join(secname, "GCPO.REMOVE.FIXED"), 1);
		mb->gcpo_param->verbose = iniparser_getboolean(ini, inla_string_join(secname, "GCPO.VERBOSE"), 0);

		char *str = NULL;
		char *str_ptr = NULL;
		char *token = NULL;
		const char *delim = " \t";
		str = GMRFLib_strdup(iniparser_getstring(ini, inla_string_join(secname, "GCPO.KEEP"), NULL));
		while ((token = GMRFLib_strtok_r(str, delim, &str_ptr))) {
			str = NULL;
			GMRFLib_str_add(&(mb->gcpo_param->keep), token);
		}

		str_ptr = NULL;
		str = GMRFLib_strdup(iniparser_getstring(ini, inla_string_join(secname, "GCPO.REMOVE"), NULL));
		while ((token = GMRFLib_strtok_r(str, delim, &str_ptr))) {
			str = NULL;
			GMRFLib_str_add(&(mb->gcpo_param->remove), token);
		}

		char *tstr = NULL;
		tstr = GMRFLib_strdup(iniparser_getstring(ini, inla_string_join(secname, "GCPO.STRATEGY"), GMRFLib_strdup("posterior")));
		if (!strcasecmp(tstr, "posterior")) {
			mb->gcpo_param->build_strategy = GMRFLib_GCPO_BUILD_STRATEGY_POSTERIOR;
		} else if (!strcasecmp(tstr, "prior")) {
			mb->gcpo_param->build_strategy = GMRFLib_GCPO_BUILD_STRATEGY_PRIOR;
		} else {
			inla_error_field_is_void(__GMRFLib_FuncName, secname, "gcpo.strategy", tstr);
			assert(0 == 1);
		}

		gfile = GMRFLib_strdup(iniparser_getstring(ini, inla_string_join(secname, "GCPO.GROUPS"), NULL));
		sfile = GMRFLib_strdup(iniparser_getstring(ini, inla_string_join(secname, "GCPO.SELECTION"), NULL));
		ffile = GMRFLib_strdup(iniparser_getstring(ini, inla_string_join(secname, "GCPO.FRIENDS"), NULL));
		assert(!(gfile && sfile && ffile));

		if (gfile) {
			FILE *fp = fopen(gfile, "rb");
			int len, total_len, glen, offset = 0;
			ret = fread((void *) &len, sizeof(int), (size_t) 1, fp);
			assert(ret == 1);
			ret = fread((void *) &total_len, sizeof(int), (size_t) 1, fp);
			assert(ret == 1);

			if (mb->gcpo_param->verbose) {
				printf("%s: read groups len %d total_len %d\n", __GMRFLib_FuncName, len, total_len);
			}
			int *buffer = Calloc(total_len, int);
			ret = fread((void *) buffer, sizeof(int), (size_t) total_len, fp);
			assert(ret == total_len);

			mb->gcpo_param->groups = GMRFLib_idxval_ncreate_x(len, IMAX(3, (total_len - len) / len));
			for (i = 0; i < len; i++) {
				glen = buffer[offset++];
				for (j = 0; j < glen; j++) {
					GMRFLib_idxval_add(&(mb->gcpo_param->groups[i]), buffer[offset++], NAN);
				}
				if (mb->gcpo_param->verbose) {
					if (mb->gcpo_param->groups[i]->n > 0) {
						char *msg;
						GMRFLib_sprintf(&msg, "group %d", i);
						GMRFLib_idxval_printf(stdout, mb->gcpo_param->groups[i], msg);
					}
				}
			}
			mb->gcpo_param->ngroups = len;
			for (i = 0; i < len; i++) {
				ngroups_eff += (mb->gcpo_param->groups[i]->n > 0 ? 1 : 0);
			}
			assert(offset == total_len);
			fclose(fp);
			Free(buffer);
		} else {
			mb->gcpo_param->ngroups = -1;
			if (sfile) {
				FILE *fp = fopen(sfile, "rb");
				int len;
				ret = fread((void *) &len, sizeof(int), (size_t) 1, fp);
				assert(ret == 1);
				if (mb->gcpo_param->verbose) {
					printf("%s: read selection len %d\n", __GMRFLib_FuncName, len);
				}
				int *buffer = Calloc(len, int);
				ret = fread((void *) buffer, sizeof(int), (size_t) len, fp);
				assert(ret == len);
				for (i = 0; i < len; i++) {
					if (mb->gcpo_param->verbose) {
						printf("%s: add idx %d\n", __GMRFLib_FuncName, buffer[i]);
					}
					GMRFLib_idx_add(&(mb->gcpo_param->selection), buffer[i]);
				}
				fclose(fp);
				Free(buffer);
			}
			if (ffile) {
				FILE *fp = fopen(ffile, "rb");
				int len;
				ret = fread((void *) &len, sizeof(int), (size_t) 1, fp);
				assert(ret == 1);
				mb->gcpo_param->friends_n = len;
				if (mb->gcpo_param->verbose) {
					printf("%s: read friends len %d\n", __GMRFLib_FuncName, len);
				}
				assert(len >= 0);

				int len_buffer = 64;
				int *buffer = Calloc(len_buffer, int);
				if (len) {
					mb->gcpo_param->friends = GMRFLib_idx_ncreate_x(len, 4);
					for (i = 0; i < len; i++) {
						if (mb->gcpo_param->verbose) {
							printf("%s: add friends for i=[%1d]: ", __GMRFLib_FuncName, i);
						}
						int local_len;
						ret = fread((void *) &local_len, sizeof(int), (size_t) 1, fp);
						assert(ret == 1);
						if (local_len > 0) {
							if (local_len > len_buffer) {
								len_buffer = local_len;
								Free(buffer);
								buffer = Calloc(len_buffer, int);
							}
							ret = fread((void *) buffer, sizeof(int), (size_t) local_len, fp);
							assert(ret == local_len);
						}
						for (j = 0; j < local_len; j++) {
							GMRFLib_idx_add(&(mb->gcpo_param->friends[i]), buffer[j]);
							if (mb->gcpo_param->verbose) {
								printf(" %1d", buffer[j]);
							}
						}
						if (mb->gcpo_param->verbose) {
							printf("\n");
						}
					}
					fclose(fp);
				}
				Free(buffer);
			}
		}
	}

	(*out)->gcpo = iniparser_getboolean(ini, inla_string_join(secname, "GCPO.ENABLE"), (*out)->gcpo);
	(*out)->cpo = iniparser_getboolean(ini, inla_string_join(secname, "CPO"), (*out)->cpo);
	(*out)->po = iniparser_getboolean(ini, inla_string_join(secname, "PO"), (*out)->po);
	(*out)->dic = iniparser_getboolean(ini, inla_string_join(secname, "DIC"), (*out)->dic);
	(*out)->summary = iniparser_getboolean(ini, inla_string_join(secname, "SUMMARY"), (*out)->summary);
	(*out)->return_marginals = iniparser_getboolean(ini, inla_string_join(secname, "RETURN.MARGINALS"), (*out)->return_marginals);
	(*out)->return_marginals_predictor =
	    iniparser_getboolean(ini, inla_string_join(secname, "RETURN.MARGINALS.PREDICTOR"), (*out)->return_marginals_predictor);
	(*out)->hyperparameters = iniparser_getboolean(ini, inla_string_join(secname, "HYPERPARAMETERS"), (*out)->hyperparameters);
	(*out)->mode = iniparser_getboolean(ini, inla_string_join(secname, "MODE"), (*out)->mode);
	(*out)->kld = iniparser_getboolean(ini, inla_string_join(secname, "KLD"), (*out)->kld);
	(*out)->mlik = iniparser_getboolean(ini, inla_string_join(secname, "MLIK"), (*out)->mlik);
	(*out)->q = iniparser_getboolean(ini, inla_string_join(secname, "Q"), (*out)->q);
	(*out)->graph = iniparser_getboolean(ini, inla_string_join(secname, "GRAPH"), (*out)->graph);
	(*out)->config = iniparser_getboolean(ini, inla_string_join(secname, "CONFIG"), (*out)->config);
	(*out)->likelihood_info = iniparser_getboolean(ini, inla_string_join(secname, "LIKELIHOOD.INFO"), (*out)->likelihood_info);
	(*out)->internal_opt = GMRFLib_internal_opt = iniparser_getboolean(ini, inla_string_join(secname, "INTERNAL.OPT"), (*out)->internal_opt);
	(*out)->save_memory = GMRFLib_save_memory = iniparser_getboolean(ini, inla_string_join(secname, "SAVE.MEMORY"), (*out)->save_memory);

	if ((*out)->likelihood_info) {
		(*out)->config = 1;
	}

	tmp = GMRFLib_strdup(iniparser_getstring(ini, inla_string_join(secname, "QUANTILES"), NULL));

	if (G.mode == INLA_MODE_HYPER) {
		/*
		 * these are the requirements for the HYPER_MODE 
		 */
		(*out)->gcpo = 0;
		(*out)->cpo = 0;
		(*out)->po = 0;
		(*out)->dic = 0;
		(*out)->mlik = 1;
	}
	if (G.mode == INLA_MODE_HYPER) {
		if (!((*out)->hyperparameters)) {
			fprintf(stderr, "*** Warning: HYPER_MODE require (*out)->hyperparameters = 1\n");
		}
		(*out)->hyperparameters = 1;
	}

	if (tmp) {
		inla_sread_doubles_q(&((*out)->quantiles), &((*out)->nquantiles), tmp);

		if ((*out)->nquantiles == 0)
			inla_error_field_is_void(__GMRFLib_FuncName, secname, "quantiles", tmp);

		for (j = 0; j < (*out)->nquantiles; j++) {
			if ((*out)->quantiles[j] <= 0.0 || (*out)->quantiles[j] >= 1.0) {
				inla_error_field_is_void(__GMRFLib_FuncName, secname, "quantiles", tmp);
			}
		}
	}
	tmp = GMRFLib_strdup(iniparser_getstring(ini, inla_string_join(secname, "CDF"), NULL));
	if (tmp) {
		inla_sread_doubles_q(&((*out)->cdf), &((*out)->ncdf), tmp);
		if ((*out)->ncdf == 0) {
			inla_error_field_is_void(__GMRFLib_FuncName, secname, "cdf", tmp);
		}
	}
	if (mb->verbose) {
		printf("\t\toutput:\n");
		if (use_defaults) {
			printf("\t\t\tgcpo=[%1d]\n", (*out)->gcpo);
			printf("\t\t\t\tnum.level.sets=[%1d]\n", mb->gcpo_param->num_level_sets);
			printf("\t\t\t\tsize.max=[%1d]\n", mb->gcpo_param->size_max);
			printf("\t\t\t\tstrategy=[%s]\n", GMRFLib_GCPO_BUILD_STRATEGY_NAME(mb->gcpo_param->build_strategy));
			printf("\t\t\t\tcorrect.hyperpar=[%1d]\n", mb->gcpo_param->correct_hyperpar);
			printf("\t\t\t\tepsilon=[%g]\n", mb->gcpo_param->epsilon);
			printf("\t\t\t\tprior.diagonal=[%g]\n", mb->gcpo_param->prior_diagonal);
			if (mb->gcpo_param->groups) {
				printf("\t\t\t\tUse user-defined gcpo-groups, ngroups.eff=[%1d]\n", ngroups_eff);
			}
			if (mb->gcpo_param->selection) {
				printf("\t\t\t\tUse user-defined selection, nselection=[%1d]\n", mb->gcpo_param->selection->n);
			}
			if (mb->gcpo_param->friends) {
				printf("\t\t\t\tUse friends-list, n=[%1d]\n", mb->gcpo_param->friends_n);
			}

			if (mb->gcpo_param->keep) {
				printf("\t\t\t\tkeep=[");
				for (i = 0; i < mb->gcpo_param->keep->n; i++) {
					if (i) {
						printf(" ");
					}
					printf("%s", mb->gcpo_param->keep->str[i]);
				}
				printf("]\n");
			} else {
				printf("\t\t\t\tkeep=[]\n");
			}

			printf("\t\t\t\tremove.fixed=[%1d]\n", mb->gcpo_param->remove_fixed);
			if (mb->gcpo_param->remove) {
				printf("\t\t\t\tremove=[");
				for (i = 0; i < mb->gcpo_param->remove->n; i++) {
					if (i) {
						printf(" ");
					}
					printf("%s", mb->gcpo_param->remove->str[i]);
				}
				printf("]\n");
			} else {
				printf("\t\t\t\tremove=[]\n");
			}

			printf("\t\t\tcpo=[%1d]\n", (*out)->cpo);
			printf("\t\t\tpo=[%1d]\n", (*out)->po);
			printf("\t\t\tdic=[%1d]\n", (*out)->dic);
			printf("\t\t\tkld=[%1d]\n", (*out)->kld);
			printf("\t\t\tmlik=[%1d]\n", (*out)->mlik);
			printf("\t\t\tq=[%1d]\n", (*out)->q);
			printf("\t\t\tgraph=[%1d]\n", (*out)->graph);
			printf("\t\t\thyperparameters=[%1d]\n", (*out)->hyperparameters);
			printf("\t\t\tconfig=[%1d]\n", (*out)->config);
			printf("\t\t\tlikelihood.info=[%1d]\n", (*out)->likelihood_info);
			printf("\t\t\tinternal.opt=[%1d]\n", (*out)->internal_opt);
			printf("\t\t\tsave.memory=[%1d]\n", (*out)->save_memory);
		}
		printf("\t\t\tsummary=[%1d]\n", (*out)->summary);
		printf("\t\t\treturn.marginals=[%1d]\n", (*out)->return_marginals);
		printf("\t\t\treturn.marginals.predictor=[%1d]\n", (*out)->return_marginals_predictor);
		printf("\t\t\tnquantiles=[%1d]  [", (*out)->nquantiles);
		for (i = 0; i < (*out)->nquantiles; i++) {
			printf(" %g", (*out)->quantiles[i]);
		}
		printf(" ]\n");
		printf("\t\t\tncdf=[%1d]  [", (*out)->ncdf);
		for (i = 0; i < (*out)->ncdf; i++) {
			printf(" %g", (*out)->cdf[i]);
		}
		printf(" ]\n");
	}
	return INLA_OK;
}

int inla_R(char **argv)
{
	while (*argv) {
		fprintf(stderr, "R> source file [%s]...\n", *argv);
		inla_R_source(*argv);
		argv++;
	}
	exit(0);

	return GMRFLib_SUCCESS;
}

