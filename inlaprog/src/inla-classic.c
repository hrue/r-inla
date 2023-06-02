
/* inla-classic.c
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

int inla_INLA(inla_tp *mb)
{
	double *c = NULL, *x = NULL, *b = NULL;
	int N, i, j, k, count, local_count;
	char *compute = NULL;
	GMRFLib_bfunc_tp **bfunc;

	if (mb->verbose) {
		printf("%s...\n", __GMRFLib_FuncName);
	}

	/*
	 * We need to determine the strategy if strategy is default 
	 */
	GMRFLib_density_storage_strategy_tp storage_scheme = GMRFLib_DENSITY_STORAGE_STRATEGY_HIGH;
	int ntot = mb->predictor_n + mb->predictor_m + mb->nlinear;

	for (i = 0; i < mb->nf; i++) {
		ntot += mb->f_graph[i]->n;
	}
	if (ntot < 50000) {
		storage_scheme = GMRFLib_DENSITY_STORAGE_STRATEGY_HIGH;
	} else {
		storage_scheme = GMRFLib_DENSITY_STORAGE_STRATEGY_LOW;
	}

	if (mb->strategy == GMRFLib_OPENMP_STRATEGY_DEFAULT) {
		/*
		 * to determine the strategy, count the size of the model 
		 */
		if (mb->verbose) {
			printf("\tStrategy = [DEFAULT]\n");
		}
		if (ntot < 500) {
			mb->strategy = GMRFLib_OPENMP_STRATEGY_SMALL;
		} else if (ntot < 2000) {
			mb->strategy = GMRFLib_OPENMP_STRATEGY_MEDIUM;
		} else if (ntot < 50000) {
			mb->strategy = GMRFLib_OPENMP_STRATEGY_LARGE;
		} else {
			mb->strategy = GMRFLib_OPENMP_STRATEGY_HUGE;
		}
	}

	GMRFLib_density_storage_strategy = storage_scheme;
	GMRFLib_openmp->strategy = mb->strategy;

	if (mb->verbose) {
		printf("\tMode..................... [%s]\n", GMRFLib_MODE_NAME());
		printf("\tSparse-matrix library.... [%s]\n", mb->smtp);
		printf("\tOpenMP strategy.......... [%s]\n", GMRFLib_OPENMP_STRATEGY_NAME(GMRFLib_openmp->strategy));
		printf("\tnum.threads.............. [%1d:%1d]\n", GMRFLib_openmp->max_threads_nested[0], GMRFLib_openmp->max_threads_nested[1]);
		if (GMRFLib_openmp->adaptive) {
			printf("\tnum.threads (adaptive)... [%1d]\n", GMRFLib_PARDISO_MAX_NUM_THREADS());
		}
		printf("\tblas.num.threads......... [%1d]\n", GMRFLib_openmp->blas_num_threads);
		printf("\tDensity-strategy......... [%s]\n",
		       (GMRFLib_density_storage_strategy == GMRFLib_DENSITY_STORAGE_STRATEGY_LOW ? "Low" : "High"));
	}

	GMRFLib_openmp_implement_strategy(GMRFLib_OPENMP_PLACES_BUILD_MODEL, NULL, NULL);

	GMRFLib_init_hgmrfm(&(mb->hgmrfm), mb->predictor_n, mb->predictor_m,
			    mb->predictor_cross_sumzero, mb->predictor_log_prec,
			    (const char *) mb->predictor_Aext_fnm, mb->predictor_Aext_precision,
			    mb->nf, mb->f_c, mb->f_weights, mb->f_graph, mb->f_Qfunc, mb->f_Qfunc_arg, mb->f_sumzero, mb->f_constr,
			    mb->ff_Qfunc, mb->ff_Qfunc_arg, mb->nlinear, mb->linear_covariate, mb->linear_precision, mb->ai_par);
	GMRFLib_openmp_implement_strategy(GMRFLib_OPENMP_PLACES_OPTIMIZE, NULL, NULL);
	N = ((GMRFLib_hgmrfm_arg_tp *) mb->hgmrfm->Qfunc_arg)->N;
	if (mb->verbose) {
		printf("\tSize of graph............ [%d]\n", N);
		printf("\tNumber of constraints.... [%d]\n", (mb->hgmrfm->constr ? mb->hgmrfm->constr->nc : 0));
	}

	mb->d = Realloc(mb->d, N, double);
	Memset(&(mb->d[mb->predictor_ndata]), 0, (N - mb->predictor_ndata) * sizeof(double));
	mb->loglikelihood = Realloc(mb->loglikelihood, N, GMRFLib_logl_tp *);
	Memset(&(mb->loglikelihood[mb->predictor_ndata]), 0, (N - mb->predictor_ndata) * sizeof(GMRFLib_logl_tp *));
	mb->loglikelihood_arg = Realloc(mb->loglikelihood_arg, N, void *);
	Memset(&(mb->loglikelihood_arg[mb->predictor_ndata]), 0, (N - mb->predictor_ndata) * sizeof(void *));

	if (0) {
		for (i = 0; i < N; i++)
			printf("d[%d]=%g\n", i, mb->d[i]);
	}

	/*
	 * add the diagonal, if any 
	 */

	c = Calloc(N, double);
	count = mb->predictor_n + mb->predictor_m;
	for (i = 0; i < mb->nf; i++) {
		for (k = 0; k < mb->f_nrep[i]; k++) {
			for (j = 0; j < mb->f_n[i]; j++) {
				c[count + j + k * mb->f_N[i]] = mb->f_diag[i];	/* yes; this is correct */
			}
		}
		count += mb->f_Ntotal[i];		       /* yes; this is correct */
	}

	/*
	 * this is an emergency option to prevent singular matrices (and is known to be >= 0) 
	 */
	if (mb->expert_diagonal_emergencey) {
		for (i = mb->predictor_n + mb->predictor_m; i < N; i++)
			c[i] += mb->expert_diagonal_emergencey;
	}

	if (0) {
		for (i = 0; i < N; i++)
			printf("c[%d]=%g\n", i, c[i]);
	}

	/*
	 * mark those we want to compute and compute the b
	 */
	compute = Calloc(N, char);
	b = Calloc(N, double);
	bfunc = Calloc(N, GMRFLib_bfunc_tp *);
	count = 0;
	if (mb->expert_cpo_manual) {
		/*
		 * if set, then only then only `linear.predictor[idx]' is set
		 */
		for (i = 0; i < mb->predictor_n + mb->predictor_m; i++) {
			compute[count] = (char) 0;
			count++;
		}

		for (i = 0; i < mb->expert_n_cpo_idx; i++) {
			compute[mb->expert_cpo_idx[i]] = (char) 1;
			mb->d[mb->expert_cpo_idx[i]] = 0.0;
		}
		mb->ai_par->cpo_manual = 1;
		mb->output->hyperparameters = GMRFLib_FALSE;

		for (i = 0; i < mb->nf; i++) {
			if (mb->f_bfunc2[i]) {
				for (j = 0; j < mb->f_Ntotal[i]; j++) {
					bfunc[count + j] = Calloc(1, GMRFLib_bfunc_tp);
					bfunc[count + j]->bdef = mb->f_bfunc2[i];
					bfunc[count + j]->idx = j;
				}
			}
			for (j = 0; j < mb->f_Ntotal[i]; j++) {
				compute[count] = (char) 0;
				count++;
			}
		}

		for (i = 0; i < mb->nlinear; i++) {
			compute[count] = (char) 0;
			b[count] = mb->linear_precision[i] * mb->linear_mean[i];
			count++;
		}
		assert(count == N);
	} else {
		/*
		 * as before 
		 */
		for (i = 0; i < mb->predictor_n + mb->predictor_m; i++) {
			compute[count++] = (char) mb->predictor_compute;
		}
		for (i = 0; i < mb->nf; i++) {
			if (mb->f_bfunc2[i]) {
				for (j = 0; j < mb->f_Ntotal[i]; j++) {
					bfunc[count + j] = Calloc(1, GMRFLib_bfunc_tp);
					bfunc[count + j]->bdef = mb->f_bfunc2[i];
					bfunc[count + j]->idx = j;
				}
			}
			for (j = 0; j < mb->f_Ntotal[i]; j++) {
				compute[count] = (char) mb->f_compute[i];
				count++;
			}
		}
		for (i = 0; i < mb->nlinear; i++) {
			compute[count] = (char) mb->linear_compute[i];
			b[count] = mb->linear_precision[i] * mb->linear_mean[i];
			count++;
		}
		if (count != N) {
			P(count);
			P(N);
			assert(count == N);
		}
	}

	// VB correct 
	char *vb_nodes = NULL;

	local_count = 0;
	if (mb->ai_par->vb_enable) {
		vb_nodes = Calloc(N, char);
		count = mb->predictor_n + mb->predictor_m;
		for (i = 0; i < mb->nf; i++) {
			GMRFLib_idx_tp *vb = mb->f_vb_correct[i];
			if ((vb->idx[0] == -1L && mb->f_Ntotal[i] <= mb->ai_par->vb_f_enable_limit_mean)) {
				for (j = 0; j < mb->f_Ntotal[i]; j++) {
					vb_nodes[count + j] = (char) 1;
					local_count++;
				}
			} else if (vb->idx[0] == -1L) {
				int len, jj;
				len = IMAX(1, mb->f_Ntotal[i] / mb->ai_par->vb_f_enable_limit_mean);	/* integer division */
				k = IMAX(1, len / 2);	       /* integer division */
				for (j = 0; j < mb->ai_par->vb_f_enable_limit_mean; j++) {
					jj = (j * len + k) % mb->f_Ntotal[i];
					vb_nodes[count + jj] = (char) 1;
					local_count++;
				}
			} else if (vb->idx[0] >= 0) {
				for (j = 0; j < vb->n; j++) {
					vb_nodes[count + IMIN(vb->idx[j], mb->f_Ntotal[i] - 1)] = (char) 1;
					local_count++;
				}
			}
			count += mb->f_Ntotal[i];
		}

		for (i = 0; i < mb->nlinear; i++) {
			vb_nodes[count++] = (char) 1;
			local_count++;
		}
		if (local_count == 0) {			       /* then there is nothting to correct for */
			Free(vb_nodes);
			vb_nodes = NULL;
		}
	}
	mb->ai_par->vb_nodes_mean = vb_nodes;

	local_count = 0;
	if (mb->ai_par->vb_enable) {
		vb_nodes = Calloc(N, char);
		count = mb->predictor_n + mb->predictor_m;
		for (i = 0; i < mb->nf; i++) {
			GMRFLib_idx_tp *vb = mb->f_vb_correct[i];
			if ((vb->idx[0] == -1L && mb->f_Ntotal[i] <= mb->ai_par->vb_f_enable_limit_variance)) {
				for (j = 0; j < mb->f_Ntotal[i]; j++) {
					vb_nodes[count + j] = (char) 1;
					local_count++;
				}
			} else if (vb->idx[0] == -1L) {
				int len, jj;
				len = IMAX(1, mb->f_Ntotal[i] / mb->ai_par->vb_f_enable_limit_variance);	/* integer division */
				k = IMAX(1, len / 2);	       /* integer division */
				for (j = 0; j < mb->ai_par->vb_f_enable_limit_variance; j++) {
					jj = (j * len + k) % mb->f_Ntotal[i];
					vb_nodes[count + jj] = (char) 1;
					local_count++;
				}
			} else if (vb->idx[0] >= 0) {
				for (j = 0; j < vb->n; j++) {
					vb_nodes[count + IMIN(vb->idx[j], mb->f_Ntotal[i] - 1)] = (char) 1;
					local_count++;
				}
			}

			count += mb->f_Ntotal[i];
		}
		for (i = 0; i < mb->nlinear; i++) {
			vb_nodes[count++] = (char) 1;
			local_count++;
		}
		if (local_count == 0) {			       /* then there is nothting to correct for */
			Free(vb_nodes);
			vb_nodes = NULL;
		}
	}
	mb->ai_par->vb_nodes_variance = vb_nodes;

	// define the adaptive strategy
	GMRFLib_ai_strategy_tp *adapt = NULL;
	if (mb->ai_par->strategy == GMRFLib_AI_STRATEGY_ADAPTIVE) {
		adapt = Calloc(N, GMRFLib_ai_strategy_tp);
		for (i = 0; i < N; i++) {
			adapt[i] = GMRFLib_AI_STRATEGY_GAUSSIAN;
		}
		count = mb->predictor_n + mb->predictor_m;
		for (i = 0; i < mb->nf; i++) {
			if (mb->f_Ntotal[i] <= mb->ai_par->adapt_max) {
				/*
				 * add also random effects with small size
				 */
				for (j = 0; j < mb->f_Ntotal[i]; j++) {
					adapt[count + j] = GMRFLib_AI_STRATEGY_MEANSKEWCORRECTED_GAUSSIAN;
				}
			}
			count += mb->f_Ntotal[i];
		}
		for (i = 0; i < mb->nlinear; i++) {
			adapt[count++] = GMRFLib_AI_STRATEGY_MEANSKEWCORRECTED_GAUSSIAN;
		}
	}
	mb->ai_par->adapt_strategy = adapt;
	mb->ai_par->adapt_len = (adapt ? N : 0);

	if (G.reorder < 0) {
		size_t nnz = 0;
		int use_g = 0;
		GMRFLib_optimize_reorder(mb->hgmrfm->graph, &nnz, &use_g, &(mb->gn));
		if (GMRFLib_smtp != GMRFLib_SMTP_PARDISO) {
			// ....
		} else {
			GMRFLib_reorder = GMRFLib_REORDER_PARDISO;
		}
		if (GMRFLib_smtp != GMRFLib_SMTP_PARDISO) {
			if (mb->verbose) {
				printf("\tFound optimal reordering=[%s] nnz(L)=[%zu] and use_global_nodes(user)=[%s]\n",
				       GMRFLib_reorder_name(GMRFLib_reorder), nnz, (use_g ? "yes" : "no"));
			}
		}
	}
	if (mb->verbose) {
		if (mb->ntheta) {
			printf("\tList of hyperparameters: \n");
			for (i = 0; i < mb->ntheta; i++) {
				printf("\t\ttheta[%1d] = [%s]\n", i, mb->theta_tag[i]);
			}
		} else {
			printf("\tNone hyperparameters\n");
		}
		printf("\n");
	}
	GMRFLib_ai_store_tp *ai_store = Calloc(1, GMRFLib_ai_store_tp);

	if (mb->output->dic) {
		mb->dic = Calloc(1, GMRFLib_ai_dic_tp);
	} else {
		mb->dic = NULL;
	}
	/*
	 * compute a 'reasonable' initial value for \eta, unless its there from before.
	 */

	int mm = mb->predictor_n + mb->predictor_m;
	Free(G_norm_const_compute);
	Free(G_norm_const);
	for (i = 0; i < G_norm_const_len; i++) {
		Free(G_norm_const_v[i]);
	}
	Free(G_norm_const_v);
	G_norm_const_len = mm;
	G_norm_const_compute = Calloc(mm, char);
	G_norm_const = Calloc(mm, double);
	G_norm_const_v = Calloc(mm, void *);
	for (i = 0; i < mm; i++) {
		G_norm_const[i] = NAN;
		G_norm_const_compute[i] = 1;
	}

	x = Calloc(N, double);
	if (mb->reuse_mode && mb->x_file) {
		if (N != mb->nx_file) {
			char *msg;
			GMRFLib_sprintf(&msg, "N = %1d but nx_file = %1d. Stop.", N, mb->nx_file);
			inla_error_general(msg);
		}
		Memcpy(x, mb->x_file, N * sizeof(double));

		/*
		 * subtract the offset 
		 */
		for (i = 0; i < mb->predictor_ndata; i++) {
			x[i] -= OFFSET3(i);
		}

	} else {
#pragma omp parallel for private(i) num_threads(GMRFLib_openmp->max_threads_outer)
		for (i = 0; i < mb->predictor_ndata; i++) {
			if (mb->d[i]) {
				x[i] = inla_compute_initial_value(i, mb->loglikelihood[i], x, (void *) mb->loglikelihood_arg[i]);
			} else {
				x[i] = 0.0;
			}
		}
	}

	/*
	 * set the flag to compute correlation-matrix or not 
	 */
	mb->misc_output = Calloc(1, GMRFLib_ai_misc_output_tp);
	if (mb->lc_derived_correlation_matrix) {
		mb->misc_output->compute_corr_lin = mb->nlc;   /* yes, pass the dimension */
	} else {
		mb->misc_output->compute_corr_lin = 0;
	}
	if (mb->output->config) {
		mb->misc_output->configs = Calloc(GMRFLib_MAX_THREADS(), GMRFLib_store_configs_tp *);
	} else {
		mb->misc_output->configs = NULL;
	}

	if (mb->fixed_mode) {
		/*
		 * then there is a request to treat the theta's as fixed and known. This little hack do the job nicely. 
		 */
		mb->ntheta = 0;
		mb->data_ntheta_all = 0;
		mb->theta = NULL;
	}

	mb->transform_funcs = Calloc(N, GMRFLib_transform_array_func_tp *);
	for (i = 0; i < mb->predictor_m + mb->predictor_n; i++) {
		/*
		 * only where we have data (ie n or m), can be a invlinkfunc different from the identity. also the compute-flag must be ON.
		 */
		if (!compute[i]) {
			mb->transform_funcs[i] = NULL;
		} else if (i < mb->predictor_ndata && mb->predictor_invlinkfunc[i]) {
			mb->transform_funcs[i] = Calloc(1, GMRFLib_transform_array_func_tp);
			mb->transform_funcs[i]->func = (GMRFLib_transform_func_tp *) mb->predictor_invlinkfunc[i];
			mb->transform_funcs[i]->arg = mb->predictor_invlinkfunc_arg[i];

			double *cov = NULL;
			if (mb->predictor_invlinkfunc_covariates && mb->predictor_invlinkfunc_covariates[i]) {
				int ncov = mb->predictor_invlinkfunc_covariates[i]->ncol;
				cov = Calloc(ncov, double);
				GMRFLib_matrix_get_row(cov, i, mb->predictor_invlinkfunc_covariates[i]);
			}
			mb->transform_funcs[i]->cov = cov;     /* yes, we store a copy here */
		} else {
			mb->transform_funcs[i] = Calloc(1, GMRFLib_transform_array_func_tp);
			mb->transform_funcs[i]->func = (GMRFLib_transform_func_tp *) link_identity;
			mb->transform_funcs[i]->arg = NULL;
			mb->transform_funcs[i]->cov = NULL;
		}
	}

	/*
	 * If Gaussian data, then force the strategy to be Gaussian  
	 */
	if (mb->gaussian_data) {
		mb->ai_par->strategy = GMRFLib_AI_STRATEGY_GAUSSIAN;
		mb->ai_par->gaussian_data = mb->gaussian_data;
	}

	/*
	 * Finally, let us do the job...
	 */
	GMRFLib_ai_INLA(&(mb->density),
			// DISABLE THIS FEATURE NOW, IT DOES NOT WORK WELL ENOGUH
			(0 ? &(mb->density_transform) : NULL), (0 ? mb->transform_funcs : NULL),
			// ....
			(mb->output->hyperparameters ? &(mb->density_hyper) : NULL),
			(mb->output->cpo || mb->expert_cpo_manual ? &(mb->cpo) : NULL),
			(mb->output->po ? &(mb->po) : NULL),
			mb->dic,
			(mb->output->mlik ? &(mb->mlik) : NULL),
			compute, mb->theta, mb->ntheta,
			extra, (void *) mb,
			x, b, c, NULL, bfunc, mb->d,
			loglikelihood_inla, (void *) mb,
			mb->hgmrfm->graph, mb->hgmrfm->Qfunc, mb->hgmrfm->Qfunc_arg, mb->hgmrfm->constr, mb->ai_par, ai_store,
			mb->nlc, mb->lc_lc, &(mb->density_lin), mb->misc_output, NULL, NULL);

	/*
	 * add the offsets to the linear predictor. Add the offsets to the 'configs' (if any), at a later stage. 
	 */
#pragma omp parallel for private(i) num_threads(GMRFLib_openmp->max_threads_outer)
	for (i = 0; i < mb->predictor_n + mb->predictor_m; i++) {
		GMRFLib_density_tp *d;
		if (mb->density[i]) {
			d = mb->density[i];
			GMRFLib_density_new_mean(&(mb->density[i]), d, d->std_mean + OFFSET3(i));
			GMRFLib_free_density(d);
		}
	}

	/*
	 * add the offset to 'x' 
	 */
	for (i = 0; i < mb->predictor_n + mb->predictor_m; i++) {
		x[i] += OFFSET3(i);
	}

	Free(mb->x_file);				       /* yes, and then */
	mb->x_file = x;					       /* just take over */
	mb->nx_file = N;

	GMRFLib_free_ai_store(ai_store);
	Free(b);
	Free(c);
	Free(compute);

	return INLA_OK;
}

int inla_INLA_preopt_stage1(inla_tp *mb, GMRFLib_preopt_res_tp *rpreopt)
{
	double *c = NULL, *x = NULL, *b = NULL;
	int N, i, j, count;
	char *compute = NULL;
	GMRFLib_bfunc_tp **bfunc;
	GMRFLib_preopt_tp *preopt = NULL;

	if (mb->verbose) {
		printf("%s...\n", __GMRFLib_FuncName);
	}

	/*
	 * We need to determine the strategy if strategy is default 
	 */
	GMRFLib_density_storage_strategy_tp storage_scheme = GMRFLib_DENSITY_STORAGE_STRATEGY_HIGH;
	int ntot = 0;

	ntot = mb->nlinear;
	for (i = 0; i < mb->nf; i++) {
		ntot += mb->f_graph[i]->n;
	}
	N = ntot;

	if (mb->strategy == GMRFLib_OPENMP_STRATEGY_DEFAULT) {
		/*
		 * to determine the strategy, count the size of the model 
		 */
		if (mb->verbose) {
			printf("\tStrategy = [DEFAULT]\n");
		}
		if (ntot < 500) {
			mb->strategy = GMRFLib_OPENMP_STRATEGY_SMALL;
		} else if (ntot < 2000) {
			mb->strategy = GMRFLib_OPENMP_STRATEGY_MEDIUM;
		} else if (ntot < 50000) {
			mb->strategy = GMRFLib_OPENMP_STRATEGY_LARGE;
		} else {
			mb->strategy = GMRFLib_OPENMP_STRATEGY_HUGE;
		}
	}

	GMRFLib_density_storage_strategy = storage_scheme;
	GMRFLib_openmp->strategy = mb->strategy;

	b = Calloc(N, double);
	bfunc = Calloc(N, GMRFLib_bfunc_tp *);
	for (count = 0, i = 0; i < mb->nf; i++) {
		if (mb->f_bfunc2[i]) {
			for (j = 0; j < mb->f_Ntotal[i]; j++) {
				bfunc[count + j] = Calloc(1, GMRFLib_bfunc_tp);
				bfunc[count + j]->bdef = mb->f_bfunc2[i];
				bfunc[count + j]->idx = j;
			}
		}
		count += mb->f_Ntotal[i];
	}

	// VB correct 
	char *vb_nodes = NULL;
	int local_count = 0;
	if (mb->ai_par->vb_enable) {
		vb_nodes = Calloc(N, char);
		count = 0;
		for (i = 0; i < mb->nf; i++) {
			GMRFLib_idx_tp *vb = mb->f_vb_correct[i];
			if ((vb->idx[0] == -1L && mb->f_Ntotal[i] <= mb->ai_par->vb_f_enable_limit_mean)) {
				for (j = 0; j < mb->f_Ntotal[i]; j++) {
					vb_nodes[count + j] = (char) 1;
					local_count++;
				}
			} else if (vb->idx[0] == -1L) {
				int len, jj;
				len = IMAX(1, mb->f_Ntotal[i] / mb->ai_par->vb_f_enable_limit_mean);	/* integer division */
				int k = IMAX(1, len / 2);      /* integer division */
				for (j = 0; j < mb->ai_par->vb_f_enable_limit_mean; j++) {
					jj = (j * len + k) % mb->f_Ntotal[i];
					vb_nodes[count + jj] = (char) 1;
					local_count++;
				}
			} else if (vb->idx[0] >= 0) {
				for (j = 0; j < vb->n; j++) {
					vb_nodes[count + IMIN(vb->idx[j], mb->f_Ntotal[i] - 1)] = (char) 1;
					local_count++;
				}
			}
			count += mb->f_Ntotal[i];
		}
		for (i = 0; i < mb->nlinear; i++) {
			vb_nodes[count++] = (char) 1;
			local_count++;
		}
		if (local_count == 0) {			       /* then there is nothting to correct for */
			Free(vb_nodes);
			vb_nodes = NULL;
		}
	}
	mb->ai_par->vb_nodes_mean = vb_nodes;

	double tref = GMRFLib_cpu();
	GMRFLib_openmp_implement_strategy(GMRFLib_OPENMP_PLACES_GCPO_BUILD, NULL, NULL);
	GMRFLib_preopt_init(&preopt,
			    mb->predictor_n, mb->nf, mb->f_c, mb->f_weights,
			    mb->f_graph, mb->f_Qfunc, mb->f_Qfunc_arg, mb->f_sumzero, mb->f_constr,
			    mb->f_diag,
			    mb->ff_Qfunc, mb->ff_Qfunc_arg,
			    mb->nlinear, mb->linear_covariate, mb->linear_precision, bfunc, mb->ai_par, mb->predictor_A_fnm, mb->global_constr);
	mb->preopt = preopt;
	assert(preopt->latent_graph->n == N);

	if (mb->verbose) {
		printf("\tMode..................... [%s]\n", GMRFLib_MODE_NAME());
		printf("\tSetup.................... [%.2fs]\n", GMRFLib_cpu() - tref);
		printf("\tSparse-matrix library.... [%s]\n", mb->smtp);
		printf("\tOpenMP strategy.......... [%s]\n", GMRFLib_OPENMP_STRATEGY_NAME(GMRFLib_openmp->strategy));
		printf("\tnum.threads.............. [%1d:%1d]\n", GMRFLib_openmp->max_threads_nested[0], GMRFLib_openmp->max_threads_nested[1]);
		if (GMRFLib_openmp->adaptive) {
			printf("\tnum.threads (adaptive)... [%1d]\n", GMRFLib_PARDISO_MAX_NUM_THREADS());
		}
		printf("\tblas.num.threads......... [%1d]\n", GMRFLib_openmp->blas_num_threads);
		printf("\tDensity-strategy......... [%s]\n",
		       (GMRFLib_density_storage_strategy == GMRFLib_DENSITY_STORAGE_STRATEGY_LOW ? "Low" : "High"));
		printf("\tSize of graph............ [%d]\n", N);
		printf("\tNumber of constraints.... [%d]\n", (preopt->latent_constr ? preopt->latent_constr->nc : 0));
	}
	GMRFLib_openmp_implement_strategy(GMRFLib_OPENMP_PLACES_OPTIMIZE, NULL, NULL);

	c = Calloc(N, double);
	if (mb->expert_diagonal_emergencey) {
		for (i = 0; i < N; i++)
			c[i] += mb->expert_diagonal_emergencey;
	}

	if (G.reorder < 0) {
		size_t nnz = 0;
		int use_g = 0;
		GMRFLib_optimize_reorder(preopt->latent_graph, &nnz, &use_g, &(mb->gn));
		if (GMRFLib_smtp != GMRFLib_SMTP_PARDISO) {
			// ....
		} else {
			GMRFLib_reorder = GMRFLib_REORDER_PARDISO;
		}
		if (GMRFLib_smtp != GMRFLib_SMTP_PARDISO) {
			if (mb->verbose) {
				printf("\tFound optimal reordering=[%s] nnz(L)=[%zu] and use_global_nodes(user)=[%s]\n",
				       GMRFLib_reorder_name(GMRFLib_reorder), nnz, (use_g ? "yes" : "no"));
			}
		}
	}
	/*
	 * mark those we want to compute and compute the b
	 */
	if (mb->verbose) {
		if (mb->ntheta) {
			printf("\tList of hyperparameters: \n");
			for (i = 0; i < mb->ntheta; i++) {
				printf("\t\ttheta[%1d] = [%s]\n", i, mb->theta_tag[i]);
			}
		} else {
			printf("\tNone hyperparameters\n");
		}
		printf("\n");
	}
	GMRFLib_ai_store_tp *ai_store = Calloc(1, GMRFLib_ai_store_tp);

	mb->dic = NULL;
	mb->misc_output = Calloc(1, GMRFLib_ai_misc_output_tp);
	x = Calloc(N, double);
	if (mb->reuse_mode && mb->x_file) {
		Memcpy(x, mb->x_file + preopt->mnpred, N * sizeof(double));
	}

	int nparam_eff = mb->ai_par->compute_nparam_eff;
	mb->ai_par->compute_nparam_eff = 0;
	compute = Calloc(N, char);

	Free(G_norm_const_compute);
	Free(G_norm_const);
	for (i = 0; i < G_norm_const_len; i++) {
		Free(G_norm_const_v[i]);
	}
	Free(G_norm_const_v);
	G_norm_const_len = preopt->Npred;
	G_norm_const_compute = Calloc(preopt->Npred, char);
	G_norm_const = Calloc(preopt->Npred, double);
	G_norm_const_v = Calloc(preopt->Npred, void *);
	for (i = 0; i < preopt->Npred; i++) {
		G_norm_const[i] = NAN;
		G_norm_const_compute[i] = 1;
	}

	GMRFLib_ai_INLA(&(mb->density),
			NULL, NULL,
			(mb->output->hyperparameters ? &(mb->density_hyper) : NULL),
			NULL,
			NULL,
			NULL,
			(mb->output->mlik ? &(mb->mlik) : NULL),
			compute, mb->theta, mb->ntheta,
			extra, (void *) mb,
			x, b, c, NULL, bfunc, mb->d,
			loglikelihood_inla, (void *) mb,
			preopt->preopt_graph, preopt->preopt_Qfunc, preopt->preopt_Qfunc_arg, preopt->latent_constr,
			mb->ai_par, ai_store, 0, NULL, NULL, mb->misc_output, preopt, rpreopt);

	if (rpreopt->int_design) {
		mb->ai_par->int_strategy = GMRFLib_AI_INT_STRATEGY_USER_EXPERT;
		GMRFLib_design_tp *design = NULL;
		GMRFLib_design_read(&design, rpreopt->int_design, 0);

		for (i = 0; i < design->nexperiments; i++) {
			design->int_weight[i] = rpreopt->adj_weights[i];
		}
		mb->ai_par->adjust_weights = 0;

		if (mb->verbose) {
			GMRFLib_design_print(stdout, design);
		}
		if (mb->verbose) {
			printf("\tPrune design: prob=0.95\n");
		}
		GMRFLib_design_prune(design, 0.95);
		if (mb->verbose) {
			GMRFLib_design_print(stdout, design);
		}
		mb->ai_par->int_design = design;
	}

	// set back
	nparam_eff = mb->ai_par->compute_nparam_eff = nparam_eff;

	GMRFLib_free_ai_store(ai_store);
	Free(x);
	Free(b);
	Free(c);
	Free(compute);
	Free(bfunc);

	return INLA_OK;
}

int inla_INLA_preopt_stage2(inla_tp *mb, GMRFLib_preopt_res_tp *rpreopt)
{
	double *c = NULL, *x = NULL, *b = NULL;
	int N, i, j, k, count, local_count;
	char *compute = NULL;
	GMRFLib_bfunc_tp **bfunc;

	if (mb->verbose) {
		printf("%s...\n", __GMRFLib_FuncName);
	}

	/*
	 * We need to determine the strategy if strategy is default 
	 */
	GMRFLib_density_storage_strategy_tp storage_scheme = GMRFLib_DENSITY_STORAGE_STRATEGY_HIGH;
	int ntot = mb->predictor_n + mb->predictor_m + mb->nlinear;

	for (i = 0; i < mb->nf; i++) {
		ntot += mb->f_graph[i]->n;
	}
	if (ntot < 50000) {
		storage_scheme = GMRFLib_DENSITY_STORAGE_STRATEGY_HIGH;
	} else {
		storage_scheme = GMRFLib_DENSITY_STORAGE_STRATEGY_LOW;
	}

	if (mb->strategy == GMRFLib_OPENMP_STRATEGY_DEFAULT) {
		/*
		 * to determine the strategy, count the size of the model 
		 */
		if (mb->verbose) {
			printf("\tStrategy = [DEFAULT]\n");
		}
		if (ntot < 500) {
			mb->strategy = GMRFLib_OPENMP_STRATEGY_SMALL;
		} else if (ntot < 2000) {
			mb->strategy = GMRFLib_OPENMP_STRATEGY_MEDIUM;
		} else if (ntot < 50000) {
			mb->strategy = GMRFLib_OPENMP_STRATEGY_LARGE;
		} else {
			mb->strategy = GMRFLib_OPENMP_STRATEGY_HUGE;
		}
	}

	GMRFLib_density_storage_strategy = storage_scheme;
	GMRFLib_openmp->strategy = mb->strategy;

	if (mb->verbose) {
		printf("\tMode..................... [%s]\n", GMRFLib_MODE_NAME());
		printf("\tSparse-matrix library.... [%s]\n", mb->smtp);
		printf("\tOpenMP strategy.......... [%s]\n", GMRFLib_OPENMP_STRATEGY_NAME(GMRFLib_openmp->strategy));
		printf("\tnum.threads.............. [%1d:%1d]\n", GMRFLib_openmp->max_threads_nested[0], GMRFLib_openmp->max_threads_nested[1]);
		if (GMRFLib_openmp->adaptive) {
			printf("\tnum.threads (adaptive)... [%1d]\n", GMRFLib_PARDISO_MAX_NUM_THREADS());
		}
		printf("\tblas.num.threads......... [%1d]\n", GMRFLib_openmp->blas_num_threads);
		printf("\tDensity-strategy......... [%s]\n",
		       (GMRFLib_density_storage_strategy == GMRFLib_DENSITY_STORAGE_STRATEGY_LOW ? "Low" : "High"));
	}

	GMRFLib_openmp_implement_strategy(GMRFLib_OPENMP_PLACES_BUILD_MODEL, NULL, NULL);

	GMRFLib_init_hgmrfm(&(mb->hgmrfm), mb->predictor_n, mb->predictor_m,
			    mb->predictor_cross_sumzero, mb->predictor_log_prec,
			    (const char *) mb->predictor_Aext_fnm, mb->predictor_Aext_precision,
			    mb->nf, mb->f_c, mb->f_weights, mb->f_graph, mb->f_Qfunc, mb->f_Qfunc_arg, mb->f_sumzero, mb->f_constr,
			    mb->ff_Qfunc, mb->ff_Qfunc_arg, mb->nlinear, mb->linear_covariate, mb->linear_precision, mb->ai_par);
	GMRFLib_openmp_implement_strategy(GMRFLib_OPENMP_PLACES_OPTIMIZE, NULL, NULL);
	N = ((GMRFLib_hgmrfm_arg_tp *) mb->hgmrfm->Qfunc_arg)->N;
	if (mb->verbose) {
		printf("\tSize of graph............ [%d]\n", N);
		printf("\tNumber of constraints.... [%d]\n", (mb->hgmrfm->constr ? mb->hgmrfm->constr->nc : 0));
	}

	mb->d = Realloc(mb->d, N, double);
	Memset(&(mb->d[mb->predictor_ndata]), 0, (N - mb->predictor_ndata) * sizeof(double));
	mb->loglikelihood = Realloc(mb->loglikelihood, N, GMRFLib_logl_tp *);
	Memset(&(mb->loglikelihood[mb->predictor_ndata]), 0, (N - mb->predictor_ndata) * sizeof(GMRFLib_logl_tp *));
	mb->loglikelihood_arg = Realloc(mb->loglikelihood_arg, N, void *);
	Memset(&(mb->loglikelihood_arg[mb->predictor_ndata]), 0, (N - mb->predictor_ndata) * sizeof(void *));

	if (0) {
		for (i = 0; i < N; i++)
			printf("d[%d]=%g\n", i, mb->d[i]);
	}

	/*
	 * add the diagonal, if any 
	 */

	c = Calloc(N, double);
	count = mb->predictor_n + mb->predictor_m;
	for (i = 0; i < mb->nf; i++) {
		for (k = 0; k < mb->f_nrep[i]; k++) {
			for (j = 0; j < mb->f_n[i]; j++) {
				c[count + j + k * mb->f_N[i]] = mb->f_diag[i];	/* yes; this is correct */
			}
		}
		count += mb->f_Ntotal[i];		       /* yes; this is correct */
	}

	/*
	 * this is an emergency option to prevent singular matrices (and is known to be >= 0) 
	 */
	if (mb->expert_diagonal_emergencey) {
		for (i = mb->predictor_n + mb->predictor_m; i < N; i++)
			c[i] += mb->expert_diagonal_emergencey;
	}

	if (0) {
		for (i = 0; i < N; i++)
			printf("c[%d]=%g\n", i, c[i]);
	}

	/*
	 * mark those we want to compute and compute the b
	 */
	compute = Calloc(N, char);
	b = Calloc(N, double);
	bfunc = Calloc(N, GMRFLib_bfunc_tp *);
	count = 0;
	if (mb->expert_cpo_manual) {
		/*
		 * if set, then only then only `linear.predictor[idx]' is set
		 */
		for (i = 0; i < mb->predictor_n + mb->predictor_m; i++) {
			compute[count] = (char) 0;
			count++;
		}

		for (i = 0; i < mb->expert_n_cpo_idx; i++) {
			compute[mb->expert_cpo_idx[i]] = (char) 1;
			mb->d[mb->expert_cpo_idx[i]] = 0.0;
		}
		mb->ai_par->cpo_manual = 1;
		mb->output->hyperparameters = GMRFLib_FALSE;

		for (i = 0; i < mb->nf; i++) {
			if (mb->f_bfunc2[i]) {
				for (j = 0; j < mb->f_Ntotal[i]; j++) {
					bfunc[count + j] = Calloc(1, GMRFLib_bfunc_tp);
					bfunc[count + j]->bdef = mb->f_bfunc2[i];
					bfunc[count + j]->idx = j;
				}
			}
			for (j = 0; j < mb->f_Ntotal[i]; j++) {
				compute[count] = (char) 0;
				count++;
			}
		}

		for (i = 0; i < mb->nlinear; i++) {
			compute[count] = (char) 0;
			b[count] = mb->linear_precision[i] * mb->linear_mean[i];
			count++;
		}
		assert(count == N);
	} else {
		/*
		 * as before 
		 */
		for (i = 0; i < mb->predictor_n + mb->predictor_m; i++) {
			compute[count++] = (char) mb->predictor_compute;
		}
		for (i = 0; i < mb->nf; i++) {
			if (mb->f_bfunc2[i]) {
				for (j = 0; j < mb->f_Ntotal[i]; j++) {
					bfunc[count + j] = Calloc(1, GMRFLib_bfunc_tp);
					bfunc[count + j]->bdef = mb->f_bfunc2[i];
					bfunc[count + j]->idx = j;
				}
			}
			for (j = 0; j < mb->f_Ntotal[i]; j++) {
				compute[count] = (char) mb->f_compute[i];
				count++;
			}
		}
		for (i = 0; i < mb->nlinear; i++) {
			compute[count] = (char) mb->linear_compute[i];
			b[count] = mb->linear_precision[i] * mb->linear_mean[i];
			count++;
		}
		if (count != N) {
			P(count);
			P(N);
			assert(count == N);
		}
	}

	// VB correct 
	char *vb_nodes = NULL;
	local_count = 0;
	if (mb->ai_par->vb_enable) {
		vb_nodes = Calloc(N, char);
		count = mb->predictor_n + mb->predictor_m;
		for (i = 0; i < mb->nf; i++) {
			GMRFLib_idx_tp *vb = mb->f_vb_correct[i];
			if ((vb->idx[0] == -1L && mb->f_Ntotal[i] <= mb->ai_par->vb_f_enable_limit_mean)) {
				for (j = 0; j < mb->f_Ntotal[i]; j++) {
					vb_nodes[count + j] = (char) 1;
					local_count++;
				}
			} else if (vb->idx[0] == -1L) {
				int len, jj;
				len = IMAX(1, mb->f_Ntotal[i] / mb->ai_par->vb_f_enable_limit_mean);	/* integer division */
				k = IMAX(1, len / 2);	       /* integer division */
				for (j = 0; j < mb->ai_par->vb_f_enable_limit_mean; j++) {
					jj = (j * len + k) % mb->f_Ntotal[i];
					vb_nodes[count + jj] = (char) 1;
					local_count++;
				}
			} else if (vb->idx[0] >= 0) {
				for (j = 0; j < vb->n; j++) {
					vb_nodes[count + IMIN(vb->idx[j], mb->f_Ntotal[i] - 1)] = (char) 1;
					local_count++;
				}
			}
			count += mb->f_Ntotal[i];
		}
		for (i = 0; i < mb->nlinear; i++) {
			vb_nodes[count++] = (char) 1;
			local_count++;
		}
		if (local_count == 0) {			       /* then there is nothting to correct for */
			Free(vb_nodes);
			vb_nodes = NULL;
		}
	}
	mb->ai_par->vb_nodes_mean = vb_nodes;

	// define the adaptive strategy
	GMRFLib_ai_strategy_tp *adapt = NULL;
	if (mb->ai_par->strategy == GMRFLib_AI_STRATEGY_ADAPTIVE) {
		adapt = Calloc(N, GMRFLib_ai_strategy_tp);
		for (i = 0; i < N; i++) {
			adapt[i] = GMRFLib_AI_STRATEGY_GAUSSIAN;
		}
		count = mb->predictor_n + mb->predictor_m;
		for (i = 0; i < mb->nf; i++) {
			if (mb->f_Ntotal[i] <= mb->ai_par->adapt_max) {
				/*
				 * add also random effects with small size
				 */
				for (j = 0; j < mb->f_Ntotal[i]; j++) {
					adapt[count + j] = GMRFLib_AI_STRATEGY_MEANSKEWCORRECTED_GAUSSIAN;
				}
			}
			count += mb->f_Ntotal[i];
		}
		for (i = 0; i < mb->nlinear; i++) {
			adapt[count++] = GMRFLib_AI_STRATEGY_MEANSKEWCORRECTED_GAUSSIAN;
		}
	}
	mb->ai_par->adapt_strategy = adapt;
	mb->ai_par->adapt_len = (adapt ? N : 0);

	if (G.reorder < 0) {
		size_t nnz = 0;
		int use_g = 0;
		GMRFLib_optimize_reorder(mb->hgmrfm->graph, &nnz, &use_g, &(mb->gn));
		if (GMRFLib_smtp != GMRFLib_SMTP_PARDISO) {
			// ....
		} else {
			GMRFLib_reorder = GMRFLib_REORDER_PARDISO;
		}
		if (GMRFLib_smtp != GMRFLib_SMTP_PARDISO) {
			if (mb->verbose) {
				printf("\tFound optimal reordering=[%s] nnz(L)=[%zu] and use_global_nodes(user)=[%s]\n",
				       GMRFLib_reorder_name(GMRFLib_reorder), nnz, (use_g ? "yes" : "no"));
			}
		}
	}
	if (mb->verbose) {
		if (mb->ntheta) {
			printf("\tList of hyperparameters: \n");
			for (i = 0; i < mb->ntheta; i++) {
				printf("\t\ttheta[%1d] = [%s]\n", i, mb->theta_tag[i]);
			}
		} else {
			printf("\tNone hyperparameters\n");
		}
		printf("\n");
	}
	GMRFLib_ai_store_tp *ai_store = Calloc(1, GMRFLib_ai_store_tp);

	if (mb->output->dic) {
		mb->dic = Calloc(1, GMRFLib_ai_dic_tp);
	} else {
		mb->dic = NULL;
	}
	/*
	 * compute a 'reasonable' initial value for \eta, unless its there from before.
	 */

	int mm = mb->predictor_n + mb->predictor_m;
	Free(G_norm_const_compute);
	Free(G_norm_const);
	for (i = 0; i < G_norm_const_len; i++) {
		Free(G_norm_const_v[i]);
	}
	Free(G_norm_const_v);
	G_norm_const_len = mm;
	G_norm_const_compute = Calloc(mm, char);
	G_norm_const = Calloc(mm, double);
	G_norm_const_v = Calloc(mm, void *);
	for (i = 0; i < mm; i++) {
		G_norm_const[i] = NAN;
		G_norm_const_compute[i] = 1;
	}

	x = Calloc(N, double);
	if (mb->reuse_mode && mb->x_file) {
		if (N != mb->nx_file) {
			char *msg;
			GMRFLib_sprintf(&msg, "N = %1d but nx_file = %1d. Stop.", N, mb->nx_file);
			inla_error_general(msg);
		}
		Memcpy(x, mb->x_file, N * sizeof(double));
		/*
		 * subtract the offset 
		 */
		for (i = 0; i < mb->predictor_ndata; i++) {
			x[i] -= OFFSET3(i);
		}

	} else {
#pragma omp parallel for private(i) num_threads(GMRFLib_openmp->max_threads_outer)
		for (i = 0; i < mb->predictor_ndata; i++) {
			if (mb->d[i]) {
				x[i] = inla_compute_initial_value(i, mb->loglikelihood[i], x, (void *) mb->loglikelihood_arg[i]);
			} else {
				x[i] = 0.0;
			}
			// printf("initial value x[%1d] = %g\n", i, x[i]);
		}
	}

	/*
	 * set the flag to compute correlation-matrix or not 
	 */
	mb->misc_output = Calloc(1, GMRFLib_ai_misc_output_tp);
	if (mb->lc_derived_correlation_matrix) {
		mb->misc_output->compute_corr_lin = mb->nlc;   /* yes, pass the dimension */
	} else {
		mb->misc_output->compute_corr_lin = 0;
	}
	if (mb->output->config) {
		mb->misc_output->configs = Calloc(GMRFLib_MAX_THREADS(), GMRFLib_store_configs_tp *);
	} else {
		mb->misc_output->configs = NULL;
	}

	if (mb->fixed_mode) {
		/*
		 * then there is a request to treat the theta's as fixed and known. This little hack do the job nicely. 
		 */
		mb->ntheta = 0;
		mb->data_ntheta_all = 0;
		mb->theta = NULL;
	}

	mb->transform_funcs = Calloc(N, GMRFLib_transform_array_func_tp *);
	for (i = 0; i < mb->predictor_m + mb->predictor_n; i++) {
		/*
		 * only where we have data (ie n or m), can be a invlinkfunc different from the identity. also the compute-flag must be ON.
		 */
		if (!compute[i]) {
			mb->transform_funcs[i] = NULL;
		} else if (i < mb->predictor_ndata && mb->predictor_invlinkfunc[i]) {
			mb->transform_funcs[i] = Calloc(1, GMRFLib_transform_array_func_tp);
			mb->transform_funcs[i]->func = (GMRFLib_transform_func_tp *) mb->predictor_invlinkfunc[i];
			mb->transform_funcs[i]->arg = mb->predictor_invlinkfunc_arg[i];

			double *cov = NULL;
			if (mb->predictor_invlinkfunc_covariates && mb->predictor_invlinkfunc_covariates[i]) {
				int ncov = mb->predictor_invlinkfunc_covariates[i]->ncol;
				cov = Calloc(ncov, double);
				GMRFLib_matrix_get_row(cov, i, mb->predictor_invlinkfunc_covariates[i]);
			}
			mb->transform_funcs[i]->cov = cov;     /* yes, we store a copy here */
		} else {
			mb->transform_funcs[i] = Calloc(1, GMRFLib_transform_array_func_tp);
			mb->transform_funcs[i]->func = (GMRFLib_transform_func_tp *) link_identity;
			mb->transform_funcs[i]->arg = NULL;
			mb->transform_funcs[i]->cov = NULL;
		}
	}

	if (mb->ai_par->int_design) {
		// make sure the dimensions are right
		if (mb->ntheta != mb->ai_par->int_design->nfactors) {
			char *msg;
			GMRFLib_sprintf(&msg, "ntheta = %1d but int.design says %1d\n", mb->ntheta, mb->ai_par->int_design->nfactors);
			inla_error_general(msg);
		}
	}

	/*
	 * If Gaussian data, then force the strategy to be Gaussian  
	 */
	if (mb->gaussian_data) {
		mb->ai_par->strategy = GMRFLib_AI_STRATEGY_GAUSSIAN;
		mb->ai_par->gaussian_data = mb->gaussian_data;
	}

	/*
	 * Finally, let us do the job...
	 */
	GMRFLib_ai_INLA(&(mb->density),
			// DISABLE THIS FEATURE NOW, IT DOES NOT WORK WELL ENOGUH
			(0 ? &(mb->density_transform) : NULL), (0 ? mb->transform_funcs : NULL),
			// ....
			NULL,
			(mb->output->cpo || mb->expert_cpo_manual ? &(mb->cpo) : NULL),
			(mb->output->po ? &(mb->po) : NULL),
			mb->dic,
			(mb->output->mlik ? &(mb->mlik) : NULL),
			compute, mb->theta, mb->ntheta,
			extra, (void *) mb,
			x, b, c, NULL, bfunc, mb->d,
			loglikelihood_inla, (void *) mb,
			mb->hgmrfm->graph, mb->hgmrfm->Qfunc, mb->hgmrfm->Qfunc_arg, mb->hgmrfm->constr, mb->ai_par, ai_store,
			mb->nlc, mb->lc_lc, &(mb->density_lin), mb->misc_output, NULL, rpreopt);

	/*
	 * add the offsets to the linear predictor. Add the offsets to the 'configs' (if any), at a later stage. 
	 */
#pragma omp parallel for private(i) num_threads(GMRFLib_openmp->max_threads_outer)
	for (i = 0; i < mb->predictor_n + mb->predictor_m; i++) {
		GMRFLib_density_tp *d;
		if (mb->density[i]) {
			d = mb->density[i];
			GMRFLib_density_new_mean(&(mb->density[i]), d, d->std_mean + OFFSET3(i));
			GMRFLib_free_density(d);
		}
	}

	/*
	 * add the offset to 'x' 
	 */
	for (i = 0; i < mb->predictor_n + mb->predictor_m; i++) {
		x[i] += OFFSET3(i);
	}

	/*
	 * need to correct configs for the integration weight, as it should not be there for _orig. only ccd has integration weights 
	 */

	Free(mb->x_file);				       /* yes, and then */
	mb->x_file = x;					       /* just take over */
	mb->nx_file = N;

	GMRFLib_free_ai_store(ai_store);
	Free(b);
	Free(c);
	Free(compute);

	return INLA_OK;
}
