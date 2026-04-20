		// compute the correction to the log-posterior for the hyperpar, based on the new mode
		double ldens_correction = 0.0;
		if (ai_par->vb_hyperpar_correct) {
			double *new_mode = Calloc(graph->n, double);
			for(i = 0; i < graph->n; i++) {
				new_mode[i] = mode[i] + gsl_vector_get(delta_mu, i);
			}
			
			for(ii = 0; ii < d_idx->n; ii++){
				double ld[2], xx[2];
				i = d_idx->idx[ii];
				xx[0] = new_mode[i];
				xx[1] = ai_store->mode[i];
				loglFunc(&ld[0], &xx[0], 1, i, mode, NULL, loglFunc_arg);
				loglFunc(&ld[1], &xx[1], 1, i, ai_store->mode, NULL, loglFunc_arg);
				ldens_correction += d[i] * (ld[0] - ld[1]);
			}
			
			// prior
			double res_new, res_old;
			GMRFLib_xQx2(&res_new, new_mode, graph, Qfunc, Qfunc_arg, c);
			GMRFLib_xQx2(&res_old, ai_store->mode, graph, Qfunc, Qfunc_arg, c);
			ldens_correction += res_new - res_old;

			// linear term and the contributions from bfunc()
			double *bnew = NULL, con = 0.0;
			GMRFLib_bnew(&bnew, &con, graph->n, b, bfunc);
			for(i = 0; i < graph->n; i++){
				ldens_correction += bnew[i] * (new_mode[i] - ai_store->mode[i]);
			}
			
			// change in the conditional distribution
			memcpy(ai_store->problem->sample, ai_store->mode, graph->n*sizeof(double));
			GMRFLib_evaluate(ai_store->problem);
			ldens_correction += ai_store->problem->sub_logdens;
			memcpy(ai_store->problem->sample, new_mode, graph->n*sizeof(double));
			GMRFLib_evaluate(ai_store->problem);
			ldens_correction -= ai_store->problem->sub_logdens;
			*ldens_hyperpar_corr = ldens_correction;
			
			Free(new_mode);
			Free(bnew);
		} else {
			if (ldens_hyperpar_corr) {
				*ldens_hyperpar_corr = 0.0;
			}
		}
