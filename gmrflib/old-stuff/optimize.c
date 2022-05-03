
int GMRFLib_optimize_set_store_flags(GMRFLib_store_tp * store)
{
	/*
	 * set flags for 'store' 
	 */
	if (store) {
		/*
		 * create logicals 
		 */
		store_store_sub_graph = (store->sub_graph ? 0 : 1);
		store_use_sub_graph = !store_store_sub_graph;
		store_store_remap = (store->remap ? 0 : 1);
		store_use_remap = !store_store_remap;

		if (GMRFLib_valid_smtp((int) store->smtp) == GMRFLib_TRUE) {
			store_smtp = store->smtp;
		} else {
			store_smtp = GMRFLib_smtp;
		}
		if (store_smtp == GMRFLib_SMTP_TAUCS) {
			store_store_symb_fact = (store->TAUCS_symb_fact ? 0 : 1);
			store_use_symb_fact = !store_store_symb_fact;
		} else {
			store_store_symb_fact = 0;
			store_use_symb_fact = 0;
		}
	} else {
		store_store_sub_graph = 0;
		store_use_sub_graph = 0;
		store_store_remap = 0;
		store_use_remap = 0;
		store_store_symb_fact = 0;
		store_use_symb_fact = 0;
		store_smtp = GMRFLib_smtp;
	}
	return GMRFLib_SUCCESS;
}

int GMRFLib_optimize(double *mode, double *b, double *c, double *mean,
		     GMRFLib_graph_tp * graph, GMRFLib_Qfunc_tp * Qfunc, void *Qfunc_args,
		     GMRFLib_constr_tp * constr, double *d, GMRFLib_logl_tp * loglFunc, void *loglFunc_arg, GMRFLib_optimize_param_tp * optpar)
{
	GMRFLib_ENTER_ROUTINE;
	GMRFLib_EWRAP1(GMRFLib_optimize_store(mode, b, c, mean, graph, Qfunc, Qfunc_args, constr, d, loglFunc, loglFunc_arg, optpar, NULL));
	GMRFLib_LEAVE_ROUTINE;
	return GMRFLib_SUCCESS;
}

int GMRFLib_optimize_store(double *mode, double *b, double *c, double *mean,
			   GMRFLib_graph_tp * graph, GMRFLib_Qfunc_tp * Qfunc, void *Qfunc_args,
			   GMRFLib_constr_tp * constr,
			   double *d, GMRFLib_logl_tp * loglFunc, void *loglFunc_arg, GMRFLib_optimize_param_tp * optpar, GMRFLib_store_tp * store)
{
	/*
	 * locate the mode starting in 'mode' using a conjugate gradient algorithm where the search direction is Q-ortogonal to 
	 * the m previous ones.
	 * 
	 * exp(-0.5(x-mean)'Q(x-mean)+b'x + \sum_i d_i loglFunc(x_i,i,logl_arg)) && fixed-flags
	 * 
	 * and linear deterministic and stochastic constaints 
	 */
	int sub_n, i, free_optpar;
	double *cc = NULL, *initial_value = NULL;
	GMRFLib_store_tp *store_ptr;
	GMRFLib_optimize_problem_tp *opt_problem = NULL;

	GMRFLib_ENTER_ROUTINE;

	GMRFLib_optimize_set_store_flags(store);

	/*
	 * our first task, is to convert the opt_problem so we get rid of the fixed_values and then write the function in the
	 * canonical form
	 * 
	 * large parts here is adapted from `problem-setup.c'!!! 
	 */

	/*
	 * create new opt_problem 
	 */
	opt_problem = Calloc(1, GMRFLib_optimize_problem_tp);

	/*
	 * get the options 
	 */
	if (optpar) {
		opt_problem->optpar = optpar;
		free_optpar = 0;
	} else {
		GMRFLib_default_optimize_param(&(opt_problem->optpar));
		free_optpar = 1;
	}

	/*
	 * first, find the new graph. 
	 */
	if (store_use_sub_graph) {
		/*
		 * copy from store 
		 */
		assert(store);
		GMRFLib_EWRAP1(GMRFLib_graph_duplicate(&(opt_problem->sub_graph), store->sub_graph));
	} else {
		/*
		 * compute it 
		 */
		GMRFLib_EWRAP1(GMRFLib_graph_comp_subgraph(&(opt_problem->sub_graph), graph, NULL, NULL));

		/*
		 * store it in store if requested 
		 */
		if (store_store_sub_graph) {
			GMRFLib_EWRAP1(GMRFLib_graph_duplicate(&(store->sub_graph), opt_problem->sub_graph));
		}
	}

	sub_n = opt_problem->sub_graph->n;
	if (sub_n == 0) {				       /* fast return if there is nothing todo */
		GMRFLib_graph_free(opt_problem->sub_graph);
		Free(opt_problem);
		GMRFLib_LEAVE_ROUTINE;
		return GMRFLib_SUCCESS;
	}

	/*
	 * setup space & misc 
	 */
	{
		{
			int ii;
			opt_problem->mode = Calloc(sub_n, double);

			for (ii = 0; ii < sub_n; ii++) {
				opt_problem->mode[ii] = mode[ii];
			}
		}
		{
			opt_problem->b = Calloc(sub_n, double);

			if (b) {
				int ii;

				for (ii = 0; ii < sub_n; ii++) {
					opt_problem->b[ii] = b[ii];
				}
			}
		}
		{
			opt_problem->d = Calloc(sub_n, double);

			if (d) {
				int ii;

				for (ii = 0; ii < sub_n; ii++) {
					opt_problem->d[ii] = d[ii];
				}
			}
		}
		{
			cc = Calloc(sub_n, double);

			if (c) {
				int ii;

				for (ii = 0; ii < sub_n; ii++) {
					cc[ii] = c[ii];
				}
			}
		}
		{
			opt_problem->x_vec = Calloc(graph->n, double);
			Memcpy(opt_problem->x_vec, mode, graph->n * sizeof(double));
		}
	}

	/*
	 * FIXME: i might want to make a wrapper of this one, then REMEMBER TO CHANGE ->map[i] in all calls to
	 * GMRFLib_2order_approx to i, or what it should be.!!! 
	 */
	opt_problem->loglFunc = loglFunc;		       /* i might want to make a wrapper of this one! */
	opt_problem->loglFunc_arg = loglFunc_arg;

	/*
	 * make the arguments to the wrapper function 
	 */
	opt_problem->sub_Qfunc = GMRFLib_Qfunc_wrapper;
	opt_problem->sub_Qfunc_arg = Calloc(1, GMRFLib_Qfunc_arg_tp);
	opt_problem->sub_Qfunc_arg->diagonal_adds = cc;

	opt_problem->sub_Qfunc_arg->user_Qfunc = Qfunc;
	opt_problem->sub_Qfunc_arg->user_Qfunc_args = Qfunc_args;

	/*
	 * now compute the new 'effective' b, and then the mean. recall to add the 'c' term manually, since we're using the
	 * original Qfunc. 
	 */
	if (mean) {
		double *tmp = NULL;
		tmp = Calloc(graph->n, double);

		GMRFLib_Qx(tmp, mean, graph, Qfunc, Qfunc_args);
		for (i = 0; i < graph->n; i++) {
			opt_problem->b[i] += tmp[i] + cc[i] * mean[i];
		}
		Free(tmp);
	}

	/*
	 * save initial value and return that unchanged if fail to converge 
	 */
	initial_value = Calloc(sub_n, double);
	Memcpy(initial_value, opt_problem->mode, sub_n * sizeof(double));

	GMRFLib_ASSERT(opt_problem->optpar->opt_type == GMRFLib_OPTTYPE_CG ||
		       opt_problem->optpar->opt_type == GMRFLib_OPTTYPE_NR ||
		       opt_problem->optpar->opt_type == GMRFLib_OPTTYPE_SAFECG
		       || opt_problem->optpar->opt_type == GMRFLib_OPTTYPE_SAFENR, GMRFLib_EPARAMETER);

	/*
	 * do stuff relating to the constraint 
	 */
	if (constr) {
		/*
		 * constraints only allowed when using the newton-raphson optimizer 
		 */
		double *b_add = NULL;

		GMRFLib_ASSERT(constr && (opt_problem->optpar->opt_type == GMRFLib_OPTTYPE_NR ||
					  opt_problem->optpar->opt_type == GMRFLib_OPTTYPE_SAFENR), GMRFLib_EPARAMETER);

		b_add = Calloc(sub_n, double);

		GMRFLib_EWRAP1(GMRFLib_recomp_constr(&(opt_problem->sub_constr), constr, mode, b_add, NULL, graph, opt_problem->sub_graph));

		for (i = 0; i < sub_n; i++) {
			opt_problem->b[i] += b_add[i];
		}
		Free(b_add);
	}

	if (opt_problem->optpar->opt_type == GMRFLib_OPTTYPE_SAFENR || opt_problem->optpar->opt_type == GMRFLib_OPTTYPE_SAFECG) {
		/*
		 * to be more safe regarding initial valus we proceed as follows: first we change the graph to a diagonal, so
		 * to optimize without any interactions using NR/CG. this optimum is then the startvalue for the NR/CG
		 * optimizer. 
		 */
		int nzero_idx = -1;
		GMRFLib_graph_tp *diag_graph = NULL, *tmp_graph;

		if (store && !(store->diag_store)) {
			/*
			 * use diag_store here 
			 */
			store->diag_store = Calloc(1, GMRFLib_store_tp);
		}
		store_ptr = (store ? store->diag_store : store);

		GMRFLib_EWRAP1(GMRFLib_graph_comp_subgraph(&diag_graph, opt_problem->sub_graph, NULL, NULL));
		for (i = 0; i < diag_graph->n; i++)
			if (diag_graph->nnbs[i]) {
				diag_graph->nnbs[i] = 0;
				if (nzero_idx < 0) {
					nzero_idx = i;	       /* the first non-zero pointer */
				}
			}

		tmp_graph = opt_problem->sub_graph;
		opt_problem->sub_graph = diag_graph;

		if (opt_problem->optpar->opt_type == GMRFLib_OPTTYPE_SAFECG) {
			GMRFLib_EWRAP1(GMRFLib_optimize2(opt_problem, store_ptr));
		} else {
			GMRFLib_EWRAP1(GMRFLib_optimize3(opt_problem, store_ptr));
		}
		opt_problem->sub_graph = tmp_graph;

		if (nzero_idx >= 0) {
			diag_graph->nnbs[nzero_idx] = 1;       /* to make _free_graph work correct */
		}
		GMRFLib_graph_free(diag_graph);
	}

	store_ptr = store;

	switch (opt_problem->optpar->opt_type) {
	case GMRFLib_OPTTYPE_CG:
	case GMRFLib_OPTTYPE_SAFECG:
		GMRFLib_EWRAP1(GMRFLib_optimize2(opt_problem, store_ptr));
		if (0) {				       /* FIXME: if (fail) then do as follows */
			Memcpy(opt_problem->mode, initial_value, sub_n * sizeof(double));
			Free(initial_value);
			GMRFLib_ERROR(GMRFLib_EOPTCG);
		}
		break;

	case GMRFLib_OPTTYPE_NR:
	case GMRFLib_OPTTYPE_SAFENR:
		GMRFLib_EWRAP1(GMRFLib_optimize3(opt_problem, store_ptr));
		if (0) {				       /* FIXME: if (fail) then do as follows */
			Memcpy(opt_problem->mode, initial_value, sub_n * sizeof(double));
			Free(initial_value);
			GMRFLib_ERROR(GMRFLib_EOPTNR);
		}
		break;
	}
	/*
	 * copy back 
	 */
	for (i = 0; i < sub_n; i++) {
		mode[i] = opt_problem->mode[i];
	}

	/*
	 * free what's malloced and return 
	 */
	Free(opt_problem->mode);
	Free(opt_problem->b);
	Free(opt_problem->d);
	Free(opt_problem->sub_Qfunc_arg);
	Free(opt_problem->x_vec);
	if (constr) {
		GMRFLib_free_constr(opt_problem->sub_constr);
	}
	if (free_optpar) {
		Free(opt_problem->optpar);
	}
	GMRFLib_graph_free(opt_problem->sub_graph);
	Free(opt_problem);
	Free(initial_value);
	Free(cc);

	GMRFLib_LEAVE_ROUTINE;
	return GMRFLib_SUCCESS;
}

int GMRFLib_optimize2(GMRFLib_optimize_problem_tp * opt_problem, GMRFLib_store_tp * UNUSED(store))
{

	/*
	 * optimize, using GMRFLib_Qfunc_wrapper as Qfunc. NOTE: the loglFunc is NOT wrapped, so the indices are in the real
	 * world.
	 * 
	 * this is conjugate gradient optimizer 
	 */

	int sub_n, i, nsdir, iter, sdir_indx, n_correct, fail, id;
	double **sdir = NULL, *grad = NULL, *omode = NULL, err, *c_orig = NULL;

	id = GMRFLib_thread_id;
	nsdir = opt_problem->optpar->nsearch_dir + 1;	       /* also contain the current dir */
	sub_n = opt_problem->sub_graph->n;
	grad = Calloc(sub_n, double);
	omode = Calloc(sub_n, double);
	c_orig = Calloc(sub_n, double);

	Memcpy(c_orig, opt_problem->sub_Qfunc_arg->diagonal_adds, sub_n * sizeof(double));

	sdir = Calloc(nsdir, double *);

	for (i = 0; i < nsdir; i++) {
		sdir[i] = Calloc(sub_n, double);
	}

	if (opt_problem->optpar->fp) {
		fprintf(opt_problem->optpar->fp, "\n%6s%22s%6s%22s\n", "Iter", "Value", "SubIt", "StepLength");
		fprintf(opt_problem->optpar->fp, " --------------------------------------------------------\n");
	}

	for (iter = 0, sdir_indx = 0, fail = 1;
	     iter < IMAX(opt_problem->optpar->max_iter, opt_problem->optpar->fixed_iter); sdir_indx = ((sdir_indx + 1) % nsdir), iter++) {
		if (opt_problem->optpar->fp) {
			fprintf(opt_problem->optpar->fp, "%6d", iter);
		}

		Memcpy(opt_problem->sub_Qfunc_arg->diagonal_adds, c_orig, sub_n * sizeof(double));

		/*
		 * first compute the gradient, -Qx + b 
		 */
		GMRFLib_Qx(grad, opt_problem->mode, opt_problem->sub_graph, opt_problem->sub_Qfunc, (void *) (opt_problem->sub_Qfunc_arg));
		for (i = 0; i < sub_n; i++) {
			grad[i] = opt_problem->b[i] - grad[i];
		}

		/*
		 * then compute contribution from the loglFunc 
		 */
		for (i = 0; i < sub_n; i++) {
			opt_problem->x_vec[i] = opt_problem->mode[i];
		}
//#pragma omp parallel for private(i)
		for (i = 0; i < sub_n; i++) {
			GMRFLib_thread_id = id;
			if (opt_problem->d[i]) {
				double bcoof, ccoof;

				GMRFLib_2order_taylor(NULL, &bcoof, &ccoof, NULL, opt_problem->d[i],
						      opt_problem->mode[i], i, opt_problem->x_vec,
						      opt_problem->loglFunc, opt_problem->loglFunc_arg, &(opt_problem->optpar->step_len),
						      &(opt_problem->optpar->stencil));
				grad[i] += bcoof;
				opt_problem->sub_Qfunc_arg->diagonal_adds[i] += DMAX(0.0, -ccoof);
			}
		}
		GMRFLib_thread_id = id;

		/*
		 * subtract terms Q-ortogonal to the previous search-directions 
		 */
		if ((iter + 1) % opt_problem->optpar->restart_interval != 0) {
			n_correct = IMIN(opt_problem->optpar->nsearch_dir, iter);
			for (i = 0; i < n_correct; i++) {
				GMRFLib_Qadjust(grad, sdir[MOD(sdir_indx - i - 1, nsdir)],
						opt_problem->sub_graph, opt_problem->sub_Qfunc, (void *) (opt_problem->sub_Qfunc_arg));
			}
		} else {
			for (i = 0; i < opt_problem->optpar->nsearch_dir; i++) {
				Memset(sdir[i], 0, sub_n * sizeof(double));
			}
		}
		Memcpy(opt_problem->sub_Qfunc_arg->diagonal_adds, c_orig, sub_n * sizeof(double));	/* go back to original */
		Memcpy(sdir[sdir_indx], grad, sub_n * sizeof(double));

		/*
		 * we got a search-direction, do a line-search in that direction and update 'mode' 
		 */
		Memcpy(omode, opt_problem->mode, sub_n * sizeof(double));
		GMRFLib_EWRAP0(GMRFLib_linesearch(opt_problem, sdir[sdir_indx]));

		for (i = 0, err = 0.0; i < sub_n; i++) {
			err += SQR(omode[i] - opt_problem->mode[i]);
		}
		err = sqrt(DMAX(0.0, err) / sub_n);
		if (opt_problem->optpar->fp) {
			fprintf(opt_problem->optpar->fp, "%22.10e\n", err);
		}

		if (opt_problem->optpar->fixed_iter > 0) {
			if (opt_problem->optpar->fixed_iter == iter + 1) {
				fail = 0;
				break;
			}
		} else {
			if (err < opt_problem->optpar->abserr_func) {
				fail = 0;
				break;
			}
		}
	}

	Memcpy(opt_problem->sub_Qfunc_arg->diagonal_adds, c_orig, sub_n * sizeof(double));

	for (i = 0; i < nsdir; i++) {
		Free(sdir[i]);
	}
	Free(sdir);
	Free(omode);
	Free(c_orig);
	Free(grad);

	return (fail ? GMRFLib_EOPTCG : GMRFLib_SUCCESS);
}

int GMRFLib_Qadjust(double *dir, double *odir, GMRFLib_graph_tp * graph, GMRFLib_Qfunc_tp * Qfunc, void *Qfunc_arg)
{
	/*
	 * compute: dir := dir - (dir,odir)/(odir,odir), where (u,v) = u^TQv 
	 */

	int i;
	double a, b, c, *v = NULL;

	v = Calloc(graph->n, double);

	GMRFLib_Qx(v, odir, graph, Qfunc, Qfunc_arg);
	for (i = 0, a = b = 0.0; i < graph->n; i++) {
		a += v[i] * dir[i];
		b += v[i] * odir[i];
	}

	if (b > 0.0 && !ISZERO(a)) {
		c = a / b;				       /* FIX? */
		for (i = 0; i < graph->n; i++) {
			dir[i] -= c * odir[i];
		}
	}

	Free(v);
	return GMRFLib_SUCCESS;
}

double GMRFLib_linesearch_func(double length, double *dir, GMRFLib_optimize_problem_tp * opt_problem)
{
	int i, sub_n, id;
	double *v = NULL, *u = NULL, fval = 0.0;

	assert(dir);
	id = GMRFLib_thread_id;
	sub_n = opt_problem->sub_graph->n;
	v = Calloc(sub_n, double);
	u = Calloc(sub_n, double);

	if (length != 0.0) {
		for (i = 0; i < sub_n; i++) {
			u[i] = opt_problem->mode[i] + length * dir[i];
		}
	} else {
		Memcpy(u, opt_problem->mode, sub_n * sizeof(double));
	}
	GMRFLib_Qx(v, u, opt_problem->sub_graph, opt_problem->sub_Qfunc, (void *) (opt_problem->sub_Qfunc_arg));

	for (i = 0; i < sub_n; i++) {
		opt_problem->x_vec[i] = u[i];
	}
#if defined(_OPENMP)
	{
		double sum = 0.0;

//#pragma omp parallel for private(i) reduction(+: sum)
		for (i = 0; i < sub_n; i++) {
			GMRFLib_thread_id = id;

			double logll;
			sum += (-0.5 * v[i] + opt_problem->b[i]) * u[i];
			if (opt_problem->d[i]) {
				(*(opt_problem->loglFunc)) (&logll, &u[i], 1, i, opt_problem->x_vec, NULL, opt_problem->loglFunc_arg);
				sum += opt_problem->d[i] * logll;
			}
		}
		GMRFLib_thread_id = id;
		fval += sum;
	}
#else
	for (i = 0, fval = 0.0; i < sub_n; i++) {
		double logll;

		GMRFLib_thread_id = id;
		fval += (-0.5 * v[i] + opt_problem->b[i]) * u[i];
		if (opt_problem->d[i]) {
			(*(opt_problem->loglFunc)) (&logll, &u[i], 1, i, opt_problem->x_vec, NULL, opt_problem->loglFunc_arg);
			fval += opt_problem->d[i] * logll;
		}
	}
	GMRFLib_thread_id = id;
#endif

	Free(u);
	Free(v);

	return fval;
}

int GMRFLib_linesearch(GMRFLib_optimize_problem_tp * opt_problem, double *dir)
{
	int i, sub_n, iter = 0, fail, id;
	double len, deriv, dderiv, *u = NULL, err, *loglikgrad = NULL, *loglikggrad = NULL;
	double *len_history, len_eps = sqrt(FLT_EPSILON), periodic_flag = 0;

	id = GMRFLib_thread_id;
	sub_n = opt_problem->sub_graph->n;
	len_history = Calloc(opt_problem->optpar->max_linesearch_iter, double);
	loglikggrad = Calloc(sub_n, double);
	loglikgrad = Calloc(sub_n, double);
	u = Calloc(sub_n, double);

	for (iter = 0, fail = 1; iter < opt_problem->optpar->max_linesearch_iter; iter++) {
		for (i = 0; i < sub_n; i++) {
			opt_problem->x_vec[i] = opt_problem->mode[i];
		}
//#pragma omp parallel for private(i)
		for (i = 0; i < sub_n; i++) {
			GMRFLib_thread_id = id;
			if (opt_problem->d[i]) {
				GMRFLib_2order_taylor(NULL, &loglikgrad[i], &loglikggrad[i], NULL, 1.0,
						      opt_problem->mode[i], i, opt_problem->x_vec,
						      opt_problem->loglFunc, opt_problem->loglFunc_arg, &(opt_problem->optpar->step_len),
						      &(opt_problem->optpar->stencil));
			} else {
				loglikgrad[i] = loglikggrad[i] = 0.0;
			}
		}
		GMRFLib_thread_id = id;

		GMRFLib_Qx(u, dir, opt_problem->sub_graph, opt_problem->sub_Qfunc, (void *) opt_problem->sub_Qfunc_arg);

		deriv = dderiv = 0.0;
		for (i = 0; i < sub_n; i++) {
			deriv += -opt_problem->mode[i] * u[i] + opt_problem->b[i] * dir[i] + opt_problem->d[i] * dir[i] * loglikgrad[i];
			dderiv += -dir[i] * u[i] + opt_problem->d[i] * SQR(dir[i]) * loglikggrad[i];
		}

		if (!ISZERO(deriv) && !ISZERO(dderiv)) {
			double len_max = 1.0;

			len = len_history[iter] = DMIN(len_max, DMAX(-len_max, -deriv / dderiv));

			/*
			 * make a robustnes check for the behaviour len = a, -a, a, -a, a, -a, etc... 
			 */
			if (iter >= 4) {
				if (ABS(len_history[iter] - len_history[iter - 2]) < len_eps &&
				    ABS(len_history[iter - 1] - len_history[iter - 3]) < len_eps
				    && ABS(len_history[iter] + len_history[iter - 1]) < len_eps) {
					len = len_history[iter] / 2.;
					len_history[iter] = len;
					periodic_flag = 1;     /* !!!stop after this step!!! */
					if (0)
						printf("\nperiodic behaviour detected, set len = %f\n", len);
				}
			}

			for (i = 0, err = 0.0; i < sub_n; i++) {
				opt_problem->mode[i] += len * dir[i];
				err += SQR(dir[i]);
			}
		} else {
			len = err = 0.0;
		}
		err = sqrt(DMAX(err * SQR(len) / sub_n, 0.0));
		if (err < opt_problem->optpar->abserr_step) {
			fail = 0;
			break;
		}
		if (periodic_flag) {
			fail = 0;
			break;
		}
		if (0) {
			printf("\n%d %.12f %.12f\n", iter, len, GMRFLib_linesearch_func(0.0, dir, opt_problem));
		}
	}
	if (fail)
		GMRFLib_ERROR(GMRFLib_EOPTCGLINE);

	if (opt_problem->optpar->fp) {
		fprintf(opt_problem->optpar->fp, "%22.10e%6d", GMRFLib_linesearch_func(0.0, dir, opt_problem), iter + 1);
	}

	Free(loglikgrad);
	Free(loglikggrad);
	Free(u);
	Free(len_history);
	return GMRFLib_SUCCESS;
}

int GMRFLib_optimize3(GMRFLib_optimize_problem_tp * opt_problem, GMRFLib_store_tp * store)
{
	/*
	 * optimize, using GMRFLib_Qfunc_wrapper as Qfunc. NOTE: the loglFunc is NOT wrapped, so the indices are in the real
	 * world.
	 * 
	 * this is the newton-raphson optimizer 
	 */

	int sub_n, i, iter, fail = 1, id;
	int *idxs = NULL, nidx = 0;
	double *bb = NULL, err, *c_orig = NULL, step_factor, f;
	GMRFLib_problem_tp *problem = NULL;

	id = GMRFLib_thread_id;
	sub_n = opt_problem->sub_graph->n;
	bb = Calloc(sub_n, double);
	c_orig = Calloc(sub_n, double);

	if (opt_problem->optpar->nr_step_factor > 0) {
		step_factor = opt_problem->optpar->nr_step_factor;
	} else {
		/*
		 * default choice
		 */
		step_factor = 1.0;
	}

	Memcpy(c_orig, opt_problem->sub_Qfunc_arg->diagonal_adds, sub_n * sizeof(double));

	if (opt_problem->optpar->fp) {
		fprintf(opt_problem->optpar->fp, "\n%6s%22s%22s\n", "Iter", "Value", "StepLength");
		fprintf(opt_problem->optpar->fp, " --------------------------------------------------\n");
	}

	for (iter = 0, fail = 1; iter < IMAX(opt_problem->optpar->max_iter, opt_problem->optpar->fixed_iter); iter++) {
		if (opt_problem->optpar->fp) {
			fprintf(opt_problem->optpar->fp, "%6d", iter);
		}
		Memcpy(opt_problem->sub_Qfunc_arg->diagonal_adds, c_orig, sub_n * sizeof(double));

		/*
		 * compute contribution from the loglFunc 
		 */
		for (i = 0; i < sub_n; i++) {
			opt_problem->x_vec[i] = opt_problem->mode[i];
		}

		if (!idxs) {
			/*
			 * be somewhat more clever in this loop to divide the work better 
			 */
			idxs = Calloc(sub_n, int);

			for (i = 0; i < sub_n; i++) {
				if (opt_problem->d[i]) {
					idxs[nidx++] = i;
				} else {
					bb[i] = opt_problem->b[i];
				}
			}
		}
//#pragma omp parallel for private(i)
		for (i = 0; i < nidx; i++) {
			GMRFLib_thread_id = id;

			int idx;
			double bcoof, ccoof;
			double cmin = 0.0;

			idx = idxs[i];
			GMRFLib_2order_approx(NULL, &bcoof, &ccoof, NULL, opt_problem->d[idx],
					      opt_problem->mode[idx], idx, opt_problem->x_vec,
					      opt_problem->loglFunc, opt_problem->loglFunc_arg, &(opt_problem->optpar->step_len),
					      &(opt_problem->optpar->stencil), &cmin);
			bb[idx] = opt_problem->b[idx] + bcoof;
			opt_problem->sub_Qfunc_arg->diagonal_adds[idx] += ccoof;
		}
		GMRFLib_thread_id = id;

		GMRFLib_EWRAP0(GMRFLib_init_problem_store(&problem, opt_problem->mode, bb, NULL, NULL, opt_problem->sub_graph,
							  opt_problem->sub_Qfunc, (void *) opt_problem->sub_Qfunc_arg,
							  opt_problem->sub_constr, store));
		f = DMIN(1.0, step_factor * (iter + 1));
		if (problem->sub_constr) {
			for (i = 0, err = 0.0; i < sub_n; i++) {
				err += SQR(problem->sub_mean_constr[i] - opt_problem->mode[i]);
				opt_problem->mode[i] = opt_problem->mode[i] + f * (problem->sub_mean_constr[i] - opt_problem->mode[i]);
			}
		} else {
			for (i = 0, err = 0.0; i < sub_n; i++) {
				err += SQR(problem->sub_mean[i] - opt_problem->mode[i]);
				opt_problem->mode[i] = opt_problem->mode[i] + f * (problem->sub_mean[i] - opt_problem->mode[i]);
			}
		}

		err = sqrt(DMAX(0.0, err) / sub_n);
		if (opt_problem->optpar->fp) {
			Memcpy(opt_problem->sub_Qfunc_arg->diagonal_adds, c_orig, sub_n * sizeof(double));
			fprintf(opt_problem->optpar->fp, "%22.10e%22.10e\n", GMRFLib_linesearch_func(0.0, NULL, opt_problem), err);
		}

		if (opt_problem->optpar->fixed_iter > 0) {
			if (opt_problem->optpar->fixed_iter == iter + 1) {
				fail = 0;
				break;
			}
		} else {
			if (err < opt_problem->optpar->abserr_step) {
				fail = 0;
				break;
			}
		}
	}

	Memcpy(opt_problem->sub_Qfunc_arg->diagonal_adds, c_orig, sub_n * sizeof(double));

	GMRFLib_free_problem(problem);
	Free(bb);
	Free(c_orig);
	Free(idxs);

	return (fail ? GMRFLib_EOPTNR : GMRFLib_SUCCESS);
}

static int store_store_sub_graph = 0, store_use_sub_graph = 0, store_store_remap = 0, store_use_remap = 0, store_store_symb_fact =
    0, store_use_symb_fact = 0, store_smtp = 0;

#pragma omp threadprivate (store_store_sub_graph, store_use_sub_graph, store_store_remap, store_use_remap, store_store_symb_fact, store_use_symb_fact, store_smtp)

