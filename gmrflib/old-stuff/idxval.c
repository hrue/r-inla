
int GMRFLib_idxval_nsort_x_OLD(GMRFLib_idxval_tp ** hold, int n, int nt, int build_groups, int merge_groups)
{
	const int limit = 8;
	if (merge_groups) {
		build_groups = 1;
	}

	int debug_details = 0;
	// int debug = GMRFLib_DEBUG_IF_TRUE();
	int debug = 0;

	int nmax = 1;
	for (int i = 0; i < n; i++) {
		GMRFLib_idxval_tp *h = hold[i];
		nmax = IMAX(nmax, h->n);
	}
	if (nmax == 1) {
		return GMRFLib_SUCCESS;
	}
#define CODE_BLOCK							\
	for (int i = 0; i < n; i++) {					\
		GMRFLib_idxval_tp *h = hold[i];				\
		if (!h) continue;					\
		double tref[6] =  {0, 0, 0, 0, 0, 0};			\
		if (debug) tref[0] -= GMRFLib_cpu();			\
		if (h->n > 1) {						\
			int is_sorted = 1;				\
			for(int j = 1; is_sorted && j < h->n; j++) {	\
				is_sorted = (h->idx[j] >= h->idx[j-1]);	\
			}						\
			if (!is_sorted) {				\
				my_sort2_id(h->idx, h->val, h->n);	\
			}						\
		}							\
		if (debug) tref[0] += GMRFLib_cpu();			\
		if (debug) tref[1] -= GMRFLib_cpu();			\
									\
		/* make them unique */					\
		if (h->n > 1) {						\
			/* do a dry-run, as its cheaper */		\
			int all_unique = 1;				\
			for(int j = 1; all_unique && j < h->n; j++) {	\
				all_unique = (h->idx[j] > h->idx[j - 1]); \
			}						\
			if (!all_unique) {				\
				int k = 0;				\
				for(int j = 1; j < h->n; j++) {		\
					if (h->idx[j] != h->idx[k]) {	\
						k++;			\
						h->idx[k] = h->idx[j];	\
						h->val[k] = h->val[j];	\
					} else {			\
						h->val[k] += h->val[j];	\
					}				\
				}					\
				if (debug && (h->n > k+1)) {		\
					printf("Make unique: reduce length from %d to %d\n", h->n, k+1); \
					fflush(stdout);			\
				}					\
				h->n = k+1;				\
			}						\
		}							\
		if (debug) tref[1] += GMRFLib_cpu();			\
		if (debug) tref[2] -= GMRFLib_cpu();			\
									\
		if (h->n <= limit || !build_groups) {			\
			h->preference = IDXVAL_SERIAL_MKL;		\
			continue;					\
		}							\
									\
		/* build basic groups with one group for each sequence and then one for each individual */ \
		int ng = 1;						\
		for (int j = 1; j < h->n; j++) {			\
			if (h->idx[j] != h->idx[j - 1] + 1) {		\
				ng++;					\
			}						\
		}							\
		int *g_istart = Calloc(ng + 1, int);			\
		int *g_len = Calloc(ng + 1, int);			\
									\
		int k = 0;						\
		g_istart[0] = 0;					\
		for (int j = 1; j < h->n; j++) {			\
			if (h->idx[j] != h->idx[j - 1] + 1) {		\
				g_len[k] = j - g_istart[k];		\
				k++;					\
				g_istart[k] = j;			\
			}						\
		}							\
		g_len[ng - 1] = h->n - g_istart[ng - 1];		\
		if (debug_details && debug) {				\
			GMRFLib_idxval_printf(stdout, h, "OLD");	\
			fflush(stdout);					\
		}							\
									\
		k = h->n;				       /* so the code below works */ \
		int kg = ng;				       /* so the code below works */ \
									\
		if (debug) tref[2] += GMRFLib_cpu();			\
		if (debug) tref[3] -= GMRFLib_cpu();			\
									\
		h->n = k;						\
		h->g_n = kg;						\
		h->g_len = g_len;					\
									\
		if (debug) tref[3] += GMRFLib_cpu();			\
		if (debug) tref[4] -= GMRFLib_cpu();			\
									\
		if (debug_details && debug) {				\
			GMRFLib_idxval_printf(stdout, h, "NEW");	\
			fflush(stdout);					\
		}							\
									\
		/*							\
		 * merge groups together by merging individual groups, and merging sequential groups, possible padding with zeros \
		 */							\
		if (merge_groups) {					\
			int gg = 0;					\
			int irregular = 0;				\
			for (int g = 0; g < h->g_n; g++) {		\
				if (h->g_len[g] == 1) {			\
					if (irregular) {		\
						h->g_len[gg]++;		\
					} else {			\
						h->g_len[gg] = 1;	\
						g_istart[gg] = g_istart[g]; \
					}				\
					irregular = 1;			\
				} else if (IABS(h->g_len[g]) > 0) {	\
					if (irregular) {		\
						gg++;			\
					}				\
					irregular = 0;			\
					h->g_len[gg] = -IABS(h->g_len[g]); \
					g_istart[gg] = g_istart[g];	\
					gg++;				\
				} else {				\
					continue;			\
				}					\
			}						\
			h->g_n = gg + irregular;			\
									\
			int fill = 8L;					\
			if (h && h->n > 1) {				\
				int len_extra = 0;			\
				for (int g = 0; g < h->g_n - 1; g++) {	\
					if (h->g_len[g] < 0 && h->g_len[g + 1] < 0) { \
						int last_this = h->idx[g_istart[g] + IABS(h->g_len[g]) - 1]; \
						int first_next = h->idx[g_istart[g + 1]]; \
						if (first_next - last_this < fill) { \
							len_extra += first_next - last_this - 1; \
							if (debug) {	\
								printf("\tmerge group %1d and %1d with dist %1d and lengths %1d and %1d\n", \
								       g, g + 1, first_next - last_this, h->g_len[g], h->g_len[g + 1]);	\
								fflush(stdout);	\
							}		\
						}			\
					}				\
				}					\
									\
				if (len_extra > 0) {			\
									\
					int n_new = h->n + len_extra;	\
					int *idx_new = Calloc(n_new, int); \
					double *val_new = Calloc(n_new, double); \
									\
					int kk = 0;			\
					int gg = 0;			\
									\
					int this_len = h->g_len[gg];	\
					int this_alen = IABS(this_len);	\
									\
					Memcpy(idx_new + kk, h->idx, this_alen * sizeof(int)); \
					Memcpy(val_new + kk, h->val, this_alen * sizeof(double)); \
									\
					kk += this_alen;		\
					g_istart[gg] = g_istart[0];	\
					h->g_len[gg] = h->g_len[0];	\
									\
					int gg_len = h->g_len[gg];	\
					int gg_alen = IABS(gg_len);	\
					int gg_first_i = g_istart[gg];	\
					int gg_last_i = gg_first_i + gg_alen - 1; \
					int gg_last_idx = idx_new[gg_last_i]; \
					int pending = 0;		\
									\
					if (debug) {			\
						printf("n=%1d len_extra=%1d n_new=%1d\n", h->n, len_extra, n_new); \
						fflush(stdout);		\
					}				\
									\
					for (int g = 1; g < h->g_n; g++) { \
									\
						if (debug) {		\
							printf("\tprocess group=%1d with len=%d kk=%1d\n", g, h->g_len[g], kk);	\
							fflush(stdout);	\
						}			\
									\
						/*			\
						 * merge this group into the current new one or just add it ? \
						 */			\
						int g_len = h->g_len[g]; \
						int g_alen = IABS(g_len); \
						int g_first_i = g_istart[g]; \
						int g_first_idx = h->idx[g_first_i]; \
									\
						int idx_diff = g_first_idx - gg_last_idx; \
						if (debug) {		\
							printf("g_len=%1d gg_len=%1d idx_diff=%1d\n", g_len, gg_len, idx_diff);	\
							fflush(stdout);	\
						}			\
									\
						if (g_len < 0 && gg_len < 0 && idx_diff < fill) { \
							if (debug) {	\
								printf("\tmerge group g=%1d into gg=%1d\n", g, gg); \
								fflush(stdout);	\
							}		\
									\
							/*		\
							 * yes, merge it \
							 */		\
							int len_fill = idx_diff - 1; \
							for (int j = 0; j < len_fill; j++) { \
								idx_new[kk + j] = idx_new[kk - 1] + 1 + j; \
							}		\
							if (debug) {	\
								printf("\t\tlen_fill=%1d\n", len_fill);	\
								fflush(stdout);	\
							}		\
									\
							h->g_len[gg] -= len_fill; \
							kk += len_fill;	\
							Memcpy(idx_new + kk, h->idx + g_first_i, g_alen * sizeof(int));	\
							Memcpy(val_new + kk, h->val + g_first_i, g_alen * sizeof(double)); \
							h->g_len[gg] -= g_alen;	\
							kk += g_alen;	\
									\
							gg_len = h->g_len[gg]; \
							gg_alen = IABS(gg_len);	\
							gg_first_i = g_istart[gg]; \
							gg_last_i = gg_first_i + gg_alen - 1; \
							gg_last_idx = idx_new[gg_last_i]; \
									\
							if (debug) {	\
								printf("\t\tappend group g=%1d with len %1d\n", g, g_len); \
								printf("\t\tgg_last_idx=%1d kk=%1d\n", gg_last_idx, kk); \
								fflush(stdout);	\
							}		\
							pending = 1;	\
						} else {		\
							/*		\
							 * no, just add a new group \
							 */		\
							gg++;		\
							Memcpy(idx_new + kk, h->idx + g_first_i, g_alen * sizeof(int));	\
							Memcpy(val_new + kk, h->val + g_first_i, g_alen * sizeof(double)); \
							h->g_len[gg] = g_len; \
							g_istart[gg] = kk; \
							kk += g_alen;	\
									\
							if (debug) {	\
								printf("\tmake a new group, gg=%1d with len=%1d\n", gg, g_len);	\
								fflush(stdout);	\
							}		\
									\
							gg_len = h->g_len[gg]; \
							gg_alen = IABS(gg_len);	\
							gg_first_i = g_istart[gg]; \
							gg_last_i = gg_first_i + gg_alen - 1; \
							gg_last_idx = idx_new[gg_last_i]; \
							pending = (g < h->g_n - 1 ? 0 : 1); \
						}			\
						if (debug) {		\
							printf("gg=%1d kk=%1d\n", gg, kk); \
							fflush(stdout);	\
						}			\
					}				\
									\
					h->g_n = gg + pending;		\
					h->g_idx = Calloc(h->g_n, int *); \
					h->g_val = Calloc(h->g_n, double *); \
					for (int g = 0; g < h->g_n; g++) { \
						h->g_idx[g] = idx_new + g_istart[g]; \
						h->g_val[g] = val_new + g_istart[g]; \
					}				\
					h->g_n_mem = 2;			\
					h->g_mem = Calloc(h->g_n_mem, void *); \
					h->g_mem[0] = (void *) idx_new;	\
					h->g_mem[1] = (void *) val_new;	\
				} else {				\
					h->g_idx = Calloc(h->g_n, int *); \
					h->g_val = Calloc(h->g_n, double *); \
					for (int g = 0; g < h->g_n; g++) { \
						h->g_idx[g] = h->idx + g_istart[g]; \
						h->g_val[g] = h->val + g_istart[g]; \
					}				\
					h->g_n_mem = 0;			\
				}					\
			}						\
		} else {						\
			h->g_idx = Calloc(h->g_n, int *);		\
			h->g_val = Calloc(h->g_n, double *);		\
			for (int g = 0; g < h->g_n; g++) {		\
				h->g_idx[g] = h->idx + g_istart[g];	\
				h->g_val[g] = h->val + g_istart[g];	\
			}						\
			h->g_n_mem = 0;					\
		}							\
		if (debug) tref[4] += GMRFLib_cpu();			\
		if (debug) tref[5] -= GMRFLib_cpu();			\
									\
		/*							\
		 * add a boolean for all(val[]==1), which makes dot-products into sums \
		 */							\
		int *g_1 = Calloc(h->g_n + 1, int);			\
		if (merge_groups || build_groups) {			\
			for (int g = 0; g < h->g_n; g++) {		\
				int is_1 = 1;				\
				double *val_ptr = h->g_val[g];		\
				int len = IABS(h->g_len[g]);		\
				for (int j = 0; is_1 && j < len; j++) {	\
					is_1 = (val_ptr[j] == 1.0);	\
				}					\
				g_1[g] = (is_1 ? 1 : 0);		\
			}						\
		}							\
		h->g_1 = g_1;						\
		Free(g_istart);						\
									\
		if (debug) tref[5] += GMRFLib_cpu();			\
									\
		if (debug_details && debug) {				\
			GMRFLib_idxval_info_printf(stdout, h, "\t");	\
			fflush(stdout);					\
		}							\
		if (debug) {						\
			double sum = tref[0] + tref[1] + tref[2] + tref[3] + tref[4] + tref[5]; \
			printf("time: [0]%.3f [1]%.3f [2]%.3f [3]%.3f [4]%.3f [5]%.3f\n", tref[0]/sum, tref[1]/sum, tref[2]/sum, tref[3]/sum, tref[4]/sum, tref[5]/sum); \
			fflush(stdout);					\
		}							\
	}

	RUN_CODE_BLOCK(nt, 0, 0);
#undef CODE_BLOCK


	if (!(build_groups && merge_groups)) {
		return GMRFLib_SUCCESS;
	}

	/*
	 * Add a tag about which option that is faster, the group or the serial algorithm, for each 'idxval'. I'm a little reluctant about doing
	 * this within the parallel loop. Usually, this is a very quick procedure, so it does not really matter...
	 */

	nmax = 8;
	for (int i = 0; i < n; i++) {
		GMRFLib_idxval_tp *h = hold[i];
		if (h->n) {
			nmax = IMAX(nmax, h->idx[h->n - 1] + 1);
		}
	}

	double *x = Calloc(nmax, double);
	assert(x);
	for (int i = 0; i < nmax; i++) {
		x[i] = GMRFLib_uniform();
	}

	double time_min = 0.0;
	double time_max = 0.0;
	int ntimes = 2;
#if defined(INLA_LINK_WITH_MKL)
	int with_mkl = 1;
#else
	int with_mkl = 0;
#endif

#define CODE_BLOCK							\
	for (int i = 0; i < n; i++) {					\
		if (hold[i]->preference != IDXVAL_UNKNOWN) {		\
			continue;					\
		}							\
		double tref[4] = { 0.0, 0.0, 0.0, 0.0 };		\
		double value[4] = { 0.0, 0.0, 0.0, 0.0 };		\
		if (debug) {						\
			printf("start testing for hold[%1d]...\n", i);	\
		}							\
		for (int time = 0; time < ntimes; time++) {		\
			int measure = (time >= 0);			\
			if (measure) {					\
				tref[0] -= GMRFLib_cpu();		\
			}						\
			value[0] = GMRFLib_dot_product_serial(hold[i], x); \
			if (measure) {					\
				tref[0] += GMRFLib_cpu();		\
			}						\
			if (with_mkl) {					\
				if (measure) {				\
					tref[1] -= GMRFLib_cpu();	\
				}					\
				value[1] = GMRFLib_dot_product_serial_mkl(hold[i], x); \
				if (measure) {				\
					tref[1] += GMRFLib_cpu();	\
				}					\
			} else {					\
				value[1] = value[0];			\
				tref[1] = tref[0];			\
			}						\
									\
			if (measure) {					\
				tref[2] -= GMRFLib_cpu();		\
			}						\
			value[2] = GMRFLib_dot_product_group(hold[i], x); \
			if (measure) {					\
				tref[2] += GMRFLib_cpu();		\
			}						\
			if (with_mkl) {					\
				if (measure) {				\
					tref[3] -= GMRFLib_cpu();	\
				}					\
				value[3] = GMRFLib_dot_product_group_mkl(hold[i], x); \
				if (measure) {				\
					tref[3] += GMRFLib_cpu();	\
				}					\
			} else {					\
				value[3] = value[2];			\
				tref[3] = tref[2];			\
			}						\
		}							\
									\
		for (int k = 1; k < 4; k++) {				\
			if (ABS(value[k] - value[0]) > FLT_EPSILON * sqrt(hold[i]->n)) { \
				P(ABS(value[k] - value[0]));		\
				P(k);					\
				P(value[0]);				\
				P(value[1]);				\
				P(value[2]);				\
				P(value[3]);				\
				if (0)assert(0 == 1);			\
			}						\
		}							\
									\
		int k = -1;						\
		double tmin = GMRFLib_min_value(tref, 4, &k);		\
		double tmax = GMRFLib_max_value(tref, 4, NULL);		\
									\
		if (debug) {						\
			double s = 1.0 / (tref[0] + tref[1] + tref[2] + tref[3]) / ntimes; \
			printf("for h[%1d] with n= %1d chose k=%1d [serial= %.3f serial.mkl= %.3f group= %.3f group.mkl= %.3f]\n", \
			       i, hold[i]->n, k, tref[0] * s, tref[1] * s, tref[2] * s, tref[3] * s); \
		}							\
									\
		switch (k) {						\
		case 0:							\
			hold[i]->preference = IDXVAL_SERIAL;		\
			break;						\
		case 1:							\
			hold[i]->preference = IDXVAL_SERIAL_MKL;	\
			break;						\
		case 2:							\
			hold[i]->preference = IDXVAL_GROUP;		\
			break;						\
		case 3:							\
			hold[i]->preference = IDXVAL_GROUP_MKL;		\
			break;						\
		default:						\
			assert(0 == 1);					\
		}							\
									\
		if (hold[i]->preference == IDXVAL_SERIAL || hold[i]->preference == IDXVAL_SERIAL_MKL) { \
			/* no need to keep the group info in the struct */ \
			hold[i]->g_n = 0;				\
			Free(hold[i]->g_idx);				\
			Free(hold[i]->g_val);				\
			Free(hold[i]->g_len);				\
			Free(hold[i]->g_1);				\
			for (int k = 0; k < hold[i]->g_n_mem; k++) {	\
				Free(hold[i]->g_mem[k]);		\
			}						\
			Free(hold[i]->g_mem);				\
			hold[i]->g_n_mem = 0;				\
		}							\
									\
		if (GMRFLib_dot_product_optim_report) {			\
			int idx;					\
			GMRFLib_CACHE_SET_ID(idx);			\
			for (int k = 0; k < 4; k++) {			\
				GMRFLib_dot_product_optim_report[idx][k] += tref[k]; \
			}						\
			GMRFLib_dot_product_optim_report[idx][4] += tmin; \
			GMRFLib_dot_product_optim_report[idx][5 + k]++;	/* count... */ \
		}							\
									\
		time_min += tmin / ntimes;				\
		time_max += tmax / ntimes;				\
	}

	RUN_CODE_BLOCK(nt, 0, 0);
#undef CODE_BLOCK

	if (debug) {
		if (time_min > 0.0) {
			printf("idxval opt: saving %.6f seconds/M.evals, %.2f%% improvement\n",
			       (time_max - time_min) * ISQR(1024), 100.0 * (1.0 - time_min / time_max));
		}
	}

	Free(x);
	return GMRFLib_SUCCESS;
}

int GMRFLib_idxval_nsort_x_core_simple(GMRFLib_idxval_tp * h, double *x)
{
	const int limit = 8L;
	const int debug = 0;
	/*
	 * need one to many here. timing only works for num.threads=1 
	 */
	static double tref[] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
	const int measure_time = 0;

#define TIME_SLOT(i_)					\
	if (measure_time) {				\
		if ((i_) > 0) {				\
			tref[(i_) -1] += GMRFLib_cpu();	\
		}					\
		tref[i_] -= GMRFLib_cpu();		\
	}						\

	if (!h || h->n == 0) {
		return GMRFLib_SUCCESS;
	}

	TIME_SLOT(0);

	// sort
	if (h->n > 1) {
		int is_sorted = 1;
		for (int j = 1; is_sorted && j < h->n; j++) {
			is_sorted = (h->idx[j] >= h->idx[j - 1]);
		}
		if (!is_sorted) {
			my_sort2_id(h->idx, h->val, h->n);
		}
	}

	TIME_SLOT(1);

	// unique
	if (h->n > 1) {
		int all_unique = 1;
		for (int j = 1; all_unique && j < h->n; j++) {
			all_unique = (h->idx[j] > h->idx[j - 1]);
		}

		if (!all_unique) {
			int k = 0;
			for (int j = 1; j < h->n; j++) {
				if (h->idx[j] != h->idx[k]) {
					k++;
					h->idx[k] = h->idx[j];
					h->val[k] = h->val[j];
				} else {
					h->val[k] += h->val[j];
				}
			}
			if (debug && (h->n > k + 1)) {
				printf("Make unique: reduce length from %d to %dn", h->n, k + 1);
			}
			h->n = k + 1;
		}
	}

	TIME_SLOT(2);

	if (h->n <= limit) {
		h->preference = IDXVAL_SERIAL_MKL;
		return GMRFLib_SUCCESS;
	}
	// an upper bound for the number of groups for memory allocation
	int ng = 3;
	int i = 1;
	while (i < h->n) {
		while (i < h->n && h->idx[i] == h->idx[i - 1] + 1)
			i++;
		while (i < h->n && h->idx[i] > h->idx[i - 1] + 1)
			i++;
		ng += 2;
	}

	int *g_istart = Calloc(ng, int);
	int *g_len = Calloc(ng, int);
	int *g_1 = Calloc(ng, int);
	int **g_idx = Calloc(ng, int *);
	double **g_val = Calloc(ng, double *);

	TIME_SLOT(3);

	// collect groups, count sequential as 1,3,5,... and then fill in the rest
	ng = 1;
	g_istart[0] = 0;
	g_len[0] = 0;
	i = 1;
	while (i < h->n) {
		if (h->idx[i] == h->idx[i - 1] + 1) {
			while (i < h->n && h->idx[i] == h->idx[i - 1] + 1)
				i++;
			g_len[ng] = -(i - g_istart[ng]);
			ng += 2;
		} else {
			while (i < h->n && h->idx[i] > h->idx[i - 1] + 1)
				i++;
			g_istart[ng] = i - 1;
		}
	}

	g_len[0] = g_istart[1];
	for (int g = 2; g < ng + 1; g += 2) {
		g_istart[g] = g_istart[g - 1] + IABS(g_len[g - 1]);
		g_len[g] = g_istart[g + 1] - g_istart[g];
	}
	g_len[ng - 1] = h->n - g_istart[ng - 1];

	if (debug) {
		for (int i = 0; i < h->n; i++) {
			printf("idx[%1d] =  %1d  val =%g\n", i, h->idx[i], h->val[i]);
		}
		printf("ng = %1d\n", ng);
		for (int g = 0; g < ng; g++) {
			printf("group %1d start %1d len %1d\n", g, g_istart[g], g_len[g]);
			for (int i = 0; i < IABS(g_len[g]); i++) {
				printf("\t\t\tidx %1d\n", h->idx[g_istart[g] + i]);
			}
		}
	}

	TIME_SLOT(4);

	// setup pointers to each sequential group
	for (int g = 0; g < ng; g++) {
		int istart = g_istart[g];
		g_idx[g] = h->idx + istart;
		g_val[g] = h->val + istart;
	}


	// set g_1's
	g_1[0] = 0;
	for (int g = 0; g < ng; g++) {
		int all_one = 1;
		double *val = g_val[g];
		for (int i = 0; all_one && i < IABS(g_len[g]); i++) {
			all_one = (val[i] == 1.0);
		}
		g_1[g] = all_one;
	}

	if (debug) {
		printf("before remove len=0\nng = %1d\n", ng);
		for (int g = 0; g < ng; g++) {
			printf("group %1d start %1d len %1d g_1 %1d\n", g, g_istart[g], g_len[g], g_1[g]);
			for (int i = 0; i < IABS(g_len[g]); i++) {
				printf("\t\t\tidx %1d val %g\n", g_idx[g][i], g_val[g][i]);
			}
		}
	}
	// remove groups with len=0
	int g = 0;
	while (1) {
		if (IABS(g_len[g] == 0)) {
			ng--;
			for (int gg = g; gg < ng; gg++) {
				g_len[gg] = g_len[gg + 1];
				g_istart[gg] = g_istart[gg + 1];
				g_1[gg] = g_1[gg + 1];
				g_idx[gg] = g_idx[gg + 1];
				g_val[gg] = g_val[gg + 1];
			}
		} else {
			g++;
		}

		if (g >= ng - 1) {
			break;
		}
	}

	if (debug) {
		printf("after removing len=0\nng = %1d\n", ng);
		for (int g = 0; g < ng; g++) {
			printf("group %1d start %1d len %1d g_1 %1d\n", g, g_istart[g], g_len[g], g_1[g]);
			for (int i = 0; i < IABS(g_len[g]); i++) {
				printf("\t\t\tidx %1d val %g\n", g_idx[g][i], g_val[g][i]);
			}
		}
	}

	h->g_n = ng;
	h->g_len = g_len;
	h->g_1 = g_1;
	h->g_idx = g_idx;
	h->g_val = g_val;
	h->g_n_mem = 0;
	h->g_mem = NULL;
	Free(g_istart);

	TIME_SLOT(5);

	double time_min = 0.0;
	double time_max = 0.0;
	int ntimes = 1;
#if defined(INLA_LINK_WITH_MKL)
	int with_mkl = 1;
#else
	int with_mkl = 0;
#endif

	if (h->preference != IDXVAL_UNKNOWN) {
		return GMRFLib_SUCCESS;
	}

	double treff[4] = { 0.0, 0.0, 0.0, 0.0 };
	double value[4] = { 0.0, 0.0, 0.0, 0.0 };

	for (int time = 0; time < ntimes; time++) {
		int measure = (time >= 0);
		if (measure) {
			treff[0] -= GMRFLib_cpu();
		}

		value[0] = GMRFLib_dot_product_serial(h, x);
		if (measure) {
			treff[0] += GMRFLib_cpu();
		}
		if (with_mkl) {
			if (measure) {
				treff[1] -= GMRFLib_cpu();
			}
			value[1] = GMRFLib_dot_product_serial_mkl(h, x);
			if (measure) {
				treff[1] += GMRFLib_cpu();
			}
		} else {
			value[1] = value[0];
			treff[1] = treff[0];
		}

		if (measure) {
			treff[2] -= GMRFLib_cpu();
		}
		value[2] = GMRFLib_dot_product_group(h, x);
		if (measure) {
			treff[2] += GMRFLib_cpu();
		}
		if (with_mkl) {
			if (measure) {
				treff[3] -= GMRFLib_cpu();
			}
			value[3] = GMRFLib_dot_product_group_mkl(h, x);
			if (measure) {
				treff[3] += GMRFLib_cpu();
			}
		} else {
			value[3] = value[2];
			treff[3] = treff[2];
		}
	}

	for (int k = 1; k < 4; k++) {
		if (ABS(value[k] - value[0]) > FLT_EPSILON * sqrt(h->n)) {
			P(ABS(value[k] - value[0]));
			P(k);
			P(value[0]);
			P(value[1]);
			P(value[2]);
			P(value[3]);

			double sum1 = 0.0, sum2 = 0.0;
			printf("n %d\n", h->n);
			for (int i = 0; i < h->n; i++) {
				printf("\tidx[%1d] =  %1d  val = %g\n", i, h->idx[i], h->val[i]);
				sum1 += h->val[i] * x[h->idx[i]];
			}
			P(sum1);
			printf("ng %d\n", h->g_n);
			for (int g = 0; g < h->g_n; g++) {
				printf("\tg = %d g_1 = %d g_len = %d\n", g, h->g_1[g], h->g_len[g]);
				for (int i = 0; i < IABS(h->g_len[g]); i++) {
					printf("\t\tidx[%1d] =  %1d  val = %g\n", i, h->g_idx[g][i], h->g_val[g][i]);
					sum2 += h->g_val[g][i] * x[h->g_idx[g][i]];
				}
			}
			P(sum1);
			P(sum2);
			assert(0 == 1);
		}
	}

	int k = -1;
	double tmin = GMRFLib_min_value(treff, 4, &k);
	double tmax = GMRFLib_max_value(treff, 4, NULL);

	if (debug) {
		double s = 1.0 / (treff[0] + treff[1] + treff[2] + treff[3]) / ntimes;
		printf("for h with n= %1d chose k=%1d [serial= %.3f serial.mkl= %.3f group= %.3f group.mkl= %.3f]\n",
		       h->n, k, treff[0] * s, treff[1] * s, treff[2] * s, treff[3] * s);
	}

	switch (k) {
	case 0:
		h->preference = IDXVAL_SERIAL;
		break;
	case 1:
		h->preference = IDXVAL_SERIAL_MKL;
		break;
	case 2:
		h->preference = IDXVAL_GROUP;
		break;
	case 3:
		h->preference = IDXVAL_GROUP_MKL;
		break;
	default:
		assert(0 == 1);
	}

	if (h->preference == IDXVAL_SERIAL || h->preference == IDXVAL_SERIAL_MKL) {
		/*
		 * no need to keep the group info in the struct 
		 */
		h->g_n = 0;
		Free(h->g_idx);
		Free(h->g_val);
		Free(h->g_len);
		Free(h->g_1);
		for (int k = 0; k < h->g_n_mem; k++) {
			Free(h->g_mem[k]);
		}
		Free(h->g_mem);
		h->g_n_mem = 0;
	}

	if (GMRFLib_dot_product_optim_report) {
		int idx;
		GMRFLib_CACHE_SET_ID(idx);
		for (int k = 0; k < 4; k++) {
			GMRFLib_dot_product_optim_report[idx][k] += treff[k];
		}
		GMRFLib_dot_product_optim_report[idx][4] += tmin;
		GMRFLib_dot_product_optim_report[idx][5 + k]++;	/* count... */
	}

	time_min += tmin / ntimes;
	time_max += tmax / ntimes;

	TIME_SLOT(6);

	if (measure_time) {
		double f = 1.0;
		printf("ACCUMULATED TIME: 0:%.3f 1:%.3f 2:%.3f 3:%.3f 4:%.3f 5:%.3f 6:%.3f\n", f * tref[0], f * tref[1],
		       f * tref[2], f * tref[3], f * tref[4], f * tref[5], f * tref[6]);
	}
#undef TIME_SLOT
	return GMRFLib_SUCCESS;
}
