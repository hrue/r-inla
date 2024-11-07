
/* inla-parse.c
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

int inla_parse_lincomb(inla_tp *mb, dictionary *ini, int sec)
{
	/*
	 * parse section = LINCOMB. Here we assume the binary files are written by Rinla, so they are index-1 based!!!!!
	 */
	int *ip = NULL, num_sections, sec_no, n, npairs, offset, i;
	const int debug = 0;
	size_t fileoffset = 0;
	char *filename = NULL, *secname = NULL, *ptr = NULL, *msg = NULL;
	GMRFLib_io_tp *io = NULL;
	GMRFLib_lc_tp *lc = NULL;

	mb->lc_tag = Realloc(mb->lc_tag, mb->nlc + 1, char *);
	mb->lc_output = Realloc(mb->lc_output, mb->nlc + 1, Output_tp *);
	mb->lc_dir = Realloc(mb->lc_dir, mb->nlc + 1, char *);
	mb->lc_order = Realloc(mb->lc_order, mb->nlc + 1, double);
	mb->lc_tag[mb->nlc] = secname = Strdup(iniparser_getsecname(ini, sec));
	mb->lc_dir[mb->nlc] = Strdup(iniparser_getstring(ini, inla_string_join(secname, "DIR"), Strdup(mb->lc_tag[mb->nlc])));

	if (mb->verbose) {
		printf("\tinla_parse_lincomb...\n\t\tsecname = [%s]\n", mb->lc_tag[mb->nlc]);
	}

	filename = Strdup(iniparser_getstring(ini, inla_string_join(secname, "FILENAME"), NULL));
	if (!filename) {
		inla_error_missing_required_field(__GMRFLib_FuncName, secname, "filename");
	}
	fileoffset = (size_t) iniparser_getdouble(ini, inla_string_join(secname, "FILE.OFFSET"), 0.0);

	if (mb->verbose) {
		printf("\t\tfilename [%s]\n", filename);
		printf("\t\tfile.offset [%zu]\n", fileoffset);
	}

	mb->lc_order[mb->nlc] = iniparser_getdouble(ini, inla_string_join(secname, "LINCOMB.ORDER"), -1.0);
	assert(mb->lc_order[mb->nlc] >= 0);

	// FORMAT:: se section.R in Rinla...

	GMRFLib_io_open(&io, filename, "rb");
	if (fileoffset > 0)
		GMRFLib_io_seek(io, fileoffset, SEEK_SET);

	if (mb->verbose) {
		printf("\t\tOpen file [%s] at location [%zu]\n", filename, fileoffset);
	}

	GMRFLib_io_read(io, &num_sections, sizeof(int));
	if (mb->verbose) {
		printf("\t\tNumber of sections [%d]\n", num_sections);
	}

	lc = Calloc(1, GMRFLib_lc_tp);
	lc->n = 0;
	lc->idx = NULL;
	lc->weight = NULL;
	lc->tinfo = Calloc(GMRFLib_CACHE_LEN(), GMRFLib_lc_tinfo_tp);
	for (i = 0; i < GMRFLib_CACHE_LEN(); i++) {
		lc->tinfo[i].first_nonzero = -1;
		lc->tinfo[i].last_nonzero = -1;
		lc->tinfo[i].first_nonzero_mapped = -1;
		lc->tinfo[i].last_nonzero_mapped = -1;
	}

	int all_weights_are_zero = 1;
	for (sec_no = 0; sec_no < num_sections; sec_no++) {

		int len;

		GMRFLib_io_read(io, &len, sizeof(int));
		ptr = Calloc(len + 1, char);
		GMRFLib_io_read(io, ptr, len + 1);	       /* includes trailing \0 */
		if (mb->verbose) {
			printf("\t\t\tSection [%1d] is named [%s]\n", sec_no, ptr);
		}
		ip = map_stri_ptr(&(mb->idx_hash), ptr);
		if (!ip) {
			GMRFLib_sprintf(&msg, "Section no [%1d] named [%s] in file [%1d] offset[%16.0g] is unknown.", sec_no, ptr,
					filename, (double) fileoffset);
			GMRFLib_io_close(io);
			inla_error_general(msg);
		}

		offset = mb->idx_start[*ip];
		n = mb->idx_n[*ip];

		if (mb->verbose)
			printf("\t\t\tSection has offset=[%1d] and n=[%1d]\n", offset, n);

		GMRFLib_io_read(io, &npairs, sizeof(int));
		assert(npairs >= 0);

		if (mb->verbose) {
			printf("\t\t\tnpairs=[%1d]\n", npairs);
		}

		int *idx = Calloc(npairs, int);
		double *w = Calloc(npairs, double);

		GMRFLib_io_read(io, idx, npairs * sizeof(int));
		lc->idx = Realloc(lc->idx, lc->n + npairs, int);
		for (i = 0; i < npairs; i++) {
			lc->idx[lc->n + i] = (idx[i] - 1) + offset;	/* `-1': convert to C-indexing */

			/*
			 * check that the index is legal
			 */
			if (!LEGAL(idx[i] - 1, n)) {
				fprintf(stderr, "\n\n");
				fprintf(stderr, "*** ERROR ***\tLincomb error for section[%s]\n", secname);
				fprintf(stderr, "*** ERROR ***\t[%s] has length %1d, but idx=%1d (R-style index) is given\n", ptr, n, idx[i]);
				GMRFLib_ASSERT(0 == 1, GMRFLib_EPARAMETER);
			}
		}

		GMRFLib_io_read(io, w, npairs * sizeof(double));
		lc->weight = Realloc(lc->weight, lc->n + npairs, double);
		for (i = 0; i < npairs; i++) {
			lc->weight[lc->n + i] = w[i];
			all_weights_are_zero &= (w[i] == 0.0);
		}

		Free(idx);
		Free(w);
		Free(ptr);

		if (debug) {
			for (i = 0; i < npairs; i++) {
				printf("\t\t\t\tC.idx+offset [%1d] weight [%g]\n", lc->idx[lc->n + i], lc->weight[lc->n + i]);
			}
		}

		lc->n += npairs;
	}
	GMRFLib_io_close(io);

	if (all_weights_are_zero) {
		fprintf(stderr, "\n\n");
		fprintf(stderr, "*** ERROR ***\tLincomb error for section[%s]\n", secname);
		fprintf(stderr, "*** ERROR ***\tAll weights are zero.\n");
		GMRFLib_ASSERT(0 == 1, GMRFLib_EPARAMETER);
	}

	/*
	 * sort them with increasing idx's (and carry the weights along) to speed things up later on. 
	 */
	// GMRFLib_qsort2((void *) lc->idx, (size_t) lc->n, sizeof(int), (void *) lc->weight, sizeof(double), NULL, 0, GMRFLib_icmp);
	my_sort2_id(lc->idx, lc->weight, lc->n);
	if (mb->verbose) {
		printf("\t\tNumber of non-zero weights [%1d]\n", lc->n);
		printf("\t\tLincomb = \tidx \tweight\n");
		for (i = 0; i < IMIN(lc->n, PREVIEW); i++) {
			printf("\t\t\t%6d \t\t%.10f\n", lc->idx[i], lc->weight[i]);
		}
	}
	mb->lc_lc = Realloc(mb->lc_lc, (mb->nlc + 1), GMRFLib_lc_tp *);
	mb->lc_lc[mb->nlc] = lc;
	inla_parse_output(mb, ini, sec, &(mb->lc_output[mb->nlc]));
	mb->nlc++;

	return INLA_OK;
}

int inla_parse_mode(inla_tp *mb, dictionary *ini, int sec)
{
	/*
	 * parse section = MODE
	 */
	int nt = 0, i;
	char *tmp = NULL, *secname = NULL;
	double *t = NULL;
	FILE *fp = NULL;
	size_t nread;

	if (mb->verbose) {
		printf("\tinla_parse_mode...\n");
	}
	secname = Strdup(iniparser_getsecname(ini, sec));
	tmp = Strdup(iniparser_getstring(ini, inla_string_join(secname, "THETA"), NULL));

	/*
	 * first try if 'tmp' is a filename, is so, read (using binary format) from that. format: NTHETA theta[0] theta[1] .... theta[ NTHETA-1 ] 
	 */

	if (tmp) {
		fp = fopen(tmp, "rb");
		if (fp) {
			nread = fread(&(mb->ntheta_file), sizeof(int), 1, fp);
			assert(nread == 1);
			mb->theta_file = Calloc(mb->ntheta_file, double);
			nread = fread(mb->theta_file, sizeof(double), mb->ntheta_file, fp);
			assert(nread == (size_t) mb->ntheta_file);
			fclose(fp);
		} else {
			inla_sread_doubles_q(&t, &nt, tmp);
			if (nt) {
				mb->ntheta_file = nt;
				mb->theta_file = t;
			} else {
				mb->ntheta_file = 0;
				mb->theta_file = NULL;
				Free(t);
			}
		}
	} else {
		mb->ntheta_file = 0;
		mb->theta_file = NULL;
	}

	if (mb->verbose) {
		if (mb->ntheta_file) {
			printf("\tUse mode in section[%s]\n", secname);
			printf("\t\ttheta = ");
			for (i = 0; i < mb->ntheta_file; i++) {
				printf(" %.4g", mb->theta_file[i]);
			}
			printf("\n");
		} else {
			printf("\t\tDid not find any theta-mode in section[%s]\n", secname);
		}
	}

	tmp = Strdup(iniparser_getstring(ini, inla_string_join(secname, "X"), NULL));
	if (tmp) {
		/*
		 * this is new code that use binary i/o 
		 */
		// format: NX x[0] x[1] .... x[ NX-1 ]
		fp = fopen(tmp, "rb");
		nread = fread(&(mb->nx_file), sizeof(int), 1, fp);
		assert(nread == 1);
		mb->x_file = Calloc(mb->nx_file, double);
		nread = fread(mb->x_file, sizeof(double), mb->nx_file, fp);
		assert(nread == (size_t) mb->nx_file);
		fclose(fp);

		if (mb->verbose) {
			printf("\t\tx = ");
			for (i = 0; i < IMIN(mb->nx_file, PREVIEW); i++) {
				printf(" %.4g", mb->x_file[i]);
			}
			printf(" ...\n");
		}
	}

	mb->mode_restart = iniparser_getboolean(ini, inla_string_join(secname, "RESTART"), 1);
	if (mb->verbose) {
		printf("\t\tMode_restart = %1d\n", mb->mode_restart);
	}

	mb->mode_fixed = iniparser_getboolean(ini, inla_string_join(secname, "FIXED"), 0);
	if (mb->verbose) {
		printf("\t\tMode_fixed = %1d\n", mb->mode_fixed);
	}

	mb->mode_use_mode = (mb->ntheta_file > 0);
	if (mb->verbose) {
		printf("\t\tMode_use_mode = %1d\n", mb->mode_use_mode);
	}

	return INLA_OK;
}

int inla_parse_problem(inla_tp *mb, dictionary *ini, int sec, int make_dir)
{
	/*
	 * parse section = PROBLEM
	 */
	int i, ok;
	char *secname = NULL, *tmp = NULL, *tmpp = NULL, *smtp = NULL, *openmp_strategy = NULL, *rinla_version = NULL, *build_date = NULL;

	if (mb->verbose) {
		printf("\tinla_parse_problem...\n");
	}
	secname = Strdup(iniparser_getsecname(ini, sec));
	mb->name = Strdup(iniparser_getstring(ini, inla_string_join(secname, "NAME"), NULL));
	if (!mb->name) {
		mb->name = Strdup(secname);
	}
	if (mb->verbose) {
		printf("\t\tname=[%s]\n", mb->name);
	}

	rinla_version = Strdup(iniparser_getstring(ini, inla_string_join(secname, "RINLA.VERSION"), Strdup("UNKNOWN")));
	build_date = Strdup(iniparser_getstring(ini, inla_string_join(secname, "RINLA.BDATE"), Strdup("UNKNOWN")));
	if (mb->verbose) {
		printf("\t\tR-INLA version = [%s]\n", rinla_version);
		printf("\t\tR-INLA build date = [%s]\n", build_date);
		printf("\t\tBuild tag = [%s]\n", INLA_TAG);
		printf("\t\tSystem memory = [%.1fGb]\n", ((double) getTotalSystemMemory()) / 1024.0);
		printf("\t\tCores = (Physical= %1d, Logical= %1d)\n", UTIL_countPhysicalCores(), UTIL_countLogicalCores());

		char a = -1;
		signed char b = -1;
		printf("\t\t'char' is %s\n", (((int) a == (int) b) ? "signed" : "unsigned"));
		printf("\t\tBUFSIZ is %1d\n", BUFSIZ);
		printf("\t\tCACHE_LINE_SIZE is %1d bytes\n", GMRFLib_cachelinesize);

#if defined(__VERSION__)
		printf("\t\tGCC/Compiler version[%s]\n", __VERSION__);
#endif
#if defined(__AVX__)
		printf("\t\tCompiler macro defined [__AVX__]\n");
#endif
#if defined(__AVX2__)
		printf("\t\tCompiler macro defined [__AVX2__]\n");
#endif
#if defined(__AVX512BW__)
		printf("\t\tCompiler macro defined [__AVX512BW__]\n");
#endif
#if defined(__AVX512CD__)
		printf("\t\tCompiler macro defined [__AVX512CD__]\n");
#endif
#if defined(__AVX512DQ__)
		printf("\t\tCompiler macro defined [__AVX512DQ__]\n");
#endif
#if defined(__AVX512F__)
		printf("\t\tCompiler macro defined [__AVX512F__]\n");
#endif
#if defined(__AVX512VL__)
		printf("\t\tCompiler macro defined [__AVX512VL__]\n");
#endif
#if defined(__AVX512VNNI__)
		printf("\t\tCompiler macro defined [__AVX512VNNI__]\n");
#endif
#if defined(__AVXVNNI__)
		printf("\t\tCompiler macro defined [__AVXVNNI__]\n");
#endif
#if defined(__MMX_WITH_SSE__)
		printf("\t\tCompiler macro defined [__MMX_WITH_SSE__]\n");
#endif
#if defined(__SSE__)
		printf("\t\tCompiler macro defined [__SSE__]\n");
#endif
#if defined(__SSE2__)
		printf("\t\tCompiler macro defined [__SSE2__]\n");
#endif
#if defined(__SSE2_MATH__)
		printf("\t\tCompiler macro defined [__SSE2_MATH__]\n");
#endif
#if defined(__SSE3__)
		printf("\t\tCompiler macro defined [__SSE3__]\n");
#endif
#if defined(__SSE4_1__)
		printf("\t\tCompiler macro defined [__SSE4_1__]\n");
#endif
#if defined(__SSE4_2__)
		printf("\t\tCompiler macro defined [__SSE4_2__]\n");
#endif
#if defined(__SSE_MATH__)
		printf("\t\tCompiler macro defined [__SSE_MATH__]\n");
#endif
#if defined(__SSSE3__)
		printf("\t\tCompiler macro defined [__SSSE3__]\n");
#endif
#if defined(INLA_WITH_PARDISO)
		printf("\t\tCompiled with -DINLA_WITH_PARDISO\n");
#endif
#if defined(INLA_WITH_PARDISO_WORKAROUND)
		printf("\t\tCompiled with -DINLA_WITH_PARDISO_WORKAROUND\n");
#endif
#if defined(INLA_WITH_LIBR)
		printf("\t\tCompiled with -DINLA_WITH_LIBR\n");
#endif
#if defined(INLA_WITH_MUPARSER)
		printf("\t\tCompiled with -DINLA_WITH_MUPARSER\n");
#endif
#if defined(INLA_WITH_SIMD)
		printf("\t\tCompiled with -DINLA_WITH_SIMD\n");
#endif
#if defined(INLA_WITH_MKL)
		printf("\t\tCompiled with -DINLA_WITH_MKL\n");
#endif
#if defined(INLA_WITH_OPENBLAS)
		printf("\t\tCompiled with -DINLA_WITH_OPENBLAS\n");
#endif
#if defined(INLA_WITH_ARMPL)
		printf("\t\tCompiled with -DINLA_WITH_ARMPL\n");
#endif
	}


	openmp_strategy = Strdup(iniparser_getstring(ini, inla_string_join(secname, "OPENMP.STRATEGY"), Strdup("DEFAULT")));
	if (mb->verbose) {
		printf("\t\topenmp.strategy=[%s]\n", openmp_strategy);
	}

	if (!strcasecmp(openmp_strategy, "DEFAULT")) {
		/*
		 * this option means that it will be determined later on. 
		 */
		mb->strategy = GMRFLib_OPENMP_STRATEGY_DEFAULT;
	} else if (!strcasecmp(openmp_strategy, "SMALL")) {
		mb->strategy = GMRFLib_OPENMP_STRATEGY_SMALL;
	} else if (!strcasecmp(openmp_strategy, "MEDIUM")) {
		mb->strategy = GMRFLib_OPENMP_STRATEGY_MEDIUM;
	} else if (!strcasecmp(openmp_strategy, "LARGE")) {
		mb->strategy = GMRFLib_OPENMP_STRATEGY_LARGE;
	} else if (!strcasecmp(openmp_strategy, "HUGE")) {
		mb->strategy = GMRFLib_OPENMP_STRATEGY_HUGE;
	} else if (!strcasecmp(openmp_strategy, "PARDISO.SERIAL")) {
		mb->strategy = GMRFLib_OPENMP_STRATEGY_PARDISO;
	} else if (!strcasecmp(openmp_strategy, "PARDISO.PARALLEL")) {
		mb->strategy = GMRFLib_OPENMP_STRATEGY_PARDISO;
	} else if (!strcasecmp(openmp_strategy, "PARDISO.NESTED")) {
		mb->strategy = GMRFLib_OPENMP_STRATEGY_PARDISO;
	} else if (!strcasecmp(openmp_strategy, "PARDISO")) {
		mb->strategy = GMRFLib_OPENMP_STRATEGY_PARDISO;
	} else {
		GMRFLib_sprintf(&tmp, "Unknown openmp.strategy [%s]\n", openmp_strategy);
		inla_error_general(tmp);
		exit(EXIT_FAILURE);
	}

	smtp = Strdup(iniparser_getstring(ini, inla_string_join(secname, "SMTP"),
					  (GMRFLib_openmp->strategy == GMRFLib_OPENMP_STRATEGY_PARDISO ? Strdup("PARDISO") : Strdup("DEFAULT"))));
	if (smtp) {
		if (!strcasecmp(smtp, "BAND")) {
			GMRFLib_smtp = GMRFLib_SMTP_BAND;
		} else if (!strcasecmp(smtp, "TAUCS")) {
			GMRFLib_smtp = GMRFLib_SMTP_TAUCS;
		} else if (!strcasecmp(smtp, "PARDISO")) {
			GMRFLib_smtp = GMRFLib_SMTP_PARDISO;
			mb->strategy = GMRFLib_OPENMP_STRATEGY_PARDISO;
		} else if (!strcasecmp(smtp, "DEFAULT")) {
			if (GMRFLib_pardiso_ok < 0) {
				GMRFLib_pardiso_ok = (GMRFLib_pardiso_check_install(0, 1) == GMRFLib_SUCCESS ? 1 : 0);
			}
			if (GMRFLib_pardiso_ok) {
				if (mb->verbose) {
					printf("\t\tpardiso-library installed and working? = [%s]\n", "yes");
				}
				mb->strategy = GMRFLib_OPENMP_STRATEGY_PARDISO;
				openmp_strategy = Strdup("pardiso");
				GMRFLib_smtp = GMRFLib_SMTP_PARDISO;
				smtp = Strdup("pardiso");
			} else {
				if (mb->verbose) {
					printf("\t\tpardiso-library installed and working? = [%s]\n", "no");
				}
				if (mb->strategy == GMRFLib_OPENMP_STRATEGY_PARDISO) {
					mb->strategy = GMRFLib_OPENMP_STRATEGY_DEFAULT;
					openmp_strategy = Strdup("default");
				}

				GMRFLib_smtp = GMRFLib_SMTP_TAUCS;
				smtp = Strdup("taucs");
			}
		} else {
			inla_error_field_is_void(__GMRFLib_FuncName, secname, "smtp", smtp);
		}
	}
	if (GMRFLib_smtp == GMRFLib_SMTP_PARDISO) {
		GMRFLib_reorder = GMRFLib_REORDER_PARDISO;
		GMRFLib_pardiso_set_parallel_reordering(1);
	}
	mb->smtp = Strdup(GMRFLib_SMTP_NAME(GMRFLib_smtp));
	if (mb->verbose) {
		printf("\t\tsmtp = [%s]\n\t\tstrategy = [%s]\n", smtp, openmp_strategy);
	}
	GMRFLib_openmp_implement_strategy(GMRFLib_OPENMP_PLACES_PARSE_MODEL, NULL, &GMRFLib_smtp);

	mb->dir = Strdup(iniparser_getstring(ini, inla_string_join(secname, "DIR"), Strdup("results-%1d")));
	ok = 0;
	int accept_argument = 0;

	if (make_dir) {
		GMRFLib_sprintf(&tmp, mb->dir, 0);
		GMRFLib_sprintf(&tmpp, mb->dir, 99);
		accept_argument = (strcmp(tmp, tmpp) == 0 ? 0 : 1);
		Free(tmp);
		Free(tmpp);
		for (i = 0; i < 9999; i++) {
			GMRFLib_sprintf(&tmp, mb->dir, i);
			if (inla_mkdir(tmp) != 0) {
				if (mb->verbose) {
					printf("\t\tfail to create directory [%s]: %s\n", tmp, strerror(errno));
				}
				if (!accept_argument) {
					fprintf(stderr, "\n\t\tFail to create directory [%s]: %s\n", tmp, strerror(errno));
					fprintf(stderr, "\t\tmb->dir=[%s] does not accept integer arguments. Cannot proceed.\n\n", mb->dir);
					exit(EXIT_FAILURE);
				}
			} else {
				if (mb->verbose) {
					printf("\t\tstore results in directory=[%s]\n", tmp);
				}
				mb->dir = tmp;
				ok = 1;
				break;
			}
			Free(tmp);
		}
		if (!ok) {
			inla_error_general("Fail to create directory. I give up.");
		}
	}
	GMRFLib_tmpdir = mb->dir;			       /* so we can easily use it elsewhere */

	inla_parse_output(mb, ini, sec, &(mb->output));
	return INLA_OK;
}

int inla_parse_predictor(inla_tp *mb, dictionary *ini, int sec)
{
	/*
	 * parse section = PREDICTOR 
	 */
	char *secname = NULL, *msg = NULL, *filename = NULL;
	int i, noffsets, nlinks_fitted_values;
	double tmp;

	if (mb->verbose) {
		printf("\tinla_parse_predictor ...\n");
	}
	secname = Strdup(iniparser_getsecname(ini, sec));
	mb->predictor_tag = Strdup("Predictor");
	GMRFLib_sprintf(&(mb->Apredictor_tag), "A%s", mb->predictor_tag);

	if (mb->verbose) {
		printf("\t\tsection=[%s]\n", secname);
	}
	mb->predictor_dir = Strdup(iniparser_getstring(ini, inla_string_join(secname, "DIR"), Strdup(mb->predictor_tag)));
	if (mb->verbose) {
		printf("\t\tdir=[%s]\n", mb->predictor_dir);
	}

	inla_read_prior(mb, ini, sec, &(mb->predictor_prior), "LOGGAMMA", NULL);

	tmp = iniparser_getdouble(ini, inla_string_join(secname, "INITIAL"), G.log_prec_initial);
	mb->predictor_fixed = iniparser_getboolean(ini, inla_string_join(secname, "FIXED"), 0);
	if (!mb->predictor_fixed && mb->mode_use_mode) {
		tmp = mb->theta_file[mb->theta_counter_file++];
		if (mb->mode_fixed)
			mb->predictor_fixed = 1;
	}
	HYPER_NEW(mb->predictor_log_prec, tmp);
	if (mb->verbose) {
		printf("\t\tinitialise log_precision[%g]\n", mb->predictor_log_prec[0][0]);
		printf("\t\tfixed=[%1d]\n", mb->predictor_fixed);
	}

	if (!mb->predictor_fixed) {
		mb->theta = Realloc(mb->theta, mb->ntheta + 1, double **);
		mb->theta_hyperid = Realloc(mb->theta_hyperid, mb->ntheta + 1, char *);
		mb->theta_hyperid[mb->ntheta] = mb->predictor_prior.hyperid;
		mb->theta[mb->ntheta] = mb->predictor_log_prec;
		mb->theta_tag = Realloc(mb->theta_tag, mb->ntheta + 1, char *);
		mb->theta_tag_userscale = Realloc(mb->theta_tag_userscale, mb->ntheta + 1, char *);
		mb->theta_dir = Realloc(mb->theta_dir, mb->ntheta + 1, char *);

		mb->theta_tag[mb->ntheta] = Strdup("Log precision for the linear predictor");
		mb->theta_tag_userscale[mb->ntheta] = Strdup("Precision for the linear predictor");
		mb->theta_dir[mb->ntheta] = Strdup(mb->predictor_dir);

		mb->theta_from = Realloc(mb->theta_from, mb->ntheta + 1, char *);
		mb->theta_to = Realloc(mb->theta_to, mb->ntheta + 1, char *);
		mb->theta_from[mb->ntheta] = Strdup(mb->predictor_prior.from_theta);
		mb->theta_to[mb->ntheta] = Strdup(mb->predictor_prior.to_theta);

		mb->theta_map = Realloc(mb->theta_map, mb->ntheta + 1, map_func_tp *);
		mb->theta_map_arg = Realloc(mb->theta_map_arg, mb->ntheta + 1, void *);
		mb->theta_map_arg[mb->ntheta] = NULL;
		mb->theta_map[mb->ntheta] = map_precision;
		mb->ntheta++;
	}
	mb->predictor_user_scale = iniparser_getboolean(ini, inla_string_join(secname, "USER.SCALE"), 1);
	if (mb->verbose) {
		printf("\t\tuser.scale=[%1d]\n", mb->predictor_user_scale);
	}

	mb->predictor_n = iniparser_getint(ini, inla_string_join(secname, "N"), -1);
	assert(mb->predictor_n > 0);
	if (mb->verbose) {
		printf("\t\tn=[%1d]\n", mb->predictor_n);
	}
	mb->predictor_m = iniparser_getint(ini, inla_string_join(secname, "M"), 0);
	assert(mb->predictor_m >= 0);
	if (mb->verbose) {
		printf("\t\tm=[%1d]\n", mb->predictor_m);
	}

	if (mb->predictor_m == 0) {
		mb->predictor_ndata = mb->predictor_n;
	} else {
		mb->predictor_ndata = mb->predictor_m;
	}
	if (mb->verbose) {
		printf("\t\tndata=[%1d]\n", mb->predictor_ndata);
	}
	mb->fl = Calloc(mb->predictor_ndata, int);

	mb->predictor_compute = iniparser_getboolean(ini, inla_string_join(secname, "COMPUTE"), 1);	// mb->output->cpo ||
	// mb->output->dic
	if (G.mode == INLA_MODE_HYPER) {
		if (mb->predictor_compute) {
			fprintf(stderr, "*** Warning: HYPER_MODE require predictor_compute = 0\n");
		}
		mb->predictor_compute = 0;
	}
	if (mb->verbose) {
		printf("\t\tcompute=[%1d]\n", mb->predictor_compute);
	}
	if ((mb->output->cpo || mb->output->dic) && !mb->predictor_compute) {
		GMRFLib_sprintf(&msg, "Illegal combination: output->cpo or dic = 1, require predictor->compute = 1, but predictor->compute = 0");
		inla_error_general(msg);
	}

	filename = Strdup(iniparser_getstring(ini, inla_string_join(secname, "OFFSET"), NULL));
	if (filename) {
		if (mb->verbose) {
			printf("\t\tread offsets from file=[%s]\n", filename);
		}
		inla_read_data_general(&(mb->offset), NULL, &noffsets, filename, mb->predictor_n + mb->predictor_m, 0, 1, mb->verbose, 0.0);
	} else {
		mb->offset = Calloc(mb->predictor_n + mb->predictor_m, double);
	}

	filename = Strdup(iniparser_getstring(ini, inla_string_join(secname, "LINK.FITTED.VALUES"), NULL));
	if (filename) {
		if (mb->verbose) {
			printf("\t\tread link.fitted.values from file=[%s]\n", filename);
		}
		inla_read_data_general(&(mb->link_fitted_values), NULL, &nlinks_fitted_values, filename,
				       mb->predictor_n + mb->predictor_m, 0, 1, mb->verbose, 0.0);
	} else {
		mb->link_fitted_values = NULL;
	}

	if (0) {
		if (mb->link_fitted_values) {
			for (i = 0; i < mb->predictor_n + mb->predictor_m; i++)
				if (gsl_isnan(mb->link_fitted_values[i]))
					fprintf(stderr, "link[%d] = NAN\n", i);
				else
					fprintf(stderr, "link[%d] = %g (%g)\n", i, mb->link_fitted_values[i], GSL_NAN);
		}
	}

	mb->predictor_cross_sumzero = NULL;
	filename = Strdup(iniparser_getstring(ini, inla_string_join(secname, "CROSS_CONSTRAINT"), NULL));
	filename = Strdup(iniparser_getstring(ini, inla_string_join(secname, "CROSS.CONSTRAINT"), filename));
	filename = Strdup(iniparser_getstring(ini, inla_string_join(secname, "CROSSCONSTRAINT"), filename));

	double *dcross = NULL;
	int *icross = NULL, len_cross = 0, nu = 0;

	if (filename) {
		inla_read_data_all(&dcross, &len_cross, filename, NULL);
		if (len_cross > 0) {
			if (len_cross != mb->predictor_n + mb->predictor_m) {
				GMRFLib_sprintf(&msg,
						"Length of cross-sum-to-zero is not equal to the TOTAL length of linear predictor: %1d != %1d\n",
						len_cross, mb->predictor_n + mb->predictor_m);
			}
			icross = Calloc(len_cross, int);
			for (i = 0; i < len_cross; i++)
				icross[i] = (int) dcross[i];
			Free(dcross);
			mb->predictor_cross_sumzero = icross;
		}
	}
	if (mb->verbose && mb->predictor_cross_sumzero) {
		GMRFLib_iuniques(&nu, NULL, mb->predictor_cross_sumzero, mb->predictor_n);
		printf("\t\tread cross-sum-to-zero from file[%s]: %1d constraints\n", filename, nu);
		for (i = 0; i < IMIN(PREVIEW, mb->predictor_n + mb->predictor_m); i++) {
			printf("\t\t\t%1d %1d\n", i, mb->predictor_cross_sumzero[i]);
		}
	}

	/*
	 * these are for the extended observational model. only valid if predictor_m > 0.
	 */
	mb->predictor_A_fnm = Strdup(iniparser_getstring(ini, inla_string_join(secname, "A"), NULL));
	mb->predictor_Aext_fnm = Strdup(iniparser_getstring(ini, inla_string_join(secname, "AEXT"), NULL));
	mb->predictor_Aext_precision = iniparser_getdouble(ini, inla_string_join(secname, "AEXTPRECISION"), 1.0e8);

	if (mb->verbose) {
		printf("\t\tA=[%s]\n", mb->predictor_A_fnm);
		printf("\t\tAext=[%s]\n", mb->predictor_Aext_fnm);
		printf("\t\tAextPrecision=[%.4g]\n", mb->predictor_Aext_precision);
	}
	if (mb->predictor_m > 0) {
		assert(mb->predictor_Aext_fnm != NULL);
	}
	if (mb->predictor_m == 0) {
		assert(mb->predictor_Aext_fnm == NULL);
	}

	inla_parse_output(mb, ini, sec, &(mb->predictor_output));

	return INLA_OK;
}

int inla_parse_data(inla_tp *mb, dictionary *ini, int sec)
{
	/*
	 * parse section = DATA 
	 */

	char *secname = NULL, *msg = NULL, *ctmp = NULL;
	int i, j, found = 0, n_data = (mb->predictor_m > 0 ? mb->predictor_m : mb->predictor_n), discrete_data = 0;
	int beta_delayed_error = 0;
	double tmp;
	Data_section_tp *ds = NULL;

	mb->nds++;
	mb->data_sections = Realloc(mb->data_sections, mb->nds, Data_section_tp);
	ds = &(mb->data_sections[mb->nds - 1]);		       /* shorthand */
	Memset(ds, 0, sizeof(Data_section_tp));

	if (mb->verbose) {
		printf("\tinla_parse_data [section %1d]...\n", mb->nds);
	}
	secname = Strdup(iniparser_getsecname(ini, sec));
	if (mb->verbose) {
		printf("\t\ttag=[%s]\n", secname);
	}
	ds->data_likelihood = Strdup(strupc(iniparser_getstring(ini, inla_string_join(secname, "LIKELIHOOD"), NULL)));
	inla_trim_family(ds->data_likelihood);

	if (mb->verbose) {
		printf("\t\tfamily=[%s]\n", ds->data_likelihood);
	}

	if (!(ds->data_likelihood)) {
		inla_error_field_is_void(__GMRFLib_FuncName, secname, "LIKELIHOOD", ds->data_likelihood);
	}
	if (!strcasecmp(ds->data_likelihood, "GAUSSIAN") || !strcasecmp(ds->data_likelihood, "NORMAL")) {
		ds->loglikelihood = (GMRFLib_logl_tp *) loglikelihood_gaussian;
		ds->data_id = L_GAUSSIAN;
	} else if (!strcasecmp(ds->data_likelihood, "STDGAUSSIAN") || !strcasecmp(ds->data_likelihood, "STDNORMAL")) {
		ds->loglikelihood = (GMRFLib_logl_tp *) loglikelihood_stdgaussian;
		ds->data_id = L_STDGAUSSIAN;
	} else if (!strcasecmp(ds->data_likelihood, "BCGAUSSIAN")) {
		ds->loglikelihood = (GMRFLib_logl_tp *) loglikelihood_bcgaussian;
		ds->data_id = L_BC_GAUSSIAN;
	} else if (!strcasecmp(ds->data_likelihood, "SEM")) {
		ds->loglikelihood = (GMRFLib_logl_tp *) loglikelihood_sem;
		ds->data_id = L_SEM;
	} else if (!strcasecmp(ds->data_likelihood, "exppower")) {
		ds->loglikelihood = (GMRFLib_logl_tp *) loglikelihood_exppower;
		ds->data_id = L_EXPPOWER;
	} else if (!strcasecmp(ds->data_likelihood, "SIMPLEX")) {
		ds->loglikelihood = (GMRFLib_logl_tp *) loglikelihood_simplex;
		ds->data_id = L_SIMPLEX;
	} else if (!strcasecmp(ds->data_likelihood, "IIDGAMMA")) {
		ds->loglikelihood = (GMRFLib_logl_tp *) loglikelihood_iid_gamma;
		ds->data_id = L_IID_GAMMA;
	} else if (!strcasecmp(ds->data_likelihood, "IIDLOGITBETA")) {
		ds->loglikelihood = (GMRFLib_logl_tp *) loglikelihood_iid_logitbeta;
		ds->data_id = L_IID_LOGITBETA;
	} else if (!strcasecmp(ds->data_likelihood, "LOGGAMMAFRAILTY")) {
		ds->loglikelihood = (GMRFLib_logl_tp *) loglikelihood_loggamma_frailty;
		ds->data_id = L_LOGGAMMA_FRAILTY;
	} else if (!strcasecmp(ds->data_likelihood, "LOGISTIC")) {
		ds->loglikelihood = (GMRFLib_logl_tp *) loglikelihood_logistic;
		ds->data_id = L_LOGISTIC;
	} else if (!strcasecmp(ds->data_likelihood, "SKEWNORMAL") || !strcasecmp(ds->data_likelihood, "SN")) {
		ds->loglikelihood = (GMRFLib_logl_tp *) loglikelihood_sn;
		ds->data_id = L_SKEWNORMAL;
	} else if (!strcasecmp(ds->data_likelihood, "GEV")) {
		ds->loglikelihood = (GMRFLib_logl_tp *) loglikelihood_gev;
		ds->data_id = L_GEV;
	} else if (!strcasecmp(ds->data_likelihood, "BGEV")) {
		ds->loglikelihood = (GMRFLib_logl_tp *) loglikelihood_bgev;
		ds->data_id = L_BGEV;
	} else if (!strcasecmp(ds->data_likelihood, "RCPOISSON")) {
		ds->loglikelihood = (GMRFLib_logl_tp *) loglikelihood_rcpoisson;
		ds->data_id = L_RCPOISSON;
	} else if (!strcasecmp(ds->data_likelihood, "TPOISSON")) {
		ds->loglikelihood = (GMRFLib_logl_tp *) loglikelihood_tpoisson;
		ds->data_id = L_TPOISSON;
	} else if (!strcasecmp(ds->data_likelihood, "GGAUSSIAN")) {
		ds->loglikelihood = (GMRFLib_logl_tp *) loglikelihood_ggaussian;
		ds->data_id = L_GGAUSSIAN;
	} else if (!strcasecmp(ds->data_likelihood, "GGAUSSIANS")) {
		ds->loglikelihood = (GMRFLib_logl_tp *) loglikelihood_ggaussianS;
		ds->data_id = L_GGAUSSIANS;
	} else if (!strcasecmp(ds->data_likelihood, "0POISSON")) {
		ds->loglikelihood = (GMRFLib_logl_tp *) loglikelihood_0poisson;
		ds->data_id = L_0POISSON;
	} else if (!strcasecmp(ds->data_likelihood, "0POISSONS")) {
		ds->loglikelihood = (GMRFLib_logl_tp *) loglikelihood_0poissonS;
		ds->data_id = L_0POISSONS;
	} else if (!strcasecmp(ds->data_likelihood, "0BINOMIAL")) {
		ds->loglikelihood = (GMRFLib_logl_tp *) loglikelihood_0binomial;
		ds->data_id = L_0BINOMIAL;
	} else if (!strcasecmp(ds->data_likelihood, "0BINOMIALS")) {
		ds->loglikelihood = (GMRFLib_logl_tp *) loglikelihood_0binomialS;
		ds->data_id = L_0BINOMIALS;
	} else if (!strcasecmp(ds->data_likelihood, "AGAUSSIAN")) {
		ds->loglikelihood = (GMRFLib_logl_tp *) loglikelihood_agaussian;
		ds->data_id = L_AGAUSSIAN;
	} else if (!strcasecmp(ds->data_likelihood, "FL")) {
		ds->loglikelihood = (GMRFLib_logl_tp *) loglikelihood_fl;
		ds->data_id = L_FL;
	} else if (!strcasecmp(ds->data_likelihood, "T")) {
		ds->loglikelihood = (GMRFLib_logl_tp *) loglikelihood_t;
		ds->data_id = L_T;
	} else if (!strcasecmp(ds->data_likelihood, "TSTRATA")) {
		ds->loglikelihood = (GMRFLib_logl_tp *) loglikelihood_tstrata;
		ds->data_id = L_TSTRATA;
	} else if (!strcasecmp(ds->data_likelihood, "BELL")) {
		ds->loglikelihood = (GMRFLib_logl_tp *) loglikelihood_bell;
		ds->data_id = L_BELL;
		discrete_data = 1;
	} else if (!strcasecmp(ds->data_likelihood, "POISSON")) {
		ds->loglikelihood = (GMRFLib_logl_tp *) loglikelihood_poisson;
		ds->data_id = L_POISSON;
		discrete_data = 1;
	} else if (!strcasecmp(ds->data_likelihood, "NPOISSON")) {
		ds->loglikelihood = (GMRFLib_logl_tp *) loglikelihood_npoisson;
		ds->data_id = L_NPOISSON;
		discrete_data = 1;
	} else if (!strcasecmp(ds->data_likelihood, "NZPOISSON")) {
		ds->loglikelihood = (GMRFLib_logl_tp *) loglikelihood_nzpoisson;
		ds->data_id = L_NZPOISSON;
		discrete_data = 1;
	} else if (!strcasecmp(ds->data_likelihood, "XPOISSON")) {
		ds->loglikelihood = (GMRFLib_logl_tp *) loglikelihood_poisson;
		ds->data_id = L_XPOISSON;
		discrete_data = 0;			       /* disable any check. Yes, this is what this one does. */
	} else if (!strcasecmp(ds->data_likelihood, "QCONTPOISSON")) {
		ds->loglikelihood = (GMRFLib_logl_tp *) loglikelihood_qcontpoisson;
		ds->data_id = L_QCONTPOISSON;
		discrete_data = 1;
	} else if (!strcasecmp(ds->data_likelihood, "CONTPOISSON")) {
		ds->loglikelihood = (GMRFLib_logl_tp *) loglikelihood_contpoisson;
		ds->data_id = L_CONTPOISSON;
	} else if (!strcasecmp(ds->data_likelihood, "CENPOISSON")) {
		ds->loglikelihood = (GMRFLib_logl_tp *) loglikelihood_cenpoisson;
		ds->data_id = L_CENPOISSON;
		discrete_data = 1;
	} else if (!strcasecmp(ds->data_likelihood, "CENPOISSON2")) {
		ds->loglikelihood = (GMRFLib_logl_tp *) loglikelihood_cenpoisson2;
		ds->data_id = L_CENPOISSON2;
		discrete_data = 1;
	} else if (!strcasecmp(ds->data_likelihood, "CENNBINOMIAL2")) {
		ds->loglikelihood = (GMRFLib_logl_tp *) loglikelihood_negative_binomial_cen2;
		ds->data_id = L_CENNBINOMIAL2;
		discrete_data = 1;
	} else if (!strcasecmp(ds->data_likelihood, "GAUSSIANJW")) {
		ds->loglikelihood = (GMRFLib_logl_tp *) loglikelihood_gaussianjw;
		ds->data_id = L_GAUSSIANJW;
		discrete_data = 0;
	} else if (!strcasecmp(ds->data_likelihood, "GPOISSON")) {
		ds->loglikelihood = (GMRFLib_logl_tp *) loglikelihood_gpoisson;
		ds->data_id = L_GPOISSON;
		discrete_data = 1;
	} else if (!strcasecmp(ds->data_likelihood, "POISSONSPECIAL1")) {
		ds->loglikelihood = (GMRFLib_logl_tp *) loglikelihood_poisson_special1;
		ds->data_id = L_POISSON_SPECIAL1;
		discrete_data = 1;
	} else if (!strcasecmp(ds->data_likelihood, "ZEROINFLATEDPOISSON0")) {
		ds->loglikelihood = (GMRFLib_logl_tp *) loglikelihood_zeroinflated_poisson0;
		ds->data_id = L_ZEROINFLATEDPOISSON0;
		discrete_data = 1;
	} else if (!strcasecmp(ds->data_likelihood, "ZEROINFLATEDPOISSON1")) {
		ds->loglikelihood = (GMRFLib_logl_tp *) loglikelihood_zeroinflated_poisson1;
		ds->data_id = L_ZEROINFLATEDPOISSON1;
		discrete_data = 1;
	} else if (!strcasecmp(ds->data_likelihood, "ZEROINFLATEDPOISSON2")) {
		ds->loglikelihood = (GMRFLib_logl_tp *) loglikelihood_zeroinflated_poisson2;
		ds->data_id = L_ZEROINFLATEDPOISSON2;
		discrete_data = 1;
	} else if (!strcasecmp(ds->data_likelihood, "ZEROINFLATEDCENPOISSON0")) {
		ds->loglikelihood = (GMRFLib_logl_tp *) loglikelihood_zeroinflated_cenpoisson0;
		ds->data_id = L_ZEROINFLATEDCENPOISSON0;
		discrete_data = 1;
	} else if (!strcasecmp(ds->data_likelihood, "ZEROINFLATEDCENPOISSON1")) {
		ds->loglikelihood = (GMRFLib_logl_tp *) loglikelihood_zeroinflated_cenpoisson1;
		ds->data_id = L_ZEROINFLATEDCENPOISSON1;
		discrete_data = 1;
	} else if (!strcasecmp(ds->data_likelihood, "BINOMIALMIX")) {
		ds->loglikelihood = (GMRFLib_logl_tp *) loglikelihood_binomialmix;
		ds->data_id = L_BINOMIALMIX;
		discrete_data = 1;
	} else if (!strcasecmp(ds->data_likelihood, "BINOMIAL")) {
		ds->loglikelihood = (GMRFLib_logl_tp *) loglikelihood_binomial;
		ds->data_id = L_BINOMIAL;
		discrete_data = 1;
	} else if (!strcasecmp(ds->data_likelihood, "XBINOMIAL")) {
		ds->loglikelihood = (GMRFLib_logl_tp *) loglikelihood_xbinomial;
		ds->data_id = L_XBINOMIAL;
		discrete_data = 0;
	} else if (!strcasecmp(ds->data_likelihood, "NBINOMIAL2")) {
		ds->loglikelihood = (GMRFLib_logl_tp *) loglikelihood_nbinomial2;
		ds->data_id = L_NBINOMIAL2;
		discrete_data = 1;
	} else if (!strcasecmp(ds->data_likelihood, "CBINOMIAL")) {
		ds->loglikelihood = (GMRFLib_logl_tp *) loglikelihood_cbinomial;
		ds->data_id = L_CBINOMIAL;
		discrete_data = 1;
	} else if (!strcasecmp(ds->data_likelihood, "POM")) {
		ds->loglikelihood = (GMRFLib_logl_tp *) loglikelihood_pom;
		ds->data_id = L_POM;
		discrete_data = 1;
	} else if (!strcasecmp(ds->data_likelihood, "GAMMA")) {
		ds->loglikelihood = (GMRFLib_logl_tp *) loglikelihood_gamma;
		ds->data_id = L_GAMMA;
	} else if (!strcasecmp(ds->data_likelihood, "MGAMMA")) {
		ds->loglikelihood = (GMRFLib_logl_tp *) loglikelihood_mgamma;
		ds->data_id = L_MGAMMA;
	} else if (!strcasecmp(ds->data_likelihood, "GAMMAJW")) {
		ds->loglikelihood = (GMRFLib_logl_tp *) loglikelihood_gammajw;
		ds->data_id = L_GAMMAJW;
	} else if (!strcasecmp(ds->data_likelihood, "GAMMACOUNT")) {
		ds->loglikelihood = (GMRFLib_logl_tp *) loglikelihood_gammacount;
		ds->data_id = L_GAMMACOUNT;
		discrete_data = 1;
	} else if (!strcasecmp(ds->data_likelihood, "QKUMAR")) {
		ds->loglikelihood = (GMRFLib_logl_tp *) loglikelihood_qkumar;
		ds->data_id = L_QKUMAR;
	} else if (!strcasecmp(ds->data_likelihood, "qloglogistic")) {
		ds->loglikelihood = (GMRFLib_logl_tp *) loglikelihood_qloglogistic;
		ds->data_id = L_QLOGLOGISTIC;
	} else if (!strcasecmp(ds->data_likelihood, "qloglogisticsurv")) {
		ds->loglikelihood = (GMRFLib_logl_tp *) loglikelihood_qloglogisticsurv;
		ds->data_id = L_QLOGLOGISTICSURV;
	} else if (!strcasecmp(ds->data_likelihood, "BETA")) {
		ds->loglikelihood = (GMRFLib_logl_tp *) loglikelihood_beta;
		ds->data_id = L_BETA;
	} else if (!strcasecmp(ds->data_likelihood, "BETABINOMIAL")) {
		ds->loglikelihood = (GMRFLib_logl_tp *) loglikelihood_betabinomial;
		ds->data_id = L_BETABINOMIAL;
		discrete_data = 1;
	} else if (!strcasecmp(ds->data_likelihood, "BETABINOMIALNA")) {
		ds->loglikelihood = (GMRFLib_logl_tp *) loglikelihood_betabinomialna;
		ds->data_id = L_BETABINOMIALNA;
	} else if (!strcasecmp(ds->data_likelihood, "ZEROINFLATEDBINOMIAL0")) {
		ds->loglikelihood = (GMRFLib_logl_tp *) loglikelihood_zeroinflated_binomial0;
		ds->data_id = L_ZEROINFLATEDBINOMIAL0;
		discrete_data = 1;
	} else if (!strcasecmp(ds->data_likelihood, "ZEROINFLATEDBINOMIAL1")) {
		ds->loglikelihood = (GMRFLib_logl_tp *) loglikelihood_zeroinflated_binomial1;
		ds->data_id = L_ZEROINFLATEDBINOMIAL1;
		discrete_data = 1;
	} else if (!strcasecmp(ds->data_likelihood, "ZEROINFLATEDBINOMIAL2")) {
		ds->loglikelihood = (GMRFLib_logl_tp *) loglikelihood_zeroinflated_binomial2;
		ds->data_id = L_ZEROINFLATEDBINOMIAL2;
		discrete_data = 1;
	} else if (!strcasecmp(ds->data_likelihood, "ZERONINFLATEDBINOMIAL2")) {
		ds->loglikelihood = (GMRFLib_logl_tp *) loglikelihood_zero_n_inflated_binomial2;
		ds->data_id = L_ZERO_N_INFLATEDBINOMIAL2;
		discrete_data = 1;
	} else if (!strcasecmp(ds->data_likelihood, "ZERONINFLATEDBINOMIAL3")) {
		ds->loglikelihood = (GMRFLib_logl_tp *) loglikelihood_zero_n_inflated_binomial3;
		ds->data_id = L_ZERO_N_INFLATEDBINOMIAL3;
		discrete_data = 1;
	} else if (!strcasecmp(ds->data_likelihood, "ZEROINFLATEDBETABINOMIAL0")) {
		ds->loglikelihood = (GMRFLib_logl_tp *) loglikelihood_zeroinflated_betabinomial0;
		ds->data_id = L_ZEROINFLATEDBETABINOMIAL0;
		discrete_data = 1;
	} else if (!strcasecmp(ds->data_likelihood, "ZEROINFLATEDBETABINOMIAL1")) {
		ds->loglikelihood = (GMRFLib_logl_tp *) loglikelihood_zeroinflated_betabinomial1;
		ds->data_id = L_ZEROINFLATEDBETABINOMIAL1;
		discrete_data = 1;
	} else if (!strcasecmp(ds->data_likelihood, "ZEROINFLATEDBETABINOMIAL2")) {
		ds->loglikelihood = (GMRFLib_logl_tp *) loglikelihood_zeroinflated_betabinomial2;
		ds->data_id = L_ZEROINFLATEDBETABINOMIAL2;
		discrete_data = 1;
	} else if (!strcasecmp(ds->data_likelihood, "NBINOMIAL")) {
		ds->loglikelihood = (GMRFLib_logl_tp *) loglikelihood_negative_binomial;
		ds->data_id = L_NBINOMIAL;
		discrete_data = 1;
	} else if (!strcasecmp(ds->data_likelihood, "ZEROINFLATEDNBINOMIAL0")) {
		ds->loglikelihood = (GMRFLib_logl_tp *) loglikelihood_zeroinflated_negative_binomial0;
		ds->data_id = L_ZEROINFLATEDNBINOMIAL0;
		discrete_data = 1;
	} else if (!strcasecmp(ds->data_likelihood, "ZEROINFLATEDNBINOMIAL1")) {
		ds->loglikelihood = (GMRFLib_logl_tp *) loglikelihood_zeroinflated_negative_binomial1;
		ds->data_id = L_ZEROINFLATEDNBINOMIAL1;
		discrete_data = 1;
	} else if (!strcasecmp(ds->data_likelihood, "ZEROINFLATEDNBINOMIAL2")) {
		ds->loglikelihood = (GMRFLib_logl_tp *) loglikelihood_zeroinflated_negative_binomial2;
		ds->data_id = L_ZEROINFLATEDNBINOMIAL2;
		discrete_data = 1;
	} else if (!strcasecmp(ds->data_likelihood, "ZEROINFLATEDNBINOMIAL1STRATA2")) {
		ds->loglikelihood = (GMRFLib_logl_tp *) loglikelihood_zeroinflated_negative_binomial1_strata2;
		ds->data_id = L_ZEROINFLATEDNBINOMIAL1STRATA2;
		discrete_data = 1;
	} else if (!strcasecmp(ds->data_likelihood, "ZEROINFLATEDNBINOMIAL1STRATA3")) {
		ds->loglikelihood = (GMRFLib_logl_tp *) loglikelihood_zeroinflated_negative_binomial1_strata3;
		ds->data_id = L_ZEROINFLATEDNBINOMIAL1STRATA3;
		discrete_data = 1;
	} else if (!strcasecmp(ds->data_likelihood, "STOCHVOL")) {
		ds->loglikelihood = (GMRFLib_logl_tp *) loglikelihood_stochvol;
		ds->data_id = L_STOCHVOL;
	} else if (!strcasecmp(ds->data_likelihood, "STOCHVOLLN")) {
		ds->loglikelihood = (GMRFLib_logl_tp *) loglikelihood_stochvolln;
		ds->data_id = L_STOCHVOL_LN;
	} else if (!strcasecmp(ds->data_likelihood, "STOCHVOLSN")) {
		ds->loglikelihood = (GMRFLib_logl_tp *) loglikelihood_stochvol_sn;
		ds->data_id = L_STOCHVOL_SN;
	} else if (!strcasecmp(ds->data_likelihood, "STOCHVOLT")) {
		ds->loglikelihood = (GMRFLib_logl_tp *) loglikelihood_stochvol_t;
		ds->data_id = L_STOCHVOL_T;
	} else if (!strcasecmp(ds->data_likelihood, "STOCHVOLNIG")) {
		ds->loglikelihood = (GMRFLib_logl_tp *) loglikelihood_stochvol_nig;
		ds->data_id = L_STOCHVOL_NIG;
	} else if (!strcasecmp(ds->data_likelihood, "LOGPERIODOGRAM")) {
		ds->loglikelihood = (GMRFLib_logl_tp *) loglikelihood_logperiodogram;
		ds->data_id = L_LOGPERIODOGRAM;
	} else if (!strcasecmp(ds->data_likelihood, "EXPONENTIAL")) {
		ds->loglikelihood = (GMRFLib_logl_tp *) loglikelihood_exp;
		ds->data_id = L_EXPONENTIAL;
	} else if (!strcasecmp(ds->data_likelihood, "EXPONENTIALSURV")) {
		ds->loglikelihood = (GMRFLib_logl_tp *) loglikelihood_expsurv;
		ds->data_id = L_EXPONENTIALSURV;
	} else if (!strcasecmp(ds->data_likelihood, "GAMMASURV")) {
		ds->loglikelihood = (GMRFLib_logl_tp *) loglikelihood_gammasurv;
		ds->data_id = L_GAMMASURV;
	} else if (!strcasecmp(ds->data_likelihood, "MGAMMASURV")) {
		ds->loglikelihood = (GMRFLib_logl_tp *) loglikelihood_mgammasurv;
		ds->data_id = L_MGAMMASURV;
	} else if (!strcasecmp(ds->data_likelihood, "GAMMAJWSURV")) {
		ds->loglikelihood = (GMRFLib_logl_tp *) loglikelihood_gammajwsurv;
		ds->data_id = L_GAMMAJWSURV;
	} else if (!strcasecmp(ds->data_likelihood, "WEIBULL")) {
		ds->loglikelihood = (GMRFLib_logl_tp *) loglikelihood_weibull;
		ds->data_id = L_WEIBULL;
	} else if (!strcasecmp(ds->data_likelihood, "WEIBULLSURV")) {
		ds->loglikelihood = (GMRFLib_logl_tp *) loglikelihood_weibullsurv;
		ds->data_id = L_WEIBULLSURV;
	} else if (!strcasecmp(ds->data_likelihood, "GOMPERTZ")) {
		ds->loglikelihood = (GMRFLib_logl_tp *) loglikelihood_gompertz;
		ds->data_id = L_GOMPERTZ;
	} else if (!strcasecmp(ds->data_likelihood, "GOMPERTZSURV")) {
		ds->loglikelihood = (GMRFLib_logl_tp *) loglikelihood_gompertzsurv;
		ds->data_id = L_GOMPERTZSURV;
	} else if (!strcasecmp(ds->data_likelihood, "LOGLOGISTIC")) {
		ds->loglikelihood = (GMRFLib_logl_tp *) loglikelihood_loglogistic;
		ds->data_id = L_LOGLOGISTIC;
	} else if (!strcasecmp(ds->data_likelihood, "LOGLOGISTICSURV")) {
		ds->loglikelihood = (GMRFLib_logl_tp *) loglikelihood_loglogisticsurv;
		ds->data_id = L_LOGLOGISTICSURV;
	} else if (!strcasecmp(ds->data_likelihood, "LOGNORMALSURV")) {
		ds->loglikelihood = (GMRFLib_logl_tp *) loglikelihood_lognormalsurv;
		ds->data_id = L_LOGNORMALSURV;
	} else if (!strcasecmp(ds->data_likelihood, "LOGNORMAL")) {
		ds->loglikelihood = (GMRFLib_logl_tp *) loglikelihood_lognormal;
		ds->data_id = L_LOGNORMAL;
	} else if (!strcasecmp(ds->data_likelihood, "CIRCULARNORMAL")) {
		ds->loglikelihood = (GMRFLib_logl_tp *) loglikelihood_circular_normal;
		ds->data_id = L_CIRCULAR_NORMAL;
	} else if (!strcasecmp(ds->data_likelihood, "WRAPPEDCAUCHY")) {
		ds->loglikelihood = (GMRFLib_logl_tp *) loglikelihood_wrapped_cauchy;
		ds->data_id = L_WRAPPED_CAUCHY;
	} else if (!strcasecmp(ds->data_likelihood, "TWEEDIE")) {
		ds->loglikelihood = (GMRFLib_logl_tp *) loglikelihood_tweedie;
		ds->data_id = L_TWEEDIE;
	} else if (!strcasecmp(ds->data_likelihood, "FMRI")) {
		ds->loglikelihood = (GMRFLib_logl_tp *) loglikelihood_fmri;
		ds->data_id = L_FMRI;
	} else if (!strcasecmp(ds->data_likelihood, "FMRISURV")) {
		ds->loglikelihood = (GMRFLib_logl_tp *) loglikelihood_fmrisurv;
		ds->data_id = L_FMRISURV;
	} else if (!strcasecmp(ds->data_likelihood, "GP")) {
		ds->loglikelihood = (GMRFLib_logl_tp *) loglikelihood_gp;
		ds->data_id = L_GP;
	} else if (!strcasecmp(ds->data_likelihood, "EGP")) {
		ds->loglikelihood = (GMRFLib_logl_tp *) loglikelihood_egp;
		ds->data_id = L_EGP;
	} else if (!strcasecmp(ds->data_likelihood, "DGP")) {
		ds->loglikelihood = (GMRFLib_logl_tp *) loglikelihood_dgp;
		ds->data_id = L_DGP;
		discrete_data = 1;
	} else if (!strcasecmp(ds->data_likelihood, "NMIX")) {
		ds->loglikelihood = (GMRFLib_logl_tp *) loglikelihood_nmix;
		ds->data_id = L_NMIX;
		discrete_data = 1;
	} else if (!strcasecmp(ds->data_likelihood, "NMIXNB")) {
		ds->loglikelihood = (GMRFLib_logl_tp *) loglikelihood_nmixnb;
		ds->data_id = L_NMIXNB;
		discrete_data = 1;
	} else if (!strcasecmp(ds->data_likelihood, "OCCUPANCY")) {
		ds->loglikelihood = (GMRFLib_logl_tp *) loglikelihood_occupancy;
		ds->data_id = L_OCCUPANCY;
		discrete_data = 1;
	} else {
		FIXME("FOUND");
		inla_error_field_is_void(__GMRFLib_FuncName, secname, "LIKELIHOOD", ds->data_likelihood);
	}
	if (mb->verbose) {
		printf("\t\tlikelihood=[%s]\n", ds->data_likelihood);
	}

	inla_read_fileinfo(mb, ini, sec, &(ds->data_file), NULL);
	inla_read_fileinfo(mb, ini, sec, &(ds->weight_file), "WEIGHTS");
	inla_read_fileinfo(mb, ini, sec, &(ds->attr_file), "ATTRIBUTES");
	inla_read_fileinfo(mb, ini, sec, &(ds->lp_scale_file), "LPSCALE");
	inla_read_data_likelihood(mb, ini, sec);

	/*
	 * validate the data 
	 */
	if (discrete_data) {
		for (i = 0; i < mb->predictor_ndata; i++) {
			if (ds->data_observations.d[i]) {
				if (ds->data_observations.y[i] != (int) ds->data_observations.y[i]) {
					GMRFLib_sprintf(&msg, "%s: %s likelihood is defined on integers, but y[%1d] = %.12g\n",
							secname, ds->data_likelihood, i, ds->data_observations.y[i]);
					inla_error_general(msg);
				}
			}
		}
	}

	switch (ds->data_id) {
	case L_SEM:
		break;

	case L_GAUSSIAN:
	case L_STDGAUSSIAN:
	{
		for (i = 0; i < mb->predictor_ndata; i++) {
			if (ds->data_observations.d[i]) {
				if (ds->data_observations.weight_gaussian[i] <= 0.0) {
					GMRFLib_sprintf(&msg, "%s: Gaussian weight[%1d] = %g is void\n", secname, i,
							ds->data_observations.weight_gaussian[i]);
					inla_error_general(msg);
				}
			}
		}
	}
		break;

	case L_EXPPOWER:
	{
		for (i = 0; i < mb->predictor_ndata; i++) {
			if (ds->data_observations.d[i]) {
				if (ds->data_observations.weight_gaussian[i] <= 0.0) {
					GMRFLib_sprintf(&msg, "%s: Exponential Power weight[%1d] = %g is void\n", secname, i,
							ds->data_observations.weight_gaussian[i]);
					inla_error_general(msg);
				}
			}
		}
	}
		break;

	case L_BC_GAUSSIAN:
	{
		for (i = 0; i < mb->predictor_ndata; i++) {
			if (ds->data_observations.d[i]) {
				if (ds->data_observations.bc_scale[i] <= 0.0 || ds->data_observations.y[i] <= 0.0
				    || ds->data_observations.bc_mean[i] <= 0.0) {
					GMRFLib_sprintf(&msg, "%s: Box-Cox Gaussian scale[%1d] = %g or y[%1d] = %g or mean[%1d] = %g is void\n",
							secname, i, ds->data_observations.bc_scale[i], i, ds->data_observations.y[i], i,
							ds->data_observations.bc_mean[i]);
					inla_error_general(msg);
				}
			}
		}
	}
		break;

	case L_GGAUSSIAN:
	{
		for (i = 0; i < mb->predictor_ndata; i++) {
			if (ds->data_observations.d[i]) {
				if (ds->data_observations.ggaussian_scale[i] <= 0.0) {
					GMRFLib_sprintf(&msg, "%s: GGaussian scale[%1d] = %g is void\n", secname, i,
							ds->data_observations.ggaussian_scale[i]);
					inla_error_general(msg);
				}
			}
		}
	}
		break;

	case L_GGAUSSIANS:
		break;

	case L_RCPOISSON:
	{
		for (i = 0; i < mb->predictor_ndata; i++) {
			if (ds->data_observations.d[i]) {
				if (ds->data_observations.y[i] < 0.0) {
					GMRFLib_sprintf(&msg, "%s: rcpoisson y[%1d]=[%g] is void\n", secname, i, ds->data_observations.y[i]);
					inla_error_general(msg);
				}
			}
		}
	}
		break;

	case L_TPOISSON:
	{
		for (i = 0; i < mb->predictor_ndata; i++) {
			if (ds->data_observations.d[i]) {
				if (ds->data_observations.y[i] < 0.0) {
					GMRFLib_sprintf(&msg, "%s: tpoisson y[%1d]=[%g] is void\n", secname, i, ds->data_observations.y[i]);
					inla_error_general(msg);
				}
			}
		}
	}
		break;

	case L_SIMPLEX:
	{
		for (i = 0; i < mb->predictor_ndata; i++) {
			if (ds->data_observations.d[i]) {
				if (ds->data_observations.weight_simplex[i] <= 0.0) {
					GMRFLib_sprintf(&msg, "%s: Simplex weight[%1d] = %g is void\n", secname, i,
							ds->data_observations.weight_simplex[i]);
					inla_error_general(msg);
				}
				if (ds->data_observations.y[i] <= 0.0 || ds->data_observations.y[i] >= 1) {
					GMRFLib_sprintf(&msg, "%s: Simplex data[%1d] (y) = (%g) is void\n", secname, i, ds->data_observations.y[i]);
					inla_error_general(msg);
				}
			}
		}
	}
		break;

	case L_CIRCULAR_NORMAL:
	{
		for (i = 0; i < mb->predictor_ndata; i++) {
			if (ds->data_observations.d[i]) {
				if (ds->data_observations.weight_circular_normal[i] <= 0.0) {
					GMRFLib_sprintf(&msg, "%s: Circular Normal weight[%1d] = %g is void\n", secname, i,
							ds->data_observations.weight_circular_normal[i]);
					inla_error_general(msg);
				}
			}
			if (ds->data_observations.d[i]) {
				if (ABS(ds->data_observations.y[i]) > 2.0 * M_PI) {
					GMRFLib_sprintf(&msg, "%s: Circular Normal observation y[%1d] = %g is void\n", secname, i,
							ds->data_observations.y[i]);
					inla_error_general(msg);
				}
			}
		}
	}
		break;

	case L_WRAPPED_CAUCHY:
	{
		for (i = 0; i < mb->predictor_ndata; i++) {
			if (ds->data_observations.d[i]) {
				if (ds->data_observations.weight_wrapped_cauchy[i] <= 0.0 || ds->data_observations.weight_wrapped_cauchy[i] > 1.0) {
					GMRFLib_sprintf(&msg, "%s: Wrapped Cauchy weight[%1d] = %g is void\n", secname, i,
							ds->data_observations.weight_wrapped_cauchy[i]);
					inla_error_general(msg);
				}
			}
			if (ds->data_observations.d[i]) {
				if (ABS(ds->data_observations.y[i]) > 2.0 * M_PI) {
					GMRFLib_sprintf(&msg, "%s: Wrapped Cauchy observation y[%1d] = %g is void\n", secname, i,
							ds->data_observations.y[i]);
					inla_error_general(msg);
				}
			}
		}
	}
		break;

	case L_TWEEDIE:
	{
		for (i = 0; i < mb->predictor_ndata; i++) {
			if (ds->data_observations.d[i]) {
				if (ds->data_observations.tweedie_w[i] <= 0.0) {
					GMRFLib_sprintf(&msg, "%s: Tweedie scale[%1d] = %g is void\n", secname, i,
							ds->data_observations.tweedie_w[i]);
					inla_error_general(msg);
				}
			}
			if (ds->data_observations.d[i]) {
				if (ABS(ds->data_observations.y[i]) < 0.0) {
					GMRFLib_sprintf(&msg, "%s: Tweedie observation y[%1d] = %g is void\n", secname, i,
							ds->data_observations.y[i]);
					inla_error_general(msg);
				}
			}
		}
	}
		break;

	case L_FMRI:
	{
		for (i = 0; i < mb->predictor_ndata; i++) {
			if (ds->data_observations.d[i]) {
				if (ds->data_observations.fmri_scale[i] <= 0.0) {
					GMRFLib_sprintf(&msg, "%s: fmri scale[%1d] = %g is void\n", secname, i,
							ds->data_observations.fmri_scale[i]);
					inla_error_general(msg);
				}
			}
			if (ds->data_observations.d[i]) {
				if (ABS(ds->data_observations.y[i]) <= 0.0) {
					GMRFLib_sprintf(&msg, "%s: fmri observation y[%1d] = %g is void\n", secname, i, ds->data_observations.y[i]);
					inla_error_general(msg);
				}
			}
		}
	}
		break;

	case L_GP:
	case L_DGP:
	case L_EGP:
	{
		for (i = 0; i < mb->predictor_ndata; i++) {
			if (ds->data_observations.d[i]) {
				if (ds->data_observations.y[i] < 0.0) {
					GMRFLib_sprintf(&msg, "%s: gp/dgp/egp observation y[%1d] = %g is void\n", secname, i,
							ds->data_observations.y[i]);
					inla_error_general(msg);
				}
			}
		}
	}
		break;

	case L_OCCUPANCY:
		// need to do this later on
		break;

	case L_IID_GAMMA:
	case L_IID_LOGITBETA:
	case L_LOGGAMMA_FRAILTY:
	{
		/*
		 * ok...
		 */
	}
		break;

	case L_LOGISTIC:
	{
		for (i = 0; i < mb->predictor_ndata; i++) {
			if (ds->data_observations.d[i]) {
				if (ds->data_observations.weight_logistic[i] <= 0.0) {
					GMRFLib_sprintf(&msg, "%s: Logistic weight[%1d] = %g is void\n", secname, i,
							ds->data_observations.weight_logistic[i]);
					inla_error_general(msg);
				}
			}
		}
	}
		break;

	case L_SKEWNORMAL:
	{
		for (i = 0; i < mb->predictor_ndata; i++) {
			if (ds->data_observations.d[i]) {
				if (ds->data_observations.sn_scale[i] <= 0.0) {
					GMRFLib_sprintf(&msg, "%s: Skewnormal scale[%1d] = %g is void\n", secname, i,
							ds->data_observations.sn_scale[i]);
					inla_error_general(msg);
				}
			}
		}
	}
		break;

	case L_GEV:
	{
		for (i = 0; i < mb->predictor_ndata; i++) {
			if (ds->data_observations.d[i]) {
				if (ds->data_observations.weight_gev[i] <= 0.0) {
					GMRFLib_sprintf(&msg, "%s: GEV weight[%1d] = %g is void\n", secname, i,
							ds->data_observations.weight_gev[i]);
					inla_error_general(msg);
				}
			}
		}
	}
		break;

	case L_BGEV:
	{
		for (i = 0; i < mb->predictor_ndata; i++) {
			if (ds->data_observations.d[i]) {
				if (ds->data_observations.bgev_scale[i] <= 0.0) {
					GMRFLib_sprintf(&msg, "%s: BGEV scale[%1d] = %g is void\n", secname, i,
							ds->data_observations.bgev_scale[i]);
					inla_error_general(msg);
				}
			}
		}
	}
		break;

	case L_0POISSON:
	case L_0POISSONS:
	{
		for (i = 0; i < mb->predictor_ndata; i++) {
			if (ds->data_observations.d[i]) {
				if (ds->data_observations.poisson0_E[i] <= 0.0 || ds->data_observations.y[i] < 0.0) {
					GMRFLib_sprintf(&msg, "%s: 0poisson(S) E[%1d] = %g y[%1d] = %g is void\n", secname, i,
							ds->data_observations.poisson0_E[i], ds->data_observations.y[i]);
					inla_error_general(msg);
				}
			}
		}
	}
		break;

	case L_0BINOMIAL:
	case L_0BINOMIALS:
	{
		for (i = 0; i < mb->predictor_ndata; i++) {
			if (ds->data_observations.d[i]) {
				if (ds->data_observations.binomial0_Ntrials[i] < ds->data_observations.y[i] || ds->data_observations.y[i] < 0.0) {
					GMRFLib_sprintf(&msg, "%s: 0binomial(S) Ntrials[%1d] = %g y[%1d] = %g is void\n", secname, i,
							ds->data_observations.binomial0_Ntrials[i], ds->data_observations.y[i]);
					inla_error_general(msg);
				}
			}
		}
	}
		break;

	case L_BINOMIALMIX:
		// later
		break;

	case L_AGAUSSIAN:
	{
		for (i = 0; i < mb->predictor_ndata; i++) {
			if (ds->data_observations.d[i]) {
				if (ds->data_observations.agaussian[0][i] < 0.0 ||
				    ds->data_observations.agaussian[2][i] <= 0.0 || ds->data_observations.agaussian[3][i] <= 0.0) {
					GMRFLib_sprintf(&msg, "%s: AGAUSSIAN data[%1d] = %g %g %g %g %g is void\n", secname, i,
							ds->data_observations.agaussian[0][i],
							ds->data_observations.agaussian[1][i],
							ds->data_observations.agaussian[2][i], ds->data_observations.agaussian[3][i]);
					inla_error_general(msg);
				}
			}
		}
	}
		break;

	case L_T:
	{
		for (i = 0; i < mb->predictor_ndata; i++) {
			if (ds->data_observations.d[i]) {
				if (ds->data_observations.weight_t[i] <= 0.0) {
					GMRFLib_sprintf(&msg, "%s: Student-t weight[%1d] = %g is void\n", secname, i,
							ds->data_observations.weight_t[i]);
					inla_error_general(msg);
				}
			}
		}
	}
		break;

	case L_TSTRATA:
	{
		for (i = 0; i < mb->predictor_ndata; i++) {
			if (ds->data_observations.d[i]) {
				if (ds->data_observations.weight_tstrata[i] <= 0.0) {
					GMRFLib_sprintf(&msg, "%s: t weight[%1d] = %g is void\n", secname, i,
							ds->data_observations.weight_tstrata[i]);
					inla_error_general(msg);
				}
				if ((int) (ds->data_observations.strata_tstrata[i]) < 0) {
					GMRFLib_sprintf(&msg, "%s: tstrata strata[%1d] = %g is void\n", secname, i,
							ds->data_observations.strata_tstrata[i]);
					inla_error_general(msg);
				}
			}
		}
	}
		break;

	case L_BELL:
	{
		// to initialize the lbell-cache
		int ymax = 0;

		for (i = 0; i < mb->predictor_ndata; i++) {
			if (ds->data_observations.d[i]) {
				if (ds->data_observations.E[i] <= 0.0 || ds->data_observations.y[i] < 0.0) {
					GMRFLib_sprintf(&msg, "%s: Bell data[%1d] (e,y) = (%g,%g) is void\n", secname, i,
							ds->data_observations.E[i], ds->data_observations.y[i]);
					inla_error_general(msg);
				}
				ymax = IMAX(ymax, (int) ds->data_observations.y[i]);
			}
		}
		// will initialize the cache, not strictly needed but its convenient to do this here 
		my_lbell(ymax);
	}
		break;

	case L_POISSON:
	case L_NPOISSON:
	case L_XPOISSON:
	{
		for (i = 0; i < mb->predictor_ndata; i++) {
			if (ds->data_observations.d[i]) {
				if (ds->data_observations.E[i] < 0.0 || ds->data_observations.y[i] < 0.0) {
					GMRFLib_sprintf(&msg, "%s: Poisson data[%1d] (e,y) = (%g,%g) is void\n", secname, i,
							ds->data_observations.E[i], ds->data_observations.y[i]);
					inla_error_general(msg);
				}
			}
		}
	}
		break;

	case L_NZPOISSON:
	{
		for (i = 0; i < mb->predictor_ndata; i++) {
			if (ds->data_observations.d[i]) {
				if (ds->data_observations.E[i] <= 0.0 || ds->data_observations.y[i] < 1.0) {
					GMRFLib_sprintf(&msg, "%s: nzPoisson data[%1d] (e,y) = (%g,%g) is void\n", secname, i,
							ds->data_observations.E[i], ds->data_observations.y[i]);
					inla_error_general(msg);
				}
			}
		}
	}
		break;

	case L_CONTPOISSON:
	case L_QCONTPOISSON:
	{
		for (i = 0; i < mb->predictor_ndata; i++) {
			if (ds->data_observations.d[i]) {
				if (ds->data_observations.E[i] < 0.0 || ds->data_observations.y[i] <= 0.0) {
					GMRFLib_sprintf(&msg, "%s: ContPoisson data[%1d] (e,y) = (%g,%g) is void\n", secname, i,
							ds->data_observations.E[i], ds->data_observations.y[i]);
					inla_error_general(msg);
				}
			}
		}
	}
		break;

	case L_CENPOISSON:
	{
		for (i = 0; i < mb->predictor_ndata; i++) {
			if (ds->data_observations.d[i]) {
				if (ds->data_observations.E[i] < 0.0 || ds->data_observations.y[i] < 0.0) {
					GMRFLib_sprintf(&msg, "%s: CenPoisson data[%1d] (e,y) = (%g,%g) is void\n", secname, i,
							ds->data_observations.E[i], ds->data_observations.y[i]);
					inla_error_general(msg);
				}
			}
		}
	}
		break;

	case L_CENPOISSON2:
	{
		for (i = 0; i < mb->predictor_ndata; i++) {
			if (ds->data_observations.d[i]) {
				if (ds->data_observations.E[i] < 0.0 || ds->data_observations.y[i] < 0.0) {
					GMRFLib_sprintf(&msg, "%s: CenPoisson2 data[%1d] (e,y) = (%g,%g) is void\n", secname, i,
							ds->data_observations.E[i], ds->data_observations.y[i]);
					inla_error_general(msg);
				}
				if ((ds->data_observations.cen_high[i] > 0 &&
				     (ds->data_observations.cen_high[i] < ds->data_observations.cen_low[i]))) {
					GMRFLib_sprintf(&msg, "%s: CenPoisson2 (idx,low,high) = (%d,%g,%g) is void\n", secname, i,
							ds->data_observations.cen_low[i], ds->data_observations.cen_high[i]);
					inla_error_general(msg);
				}
			}
		}
	}
		break;

	case L_CENNBINOMIAL2:
	{
		for (i = 0; i < mb->predictor_ndata; i++) {
			if (ds->data_observations.d[i]) {
				if (ds->data_observations.E[i] < 0.0 || ds->data_observations.y[i] < 0.0) {
					GMRFLib_sprintf(&msg, "%s: CenNBinomial2 data[%1d] (e,y) = (%g,%g) is void\n", secname, i,
							ds->data_observations.E[i], ds->data_observations.y[i]);
					inla_error_general(msg);
				}
				if ((ds->data_observations.cen_high[i] > 0 &&
				     (ds->data_observations.cen_high[i] < ds->data_observations.cen_low[i]))) {
					GMRFLib_sprintf(&msg, "%s: CenNBinomial2 (idx,low,high) = (%d,%g,%g) is void\n", secname, i,
							ds->data_observations.cen_low[i], ds->data_observations.cen_high[i]);
					inla_error_general(msg);
				}
				if (ds->data_observations.S[i] <= 0.0) {
					GMRFLib_sprintf(&msg, "%s: CenNBinomial2 data[%1d] S = %g is void\n", secname, i,
							ds->data_observations.S[i]);
					inla_error_general(msg);
				}
			}
		}
	}
		break;

	case L_GAUSSIANJW:
	{
		for (i = 0; i < mb->predictor_ndata; i++) {
			if (ds->data_observations.d[i]) {
				if (ds->data_observations.gjw_n[i] <= 0.0 || ds->data_observations.gjw_df[i] <= 0.0 ||
				    ds->data_observations.gjw_var[i] <= 0.0) {
					GMRFLib_sprintf(&msg, "%s: GaussianJW data[%1d] (n,df,var) = (%g,%g,%g) is void\n", secname, i,
							ds->data_observations.gjw_n[i],
							ds->data_observations.gjw_df[i], ds->data_observations.gjw_var[i]);
					inla_error_general(msg);
				}
			}
		}
	}
		break;

	case L_POISSON_SPECIAL1:
	{
		for (i = 0; i < mb->predictor_ndata; i++) {
			if (ds->data_observations.d[i]) {
				if (ds->data_observations.E[i] < 0.0 || ds->data_observations.y[i] < 1.0) {
					GMRFLib_sprintf(&msg, "%s: Poisson.special1 data[%1d] (e,y) = (%g,%g) is void\n", secname, i,
							ds->data_observations.E[i], ds->data_observations.y[i]);
					inla_error_general(msg);
				}
			}
		}
	}
		break;

	case L_POM:
	{
		int nclasses = -1, iy;
		for (i = 0; i < mb->predictor_ndata; i++) {
			if (ds->data_observations.d[i]) {
				iy = (int) ds->data_observations.y[i];
				nclasses = IMAX(nclasses, iy);
				if (iy <= 0) {
					GMRFLib_sprintf(&msg, "%s: POM data[%1d] (y) = %g is void\n", secname, i, ds->data_observations.y[i]);
					inla_error_general(msg);
				}
				if (iy != (int) ds->data_observations.y[i]) {
					GMRFLib_sprintf(&msg, "%s: POM data[%1d] (y) = (%g) is not an integer\n", secname, i,
							ds->data_observations.y[i]);
					inla_error_general(msg);
				}
			}
		}
		GMRFLib_ASSERT(nclasses <= POM_MAXTHETA, GMRFLib_EPARAMETER);
		GMRFLib_ASSERT(nclasses > 1, GMRFLib_EPARAMETER);
		ds->data_observations.pom_nclasses = nclasses;
		assert(nclasses > 0);

		int *check = Calloc(nclasses + 1, int);
		for (i = 0; i < mb->predictor_ndata; i++) {
			if (ds->data_observations.d[i]) {
				iy = (int) ds->data_observations.y[i];
				check[iy] = 1;
			}
		}
		for (int k = 1; k <= nclasses; k++) {
			if (check[k] == 0) {
				GMRFLib_sprintf(&msg, "%s: POM: There are no observations for class [%1d]. Not allowed\n", secname, k);
				inla_error_general(msg);
			}
		}
		Free(check);
	}
		break;

	case L_GAMMACOUNT:
	{
		for (i = 0; i < mb->predictor_ndata; i++) {
			if (ds->data_observations.d[i]) {
				if (ds->data_observations.E[i] < 0.0 || ds->data_observations.y[i] < 0.0) {
					GMRFLib_sprintf(&msg, "%s: Gammacount data[%1d] (E,y) = (%g, %g) is void\n", secname, i,
							ds->data_observations.E[i], ds->data_observations.y[i]);
					inla_error_general(msg);
				}
			}
		}
	}
		break;

	case L_EXPONENTIAL:
	{
		for (i = 0; i < mb->predictor_ndata; i++) {
			if (ds->data_observations.d[i]) {
				if (ds->data_observations.y[i] < 0.0) {
					GMRFLib_sprintf(&msg, "%s: Exponential data[%1d] (y) = (%g) is void\n", secname, i,
							ds->data_observations.y[i]);
					inla_error_general(msg);
				}
			}
		}
	}
		break;

	case L_WEIBULL:
	{
		for (i = 0; i < mb->predictor_ndata; i++) {
			if (ds->data_observations.d[i]) {
				if (ds->data_observations.y[i] <= 0.0) {
					GMRFLib_sprintf(&msg, "%s: Weibull data[%1d] (y) = (%g) is void\n", secname, i, ds->data_observations.y[i]);
					inla_error_general(msg);
				}
			}
		}
	}
		break;

	case L_GOMPERTZ:
	{
		for (i = 0; i < mb->predictor_ndata; i++) {
			if (ds->data_observations.d[i]) {
				if (ds->data_observations.y[i] < 0.0) {
					GMRFLib_sprintf(&msg, "%s: GOMPERTZ data[%1d] (y) = (%g) is void\n", secname, i,
							ds->data_observations.y[i]);
					inla_error_general(msg);
				}
			}
		}
	}
		break;

	case L_GPOISSON:
	case L_ZEROINFLATEDPOISSON0:
	case L_ZEROINFLATEDPOISSON1:
	case L_ZEROINFLATEDPOISSON2:
	case L_ZEROINFLATEDCENPOISSON0:
	case L_ZEROINFLATEDCENPOISSON1:
	case L_ZEROINFLATEDNBINOMIAL0:
	case L_ZEROINFLATEDNBINOMIAL1:
	case L_ZEROINFLATEDNBINOMIAL2:
	{
		for (i = 0; i < mb->predictor_ndata; i++) {
			if (ds->data_observations.d[i]) {
				if (ds->data_observations.E[i] <= 0.0 || ds->data_observations.y[i] < 0.0) {
					GMRFLib_sprintf(&msg, "%s: Poisson-like data[%1d] (e,y) = (%g,%g) is void\n", secname, i,
							ds->data_observations.E[i], ds->data_observations.y[i]);
					inla_error_general(msg);
				}
			}
		}
	}
		break;

	case L_NBINOMIAL:
	{
		for (i = 0; i < mb->predictor_ndata; i++) {
			if (ds->data_observations.d[i]) {
				if (ds->data_observations.E[i] <= 0.0 || ds->data_observations.y[i] < 0.0) {
					GMRFLib_sprintf(&msg, "%s: Poisson-like data[%1d] (E,y) = (%g,%g) is void\n", secname, i,
							ds->data_observations.E[i], ds->data_observations.y[i]);
					inla_error_general(msg);
				}
				if (ds->data_observations.S[i] <= 0.0) {
					GMRFLib_sprintf(&msg, "%s: Poisson-like data[%1d] S = %g is void\n", secname, i,
							ds->data_observations.S[i]);
					inla_error_general(msg);
				}
			}
		}
	}
		break;

	case L_ZEROINFLATEDNBINOMIAL1STRATA2:
	case L_ZEROINFLATEDNBINOMIAL1STRATA3:
	{
		for (i = 0; i < mb->predictor_ndata; i++) {
			if (ds->data_observations.d[i]) {
				if (ds->data_observations.E[i] < 0.0 || ds->data_observations.y[i] < 0.0) {
					GMRFLib_sprintf(&msg, "%s: Poisson data[%1d] (e,y) = (%g,%g) is void\n", secname, i,
							ds->data_observations.E[i], ds->data_observations.y[i]);
					inla_error_general(msg);
				}
			}
		}
	}
		break;

	case L_GAMMA:
	case L_MGAMMA:
	{
		for (i = 0; i < mb->predictor_ndata; i++) {
			if (ds->data_observations.d[i]) {
				if (ds->data_observations.y[i] <= 0.0 || ds->data_observations.gamma_scale[i] <= 0.0) {
					GMRFLib_sprintf(&msg, "%s: (m)Gamma data[%1d] (y) = %g or weight %g is void\n", secname, i,
							ds->data_observations.y[i], ds->data_observations.gamma_scale[i]);
					inla_error_general(msg);
				}
			}
		}
	}
		break;

	case L_GAMMAJW:
	{
		for (i = 0; i < mb->predictor_ndata; i++) {
			if (ds->data_observations.d[i]) {
				if (ds->data_observations.y[i] <= 0.0) {
					GMRFLib_sprintf(&msg, "%s: Gammajw data[%1d] (y) = %g or weight %g is void\n", secname, i,
							ds->data_observations.y[i], ds->data_observations.gamma_scale[i]);
					inla_error_general(msg);
				}
			}
		}
	}
		break;

	case L_QKUMAR:
	{
		for (i = 0; i < mb->predictor_ndata; i++) {
			if (ds->data_observations.d[i]) {
				if (ds->data_observations.y[i] <= 0.0 || ds->data_observations.y[i] >= 1.0) {
					GMRFLib_sprintf(&msg, "%s: qKumar data[%1d] (y) = (%g) is void\n", secname, i, ds->data_observations.y[i]);
					inla_error_general(msg);
				}
			}
		}
	}
		break;

	case L_LOGLOGISTIC:
	case L_QLOGLOGISTIC:
	{
		for (i = 0; i < mb->predictor_ndata; i++) {
			if (ds->data_observations.d[i]) {
				if (ds->data_observations.y[i] <= 0.0) {
					GMRFLib_sprintf(&msg, "%s: log-logistic/normal data[%1d] (y) = (%g) is void\n", secname, i,
							ds->data_observations.y[i]);
					inla_error_general(msg);
				}
			}
		}
	}
		break;

	case L_BETA:
	{
		for (i = 0; i < mb->predictor_ndata; i++) {
			if (ds->data_observations.d[i]) {
				if (ds->data_observations.y[i] == 0.0 || ds->data_observations.y[i] == 1.0) {
					beta_delayed_error++;
				} else if (ds->data_observations.y[i] < 0.0 || ds->data_observations.y[i] > 1.0
					   || ds->data_observations.beta_weight[i] <= 0.0) {
					GMRFLib_sprintf(&msg, "%s: Beta data[%1d] (y) = (%g) or weight (%g)is void\n", secname, i,
							ds->data_observations.y[i], ds->data_observations.beta_weight[i]);
					inla_error_general(msg);
				}
			}
		}
	}
		break;

	case L_ZERO_N_INFLATEDBINOMIAL2:
	case L_ZERO_N_INFLATEDBINOMIAL3:
	{
		for (i = 0; i < mb->predictor_ndata; i++) {
			if (ds->data_observations.d[i]) {
				if (ds->data_observations.nb[i] <= 0.0 ||
				    ds->data_observations.y[i] > ds->data_observations.nb[i] || ds->data_observations.y[i] < 0.0) {
					// we allow for binomial(0,p) if y = 0
					if (!(ISZERO(ds->data_observations.nb[i]) && ISZERO(ds->data_observations.y[i]))) {
						GMRFLib_sprintf(&msg, "%s: Binomial data[%1d] (nb,y) = (%g,%g) is void\n", secname,
								i, ds->data_observations.nb[i], ds->data_observations.y[i]);
						inla_error_general(msg);
					}
				}
			}
		}
	}
		break;

	case L_XBINOMIAL:
	{
		for (i = 0; i < mb->predictor_ndata; i++) {
			if (ds->data_observations.d[i]) {
				if (ds->data_observations.nb[i] <= 0.0 ||
				    ds->data_observations.y[i] > ds->data_observations.nb[i] || ds->data_observations.y[i] < 0.0 ||
				    ds->data_observations.p_scale[i] <= 0.0 || ds->data_observations.p_scale[i] > 1.0) {
					GMRFLib_sprintf(&msg, "%s: xBinomial data[%1d] (nb,p.scale,y) = (%g,%g,%g) is void\n", secname,
							i, ds->data_observations.nb[i], ds->data_observations.p_scale[i],
							ds->data_observations.y[i]);
					inla_error_general(msg);
				}
			}
		}
	}
		break;

	case L_BINOMIAL:
	{
		if (ds->variant == 0) {
			for (i = 0; i < mb->predictor_ndata; i++) {
				if (ds->data_observations.d[i]) {
					if (ds->data_observations.nb[i] <= 0.0 ||
					    ds->data_observations.y[i] > ds->data_observations.nb[i] || ds->data_observations.y[i] < 0.0) {
						// we allow for binomial(0,p) if y = 0
						if (!(ISZERO(ds->data_observations.nb[i]) && ISZERO(ds->data_observations.y[i]))) {
							GMRFLib_sprintf(&msg, "%s: Binomial data[%1d] (nb,y) = (%g,%g) is void\n", secname,
									i, ds->data_observations.nb[i], ds->data_observations.y[i]);
							inla_error_general(msg);
						}
					}
				}
			}
		} else {
			// neg binomial
			for (i = 0; i < mb->predictor_ndata; i++) {
				if (ds->data_observations.d[i]) {
					if (!((ds->data_observations.nb[i] - ds->data_observations.y[i]) >= 0.0 &&
					      ds->data_observations.y[i] >= 1.0)) {
						GMRFLib_sprintf(&msg, "%s: Binomial data[%1d] variant=%1d, (nb,y) = (%g,%g) is void\n", secname,
								i, ds->variant, ds->data_observations.nb[i], ds->data_observations.y[i]);
						inla_error_general(msg);
					}
				}
			}
		}
	}
		break;

	case L_NBINOMIAL2:
	{
		for (i = 0; i < mb->predictor_ndata; i++) {
			if (ds->data_observations.d[i]) {
				if (ds->data_observations.nb[i] <= 0.0 || ds->data_observations.y[i] < 0.0) {
					GMRFLib_sprintf(&msg, "%s: NBinomial2 data[%1d] (nb,y) = (%g,%g) is void\n", secname,
							i, ds->data_observations.nb[i], ds->data_observations.y[i]);
					inla_error_general(msg);
				}
			}
		}
	}
		break;

	case L_NMIX:
	case L_NMIXNB:
	{
		for (i = 0; i < mb->predictor_ndata; i++) {
			if (ds->data_observations.d[i]) {
				for (j = 0; j > -1; j++) {
					// printf("%d %d %g\n", i, j, ds->data_observations.nmix_y[j][i]);
					if (gsl_isnan(ds->data_observations.nmix_y[j][i]))
						break;
					if (ds->data_observations.nmix_y[j][i] < 0) {
						GMRFLib_sprintf(&msg, "%s: NMix data[%1d][%1d] (y) = (%g) is void\n", secname,
								i, j, ds->data_observations.nmix_y[j][i]);
						inla_error_general(msg);
					}
				}
				assert(ds->data_observations.y[i] < 0);	/* have to be void */
			}
		}
	}
		break;

	case L_ZEROINFLATEDBINOMIAL0:
	case L_ZEROINFLATEDBINOMIAL1:
	case L_ZEROINFLATEDBINOMIAL2:
	case L_ZEROINFLATEDBETABINOMIAL0:
	case L_ZEROINFLATEDBETABINOMIAL1:
	case L_ZEROINFLATEDBETABINOMIAL2:
	case L_BETABINOMIAL:
	{
		for (i = 0; i < mb->predictor_ndata; i++) {
			if (ds->data_observations.d[i]) {
				if (ds->data_observations.nb[i] <= 0.0 ||
				    ds->data_observations.y[i] > ds->data_observations.nb[i] || ds->data_observations.y[i] < 0.0) {
					GMRFLib_sprintf(&msg, "%s: Binomial data[%1d] (nb,y) = (%g,%g) is void\n", secname,
							i, ds->data_observations.nb[i], ds->data_observations.y[i]);
					inla_error_general(msg);
				}
			}
		}
	}
		break;

	case L_BETABINOMIALNA:
	{
		// Since we're using a normal approximation, the data can be negative and also larger then nb
		for (i = 0; i < mb->predictor_ndata; i++) {
			if (ds->data_observations.d[i]) {
				if (ds->data_observations.nb[i] <= 0.0 || ds->data_observations.betabinomialnb_scale[i] < 0.0) {
					GMRFLib_sprintf(&msg, "%s: BetaBinomialNA data[%1d] (nb,scale,y) = (%g,%g,%g) is void\n", secname,
							i, ds->data_observations.nb[i],
							ds->data_observations.betabinomialnb_scale[i], ds->data_observations.y[i]);
					inla_error_general(msg);
				}
			}
		}
	}
		break;

	case L_CBINOMIAL:
	{
		for (i = 0; i < mb->predictor_ndata; i++) {
			if (ds->data_observations.d[i]) {
				if (ds->data_observations.cbinomial_n[i] <= 0.0 ||
				    ds->data_observations.cbinomial_k[i] <= 0.0 || ds->data_observations.y[i] > ds->data_observations.cbinomial_k[i]
				    || ds->data_observations.y[i] < 0.0) {
					GMRFLib_sprintf(&msg, "%s: CBinomial data[%1d] (k,n,y) = (%g,%g,%g) is void\n", secname, i,
							ds->data_observations.cbinomial_k[i], ds->data_observations.cbinomial_n[i],
							ds->data_observations.y[i]);
					inla_error_general(msg);
				}
			}
		}
	}
		break;


	case L_EXPONENTIALSURV:
	case L_GAMMASURV:
	case L_MGAMMASURV:
	case L_GAMMAJWSURV:
	case L_WEIBULLSURV:
	case L_LOGLOGISTICSURV:
	case L_QLOGLOGISTICSURV:
	case L_LOGNORMALSURV:
	case L_FMRISURV:
	case L_GOMPERTZSURV:
	{
		switch (ds->data_id) {
		case L_WEIBULLSURV:
		case L_GAMMASURV:
		case L_MGAMMASURV:
		case L_LOGNORMAL:
		{
			// those who cannot take y=0

			for (i = 0; i < mb->predictor_ndata; i++) {
				if (ds->data_observations.d[i] && ((int) ds->data_observations.event[i] == SURV_EVENT_FAILURE ||
								   (int) ds->data_observations.event[i] == SURV_EVENT_ININTERVAL)) {
					if (ds->data_observations.y[i] <= 0.0) {
						GMRFLib_sprintf(&msg, "%s: Weibull/Gamma/logNormal data[%1d] (y) = (%.12g) is void\n",
								secname, i, ds->data_observations.y[i]);
						inla_error_general(msg);
					}
				}
			}
			break;
		}
		default:
			break;
		}

		for (i = 0; i < mb->predictor_ndata; i++) {
			if (ds->data_observations.d[i]) {
				int event;
				double truncation, lower, upper, ttime;

				truncation = ds->data_observations.truncation[i];
				lower = ds->data_observations.lower[i];
				upper = ds->data_observations.upper[i];
				ttime = ds->data_observations.y[i];
				event = (int) (ds->data_observations.event[i]);

#define _SERR { GMRFLib_sprintf(&msg, "%s: survival data[%1d] (event,trunc,lower,upper,y) = (%g,%g,%g,%g,%g) is void\n", \
				secname, i, (double)event, truncation, lower,  upper,  ttime); inla_error_general(msg); }

				if (truncation < 0.0 || lower < 0.0 || upper < 0.0 || ttime < 0.0)
					_SERR;

				switch (event) {
				case SURV_EVENT_FAILURE:
				{
					if (ttime < truncation)
						_SERR;
				}
					break;
				case SURV_EVENT_RIGHT:
					if (lower < truncation)
						_SERR;
					break;
				case SURV_EVENT_LEFT:
				{
					if (upper < truncation)
						_SERR;
				}
					break;
				case SURV_EVENT_INTERVAL:
				{
					if (DMIN(lower, upper) < truncation || upper < lower)
						_SERR;
				}
					break;
				case SURV_EVENT_ININTERVAL:
				{
					if (DMIN(lower, upper) < truncation || upper < lower || ttime < lower || ttime > upper)
						_SERR;
				}
					break;
				default:
					_SERR;
					break;
				}
#undef _SERR
			}
		}
	}
		break;

	case L_FL:
		break;

	default:
		break;
	}

	/*
	 * common for all 
	 */
	ds->variant = iniparser_getint(ini, inla_string_join(secname, "VARIANT"), 0);
	if (mb->verbose) {
		printf("\t\tlikelihood.variant=[%1d]\n", ds->variant);
	}
	ds->data_observations.quantile = iniparser_getdouble(ini, inla_string_join(secname, "QUANTILE"), -1.0);
	if (mb->verbose) {
		if (ds->data_observations.quantile > 0.0 && ds->data_observations.quantile < 1.0) {
			printf("\t\tlikelihood.quantile=[%g]\n", ds->data_observations.quantile);
		}
	}

	/*
	 * read spesific options and define hyperparameters, if any.
	 */
	switch (ds->data_id) {
	case L_GAUSSIAN:
	{
		/*
		 * get options related to the gaussian 
		 */
		tmp = iniparser_getdouble(ini, inla_string_join(secname, "INITIAL0"), G.log_prec_initial);
		ds->data_fixed0 = iniparser_getboolean(ini, inla_string_join(secname, "FIXED0"), 0);
		if (!ds->data_fixed0 && mb->mode_use_mode) {
			tmp = mb->theta_file[mb->theta_counter_file++];
			if (mb->mode_fixed)
				ds->data_fixed0 = 1;
		}
		HYPER_NEW(ds->data_observations.log_prec_gaussian, tmp);
		if (mb->verbose) {
			printf("\t\tinitialise log_precision[%g]\n", ds->data_observations.log_prec_gaussian[0][0]);
			printf("\t\tfixed0=[%1d]\n", ds->data_fixed0);
		}
		inla_read_prior0(mb, ini, sec, &(ds->data_prior0), "LOGGAMMA", NULL);

		/*
		 * add theta 
		 */
		if (!ds->data_fixed0) {
			mb->theta = Realloc(mb->theta, mb->ntheta + 1, double **);
			mb->theta_hyperid = Realloc(mb->theta_hyperid, mb->ntheta + 1, char *);
			mb->theta_hyperid[mb->ntheta] = ds->data_prior0.hyperid;
			mb->theta_tag = Realloc(mb->theta_tag, mb->ntheta + 1, char *);
			mb->theta_tag_userscale = Realloc(mb->theta_tag_userscale, mb->ntheta + 1, char *);
			mb->theta_dir = Realloc(mb->theta_dir, mb->ntheta + 1, char *);
			mb->theta_tag[mb->ntheta] = inla_make_tag("Log precision for the Gaussian observations", mb->ds);
			mb->theta_tag_userscale[mb->ntheta] = inla_make_tag("Precision for the Gaussian observations", mb->ds);
			GMRFLib_sprintf(&msg, "%s-parameter0", secname);
			mb->theta_dir[mb->ntheta] = msg;

			mb->theta_from = Realloc(mb->theta_from, mb->ntheta + 1, char *);
			mb->theta_to = Realloc(mb->theta_to, mb->ntheta + 1, char *);
			mb->theta_from[mb->ntheta] = Strdup(ds->data_prior0.from_theta);
			mb->theta_to[mb->ntheta] = Strdup(ds->data_prior0.to_theta);

			mb->theta[mb->ntheta] = ds->data_observations.log_prec_gaussian;
			mb->theta_map = Realloc(mb->theta_map, mb->ntheta + 1, map_func_tp *);
			mb->theta_map[mb->ntheta] = map_precision;
			mb->theta_map_arg = Realloc(mb->theta_map_arg, mb->ntheta + 1, void *);
			mb->theta_map_arg[mb->ntheta] = NULL;
			mb->ntheta++;
			ds->data_ntheta++;
		}
		// this one must be fixed.
		tmp = iniparser_getdouble(ini, inla_string_join(secname, "INITIAL1"), 999999.9);
		ds->data_fixed1 = iniparser_getboolean(ini, inla_string_join(secname, "FIXED1"), 1);
		HYPER_NEW(ds->data_observations.log_prec_gaussian_offset, tmp);
		if (mb->verbose) {
			printf("\t\tinitialise log_precision offset[%g]\n", ds->data_observations.log_prec_gaussian_offset[0][0]);
			printf("\t\tfixed1=[%1d]\n", ds->data_fixed1);
		}
		inla_read_prior1(mb, ini, sec, &(ds->data_prior1), "NONE", NULL);
		if (!ds->data_fixed1) {
			inla_error_field_is_void(__GMRFLib_FuncName, secname, "FIXED1", "1");
		}
	}
		break;

	case L_AGAUSSIAN:
	{
		/*
		 * get options related to the agaussian 
		 */
		tmp = iniparser_getdouble(ini, inla_string_join(secname, "INITIAL"), G.log_prec_initial);
		ds->data_fixed = iniparser_getboolean(ini, inla_string_join(secname, "FIXED"), 0);
		if (!ds->data_fixed && mb->mode_use_mode) {
			tmp = mb->theta_file[mb->theta_counter_file++];
			if (mb->mode_fixed)
				ds->data_fixed = 1;
		}
		HYPER_NEW(ds->data_observations.log_prec_gaussian, tmp);
		if (mb->verbose) {
			printf("\t\tinitialise log_precision[%g]\n", ds->data_observations.log_prec_gaussian[0][0]);
			printf("\t\tfixed=[%1d]\n", ds->data_fixed);
		}
		inla_read_prior(mb, ini, sec, &(ds->data_prior), "LOGGAMMA", NULL);

		/*
		 * add theta 
		 */
		if (!ds->data_fixed) {
			mb->theta = Realloc(mb->theta, mb->ntheta + 1, double **);
			mb->theta_hyperid = Realloc(mb->theta_hyperid, mb->ntheta + 1, char *);
			mb->theta_hyperid[mb->ntheta] = ds->data_prior.hyperid;
			mb->theta_tag = Realloc(mb->theta_tag, mb->ntheta + 1, char *);
			mb->theta_tag_userscale = Realloc(mb->theta_tag_userscale, mb->ntheta + 1, char *);
			mb->theta_dir = Realloc(mb->theta_dir, mb->ntheta + 1, char *);
			mb->theta_tag[mb->ntheta] = inla_make_tag("Log precision for the AggGaussian observations", mb->ds);
			mb->theta_tag_userscale[mb->ntheta] = inla_make_tag("Precision for the AggGaussian observations", mb->ds);
			GMRFLib_sprintf(&msg, "%s-parameter", secname);
			mb->theta_dir[mb->ntheta] = msg;

			mb->theta_from = Realloc(mb->theta_from, mb->ntheta + 1, char *);
			mb->theta_to = Realloc(mb->theta_to, mb->ntheta + 1, char *);
			mb->theta_from[mb->ntheta] = Strdup(ds->data_prior.from_theta);
			mb->theta_to[mb->ntheta] = Strdup(ds->data_prior.to_theta);

			mb->theta[mb->ntheta] = ds->data_observations.log_prec_gaussian;
			mb->theta_map = Realloc(mb->theta_map, mb->ntheta + 1, map_func_tp *);
			mb->theta_map[mb->ntheta] = map_precision;
			mb->theta_map_arg = Realloc(mb->theta_map_arg, mb->ntheta + 1, void *);
			mb->theta_map_arg[mb->ntheta] = NULL;
			mb->ntheta++;
			ds->data_ntheta++;
		}
	}
		break;

	case L_LOGNORMAL:
	{
		/*
		 * get options related to the lognormal
		 */
		tmp = iniparser_getdouble(ini, inla_string_join(secname, "INITIAL"), G.log_prec_initial);
		ds->data_fixed = iniparser_getboolean(ini, inla_string_join(secname, "FIXED"), 0);
		if (!ds->data_fixed && mb->mode_use_mode) {
			tmp = mb->theta_file[mb->theta_counter_file++];
			if (mb->mode_fixed)
				ds->data_fixed = 1;
		}
		HYPER_NEW(ds->data_observations.log_prec_gaussian, tmp);
		if (mb->verbose) {
			printf("\t\tinitialise log_precision[%g]\n", ds->data_observations.log_prec_gaussian[0][0]);
			printf("\t\tfixed=[%1d]\n", ds->data_fixed);
		}
		inla_read_prior(mb, ini, sec, &(ds->data_prior), "LOGGAMMA", NULL);

		/*
		 * add theta 
		 */
		if (!ds->data_fixed) {
			mb->theta = Realloc(mb->theta, mb->ntheta + 1, double **);
			mb->theta_hyperid = Realloc(mb->theta_hyperid, mb->ntheta + 1, char *);
			mb->theta_hyperid[mb->ntheta] = ds->data_prior.hyperid;
			mb->theta_tag = Realloc(mb->theta_tag, mb->ntheta + 1, char *);
			mb->theta_tag_userscale = Realloc(mb->theta_tag_userscale, mb->ntheta + 1, char *);
			mb->theta_dir = Realloc(mb->theta_dir, mb->ntheta + 1, char *);

			if (ds->data_id == L_LOGNORMAL) {
				mb->theta_tag[mb->ntheta] = inla_make_tag("Log precision for the lognormal observations", mb->ds);
				mb->theta_tag_userscale[mb->ntheta] = inla_make_tag("Precision for the lognormal observations", mb->ds);
			} else {
				mb->theta_tag[mb->ntheta] = inla_make_tag("Log precision for the lognormalsurv observations", mb->ds);
				mb->theta_tag_userscale[mb->ntheta] = inla_make_tag("Precision for the lognormalsurv observations", mb->ds);
			}
			GMRFLib_sprintf(&msg, "%s-parameter", secname);
			mb->theta_dir[mb->ntheta] = msg;

			mb->theta_from = Realloc(mb->theta_from, mb->ntheta + 1, char *);
			mb->theta_to = Realloc(mb->theta_to, mb->ntheta + 1, char *);
			mb->theta_from[mb->ntheta] = Strdup(ds->data_prior.from_theta);
			mb->theta_to[mb->ntheta] = Strdup(ds->data_prior.to_theta);

			mb->theta[mb->ntheta] = ds->data_observations.log_prec_gaussian;
			mb->theta_map = Realloc(mb->theta_map, mb->ntheta + 1, map_func_tp *);
			mb->theta_map[mb->ntheta] = map_precision;
			mb->theta_map_arg = Realloc(mb->theta_map_arg, mb->ntheta + 1, void *);
			mb->theta_map_arg[mb->ntheta] = NULL;
			mb->ntheta++;
			ds->data_ntheta++;
		}
	}
		break;

	case L_LOGNORMALSURV:
	{
		for (i = 0; i < CURE_MAXTHETA + 1; i++) {
			GMRFLib_sprintf(&ctmp, "FIXED%1d", i);
			iniparser_getstring(ini, inla_string_join(secname, ctmp), NULL);

			GMRFLib_sprintf(&ctmp, "INITIAL%1d", i);
			iniparser_getstring(ini, inla_string_join(secname, ctmp), NULL);

			GMRFLib_sprintf(&ctmp, "PRIOR%1d", i);
			iniparser_getstring(ini, inla_string_join(secname, ctmp), NULL);

			GMRFLib_sprintf(&ctmp, "HYPERID%1d", i);
			iniparser_getstring(ini, inla_string_join(secname, ctmp), NULL);

			GMRFLib_sprintf(&ctmp, "PARAMETERS%1d", i);
			iniparser_getstring(ini, inla_string_join(secname, ctmp), NULL);

			GMRFLib_sprintf(&ctmp, "to.theta%1d", i);
			iniparser_getstring(ini, inla_string_join(secname, ctmp), NULL);

			GMRFLib_sprintf(&ctmp, "from.theta%1d", i);
			iniparser_getstring(ini, inla_string_join(secname, ctmp), NULL);
		}

		ds->data_nfixed = Calloc(CURE_MAXTHETA + 1, int);
		ds->data_nprior = Calloc(CURE_MAXTHETA + 1, Prior_tp);
		ds->data_observations.cure_beta = Calloc(CURE_MAXTHETA, double **);

		tmp = iniparser_getdouble(ini, inla_string_join(secname, "INITIAL0"), G.log_prec_initial);
		ds->data_nfixed[0] = iniparser_getboolean(ini, inla_string_join(secname, "FIXED0"), 0);
		if (!ds->data_nfixed[0] && mb->mode_use_mode) {
			tmp = mb->theta_file[mb->theta_counter_file++];
			if (mb->mode_fixed)
				ds->data_nfixed[0] = 1;
		}
		HYPER_NEW(ds->data_observations.log_prec_gaussian, tmp);
		if (mb->verbose) {
			printf("\t\tinitialise log_precision[%g]\n", ds->data_observations.log_prec_gaussian[0][0]);
			printf("\t\tfixed=[%1d]\n", ds->data_nfixed[0]);
		}
		inla_read_priorN(mb, ini, sec, &(ds->data_nprior[0]), "LOGGAMMA", 0, NULL);

		if (!ds->data_nfixed[0]) {
			mb->theta = Realloc(mb->theta, mb->ntheta + 1, double **);
			mb->theta_hyperid = Realloc(mb->theta_hyperid, mb->ntheta + 1, char *);
			mb->theta_hyperid[mb->ntheta] = ds->data_nprior[0].hyperid;
			mb->theta_tag = Realloc(mb->theta_tag, mb->ntheta + 1, char *);
			mb->theta_tag_userscale = Realloc(mb->theta_tag_userscale, mb->ntheta + 1, char *);
			mb->theta_dir = Realloc(mb->theta_dir, mb->ntheta + 1, char *);
			mb->theta_tag[mb->ntheta] = inla_make_tag("Log precision for the lognormalsurv observations", mb->ds);
			mb->theta_tag_userscale[mb->ntheta] = inla_make_tag("Precision for the lognormalsurv observations", mb->ds);
			GMRFLib_sprintf(&msg, "%s-parameter", secname);
			mb->theta_dir[mb->ntheta] = msg;

			mb->theta_from = Realloc(mb->theta_from, mb->ntheta + 1, char *);
			mb->theta_to = Realloc(mb->theta_to, mb->ntheta + 1, char *);
			mb->theta_from[mb->ntheta] = Strdup(ds->data_nprior[0].from_theta);
			mb->theta_to[mb->ntheta] = Strdup(ds->data_nprior[0].to_theta);

			mb->theta[mb->ntheta] = ds->data_observations.log_prec_gaussian;
			mb->theta_map = Realloc(mb->theta_map, mb->ntheta + 1, map_func_tp *);
			mb->theta_map[mb->ntheta] = map_precision;
			mb->theta_map_arg = Realloc(mb->theta_map_arg, mb->ntheta + 1, void *);
			mb->theta_map_arg[mb->ntheta] = NULL;
			mb->ntheta++;
			ds->data_ntheta++;
		}

		for (i = 1; i < ds->data_observations.cure_ncov + 1; i++) {
			GMRFLib_sprintf(&ctmp, "INITIAL%1d", i);
			tmp = iniparser_getdouble(ini, inla_string_join(secname, ctmp), 0.0);

			GMRFLib_sprintf(&ctmp, "FIXED%1d", i);
			ds->data_nfixed[i] = iniparser_getboolean(ini, inla_string_join(secname, ctmp), 0);
			if (!ds->data_nfixed[i] && mb->mode_use_mode) {
				tmp = mb->theta_file[mb->theta_counter_file++];
				if (mb->mode_fixed)
					ds->data_nfixed[i] = 1;
			}

			HYPER_NEW(ds->data_observations.cure_beta[i - 1], tmp);
			if (mb->verbose) {
				printf("\t\tbeta[%1d] = %g\n", i, ds->data_observations.cure_beta[i - 1][0][0]);
				printf("\t\tfixed[%1d] = %1d\n", i, ds->data_nfixed[i]);
			}
			inla_read_priorN(mb, ini, sec, &(ds->data_nprior[i]), "GAUSSIAN-std", i, NULL);

			if (!ds->data_nfixed[i]) {
				mb->theta = Realloc(mb->theta, mb->ntheta + 1, double **);
				mb->theta_hyperid = Realloc(mb->theta_hyperid, mb->ntheta + 1, char *);
				mb->theta_hyperid[mb->ntheta] = ds->data_nprior[i].hyperid;
				mb->theta_tag = Realloc(mb->theta_tag, mb->ntheta + 1, char *);
				mb->theta_tag_userscale = Realloc(mb->theta_tag_userscale, mb->ntheta + 1, char *);
				mb->theta_dir = Realloc(mb->theta_dir, mb->ntheta + 1, char *);
				GMRFLib_sprintf(&ctmp, "beta%1d for logNormal-Cure", i);

				mb->theta_tag[mb->ntheta] = inla_make_tag(ctmp, mb->ds);
				mb->theta_tag_userscale[mb->ntheta] = inla_make_tag(ctmp, mb->ds);
				GMRFLib_sprintf(&msg, "%s-parameter%1d", secname, i);
				mb->theta_dir[mb->ntheta] = msg;

				mb->theta_from = Realloc(mb->theta_from, mb->ntheta + 1, char *);
				mb->theta_to = Realloc(mb->theta_to, mb->ntheta + 1, char *);
				mb->theta_from[mb->ntheta] = Strdup(ds->data_nprior[i].from_theta);
				mb->theta_to[mb->ntheta] = Strdup(ds->data_nprior[i].to_theta);

				mb->theta[mb->ntheta] = ds->data_observations.cure_beta[i - 1];
				mb->theta_map = Realloc(mb->theta_map, mb->ntheta + 1, map_func_tp *);

				mb->theta_map[mb->ntheta] = map_identity;
				mb->theta_map_arg = Realloc(mb->theta_map_arg, mb->ntheta + 1, void *);
				mb->theta_map_arg[mb->ntheta] = NULL;

				mb->ntheta++;
				ds->data_ntheta++;
			}
		}

		// all the remaining ones are fixed
		for (i = 1 + ds->data_observations.cure_ncov; i < 1 + CURE_MAXTHETA; i++) {
			ds->data_nfixed[i] = 1;
		}
	}
		break;

	case L_GAUSSIANJW:
	{
		for (i = 0; i < 3; i++) {
			GMRFLib_sprintf(&ctmp, "FIXED%1d", i);
			iniparser_getstring(ini, inla_string_join(secname, ctmp), NULL);

			GMRFLib_sprintf(&ctmp, "INITIAL%1d", i);
			iniparser_getstring(ini, inla_string_join(secname, ctmp), NULL);

			GMRFLib_sprintf(&ctmp, "PRIOR%1d", i);
			iniparser_getstring(ini, inla_string_join(secname, ctmp), NULL);

			GMRFLib_sprintf(&ctmp, "HYPERID%1d", i);
			iniparser_getstring(ini, inla_string_join(secname, ctmp), NULL);

			GMRFLib_sprintf(&ctmp, "PARAMETERS%1d", i);
			iniparser_getstring(ini, inla_string_join(secname, ctmp), NULL);

			GMRFLib_sprintf(&ctmp, "to.theta%1d", i);
			iniparser_getstring(ini, inla_string_join(secname, ctmp), NULL);

			GMRFLib_sprintf(&ctmp, "from.theta%1d", i);
			iniparser_getstring(ini, inla_string_join(secname, ctmp), NULL);
		}

		ds->data_observations.gjw_beta = Calloc(3, double **);
		ds->data_nfixed = Calloc(3, int);
		ds->data_nprior = Calloc(3, Prior_tp);
		for (i = 0; i < 3; i++) {
			GMRFLib_sprintf(&ctmp, "INITIAL%1d", i);
			tmp = iniparser_getdouble(ini, inla_string_join(secname, ctmp), 0.0);	/* YES! */

			GMRFLib_sprintf(&ctmp, "FIXED%1d", i);
			ds->data_nfixed[i] = iniparser_getboolean(ini, inla_string_join(secname, ctmp), 0);
			if (!ds->data_nfixed[i] && mb->mode_use_mode) {
				tmp = mb->theta_file[mb->theta_counter_file++];
				if (mb->mode_fixed)
					ds->data_nfixed[i] = 1;
			}

			HYPER_NEW(ds->data_observations.gjw_beta[i], tmp);
			if (mb->verbose) {
				printf("\t\tbeta[%1d] = %g\n", i + 1, ds->data_observations.gjw_beta[i][0][0]);
				printf("\t\tfixed[%1d] = %1d\n", i, ds->data_nfixed[i]);
			}
			inla_read_priorN(mb, ini, sec, &(ds->data_nprior[i]), "GAUSSIAN-std", i, NULL);

			if (!ds->data_nfixed[i]) {
				mb->theta = Realloc(mb->theta, mb->ntheta + 1, double **);
				mb->theta_hyperid = Realloc(mb->theta_hyperid, mb->ntheta + 1, char *);
				mb->theta_hyperid[mb->ntheta] = ds->data_nprior[i].hyperid;
				mb->theta_tag = Realloc(mb->theta_tag, mb->ntheta + 1, char *);
				mb->theta_tag_userscale = Realloc(mb->theta_tag_userscale, mb->ntheta + 1, char *);
				mb->theta_dir = Realloc(mb->theta_dir, mb->ntheta + 1, char *);
				GMRFLib_sprintf(&ctmp, "beta%1d for GaussianJW observations", i);

				mb->theta_tag[mb->ntheta] = inla_make_tag(ctmp, mb->ds);
				mb->theta_tag_userscale[mb->ntheta] = inla_make_tag(ctmp, mb->ds);
				GMRFLib_sprintf(&msg, "%s-parameter%1d", secname, i);
				mb->theta_dir[mb->ntheta] = msg;

				mb->theta_from = Realloc(mb->theta_from, mb->ntheta + 1, char *);
				mb->theta_to = Realloc(mb->theta_to, mb->ntheta + 1, char *);
				mb->theta_from[mb->ntheta] = Strdup(ds->data_nprior[i].from_theta);
				mb->theta_to[mb->ntheta] = Strdup(ds->data_nprior[i].to_theta);

				mb->theta[mb->ntheta] = ds->data_observations.gjw_beta[i];
				mb->theta_map = Realloc(mb->theta_map, mb->ntheta + 1, map_func_tp *);

				mb->theta_map[mb->ntheta] = map_identity;
				mb->theta_map_arg = Realloc(mb->theta_map_arg, mb->ntheta + 1, void *);
				mb->theta_map_arg[mb->ntheta] = NULL;

				mb->ntheta++;
				ds->data_ntheta++;
			}
		}
	}
		break;

	case L_BC_GAUSSIAN:
	{
		tmp = iniparser_getdouble(ini, inla_string_join(secname, "INITIAL0"), G.log_prec_initial);
		ds->data_fixed0 = iniparser_getboolean(ini, inla_string_join(secname, "FIXED0"), 0);
		if (!ds->data_fixed0 && mb->mode_use_mode) {
			tmp = mb->theta_file[mb->theta_counter_file++];
			if (mb->mode_fixed)
				ds->data_fixed0 = 1;
		}
		HYPER_NEW(ds->data_observations.log_prec_gaussian, tmp);
		if (mb->verbose) {
			printf("\t\tinitialise log_precision[%g]\n", ds->data_observations.log_prec_gaussian[0][0]);
			printf("\t\tfixed0=[%1d]\n", ds->data_fixed0);
		}
		inla_read_prior0(mb, ini, sec, &(ds->data_prior0), "LOGGAMMA", NULL);

		/*
		 * add theta 
		 */
		if (!ds->data_fixed0) {
			mb->theta = Realloc(mb->theta, mb->ntheta + 1, double **);
			mb->theta_hyperid = Realloc(mb->theta_hyperid, mb->ntheta + 1, char *);
			mb->theta_hyperid[mb->ntheta] = ds->data_prior0.hyperid;
			mb->theta_tag = Realloc(mb->theta_tag, mb->ntheta + 1, char *);
			mb->theta_tag_userscale = Realloc(mb->theta_tag_userscale, mb->ntheta + 1, char *);
			mb->theta_dir = Realloc(mb->theta_dir, mb->ntheta + 1, char *);
			mb->theta_tag[mb->ntheta] = inla_make_tag("Log precision for the Box-Cox Gaussian observations", mb->ds);
			mb->theta_tag_userscale[mb->ntheta] = inla_make_tag("Precision for the Box-Cox Gaussian observations", mb->ds);
			GMRFLib_sprintf(&msg, "%s-parameter0", secname);
			mb->theta_dir[mb->ntheta] = msg;

			mb->theta_from = Realloc(mb->theta_from, mb->ntheta + 1, char *);
			mb->theta_to = Realloc(mb->theta_to, mb->ntheta + 1, char *);
			mb->theta_from[mb->ntheta] = Strdup(ds->data_prior0.from_theta);
			mb->theta_to[mb->ntheta] = Strdup(ds->data_prior0.to_theta);

			mb->theta[mb->ntheta] = ds->data_observations.log_prec_gaussian;
			mb->theta_map = Realloc(mb->theta_map, mb->ntheta + 1, map_func_tp *);
			mb->theta_map[mb->ntheta] = map_precision;
			mb->theta_map_arg = Realloc(mb->theta_map_arg, mb->ntheta + 1, void *);
			mb->theta_map_arg[mb->ntheta] = NULL;
			mb->ntheta++;
			ds->data_ntheta++;
		}

		tmp = iniparser_getdouble(ini, inla_string_join(secname, "INITIAL1"), 1.0);
		ds->data_fixed1 = iniparser_getboolean(ini, inla_string_join(secname, "FIXED1"), 0);
		if (!ds->data_fixed1 && mb->mode_use_mode) {
			tmp = mb->theta_file[mb->theta_counter_file++];
			if (mb->mode_fixed)
				ds->data_fixed1 = 1;
		}
		HYPER_NEW(ds->data_observations.bc_lambda, tmp);
		if (mb->verbose) {
			printf("\t\tinitialise Box-Cox parameter[%g]\n", ds->data_observations.bc_lambda[0][0]);
			printf("\t\tfixed1=[%1d]\n", ds->data_fixed1);
		}
		inla_read_prior1(mb, ini, sec, &(ds->data_prior1), "GAUSSIAN", NULL);

		if (!ds->data_fixed1) {
			mb->theta = Realloc(mb->theta, mb->ntheta + 1, double **);
			mb->theta_hyperid = Realloc(mb->theta_hyperid, mb->ntheta + 1, char *);
			mb->theta_hyperid[mb->ntheta] = ds->data_prior1.hyperid;
			mb->theta_tag = Realloc(mb->theta_tag, mb->ntheta + 1, char *);
			mb->theta_tag_userscale = Realloc(mb->theta_tag_userscale, mb->ntheta + 1, char *);
			mb->theta_dir = Realloc(mb->theta_dir, mb->ntheta + 1, char *);
			mb->theta_tag[mb->ntheta] = inla_make_tag("Box-Cox parameter", mb->ds);
			mb->theta_tag_userscale[mb->ntheta] = inla_make_tag("Box-Cox parameter", mb->ds);
			GMRFLib_sprintf(&msg, "%s-parameter1", secname);
			mb->theta_dir[mb->ntheta] = msg;

			mb->theta_from = Realloc(mb->theta_from, mb->ntheta + 1, char *);
			mb->theta_to = Realloc(mb->theta_to, mb->ntheta + 1, char *);
			mb->theta_from[mb->ntheta] = Strdup(ds->data_prior1.from_theta);
			mb->theta_to[mb->ntheta] = Strdup(ds->data_prior1.to_theta);

			mb->theta[mb->ntheta] = ds->data_observations.bc_lambda;
			mb->theta_map = Realloc(mb->theta_map, mb->ntheta + 1, map_func_tp *);
			mb->theta_map[mb->ntheta] = map_identity;
			mb->theta_map_arg = Realloc(mb->theta_map_arg, mb->ntheta + 1, void *);
			mb->theta_map_arg[mb->ntheta] = NULL;
			mb->ntheta++;
			ds->data_ntheta++;
		}
	}
		break;

	case L_EXPPOWER:
	{
		/*
		 * get options related to the generalized gaussian 
		 */
		tmp = iniparser_getdouble(ini, inla_string_join(secname, "INITIAL0"), G.log_prec_initial);
		ds->data_fixed0 = iniparser_getboolean(ini, inla_string_join(secname, "FIXED0"), 0);
		if (!ds->data_fixed0 && mb->mode_use_mode) {
			tmp = mb->theta_file[mb->theta_counter_file++];
			if (mb->mode_fixed)
				ds->data_fixed0 = 1;
		}
		HYPER_NEW(ds->data_observations.log_prec_gaussian, tmp);
		if (mb->verbose) {
			printf("\t\tinitialise log_precision[%g]\n", ds->data_observations.log_prec_gaussian[0][0]);
			printf("\t\tfixed0=[%1d]\n", ds->data_fixed0);
		}
		inla_read_prior0(mb, ini, sec, &(ds->data_prior0), "LOGGAMMA", NULL);

		/*
		 * add theta 
		 */
		if (!ds->data_fixed0) {
			mb->theta = Realloc(mb->theta, mb->ntheta + 1, double **);
			mb->theta_hyperid = Realloc(mb->theta_hyperid, mb->ntheta + 1, char *);
			mb->theta_hyperid[mb->ntheta] = ds->data_prior0.hyperid;
			mb->theta_tag = Realloc(mb->theta_tag, mb->ntheta + 1, char *);
			mb->theta_tag_userscale = Realloc(mb->theta_tag_userscale, mb->ntheta + 1, char *);
			mb->theta_dir = Realloc(mb->theta_dir, mb->ntheta + 1, char *);
			mb->theta_tag[mb->ntheta] = inla_make_tag("Log precision for exponential power observations", mb->ds);
			mb->theta_tag_userscale[mb->ntheta] = inla_make_tag("Precision for exponential power observations", mb->ds);
			GMRFLib_sprintf(&msg, "%s-parameter0", secname);
			mb->theta_dir[mb->ntheta] = msg;

			mb->theta_from = Realloc(mb->theta_from, mb->ntheta + 1, char *);
			mb->theta_to = Realloc(mb->theta_to, mb->ntheta + 1, char *);
			mb->theta_from[mb->ntheta] = Strdup(ds->data_prior0.from_theta);
			mb->theta_to[mb->ntheta] = Strdup(ds->data_prior0.to_theta);

			mb->theta[mb->ntheta] = ds->data_observations.log_prec_gaussian;
			mb->theta_map = Realloc(mb->theta_map, mb->ntheta + 1, map_func_tp *);
			mb->theta_map[mb->ntheta] = map_precision;
			mb->theta_map_arg = Realloc(mb->theta_map_arg, mb->ntheta + 1, void *);
			mb->theta_map_arg[mb->ntheta] = NULL;
			mb->ntheta++;
			ds->data_ntheta++;
		}

		tmp = iniparser_getdouble(ini, inla_string_join(secname, "INITIAL1"), 0.0);
		ds->data_fixed1 = iniparser_getboolean(ini, inla_string_join(secname, "FIXED1"), 0);
		if (!ds->data_fixed1 && mb->mode_use_mode) {
			tmp = mb->theta_file[mb->theta_counter_file++];
			if (mb->mode_fixed)
				ds->data_fixed1 = 1;
		}
		HYPER_NEW(ds->data_observations.log_power, tmp);
		if (mb->verbose) {
			printf("\t\tinitialise log_power[%g]\n", ds->data_observations.log_power[0][0]);
			printf("\t\tfixed1=[%1d]\n", ds->data_fixed1);
		}
		inla_read_prior1(mb, ini, sec, &(ds->data_prior1), "NORMAL", NULL);

		/*
		 * add theta 
		 */
		if (!ds->data_fixed1) {
			mb->theta = Realloc(mb->theta, mb->ntheta + 1, double **);
			mb->theta_hyperid = Realloc(mb->theta_hyperid, mb->ntheta + 1, char *);
			mb->theta_hyperid[mb->ntheta] = ds->data_prior1.hyperid;
			mb->theta_tag = Realloc(mb->theta_tag, mb->ntheta + 1, char *);
			mb->theta_tag_userscale = Realloc(mb->theta_tag_userscale, mb->ntheta + 1, char *);
			mb->theta_dir = Realloc(mb->theta_dir, mb->ntheta + 1, char *);
			mb->theta_tag[mb->ntheta] = inla_make_tag("Log power for exponential power observations", mb->ds);
			mb->theta_tag_userscale[mb->ntheta] = inla_make_tag("Power for exponential power observations", mb->ds);
			GMRFLib_sprintf(&msg, "%s-parameter1", secname);
			mb->theta_dir[mb->ntheta] = msg;

			mb->theta_from = Realloc(mb->theta_from, mb->ntheta + 1, char *);
			mb->theta_to = Realloc(mb->theta_to, mb->ntheta + 1, char *);
			mb->theta_from[mb->ntheta] = Strdup(ds->data_prior1.from_theta);
			mb->theta_to[mb->ntheta] = Strdup(ds->data_prior1.to_theta);

			mb->theta[mb->ntheta] = ds->data_observations.log_power;
			mb->theta_map = Realloc(mb->theta_map, mb->ntheta + 1, map_func_tp *);
			mb->theta_map[mb->ntheta] = map_one_plus_exp;	// this is power = 1 + exp(theta)
			mb->theta_map_arg = Realloc(mb->theta_map_arg, mb->ntheta + 1, void *);
			mb->theta_map_arg[mb->ntheta] = NULL;
			mb->ntheta++;
			ds->data_ntheta++;
		}
	}
		break;

	case L_EXPONENTIALSURV:
	{
		for (i = 0; i < CURE_MAXTHETA; i++) {
			GMRFLib_sprintf(&ctmp, "FIXED%1d", i);
			iniparser_getstring(ini, inla_string_join(secname, ctmp), NULL);

			GMRFLib_sprintf(&ctmp, "INITIAL%1d", i);
			iniparser_getstring(ini, inla_string_join(secname, ctmp), NULL);

			GMRFLib_sprintf(&ctmp, "PRIOR%1d", i);
			iniparser_getstring(ini, inla_string_join(secname, ctmp), NULL);

			GMRFLib_sprintf(&ctmp, "HYPERID%1d", i);
			iniparser_getstring(ini, inla_string_join(secname, ctmp), NULL);

			GMRFLib_sprintf(&ctmp, "PARAMETERS%1d", i);
			iniparser_getstring(ini, inla_string_join(secname, ctmp), NULL);

			GMRFLib_sprintf(&ctmp, "to.theta%1d", i);
			iniparser_getstring(ini, inla_string_join(secname, ctmp), NULL);

			GMRFLib_sprintf(&ctmp, "from.theta%1d", i);
			iniparser_getstring(ini, inla_string_join(secname, ctmp), NULL);
		}

		ds->data_nfixed = Calloc(CURE_MAXTHETA, int);
		ds->data_nprior = Calloc(CURE_MAXTHETA, Prior_tp);
		ds->data_observations.cure_beta = Calloc(CURE_MAXTHETA, double **);

		for (i = 0; i < ds->data_observations.cure_ncov; i++) {
			GMRFLib_sprintf(&ctmp, "INITIAL%1d", i);
			tmp = iniparser_getdouble(ini, inla_string_join(secname, ctmp), 0.0);

			GMRFLib_sprintf(&ctmp, "FIXED%1d", i);
			ds->data_nfixed[i] = iniparser_getboolean(ini, inla_string_join(secname, ctmp), 0);
			if (!ds->data_nfixed[i] && mb->mode_use_mode) {
				tmp = mb->theta_file[mb->theta_counter_file++];
				if (mb->mode_fixed)
					ds->data_nfixed[i] = 1;
			}

			HYPER_NEW(ds->data_observations.cure_beta[i], tmp);
			if (mb->verbose) {
				printf("\t\tbeta[%1d] = %g\n", i, ds->data_observations.cure_beta[i][0][0]);
				printf("\t\tfixed[%1d] = %1d\n", i, ds->data_nfixed[i]);
			}
			inla_read_priorN(mb, ini, sec, &(ds->data_nprior[i]), "GAUSSIAN-std", i, NULL);

			if (!ds->data_nfixed[i]) {
				mb->theta = Realloc(mb->theta, mb->ntheta + 1, double **);
				mb->theta_hyperid = Realloc(mb->theta_hyperid, mb->ntheta + 1, char *);
				mb->theta_hyperid[mb->ntheta] = ds->data_nprior[i].hyperid;
				mb->theta_tag = Realloc(mb->theta_tag, mb->ntheta + 1, char *);
				mb->theta_tag_userscale = Realloc(mb->theta_tag_userscale, mb->ntheta + 1, char *);
				mb->theta_dir = Realloc(mb->theta_dir, mb->ntheta + 1, char *);
				GMRFLib_sprintf(&ctmp, "beta%1d for Exp-Cure", i + 1);

				mb->theta_tag[mb->ntheta] = inla_make_tag(ctmp, mb->ds);
				mb->theta_tag_userscale[mb->ntheta] = inla_make_tag(ctmp, mb->ds);
				GMRFLib_sprintf(&msg, "%s-parameter%1d", secname, i);
				mb->theta_dir[mb->ntheta] = msg;

				mb->theta_from = Realloc(mb->theta_from, mb->ntheta + 1, char *);
				mb->theta_to = Realloc(mb->theta_to, mb->ntheta + 1, char *);
				mb->theta_from[mb->ntheta] = Strdup(ds->data_nprior[i].from_theta);
				mb->theta_to[mb->ntheta] = Strdup(ds->data_nprior[i].to_theta);

				mb->theta[mb->ntheta] = ds->data_observations.cure_beta[i];
				mb->theta_map = Realloc(mb->theta_map, mb->ntheta + 1, map_func_tp *);

				mb->theta_map[mb->ntheta] = map_identity;
				mb->theta_map_arg = Realloc(mb->theta_map_arg, mb->ntheta + 1, void *);
				mb->theta_map_arg[mb->ntheta] = NULL;

				mb->ntheta++;
				ds->data_ntheta++;
			}
		}
		// all the remaining ones are fixed
		for (i = ds->data_observations.cure_ncov; i < CURE_MAXTHETA; i++) {
			ds->data_nfixed[i] = 1;
		}

	}
		break;

	case L_SIMPLEX:
	{
		/*
		 * get options related to the gaussian 
		 */
		tmp = iniparser_getdouble(ini, inla_string_join(secname, "INITIAL"), G.log_prec_initial);
		ds->data_fixed = iniparser_getboolean(ini, inla_string_join(secname, "FIXED"), 0);
		if (!ds->data_fixed && mb->mode_use_mode) {
			tmp = mb->theta_file[mb->theta_counter_file++];
			if (mb->mode_fixed)
				ds->data_fixed = 1;
		}
		HYPER_NEW(ds->data_observations.log_prec_simplex, tmp);
		if (mb->verbose) {
			printf("\t\tinitialise log_prec[%g]\n", ds->data_observations.log_prec_simplex[0][0]);
			printf("\t\tfixed=[%1d]\n", ds->data_fixed);
		}
		inla_read_prior(mb, ini, sec, &(ds->data_prior), "LOGGAMMA", NULL);

		/*
		 * add theta 
		 */
		if (!ds->data_fixed) {
			mb->theta = Realloc(mb->theta, mb->ntheta + 1, double **);
			mb->theta_hyperid = Realloc(mb->theta_hyperid, mb->ntheta + 1, char *);
			mb->theta_hyperid[mb->ntheta] = ds->data_prior.hyperid;
			mb->theta_tag = Realloc(mb->theta_tag, mb->ntheta + 1, char *);
			mb->theta_tag_userscale = Realloc(mb->theta_tag_userscale, mb->ntheta + 1, char *);
			mb->theta_dir = Realloc(mb->theta_dir, mb->ntheta + 1, char *);
			mb->theta_tag[mb->ntheta] = inla_make_tag("Log precision for the Simplex observations", mb->ds);
			mb->theta_tag_userscale[mb->ntheta] = inla_make_tag("Precision for the Simplex observations", mb->ds);
			GMRFLib_sprintf(&msg, "%s-parameter", secname);
			mb->theta_dir[mb->ntheta] = msg;

			mb->theta_from = Realloc(mb->theta_from, mb->ntheta + 1, char *);
			mb->theta_to = Realloc(mb->theta_to, mb->ntheta + 1, char *);
			mb->theta_from[mb->ntheta] = Strdup(ds->data_prior.from_theta);
			mb->theta_to[mb->ntheta] = Strdup(ds->data_prior.to_theta);

			mb->theta[mb->ntheta] = ds->data_observations.log_prec_simplex;
			mb->theta_map = Realloc(mb->theta_map, mb->ntheta + 1, map_func_tp *);
			mb->theta_map[mb->ntheta] = map_precision;
			mb->theta_map_arg = Realloc(mb->theta_map_arg, mb->ntheta + 1, void *);
			mb->theta_map_arg[mb->ntheta] = NULL;
			mb->ntheta++;
			ds->data_ntheta++;
		}
	}
		break;

	case L_SEM:
	{
		char *Bfile = iniparser_getstring(ini, inla_string_join(secname, "CONTROL.SEM.B"), NULL);
		assert(Bfile);

		FILE *fp = fopen(Bfile, "r");
		assert(fp);

		int idx = iniparser_getint(ini, inla_string_join(secname, "CONTROL.SEM.IDX"), -1);
		int dim = 0;
		assert(fscanf(fp, "%d", &dim) == 1);
		assert(idx >= 0 && idx < dim);

		double *A = Calloc(ISQR(dim), double);
		char **B = Calloc(ISQR(dim), char *);
		size_t len = 4096L;
		char *cache = Calloc(len + 1, char);
		for (int ii = 0; ii < ISQR(dim); ii++) {
			int ret = fscanf(fp, "%lf", &(A[ii]));
			assert(ret == 1);
			ret = fscanf(fp, "%s\n", cache);
			assert(ret == 1);
			size_t len0 = strlen(cache);
			size_t len1 = len0 + 1L;
			assert(len0 <= len);
			B[ii] = Calloc(len1, char);
			Memcpy(B[ii], cache, len1);
		}
		Free(cache);
		fclose(fp);

		if (mb->verbose) {
			printf("\t\tsem.B[%1d x %1d] =\n", dim, dim);
			for (int ii = 0; ii < dim; ii++) {
				for (int jj = 0; jj < dim; jj++) {
					int kk = ii + jj * dim;
					printf("\t\t\tB[%1d, %1d] = [%8.4f] x [%s]\n", ii, jj, A[kk], B[kk]);
				}
			}
		}
		ds->data_observations.sem_dim = dim;
		ds->data_observations.sem_A = A;
		ds->data_observations.sem_B = B;
		ds->data_observations.sem_idx = idx;
		ds->data_observations.sem_cache = NULL;
	}
		break;

	case L_BELL:
	case L_STDGAUSSIAN:
	case L_POISSON:
	case L_NPOISSON:
	case L_XPOISSON:
	case L_CONTPOISSON:
		break;

	case L_QCONTPOISSON:
	{
		GMRFLib_ASSERT(ds->data_observations.quantile > 0.0 && ds->data_observations.quantile < 1.0, GMRFLib_EPARAMETER);
		ds->data_observations.qcontpoisson_func = inla_qcontpois_func(ds->data_observations.quantile, GMRFLib_CACHE_LEN());
	}
		break;

	case L_CENPOISSON:
	{
		/*
		 * get options related to the cenpoisson 
		 */
		ds->data_observations.cenpoisson_interval = Calloc(2, double);
		ctmp = Strdup(iniparser_getstring(ini, inla_string_join(secname, "CENPOISSON.I"), NULL));
		if (inla_sread_doubles(ds->data_observations.cenpoisson_interval, 2, ctmp) == INLA_FAIL) {
			inla_error_field_is_void(__GMRFLib_FuncName, secname, "CENPOISSON.I", ctmp);
		}
		if (mb->verbose) {
			printf("\t\tcenpoisson censor-interval = [%g, %g]\n",
			       ds->data_observations.cenpoisson_interval[0], ds->data_observations.cenpoisson_interval[1]);
		}
		Free(ctmp);
	}
		break;

	case L_CENPOISSON2:
	{
		/*
		 * get options related to the cenpoisson2
		 */
	}
		break;

	case L_GPOISSON:
	{
		/*
		 * get options related to the gpoisson 
		 */
		tmp = iniparser_getdouble(ini, inla_string_join(secname, "INITIAL0"), -5.0);
		ds->data_fixed0 = iniparser_getboolean(ini, inla_string_join(secname, "FIXED0"), 0);
		if (!ds->data_fixed0 && mb->mode_use_mode) {
			tmp = mb->theta_file[mb->theta_counter_file++];
			if (mb->mode_fixed)
				ds->data_fixed0 = 1;
		}
		HYPER_NEW(ds->data_observations.gpoisson_overdispersion, tmp);
		if (mb->verbose) {
			printf("\t\tinitialise log_overdispersion[%g]\n", ds->data_observations.gpoisson_overdispersion[0][0]);
			printf("\t\tfixed=[%1d]\n", ds->data_fixed0);
		}
		inla_read_prior0(mb, ini, sec, &(ds->data_prior0), "LOGGAMMA", NULL);

		/*
		 * add theta 
		 */
		if (!ds->data_fixed0) {
			mb->theta = Realloc(mb->theta, mb->ntheta + 1, double **);
			mb->theta_hyperid = Realloc(mb->theta_hyperid, mb->ntheta + 1, char *);
			mb->theta_hyperid[mb->ntheta] = ds->data_prior0.hyperid;
			mb->theta_tag = Realloc(mb->theta_tag, mb->ntheta + 1, char *);
			mb->theta_tag_userscale = Realloc(mb->theta_tag_userscale, mb->ntheta + 1, char *);
			mb->theta_dir = Realloc(mb->theta_dir, mb->ntheta + 1, char *);
			mb->theta_tag[mb->ntheta] = inla_make_tag("Log overdispersion for gpoisson", mb->ds);
			mb->theta_tag_userscale[mb->ntheta] = inla_make_tag("Overdispersion for gpoisson", mb->ds);
			GMRFLib_sprintf(&msg, "%s-parameter0", secname);
			mb->theta_dir[mb->ntheta] = msg;

			mb->theta_from = Realloc(mb->theta_from, mb->ntheta + 1, char *);
			mb->theta_to = Realloc(mb->theta_to, mb->ntheta + 1, char *);
			mb->theta_from[mb->ntheta] = Strdup(ds->data_prior0.from_theta);
			mb->theta_to[mb->ntheta] = Strdup(ds->data_prior0.to_theta);

			mb->theta[mb->ntheta] = ds->data_observations.gpoisson_overdispersion;
			mb->theta_map = Realloc(mb->theta_map, mb->ntheta + 1, map_func_tp *);
			mb->theta_map[mb->ntheta] = map_exp;
			mb->theta_map_arg = Realloc(mb->theta_map_arg, mb->ntheta + 1, void *);
			mb->theta_map_arg[mb->ntheta] = NULL;
			mb->ntheta++;
			ds->data_ntheta++;
		}

		/*
		 * the 'p' parameter
		 */
		tmp = iniparser_getdouble(ini, inla_string_join(secname, "INITIAL1"), 1);
		ds->data_fixed1 = iniparser_getboolean(ini, inla_string_join(secname, "FIXED1"), 0);
		if (!ds->data_fixed1 && mb->mode_use_mode) {
			tmp = mb->theta_file[mb->theta_counter_file++];
			if (mb->mode_fixed)
				ds->data_fixed1 = 1;
		}
		HYPER_NEW(ds->data_observations.gpoisson_p, tmp);
		if (mb->verbose) {
			printf("\t\tinitialise p[%g]\n", ds->data_observations.gpoisson_p[0][0]);
			printf("\t\tfixed=[%1d]\n", ds->data_fixed1);
		}
		inla_read_prior1(mb, ini, sec, &(ds->data_prior1), "normal", NULL);

		/*
		 * add theta 
		 */
		if (!ds->data_fixed1) {
			mb->theta = Realloc(mb->theta, mb->ntheta + 1, double **);
			mb->theta_hyperid = Realloc(mb->theta_hyperid, mb->ntheta + 1, char *);
			mb->theta_hyperid[mb->ntheta] = ds->data_prior1.hyperid;
			mb->theta_tag = Realloc(mb->theta_tag, mb->ntheta + 1, char *);
			mb->theta_tag_userscale = Realloc(mb->theta_tag_userscale, mb->ntheta + 1, char *);
			mb->theta_dir = Realloc(mb->theta_dir, mb->ntheta + 1, char *);
			mb->theta_tag[mb->ntheta] = inla_make_tag("Parameter p for gpoisson", mb->ds);
			mb->theta_tag_userscale[mb->ntheta] = inla_make_tag("Parameter p_intern for gpoisson", mb->ds);
			GMRFLib_sprintf(&msg, "%s-parameter1", secname);
			mb->theta_dir[mb->ntheta] = msg;

			mb->theta_from = Realloc(mb->theta_from, mb->ntheta + 1, char *);
			mb->theta_to = Realloc(mb->theta_to, mb->ntheta + 1, char *);
			mb->theta_from[mb->ntheta] = Strdup(ds->data_prior1.from_theta);
			mb->theta_to[mb->ntheta] = Strdup(ds->data_prior1.to_theta);

			mb->theta[mb->ntheta] = ds->data_observations.gpoisson_p;
			mb->theta_map = Realloc(mb->theta_map, mb->ntheta + 1, map_func_tp *);
			mb->theta_map[mb->ntheta] = map_identity;
			mb->theta_map_arg = Realloc(mb->theta_map_arg, mb->ntheta + 1, void *);
			mb->theta_map_arg[mb->ntheta] = NULL;
			mb->ntheta++;
			ds->data_ntheta++;
		}
	}
		break;

	case L_FMRI:
	case L_FMRISURV:
	{
		/*
		 */
		char *lname = (ds->data_id == L_FMRI ? Strdup(" fmri") : Strdup(" fmrisurv"));

		tmp = iniparser_getdouble(ini, inla_string_join(secname, "INITIAL0"), -4.0);
		ds->data_fixed0 = iniparser_getboolean(ini, inla_string_join(secname, "FIXED0"), 0);
		if (!ds->data_fixed0 && mb->mode_use_mode) {
			tmp = mb->theta_file[mb->theta_counter_file++];
			if (mb->mode_fixed)
				ds->data_fixed0 = 1;
		}
		HYPER_NEW(ds->data_observations.fmri_lprec, tmp);
		if (mb->verbose) {
			printf("\t\tinitialise log_precision[%g]\n", ds->data_observations.fmri_lprec[0][0]);
			printf("\t\tfixed=[%1d]\n", ds->data_fixed0);
		}
		inla_read_prior0(mb, ini, sec, &(ds->data_prior0), "LOGGAMMA", NULL);

		/*
		 * add theta 
		 */
		if (!ds->data_fixed0) {
			mb->theta = Realloc(mb->theta, mb->ntheta + 1, double **);
			mb->theta_hyperid = Realloc(mb->theta_hyperid, mb->ntheta + 1, char *);
			mb->theta_hyperid[mb->ntheta] = ds->data_prior0.hyperid;
			mb->theta_tag = Realloc(mb->theta_tag, mb->ntheta + 1, char *);
			mb->theta_tag_userscale = Realloc(mb->theta_tag_userscale, mb->ntheta + 1, char *);
			mb->theta_dir = Realloc(mb->theta_dir, mb->ntheta + 1, char *);
			mb->theta_tag[mb->ntheta] = inla_make_tag2("Log precision for", mb->ds, lname);
			mb->theta_tag_userscale[mb->ntheta] = inla_make_tag2("Precision for", mb->ds, lname);
			GMRFLib_sprintf(&msg, "%s-parameter0", secname);
			mb->theta_dir[mb->ntheta] = msg;

			mb->theta_from = Realloc(mb->theta_from, mb->ntheta + 1, char *);
			mb->theta_to = Realloc(mb->theta_to, mb->ntheta + 1, char *);
			mb->theta_from[mb->ntheta] = Strdup(ds->data_prior0.from_theta);
			mb->theta_to[mb->ntheta] = Strdup(ds->data_prior0.to_theta);

			mb->theta[mb->ntheta] = ds->data_observations.fmri_lprec;
			mb->theta_map = Realloc(mb->theta_map, mb->ntheta + 1, map_func_tp *);
			mb->theta_map[mb->ntheta] = map_exp;
			mb->theta_map_arg = Realloc(mb->theta_map_arg, mb->ntheta + 1, void *);
			mb->theta_map_arg[mb->ntheta] = NULL;
			mb->ntheta++;
			ds->data_ntheta++;
		}

		/*
		 * the 'dof' parameter. must be fixed
		 */
		tmp = iniparser_getdouble(ini, inla_string_join(secname, "INITIAL1"), 1);
		ds->data_fixed1 = iniparser_getboolean(ini, inla_string_join(secname, "FIXED1"), 0);
		if (!ds->data_fixed1 && mb->mode_use_mode) {
			tmp = mb->theta_file[mb->theta_counter_file++];
			if (mb->mode_fixed)
				ds->data_fixed1 = 1;
		}
		HYPER_NEW(ds->data_observations.fmri_ldof, tmp);
		if (mb->verbose) {
			printf("\t\tinitialise dof[%g]\n", ds->data_observations.fmri_ldof[0][0]);
			printf("\t\tfixed=[%1d]\n", ds->data_fixed1);
		}
		GMRFLib_ASSERT(ds->data_fixed1, GMRFLib_EPARAMETER);
		GMRFLib_ASSERT(ds->data_observations.fmri_ldof[0][0] == round(ds->data_observations.fmri_ldof[0][0]), GMRFLib_EPARAMETER);

		inla_read_prior1(mb, ini, sec, &(ds->data_prior1), "normal", NULL);
		/*
		 * add theta 
		 */
		if (!ds->data_fixed1) {
			mb->theta = Realloc(mb->theta, mb->ntheta + 1, double **);
			mb->theta_hyperid = Realloc(mb->theta_hyperid, mb->ntheta + 1, char *);
			mb->theta_hyperid[mb->ntheta] = ds->data_prior1.hyperid;
			mb->theta_tag = Realloc(mb->theta_tag, mb->ntheta + 1, char *);
			mb->theta_tag_userscale = Realloc(mb->theta_tag_userscale, mb->ntheta + 1, char *);
			mb->theta_dir = Realloc(mb->theta_dir, mb->ntheta + 1, char *);
			mb->theta_tag[mb->ntheta] = inla_make_tag2("DOF for", mb->ds, lname);
			mb->theta_tag_userscale[mb->ntheta] = inla_make_tag2("DOF for", mb->ds, lname);
			GMRFLib_sprintf(&msg, "%s-parameter1", secname);
			mb->theta_dir[mb->ntheta] = msg;

			mb->theta_from = Realloc(mb->theta_from, mb->ntheta + 1, char *);
			mb->theta_to = Realloc(mb->theta_to, mb->ntheta + 1, char *);
			mb->theta_from[mb->ntheta] = Strdup(ds->data_prior1.from_theta);
			mb->theta_to[mb->ntheta] = Strdup(ds->data_prior1.to_theta);

			mb->theta[mb->ntheta] = ds->data_observations.fmri_ldof;
			mb->theta_map = Realloc(mb->theta_map, mb->ntheta + 1, map_func_tp *);
			mb->theta_map[mb->ntheta] = map_identity;
			mb->theta_map_arg = Realloc(mb->theta_map_arg, mb->ntheta + 1, void *);
			mb->theta_map_arg[mb->ntheta] = NULL;
			mb->ntheta++;
			ds->data_ntheta++;
		}
	}
		break;

	case L_POM:
	{
#define _Q(_x) log((_x)/(1.0-(_x)))			       /* we can still use this one for default initial values */
		/*
		 * get options for the POM model. note that all theta`K' for K > nclasses-1 are not used and must be fixed no
		 * matter their input.
		 */
		int nclasses = ds->data_observations.pom_nclasses;
		ds->data_observations.pom_theta = Calloc(POM_MAXTHETA, double **);
		ds->data_nfixed = Calloc(POM_MAXTHETA, int);
		ds->data_nprior = Calloc(POM_MAXTHETA, Prior_tp);

		ds->data_observations.pom_fast_probit = iniparser_getboolean(ini, inla_string_join(secname, "POM.FAST.PROBIT"), 0);
		ctmp = iniparser_getstring(ini, inla_string_join(secname, "pom.cdf"), Strdup("DEFAULT"));
		if (!strcasecmp(ctmp, "DEFAULT")) {
			ds->data_observations.pom_cdf = POM_CDF_DEFAULT;
		} else if (!strcasecmp(ctmp, "LOGIT")) {
			ds->data_observations.pom_cdf = POM_CDF_LOGIT;
		} else if (!strcasecmp(ctmp, "PROBIT")) {
			ds->data_observations.pom_cdf = POM_CDF_PROBIT;
		} else {
			GMRFLib_ASSERT(ds->data_observations.pom_cdf == POM_CDF_DEFAULT ||
				       ds->data_observations.pom_cdf == POM_CDF_LOGIT ||
				       ds->data_observations.pom_cdf == POM_CDF_PROBIT, GMRFLib_EPARAMETER);
		}
		if (mb->verbose) {
			printf("\t\tPOM cdf = [%s]\n", (ds->data_observations.pom_cdf == POM_CDF_DEFAULT ? "default" :
							(ds->data_observations.pom_cdf == POM_CDF_LOGIT ? "logit" :
							 (ds->data_observations.pom_cdf == POM_CDF_PROBIT ? "probit" : "UNKNOWN"))));
			printf("\t\tPOM fast.probit = [%s]\n", (ds->data_observations.pom_fast_probit ? "Yes" : "No"));
		}

		if (mb->verbose) {
			printf("\t\tPOM nclasses = [%1d]\n", nclasses);
		}

		for (int count = 0; count < POM_MAXTHETA; count++) {
			GMRFLib_sprintf(&ctmp, "INITIAL%1d", count);
			tmp = iniparser_getdouble(ini, inla_string_join(secname, ctmp), NAN);
			if (ISNAN(tmp)) {
				if (count == 0) {
					tmp = _Q(1.0 / (nclasses + 1.0));
				} else {
					if (count < nclasses) {
						tmp = log(_Q((count + 1.0) / (nclasses + 1.0)) - _Q(count / (nclasses + 1.0)));
					} else {
						tmp = 0.0;     /* not in use */
					}
				}
			}

			GMRFLib_sprintf(&ctmp, "FIXED%1d", count);
			ds->data_nfixed[count] = iniparser_getboolean(ini, inla_string_join(secname, ctmp), 0);
			if (count + 1 >= ds->data_observations.pom_nclasses) {
				ds->data_nfixed[count] = 1;
			}
			if (!ds->data_nfixed[count] && mb->mode_use_mode) {
				tmp = mb->theta_file[mb->theta_counter_file++];
				if (mb->mode_fixed)
					ds->data_nfixed[count] = 1;
			}
			HYPER_NEW(ds->data_observations.pom_theta[count], tmp);

			if (mb->verbose) {
				printf("\t\tinitialise theta%1d[%g]\n", count, ds->data_observations.pom_theta[count][0][0]);
				printf("\t\tfixed%1d=[%1d]\n", count, ds->data_nfixed[count]);
			}

			if (count == 0) {
				inla_read_priorN(mb, ini, sec, &(ds->data_nprior[count]), "DIRICHLET", count, NULL);
				assert(ds->data_nprior[count].id == P_DIRICHLET);
				ds->data_nprior[count].parameters[1] = nclasses;
				ds->data_nprior[count].parameters[2] = ds->data_observations.pom_cdf;
			} else {
				inla_read_priorN(mb, ini, sec, &(ds->data_nprior[count]), "NONE", count, NULL);
				assert(ds->data_nprior[count].id == P_NONE);
			}

			/*
			 * add theta 
			 */
			if (!ds->data_nfixed[count]) {
				mb->theta = Realloc(mb->theta, mb->ntheta + 1, double **);
				mb->theta_hyperid = Realloc(mb->theta_hyperid, mb->ntheta + 1, char *);
				mb->theta_hyperid[mb->ntheta] = ds->data_nprior[count].hyperid;
				mb->theta_tag = Realloc(mb->theta_tag, mb->ntheta + 1, char *);
				mb->theta_tag_userscale = Realloc(mb->theta_tag_userscale, mb->ntheta + 1, char *);
				mb->theta_dir = Realloc(mb->theta_dir, mb->ntheta + 1, char *);

				GMRFLib_sprintf(&ctmp, "theta%1d for POM", count + 1);
				mb->theta_tag[mb->ntheta] = inla_make_tag(ctmp, mb->ds);
				GMRFLib_sprintf(&ctmp, "theta%1d for POM", count + 1);
				mb->theta_tag_userscale[mb->ntheta] = inla_make_tag(ctmp, mb->ds);
				GMRFLib_sprintf(&msg, "%s-parameter%1d", secname, count + 1);
				mb->theta_dir[mb->ntheta] = msg;

				mb->theta_from = Realloc(mb->theta_from, mb->ntheta + 1, char *);
				mb->theta_to = Realloc(mb->theta_to, mb->ntheta + 1, char *);
				mb->theta_from[mb->ntheta] = Strdup(ds->data_nprior[count].from_theta);
				mb->theta_to[mb->ntheta] = Strdup(ds->data_nprior[count].to_theta);

				mb->theta[mb->ntheta] = ds->data_observations.pom_theta[count];
				mb->theta_map = Realloc(mb->theta_map, mb->ntheta + 1, map_func_tp *);
				mb->theta_map[mb->ntheta] = map_identity;
				mb->theta_map_arg = Realloc(mb->theta_map_arg, mb->ntheta + 1, void *);
				mb->theta_map_arg[mb->ntheta] = NULL;
				mb->ntheta++;
				ds->data_ntheta++;
			}
		}

		int all_fixed = 1, all_nonfixed = 1;
		for (int count = 0; count < nclasses - 1; count++) {
			all_fixed = (all_fixed && (ds->data_nfixed[count] == 1));
			all_nonfixed = (all_nonfixed && (ds->data_nfixed[count] == 0));
		}
		if (all_fixed == 0 && all_nonfixed == 0) {
			GMRFLib_sprintf(&msg, "Hyperparameters of the POM model must either be all fixed, or all non-fixed.");
			inla_error_general(msg);
			exit(1);
		}
#undef _Q
	}
		break;

	case L_CIRCULAR_NORMAL:
	{
		/*
		 * get options related to the circular normal
		 */
		tmp = iniparser_getdouble(ini, inla_string_join(secname, "INITIAL"), G.log_prec_initial);
		ds->data_fixed = iniparser_getboolean(ini, inla_string_join(secname, "FIXED"), 0);
		if (!ds->data_fixed && mb->mode_use_mode) {
			tmp = mb->theta_file[mb->theta_counter_file++];
			if (mb->mode_fixed)
				ds->data_fixed = 1;
		}
		HYPER_NEW(ds->data_observations.log_prec_circular_normal, tmp);
		if (mb->verbose) {
			printf("\t\tinitialise log_precision parameter[%g]\n", ds->data_observations.log_prec_circular_normal[0][0]);
			printf("\t\tfixed=[%1d]\n", ds->data_fixed);
		}
		inla_read_prior(mb, ini, sec, &(ds->data_prior), "LOGGAMMA", NULL);

		/*
		 * add theta 
		 */
		if (!ds->data_fixed) {
			mb->theta = Realloc(mb->theta, mb->ntheta + 1, double **);
			mb->theta_hyperid = Realloc(mb->theta_hyperid, mb->ntheta + 1, char *);
			mb->theta_hyperid[mb->ntheta] = ds->data_prior.hyperid;
			mb->theta_tag = Realloc(mb->theta_tag, mb->ntheta + 1, char *);
			mb->theta_tag_userscale = Realloc(mb->theta_tag_userscale, mb->ntheta + 1, char *);
			mb->theta_dir = Realloc(mb->theta_dir, mb->ntheta + 1, char *);
			mb->theta_tag[mb->ntheta] = inla_make_tag("Log precision parameter for the Circular Normal observations", mb->ds);
			mb->theta_tag_userscale[mb->ntheta] = inla_make_tag("Precision parameter for the Circular Normal observations", mb->ds);
			GMRFLib_sprintf(&msg, "%s-parameter", secname);
			mb->theta_dir[mb->ntheta] = msg;

			mb->theta_from = Realloc(mb->theta_from, mb->ntheta + 1, char *);
			mb->theta_to = Realloc(mb->theta_to, mb->ntheta + 1, char *);
			mb->theta_from[mb->ntheta] = Strdup(ds->data_prior.from_theta);
			mb->theta_to[mb->ntheta] = Strdup(ds->data_prior.to_theta);

			mb->theta[mb->ntheta] = ds->data_observations.log_prec_circular_normal;
			mb->theta_map = Realloc(mb->theta_map, mb->ntheta + 1, map_func_tp *);
			mb->theta_map[mb->ntheta] = map_precision;
			mb->theta_map_arg = Realloc(mb->theta_map_arg, mb->ntheta + 1, void *);
			mb->theta_map_arg[mb->ntheta] = NULL;
			mb->ntheta++;
			ds->data_ntheta++;
		}
	}
		break;

	case L_WRAPPED_CAUCHY:
	{
		/*
		 * get options related to the circular cauchy
		 */
		tmp = iniparser_getdouble(ini, inla_string_join(secname, "INITIAL"), G.log_prec_initial);
		ds->data_fixed = iniparser_getboolean(ini, inla_string_join(secname, "FIXED"), 0);
		if (!ds->data_fixed && mb->mode_use_mode) {
			tmp = mb->theta_file[mb->theta_counter_file++];
			if (mb->mode_fixed)
				ds->data_fixed = 1;
		}
		HYPER_NEW(ds->data_observations.log_prec_wrapped_cauchy, tmp);
		if (mb->verbose) {
			printf("\t\tinitialise log_precision parameter[%g]\n", ds->data_observations.log_prec_wrapped_cauchy[0][0]);
			printf("\t\tfixed=[%1d]\n", ds->data_fixed);
		}
		inla_read_prior(mb, ini, sec, &(ds->data_prior), "LOGGAMMA", NULL);

		/*
		 * add theta 
		 */
		if (!ds->data_fixed) {
			mb->theta = Realloc(mb->theta, mb->ntheta + 1, double **);
			mb->theta_hyperid = Realloc(mb->theta_hyperid, mb->ntheta + 1, char *);
			mb->theta_hyperid[mb->ntheta] = ds->data_prior.hyperid;
			mb->theta_tag = Realloc(mb->theta_tag, mb->ntheta + 1, char *);
			mb->theta_tag_userscale = Realloc(mb->theta_tag_userscale, mb->ntheta + 1, char *);
			mb->theta_dir = Realloc(mb->theta_dir, mb->ntheta + 1, char *);
			mb->theta_tag[mb->ntheta] = inla_make_tag("Log precision parameter for the Wrapped Cauchy observations", mb->ds);
			mb->theta_tag_userscale[mb->ntheta] = inla_make_tag("Precision parameter for the Wrapped Cauchy observations", mb->ds);
			GMRFLib_sprintf(&msg, "%s-parameter", secname);
			mb->theta_dir[mb->ntheta] = msg;

			mb->theta_from = Realloc(mb->theta_from, mb->ntheta + 1, char *);
			mb->theta_to = Realloc(mb->theta_to, mb->ntheta + 1, char *);
			mb->theta_from[mb->ntheta] = Strdup(ds->data_prior.from_theta);
			mb->theta_to[mb->ntheta] = Strdup(ds->data_prior.to_theta);

			mb->theta[mb->ntheta] = ds->data_observations.log_prec_wrapped_cauchy;
			mb->theta_map = Realloc(mb->theta_map, mb->ntheta + 1, map_func_tp *);
			mb->theta_map[mb->ntheta] = map_probability;
			mb->theta_map_arg = Realloc(mb->theta_map_arg, mb->ntheta + 1, void *);
			mb->theta_map_arg[mb->ntheta] = NULL;
			mb->ntheta++;
			ds->data_ntheta++;
		}
	}
		break;

	case L_TWEEDIE:
	{
		/*
		 * get options related to the tweedie
		 */
		dtweedie_init_cache();			       // will only initialize once
		tmp = iniparser_getdouble(ini, inla_string_join(secname, "INITIAL0"), 0.0);	/* yes! */
		ds->data_fixed0 = iniparser_getboolean(ini, inla_string_join(secname, "FIXED0"), 0);
		if (!ds->data_fixed0 && mb->mode_use_mode) {
			tmp = mb->theta_file[mb->theta_counter_file++];
			if (mb->mode_fixed)
				ds->data_fixed0 = 1;
		}
		HYPER_NEW(ds->data_observations.tweedie_p_intern, tmp);
		if (mb->verbose) {
			printf("\t\tinitialise p_intern[%g]\n", ds->data_observations.tweedie_p_intern[0][0]);
			printf("\t\tfixed=[%1d]\n", ds->data_fixed0);
		}
		inla_read_prior0(mb, ini, sec, &(ds->data_prior0), "LOGGAMMA", NULL);

		/*
		 * add theta 
		 */
		if (!ds->data_fixed0) {
			mb->theta = Realloc(mb->theta, mb->ntheta + 1, double **);
			mb->theta_hyperid = Realloc(mb->theta_hyperid, mb->ntheta + 1, char *);
			mb->theta_hyperid[mb->ntheta] = ds->data_prior0.hyperid;
			mb->theta_tag = Realloc(mb->theta_tag, mb->ntheta + 1, char *);
			mb->theta_tag_userscale = Realloc(mb->theta_tag_userscale, mb->ntheta + 1, char *);
			mb->theta_dir = Realloc(mb->theta_dir, mb->ntheta + 1, char *);
			mb->theta_tag[mb->ntheta] = inla_make_tag("p_intern parameter for Tweedie", mb->ds);
			mb->theta_tag_userscale[mb->ntheta] = inla_make_tag("p parameter for Tweedie", mb->ds);
			GMRFLib_sprintf(&msg, "%s-parameter0", secname);
			mb->theta_dir[mb->ntheta] = msg;

			mb->theta_from = Realloc(mb->theta_from, mb->ntheta + 1, char *);
			mb->theta_to = Realloc(mb->theta_to, mb->ntheta + 1, char *);
			mb->theta_from[mb->ntheta] = Strdup(ds->data_prior0.from_theta);
			mb->theta_to[mb->ntheta] = Strdup(ds->data_prior0.to_theta);

			mb->theta[mb->ntheta] = ds->data_observations.tweedie_p_intern;
			mb->theta_map = Realloc(mb->theta_map, mb->ntheta + 1, map_func_tp *);
			mb->theta_map[mb->ntheta] = map_interval;
			mb->theta_map_arg = Realloc(mb->theta_map_arg, mb->ntheta + 1, void *);
			double *AB = Calloc(2, double);
			AB[0] = 1.0;
			AB[1] = 2.0;
			mb->theta_map_arg[mb->ntheta] = (void *) AB;
			mb->ntheta++;
			ds->data_ntheta++;
		}

		/*
		 * the 'dispersion' parameter
		 */
		tmp = iniparser_getdouble(ini, inla_string_join(secname, "INITIAL1"), 0.0);
		ds->data_fixed1 = iniparser_getboolean(ini, inla_string_join(secname, "FIXED1"), 0);
		if (!ds->data_fixed1 && mb->mode_use_mode) {
			tmp = mb->theta_file[mb->theta_counter_file++];
			if (mb->mode_fixed)
				ds->data_fixed1 = 1;
		}
		HYPER_NEW(ds->data_observations.tweedie_phi_intern, tmp);
		if (mb->verbose) {
			printf("\t\tinitialise phi_intern[%g]\n", ds->data_observations.tweedie_phi_intern[0][0]);
			printf("\t\tfixed=[%1d]\n", ds->data_fixed1);
		}
		inla_read_prior1(mb, ini, sec, &(ds->data_prior1), "loggamma", NULL);

		/*
		 * add theta 
		 */
		if (!ds->data_fixed1) {
			mb->theta = Realloc(mb->theta, mb->ntheta + 1, double **);
			mb->theta_hyperid = Realloc(mb->theta_hyperid, mb->ntheta + 1, char *);
			mb->theta_hyperid[mb->ntheta] = ds->data_prior1.hyperid;
			mb->theta_tag = Realloc(mb->theta_tag, mb->ntheta + 1, char *);
			mb->theta_tag_userscale = Realloc(mb->theta_tag_userscale, mb->ntheta + 1, char *);
			mb->theta_dir = Realloc(mb->theta_dir, mb->ntheta + 1, char *);
			mb->theta_tag[mb->ntheta] = inla_make_tag("Log dispersion parameter for Tweedie", mb->ds);
			mb->theta_tag_userscale[mb->ntheta] = inla_make_tag("Dispersion parameter for Tweedie", mb->ds);
			GMRFLib_sprintf(&msg, "%s-parameter1", secname);
			mb->theta_dir[mb->ntheta] = msg;

			mb->theta_from = Realloc(mb->theta_from, mb->ntheta + 1, char *);
			mb->theta_to = Realloc(mb->theta_to, mb->ntheta + 1, char *);
			mb->theta_from[mb->ntheta] = Strdup(ds->data_prior1.from_theta);
			mb->theta_to[mb->ntheta] = Strdup(ds->data_prior1.to_theta);

			mb->theta[mb->ntheta] = ds->data_observations.tweedie_phi_intern;
			mb->theta_map = Realloc(mb->theta_map, mb->ntheta + 1, map_func_tp *);
			mb->theta_map[mb->ntheta] = map_exp;
			mb->theta_map_arg = Realloc(mb->theta_map_arg, mb->ntheta + 1, void *);
			mb->theta_map_arg[mb->ntheta] = NULL;
			mb->ntheta++;
			ds->data_ntheta++;
		}
	}
		break;

	case L_GP:
	case L_DGP:
	{
		/*
		 * get options related to the gp/dgp
		 */
		GMRFLib_ASSERT(ds->data_observations.quantile > 0.0 && ds->data_observations.quantile < 1.0, GMRFLib_EPARAMETER);
		tmp = iniparser_getdouble(ini, inla_string_join(secname, "INITIAL"), -3.0);
		ds->data_fixed = iniparser_getboolean(ini, inla_string_join(secname, "FIXED"), 0);
		if (!ds->data_fixed && mb->mode_use_mode) {
			tmp = mb->theta_file[mb->theta_counter_file++];
			if (mb->mode_fixed)
				ds->data_fixed = 1;
		}
		HYPER_NEW(ds->data_observations.gp_intern_tail, tmp);
		if (mb->verbose) {
			printf("\t\tinitialise internal_tail[%g]\n", ds->data_observations.gp_intern_tail[0][0]);
			printf("\t\tfixed=[%1d]\n", ds->data_fixed);
		}
		inla_read_prior(mb, ini, sec, &(ds->data_prior), "PCGEVTAIL", NULL);
		ds->data_observations.gp_tail_interval = Calloc(2, double);
		if (ds->data_prior.id == P_PC_GEVTAIL) {
			ds->data_observations.gp_tail_interval[0] = ds->data_prior.parameters[1];
			ds->data_observations.gp_tail_interval[1] = ds->data_prior.parameters[2];
		} else {
			// used a fixed interval then
			ds->data_observations.gp_tail_interval[0] = 0.0;
			ds->data_observations.gp_tail_interval[1] = 0.5;
		}

		if (DMIN(ds->data_observations.gp_tail_interval[0], ds->data_observations.gp_tail_interval[1]) < 0.0 ||
		    DMAX(ds->data_observations.gp_tail_interval[0], ds->data_observations.gp_tail_interval[1]) >= 1.0 ||
		    ds->data_observations.gp_tail_interval[0] >= ds->data_observations.gp_tail_interval[1]) {
			inla_error_field_is_void(__GMRFLib_FuncName, secname, "TAIL.INTERVAL", ctmp);
		}
		if (mb->verbose) {
			printf("\t\tgp.tail.interval [%g %g]\n", ds->data_observations.gp_tail_interval[0],
			       ds->data_observations.gp_tail_interval[1]);
		}

		/*
		 * add theta 
		 */
		if (!ds->data_fixed) {
			mb->theta = Realloc(mb->theta, mb->ntheta + 1, double **);
			mb->theta_hyperid = Realloc(mb->theta_hyperid, mb->ntheta + 1, char *);
			mb->theta_hyperid[mb->ntheta] = ds->data_prior.hyperid;
			mb->theta_tag = Realloc(mb->theta_tag, mb->ntheta + 1, char *);
			mb->theta_tag_userscale = Realloc(mb->theta_tag_userscale, mb->ntheta + 1, char *);
			mb->theta_dir = Realloc(mb->theta_dir, mb->ntheta + 1, char *);

			if (ds->data_id == L_GP) {
				mb->theta_tag[mb->ntheta] = inla_make_tag("Intern tail parameter for the gp observations", mb->ds);
				mb->theta_tag_userscale[mb->ntheta] = inla_make_tag("Tail parameter for the gp observations", mb->ds);
			} else {
				mb->theta_tag[mb->ntheta] = inla_make_tag("Intern tail parameter for the dgp observations", mb->ds);
				mb->theta_tag_userscale[mb->ntheta] = inla_make_tag("Tail parameter for the dgp observations", mb->ds);
			}

			GMRFLib_sprintf(&msg, "%s-parameter", secname);
			mb->theta_dir[mb->ntheta] = msg;

			mb->theta_from = Realloc(mb->theta_from, mb->ntheta + 1, char *);
			mb->theta_to = Realloc(mb->theta_to, mb->ntheta + 1, char *);
			mb->theta_from[mb->ntheta] = Strdup(ds->data_prior.from_theta);
			mb->theta_to[mb->ntheta] = Strdup(ds->data_prior.to_theta);

			mb->theta[mb->ntheta] = ds->data_observations.gp_intern_tail;
			mb->theta_map = Realloc(mb->theta_map, mb->ntheta + 1, map_func_tp *);
			mb->theta_map[mb->ntheta] = map_interval;
			mb->theta_map_arg = Realloc(mb->theta_map_arg, mb->ntheta + 1, void *);
			mb->theta_map_arg[mb->ntheta] = (void *) (ds->data_observations.gp_tail_interval);
			mb->ntheta++;
			ds->data_ntheta++;
		}
	}
		break;

	case L_EGP:
	{
		/*
		 * get options related to the egp
		 */
		GMRFLib_ASSERT(ds->data_observations.quantile > 0.0 && ds->data_observations.quantile < 1.0, GMRFLib_EPARAMETER);

		tmp = iniparser_getdouble(ini, inla_string_join(secname, "INITIAL0"), 0.0);
		ds->data_fixed0 = iniparser_getboolean(ini, inla_string_join(secname, "FIXED0"), 0);
		if (!ds->data_fixed0 && mb->mode_use_mode) {
			tmp = mb->theta_file[mb->theta_counter_file++];
			if (mb->mode_fixed)
				ds->data_fixed0 = 1;
		}

		HYPER_NEW(ds->data_observations.egp_intern_tail, tmp);
		if (mb->verbose) {
			printf("\t\tinitialise internal_tail[%g]\n", ds->data_observations.egp_intern_tail[0][0]);
			printf("\t\tfixed0=[%1d]\n", ds->data_fixed0);
		}

		inla_read_prior0(mb, ini, sec, &(ds->data_prior0), "PCGEVTAIL", NULL);
		ds->data_observations.egp_tail_interval = Calloc(2, double);
		if (ds->data_prior0.id == P_PC_GEVTAIL) {
			ds->data_observations.egp_tail_interval[0] = ds->data_prior0.parameters[1];
			ds->data_observations.egp_tail_interval[1] = ds->data_prior0.parameters[2];
		} else {
			// used a fixed interval then
			ds->data_observations.egp_tail_interval[0] = -0.5;
			ds->data_observations.egp_tail_interval[1] = 0.5;
		}

		if (DMAX(ds->data_observations.egp_tail_interval[0], ds->data_observations.egp_tail_interval[1]) >= 1.0 ||
		    ds->data_observations.egp_tail_interval[0] >= ds->data_observations.egp_tail_interval[1]) {
			inla_error_field_is_void(__GMRFLib_FuncName, secname, "TAIL.INTERVAL", ctmp);
		}
		if (mb->verbose) {
			printf("\t\tegp.tail.interval [%g , %g]\n", ds->data_observations.egp_tail_interval[0],
			       ds->data_observations.egp_tail_interval[1]);
		}

		if (!ds->data_fixed0) {
			mb->theta = Realloc(mb->theta, mb->ntheta + 1, double **);
			mb->theta_hyperid = Realloc(mb->theta_hyperid, mb->ntheta + 1, char *);
			mb->theta_hyperid[mb->ntheta] = ds->data_prior0.hyperid;
			mb->theta_tag = Realloc(mb->theta_tag, mb->ntheta + 1, char *);
			mb->theta_tag_userscale = Realloc(mb->theta_tag_userscale, mb->ntheta + 1, char *);
			mb->theta_dir = Realloc(mb->theta_dir, mb->ntheta + 1, char *);

			mb->theta_tag[mb->ntheta] = inla_make_tag("Intern tail parameter for egp observations", mb->ds);
			mb->theta_tag_userscale[mb->ntheta] = inla_make_tag("Tail parameter for egp observations", mb->ds);

			GMRFLib_sprintf(&msg, "%s-parameter", secname);
			mb->theta_dir[mb->ntheta] = msg;

			mb->theta_from = Realloc(mb->theta_from, mb->ntheta + 1, char *);
			mb->theta_to = Realloc(mb->theta_to, mb->ntheta + 1, char *);
			mb->theta_from[mb->ntheta] = Strdup(ds->data_prior0.from_theta);
			mb->theta_to[mb->ntheta] = Strdup(ds->data_prior0.to_theta);

			mb->theta[mb->ntheta] = ds->data_observations.egp_intern_tail;
			mb->theta_map = Realloc(mb->theta_map, mb->ntheta + 1, map_func_tp *);
			mb->theta_map[mb->ntheta] = map_interval;
			mb->theta_map_arg = Realloc(mb->theta_map_arg, mb->ntheta + 1, void *);
			mb->theta_map_arg[mb->ntheta] = (void *) (ds->data_observations.egp_tail_interval);
			mb->ntheta++;
			ds->data_ntheta++;
		}

		tmp = iniparser_getdouble(ini, inla_string_join(secname, "INITIAL1"), 0);
		ds->data_fixed1 = iniparser_getboolean(ini, inla_string_join(secname, "FIXED1"), 0);
		if (!ds->data_fixed1 && mb->mode_use_mode) {
			tmp = mb->theta_file[mb->theta_counter_file++];
			if (mb->mode_fixed)
				ds->data_fixed1 = 1;
		}

		HYPER_NEW(ds->data_observations.egp_intern_shape, tmp);
		if (mb->verbose) {
			printf("\t\tinitialise internal_shape[%g]\n", ds->data_observations.egp_intern_shape[0][0]);
			printf("\t\tfixed1=[%1d]\n", ds->data_fixed1);
		}
		inla_read_prior1(mb, ini, sec, &(ds->data_prior1), "Gaussian", NULL);

		if (!ds->data_fixed1) {
			mb->theta = Realloc(mb->theta, mb->ntheta + 1, double **);
			mb->theta_hyperid = Realloc(mb->theta_hyperid, mb->ntheta + 1, char *);
			mb->theta_hyperid[mb->ntheta] = ds->data_prior1.hyperid;
			mb->theta_tag = Realloc(mb->theta_tag, mb->ntheta + 1, char *);
			mb->theta_tag_userscale = Realloc(mb->theta_tag_userscale, mb->ntheta + 1, char *);
			mb->theta_dir = Realloc(mb->theta_dir, mb->ntheta + 1, char *);

			mb->theta_tag[mb->ntheta] = inla_make_tag("Intern shape parameter for egp observations", mb->ds);
			mb->theta_tag_userscale[mb->ntheta] = inla_make_tag("Shape parameter for egp observations", mb->ds);

			GMRFLib_sprintf(&msg, "%s-parameter", secname);
			mb->theta_dir[mb->ntheta] = msg;

			mb->theta_from = Realloc(mb->theta_from, mb->ntheta + 1, char *);
			mb->theta_to = Realloc(mb->theta_to, mb->ntheta + 1, char *);
			mb->theta_from[mb->ntheta] = Strdup(ds->data_prior1.from_theta);
			mb->theta_to[mb->ntheta] = Strdup(ds->data_prior1.to_theta);

			mb->theta[mb->ntheta] = ds->data_observations.egp_intern_shape;
			mb->theta_map = Realloc(mb->theta_map, mb->ntheta + 1, map_func_tp *);
			mb->theta_map[mb->ntheta] = map_exp;
			mb->theta_map_arg = Realloc(mb->theta_map_arg, mb->ntheta + 1, void *);
			mb->theta_map_arg[mb->ntheta] = NULL;
			mb->ntheta++;
			ds->data_ntheta++;
		}
	}
		break;

	case L_IID_GAMMA:
	{
		/*
		 * get options related to the iid_gamma
		 */
		tmp = iniparser_getdouble(ini, inla_string_join(secname, "INITIAL0"), 0.0);	/* yes! */
		ds->data_fixed0 = iniparser_getboolean(ini, inla_string_join(secname, "FIXED0"), 0);
		if (!ds->data_fixed0 && mb->mode_use_mode) {
			tmp = mb->theta_file[mb->theta_counter_file++];
			if (mb->mode_fixed)
				ds->data_fixed0 = 1;
		}
		HYPER_NEW(ds->data_observations.iid_gamma_log_shape, tmp);
		if (mb->verbose) {
			printf("\t\tinitialise log_shape[%g]\n", ds->data_observations.iid_gamma_log_shape[0][0]);
			printf("\t\tfixed=[%1d]\n", ds->data_fixed0);
		}
		inla_read_prior0(mb, ini, sec, &(ds->data_prior0), "LOGGAMMA", NULL);

		/*
		 * add theta 
		 */
		if (!ds->data_fixed0) {
			mb->theta = Realloc(mb->theta, mb->ntheta + 1, double **);
			mb->theta_hyperid = Realloc(mb->theta_hyperid, mb->ntheta + 1, char *);
			mb->theta_hyperid[mb->ntheta] = ds->data_prior0.hyperid;
			mb->theta_tag = Realloc(mb->theta_tag, mb->ntheta + 1, char *);
			mb->theta_tag_userscale = Realloc(mb->theta_tag_userscale, mb->ntheta + 1, char *);
			mb->theta_dir = Realloc(mb->theta_dir, mb->ntheta + 1, char *);
			mb->theta_tag[mb->ntheta] = inla_make_tag("Log shape parameter for iid-gamma", mb->ds);
			mb->theta_tag_userscale[mb->ntheta] = inla_make_tag("Shape parameter for iid-gamma", mb->ds);
			GMRFLib_sprintf(&msg, "%s-parameter0", secname);
			mb->theta_dir[mb->ntheta] = msg;

			mb->theta_from = Realloc(mb->theta_from, mb->ntheta + 1, char *);
			mb->theta_to = Realloc(mb->theta_to, mb->ntheta + 1, char *);
			mb->theta_from[mb->ntheta] = Strdup(ds->data_prior0.from_theta);
			mb->theta_to[mb->ntheta] = Strdup(ds->data_prior0.to_theta);

			mb->theta[mb->ntheta] = ds->data_observations.iid_gamma_log_shape;
			mb->theta_map = Realloc(mb->theta_map, mb->ntheta + 1, map_func_tp *);
			mb->theta_map[mb->ntheta] = map_exp;
			mb->theta_map_arg = Realloc(mb->theta_map_arg, mb->ntheta + 1, void *);
			mb->theta_map_arg[mb->ntheta] = NULL;
			mb->ntheta++;
			ds->data_ntheta++;
		}

		/*
		 * the 'rate' parameter
		 */
		tmp = iniparser_getdouble(ini, inla_string_join(secname, "INITIAL1"), 0.0);
		ds->data_fixed1 = iniparser_getboolean(ini, inla_string_join(secname, "FIXED1"), 0);
		if (!ds->data_fixed1 && mb->mode_use_mode) {
			tmp = mb->theta_file[mb->theta_counter_file++];
			if (mb->mode_fixed)
				ds->data_fixed1 = 1;
		}
		HYPER_NEW(ds->data_observations.iid_gamma_log_rate, tmp);
		if (mb->verbose) {
			printf("\t\tinitialise log_rate[%g]\n", ds->data_observations.iid_gamma_log_rate[0][0]);
			printf("\t\tfixed=[%1d]\n", ds->data_fixed1);
		}
		inla_read_prior1(mb, ini, sec, &(ds->data_prior1), "loggamma", NULL);

		/*
		 * add theta 
		 */
		if (!ds->data_fixed1) {
			mb->theta = Realloc(mb->theta, mb->ntheta + 1, double **);
			mb->theta_hyperid = Realloc(mb->theta_hyperid, mb->ntheta + 1, char *);
			mb->theta_hyperid[mb->ntheta] = ds->data_prior1.hyperid;
			mb->theta_tag = Realloc(mb->theta_tag, mb->ntheta + 1, char *);
			mb->theta_tag_userscale = Realloc(mb->theta_tag_userscale, mb->ntheta + 1, char *);
			mb->theta_dir = Realloc(mb->theta_dir, mb->ntheta + 1, char *);
			mb->theta_tag[mb->ntheta] = inla_make_tag("Log rate parameter for iid-gamma", mb->ds);
			mb->theta_tag_userscale[mb->ntheta] = inla_make_tag("Rate parameter for iid-gamma", mb->ds);
			GMRFLib_sprintf(&msg, "%s-parameter1", secname);
			mb->theta_dir[mb->ntheta] = msg;

			mb->theta_from = Realloc(mb->theta_from, mb->ntheta + 1, char *);
			mb->theta_to = Realloc(mb->theta_to, mb->ntheta + 1, char *);
			mb->theta_from[mb->ntheta] = Strdup(ds->data_prior1.from_theta);
			mb->theta_to[mb->ntheta] = Strdup(ds->data_prior1.to_theta);

			mb->theta[mb->ntheta] = ds->data_observations.iid_gamma_log_rate;
			mb->theta_map = Realloc(mb->theta_map, mb->ntheta + 1, map_func_tp *);
			mb->theta_map[mb->ntheta] = map_exp;
			mb->theta_map_arg = Realloc(mb->theta_map_arg, mb->ntheta + 1, void *);
			mb->theta_map_arg[mb->ntheta] = NULL;
			mb->ntheta++;
			ds->data_ntheta++;
		}
	}
		break;

	case L_IID_LOGITBETA:
	{
		/*
		 * get options related to the iid_logitbeta. first log(a)
		 */
		tmp = iniparser_getdouble(ini, inla_string_join(secname, "INITIAL0"), 0.0);	/* yes! */
		ds->data_fixed0 = iniparser_getboolean(ini, inla_string_join(secname, "FIXED0"), 0);
		if (!ds->data_fixed0 && mb->mode_use_mode) {
			tmp = mb->theta_file[mb->theta_counter_file++];
			if (mb->mode_fixed)
				ds->data_fixed0 = 1;
		}
		HYPER_NEW(ds->data_observations.iid_logitbeta_log_a, tmp);
		if (mb->verbose) {
			printf("\t\tinitialise log_a[%g]\n", ds->data_observations.iid_logitbeta_log_a[0][0]);
			printf("\t\tfixed=[%1d]\n", ds->data_fixed0);
		}
		inla_read_prior0(mb, ini, sec, &(ds->data_prior0), "LOGGAMMA", NULL);

		/*
		 * add theta 
		 */
		if (!ds->data_fixed0) {
			mb->theta = Realloc(mb->theta, mb->ntheta + 1, double **);
			mb->theta_hyperid = Realloc(mb->theta_hyperid, mb->ntheta + 1, char *);
			mb->theta_hyperid[mb->ntheta] = ds->data_prior0.hyperid;
			mb->theta_tag = Realloc(mb->theta_tag, mb->ntheta + 1, char *);
			mb->theta_tag_userscale = Realloc(mb->theta_tag_userscale, mb->ntheta + 1, char *);
			mb->theta_dir = Realloc(mb->theta_dir, mb->ntheta + 1, char *);
			mb->theta_tag[mb->ntheta] = inla_make_tag("Log a parameter for iid-beta", mb->ds);
			mb->theta_tag_userscale[mb->ntheta] = inla_make_tag("a parameter for iid-beta", mb->ds);
			GMRFLib_sprintf(&msg, "%s-parameter0", secname);
			mb->theta_dir[mb->ntheta] = msg;

			mb->theta_from = Realloc(mb->theta_from, mb->ntheta + 1, char *);
			mb->theta_to = Realloc(mb->theta_to, mb->ntheta + 1, char *);
			mb->theta_from[mb->ntheta] = Strdup(ds->data_prior0.from_theta);
			mb->theta_to[mb->ntheta] = Strdup(ds->data_prior0.to_theta);

			mb->theta[mb->ntheta] = ds->data_observations.iid_logitbeta_log_a;
			mb->theta_map = Realloc(mb->theta_map, mb->ntheta + 1, map_func_tp *);
			mb->theta_map[mb->ntheta] = map_exp;
			mb->theta_map_arg = Realloc(mb->theta_map_arg, mb->ntheta + 1, void *);
			mb->theta_map_arg[mb->ntheta] = NULL;
			mb->ntheta++;
			ds->data_ntheta++;
		}

		/*
		 * log(b)
		 */
		tmp = iniparser_getdouble(ini, inla_string_join(secname, "INITIAL1"), 0.0);
		ds->data_fixed1 = iniparser_getboolean(ini, inla_string_join(secname, "FIXED1"), 0);
		if (!ds->data_fixed1 && mb->mode_use_mode) {
			tmp = mb->theta_file[mb->theta_counter_file++];
			if (mb->mode_fixed)
				ds->data_fixed1 = 1;
		}
		HYPER_NEW(ds->data_observations.iid_logitbeta_log_b, tmp);
		if (mb->verbose) {
			printf("\t\tinitialise log_b[%g]\n", ds->data_observations.iid_logitbeta_log_b[0][0]);
			printf("\t\tfixed=[%1d]\n", ds->data_fixed1);
		}
		inla_read_prior1(mb, ini, sec, &(ds->data_prior1), "loggamma", NULL);

		/*
		 * add theta 
		 */
		if (!ds->data_fixed1) {
			mb->theta = Realloc(mb->theta, mb->ntheta + 1, double **);
			mb->theta_hyperid = Realloc(mb->theta_hyperid, mb->ntheta + 1, char *);
			mb->theta_hyperid[mb->ntheta] = ds->data_prior1.hyperid;
			mb->theta_tag = Realloc(mb->theta_tag, mb->ntheta + 1, char *);
			mb->theta_tag_userscale = Realloc(mb->theta_tag_userscale, mb->ntheta + 1, char *);
			mb->theta_dir = Realloc(mb->theta_dir, mb->ntheta + 1, char *);
			mb->theta_tag[mb->ntheta] = inla_make_tag("Log b parameter for iid-beta", mb->ds);
			mb->theta_tag_userscale[mb->ntheta] = inla_make_tag("b parameter for iid-beta", mb->ds);
			GMRFLib_sprintf(&msg, "%s-parameter1", secname);
			mb->theta_dir[mb->ntheta] = msg;

			mb->theta_from = Realloc(mb->theta_from, mb->ntheta + 1, char *);
			mb->theta_to = Realloc(mb->theta_to, mb->ntheta + 1, char *);
			mb->theta_from[mb->ntheta] = Strdup(ds->data_prior1.from_theta);
			mb->theta_to[mb->ntheta] = Strdup(ds->data_prior1.to_theta);

			mb->theta[mb->ntheta] = ds->data_observations.iid_logitbeta_log_b;
			mb->theta_map = Realloc(mb->theta_map, mb->ntheta + 1, map_func_tp *);
			mb->theta_map[mb->ntheta] = map_exp;
			mb->theta_map_arg = Realloc(mb->theta_map_arg, mb->ntheta + 1, void *);
			mb->theta_map_arg[mb->ntheta] = NULL;
			mb->ntheta++;
			ds->data_ntheta++;
		}
	}
		break;

	case L_LOGGAMMA_FRAILTY:
	{
		/*
		 * get options related to the loggammafrailty
		 */

		tmp = iniparser_getdouble(ini, inla_string_join(secname, "INITIAL"), G.log_prec_initial);
		ds->data_fixed = iniparser_getboolean(ini, inla_string_join(secname, "FIXED"), 0);
		if (!ds->data_fixed && mb->mode_use_mode) {
			tmp = mb->theta_file[mb->theta_counter_file++];
			if (mb->mode_fixed)
				ds->data_fixed = 1;
		}
		HYPER_NEW(ds->data_observations.log_prec_loggamma_frailty, tmp);
		if (mb->verbose) {
			printf("\t\tinitialise log_precision[%g]\n", ds->data_observations.log_prec_loggamma_frailty[0][0]);
			printf("\t\tfixed=[%1d]\n", ds->data_fixed);
		}
		inla_read_prior(mb, ini, sec, &(ds->data_prior), "LOGGAMMA", NULL);

		/*
		 * add theta 
		 */
		if (!ds->data_fixed) {
			mb->theta = Realloc(mb->theta, mb->ntheta + 1, double **);
			mb->theta_hyperid = Realloc(mb->theta_hyperid, mb->ntheta + 1, char *);
			mb->theta_hyperid[mb->ntheta] = ds->data_prior.hyperid;
			mb->theta_tag = Realloc(mb->theta_tag, mb->ntheta + 1, char *);
			mb->theta_tag_userscale = Realloc(mb->theta_tag_userscale, mb->ntheta + 1, char *);
			mb->theta_dir = Realloc(mb->theta_dir, mb->ntheta + 1, char *);
			mb->theta_tag[mb->ntheta] = inla_make_tag("log precision for the gamma frailty", mb->ds);
			mb->theta_tag_userscale[mb->ntheta] = inla_make_tag("precision for the gamma frailty", mb->ds);
			GMRFLib_sprintf(&msg, "%s-parameter", secname);
			mb->theta_dir[mb->ntheta] = msg;

			mb->theta_from = Realloc(mb->theta_from, mb->ntheta + 1, char *);
			mb->theta_to = Realloc(mb->theta_to, mb->ntheta + 1, char *);
			mb->theta_from[mb->ntheta] = Strdup(ds->data_prior.from_theta);
			mb->theta_to[mb->ntheta] = Strdup(ds->data_prior.to_theta);

			mb->theta[mb->ntheta] = ds->data_observations.log_prec_loggamma_frailty;
			mb->theta_map = Realloc(mb->theta_map, mb->ntheta + 1, map_func_tp *);
			mb->theta_map[mb->ntheta] = map_precision;
			mb->theta_map_arg = Realloc(mb->theta_map_arg, mb->ntheta + 1, void *);
			mb->theta_map_arg[mb->ntheta] = NULL;
			mb->ntheta++;
			ds->data_ntheta++;
		}
	}
		break;

	case L_LOGISTIC:
	{
		/*
		 * get options related to the logistic 
		 */
		tmp = iniparser_getdouble(ini, inla_string_join(secname, "INITIAL"), G.log_prec_initial);
		ds->data_fixed = iniparser_getboolean(ini, inla_string_join(secname, "FIXED"), 0);
		if (!ds->data_fixed && mb->mode_use_mode) {
			tmp = mb->theta_file[mb->theta_counter_file++];
			if (mb->mode_fixed)
				ds->data_fixed = 1;
		}
		HYPER_NEW(ds->data_observations.log_prec_logistic, tmp);
		if (mb->verbose) {
			printf("\t\tinitialise log_precision[%g]\n", ds->data_observations.log_prec_logistic[0][0]);
			printf("\t\tfixed=[%1d]\n", ds->data_fixed);
		}
		inla_read_prior(mb, ini, sec, &(ds->data_prior), "LOGGAMMA", NULL);

		/*
		 * add theta 
		 */
		if (!ds->data_fixed) {
			mb->theta = Realloc(mb->theta, mb->ntheta + 1, double **);
			mb->theta_hyperid = Realloc(mb->theta_hyperid, mb->ntheta + 1, char *);
			mb->theta_hyperid[mb->ntheta] = ds->data_prior.hyperid;
			mb->theta_tag = Realloc(mb->theta_tag, mb->ntheta + 1, char *);
			mb->theta_tag_userscale = Realloc(mb->theta_tag_userscale, mb->ntheta + 1, char *);
			mb->theta_dir = Realloc(mb->theta_dir, mb->ntheta + 1, char *);
			mb->theta_tag[mb->ntheta] = inla_make_tag("log precision for the logistic observations", mb->ds);
			mb->theta_tag_userscale[mb->ntheta] = inla_make_tag("precision for the logistic observations", mb->ds);
			GMRFLib_sprintf(&msg, "%s-parameter", secname);
			mb->theta_dir[mb->ntheta] = msg;

			mb->theta_from = Realloc(mb->theta_from, mb->ntheta + 1, char *);
			mb->theta_to = Realloc(mb->theta_to, mb->ntheta + 1, char *);
			mb->theta_from[mb->ntheta] = Strdup(ds->data_prior.from_theta);
			mb->theta_to[mb->ntheta] = Strdup(ds->data_prior.to_theta);

			mb->theta[mb->ntheta] = ds->data_observations.log_prec_logistic;
			mb->theta_map = Realloc(mb->theta_map, mb->ntheta + 1, map_func_tp *);
			mb->theta_map[mb->ntheta] = map_precision;
			mb->theta_map_arg = Realloc(mb->theta_map_arg, mb->ntheta + 1, void *);
			mb->theta_map_arg[mb->ntheta] = NULL;
			mb->ntheta++;
			ds->data_ntheta++;
		}
	}
		break;

	case L_SKEWNORMAL:
	{
		/*
		 * get options related to the skew-normal
		 */
		tmp = iniparser_getdouble(ini, inla_string_join(secname, "INITIAL0"), 0.0);	/* yes! */
		ds->data_fixed0 = iniparser_getboolean(ini, inla_string_join(secname, "FIXED0"), 0);
		if (!ds->data_fixed0 && mb->mode_use_mode) {
			tmp = mb->theta_file[mb->theta_counter_file++];
			if (mb->mode_fixed)
				ds->data_fixed0 = 1;
		}
		HYPER_NEW(ds->data_observations.sn_lprec, tmp);
		if (mb->verbose) {
			printf("\t\tinitialise log_precision[%g]\n", ds->data_observations.sn_lprec[0][0]);
			printf("\t\tfixed=[%1d]\n", ds->data_fixed0);
		}
		inla_read_prior0(mb, ini, sec, &(ds->data_prior0), "LOGGAMMA", NULL);

		/*
		 * add theta 
		 */
		if (!ds->data_fixed0) {
			mb->theta = Realloc(mb->theta, mb->ntheta + 1, double **);
			mb->theta_hyperid = Realloc(mb->theta_hyperid, mb->ntheta + 1, char *);
			mb->theta_hyperid[mb->ntheta] = ds->data_prior0.hyperid;
			mb->theta_tag = Realloc(mb->theta_tag, mb->ntheta + 1, char *);
			mb->theta_tag_userscale = Realloc(mb->theta_tag_userscale, mb->ntheta + 1, char *);
			mb->theta_dir = Realloc(mb->theta_dir, mb->ntheta + 1, char *);
			mb->theta_tag[mb->ntheta] = inla_make_tag("log precision for skew-normal observations", mb->ds);
			mb->theta_tag_userscale[mb->ntheta] = inla_make_tag("precision for skew-normal observations", mb->ds);
			GMRFLib_sprintf(&msg, "%s-parameter0", secname);
			mb->theta_dir[mb->ntheta] = msg;

			mb->theta_from = Realloc(mb->theta_from, mb->ntheta + 1, char *);
			mb->theta_to = Realloc(mb->theta_to, mb->ntheta + 1, char *);
			mb->theta_from[mb->ntheta] = Strdup(ds->data_prior0.from_theta);
			mb->theta_to[mb->ntheta] = Strdup(ds->data_prior0.to_theta);

			mb->theta[mb->ntheta] = ds->data_observations.sn_lprec;
			mb->theta_map = Realloc(mb->theta_map, mb->ntheta + 1, map_func_tp *);
			mb->theta_map[mb->ntheta] = map_precision;
			mb->theta_map_arg = Realloc(mb->theta_map_arg, mb->ntheta + 1, void *);
			mb->theta_map_arg[mb->ntheta] = NULL;
			mb->ntheta++;
			ds->data_ntheta++;
		}

		tmp = iniparser_getdouble(ini, inla_string_join(secname, "INITIAL1"), 0.0);
		ds->data_fixed1 = iniparser_getboolean(ini, inla_string_join(secname, "FIXED1"), 0);
		if (!ds->data_fixed1 && mb->mode_use_mode) {
			tmp = mb->theta_file[mb->theta_counter_file++];
			if (mb->mode_fixed)
				ds->data_fixed1 = 1;
		}
		HYPER_NEW(ds->data_observations.sn_skew, tmp);
		if (mb->verbose) {
			printf("\t\tinitialise intern_skewness[%g]\n", ds->data_observations.sn_skew[0][0]);
			printf("\t\tfixed=[%1d]\n", ds->data_fixed1);
		}
		inla_read_prior1(mb, ini, sec, &(ds->data_prior1), "PCSN", NULL);

		/*
		 * add theta 
		 */
		if (!ds->data_fixed1) {
			mb->theta = Realloc(mb->theta, mb->ntheta + 1, double **);
			mb->theta_hyperid = Realloc(mb->theta_hyperid, mb->ntheta + 1, char *);
			mb->theta_hyperid[mb->ntheta] = ds->data_prior1.hyperid;
			mb->theta_tag = Realloc(mb->theta_tag, mb->ntheta + 1, char *);
			mb->theta_tag_userscale = Realloc(mb->theta_tag_userscale, mb->ntheta + 1, char *);
			mb->theta_dir = Realloc(mb->theta_dir, mb->ntheta + 1, char *);
			mb->theta_tag[mb->ntheta] = inla_make_tag("Intern skewness for skew-normal observations", mb->ds);
			mb->theta_tag_userscale[mb->ntheta] = inla_make_tag("Skewness for skew-normal observations", mb->ds);
			GMRFLib_sprintf(&msg, "%s-parameter1", secname);
			mb->theta_dir[mb->ntheta] = msg;

			mb->theta_from = Realloc(mb->theta_from, mb->ntheta + 1, char *);
			mb->theta_to = Realloc(mb->theta_to, mb->ntheta + 1, char *);
			mb->theta_from[mb->ntheta] = Strdup(ds->data_prior1.from_theta);
			mb->theta_to[mb->ntheta] = Strdup(ds->data_prior1.to_theta);
			mb->theta[mb->ntheta] = ds->data_observations.sn_skew;

			double *skewmax = Calloc(1, double);
			*skewmax = GMRFLib_SN_SKEWMAX;	       /* yes, this is correct */

			mb->theta_map = Realloc(mb->theta_map, mb->ntheta + 1, map_func_tp *);
			mb->theta_map[mb->ntheta] = map_phi;
			mb->theta_map_arg = Realloc(mb->theta_map_arg, mb->ntheta + 1, void *);
			mb->theta_map_arg[mb->ntheta] = (void *) skewmax;
			mb->ntheta++;
			ds->data_ntheta++;
		}
	}
		break;

	case L_GEV:
	{
		/*
		 * get options related to the gev
		 */
		tmp = iniparser_getdouble(ini, inla_string_join(secname, "gev.scale.xi"), 0.01);	/* yes! */
		assert(tmp > 0.0);
		ds->data_observations.gev_scale_xi = tmp;
		if (mb->verbose) {
			printf("\t\tgev.scale.xi [%g]\n", ds->data_observations.gev_scale_xi);
		}

		tmp = iniparser_getdouble(ini, inla_string_join(secname, "INITIAL0"), 0.0);	/* YES! */
		ds->data_fixed0 = iniparser_getboolean(ini, inla_string_join(secname, "FIXED0"), 0);
		if (!ds->data_fixed0 && mb->mode_use_mode) {
			tmp = mb->theta_file[mb->theta_counter_file++];
			if (mb->mode_fixed)
				ds->data_fixed0 = 1;
		}
		HYPER_NEW(ds->data_observations.log_prec_gev, tmp);
		if (mb->verbose) {
			printf("\t\tinitialise log_precision[%g]\n", ds->data_observations.log_prec_gev[0][0]);
			printf("\t\tfixed=[%1d]\n", ds->data_fixed0);
		}
		inla_read_prior0(mb, ini, sec, &(ds->data_prior0), "LOGGAMMA", NULL);

		/*
		 * add theta 
		 */
		if (!ds->data_fixed0) {
			mb->theta = Realloc(mb->theta, mb->ntheta + 1, double **);
			mb->theta_hyperid = Realloc(mb->theta_hyperid, mb->ntheta + 1, char *);
			mb->theta_hyperid[mb->ntheta] = ds->data_prior0.hyperid;
			mb->theta_tag = Realloc(mb->theta_tag, mb->ntheta + 1, char *);
			mb->theta_tag_userscale = Realloc(mb->theta_tag_userscale, mb->ntheta + 1, char *);
			mb->theta_dir = Realloc(mb->theta_dir, mb->ntheta + 1, char *);
			mb->theta_tag[mb->ntheta] = inla_make_tag("log precision for GEV observations", mb->ds);
			mb->theta_tag_userscale[mb->ntheta] = inla_make_tag("precision for GEV observations", mb->ds);
			GMRFLib_sprintf(&msg, "%s-parameter0", secname);
			mb->theta_dir[mb->ntheta] = msg;

			mb->theta_from = Realloc(mb->theta_from, mb->ntheta + 1, char *);
			mb->theta_to = Realloc(mb->theta_to, mb->ntheta + 1, char *);
			mb->theta_from[mb->ntheta] = Strdup(ds->data_prior0.from_theta);
			mb->theta_to[mb->ntheta] = Strdup(ds->data_prior0.to_theta);

			mb->theta[mb->ntheta] = ds->data_observations.log_prec_gev;
			mb->theta_map = Realloc(mb->theta_map, mb->ntheta + 1, map_func_tp *);
			mb->theta_map[mb->ntheta] = map_precision;
			mb->theta_map_arg = Realloc(mb->theta_map_arg, mb->ntheta + 1, void *);
			mb->theta_map_arg[mb->ntheta] = NULL;
			mb->ntheta++;
			ds->data_ntheta++;
		}

		/*
		 * the 'xi' parameter/ the gev-parameter. we need a little care, as the user see 'xi' while we work internally with 'xi/xi.scale'
		 */
		tmp = iniparser_getdouble(ini, inla_string_join(secname, "INITIAL1"), 0.0) / ds->data_observations.gev_scale_xi;	/* scale
																	 * here */
		ds->data_fixed1 = iniparser_getboolean(ini, inla_string_join(secname, "FIXED1"), 0);
		if (!ds->data_fixed1 && mb->mode_use_mode) {
			// tmp = mb->theta_file[mb->theta_counter_file++]/ds->data_observations.gev_scale_xi; /* scale here */
			if (mb->mode_fixed)
				ds->data_fixed1 = 1;
			tmp = mb->theta_file[mb->theta_counter_file++];	/* DO NOT scale here */
		}
		HYPER_NEW(ds->data_observations.xi_gev, tmp);
		if (mb->verbose) {
			printf("\t\tinitialise gev-parameter[%g] (scaled with scale.xi=[%g])\n", ds->data_observations.xi_gev[0][0],
			       ds->data_observations.gev_scale_xi);
			printf("\t\tfixed=[%1d]\n", ds->data_fixed1);
		}
		inla_read_prior1(mb, ini, sec, &(ds->data_prior1), "GAUSSIAN-a", NULL);

		/*
		 * add theta 
		 */
		if (!ds->data_fixed1) {
			mb->theta = Realloc(mb->theta, mb->ntheta + 1, double **);
			mb->theta_hyperid = Realloc(mb->theta_hyperid, mb->ntheta + 1, char *);
			mb->theta_hyperid[mb->ntheta] = ds->data_prior1.hyperid;
			mb->theta_tag = Realloc(mb->theta_tag, mb->ntheta + 1, char *);
			mb->theta_tag_userscale = Realloc(mb->theta_tag_userscale, mb->ntheta + 1, char *);
			mb->theta_dir = Realloc(mb->theta_dir, mb->ntheta + 1, char *);
			mb->theta_tag[mb->ntheta] = inla_make_tag("tail parameter for GEV observations", mb->ds);
			mb->theta_tag_userscale[mb->ntheta] = inla_make_tag("tail parameter for GEV observations", mb->ds);
			GMRFLib_sprintf(&msg, "%s-parameter1", secname);
			mb->theta_dir[mb->ntheta] = msg;

			mb->theta_from = Realloc(mb->theta_from, mb->ntheta + 1, char *);
			mb->theta_to = Realloc(mb->theta_to, mb->ntheta + 1, char *);
			mb->theta_from[mb->ntheta] = Strdup(ds->data_prior1.from_theta);
			mb->theta_to[mb->ntheta] = Strdup(ds->data_prior1.to_theta);

			mb->theta[mb->ntheta] = ds->data_observations.xi_gev;
			mb->theta_map = Realloc(mb->theta_map, mb->ntheta + 1, map_func_tp *);
			mb->theta_map[mb->ntheta] = map_identity_scale;
			mb->theta_map_arg = Realloc(mb->theta_map_arg, mb->ntheta + 1, void *);
			mb->theta_map_arg[mb->ntheta] = (void *) &(ds->data_observations.gev_scale_xi);
			mb->ntheta++;
			ds->data_ntheta++;
		}
	}
		break;

	case L_BGEV:
	{
		/*
		 * get options related to the bgev
		 */
		tmp = iniparser_getdouble(ini, inla_string_join(secname, "bgev.q.location"), 0.5);
		assert(tmp > 0.0 && tmp < 1.0);
		ds->data_observations.bgev_qlocation = tmp;
		if (mb->verbose) {
			printf("\t\tbgev.q.location [%g]\n", ds->data_observations.bgev_qlocation);
		}

		tmp = iniparser_getdouble(ini, inla_string_join(secname, "bgev.q.spread"), 0.25);
		assert(tmp > 0.0 && tmp < 1.0);
		ds->data_observations.bgev_qspread = tmp;
		if (mb->verbose) {
			printf("\t\tbgev.q.spread [%g]\n", ds->data_observations.bgev_qspread);
		}

		ctmp = iniparser_getstring(ini, inla_string_join(secname, "bgev.q.mix"), Strdup("0.10 0.20"));
		ds->data_observations.bgev_qmix = Calloc(2, double);
		if (inla_sread_doubles(ds->data_observations.bgev_qmix, 2, ctmp) == INLA_FAIL ||
		    DMIN(ds->data_observations.bgev_qmix[0], ds->data_observations.bgev_qmix[1]) <= 0.0 ||
		    DMAX(ds->data_observations.bgev_qmix[0], ds->data_observations.bgev_qmix[1]) >= 1.0 ||
		    ds->data_observations.bgev_qmix[0] >= ds->data_observations.bgev_qmix[1]) {
			inla_error_field_is_void(__GMRFLib_FuncName, secname, "BGEV.Q.MIX", ctmp);
		}
		if (mb->verbose) {
			printf("\t\tbgev.q.mix [%g %g]\n", ds->data_observations.bgev_qmix[0], ds->data_observations.bgev_qmix[1]);
		}

		tmp = iniparser_getdouble(ini, inla_string_join(secname, "bgev.beta.ab"), 5);
		ds->data_observations.bgev_beta_ab = tmp;
		if (mb->verbose) {
			printf("\t\tbgev.beta.ab [%g]\n", ds->data_observations.bgev_beta_ab);
		}

		/*
		 * mark all as read 
		 */
		for (i = 0; i < BGEV_MAXTHETA + 2; i++) {
			GMRFLib_sprintf(&ctmp, "FIXED%1d", i);
			iniparser_getstring(ini, inla_string_join(secname, ctmp), NULL);

			GMRFLib_sprintf(&ctmp, "INITIAL%1d", i);
			iniparser_getstring(ini, inla_string_join(secname, ctmp), NULL);

			GMRFLib_sprintf(&ctmp, "PRIOR%1d", i);
			iniparser_getstring(ini, inla_string_join(secname, ctmp), NULL);

			GMRFLib_sprintf(&ctmp, "HYPERID%1d", i);
			iniparser_getstring(ini, inla_string_join(secname, ctmp), NULL);

			GMRFLib_sprintf(&ctmp, "PARAMETERS%1d", i);
			iniparser_getstring(ini, inla_string_join(secname, ctmp), NULL);

			GMRFLib_sprintf(&ctmp, "to.theta%1d", i);
			iniparser_getstring(ini, inla_string_join(secname, ctmp), NULL);

			GMRFLib_sprintf(&ctmp, "from.theta%1d", i);
			iniparser_getstring(ini, inla_string_join(secname, ctmp), NULL);
		}

		/*
		 * this gives the number of betas in the spread and tail
		 */
		assert(BGEV_MAXTHETA == 10);		       /* otherwise change below... */
		ds->data_observations.bgev_nbetas[0] = ds->data_observations.attr[0];
		ds->data_observations.bgev_nbetas[1] = ds->data_observations.attr[1];
		int nbetas = ds->data_observations.bgev_nbetas[0] + ds->data_observations.bgev_nbetas[1];

		ds->data_observations.bgev_betas = Calloc(BGEV_MAXTHETA, double **);
		ds->data_nfixed = Calloc(BGEV_MAXTHETA + 2, int);	/* +2 for spread and tail */
		ds->data_nprior = Calloc(BGEV_MAXTHETA + 2, Prior_tp);

		tmp = iniparser_getdouble(ini, inla_string_join(secname, "INITIAL0"), 0.0);
		ds->data_nfixed[0] = iniparser_getboolean(ini, inla_string_join(secname, "FIXED0"), 0);
		if (!ds->data_nfixed[0] && mb->mode_use_mode) {
			tmp = mb->theta_file[mb->theta_counter_file++];
			if (mb->mode_fixed)
				ds->data_nfixed[0] = 1;
		}
		HYPER_NEW(ds->data_observations.bgev_log_spread, tmp);
		if (mb->verbose) {
			printf("\t\tinitialise log.spread[%g]\n", ds->data_observations.bgev_log_spread[0][0]);
			printf("\t\tfixed=[%1d]\n", ds->data_nfixed[0]);
		}
		inla_read_priorN(mb, ini, sec, &(ds->data_nprior[0]), "LOGGAMMA", 0, NULL);

		/*
		 * add theta 
		 */
		if (!ds->data_nfixed[0]) {
			mb->theta = Realloc(mb->theta, mb->ntheta + 1, double **);
			mb->theta_hyperid = Realloc(mb->theta_hyperid, mb->ntheta + 1, char *);
			mb->theta_hyperid[mb->ntheta] = ds->data_nprior[0].hyperid;
			mb->theta_tag = Realloc(mb->theta_tag, mb->ntheta + 1, char *);
			mb->theta_tag_userscale = Realloc(mb->theta_tag_userscale, mb->ntheta + 1, char *);
			mb->theta_dir = Realloc(mb->theta_dir, mb->ntheta + 1, char *);
			mb->theta_tag[mb->ntheta] = inla_make_tag("log spread for BGEV observations", mb->ds);
			mb->theta_tag_userscale[mb->ntheta] = inla_make_tag("spread for BGEV observations", mb->ds);
			GMRFLib_sprintf(&msg, "%s-parameter0", secname);
			mb->theta_dir[mb->ntheta] = msg;

			mb->theta_from = Realloc(mb->theta_from, mb->ntheta + 1, char *);
			mb->theta_to = Realloc(mb->theta_to, mb->ntheta + 1, char *);
			mb->theta_from[mb->ntheta] = Strdup(ds->data_nprior[0].from_theta);
			mb->theta_to[mb->ntheta] = Strdup(ds->data_nprior[0].to_theta);

			mb->theta[mb->ntheta] = ds->data_observations.bgev_log_spread;
			mb->theta_map = Realloc(mb->theta_map, mb->ntheta + 1, map_func_tp *);
			mb->theta_map[mb->ntheta] = map_exp;
			mb->theta_map_arg = Realloc(mb->theta_map_arg, mb->ntheta + 1, void *);
			mb->theta_map_arg[mb->ntheta] = NULL;
			mb->ntheta++;
			ds->data_ntheta++;
		}

		tmp = iniparser_getdouble(ini, inla_string_join(secname, "INITIAL1"), 0.0);
		ds->data_nfixed[1] = iniparser_getboolean(ini, inla_string_join(secname, "FIXED1"), 0);
		if (!ds->data_nfixed[1] && mb->mode_use_mode) {
			tmp = mb->theta_file[mb->theta_counter_file++];
			if (mb->mode_fixed)
				ds->data_nfixed[1] = 1;
		}
		HYPER_NEW(ds->data_observations.bgev_intern_tail, tmp);
		if (mb->verbose) {
			printf("\t\tinitialise internal tail[%g]\n", ds->data_observations.bgev_intern_tail[0][0]);
			printf("\t\tfixed=[%1d]\n", ds->data_nfixed[1]);
		}
		inla_read_priorN(mb, ini, sec, &(ds->data_nprior[1]), "PCGEVTAIL", 1, NULL);
		ds->data_observations.bgev_tail_interval = Calloc(2, double);
		if (ds->data_nprior[1].id == P_PC_GEVTAIL) {
			ds->data_observations.bgev_tail_interval[0] = ds->data_nprior[1].parameters[1];
			ds->data_observations.bgev_tail_interval[1] = ds->data_nprior[1].parameters[2];
		} else {
			// used a fixed interval then
			ds->data_observations.bgev_tail_interval[0] = 0.0;
			ds->data_observations.bgev_tail_interval[1] = 0.5;
		}

		if (DMIN(ds->data_observations.bgev_tail_interval[0], ds->data_observations.bgev_tail_interval[1]) < 0.0 ||
		    DMAX(ds->data_observations.bgev_tail_interval[0], ds->data_observations.bgev_tail_interval[1]) >= 1.0 ||
		    ds->data_observations.bgev_tail_interval[0] >= ds->data_observations.bgev_tail_interval[1]) {
			inla_error_field_is_void(__GMRFLib_FuncName, secname, "BGEV.TAIL.INTERVAL", ctmp);
		}
		if (mb->verbose) {
			printf("\t\tbgev.tail.interval [%g %g]\n", ds->data_observations.bgev_tail_interval[0],
			       ds->data_observations.bgev_tail_interval[1]);
		}

		/*
		 * add theta 
		 */
		if (!ds->data_nfixed[1]) {
			mb->theta = Realloc(mb->theta, mb->ntheta + 1, double **);
			mb->theta_hyperid = Realloc(mb->theta_hyperid, mb->ntheta + 1, char *);
			mb->theta_hyperid[mb->ntheta] = ds->data_nprior[1].hyperid;
			mb->theta_tag = Realloc(mb->theta_tag, mb->ntheta + 1, char *);
			mb->theta_tag_userscale = Realloc(mb->theta_tag_userscale, mb->ntheta + 1, char *);
			mb->theta_dir = Realloc(mb->theta_dir, mb->ntheta + 1, char *);
			mb->theta_tag[mb->ntheta] = inla_make_tag("intern tail for BGEV observations", mb->ds);
			mb->theta_tag_userscale[mb->ntheta] = inla_make_tag("tail for BGEV observations", mb->ds);
			GMRFLib_sprintf(&msg, "%s-parameter0", secname);
			mb->theta_dir[mb->ntheta] = msg;

			mb->theta_from = Realloc(mb->theta_from, mb->ntheta + 1, char *);
			mb->theta_to = Realloc(mb->theta_to, mb->ntheta + 1, char *);
			mb->theta_from[mb->ntheta] = Strdup(ds->data_nprior[1].from_theta);
			mb->theta_to[mb->ntheta] = Strdup(ds->data_nprior[1].to_theta);

			mb->theta[mb->ntheta] = ds->data_observations.bgev_intern_tail;
			mb->theta_map = Realloc(mb->theta_map, mb->ntheta + 1, map_func_tp *);
			mb->theta_map[mb->ntheta] = map_interval;
			mb->theta_map_arg = Realloc(mb->theta_map_arg, mb->ntheta + 1, void *);
			mb->theta_map_arg[mb->ntheta] = (void *) (ds->data_observations.bgev_tail_interval);
			mb->ntheta++;
			ds->data_ntheta++;
		}

		for (i = 0; i < nbetas; i++) {
			int ii, idx = i + 2;

			GMRFLib_sprintf(&ctmp, "INITIAL%1d", idx);
			tmp = iniparser_getdouble(ini, inla_string_join(secname, ctmp), 0.0);	/* YES! */

			GMRFLib_sprintf(&ctmp, "FIXED%1d", idx);
			ds->data_nfixed[idx] = iniparser_getboolean(ini, inla_string_join(secname, ctmp), 0);
			if (!ds->data_nfixed[idx] && mb->mode_use_mode) {
				tmp = mb->theta_file[mb->theta_counter_file++];
				if (mb->mode_fixed)
					ds->data_nfixed[idx] = 1;
			}

			if (i < ds->data_observations.bgev_nbetas[0]) {
				// log-spread part, as normal
				HYPER_NEW(ds->data_observations.bgev_betas[i], tmp);
				if (mb->verbose) {
					printf("\t\tbetas[%1d] (spread) = %g\n", i, ds->data_observations.bgev_betas[i][0][0]);
					printf("\t\tfixed[%1d] = %1d\n", i, ds->data_nfixed[i]);
				}
			} else {
				// tail part
				HYPER_NEW(ds->data_observations.bgev_betas[i], tmp);
				if (mb->verbose) {
					printf("\t\tbetas[%1d] (tail) = %g\n", i, ds->data_observations.bgev_betas[i][0][0]);
					printf("\t\tfixed[%1d]= %1d\n", i, ds->data_nfixed[i]);
				}
			}

			inla_read_priorN(mb, ini, sec, &(ds->data_nprior[idx]), "GAUSSIAN-std", idx, NULL);

			if (!ds->data_nfixed[idx]) {
				mb->theta = Realloc(mb->theta, mb->ntheta + 1, double **);
				mb->theta_hyperid = Realloc(mb->theta_hyperid, mb->ntheta + 1, char *);
				mb->theta_hyperid[mb->ntheta] = ds->data_nprior[idx].hyperid;
				mb->theta_tag = Realloc(mb->theta_tag, mb->ntheta + 1, char *);
				mb->theta_tag_userscale = Realloc(mb->theta_tag_userscale, mb->ntheta + 1, char *);
				mb->theta_dir = Realloc(mb->theta_dir, mb->ntheta + 1, char *);
				if (i < ds->data_observations.bgev_nbetas[0]) {
					ii = i + 1;
					GMRFLib_sprintf(&ctmp, "beta%1d (spread) for BGEV observations", ii);
				} else {
					ii = i + 1;
					GMRFLib_sprintf(&ctmp, "beta%1d (tail) for BGEV observations", ii);
				}

				mb->theta_tag[mb->ntheta] = inla_make_tag(ctmp, mb->ds);
				mb->theta_tag_userscale[mb->ntheta] = inla_make_tag(ctmp, mb->ds);
				GMRFLib_sprintf(&msg, "%s-parameter%1d", secname, i);
				mb->theta_dir[mb->ntheta] = msg;

				mb->theta_from = Realloc(mb->theta_from, mb->ntheta + 1, char *);
				mb->theta_to = Realloc(mb->theta_to, mb->ntheta + 1, char *);
				mb->theta_from[mb->ntheta] = Strdup(ds->data_nprior[idx].from_theta);
				mb->theta_to[mb->ntheta] = Strdup(ds->data_nprior[idx].to_theta);

				mb->theta[mb->ntheta] = ds->data_observations.bgev_betas[i];
				mb->theta_map = Realloc(mb->theta_map, mb->ntheta + 1, map_func_tp *);

				if (i < ds->data_observations.bgev_nbetas[0]) {
					// spread
					mb->theta_map[mb->ntheta] = map_identity;
					mb->theta_map_arg = Realloc(mb->theta_map_arg, mb->ntheta + 1, void *);
					mb->theta_map_arg[mb->ntheta] = NULL;
				} else {
					// tail
					mb->theta_map[mb->ntheta] = map_identity;
					mb->theta_map_arg = Realloc(mb->theta_map_arg, mb->ntheta + 1, void *);
					mb->theta_map_arg[mb->ntheta] = (void *) NULL;
				}

				mb->ntheta++;
				ds->data_ntheta++;
			}
		}

	}
		break;

	case L_GGAUSSIAN:
	case L_GGAUSSIANS:
	{
		for (i = 0; i < GGAUSSIAN_MAXTHETA; i++) {
			GMRFLib_sprintf(&ctmp, "FIXED%1d", i);
			iniparser_getstring(ini, inla_string_join(secname, ctmp), NULL);

			GMRFLib_sprintf(&ctmp, "INITIAL%1d", i);
			iniparser_getstring(ini, inla_string_join(secname, ctmp), NULL);

			GMRFLib_sprintf(&ctmp, "PRIOR%1d", i);
			iniparser_getstring(ini, inla_string_join(secname, ctmp), NULL);

			GMRFLib_sprintf(&ctmp, "HYPERID%1d", i);
			iniparser_getstring(ini, inla_string_join(secname, ctmp), NULL);

			GMRFLib_sprintf(&ctmp, "PARAMETERS%1d", i);
			iniparser_getstring(ini, inla_string_join(secname, ctmp), NULL);

			GMRFLib_sprintf(&ctmp, "to.theta%1d", i);
			iniparser_getstring(ini, inla_string_join(secname, ctmp), NULL);

			GMRFLib_sprintf(&ctmp, "from.theta%1d", i);
			iniparser_getstring(ini, inla_string_join(secname, ctmp), NULL);
		}

		char *link_simple = iniparser_getstring(ini, inla_string_join(secname, "LINK.SIMPLE"), NULL);
		ds->data_observations.link_simple_name = link_simple;
		if (!strcasecmp(link_simple, "IDENTITY")) {
			ds->data_observations.link_simple_invlinkfunc = link_identity;
		} else if (!strcasecmp(link_simple, "LOG")) {
			ds->data_observations.link_simple_invlinkfunc = link_log;
		} else if (!strcasecmp(link_simple, "PROBIT")) {
			ds->data_observations.link_simple_invlinkfunc = link_probit;
		} else if (!strcasecmp(link_simple, "CLOGLOG")) {
			ds->data_observations.link_simple_invlinkfunc = link_cloglog;
		} else if (!strcasecmp(link_simple, "LOGIT")) {
			ds->data_observations.link_simple_invlinkfunc = link_logit;
		} else {
			GMRFLib_sprintf(&msg, "%s: ggaussian(S) likelihood: no valid link.simple[%s]", secname, link_simple);
			inla_error_general(msg);
			exit(1);
		}
		if (mb->verbose) {
			printf("\t\tlink.simple[%s]\n", ds->data_observations.link_simple_name);
		}

		const char *suff = Strdup((ds->data_id == L_GGAUSSIAN ? "" : "S"));
		ds->data_nfixed = Calloc(GGAUSSIAN_MAXTHETA + 1, int);
		ds->data_nprior = Calloc(GGAUSSIAN_MAXTHETA, Prior_tp);
		ds->data_observations.ggaussian_beta = Calloc(GGAUSSIAN_MAXTHETA, double **);

		for (i = 0; i < ds->data_observations.ggaussian_nbeta; i++) {
			GMRFLib_sprintf(&ctmp, "INITIAL%1d", i);
			tmp = iniparser_getdouble(ini, inla_string_join(secname, ctmp), 0.0);	/* YES! */

			GMRFLib_sprintf(&ctmp, "FIXED%1d", i);
			ds->data_nfixed[i] = iniparser_getboolean(ini, inla_string_join(secname, ctmp), 0);
			if (!ds->data_nfixed[i] && mb->mode_use_mode) {
				tmp = mb->theta_file[mb->theta_counter_file++];
				if (mb->mode_fixed)
					ds->data_nfixed[i] = 1;
			}

			HYPER_NEW(ds->data_observations.ggaussian_beta[i], tmp);
			if (mb->verbose) {
				printf("\t\tbeta[%1d] = %g\n", i, ds->data_observations.ggaussian_beta[i][0][0]);
				printf("\t\tfixed[%1d] = %1d\n", i, ds->data_nfixed[i]);
			}

			inla_read_priorN(mb, ini, sec, &(ds->data_nprior[i]), "GAUSSIAN-std", i, NULL);

			if (!ds->data_nfixed[i]) {
				mb->theta = Realloc(mb->theta, mb->ntheta + 1, double **);
				mb->theta_hyperid = Realloc(mb->theta_hyperid, mb->ntheta + 1, char *);
				mb->theta_hyperid[mb->ntheta] = ds->data_nprior[i].hyperid;
				mb->theta_tag = Realloc(mb->theta_tag, mb->ntheta + 1, char *);
				mb->theta_tag_userscale = Realloc(mb->theta_tag_userscale, mb->ntheta + 1, char *);
				mb->theta_dir = Realloc(mb->theta_dir, mb->ntheta + 1, char *);

				GMRFLib_sprintf(&ctmp, "beta%1d for ggaussian%1s observations", i + 1, suff);
				mb->theta_tag[mb->ntheta] = inla_make_tag(ctmp, mb->ds);
				mb->theta_tag_userscale[mb->ntheta] = inla_make_tag(ctmp, mb->ds);
				GMRFLib_sprintf(&msg, "%s-parameter%1d", secname, i);
				mb->theta_dir[mb->ntheta] = msg;

				mb->theta_from = Realloc(mb->theta_from, mb->ntheta + 1, char *);
				mb->theta_to = Realloc(mb->theta_to, mb->ntheta + 1, char *);
				mb->theta_from[mb->ntheta] = Strdup(ds->data_nprior[i].from_theta);
				mb->theta_to[mb->ntheta] = Strdup(ds->data_nprior[i].to_theta);

				mb->theta[mb->ntheta] = ds->data_observations.ggaussian_beta[i];
				mb->theta_map = Realloc(mb->theta_map, mb->ntheta + 1, map_func_tp *);

				mb->theta_map[mb->ntheta] = map_identity;
				mb->theta_map_arg = Realloc(mb->theta_map_arg, mb->ntheta + 1, void *);
				mb->theta_map_arg[mb->ntheta] = NULL;
				mb->ntheta++;
				ds->data_ntheta++;
			}
		}
	}
		break;

	case L_RCPOISSON:
	{
		for (i = 0; i < RCPOISSON_MAXTHETA; i++) {
			GMRFLib_sprintf(&ctmp, "FIXED%1d", i);
			iniparser_getstring(ini, inla_string_join(secname, ctmp), NULL);

			GMRFLib_sprintf(&ctmp, "INITIAL%1d", i);
			iniparser_getstring(ini, inla_string_join(secname, ctmp), NULL);

			GMRFLib_sprintf(&ctmp, "PRIOR%1d", i);
			iniparser_getstring(ini, inla_string_join(secname, ctmp), NULL);

			GMRFLib_sprintf(&ctmp, "HYPERID%1d", i);
			iniparser_getstring(ini, inla_string_join(secname, ctmp), NULL);

			GMRFLib_sprintf(&ctmp, "PARAMETERS%1d", i);
			iniparser_getstring(ini, inla_string_join(secname, ctmp), NULL);

			GMRFLib_sprintf(&ctmp, "to.theta%1d", i);
			iniparser_getstring(ini, inla_string_join(secname, ctmp), NULL);

			GMRFLib_sprintf(&ctmp, "from.theta%1d", i);
			iniparser_getstring(ini, inla_string_join(secname, ctmp), NULL);
		}

		ds->data_nfixed = Calloc(RCPOISSON_MAXTHETA, int);
		ds->data_nprior = Calloc(RCPOISSON_MAXTHETA, Prior_tp);
		ds->data_observations.rcp_beta = Calloc(RCPOISSON_MAXTHETA, double **);

		for (i = ds->data_observations.rcp_nbeta; i < RCPOISSON_MAXTHETA; i++) {
			ds->data_nfixed[i] = 1;
		}

		for (i = 0; i < ds->data_observations.rcp_nbeta; i++) {
			GMRFLib_sprintf(&ctmp, "INITIAL%1d", i);
			tmp = iniparser_getdouble(ini, inla_string_join(secname, ctmp), 0.0);	/* YES! */

			GMRFLib_sprintf(&ctmp, "FIXED%1d", i);
			ds->data_nfixed[i] = iniparser_getboolean(ini, inla_string_join(secname, ctmp), 0);
			if (!ds->data_nfixed[i] && mb->mode_use_mode) {
				tmp = mb->theta_file[mb->theta_counter_file++];
				if (mb->mode_fixed)
					ds->data_nfixed[i] = 1;
			}

			HYPER_NEW(ds->data_observations.rcp_beta[i], tmp);
			if (mb->verbose) {
				printf("\t\tbeta[%1d] = %g\n", i, ds->data_observations.rcp_beta[i][0][0]);
				printf("\t\tfixed[%1d] = %1d\n", i, ds->data_nfixed[i]);
			}

			inla_read_priorN(mb, ini, sec, &(ds->data_nprior[i]), "GAUSSIAN-std", i, NULL);

			if (!ds->data_nfixed[i]) {
				mb->theta = Realloc(mb->theta, mb->ntheta + 1, double **);
				mb->theta_hyperid = Realloc(mb->theta_hyperid, mb->ntheta + 1, char *);
				mb->theta_hyperid[mb->ntheta] = ds->data_nprior[i].hyperid;
				mb->theta_tag = Realloc(mb->theta_tag, mb->ntheta + 1, char *);
				mb->theta_tag_userscale = Realloc(mb->theta_tag_userscale, mb->ntheta + 1, char *);
				mb->theta_dir = Realloc(mb->theta_dir, mb->ntheta + 1, char *);

				GMRFLib_sprintf(&ctmp, "beta%1d for rcpoisson observations", i + 1);
				mb->theta_tag[mb->ntheta] = inla_make_tag(ctmp, mb->ds);
				mb->theta_tag_userscale[mb->ntheta] = inla_make_tag(ctmp, mb->ds);
				GMRFLib_sprintf(&msg, "%s-parameter%1d", secname, i);
				mb->theta_dir[mb->ntheta] = msg;

				mb->theta_from = Realloc(mb->theta_from, mb->ntheta + 1, char *);
				mb->theta_to = Realloc(mb->theta_to, mb->ntheta + 1, char *);
				mb->theta_from[mb->ntheta] = Strdup(ds->data_nprior[i].from_theta);
				mb->theta_to[mb->ntheta] = Strdup(ds->data_nprior[i].to_theta);

				mb->theta[mb->ntheta] = ds->data_observations.rcp_beta[i];
				mb->theta_map = Realloc(mb->theta_map, mb->ntheta + 1, map_func_tp *);

				mb->theta_map[mb->ntheta] = map_identity;
				mb->theta_map_arg = Realloc(mb->theta_map_arg, mb->ntheta + 1, void *);
				mb->theta_map_arg[mb->ntheta] = NULL;
				mb->ntheta++;
				ds->data_ntheta++;
			}
		}
	}
		break;

	case L_TPOISSON:
	{
		for (i = 0; i < TPOISSON_MAXTHETA; i++) {
			GMRFLib_sprintf(&ctmp, "FIXED%1d", i);
			iniparser_getstring(ini, inla_string_join(secname, ctmp), NULL);

			GMRFLib_sprintf(&ctmp, "INITIAL%1d", i);
			iniparser_getstring(ini, inla_string_join(secname, ctmp), NULL);

			GMRFLib_sprintf(&ctmp, "PRIOR%1d", i);
			iniparser_getstring(ini, inla_string_join(secname, ctmp), NULL);

			GMRFLib_sprintf(&ctmp, "HYPERID%1d", i);
			iniparser_getstring(ini, inla_string_join(secname, ctmp), NULL);

			GMRFLib_sprintf(&ctmp, "PARAMETERS%1d", i);
			iniparser_getstring(ini, inla_string_join(secname, ctmp), NULL);

			GMRFLib_sprintf(&ctmp, "to.theta%1d", i);
			iniparser_getstring(ini, inla_string_join(secname, ctmp), NULL);

			GMRFLib_sprintf(&ctmp, "from.theta%1d", i);
			iniparser_getstring(ini, inla_string_join(secname, ctmp), NULL);
		}

		ds->data_nfixed = Calloc(TPOISSON_MAXTHETA, int);
		ds->data_nprior = Calloc(TPOISSON_MAXTHETA, Prior_tp);
		ds->data_observations.tp_beta = Calloc(TPOISSON_MAXTHETA, double **);

		for (i = ds->data_observations.tp_nbeta; i < TPOISSON_MAXTHETA; i++) {
			ds->data_nfixed[i] = 1;
		}

		for (i = 0; i < ds->data_observations.tp_nbeta; i++) {
			GMRFLib_sprintf(&ctmp, "INITIAL%1d", i);
			tmp = iniparser_getdouble(ini, inla_string_join(secname, ctmp), 0.0);	/* YES! */

			GMRFLib_sprintf(&ctmp, "FIXED%1d", i);
			ds->data_nfixed[i] = iniparser_getboolean(ini, inla_string_join(secname, ctmp), 0);
			if (!ds->data_nfixed[i] && mb->mode_use_mode) {
				tmp = mb->theta_file[mb->theta_counter_file++];
				if (mb->mode_fixed)
					ds->data_nfixed[i] = 1;
			}

			HYPER_NEW(ds->data_observations.tp_beta[i], tmp);
			if (mb->verbose) {
				printf("\t\tbeta[%1d] = %g\n", i, ds->data_observations.tp_beta[i][0][0]);
				printf("\t\tfixed[%1d] = %1d\n", i, ds->data_nfixed[i]);
			}

			inla_read_priorN(mb, ini, sec, &(ds->data_nprior[i]), "GAUSSIAN-std", i, NULL);

			if (!ds->data_nfixed[i]) {
				mb->theta = Realloc(mb->theta, mb->ntheta + 1, double **);
				mb->theta_hyperid = Realloc(mb->theta_hyperid, mb->ntheta + 1, char *);
				mb->theta_hyperid[mb->ntheta] = ds->data_nprior[i].hyperid;
				mb->theta_tag = Realloc(mb->theta_tag, mb->ntheta + 1, char *);
				mb->theta_tag_userscale = Realloc(mb->theta_tag_userscale, mb->ntheta + 1, char *);
				mb->theta_dir = Realloc(mb->theta_dir, mb->ntheta + 1, char *);

				GMRFLib_sprintf(&ctmp, "beta%1d for tpoisson observations", i + 1);
				mb->theta_tag[mb->ntheta] = inla_make_tag(ctmp, mb->ds);
				mb->theta_tag_userscale[mb->ntheta] = inla_make_tag(ctmp, mb->ds);
				GMRFLib_sprintf(&msg, "%s-parameter%1d", secname, i);
				mb->theta_dir[mb->ntheta] = msg;

				mb->theta_from = Realloc(mb->theta_from, mb->ntheta + 1, char *);
				mb->theta_to = Realloc(mb->theta_to, mb->ntheta + 1, char *);
				mb->theta_from[mb->ntheta] = Strdup(ds->data_nprior[i].from_theta);
				mb->theta_to[mb->ntheta] = Strdup(ds->data_nprior[i].to_theta);

				mb->theta[mb->ntheta] = ds->data_observations.tp_beta[i];
				mb->theta_map = Realloc(mb->theta_map, mb->ntheta + 1, map_func_tp *);

				mb->theta_map[mb->ntheta] = map_identity;
				mb->theta_map_arg = Realloc(mb->theta_map_arg, mb->ntheta + 1, void *);
				mb->theta_map_arg[mb->ntheta] = NULL;
				mb->ntheta++;
				ds->data_ntheta++;
			}
		}
	}
		break;

	case L_0POISSON:
	case L_0POISSONS:
	{
		for (i = 0; i < POISSON0_MAXTHETA; i++) {
			GMRFLib_sprintf(&ctmp, "FIXED%1d", i);
			iniparser_getstring(ini, inla_string_join(secname, ctmp), NULL);

			GMRFLib_sprintf(&ctmp, "INITIAL%1d", i);
			iniparser_getstring(ini, inla_string_join(secname, ctmp), NULL);

			GMRFLib_sprintf(&ctmp, "PRIOR%1d", i);
			iniparser_getstring(ini, inla_string_join(secname, ctmp), NULL);

			GMRFLib_sprintf(&ctmp, "HYPERID%1d", i);
			iniparser_getstring(ini, inla_string_join(secname, ctmp), NULL);

			GMRFLib_sprintf(&ctmp, "PARAMETERS%1d", i);
			iniparser_getstring(ini, inla_string_join(secname, ctmp), NULL);

			GMRFLib_sprintf(&ctmp, "to.theta%1d", i);
			iniparser_getstring(ini, inla_string_join(secname, ctmp), NULL);

			GMRFLib_sprintf(&ctmp, "from.theta%1d", i);
			iniparser_getstring(ini, inla_string_join(secname, ctmp), NULL);
		}

		char *link_simple = iniparser_getstring(ini, inla_string_join(secname, "LINK.SIMPLE"), NULL);
		ds->data_observations.link_simple_name = link_simple;
		if (!strcasecmp(link_simple, "IDENTITY")) {
			ds->data_observations.link_simple_invlinkfunc = link_identity;
		} else if (!strcasecmp(link_simple, "LOG")) {
			ds->data_observations.link_simple_invlinkfunc = link_log;
		} else if (!strcasecmp(link_simple, "PROBIT")) {
			ds->data_observations.link_simple_invlinkfunc = link_probit;
		} else if (!strcasecmp(link_simple, "CLOGLOG")) {
			ds->data_observations.link_simple_invlinkfunc = link_cloglog;
		} else if (!strcasecmp(link_simple, "LOGIT")) {
			ds->data_observations.link_simple_invlinkfunc = link_logit;
		} else {
			GMRFLib_sprintf(&msg, "%s: 0poisson(S) likelihood: no valid link.simple[%s]", secname, link_simple);
			inla_error_general(msg);
			exit(1);
		}
		if (mb->verbose) {
			printf("\t\tlink.simple[%s]\n", ds->data_observations.link_simple_name);
		}

		const char *suff = Strdup((ds->data_id == L_0POISSON ? "" : "S"));
		ds->data_nfixed = Calloc(POISSON0_MAXTHETA + 1, int);
		ds->data_nprior = Calloc(POISSON0_MAXTHETA, Prior_tp);
		ds->data_observations.poisson0_beta = Calloc(POISSON0_MAXTHETA, double **);

		for (i = 0; i < ds->data_observations.poisson0_nbeta; i++) {
			GMRFLib_sprintf(&ctmp, "INITIAL%1d", i);
			tmp = iniparser_getdouble(ini, inla_string_join(secname, ctmp), 0.0);	/* YES! */

			GMRFLib_sprintf(&ctmp, "FIXED%1d", i);
			ds->data_nfixed[i] = iniparser_getboolean(ini, inla_string_join(secname, ctmp), 0);
			if (!ds->data_nfixed[i] && mb->mode_use_mode) {
				tmp = mb->theta_file[mb->theta_counter_file++];
				if (mb->mode_fixed)
					ds->data_nfixed[i] = 1;
			}

			HYPER_NEW(ds->data_observations.poisson0_beta[i], tmp);
			if (mb->verbose) {
				printf("\t\tbeta[%1d] = %g\n", i, ds->data_observations.poisson0_beta[i][0][0]);
				printf("\t\tfixed[%1d] = %1d\n", i, ds->data_nfixed[i]);
			}

			inla_read_priorN(mb, ini, sec, &(ds->data_nprior[i]), "GAUSSIAN-std", i, NULL);

			if (!ds->data_nfixed[i]) {
				mb->theta = Realloc(mb->theta, mb->ntheta + 1, double **);
				mb->theta_hyperid = Realloc(mb->theta_hyperid, mb->ntheta + 1, char *);
				mb->theta_hyperid[mb->ntheta] = ds->data_nprior[i].hyperid;
				mb->theta_tag = Realloc(mb->theta_tag, mb->ntheta + 1, char *);
				mb->theta_tag_userscale = Realloc(mb->theta_tag_userscale, mb->ntheta + 1, char *);
				mb->theta_dir = Realloc(mb->theta_dir, mb->ntheta + 1, char *);

				GMRFLib_sprintf(&ctmp, "beta%1d for 0poisson%1s observations", i + 1, suff);
				mb->theta_tag[mb->ntheta] = inla_make_tag(ctmp, mb->ds);
				mb->theta_tag_userscale[mb->ntheta] = inla_make_tag(ctmp, mb->ds);
				GMRFLib_sprintf(&msg, "%s-parameter%1d", secname, i);
				mb->theta_dir[mb->ntheta] = msg;

				mb->theta_from = Realloc(mb->theta_from, mb->ntheta + 1, char *);
				mb->theta_to = Realloc(mb->theta_to, mb->ntheta + 1, char *);
				mb->theta_from[mb->ntheta] = Strdup(ds->data_nprior[i].from_theta);
				mb->theta_to[mb->ntheta] = Strdup(ds->data_nprior[i].to_theta);

				mb->theta[mb->ntheta] = ds->data_observations.poisson0_beta[i];
				mb->theta_map = Realloc(mb->theta_map, mb->ntheta + 1, map_func_tp *);

				mb->theta_map[mb->ntheta] = map_identity;
				mb->theta_map_arg = Realloc(mb->theta_map_arg, mb->ntheta + 1, void *);
				mb->theta_map_arg[mb->ntheta] = NULL;
				mb->ntheta++;
				ds->data_ntheta++;
			}
		}
	}
		break;

	case L_0BINOMIAL:
	case L_0BINOMIALS:
	{
		for (i = 0; i < BINOMIAL0_MAXTHETA; i++) {
			GMRFLib_sprintf(&ctmp, "FIXED%1d", i);
			iniparser_getstring(ini, inla_string_join(secname, ctmp), NULL);

			GMRFLib_sprintf(&ctmp, "INITIAL%1d", i);
			iniparser_getstring(ini, inla_string_join(secname, ctmp), NULL);

			GMRFLib_sprintf(&ctmp, "PRIOR%1d", i);
			iniparser_getstring(ini, inla_string_join(secname, ctmp), NULL);

			GMRFLib_sprintf(&ctmp, "HYPERID%1d", i);
			iniparser_getstring(ini, inla_string_join(secname, ctmp), NULL);

			GMRFLib_sprintf(&ctmp, "PARAMETERS%1d", i);
			iniparser_getstring(ini, inla_string_join(secname, ctmp), NULL);

			GMRFLib_sprintf(&ctmp, "to.theta%1d", i);
			iniparser_getstring(ini, inla_string_join(secname, ctmp), NULL);

			GMRFLib_sprintf(&ctmp, "from.theta%1d", i);
			iniparser_getstring(ini, inla_string_join(secname, ctmp), NULL);
		}

		char *link_simple = iniparser_getstring(ini, inla_string_join(secname, "LINK.SIMPLE"), NULL);
		ds->data_observations.link_simple_name = link_simple;
		if (!strcasecmp(link_simple, "IDENTITY")) {
			ds->data_observations.link_simple_invlinkfunc = link_identity;
		} else if (!strcasecmp(link_simple, "LOG")) {
			ds->data_observations.link_simple_invlinkfunc = link_log;
		} else if (!strcasecmp(link_simple, "PROBIT")) {
			ds->data_observations.link_simple_invlinkfunc = link_probit;
		} else if (!strcasecmp(link_simple, "CLOGLOG")) {
			ds->data_observations.link_simple_invlinkfunc = link_cloglog;
		} else if (!strcasecmp(link_simple, "LOGIT")) {
			ds->data_observations.link_simple_invlinkfunc = link_logit;
		} else {
			GMRFLib_sprintf(&msg, "%s: 0binomial(S) likelihood: no valid link.simple[%s]", secname, link_simple);
			inla_error_general(msg);
			exit(1);
		}
		if (mb->verbose) {
			printf("\t\tlink.simple[%s]\n", ds->data_observations.link_simple_name);
		}

		const char *suff = Strdup((ds->data_id == L_0BINOMIAL ? "" : "S"));
		ds->data_nfixed = Calloc(BINOMIAL0_MAXTHETA + 1, int);
		ds->data_nprior = Calloc(BINOMIAL0_MAXTHETA, Prior_tp);
		ds->data_observations.binomial0_beta = Calloc(BINOMIAL0_MAXTHETA, double **);

		for (i = 0; i < ds->data_observations.binomial0_nbeta; i++) {
			GMRFLib_sprintf(&ctmp, "INITIAL%1d", i);
			tmp = iniparser_getdouble(ini, inla_string_join(secname, ctmp), 0.0);	/* YES! */

			GMRFLib_sprintf(&ctmp, "FIXED%1d", i);
			ds->data_nfixed[i] = iniparser_getboolean(ini, inla_string_join(secname, ctmp), 0);
			if (!ds->data_nfixed[i] && mb->mode_use_mode) {
				tmp = mb->theta_file[mb->theta_counter_file++];
				if (mb->mode_fixed)
					ds->data_nfixed[i] = 1;
			}

			HYPER_NEW(ds->data_observations.binomial0_beta[i], tmp);
			if (mb->verbose) {
				printf("\t\tbeta[%1d] = %g\n", i, ds->data_observations.binomial0_beta[i][0][0]);
				printf("\t\tfixed[%1d] = %1d\n", i, ds->data_nfixed[i]);
			}

			inla_read_priorN(mb, ini, sec, &(ds->data_nprior[i]), "GAUSSIAN-std", i, NULL);

			if (!ds->data_nfixed[i]) {
				mb->theta = Realloc(mb->theta, mb->ntheta + 1, double **);
				mb->theta_hyperid = Realloc(mb->theta_hyperid, mb->ntheta + 1, char *);
				mb->theta_hyperid[mb->ntheta] = ds->data_nprior[i].hyperid;
				mb->theta_tag = Realloc(mb->theta_tag, mb->ntheta + 1, char *);
				mb->theta_tag_userscale = Realloc(mb->theta_tag_userscale, mb->ntheta + 1, char *);
				mb->theta_dir = Realloc(mb->theta_dir, mb->ntheta + 1, char *);

				GMRFLib_sprintf(&ctmp, "beta%1d for 0binomial%1s observations", i + 1, suff);
				mb->theta_tag[mb->ntheta] = inla_make_tag(ctmp, mb->ds);
				mb->theta_tag_userscale[mb->ntheta] = inla_make_tag(ctmp, mb->ds);
				GMRFLib_sprintf(&msg, "%s-parameter%1d", secname, i);
				mb->theta_dir[mb->ntheta] = msg;

				mb->theta_from = Realloc(mb->theta_from, mb->ntheta + 1, char *);
				mb->theta_to = Realloc(mb->theta_to, mb->ntheta + 1, char *);
				mb->theta_from[mb->ntheta] = Strdup(ds->data_nprior[i].from_theta);
				mb->theta_to[mb->ntheta] = Strdup(ds->data_nprior[i].to_theta);

				mb->theta[mb->ntheta] = ds->data_observations.binomial0_beta[i];
				mb->theta_map = Realloc(mb->theta_map, mb->ntheta + 1, map_func_tp *);

				mb->theta_map[mb->ntheta] = map_identity;
				mb->theta_map_arg = Realloc(mb->theta_map_arg, mb->ntheta + 1, void *);
				mb->theta_map_arg[mb->ntheta] = NULL;
				mb->ntheta++;
				ds->data_ntheta++;
			}
		}
	}
		break;

	case L_BINOMIALMIX:
	{
		const int nbeta = BINOMIALMIX_NBETA;
		for (i = 0; i < nbeta; i++) {
			GMRFLib_sprintf(&ctmp, "FIXED%1d", i);
			iniparser_getstring(ini, inla_string_join(secname, ctmp), NULL);

			GMRFLib_sprintf(&ctmp, "INITIAL%1d", i);
			iniparser_getstring(ini, inla_string_join(secname, ctmp), NULL);

			GMRFLib_sprintf(&ctmp, "PRIOR%1d", i);
			iniparser_getstring(ini, inla_string_join(secname, ctmp), NULL);

			GMRFLib_sprintf(&ctmp, "HYPERID%1d", i);
			iniparser_getstring(ini, inla_string_join(secname, ctmp), NULL);

			GMRFLib_sprintf(&ctmp, "PARAMETERS%1d", i);
			iniparser_getstring(ini, inla_string_join(secname, ctmp), NULL);

			GMRFLib_sprintf(&ctmp, "to.theta%1d", i);
			iniparser_getstring(ini, inla_string_join(secname, ctmp), NULL);

			GMRFLib_sprintf(&ctmp, "from.theta%1d", i);
			iniparser_getstring(ini, inla_string_join(secname, ctmp), NULL);
		}

		ds->data_nfixed = Calloc(nbeta, int);
		ds->data_nprior = Calloc(nbeta, Prior_tp);
		ds->data_observations.binmix_beta = Calloc(nbeta, double **);

		for (i = 0; i < nbeta; i++) {
			GMRFLib_sprintf(&ctmp, "INITIAL%1d", i);
			tmp = iniparser_getdouble(ini, inla_string_join(secname, ctmp), 0.0);	/* YES! */

			GMRFLib_sprintf(&ctmp, "FIXED%1d", i);
			ds->data_nfixed[i] = iniparser_getboolean(ini, inla_string_join(secname, ctmp), 0);
			if (!ds->data_nfixed[i] && mb->mode_use_mode) {
				tmp = mb->theta_file[mb->theta_counter_file++];
				if (mb->mode_fixed)
					ds->data_nfixed[i] = 1;
			}

			HYPER_NEW(ds->data_observations.binmix_beta[i], tmp);
			if (mb->verbose) {
				printf("\t\tbeta[%1d] = %g\n", i, ds->data_observations.binmix_beta[i][0][0]);
				printf("\t\tfixed[%1d] = %1d\n", i, ds->data_nfixed[i]);
			}

			inla_read_priorN(mb, ini, sec, &(ds->data_nprior[i]), "GAUSSIAN-std", i, NULL);

			if (!ds->data_nfixed[i]) {
				mb->theta = Realloc(mb->theta, mb->ntheta + 1, double **);
				mb->theta_hyperid = Realloc(mb->theta_hyperid, mb->ntheta + 1, char *);
				mb->theta_hyperid[mb->ntheta] = ds->data_nprior[i].hyperid;
				mb->theta_tag = Realloc(mb->theta_tag, mb->ntheta + 1, char *);
				mb->theta_tag_userscale = Realloc(mb->theta_tag_userscale, mb->ntheta + 1, char *);
				mb->theta_dir = Realloc(mb->theta_dir, mb->ntheta + 1, char *);

				GMRFLib_sprintf(&ctmp, "beta%1d for binomialmix observations", i + 1);
				mb->theta_tag[mb->ntheta] = inla_make_tag(ctmp, mb->ds);
				mb->theta_tag_userscale[mb->ntheta] = inla_make_tag(ctmp, mb->ds);
				GMRFLib_sprintf(&msg, "%s-parameter%1d", secname, i);
				mb->theta_dir[mb->ntheta] = msg;

				mb->theta_from = Realloc(mb->theta_from, mb->ntheta + 1, char *);
				mb->theta_to = Realloc(mb->theta_to, mb->ntheta + 1, char *);
				mb->theta_from[mb->ntheta] = Strdup(ds->data_nprior[i].from_theta);
				mb->theta_to[mb->ntheta] = Strdup(ds->data_nprior[i].to_theta);

				mb->theta[mb->ntheta] = ds->data_observations.binmix_beta[i];
				mb->theta_map = Realloc(mb->theta_map, mb->ntheta + 1, map_func_tp *);

				mb->theta_map[mb->ntheta] = map_identity;
				mb->theta_map_arg = Realloc(mb->theta_map_arg, mb->ntheta + 1, void *);
				mb->theta_map_arg[mb->ntheta] = NULL;
				mb->ntheta++;
				ds->data_ntheta++;
			}
		}
	}
		break;

	case L_GAMMA:
	case L_MGAMMA:
	{
		/*
		 * get options related to the gamma
		 */
		char *nm = (ds->data_id == L_GAMMA ? strdup("Gamma") : strdup("mGamma"));
		tmp = iniparser_getdouble(ini, inla_string_join(secname, "INITIAL"), 0.0);
		ds->data_fixed = iniparser_getboolean(ini, inla_string_join(secname, "FIXED"), 0);
		if (!ds->data_fixed && mb->mode_use_mode) {
			tmp = mb->theta_file[mb->theta_counter_file++];
			if (mb->mode_fixed)
				ds->data_fixed = 1;
		}
		HYPER_NEW(ds->data_observations.gamma_log_prec, tmp);
		if (mb->verbose) {
			printf("\t\tinitialise precision_intern[%g]\n", ds->data_observations.gamma_log_prec[0][0]);
			printf("\t\tfixed=[%1d]\n", ds->data_fixed);
		}
		inla_read_prior(mb, ini, sec, &(ds->data_prior), "LOGGAMMA", NULL);

		/*
		 * add theta 
		 */
		if (!ds->data_fixed) {
			mb->theta = Realloc(mb->theta, mb->ntheta + 1, double **);
			mb->theta_hyperid = Realloc(mb->theta_hyperid, mb->ntheta + 1, char *);
			mb->theta_hyperid[mb->ntheta] = ds->data_prior.hyperid;
			mb->theta_tag = Realloc(mb->theta_tag, mb->ntheta + 1, char *);
			mb->theta_tag_userscale = Realloc(mb->theta_tag_userscale, mb->ntheta + 1, char *);
			mb->theta_dir = Realloc(mb->theta_dir, mb->ntheta + 1, char *);
			char *nnm = NULL;
			GMRFLib_sprintf(&nnm, "Intern precision-parameter for the %s observations", nm);
			mb->theta_tag[mb->ntheta] = inla_make_tag(nnm, mb->ds);
			Free(nnm);
			GMRFLib_sprintf(&nnm, "Precision-parameter for the %s observations", nm);
			mb->theta_tag_userscale[mb->ntheta] = inla_make_tag(nnm, mb->ds);
			Free(nnm);

			GMRFLib_sprintf(&msg, "%s-parameter", secname);
			mb->theta_dir[mb->ntheta] = msg;

			mb->theta_from = Realloc(mb->theta_from, mb->ntheta + 1, char *);
			mb->theta_to = Realloc(mb->theta_to, mb->ntheta + 1, char *);
			mb->theta_from[mb->ntheta] = Strdup(ds->data_prior.from_theta);
			mb->theta_to[mb->ntheta] = Strdup(ds->data_prior.to_theta);

			mb->theta[mb->ntheta] = ds->data_observations.gamma_log_prec;
			mb->theta_map = Realloc(mb->theta_map, mb->ntheta + 1, map_func_tp *);
			mb->theta_map[mb->ntheta] = map_exp;
			mb->theta_map_arg = Realloc(mb->theta_map_arg, mb->ntheta + 1, void *);
			mb->theta_map_arg[mb->ntheta] = NULL;
			mb->ntheta++;
			ds->data_ntheta++;
		}
		Free(nm);
	}
		break;

	case L_GAMMASURV:
	case L_MGAMMASURV:
	{
		/*
		 * get options related to the gamma
		 */
		char *nm = (ds->data_id == L_GAMMASURV ? strdup("Gamma") : strdup("mGamma"));

		for (i = 0; i < CURE_MAXTHETA + 1; i++) {
			GMRFLib_sprintf(&ctmp, "FIXED%1d", i);
			iniparser_getstring(ini, inla_string_join(secname, ctmp), NULL);

			GMRFLib_sprintf(&ctmp, "INITIAL%1d", i);
			iniparser_getstring(ini, inla_string_join(secname, ctmp), NULL);

			GMRFLib_sprintf(&ctmp, "PRIOR%1d", i);
			iniparser_getstring(ini, inla_string_join(secname, ctmp), NULL);

			GMRFLib_sprintf(&ctmp, "HYPERID%1d", i);
			iniparser_getstring(ini, inla_string_join(secname, ctmp), NULL);

			GMRFLib_sprintf(&ctmp, "PARAMETERS%1d", i);
			iniparser_getstring(ini, inla_string_join(secname, ctmp), NULL);

			GMRFLib_sprintf(&ctmp, "to.theta%1d", i);
			iniparser_getstring(ini, inla_string_join(secname, ctmp), NULL);

			GMRFLib_sprintf(&ctmp, "from.theta%1d", i);
			iniparser_getstring(ini, inla_string_join(secname, ctmp), NULL);
		}

		ds->data_nfixed = Calloc(CURE_MAXTHETA + 1, int);
		ds->data_nprior = Calloc(CURE_MAXTHETA + 1, Prior_tp);
		ds->data_observations.cure_beta = Calloc(CURE_MAXTHETA, double **);

		tmp = iniparser_getdouble(ini, inla_string_join(secname, "INITIAL0"), 0.0);
		ds->data_nfixed[0] = iniparser_getboolean(ini, inla_string_join(secname, "FIXED0"), 0);
		if (!ds->data_nfixed[0] && mb->mode_use_mode) {
			tmp = mb->theta_file[mb->theta_counter_file++];
			if (mb->mode_fixed)
				ds->data_nfixed[0] = 1;
		}
		HYPER_NEW(ds->data_observations.gamma_log_prec, tmp);
		if (mb->verbose) {
			printf("\t\tinitialise precision_intern[%g]\n", ds->data_observations.gamma_log_prec[0][0]);
			printf("\t\tfixed=[%1d]\n", ds->data_fixed);
		}
		inla_read_priorN(mb, ini, sec, &(ds->data_nprior[0]), "LOGGAMMA", 0, NULL);

		if (!ds->data_nfixed[0]) {
			mb->theta = Realloc(mb->theta, mb->ntheta + 1, double **);
			mb->theta_hyperid = Realloc(mb->theta_hyperid, mb->ntheta + 1, char *);
			mb->theta_hyperid[mb->ntheta] = ds->data_nprior[0].hyperid;
			mb->theta_tag = Realloc(mb->theta_tag, mb->ntheta + 1, char *);
			mb->theta_tag_userscale = Realloc(mb->theta_tag_userscale, mb->ntheta + 1, char *);
			mb->theta_dir = Realloc(mb->theta_dir, mb->ntheta + 1, char *);

			char *nnm = NULL;
			GMRFLib_sprintf(&nnm, "Intern precision-parameter for the %s surv", nm);
			mb->theta_tag[mb->ntheta] = inla_make_tag(nnm, mb->ds);
			Free(nnm);

			GMRFLib_sprintf(&nnm, "Precision-parameter for the %s surv", nm);
			mb->theta_tag_userscale[mb->ntheta] = inla_make_tag(nnm, mb->ds);
			Free(nnm);

			GMRFLib_sprintf(&msg, "%s-parameter", secname);
			mb->theta_dir[mb->ntheta] = msg;

			mb->theta_from = Realloc(mb->theta_from, mb->ntheta + 1, char *);
			mb->theta_to = Realloc(mb->theta_to, mb->ntheta + 1, char *);
			mb->theta_from[mb->ntheta] = Strdup(ds->data_nprior[0].from_theta);
			mb->theta_to[mb->ntheta] = Strdup(ds->data_nprior[0].to_theta);

			mb->theta[mb->ntheta] = ds->data_observations.gamma_log_prec;
			mb->theta_map = Realloc(mb->theta_map, mb->ntheta + 1, map_func_tp *);
			mb->theta_map[mb->ntheta] = map_exp;
			mb->theta_map_arg = Realloc(mb->theta_map_arg, mb->ntheta + 1, void *);
			mb->theta_map_arg[mb->ntheta] = NULL;
			mb->ntheta++;
			ds->data_ntheta++;
		}

		for (i = 1; i < ds->data_observations.cure_ncov + 1; i++) {
			GMRFLib_sprintf(&ctmp, "INITIAL%1d", i);
			tmp = iniparser_getdouble(ini, inla_string_join(secname, ctmp), 0.0);	/* YES! */

			GMRFLib_sprintf(&ctmp, "FIXED%1d", i);
			ds->data_nfixed[i] = iniparser_getboolean(ini, inla_string_join(secname, ctmp), 0);
			if (!ds->data_nfixed[i] && mb->mode_use_mode) {
				tmp = mb->theta_file[mb->theta_counter_file++];
				if (mb->mode_fixed)
					ds->data_nfixed[i] = 1;
			}

			HYPER_NEW(ds->data_observations.cure_beta[i - 1], tmp);
			if (mb->verbose) {
				printf("\t\tbeta[%1d] = %g\n", i, ds->data_observations.cure_beta[i - 1][0][0]);
				printf("\t\tfixed[%1d] = %1d\n", i, ds->data_nfixed[i]);
			}
			inla_read_priorN(mb, ini, sec, &(ds->data_nprior[i]), "GAUSSIAN-std", i, NULL);

			if (!ds->data_nfixed[i]) {
				mb->theta = Realloc(mb->theta, mb->ntheta + 1, double **);
				mb->theta_hyperid = Realloc(mb->theta_hyperid, mb->ntheta + 1, char *);
				mb->theta_hyperid[mb->ntheta] = ds->data_nprior[i].hyperid;
				mb->theta_tag = Realloc(mb->theta_tag, mb->ntheta + 1, char *);
				mb->theta_tag_userscale = Realloc(mb->theta_tag_userscale, mb->ntheta + 1, char *);
				mb->theta_dir = Realloc(mb->theta_dir, mb->ntheta + 1, char *);

				if (ds->data_id == L_GAMMASURV) {
					GMRFLib_sprintf(&ctmp, "beta%1d for Gamma-Cure", i);
				} else {
					GMRFLib_sprintf(&ctmp, "beta%1d for mGamma-Cure", i);
				}
				mb->theta_tag[mb->ntheta] = inla_make_tag(ctmp, mb->ds);
				mb->theta_tag_userscale[mb->ntheta] = inla_make_tag(ctmp, mb->ds);
				GMRFLib_sprintf(&msg, "%s-parameter%1d", secname, i);
				mb->theta_dir[mb->ntheta] = msg;

				mb->theta_from = Realloc(mb->theta_from, mb->ntheta + 1, char *);
				mb->theta_to = Realloc(mb->theta_to, mb->ntheta + 1, char *);
				mb->theta_from[mb->ntheta] = Strdup(ds->data_nprior[i].from_theta);
				mb->theta_to[mb->ntheta] = Strdup(ds->data_nprior[i].to_theta);

				mb->theta[mb->ntheta] = ds->data_observations.cure_beta[i - 1];
				mb->theta_map = Realloc(mb->theta_map, mb->ntheta + 1, map_func_tp *);

				mb->theta_map[mb->ntheta] = map_identity;
				mb->theta_map_arg = Realloc(mb->theta_map_arg, mb->ntheta + 1, void *);
				mb->theta_map_arg[mb->ntheta] = NULL;

				mb->ntheta++;
				ds->data_ntheta++;
			}
		}
		// all the remaining ones are fixed
		for (i = 1 + ds->data_observations.cure_ncov; i < 1 + CURE_MAXTHETA; i++) {
			ds->data_nfixed[i] = 1;
		}
		Free(nm);
	}
		break;

	case L_GAMMAJW:
	{
		// nothing
	}
		break;

	case L_GAMMAJWSURV:
	{
		/*
		 * get options related to the gammajw
		 */
		for (i = 0; i < CURE_MAXTHETA; i++) {
			GMRFLib_sprintf(&ctmp, "FIXED%1d", i);
			iniparser_getstring(ini, inla_string_join(secname, ctmp), NULL);

			GMRFLib_sprintf(&ctmp, "INITIAL%1d", i);
			iniparser_getstring(ini, inla_string_join(secname, ctmp), NULL);

			GMRFLib_sprintf(&ctmp, "PRIOR%1d", i);
			iniparser_getstring(ini, inla_string_join(secname, ctmp), NULL);

			GMRFLib_sprintf(&ctmp, "HYPERID%1d", i);
			iniparser_getstring(ini, inla_string_join(secname, ctmp), NULL);

			GMRFLib_sprintf(&ctmp, "PARAMETERS%1d", i);
			iniparser_getstring(ini, inla_string_join(secname, ctmp), NULL);

			GMRFLib_sprintf(&ctmp, "to.theta%1d", i);
			iniparser_getstring(ini, inla_string_join(secname, ctmp), NULL);

			GMRFLib_sprintf(&ctmp, "from.theta%1d", i);
			iniparser_getstring(ini, inla_string_join(secname, ctmp), NULL);
		}

		ds->data_nfixed = Calloc(CURE_MAXTHETA, int);
		ds->data_nprior = Calloc(CURE_MAXTHETA, Prior_tp);
		ds->data_observations.cure_beta = Calloc(CURE_MAXTHETA, double **);

		for (i = 0; i < ds->data_observations.cure_ncov; i++) {
			GMRFLib_sprintf(&ctmp, "INITIAL%1d", i);
			tmp = iniparser_getdouble(ini, inla_string_join(secname, ctmp), 0.0);	/* YES! */

			GMRFLib_sprintf(&ctmp, "FIXED%1d", i);
			ds->data_nfixed[i] = iniparser_getboolean(ini, inla_string_join(secname, ctmp), 0);
			if (!ds->data_nfixed[i] && mb->mode_use_mode) {
				tmp = mb->theta_file[mb->theta_counter_file++];
				if (mb->mode_fixed)
					ds->data_nfixed[i] = 1;
			}

			HYPER_NEW(ds->data_observations.cure_beta[i], tmp);
			if (mb->verbose) {
				printf("\t\tbeta[%1d] = %g\n", i, ds->data_observations.cure_beta[i][0][0]);
				printf("\t\tfixed[%1d] = %1d\n", i, ds->data_nfixed[i]);
			}
			inla_read_priorN(mb, ini, sec, &(ds->data_nprior[i]), "GAUSSIAN-std", i, NULL);

			if (!ds->data_nfixed[i]) {
				mb->theta = Realloc(mb->theta, mb->ntheta + 1, double **);
				mb->theta_hyperid = Realloc(mb->theta_hyperid, mb->ntheta + 1, char *);
				mb->theta_hyperid[mb->ntheta] = ds->data_nprior[i].hyperid;
				mb->theta_tag = Realloc(mb->theta_tag, mb->ntheta + 1, char *);
				mb->theta_tag_userscale = Realloc(mb->theta_tag_userscale, mb->ntheta + 1, char *);
				mb->theta_dir = Realloc(mb->theta_dir, mb->ntheta + 1, char *);
				GMRFLib_sprintf(&ctmp, "beta%1d for GammaJW-Cure", i + 1);

				mb->theta_tag[mb->ntheta] = inla_make_tag(ctmp, mb->ds);
				mb->theta_tag_userscale[mb->ntheta] = inla_make_tag(ctmp, mb->ds);
				GMRFLib_sprintf(&msg, "%s-parameter%1d", secname, i);
				mb->theta_dir[mb->ntheta] = msg;

				mb->theta_from = Realloc(mb->theta_from, mb->ntheta + 1, char *);
				mb->theta_to = Realloc(mb->theta_to, mb->ntheta + 1, char *);
				mb->theta_from[mb->ntheta] = Strdup(ds->data_nprior[i].from_theta);
				mb->theta_to[mb->ntheta] = Strdup(ds->data_nprior[i].to_theta);

				mb->theta[mb->ntheta] = ds->data_observations.cure_beta[i];
				mb->theta_map = Realloc(mb->theta_map, mb->ntheta + 1, map_func_tp *);

				mb->theta_map[mb->ntheta] = map_identity;
				mb->theta_map_arg = Realloc(mb->theta_map_arg, mb->ntheta + 1, void *);
				mb->theta_map_arg[mb->ntheta] = NULL;

				mb->ntheta++;
				ds->data_ntheta++;
			}
		}
		// all the remaining ones are fixed
		for (i = ds->data_observations.cure_ncov; i < CURE_MAXTHETA; i++) {
			ds->data_nfixed[i] = 1;
		}
	}
		break;

	case L_GAMMACOUNT:
	{
		/*
		 * get options related to the gammacount
		 */
		tmp = iniparser_getdouble(ini, inla_string_join(secname, "INITIAL"), 0.0);
		ds->data_fixed = iniparser_getboolean(ini, inla_string_join(secname, "FIXED"), 0);
		if (!ds->data_fixed && mb->mode_use_mode) {
			tmp = mb->theta_file[mb->theta_counter_file++];
			if (mb->mode_fixed)
				ds->data_fixed = 1;
		}
		HYPER_NEW(ds->data_observations.gammacount_log_alpha, tmp);
		if (mb->verbose) {
			printf("\t\tinitialise log_alpha[%g]\n", ds->data_observations.gammacount_log_alpha[0][0]);
			printf("\t\tfixed=[%1d]\n", ds->data_fixed);
		}
		inla_read_prior(mb, ini, sec, &(ds->data_prior), "LOGGAMMA", NULL);

		/*
		 * add theta 
		 */
		if (!ds->data_fixed) {
			mb->theta = Realloc(mb->theta, mb->ntheta + 1, double **);
			mb->theta_hyperid = Realloc(mb->theta_hyperid, mb->ntheta + 1, char *);
			mb->theta_hyperid[mb->ntheta] = ds->data_prior.hyperid;
			mb->theta_tag = Realloc(mb->theta_tag, mb->ntheta + 1, char *);
			mb->theta_tag_userscale = Realloc(mb->theta_tag_userscale, mb->ntheta + 1, char *);
			mb->theta_dir = Realloc(mb->theta_dir, mb->ntheta + 1, char *);
			mb->theta_tag[mb->ntheta] = inla_make_tag("Log-alpha parameter for Gammacount observations", mb->ds);
			mb->theta_tag_userscale[mb->ntheta] = inla_make_tag("Alpha parameter for Gammacount observations", mb->ds);
			GMRFLib_sprintf(&msg, "%s-parameter", secname);
			mb->theta_dir[mb->ntheta] = msg;

			mb->theta_from = Realloc(mb->theta_from, mb->ntheta + 1, char *);
			mb->theta_to = Realloc(mb->theta_to, mb->ntheta + 1, char *);
			mb->theta_from[mb->ntheta] = Strdup(ds->data_prior.from_theta);
			mb->theta_to[mb->ntheta] = Strdup(ds->data_prior.to_theta);

			mb->theta[mb->ntheta] = ds->data_observations.gammacount_log_alpha;
			mb->theta_map = Realloc(mb->theta_map, mb->ntheta + 1, map_func_tp *);
			mb->theta_map[mb->ntheta] = map_exp;
			mb->theta_map_arg = Realloc(mb->theta_map_arg, mb->ntheta + 1, void *);
			mb->theta_map_arg[mb->ntheta] = NULL;
			mb->ntheta++;
			ds->data_ntheta++;
		}
	}
		break;

	case L_QKUMAR:
	{
		/*
		 * get options related to the qkumar-distribution
		 */
		GMRFLib_ASSERT(ds->data_observations.quantile > 0.0 && ds->data_observations.quantile < 1.0, GMRFLib_EPARAMETER);

		tmp = iniparser_getdouble(ini, inla_string_join(secname, "INITIAL"), 0.0);	/* yes! */
		ds->data_fixed = iniparser_getboolean(ini, inla_string_join(secname, "FIXED"), 0);
		if (!ds->data_fixed && mb->mode_use_mode) {
			tmp = mb->theta_file[mb->theta_counter_file++];
			if (mb->mode_fixed)
				ds->data_fixed = 1;
		}
		HYPER_NEW(ds->data_observations.qkumar_log_prec, tmp);
		if (mb->verbose) {
			printf("\t\tinitialise log_precision[%g]\n", ds->data_observations.qkumar_log_prec[0][0]);
			printf("\t\tfixed=[%1d]\n", ds->data_fixed);
		}
		inla_read_prior(mb, ini, sec, &(ds->data_prior), "LOGGAMMA", NULL);

		/*
		 * add theta 
		 */
		if (!ds->data_fixed) {
			mb->theta = Realloc(mb->theta, mb->ntheta + 1, double **);
			mb->theta_hyperid = Realloc(mb->theta_hyperid, mb->ntheta + 1, char *);
			mb->theta_hyperid[mb->ntheta] = ds->data_prior.hyperid;
			mb->theta_tag = Realloc(mb->theta_tag, mb->ntheta + 1, char *);
			mb->theta_tag_userscale = Realloc(mb->theta_tag_userscale, mb->ntheta + 1, char *);
			mb->theta_dir = Realloc(mb->theta_dir, mb->ntheta + 1, char *);
			mb->theta_tag[mb->ntheta] = inla_make_tag("log precision for qkumar observations", mb->ds);
			mb->theta_tag_userscale[mb->ntheta] = inla_make_tag("precision for qkumar observations", mb->ds);
			GMRFLib_sprintf(&msg, "%s-parameter", secname);
			mb->theta_dir[mb->ntheta] = msg;

			mb->theta_from = Realloc(mb->theta_from, mb->ntheta + 1, char *);
			mb->theta_to = Realloc(mb->theta_to, mb->ntheta + 1, char *);
			mb->theta_from[mb->ntheta] = Strdup(ds->data_prior.from_theta);
			mb->theta_to[mb->ntheta] = Strdup(ds->data_prior.to_theta);

			mb->theta[mb->ntheta] = ds->data_observations.qkumar_log_prec;
			mb->theta_map = Realloc(mb->theta_map, mb->ntheta + 1, map_func_tp *);
			mb->theta_map[mb->ntheta] = map_prec_qkumar;
			mb->theta_map_arg = Realloc(mb->theta_map_arg, mb->ntheta + 1, void *);
			mb->theta_map_arg[mb->ntheta] = NULL;
			mb->ntheta++;
			ds->data_ntheta++;
		}
	}
		break;

	case L_QLOGLOGISTIC:
	case L_LOGLOGISTIC:
	{

		if (ds->data_id == L_QLOGLOGISTIC) {
			GMRFLib_ASSERT(ds->data_observations.quantile > 0.0 && ds->data_observations.quantile < 1.0, GMRFLib_EPARAMETER);
		}
		assert(ds->variant == 0 || ds->variant == 1);
		tmp = iniparser_getdouble(ini, inla_string_join(secname, "INITIAL"), 0.0);	/* yes! */
		ds->data_fixed = iniparser_getboolean(ini, inla_string_join(secname, "FIXED"), 0);
		if (!ds->data_fixed && mb->mode_use_mode) {
			tmp = mb->theta_file[mb->theta_counter_file++];
			if (mb->mode_fixed)
				ds->data_fixed = 1;
		}
		HYPER_NEW(ds->data_observations.alpha_intern, tmp);
		if (mb->verbose) {
			printf("\t\tinitialise log(alpha) [%g]\n", ds->data_observations.alpha_intern[0][0]);
			printf("\t\tfixed=[%1d]\n", ds->data_fixed);
		}
		inla_read_prior(mb, ini, sec, &(ds->data_prior), "LOGGAMMA", NULL);
		/*
		 * add theta 
		 */
		if (!ds->data_fixed) {
			mb->theta = Realloc(mb->theta, mb->ntheta + 1, double **);
			mb->theta_hyperid = Realloc(mb->theta_hyperid, mb->ntheta + 1, char *);
			mb->theta_hyperid[mb->ntheta] = ds->data_prior.hyperid;
			mb->theta_tag = Realloc(mb->theta_tag, mb->ntheta + 1, char *);
			mb->theta_tag_userscale = Realloc(mb->theta_tag_userscale, mb->ntheta + 1, char *);
			mb->theta_dir = Realloc(mb->theta_dir, mb->ntheta + 1, char *);
			if (ds->data_id == L_QLOGLOGISTIC) {
				mb->theta_tag[mb->ntheta] = inla_make_tag("log alpha for qloglogistic observations", mb->ds);
				mb->theta_tag_userscale[mb->ntheta] = inla_make_tag("alpha for qloglogistic observations", mb->ds);
			} else {
				mb->theta_tag[mb->ntheta] = inla_make_tag("log alpha for loglogistic observations", mb->ds);
				mb->theta_tag_userscale[mb->ntheta] = inla_make_tag("alpha for loglogistic observations", mb->ds);
			}
			GMRFLib_sprintf(&msg, "%s-parameter", secname);
			mb->theta_dir[mb->ntheta] = msg;

			mb->theta_from = Realloc(mb->theta_from, mb->ntheta + 1, char *);
			mb->theta_to = Realloc(mb->theta_to, mb->ntheta + 1, char *);
			mb->theta_from[mb->ntheta] = Strdup(ds->data_prior.from_theta);
			mb->theta_to[mb->ntheta] = Strdup(ds->data_prior.to_theta);

			mb->theta[mb->ntheta] = ds->data_observations.alpha_intern;
			mb->theta_map = Realloc(mb->theta_map, mb->ntheta + 1, map_func_tp *);
			mb->theta_map[mb->ntheta] = map_exp;
			mb->theta_map_arg = Realloc(mb->theta_map_arg, mb->ntheta + 1, void *);
			mb->theta_map_arg[mb->ntheta] = NULL;
			mb->ntheta++;
			ds->data_ntheta++;
		}
	}
		break;

	case L_QLOGLOGISTICSURV:
	case L_LOGLOGISTICSURV:
	{
		for (i = 0; i < CURE_MAXTHETA + 1; i++) {
			GMRFLib_sprintf(&ctmp, "FIXED%1d", i);
			iniparser_getstring(ini, inla_string_join(secname, ctmp), NULL);

			GMRFLib_sprintf(&ctmp, "INITIAL%1d", i);
			iniparser_getstring(ini, inla_string_join(secname, ctmp), NULL);

			GMRFLib_sprintf(&ctmp, "PRIOR%1d", i);
			iniparser_getstring(ini, inla_string_join(secname, ctmp), NULL);

			GMRFLib_sprintf(&ctmp, "HYPERID%1d", i);
			iniparser_getstring(ini, inla_string_join(secname, ctmp), NULL);

			GMRFLib_sprintf(&ctmp, "PARAMETERS%1d", i);
			iniparser_getstring(ini, inla_string_join(secname, ctmp), NULL);

			GMRFLib_sprintf(&ctmp, "to.theta%1d", i);
			iniparser_getstring(ini, inla_string_join(secname, ctmp), NULL);

			GMRFLib_sprintf(&ctmp, "from.theta%1d", i);
			iniparser_getstring(ini, inla_string_join(secname, ctmp), NULL);
		}

		ds->data_nfixed = Calloc(CURE_MAXTHETA + 1, int);
		ds->data_nprior = Calloc(CURE_MAXTHETA + 1, Prior_tp);
		ds->data_observations.cure_beta = Calloc(CURE_MAXTHETA, double **);

		if (ds->data_id == L_QLOGLOGISTICSURV) {
			GMRFLib_ASSERT(ds->data_observations.quantile > 0.0 && ds->data_observations.quantile < 1.0, GMRFLib_EPARAMETER);
		}
		assert(ds->variant == 0 || ds->variant == 1);
		tmp = iniparser_getdouble(ini, inla_string_join(secname, "INITIAL0"), 0.0);	/* yes! */
		ds->data_nfixed[0] = iniparser_getboolean(ini, inla_string_join(secname, "FIXED0"), 0);
		if (!ds->data_nfixed[0] && mb->mode_use_mode) {
			tmp = mb->theta_file[mb->theta_counter_file++];
			if (mb->mode_fixed)
				ds->data_nfixed[0] = 1;
		}
		HYPER_NEW(ds->data_observations.alpha_intern, tmp);
		if (mb->verbose) {
			printf("\t\tinitialise log(alpha) [%g]\n", ds->data_observations.alpha_intern[0][0]);
			printf("\t\tfixed=[%1d]\n", ds->data_fixed);
		}
		inla_read_priorN(mb, ini, sec, &(ds->data_nprior[0]), "LOGGAMMA", 0, NULL);
		/*
		 * add theta 
		 */
		if (!ds->data_nfixed[0]) {
			mb->theta = Realloc(mb->theta, mb->ntheta + 1, double **);
			mb->theta_hyperid = Realloc(mb->theta_hyperid, mb->ntheta + 1, char *);
			mb->theta_hyperid[mb->ntheta] = ds->data_nprior[0].hyperid;
			mb->theta_tag = Realloc(mb->theta_tag, mb->ntheta + 1, char *);
			mb->theta_tag_userscale = Realloc(mb->theta_tag_userscale, mb->ntheta + 1, char *);
			mb->theta_dir = Realloc(mb->theta_dir, mb->ntheta + 1, char *);
			if (ds->data_id == L_QLOGLOGISTICSURV) {
				mb->theta_tag[mb->ntheta] = inla_make_tag("log alpha for qloglogisticsurv observations", mb->ds);
				mb->theta_tag_userscale[mb->ntheta] = inla_make_tag("alpha for qloglogisticsurv observations", mb->ds);
			} else {
				mb->theta_tag[mb->ntheta] = inla_make_tag("log alpha for qloglogisticsurv observations", mb->ds);
				mb->theta_tag_userscale[mb->ntheta] = inla_make_tag("alpha for qloglogisticsurv observations", mb->ds);
			}
			GMRFLib_sprintf(&msg, "%s-parameter", secname);
			mb->theta_dir[mb->ntheta] = msg;

			mb->theta_from = Realloc(mb->theta_from, mb->ntheta + 1, char *);
			mb->theta_to = Realloc(mb->theta_to, mb->ntheta + 1, char *);
			mb->theta_from[mb->ntheta] = Strdup(ds->data_nprior[0].from_theta);
			mb->theta_to[mb->ntheta] = Strdup(ds->data_nprior[0].to_theta);

			mb->theta[mb->ntheta] = ds->data_observations.alpha_intern;
			mb->theta_map = Realloc(mb->theta_map, mb->ntheta + 1, map_func_tp *);
			mb->theta_map[mb->ntheta] = map_exp;
			mb->theta_map_arg = Realloc(mb->theta_map_arg, mb->ntheta + 1, void *);
			mb->theta_map_arg[mb->ntheta] = NULL;
			mb->ntheta++;
			ds->data_ntheta++;
		}

		for (i = 1; i < ds->data_observations.cure_ncov + 1; i++) {
			GMRFLib_sprintf(&ctmp, "INITIAL%1d", i);
			tmp = iniparser_getdouble(ini, inla_string_join(secname, ctmp), 0.0);	/* YES! */

			GMRFLib_sprintf(&ctmp, "FIXED%1d", i);
			ds->data_nfixed[i] = iniparser_getboolean(ini, inla_string_join(secname, ctmp), 0);
			if (!ds->data_nfixed[i] && mb->mode_use_mode) {
				tmp = mb->theta_file[mb->theta_counter_file++];
				if (mb->mode_fixed)
					ds->data_nfixed[i] = 1;
			}

			HYPER_NEW(ds->data_observations.cure_beta[i - 1], tmp);
			if (mb->verbose) {
				printf("\t\tbeta[%1d] = %g\n", i, ds->data_observations.cure_beta[i - 1][0][0]);
				printf("\t\tfixed[%1d] = %1d\n", i, ds->data_nfixed[i]);
			}
			inla_read_priorN(mb, ini, sec, &(ds->data_nprior[i]), "GAUSSIAN-std", i, NULL);

			if (!ds->data_nfixed[i]) {
				mb->theta = Realloc(mb->theta, mb->ntheta + 1, double **);
				mb->theta_hyperid = Realloc(mb->theta_hyperid, mb->ntheta + 1, char *);
				mb->theta_hyperid[mb->ntheta] = ds->data_nprior[i].hyperid;
				mb->theta_tag = Realloc(mb->theta_tag, mb->ntheta + 1, char *);
				mb->theta_tag_userscale = Realloc(mb->theta_tag_userscale, mb->ntheta + 1, char *);
				mb->theta_dir = Realloc(mb->theta_dir, mb->ntheta + 1, char *);
				if (ds->data_id == L_QLOGLOGISTICSURV) {
					GMRFLib_sprintf(&ctmp, "beta%1d for qlogLogistic-Cure", i);
				} else {
					GMRFLib_sprintf(&ctmp, "beta%1d for logLogistic-Cure", i);
				}
				mb->theta_tag[mb->ntheta] = inla_make_tag(ctmp, mb->ds);
				mb->theta_tag_userscale[mb->ntheta] = inla_make_tag(ctmp, mb->ds);
				GMRFLib_sprintf(&msg, "%s-parameter%1d", secname, i);
				mb->theta_dir[mb->ntheta] = msg;

				mb->theta_from = Realloc(mb->theta_from, mb->ntheta + 1, char *);
				mb->theta_to = Realloc(mb->theta_to, mb->ntheta + 1, char *);
				mb->theta_from[mb->ntheta] = Strdup(ds->data_nprior[i].from_theta);
				mb->theta_to[mb->ntheta] = Strdup(ds->data_nprior[i].to_theta);

				mb->theta[mb->ntheta] = ds->data_observations.cure_beta[i - 1];
				mb->theta_map = Realloc(mb->theta_map, mb->ntheta + 1, map_func_tp *);

				mb->theta_map[mb->ntheta] = map_identity;
				mb->theta_map_arg = Realloc(mb->theta_map_arg, mb->ntheta + 1, void *);
				mb->theta_map_arg[mb->ntheta] = NULL;

				mb->ntheta++;
				ds->data_ntheta++;
			}
		}
	}
		break;

	case L_BETA:
	{
		/*
		 * get options related to the beta
		 */
		ds->data_observations.beta_censor_value = iniparser_getdouble(ini, inla_string_join(secname, "BETA.CENSOR.VALUE"), 0.0);
		if (mb->verbose) {
			printf("\t\tbcensor.value [%g]\n", ds->data_observations.beta_censor_value);
		}

		if (beta_delayed_error && ds->data_observations.beta_censor_value == 0.0) {
			GMRFLib_sprintf(&msg,
					"%s: Beta data: %1d observations are either 0 or 1,\n\t\tbut then you need to enable censoring, see inla.doc('^beta$')\n",
					secname, beta_delayed_error);
			inla_error_general(msg);
		}

		tmp = iniparser_getdouble(ini, inla_string_join(secname, "INITIAL"), 0.0);
		ds->data_fixed = iniparser_getboolean(ini, inla_string_join(secname, "FIXED"), 0);
		if (!ds->data_fixed && mb->mode_use_mode) {
			tmp = mb->theta_file[mb->theta_counter_file++];
			if (mb->mode_fixed)
				ds->data_fixed = 1;
		}
		HYPER_NEW(ds->data_observations.beta_precision_intern, tmp);
		if (mb->verbose) {
			printf("\t\tinitialise precision_intern[%g]\n", ds->data_observations.beta_precision_intern[0][0]);
			printf("\t\tfixed=[%1d]\n", ds->data_fixed);
		}
		inla_read_prior(mb, ini, sec, &(ds->data_prior), "LOGGAMMA", NULL);

		/*
		 * add theta 
		 */
		if (!ds->data_fixed) {
			mb->theta = Realloc(mb->theta, mb->ntheta + 1, double **);
			mb->theta_hyperid = Realloc(mb->theta_hyperid, mb->ntheta + 1, char *);
			mb->theta_hyperid[mb->ntheta] = ds->data_prior.hyperid;
			mb->theta_tag = Realloc(mb->theta_tag, mb->ntheta + 1, char *);
			mb->theta_tag_userscale = Realloc(mb->theta_tag_userscale, mb->ntheta + 1, char *);
			mb->theta_dir = Realloc(mb->theta_dir, mb->ntheta + 1, char *);
			mb->theta_tag[mb->ntheta] = inla_make_tag("intern precision-parameter for the beta observations", mb->ds);
			mb->theta_tag_userscale[mb->ntheta] = inla_make_tag("precision parameter for the beta observations", mb->ds);
			GMRFLib_sprintf(&msg, "%s-parameter", secname);
			mb->theta_dir[mb->ntheta] = msg;

			mb->theta_from = Realloc(mb->theta_from, mb->ntheta + 1, char *);
			mb->theta_to = Realloc(mb->theta_to, mb->ntheta + 1, char *);
			mb->theta_from[mb->ntheta] = Strdup(ds->data_prior.from_theta);
			mb->theta_to[mb->ntheta] = Strdup(ds->data_prior.to_theta);

			mb->theta[mb->ntheta] = ds->data_observations.beta_precision_intern;
			mb->theta_map = Realloc(mb->theta_map, mb->ntheta + 1, map_func_tp *);
			mb->theta_map[mb->ntheta] = map_exp;
			mb->theta_map_arg = Realloc(mb->theta_map_arg, mb->ntheta + 1, void *);
			mb->theta_map_arg[mb->ntheta] = NULL;
			mb->ntheta++;
			ds->data_ntheta++;
		}
	}
		break;

	case L_BETABINOMIAL:
	{
		/*
		 * get options related to the betabinomial
		 */
		tmp = iniparser_getdouble(ini, inla_string_join(secname, "INITIAL"), 0.0);
		ds->data_fixed = iniparser_getboolean(ini, inla_string_join(secname, "FIXED"), 0);
		if (!ds->data_fixed && mb->mode_use_mode) {
			tmp = mb->theta_file[mb->theta_counter_file++];
			if (mb->mode_fixed)
				ds->data_fixed = 1;
		}
		HYPER_NEW(ds->data_observations.betabinomial_overdispersion_intern, tmp);
		if (mb->verbose) {
			printf("\t\tinitialise overdispersion_intern[%g]\n", ds->data_observations.betabinomial_overdispersion_intern[0][0]);
			printf("\t\tfixed=[%1d]\n", ds->data_fixed);
		}
		inla_read_prior(mb, ini, sec, &(ds->data_prior), "LOGGAMMA", NULL);

		/*
		 * add theta 
		 */
		if (!ds->data_fixed) {
			mb->theta = Realloc(mb->theta, mb->ntheta + 1, double **);
			mb->theta_hyperid = Realloc(mb->theta_hyperid, mb->ntheta + 1, char *);
			mb->theta_hyperid[mb->ntheta] = ds->data_prior.hyperid;
			mb->theta_tag = Realloc(mb->theta_tag, mb->ntheta + 1, char *);
			mb->theta_tag_userscale = Realloc(mb->theta_tag_userscale, mb->ntheta + 1, char *);
			mb->theta_dir = Realloc(mb->theta_dir, mb->ntheta + 1, char *);
			mb->theta_tag[mb->ntheta] = inla_make_tag("intern overdispersion for the betabinomial observations", mb->ds);
			mb->theta_tag_userscale[mb->ntheta] = inla_make_tag("overdispersion for the betabinomial observations", mb->ds);
			GMRFLib_sprintf(&msg, "%s-parameter", secname);
			mb->theta_dir[mb->ntheta] = msg;

			mb->theta_from = Realloc(mb->theta_from, mb->ntheta + 1, char *);
			mb->theta_to = Realloc(mb->theta_to, mb->ntheta + 1, char *);
			mb->theta_from[mb->ntheta] = Strdup(ds->data_prior.from_theta);
			mb->theta_to[mb->ntheta] = Strdup(ds->data_prior.to_theta);

			mb->theta[mb->ntheta] = ds->data_observations.betabinomial_overdispersion_intern;
			mb->theta_map = Realloc(mb->theta_map, mb->ntheta + 1, map_func_tp *);
			mb->theta_map[mb->ntheta] = map_probability;
			mb->theta_map_arg = Realloc(mb->theta_map_arg, mb->ntheta + 1, void *);
			mb->theta_map_arg[mb->ntheta] = NULL;
			mb->ntheta++;
			ds->data_ntheta++;
		}
	}
		break;

	case L_BETABINOMIALNA:
	{
		/*
		 * get options related to the betabinomialna
		 */
		tmp = iniparser_getdouble(ini, inla_string_join(secname, "INITIAL"), 0.0);
		ds->data_fixed = iniparser_getboolean(ini, inla_string_join(secname, "FIXED"), 0);
		if (!ds->data_fixed && mb->mode_use_mode) {
			tmp = mb->theta_file[mb->theta_counter_file++];
			if (mb->mode_fixed)
				ds->data_fixed = 1;
		}
		HYPER_NEW(ds->data_observations.betabinomial_overdispersion_intern, tmp);
		if (mb->verbose) {
			printf("\t\tinitialise overdispersion_intern[%g]\n", ds->data_observations.betabinomial_overdispersion_intern[0][0]);
			printf("\t\tfixed=[%1d]\n", ds->data_fixed);
		}
		inla_read_prior(mb, ini, sec, &(ds->data_prior), "LOGGAMMA", NULL);

		/*
		 * add theta 
		 */
		if (!ds->data_fixed) {
			mb->theta = Realloc(mb->theta, mb->ntheta + 1, double **);
			mb->theta_hyperid = Realloc(mb->theta_hyperid, mb->ntheta + 1, char *);
			mb->theta_hyperid[mb->ntheta] = ds->data_prior.hyperid;
			mb->theta_tag = Realloc(mb->theta_tag, mb->ntheta + 1, char *);
			mb->theta_tag_userscale = Realloc(mb->theta_tag_userscale, mb->ntheta + 1, char *);
			mb->theta_dir = Realloc(mb->theta_dir, mb->ntheta + 1, char *);
			mb->theta_tag[mb->ntheta] = inla_make_tag("intern overdispersion for the betabinomialna observations", mb->ds);
			mb->theta_tag_userscale[mb->ntheta] = inla_make_tag("overdispersion for the betabinomialna observations", mb->ds);
			GMRFLib_sprintf(&msg, "%s-parameter", secname);
			mb->theta_dir[mb->ntheta] = msg;

			mb->theta_from = Realloc(mb->theta_from, mb->ntheta + 1, char *);
			mb->theta_to = Realloc(mb->theta_to, mb->ntheta + 1, char *);
			mb->theta_from[mb->ntheta] = Strdup(ds->data_prior.from_theta);
			mb->theta_to[mb->ntheta] = Strdup(ds->data_prior.to_theta);

			mb->theta[mb->ntheta] = ds->data_observations.betabinomial_overdispersion_intern;
			mb->theta_map = Realloc(mb->theta_map, mb->ntheta + 1, map_func_tp *);
			mb->theta_map[mb->ntheta] = map_probability;
			mb->theta_map_arg = Realloc(mb->theta_map_arg, mb->ntheta + 1, void *);
			mb->theta_map_arg[mb->ntheta] = NULL;
			mb->ntheta++;
			ds->data_ntheta++;
		}
	}
		break;

	case L_NBINOMIAL:
	{
		/*
		 * get options related to the negative binomial
		 */
		tmp = iniparser_getdouble(ini, inla_string_join(secname, "INITIAL"), log(10.0));
		ds->data_fixed = iniparser_getboolean(ini, inla_string_join(secname, "FIXED"), 0);
		if (!ds->data_fixed && mb->mode_use_mode) {
			tmp = mb->theta_file[mb->theta_counter_file++];
			if (mb->mode_fixed)
				ds->data_fixed = 1;
		}
		HYPER_NEW(ds->data_observations.log_size, tmp);
		assert(ds->variant == 0 || ds->variant == 1 || ds->variant == 2);
		if (mb->verbose) {
			printf("\t\tinitialise log_size[%g]\n", ds->data_observations.log_size[0][0]);
			printf("\t\tfixed=[%1d]\n", ds->data_fixed);
			printf("\t\tuse parameterization variant=[%1d]; see doc for details\n", ds->variant);
		}
		inla_read_prior(mb, ini, sec, &(ds->data_prior), "LOGGAMMA", NULL);

		/*
		 * add theta 
		 */
		if (!ds->data_fixed) {
			mb->theta = Realloc(mb->theta, mb->ntheta + 1, double **);
			mb->theta_hyperid = Realloc(mb->theta_hyperid, mb->ntheta + 1, char *);
			mb->theta_hyperid[mb->ntheta] = ds->data_prior.hyperid;
			mb->theta_tag = Realloc(mb->theta_tag, mb->ntheta + 1, char *);
			mb->theta_tag_userscale = Realloc(mb->theta_tag_userscale, mb->ntheta + 1, char *);
			mb->theta_dir = Realloc(mb->theta_dir, mb->ntheta + 1, char *);
			if (ds->variant == 0 || ds->variant == 1 || ds->variant == 2) {
				mb->theta_tag[mb->ntheta] = inla_make_tag("log size for the nbinomial observations (1/overdispersion)", mb->ds);
				mb->theta_tag_userscale[mb->ntheta] =
				    inla_make_tag("size for the nbinomial observations (1/overdispersion)", mb->ds);
			} else {
				assert(0 == 1);
			}
			GMRFLib_sprintf(&msg, "%s-parameter", secname);
			mb->theta_dir[mb->ntheta] = msg;

			mb->theta_from = Realloc(mb->theta_from, mb->ntheta + 1, char *);
			mb->theta_to = Realloc(mb->theta_to, mb->ntheta + 1, char *);
			mb->theta_from[mb->ntheta] = Strdup(ds->data_prior.from_theta);
			mb->theta_to[mb->ntheta] = Strdup(ds->data_prior.to_theta);

			mb->theta[mb->ntheta] = ds->data_observations.log_size;
			mb->theta_map = Realloc(mb->theta_map, mb->ntheta + 1, map_func_tp *);
			mb->theta_map[mb->ntheta] = map_exp;
			mb->theta_map_arg = Realloc(mb->theta_map_arg, mb->ntheta + 1, void *);
			mb->theta_map_arg[mb->ntheta] = NULL;
			mb->ntheta++;
			ds->data_ntheta++;
		}
	}
		break;

	case L_CENNBINOMIAL2:
	{
		/*
		 * get options related to the cen negative binomial 2
		 */
		tmp = iniparser_getdouble(ini, inla_string_join(secname, "INITIAL"), log(10.0));
		ds->data_fixed = iniparser_getboolean(ini, inla_string_join(secname, "FIXED"), 0);
		if (!ds->data_fixed && mb->mode_use_mode) {
			tmp = mb->theta_file[mb->theta_counter_file++];
			if (mb->mode_fixed)
				ds->data_fixed = 1;
		}
		HYPER_NEW(ds->data_observations.log_size, tmp);
		assert(ds->variant == 0 || ds->variant == 1 || ds->variant == 2);
		if (mb->verbose) {
			printf("\t\tinitialise log_size[%g]\n", ds->data_observations.log_size[0][0]);
			printf("\t\tfixed=[%1d]\n", ds->data_fixed);
			printf("\t\tuse parameterization variant=[%1d]; see doc for details\n", ds->variant);
		}
		inla_read_prior(mb, ini, sec, &(ds->data_prior), "LOGGAMMA", NULL);

		/*
		 * add theta 
		 */
		if (!ds->data_fixed) {
			mb->theta = Realloc(mb->theta, mb->ntheta + 1, double **);
			mb->theta_hyperid = Realloc(mb->theta_hyperid, mb->ntheta + 1, char *);
			mb->theta_hyperid[mb->ntheta] = ds->data_prior.hyperid;
			mb->theta_tag = Realloc(mb->theta_tag, mb->ntheta + 1, char *);
			mb->theta_tag_userscale = Realloc(mb->theta_tag_userscale, mb->ntheta + 1, char *);
			mb->theta_dir = Realloc(mb->theta_dir, mb->ntheta + 1, char *);
			if (ds->variant == 0 || ds->variant == 1 || ds->variant == 2) {
				mb->theta_tag[mb->ntheta] = inla_make_tag("log size for the cennbinomial2 observations (1/overdispersion)", mb->ds);
				mb->theta_tag_userscale[mb->ntheta] =
				    inla_make_tag("size for the cennbinomial2 observations (1/overdispersion)", mb->ds);
			} else {
				assert(0 == 1);
			}
			GMRFLib_sprintf(&msg, "%s-parameter", secname);
			mb->theta_dir[mb->ntheta] = msg;

			mb->theta_from = Realloc(mb->theta_from, mb->ntheta + 1, char *);
			mb->theta_to = Realloc(mb->theta_to, mb->ntheta + 1, char *);
			mb->theta_from[mb->ntheta] = Strdup(ds->data_prior.from_theta);
			mb->theta_to[mb->ntheta] = Strdup(ds->data_prior.to_theta);

			mb->theta[mb->ntheta] = ds->data_observations.log_size;
			mb->theta_map = Realloc(mb->theta_map, mb->ntheta + 1, map_func_tp *);
			mb->theta_map[mb->ntheta] = map_exp;
			mb->theta_map_arg = Realloc(mb->theta_map_arg, mb->ntheta + 1, void *);
			mb->theta_map_arg[mb->ntheta] = NULL;
			mb->ntheta++;
			ds->data_ntheta++;
		}
	}
		break;

	case L_ZEROINFLATEDNBINOMIAL0:
	case L_ZEROINFLATEDNBINOMIAL1:
	{
		/*
		 * get options related to the zeroinflated negative binomial_0/1
		 */
		tmp = iniparser_getdouble(ini, inla_string_join(secname, "INITIAL0"), log(10.0));
		ds->data_fixed0 = iniparser_getboolean(ini, inla_string_join(secname, "FIXED0"), 0);
		if (!ds->data_fixed0 && mb->mode_use_mode) {
			tmp = mb->theta_file[mb->theta_counter_file++];
			if (mb->mode_fixed)
				ds->data_fixed0 = 1;
		}
		HYPER_NEW(ds->data_observations.log_size, tmp);
		if (mb->verbose) {
			printf("\t\tinitialise log_size[%g]\n", ds->data_observations.log_size[0][0]);
			printf("\t\tfixed=[%1d]\n", ds->data_fixed0);
		}
		inla_read_prior0(mb, ini, sec, &(ds->data_prior0), "LOGGAMMA", NULL);

		/*
		 * add theta 
		 */
		if (!ds->data_fixed0) {
			mb->theta = Realloc(mb->theta, mb->ntheta + 1, double **);
			mb->theta_hyperid = Realloc(mb->theta_hyperid, mb->ntheta + 1, char *);
			mb->theta_hyperid[mb->ntheta] = ds->data_prior0.hyperid;
			mb->theta_tag = Realloc(mb->theta_tag, mb->ntheta + 1, char *);
			mb->theta_tag_userscale = Realloc(mb->theta_tag_userscale, mb->ntheta + 1, char *);
			mb->theta_dir = Realloc(mb->theta_dir, mb->ntheta + 1, char *);
			if (ds->data_id == L_ZEROINFLATEDNBINOMIAL0) {
				mb->theta_tag[mb->ntheta] = inla_make_tag("log size for nbinomial_0 zero-inflated observations", mb->ds);
				mb->theta_tag_userscale[mb->ntheta] = inla_make_tag("size for nbinomial_0 zero-inflated observations", mb->ds);
			} else {
				mb->theta_tag[mb->ntheta] = inla_make_tag("log size for nbinomial_1 zero-inflated observations", mb->ds);
				mb->theta_tag_userscale[mb->ntheta] = inla_make_tag("size for nbinomial_1 zero-inflated observations", mb->ds);
			}
			GMRFLib_sprintf(&msg, "%s-parameter0", secname);
			mb->theta_dir[mb->ntheta] = msg;

			mb->theta_from = Realloc(mb->theta_from, mb->ntheta + 1, char *);
			mb->theta_to = Realloc(mb->theta_to, mb->ntheta + 1, char *);
			mb->theta_from[mb->ntheta] = Strdup(ds->data_prior0.from_theta);
			mb->theta_to[mb->ntheta] = Strdup(ds->data_prior0.to_theta);

			mb->theta[mb->ntheta] = ds->data_observations.log_size;
			mb->theta_map = Realloc(mb->theta_map, mb->ntheta + 1, map_func_tp *);
			mb->theta_map[mb->ntheta] = map_exp;
			mb->theta_map_arg = Realloc(mb->theta_map_arg, mb->ntheta + 1, void *);
			mb->theta_map_arg[mb->ntheta] = NULL;
			mb->ntheta++;
			ds->data_ntheta++;
		}

		/*
		 * the zeroinflation parameter 
		 */
		tmp = iniparser_getdouble(ini, inla_string_join(secname, "INITIAL1"), -1.0);
		ds->data_fixed1 = iniparser_getboolean(ini, inla_string_join(secname, "FIXED1"), 0);
		if (!ds->data_fixed1 && mb->mode_use_mode) {
			tmp = mb->theta_file[mb->theta_counter_file++];
			if (mb->mode_fixed)
				ds->data_fixed1 = 1;
		}
		HYPER_NEW(ds->data_observations.prob_intern, tmp);
		if (mb->verbose) {
			printf("\t\tinitialise prob_intern[%g]\n", ds->data_observations.prob_intern[0][0]);
			printf("\t\tfixed=[%1d]\n", ds->data_fixed1);
		}
		inla_read_prior1(mb, ini, sec, &(ds->data_prior1), "GAUSSIAN-std", NULL);

		/*
		 * add theta 
		 */
		if (!ds->data_fixed1) {
			mb->theta = Realloc(mb->theta, mb->ntheta + 1, double **);
			mb->theta_hyperid = Realloc(mb->theta_hyperid, mb->ntheta + 1, char *);
			mb->theta_hyperid[mb->ntheta] = ds->data_prior1.hyperid;
			mb->theta_tag = Realloc(mb->theta_tag, mb->ntheta + 1, char *);
			mb->theta_tag_userscale = Realloc(mb->theta_tag_userscale, mb->ntheta + 1, char *);
			mb->theta_dir = Realloc(mb->theta_dir, mb->ntheta + 1, char *);
			if (ds->data_id == L_ZEROINFLATEDNBINOMIAL0) {
				mb->theta_tag[mb->ntheta] =
				    inla_make_tag("intern zero-probability parameter for zero-inflated nbinomial_0", mb->ds);
				mb->theta_tag_userscale[mb->ntheta] =
				    inla_make_tag("zero-probability parameter for zero-inflated nbinomial_0", mb->ds);
			} else {
				mb->theta_tag[mb->ntheta] =
				    inla_make_tag("intern zero-probability parameter for zero-inflated nbinomial_1", mb->ds);
				mb->theta_tag_userscale[mb->ntheta] =
				    inla_make_tag("zero-probability parameter for zero-inflated nbinomial_1", mb->ds);
			}
			GMRFLib_sprintf(&msg, "%s-parameter1", secname);
			mb->theta_dir[mb->ntheta] = msg;

			mb->theta_from = Realloc(mb->theta_from, mb->ntheta + 1, char *);
			mb->theta_to = Realloc(mb->theta_to, mb->ntheta + 1, char *);
			mb->theta_from[mb->ntheta] = Strdup(ds->data_prior1.from_theta);
			mb->theta_to[mb->ntheta] = Strdup(ds->data_prior1.to_theta);

			mb->theta[mb->ntheta] = ds->data_observations.prob_intern;
			mb->theta_map = Realloc(mb->theta_map, mb->ntheta + 1, map_func_tp *);
			mb->theta_map[mb->ntheta] = map_probability;
			mb->theta_map_arg = Realloc(mb->theta_map_arg, mb->ntheta + 1, void *);
			mb->theta_map_arg[mb->ntheta] = NULL;
			mb->ntheta++;
			ds->data_ntheta++;
		}
	}
		break;

	case L_ZEROINFLATEDBETABINOMIAL0:
	case L_ZEROINFLATEDBETABINOMIAL1:
	{
		/*
		 * get options related to the zeroinflated beta binomial_0/1
		 */
		tmp = iniparser_getdouble(ini, inla_string_join(secname, "INITIAL0"), log(10.0));
		ds->data_fixed0 = iniparser_getboolean(ini, inla_string_join(secname, "FIXED0"), 0);
		if (!ds->data_fixed0 && mb->mode_use_mode) {
			tmp = mb->theta_file[mb->theta_counter_file++];
			if (mb->mode_fixed)
				ds->data_fixed0 = 1;
		}
		HYPER_NEW(ds->data_observations.zeroinflated_rho_intern, tmp);
		if (mb->verbose) {
			printf("\t\tinitialise rho_intern[%g]\n", ds->data_observations.zeroinflated_rho_intern[0][0]);
			printf("\t\tfixed=[%1d]\n", ds->data_fixed0);
		}
		inla_read_prior0(mb, ini, sec, &(ds->data_prior0), "GAUSSIAN-std", NULL);

		/*
		 * add theta 
		 */
		if (!ds->data_fixed0) {
			mb->theta = Realloc(mb->theta, mb->ntheta + 1, double **);
			mb->theta_hyperid = Realloc(mb->theta_hyperid, mb->ntheta + 1, char *);
			mb->theta_hyperid[mb->ntheta] = ds->data_prior0.hyperid;
			mb->theta_tag = Realloc(mb->theta_tag, mb->ntheta + 1, char *);
			mb->theta_tag_userscale = Realloc(mb->theta_tag_userscale, mb->ntheta + 1, char *);
			mb->theta_dir = Realloc(mb->theta_dir, mb->ntheta + 1, char *);

			if (ds->data_id == L_ZEROINFLATEDBETABINOMIAL0) {
				mb->theta_tag[mb->ntheta] = inla_make_tag("rho_intern for zero-inflated betabinomial_0", mb->ds);
				mb->theta_tag_userscale[mb->ntheta] = inla_make_tag("rho for zero-inflated betabinomial_0", mb->ds);
			} else {
				mb->theta_tag[mb->ntheta] = inla_make_tag("rho_intern for zero-inflated betabinomial_1", mb->ds);
				mb->theta_tag_userscale[mb->ntheta] = inla_make_tag("rho for zero-inflated betabinomial_1", mb->ds);
			}
			GMRFLib_sprintf(&msg, "%s-parameter0", secname);
			mb->theta_dir[mb->ntheta] = msg;

			mb->theta_from = Realloc(mb->theta_from, mb->ntheta + 1, char *);
			mb->theta_to = Realloc(mb->theta_to, mb->ntheta + 1, char *);
			mb->theta_from[mb->ntheta] = Strdup(ds->data_prior0.from_theta);
			mb->theta_to[mb->ntheta] = Strdup(ds->data_prior0.to_theta);

			mb->theta[mb->ntheta] = ds->data_observations.zeroinflated_rho_intern;
			mb->theta_map = Realloc(mb->theta_map, mb->ntheta + 1, map_func_tp *);
			mb->theta_map[mb->ntheta] = map_probability;
			mb->theta_map_arg = Realloc(mb->theta_map_arg, mb->ntheta + 1, void *);
			mb->theta_map_arg[mb->ntheta] = NULL;
			mb->ntheta++;
			ds->data_ntheta++;
		}

		/*
		 * the zeroinflation parameter 
		 */
		tmp = iniparser_getdouble(ini, inla_string_join(secname, "INITIAL1"), -1.0);
		ds->data_fixed1 = iniparser_getboolean(ini, inla_string_join(secname, "FIXED1"), 0);
		if (!ds->data_fixed1 && mb->mode_use_mode) {
			tmp = mb->theta_file[mb->theta_counter_file++];
			if (mb->mode_fixed)
				ds->data_fixed1 = 1;
		}
		HYPER_NEW(ds->data_observations.prob_intern, tmp);
		if (mb->verbose) {
			printf("\t\tinitialise prob_intern[%g]\n", ds->data_observations.prob_intern[0][0]);
			printf("\t\tfixed=[%1d]\n", ds->data_fixed1);
		}
		inla_read_prior1(mb, ini, sec, &(ds->data_prior1), "GAUSSIAN-std", NULL);

		/*
		 * add theta 
		 */
		if (!ds->data_fixed1) {
			mb->theta = Realloc(mb->theta, mb->ntheta + 1, double **);
			mb->theta_hyperid = Realloc(mb->theta_hyperid, mb->ntheta + 1, char *);
			mb->theta_hyperid[mb->ntheta] = ds->data_prior1.hyperid;
			mb->theta_tag = Realloc(mb->theta_tag, mb->ntheta + 1, char *);
			mb->theta_tag_userscale = Realloc(mb->theta_tag_userscale, mb->ntheta + 1, char *);
			mb->theta_dir = Realloc(mb->theta_dir, mb->ntheta + 1, char *);
			if (ds->data_id == L_ZEROINFLATEDBETABINOMIAL0) {
				mb->theta_tag[mb->ntheta] =
				    inla_make_tag("intern zero-probability parameter for zero-inflated betabinomial_0", mb->ds);
				mb->theta_tag_userscale[mb->ntheta] =
				    inla_make_tag("zero-probability parameter for zero-inflated betabinomial_0", mb->ds);
			} else {
				mb->theta_tag[mb->ntheta] =
				    inla_make_tag("intern zero-probability parameter for zero-inflated betabinomial_1", mb->ds);
				mb->theta_tag_userscale[mb->ntheta] =
				    inla_make_tag("zero-probability parameter for zero-inflated betabinomial_1", mb->ds);
			}
			GMRFLib_sprintf(&msg, "%s-parameter1", secname);
			mb->theta_dir[mb->ntheta] = msg;

			mb->theta_from = Realloc(mb->theta_from, mb->ntheta + 1, char *);
			mb->theta_to = Realloc(mb->theta_to, mb->ntheta + 1, char *);
			mb->theta_from[mb->ntheta] = Strdup(ds->data_prior1.from_theta);
			mb->theta_to[mb->ntheta] = Strdup(ds->data_prior1.to_theta);

			mb->theta[mb->ntheta] = ds->data_observations.prob_intern;
			mb->theta_map = Realloc(mb->theta_map, mb->ntheta + 1, map_func_tp *);
			mb->theta_map[mb->ntheta] = map_probability;
			mb->theta_map_arg = Realloc(mb->theta_map_arg, mb->ntheta + 1, void *);
			mb->theta_map_arg[mb->ntheta] = NULL;
			mb->ntheta++;
			ds->data_ntheta++;
		}
	}
		break;

	case L_ZEROINFLATEDNBINOMIAL1STRATA2:
	{
		/*
		 * get options related to the zeroinflated negative binomial_0/1, strata2
		 */
		tmp = iniparser_getdouble(ini, inla_string_join(secname, "INITIAL0"), log(10.0));
		ds->data_fixed0 = iniparser_getboolean(ini, inla_string_join(secname, "FIXED0"), 0);
		if (!ds->data_fixed0 && mb->mode_use_mode) {
			tmp = mb->theta_file[mb->theta_counter_file++];
			if (mb->mode_fixed)
				ds->data_fixed0 = 1;
		}
		HYPER_NEW(ds->data_observations.log_size, tmp);
		if (mb->verbose) {
			printf("\t\tinitialise log_size[%g]\n", ds->data_observations.log_size[0][0]);
			printf("\t\tfixed=[%1d]\n", ds->data_fixed0);
		}
		inla_read_prior0(mb, ini, sec, &(ds->data_prior0), "LOGGAMMA", NULL);

		/*
		 * add theta 
		 */
		if (!ds->data_fixed0) {
			mb->theta = Realloc(mb->theta, mb->ntheta + 1, double **);
			mb->theta_hyperid = Realloc(mb->theta_hyperid, mb->ntheta + 1, char *);
			mb->theta_hyperid[mb->ntheta] = ds->data_prior0.hyperid;
			mb->theta_tag = Realloc(mb->theta_tag, mb->ntheta + 1, char *);
			mb->theta_tag_userscale = Realloc(mb->theta_tag_userscale, mb->ntheta + 1, char *);
			mb->theta_dir = Realloc(mb->theta_dir, mb->ntheta + 1, char *);
			mb->theta_tag[mb->ntheta] = inla_make_tag("log size for zero-inflated nbinomial_1_strata2", mb->ds);
			mb->theta_tag_userscale[mb->ntheta] = inla_make_tag("size for zero-inflated nbinomial_1_strata2", mb->ds);
			GMRFLib_sprintf(&msg, "%s-parameter0", secname);
			mb->theta_dir[mb->ntheta] = msg;

			mb->theta_from = Realloc(mb->theta_from, mb->ntheta + 1, char *);
			mb->theta_to = Realloc(mb->theta_to, mb->ntheta + 1, char *);
			mb->theta_from[mb->ntheta] = Strdup(ds->data_prior0.from_theta);
			mb->theta_to[mb->ntheta] = Strdup(ds->data_prior0.to_theta);

			mb->theta[mb->ntheta] = ds->data_observations.log_size;
			mb->theta_map = Realloc(mb->theta_map, mb->ntheta + 1, map_func_tp *);
			mb->theta_map[mb->ntheta] = map_exp;
			mb->theta_map_arg = Realloc(mb->theta_map_arg, mb->ntheta + 1, void *);
			mb->theta_map_arg[mb->ntheta] = NULL;
			mb->ntheta++;
			ds->data_ntheta++;
		}

		/*
		 * THERE are up to STRATA_MAXTHETA of the probs, called prob 1... 10 
		 */
		ds->data_observations.probN_intern = Calloc(STRATA_MAXTHETA, double **);
		ds->data_nfixed = Calloc(STRATA_MAXTHETA, int);
		ds->data_nprior = Calloc(STRATA_MAXTHETA, Prior_tp);

		for (int count = 0; count < STRATA_MAXTHETA; count++) {
			/*
			 * the zeroinflation prob-parameter
			 */
			GMRFLib_sprintf(&ctmp, "INITIAL%1d", count + 1);
			tmp = iniparser_getdouble(ini, inla_string_join(secname, ctmp), -1.0);

			GMRFLib_sprintf(&ctmp, "FIXED%1d", count + 1);
			ds->data_nfixed[count] = iniparser_getboolean(ini, inla_string_join(secname, ctmp), 0);
			if (!ds->data_nfixed[count] && mb->mode_use_mode) {
				tmp = mb->theta_file[mb->theta_counter_file++];
				if (mb->mode_fixed)
					ds->data_nfixed[count] = 1;
			}
			HYPER_NEW(ds->data_observations.probN_intern[count], tmp);

			if (mb->verbose) {
				printf("\t\tinitialise prob%1d_intern[%g]\n", count + 1, ds->data_observations.probN_intern[count][0][0]);
				printf("\t\tfixed%1d=[%1d]\n", count + 1, ds->data_nfixed[count]);
			}
			inla_read_priorN(mb, ini, sec, &(ds->data_nprior[count]), "GAUSSIAN-std", count + 1, NULL);

			/*
			 * add theta 
			 */
			if (!ds->data_nfixed[count]) {
				mb->theta = Realloc(mb->theta, mb->ntheta + 1, double **);
				mb->theta_hyperid = Realloc(mb->theta_hyperid, mb->ntheta + 1, char *);
				mb->theta_hyperid[mb->ntheta] = ds->data_nprior[count].hyperid;
				mb->theta_tag = Realloc(mb->theta_tag, mb->ntheta + 1, char *);
				mb->theta_tag_userscale = Realloc(mb->theta_tag_userscale, mb->ntheta + 1, char *);
				mb->theta_dir = Realloc(mb->theta_dir, mb->ntheta + 1, char *);

				GMRFLib_sprintf(&ctmp, "intern zero-probability%1d for zero-inflated nbinomial_1_strata2", count + 1);
				mb->theta_tag[mb->ntheta] = inla_make_tag(ctmp, mb->ds);
				GMRFLib_sprintf(&ctmp, "zero-probability%1d for zero-inflated nbinomial_1_strata2", count + 1);
				mb->theta_tag_userscale[mb->ntheta] = inla_make_tag(ctmp, mb->ds);
				GMRFLib_sprintf(&msg, "%s-parameter%1d", secname, count + 1);
				mb->theta_dir[mb->ntheta] = msg;

				mb->theta_from = Realloc(mb->theta_from, mb->ntheta + 1, char *);
				mb->theta_to = Realloc(mb->theta_to, mb->ntheta + 1, char *);
				mb->theta_from[mb->ntheta] = Strdup(ds->data_nprior[count].from_theta);
				mb->theta_to[mb->ntheta] = Strdup(ds->data_nprior[count].to_theta);

				mb->theta[mb->ntheta] = ds->data_observations.probN_intern[count];
				mb->theta_map = Realloc(mb->theta_map, mb->ntheta + 1, map_func_tp *);
				mb->theta_map[mb->ntheta] = map_probability;
				mb->theta_map_arg = Realloc(mb->theta_map_arg, mb->ntheta + 1, void *);
				mb->theta_map_arg[mb->ntheta] = NULL;
				mb->ntheta++;
				ds->data_ntheta++;
			}
		}
	}
		break;

	case L_ZEROINFLATEDNBINOMIAL1STRATA3:
	{
		/*
		 * get options related to the zeroinflated negative binomial_0/1, strata3
		 */
		tmp = iniparser_getdouble(ini, inla_string_join(secname, "INITIAL0"), log(10.0));
		ds->data_fixed0 = iniparser_getboolean(ini, inla_string_join(secname, "FIXED0"), 0);
		if (!ds->data_fixed0 && mb->mode_use_mode) {
			tmp = mb->theta_file[mb->theta_counter_file++];
			if (mb->mode_fixed)
				ds->data_fixed0 = 1;
		}
		HYPER_NEW(ds->data_observations.prob_intern, tmp);
		if (mb->verbose) {
			printf("\t\tinitialise prob_intern[%g]\n", ds->data_observations.prob_intern[0][0]);
			printf("\t\tfixed0=[%1d]\n", ds->data_fixed0);
		}
		inla_read_prior0(mb, ini, sec, &(ds->data_prior0), "GAUSSIAN-std", NULL);

		/*
		 * add theta 
		 */
		if (!ds->data_fixed0) {
			mb->theta = Realloc(mb->theta, mb->ntheta + 1, double **);
			mb->theta_hyperid = Realloc(mb->theta_hyperid, mb->ntheta + 1, char *);
			mb->theta_hyperid[mb->ntheta] = ds->data_prior0.hyperid;
			mb->theta_tag = Realloc(mb->theta_tag, mb->ntheta + 1, char *);
			mb->theta_tag_userscale = Realloc(mb->theta_tag_userscale, mb->ntheta + 1, char *);
			mb->theta_dir = Realloc(mb->theta_dir, mb->ntheta + 1, char *);
			mb->theta_tag[mb->ntheta] = inla_make_tag("intern zero-probability for zero-inflated nbinomial_1_strata3", mb->ds);
			mb->theta_tag_userscale[mb->ntheta] = inla_make_tag("zero-probability for zero-inflated nbinomial_1_strata3", mb->ds);
			GMRFLib_sprintf(&msg, "%s-parameter0", secname);
			mb->theta_dir[mb->ntheta] = msg;

			mb->theta_from = Realloc(mb->theta_from, mb->ntheta + 1, char *);
			mb->theta_to = Realloc(mb->theta_to, mb->ntheta + 1, char *);
			mb->theta_from[mb->ntheta] = Strdup(ds->data_prior0.from_theta);
			mb->theta_to[mb->ntheta] = Strdup(ds->data_prior0.to_theta);

			mb->theta[mb->ntheta] = ds->data_observations.prob_intern;
			mb->theta_map = Realloc(mb->theta_map, mb->ntheta + 1, map_func_tp *);
			mb->theta_map[mb->ntheta] = map_probability;
			mb->theta_map_arg = Realloc(mb->theta_map_arg, mb->ntheta + 1, void *);
			mb->theta_map_arg[mb->ntheta] = NULL;
			mb->ntheta++;
			ds->data_ntheta++;
		}

		/*
		 * THERE are up to STRATA_MAXTHETA of the size, called size 1... 10 
		 */
		ds->data_observations.log_sizes = Calloc(STRATA_MAXTHETA, double **);
		ds->data_nfixed = Calloc(STRATA_MAXTHETA, int);
		ds->data_nprior = Calloc(STRATA_MAXTHETA, Prior_tp);

		for (int count = 0; count < STRATA_MAXTHETA; count++) {
			/*
			 * the zeroinflation prob-parameter
			 */
			GMRFLib_sprintf(&ctmp, "INITIAL%1d", count + 1);
			tmp = iniparser_getdouble(ini, inla_string_join(secname, ctmp), -1.0);

			GMRFLib_sprintf(&ctmp, "FIXED%1d", count + 1);
			ds->data_nfixed[count] = iniparser_getboolean(ini, inla_string_join(secname, ctmp), 0);
			if (!ds->data_nfixed[count] && mb->mode_use_mode) {
				tmp = mb->theta_file[mb->theta_counter_file++];
				if (mb->mode_fixed)
					ds->data_nfixed[count] = 1;
			}
			HYPER_NEW(ds->data_observations.log_sizes[count], tmp);

			if (mb->verbose) {
				printf("\t\tinitialise log_size%1d[%g]\n", count + 1, ds->data_observations.log_sizes[count][0][0]);
				printf("\t\tfixed%1d=[%1d]\n", count + 1, ds->data_nfixed[count]);
			}
			inla_read_priorN(mb, ini, sec, &(ds->data_nprior[count]), "LOGGAMMA", count + 1, NULL);

			/*
			 * add theta 
			 */
			if (!ds->data_nfixed[count]) {
				mb->theta = Realloc(mb->theta, mb->ntheta + 1, double **);
				mb->theta_hyperid = Realloc(mb->theta_hyperid, mb->ntheta + 1, char *);
				mb->theta_hyperid[mb->ntheta] = ds->data_nprior[count].hyperid;
				mb->theta_tag = Realloc(mb->theta_tag, mb->ntheta + 1, char *);
				mb->theta_tag_userscale = Realloc(mb->theta_tag_userscale, mb->ntheta + 1, char *);
				mb->theta_dir = Realloc(mb->theta_dir, mb->ntheta + 1, char *);

				GMRFLib_sprintf(&ctmp, "log_size%1d for zero-inflated nbinomial_strata3", count + 1);
				mb->theta_tag[mb->ntheta] = inla_make_tag(ctmp, mb->ds);
				GMRFLib_sprintf(&ctmp, "size%1d for zero-inflated nbinomial_strata3", count + 1);
				mb->theta_tag_userscale[mb->ntheta] = inla_make_tag(ctmp, mb->ds);
				GMRFLib_sprintf(&msg, "%s-parameter%1d", secname, count + 1);
				mb->theta_dir[mb->ntheta] = msg;

				mb->theta_from = Realloc(mb->theta_from, mb->ntheta + 1, char *);
				mb->theta_to = Realloc(mb->theta_to, mb->ntheta + 1, char *);
				mb->theta_from[mb->ntheta] = Strdup(ds->data_nprior[count].from_theta);
				mb->theta_to[mb->ntheta] = Strdup(ds->data_nprior[count].to_theta);

				mb->theta[mb->ntheta] = ds->data_observations.log_sizes[count];
				mb->theta_map = Realloc(mb->theta_map, mb->ntheta + 1, map_func_tp *);
				mb->theta_map[mb->ntheta] = map_exp;
				mb->theta_map_arg = Realloc(mb->theta_map_arg, mb->ntheta + 1, void *);
				mb->theta_map_arg[mb->ntheta] = NULL;
				mb->ntheta++;
				ds->data_ntheta++;
			}
		}
	}
		break;

	case L_ZEROINFLATEDNBINOMIAL2:
	{
		/*
		 * get options related to the zeroinflated negative binomial_2
		 */
		tmp = iniparser_getdouble(ini, inla_string_join(secname, "INITIAL0"), log(10.0));
		ds->data_fixed0 = iniparser_getboolean(ini, inla_string_join(secname, "FIXED0"), 0);
		if (!ds->data_fixed0 && mb->mode_use_mode) {
			tmp = mb->theta_file[mb->theta_counter_file++];
			if (mb->mode_fixed)
				ds->data_fixed0 = 1;
		}
		HYPER_NEW(ds->data_observations.log_size, tmp);
		if (mb->verbose) {
			printf("\t\tinitialise log_size[%g]\n", ds->data_observations.log_size[0][0]);
			printf("\t\tfixed=[%1d]\n", ds->data_fixed0);
		}
		inla_read_prior0(mb, ini, sec, &(ds->data_prior0), "LOGGAMMA", NULL);

		/*
		 * add theta 
		 */
		if (!ds->data_fixed0) {
			mb->theta = Realloc(mb->theta, mb->ntheta + 1, double **);
			mb->theta_hyperid = Realloc(mb->theta_hyperid, mb->ntheta + 1, char *);
			mb->theta_hyperid[mb->ntheta] = ds->data_prior0.hyperid;
			mb->theta_tag = Realloc(mb->theta_tag, mb->ntheta + 1, char *);
			mb->theta_tag_userscale = Realloc(mb->theta_tag_userscale, mb->ntheta + 1, char *);
			mb->theta_dir = Realloc(mb->theta_dir, mb->ntheta + 1, char *);
			mb->theta_tag[mb->ntheta] = inla_make_tag("log size for nbinomial zero-inflated observations", mb->ds);
			mb->theta_tag_userscale[mb->ntheta] = inla_make_tag("size for nbinomial zero-inflated observations", mb->ds);
			GMRFLib_sprintf(&msg, "%s-parameter0", secname);
			mb->theta_dir[mb->ntheta] = msg;

			mb->theta_from = Realloc(mb->theta_from, mb->ntheta + 1, char *);
			mb->theta_to = Realloc(mb->theta_to, mb->ntheta + 1, char *);
			mb->theta_from[mb->ntheta] = Strdup(ds->data_prior0.from_theta);
			mb->theta_to[mb->ntheta] = Strdup(ds->data_prior0.to_theta);

			mb->theta[mb->ntheta] = ds->data_observations.log_size;
			mb->theta_map = Realloc(mb->theta_map, mb->ntheta + 1, map_func_tp *);
			mb->theta_map[mb->ntheta] = map_exp;
			mb->theta_map_arg = Realloc(mb->theta_map_arg, mb->ntheta + 1, void *);
			mb->theta_map_arg[mb->ntheta] = NULL;
			mb->ntheta++;
			ds->data_ntheta++;
		}

		/*
		 * the zeroinflation parameter; the parameter alpha (see the documentation) 
		 */
		tmp = iniparser_getdouble(ini, inla_string_join(secname, "INITIAL1"), log(2.0));
		ds->data_fixed1 = iniparser_getboolean(ini, inla_string_join(secname, "FIXED1"), 0);
		if (!ds->data_fixed1 && mb->mode_use_mode) {
			tmp = mb->theta_file[mb->theta_counter_file++];
			if (mb->mode_fixed)
				ds->data_fixed1 = 1;
		}
		HYPER_NEW(ds->data_observations.zeroinflated_alpha_intern, tmp);
		if (mb->verbose) {
			printf("\t\tinitialise alpha_intern[%g]\n", ds->data_observations.zeroinflated_alpha_intern[0][0]);
			printf("\t\tfixed=[%1d]\n", ds->data_fixed1);
		}
		inla_read_prior1(mb, ini, sec, &(ds->data_prior1), "gaussian-std", NULL);

		/*
		 * add theta 
		 */
		if (!ds->data_fixed1) {
			mb->theta = Realloc(mb->theta, mb->ntheta + 1, double **);
			mb->theta_hyperid = Realloc(mb->theta_hyperid, mb->ntheta + 1, char *);
			mb->theta_hyperid[mb->ntheta] = ds->data_prior1.hyperid;
			mb->theta_tag = Realloc(mb->theta_tag, mb->ntheta + 1, char *);
			mb->theta_tag_userscale = Realloc(mb->theta_tag_userscale, mb->ntheta + 1, char *);
			mb->theta_dir = Realloc(mb->theta_dir, mb->ntheta + 1, char *);
			mb->theta_tag[mb->ntheta] = inla_make_tag("parameter alpha.intern for zero-inflated nbinomial2", mb->ds);
			mb->theta_tag_userscale[mb->ntheta] = inla_make_tag("parameter alpha for zero-inflated nbinomial2", mb->ds);
			GMRFLib_sprintf(&msg, "%s-parameter1", secname);
			mb->theta_dir[mb->ntheta] = msg;

			mb->theta_from = Realloc(mb->theta_from, mb->ntheta + 1, char *);
			mb->theta_to = Realloc(mb->theta_to, mb->ntheta + 1, char *);
			mb->theta_from[mb->ntheta] = Strdup(ds->data_prior1.from_theta);
			mb->theta_to[mb->ntheta] = Strdup(ds->data_prior1.to_theta);

			mb->theta[mb->ntheta] = ds->data_observations.zeroinflated_alpha_intern;
			mb->theta_map = Realloc(mb->theta_map, mb->ntheta + 1, map_func_tp *);
			mb->theta_map[mb->ntheta] = map_exp;
			mb->theta_map_arg = Realloc(mb->theta_map_arg, mb->ntheta + 1, void *);
			mb->theta_map_arg[mb->ntheta] = NULL;
			mb->ntheta++;
			ds->data_ntheta++;
		}
	}
		break;

	case L_T:
	{
		/*
		 * get options related to the t
		 */

		tmp = iniparser_getdouble(ini, inla_string_join(secname, "INITIAL0"), G.log_prec_initial);
		ds->data_fixed0 = iniparser_getboolean(ini, inla_string_join(secname, "FIXED0"), 0);
		if (!ds->data_fixed0 && mb->mode_use_mode) {
			tmp = mb->theta_file[mb->theta_counter_file++];
			if (mb->mode_fixed)
				ds->data_fixed0 = 1;
		}
		HYPER_NEW(ds->data_observations.log_prec_t, tmp);
		if (mb->verbose) {
			printf("\t\tinitialise log_precision[%g]\n", ds->data_observations.log_prec_t[0][0]);
			printf("\t\tfixed0=[%1d]\n", ds->data_fixed0);
		}
		inla_read_prior0(mb, ini, sec, &(ds->data_prior0), "LOGGAMMA", NULL);

		if (!ds->data_fixed0) {
			mb->theta = Realloc(mb->theta, mb->ntheta + 1, double **);
			mb->theta_hyperid = Realloc(mb->theta_hyperid, mb->ntheta + 1, char *);
			mb->theta_hyperid[mb->ntheta] = ds->data_prior0.hyperid;
			mb->theta_tag = Realloc(mb->theta_tag, mb->ntheta + 1, char *);
			mb->theta_tag_userscale = Realloc(mb->theta_tag_userscale, mb->ntheta + 1, char *);
			mb->theta_dir = Realloc(mb->theta_dir, mb->ntheta + 1, char *);
			mb->theta_tag[mb->ntheta] = inla_make_tag("log precision for the student-t observations", mb->ds);
			mb->theta_tag_userscale[mb->ntheta] = inla_make_tag("precision for the student-t observations", mb->ds);
			GMRFLib_sprintf(&msg, "%s-parameter0", secname);
			mb->theta_dir[mb->ntheta] = msg;

			mb->theta_from = Realloc(mb->theta_from, mb->ntheta + 1, char *);
			mb->theta_to = Realloc(mb->theta_to, mb->ntheta + 1, char *);
			mb->theta_from[mb->ntheta] = Strdup(ds->data_prior0.from_theta);
			mb->theta_to[mb->ntheta] = Strdup(ds->data_prior0.to_theta);

			mb->theta[mb->ntheta] = ds->data_observations.log_prec_t;
			mb->theta_map = Realloc(mb->theta_map, mb->ntheta + 1, map_func_tp *);
			mb->theta_map[mb->ntheta] = map_precision;
			mb->theta_map_arg = Realloc(mb->theta_map_arg, mb->ntheta + 1, void *);
			mb->theta_map_arg[mb->ntheta] = NULL;
			mb->ntheta++;
			ds->data_ntheta++;
		}
		/*
		 * dof 
		 */
		tmp = iniparser_getdouble(ini, inla_string_join(secname, "INITIAL1"), 3.0);
		if (ISINF(tmp) || ISNAN(tmp)) {
			GMRFLib_sprintf(&msg, "Initial value for dof_internal is void [%g]\n", tmp);
			inla_error_general(msg);
			exit(1);
		}
		ds->data_fixed1 = iniparser_getboolean(ini, inla_string_join(secname, "FIXED1"), 0);
		if (!ds->data_fixed1 && mb->mode_use_mode) {
			tmp = mb->theta_file[mb->theta_counter_file++];
			if (mb->mode_fixed)
				ds->data_fixed1 = 1;
		}
		HYPER_NEW(ds->data_observations.dof_intern_t, tmp);
		if (mb->verbose) {
			printf("\t\tinitialise dof_intern_t[%g]\n", ds->data_observations.dof_intern_t[0][0]);
			printf("\t\tfixed1=[%1d]\n", ds->data_fixed1);
		}
		inla_read_prior1(mb, ini, sec, &(ds->data_prior1), "normal", NULL);

		if (!ds->data_fixed1) {
			mb->theta = Realloc(mb->theta, mb->ntheta + 1, double **);
			mb->theta_hyperid = Realloc(mb->theta_hyperid, mb->ntheta + 1, char *);
			mb->theta_hyperid[mb->ntheta] = ds->data_prior1.hyperid;
			mb->theta_tag = Realloc(mb->theta_tag, mb->ntheta + 1, char *);
			mb->theta_tag_userscale = Realloc(mb->theta_tag_userscale, mb->ntheta + 1, char *);
			mb->theta_dir = Realloc(mb->theta_dir, mb->ntheta + 1, char *);
			mb->theta_tag[mb->ntheta] = inla_make_tag("dof_intern for student-t", mb->ds);
			mb->theta_tag_userscale[mb->ntheta] = inla_make_tag("degrees of freedom for student-t", mb->ds);
			GMRFLib_sprintf(&msg, "%s-parameter1", secname);
			mb->theta_dir[mb->ntheta] = msg;

			mb->theta_from = Realloc(mb->theta_from, mb->ntheta + 1, char *);
			mb->theta_to = Realloc(mb->theta_to, mb->ntheta + 1, char *);
			mb->theta_from[mb->ntheta] = Strdup(ds->data_prior1.from_theta);
			mb->theta_to[mb->ntheta] = Strdup(ds->data_prior1.to_theta);

			mb->theta[mb->ntheta] = ds->data_observations.dof_intern_t;
			mb->theta_map = Realloc(mb->theta_map, mb->ntheta + 1, map_func_tp *);
			mb->theta_map[mb->ntheta] = map_dof;
			mb->theta_map_arg = Realloc(mb->theta_map_arg, mb->ntheta + 1, void *);
			mb->theta_map_arg[mb->ntheta] = NULL;
			mb->ntheta++;
			ds->data_ntheta++;
		}
	}
		break;

	case L_TSTRATA:
	{
		/*
		 * get options related to the tstrata
		 */
		int k;

		ds->data_nprior = Calloc(TSTRATA_MAXTHETA, Prior_tp);
		ds->data_nfixed = Calloc(TSTRATA_MAXTHETA, int);

		/*
		 * check how many strata we have 
		 */
		int nstrata = 0;
		for (k = 0; k < mb->predictor_ndata; k++) {
			if (ds->data_observations.d[k]) {
				nstrata = IMAX(nstrata, (int) ds->data_observations.strata_tstrata[k]);
			}
		}
		nstrata++;
		assert(nstrata <= TSTRATA_MAXTHETA - 1);

		/*
		 * dof = theta0
		 */
		tmp = iniparser_getdouble(ini, inla_string_join(secname, "INITIAL0"), 3.0);
		if (ISINF(tmp) || ISNAN(tmp)) {
			GMRFLib_sprintf(&msg, "Initial value for dof_internal is void [%g]\n", tmp);
			inla_error_general(msg);
			exit(1);
		}
		ds->data_nfixed[0] = iniparser_getboolean(ini, inla_string_join(secname, "FIXED0"), 0);
		if (!ds->data_nfixed[0] && mb->mode_use_mode) {
			tmp = mb->theta_file[mb->theta_counter_file++];
			if (mb->mode_fixed)
				ds->data_nfixed[0] = 1;
		}
		HYPER_NEW(ds->data_observations.dof_intern_tstrata, tmp);
		if (mb->verbose) {
			printf("\t\tinitialise dof_intern_tstrata[%g]\n", ds->data_observations.dof_intern_tstrata[0][0]);
			printf("\t\tfixed0=[%1d]\n", ds->data_nfixed[0]);
		}
		inla_read_prior0(mb, ini, sec, &(ds->data_nprior[0]), "normal", NULL);

		if (!ds->data_nfixed[0]) {
			mb->theta = Realloc(mb->theta, mb->ntheta + 1, double **);
			mb->theta_hyperid = Realloc(mb->theta_hyperid, mb->ntheta + 1, char *);
			mb->theta_hyperid[mb->ntheta] = ds->data_nprior[0].hyperid;
			mb->theta_tag = Realloc(mb->theta_tag, mb->ntheta + 1, char *);
			mb->theta_tag_userscale = Realloc(mb->theta_tag_userscale, mb->ntheta + 1, char *);
			mb->theta_dir = Realloc(mb->theta_dir, mb->ntheta + 1, char *);
			mb->theta_tag[mb->ntheta] = inla_make_tag("dof_intern for tstrata", mb->ds);
			mb->theta_tag_userscale[mb->ntheta] = inla_make_tag("degrees of freedom for tstrata", mb->ds);
			GMRFLib_sprintf(&msg, "%s-parameter0", secname);
			mb->theta_dir[mb->ntheta] = msg;

			mb->theta_from = Realloc(mb->theta_from, mb->ntheta + 1, char *);
			mb->theta_to = Realloc(mb->theta_to, mb->ntheta + 1, char *);
			mb->theta_from[mb->ntheta] = Strdup(ds->data_nprior[0].from_theta);
			mb->theta_to[mb->ntheta] = Strdup(ds->data_nprior[0].to_theta);

			mb->theta[mb->ntheta] = ds->data_observations.dof_intern_tstrata;
			mb->theta_map = Realloc(mb->theta_map, mb->ntheta + 1, map_func_tp *);
			mb->theta_map[mb->ntheta] = map_dof;
			mb->theta_map_arg = Realloc(mb->theta_map_arg, mb->ntheta + 1, void *);
			mb->theta_map_arg[mb->ntheta] = NULL;
			mb->ntheta++;
			ds->data_ntheta++;
		}

		ds->data_observations.log_prec_tstrata = Calloc(TSTRATA_MAXTHETA - 1, double **);

		for (k = 1; k < TSTRATA_MAXTHETA; k++) {

			char *pri = NULL, *par = NULL, *from_theta = NULL, *to_theta = NULL, *hyperid = NULL;

			GMRFLib_sprintf(&pri, "PRIOR%1d", k);
			GMRFLib_sprintf(&par, "PARAMETERS%1d", k);
			GMRFLib_sprintf(&from_theta, "FROM.THETA%1d", k);
			GMRFLib_sprintf(&to_theta, "TO.THETA%1d", k);
			GMRFLib_sprintf(&hyperid, "HYPERID%1d", k);
			inla_read_prior_generic(mb, ini, sec, &(ds->data_nprior[k]), pri, par, from_theta, to_theta, hyperid, "normal", NULL);

			GMRFLib_sprintf(&ctmp, "FIXED%1d", k);
			ds->data_nfixed[k] = iniparser_getboolean(ini, inla_string_join(secname, ctmp), 0);
			/*
			 * if above number of stata, then its fixed for sure! 
			 */
			if (k > nstrata) {
				ds->data_nfixed[k] = 1;
			}
			GMRFLib_sprintf(&ctmp, "INITIAL%1d", k);
			double initial = iniparser_getdouble(ini, inla_string_join(secname, ctmp), G.log_prec_initial);

			if (!ds->data_nfixed[k] && mb->mode_use_mode) {
				initial = mb->theta_file[mb->theta_counter_file++];
				if (mb->mode_fixed)
					ds->data_nfixed[k] = 1;
			}
			HYPER_NEW(ds->data_observations.log_prec_tstrata[k - 1], initial);	/* yes, its a -1, prec0, prec1, etc... */
			if (mb->verbose) {
				printf("\t\tinitialise log_prec_tstrata[%1d][%g]\n", ds->data_nfixed[k],
				       ds->data_observations.log_prec_tstrata[k - 1][0][0]);
				printf("\t\tfixed%1d=[%1d]\n", k, ds->data_nfixed[k]);
			}

			if (!ds->data_nfixed[k]) {

				mb->theta = Realloc(mb->theta, mb->ntheta + 1, double **);
				mb->theta_hyperid = Realloc(mb->theta_hyperid, mb->ntheta + 1, char *);
				mb->theta_hyperid[mb->ntheta] = ds->data_nprior[k].hyperid;
				mb->theta_tag = Realloc(mb->theta_tag, mb->ntheta + 1, char *);
				mb->theta_tag_userscale = Realloc(mb->theta_tag_userscale, mb->ntheta + 1, char *);
				mb->theta_dir = Realloc(mb->theta_dir, mb->ntheta + 1, char *);

				mb->theta_tag[mb->ntheta] = inla_make_tag("Log prec for tstrata strata", k - 1);
				mb->theta_tag_userscale[mb->ntheta] = inla_make_tag("Prec tstrata strata", k - 1);
				GMRFLib_sprintf(&msg, "%s-parameter%1d", secname, k);
				mb->theta_dir[mb->ntheta] = msg;

				mb->theta_from = Realloc(mb->theta_from, mb->ntheta + 1, char *);
				mb->theta_to = Realloc(mb->theta_to, mb->ntheta + 1, char *);
				mb->theta_from[mb->ntheta] = Strdup(ds->data_nprior[k].from_theta);
				mb->theta_to[mb->ntheta] = Strdup(ds->data_nprior[k].to_theta);

				mb->theta[mb->ntheta] = ds->data_observations.log_prec_tstrata[k - 1];	/* yes its a -1 */
				mb->theta_map = Realloc(mb->theta_map, mb->ntheta + 1, map_func_tp *);
				mb->theta_map[mb->ntheta] = map_precision;
				mb->theta_map_arg = Realloc(mb->theta_map_arg, mb->ntheta + 1, void *);
				mb->theta_map_arg[mb->ntheta] = NULL;
				mb->ntheta++;
				ds->data_ntheta++;
			}

			Free(pri);
			Free(par);
			Free(from_theta);
			Free(to_theta);
			Free(hyperid);
			Free(ctmp);
		}
	}
		break;

	case L_STOCHVOL:
	{
		/*
		 * get options related to the stochvol; the log-offset in the variance
		 */
		double initial_value = 500.0;

		tmp = iniparser_getdouble(ini, inla_string_join(secname, "INITIAL"), initial_value);
		ds->data_fixed = iniparser_getboolean(ini, inla_string_join(secname, "FIXED"), 1);	/* yes, default fixed */
		if (!ds->data_fixed && mb->mode_use_mode) {
			tmp = mb->theta_file[mb->theta_counter_file++];
			if (mb->mode_fixed)
				ds->data_fixed = 1;
		}
		HYPER_NEW(ds->data_observations.log_offset_prec, tmp);
		if (mb->verbose) {
			printf("\t\tinitialise log.offset[%g]\n", ds->data_observations.log_offset_prec[0][0]);
			printf("\t\tfixed=[%1d]\n", ds->data_fixed);
		}
		inla_read_prior(mb, ini, sec, &(ds->data_prior), "LOGGAMMA", NULL);

		/*
		 * add theta 
		 */
		if (!ds->data_fixed) {
			mb->theta = Realloc(mb->theta, mb->ntheta + 1, double **);
			mb->theta_hyperid = Realloc(mb->theta_hyperid, mb->ntheta + 1, char *);
			mb->theta_hyperid[mb->ntheta] = ds->data_prior.hyperid;
			mb->theta_tag = Realloc(mb->theta_tag, mb->ntheta + 1, char *);
			mb->theta_tag_userscale = Realloc(mb->theta_tag_userscale, mb->ntheta + 1, char *);
			mb->theta_dir = Realloc(mb->theta_dir, mb->ntheta + 1, char *);
			mb->theta_tag[mb->ntheta] = inla_make_tag("Log offset precision for stochvol", mb->ds);
			mb->theta_tag_userscale[mb->ntheta] = inla_make_tag("Offset precision for stochvol", mb->ds);
			GMRFLib_sprintf(&msg, "%s-parameter", secname);
			mb->theta_dir[mb->ntheta] = msg;

			mb->theta_from = Realloc(mb->theta_from, mb->ntheta + 1, char *);
			mb->theta_to = Realloc(mb->theta_to, mb->ntheta + 1, char *);
			mb->theta_from[mb->ntheta] = Strdup(ds->data_prior.from_theta);
			mb->theta_to[mb->ntheta] = Strdup(ds->data_prior.to_theta);

			mb->theta[mb->ntheta] = ds->data_observations.log_offset_prec;
			mb->theta_map = Realloc(mb->theta_map, mb->ntheta + 1, map_func_tp *);
			mb->theta_map[mb->ntheta] = map_precision;
			mb->theta_map_arg = Realloc(mb->theta_map_arg, mb->ntheta + 1, void *);
			mb->theta_map_arg[mb->ntheta] = NULL;
			mb->ntheta++;
			ds->data_ntheta++;
		}
	}
		break;

	case L_STOCHVOL_LN:
	{
		/*
		 * get options related to the stochvol; the log-offset in the variance
		 */
		tmp = iniparser_getdouble(ini, inla_string_join(secname, "INITIAL"), 0.0);
		ds->data_fixed = iniparser_getboolean(ini, inla_string_join(secname, "FIXED"), 1);	/* yes, default fixed */
		if (!ds->data_fixed && mb->mode_use_mode) {
			tmp = mb->theta_file[mb->theta_counter_file++];
			if (mb->mode_fixed)
				ds->data_fixed = 1;
		}
		HYPER_NEW(ds->data_observations.stochvolln_c, tmp);
		if (mb->verbose) {
			printf("\t\tinitialise mean.offset[%g]\n", ds->data_observations.stochvolln_c[0][0]);
			printf("\t\tfixed=[%1d]\n", ds->data_fixed);
		}

		inla_read_prior(mb, ini, sec, &(ds->data_prior), "NORMAL", NULL);

		/*
		 * add theta 
		 */
		if (!ds->data_fixed) {
			mb->theta = Realloc(mb->theta, mb->ntheta + 1, double **);
			mb->theta_hyperid = Realloc(mb->theta_hyperid, mb->ntheta + 1, char *);
			mb->theta_hyperid[mb->ntheta] = ds->data_prior.hyperid;
			mb->theta_tag = Realloc(mb->theta_tag, mb->ntheta + 1, char *);
			mb->theta_tag_userscale = Realloc(mb->theta_tag_userscale, mb->ntheta + 1, char *);
			mb->theta_dir = Realloc(mb->theta_dir, mb->ntheta + 1, char *);
			mb->theta_tag[mb->ntheta] = inla_make_tag("Mean offset for stochvolln", mb->ds);
			mb->theta_tag_userscale[mb->ntheta] = inla_make_tag("Mean offset for stochvolln", mb->ds);
			GMRFLib_sprintf(&msg, "%s-parameter", secname);
			mb->theta_dir[mb->ntheta] = msg;

			mb->theta_from = Realloc(mb->theta_from, mb->ntheta + 1, char *);
			mb->theta_to = Realloc(mb->theta_to, mb->ntheta + 1, char *);
			mb->theta_from[mb->ntheta] = Strdup(ds->data_prior.from_theta);
			mb->theta_to[mb->ntheta] = Strdup(ds->data_prior.to_theta);

			mb->theta[mb->ntheta] = ds->data_observations.stochvolln_c;
			mb->theta_map = Realloc(mb->theta_map, mb->ntheta + 1, map_func_tp *);
			mb->theta_map[mb->ntheta] = map_identity;
			mb->theta_map_arg = Realloc(mb->theta_map_arg, mb->ntheta + 1, void *);
			mb->theta_map_arg[mb->ntheta] = NULL;
			mb->ntheta++;
			ds->data_ntheta++;
		}
	}
		break;

	case L_STOCHVOL_SN:
	{
		/*
		 * get options related to the stochvol_sn
		 */
		tmp = iniparser_getdouble(ini, inla_string_join(secname, "INITIAL0"), 0.00123456789);	/* yes! */
		ds->data_fixed0 = iniparser_getboolean(ini, inla_string_join(secname, "FIXED0"), 0);
		if (!ds->data_fixed0 && mb->mode_use_mode) {
			tmp = mb->theta_file[mb->theta_counter_file++];
			if (mb->mode_fixed)
				ds->data_fixed0 = 1;
		}
		HYPER_NEW(ds->data_observations.sn_skew, tmp);
		if (mb->verbose) {
			printf("\t\tinitialise intern skewness[%g]\n", ds->data_observations.sn_skew[0][0]);
			printf("\t\tfixed=[%1d]\n", ds->data_fixed0);
		}
		inla_read_prior0(mb, ini, sec, &(ds->data_prior0), "PCSN", NULL);

		/*
		 * add theta 
		 */
		if (!ds->data_fixed0) {
			mb->theta = Realloc(mb->theta, mb->ntheta + 1, double **);
			mb->theta_hyperid = Realloc(mb->theta_hyperid, mb->ntheta + 1, char *);
			mb->theta_hyperid[mb->ntheta] = ds->data_prior0.hyperid;
			mb->theta_tag = Realloc(mb->theta_tag, mb->ntheta + 1, char *);
			mb->theta_tag_userscale = Realloc(mb->theta_tag_userscale, mb->ntheta + 1, char *);
			mb->theta_dir = Realloc(mb->theta_dir, mb->ntheta + 1, char *);
			mb->theta_tag[mb->ntheta] = inla_make_tag("Intern skewness for stochvol_sn observations", mb->ds);
			mb->theta_tag_userscale[mb->ntheta] = inla_make_tag("Skewness for stochvol_sn observations", mb->ds);
			GMRFLib_sprintf(&msg, "%s-parameter0", secname);
			mb->theta_dir[mb->ntheta] = msg;

			mb->theta_from = Realloc(mb->theta_from, mb->ntheta + 1, char *);
			mb->theta_to = Realloc(mb->theta_to, mb->ntheta + 1, char *);
			mb->theta_from[mb->ntheta] = Strdup(ds->data_prior0.from_theta);
			mb->theta_to[mb->ntheta] = Strdup(ds->data_prior0.to_theta);
			mb->theta[mb->ntheta] = ds->data_observations.sn_skew;

			double *skewmax = Calloc(1, double);
			*skewmax = GMRFLib_SN_SKEWMAX;	       /* yes, this is correct */

			mb->theta_map = Realloc(mb->theta_map, mb->ntheta + 1, map_func_tp *);
			mb->theta_map[mb->ntheta] = map_phi;
			mb->theta_map_arg = Realloc(mb->theta_map_arg, mb->ntheta + 1, void *);
			mb->theta_map_arg[mb->ntheta] = (void *) skewmax;
			mb->ntheta++;
			ds->data_ntheta++;
		}

		double initial_value = 500.0;
		tmp = iniparser_getdouble(ini, inla_string_join(secname, "INITIAL1"), initial_value);
		ds->data_fixed1 = iniparser_getboolean(ini, inla_string_join(secname, "FIXED1"), 1);	/* yes, default fixed */
		if (!ds->data_fixed1 && mb->mode_use_mode) {
			tmp = mb->theta_file[mb->theta_counter_file++];
			if (mb->mode_fixed)
				ds->data_fixed1 = 1;
		}
		HYPER_NEW(ds->data_observations.log_offset_prec, tmp);
		if (mb->verbose) {
			printf("\t\tinitialise log_offset_prec[%g]\n", ds->data_observations.log_offset_prec[0][0]);
			printf("\t\tfixed=[%1d]\n", ds->data_fixed1);
		}
		inla_read_prior1(mb, ini, sec, &(ds->data_prior1), "LOGGAMMA", NULL);

		/*
		 * add theta 
		 */
		if (!ds->data_fixed1) {
			mb->theta = Realloc(mb->theta, mb->ntheta + 1, double **);
			mb->theta_hyperid = Realloc(mb->theta_hyperid, mb->ntheta + 1, char *);
			mb->theta_hyperid[mb->ntheta] = ds->data_prior1.hyperid;
			mb->theta_tag = Realloc(mb->theta_tag, mb->ntheta + 1, char *);
			mb->theta_tag_userscale = Realloc(mb->theta_tag_userscale, mb->ntheta + 1, char *);
			mb->theta_dir = Realloc(mb->theta_dir, mb->ntheta + 1, char *);
			mb->theta_tag[mb->ntheta] = inla_make_tag("Log offset precision for stochvol_sn", mb->ds);
			mb->theta_tag_userscale[mb->ntheta] = inla_make_tag("Offset precision for stochvol_sn", mb->ds);
			GMRFLib_sprintf(&msg, "%s-parameter1", secname);
			mb->theta_dir[mb->ntheta] = msg;

			mb->theta_from = Realloc(mb->theta_from, mb->ntheta + 1, char *);
			mb->theta_to = Realloc(mb->theta_to, mb->ntheta + 1, char *);
			mb->theta_from[mb->ntheta] = Strdup(ds->data_prior1.from_theta);
			mb->theta_to[mb->ntheta] = Strdup(ds->data_prior1.to_theta);

			mb->theta[mb->ntheta] = ds->data_observations.log_offset_prec;
			mb->theta_map = Realloc(mb->theta_map, mb->ntheta + 1, map_func_tp *);
			mb->theta_map[mb->ntheta] = map_precision;
			mb->theta_map_arg = Realloc(mb->theta_map_arg, mb->ntheta + 1, void *);
			mb->theta_map_arg[mb->ntheta] = NULL;
			mb->ntheta++;
			ds->data_ntheta++;
		}
	}
		break;

	case L_STOCHVOL_T:
	{
		/*
		 * get options related to the stochvol_t
		 */
		double initial_value = 3.0;

		tmp = iniparser_getdouble(ini, inla_string_join(secname, "INITIAL"), initial_value);
		ds->data_fixed = iniparser_getboolean(ini, inla_string_join(secname, "FIXED"), 0);
		if (!ds->data_fixed && mb->mode_use_mode) {
			tmp = mb->theta_file[mb->theta_counter_file++];
			if (mb->mode_fixed)
				ds->data_fixed = 1;
		}
		HYPER_NEW(ds->data_observations.dof_intern_svt, tmp);
		if (mb->verbose) {
			printf("\t\tinitialise dof_intern[%g]\n", ds->data_observations.dof_intern_svt[0][0]);
			printf("\t\tfixed=[%1d]\n", ds->data_fixed);
		}
		inla_read_prior(mb, ini, sec, &(ds->data_prior), "GAUSSIAN", NULL);

		/*
		 * add theta 
		 */
		if (!ds->data_fixed) {
			mb->theta = Realloc(mb->theta, mb->ntheta + 1, double **);
			mb->theta_hyperid = Realloc(mb->theta_hyperid, mb->ntheta + 1, char *);
			mb->theta_hyperid[mb->ntheta] = ds->data_prior.hyperid;
			mb->theta_tag = Realloc(mb->theta_tag, mb->ntheta + 1, char *);
			mb->theta_tag_userscale = Realloc(mb->theta_tag_userscale, mb->ntheta + 1, char *);
			mb->theta_dir = Realloc(mb->theta_dir, mb->ntheta + 1, char *);
			mb->theta_tag[mb->ntheta] = inla_make_tag("dof_intern for stochvol student-t", mb->ds);
			mb->theta_tag_userscale[mb->ntheta] = inla_make_tag("degrees of freedom for stochvol student-t", mb->ds);
			GMRFLib_sprintf(&msg, "%s-parameter", secname);
			mb->theta_dir[mb->ntheta] = msg;

			mb->theta_from = Realloc(mb->theta_from, mb->ntheta + 1, char *);
			mb->theta_to = Realloc(mb->theta_to, mb->ntheta + 1, char *);
			mb->theta_from[mb->ntheta] = Strdup(ds->data_prior.from_theta);
			mb->theta_to[mb->ntheta] = Strdup(ds->data_prior.to_theta);

			mb->theta[mb->ntheta] = ds->data_observations.dof_intern_svt;
			mb->theta_map = Realloc(mb->theta_map, mb->ntheta + 1, map_func_tp *);
			mb->theta_map[mb->ntheta] = map_dof;
			mb->theta_map_arg = Realloc(mb->theta_map_arg, mb->ntheta + 1, void *);
			mb->theta_map_arg[mb->ntheta] = NULL;
			mb->ntheta++;
			ds->data_ntheta++;
		}
	}
		break;

	case L_STOCHVOL_NIG:
	{
		/*
		 * get options related to the stochvol_nig
		 *
		 * first parameter is skew, second is shape
		 */
		double initial0 = 0.0, initial1 = 1.0;

		tmp = iniparser_getdouble(ini, inla_string_join(secname, "INITIAL0"), initial0);
		ds->data_fixed0 = iniparser_getboolean(ini, inla_string_join(secname, "FIXED0"), 0);
		if (!ds->data_fixed0 && mb->mode_use_mode) {
			tmp = mb->theta_file[mb->theta_counter_file++];
			if (mb->mode_fixed)
				ds->data_fixed0 = 1;
		}
		HYPER_NEW(ds->data_observations.skew_intern_svnig, tmp);
		if (mb->verbose) {
			printf("\t\tinitialise skew_intern_svnig[%g]\n", ds->data_observations.skew_intern_svnig[0][0]);
			printf("\t\tfixed=[%1d]\n", ds->data_fixed0);
		}
		inla_read_prior0(mb, ini, sec, &(ds->data_prior0), "GAUSSIAN", NULL);

		tmp = iniparser_getdouble(ini, inla_string_join(secname, "INITIAL1"), initial1);
		ds->data_fixed1 = iniparser_getboolean(ini, inla_string_join(secname, "FIXED1"), 0);
		if (!ds->data_fixed1 && mb->mode_use_mode) {
			tmp = mb->theta_file[mb->theta_counter_file++];
			if (mb->mode_fixed)
				ds->data_fixed1 = 1;
		}
		HYPER_NEW(ds->data_observations.shape_intern_svnig, tmp);
		if (mb->verbose) {
			printf("\t\tinitialise shape_intern_snvig[%g]\n", tmp);
			printf("\t\tfixed=[%1d]\n", ds->data_fixed1);
		}
		inla_read_prior1(mb, ini, sec, &(ds->data_prior1), "GAUSSIAN", NULL);

		/*
		 * add theta 
		 */
		if (!ds->data_fixed0) {
			mb->theta = Realloc(mb->theta, mb->ntheta + 1, double **);
			mb->theta_hyperid = Realloc(mb->theta_hyperid, mb->ntheta + 1, char *);
			mb->theta_hyperid[mb->ntheta] = ds->data_prior0.hyperid;
			mb->theta_tag = Realloc(mb->theta_tag, mb->ntheta + 1, char *);
			mb->theta_tag_userscale = Realloc(mb->theta_tag_userscale, mb->ntheta + 1, char *);
			mb->theta_dir = Realloc(mb->theta_dir, mb->ntheta + 1, char *);
			mb->theta_tag[mb->ntheta] = inla_make_tag("skewness_param_intern for stochvol-nig", mb->ds);
			mb->theta_tag_userscale[mb->ntheta] = inla_make_tag("skewness parameter for stochvol-nig", mb->ds);
			GMRFLib_sprintf(&msg, "%s-parameter0", secname);
			mb->theta_dir[mb->ntheta] = msg;

			mb->theta_from = Realloc(mb->theta_from, mb->ntheta + 1, char *);
			mb->theta_to = Realloc(mb->theta_to, mb->ntheta + 1, char *);
			mb->theta_from[mb->ntheta] = Strdup(ds->data_prior0.from_theta);
			mb->theta_to[mb->ntheta] = Strdup(ds->data_prior0.to_theta);

			mb->theta[mb->ntheta] = ds->data_observations.skew_intern_svnig;
			mb->theta_map = Realloc(mb->theta_map, mb->ntheta + 1, map_func_tp *);
			mb->theta_map[mb->ntheta] = map_identity;
			mb->theta_map_arg = Realloc(mb->theta_map_arg, mb->ntheta + 1, void *);
			mb->theta_map_arg[mb->ntheta] = NULL;
			mb->ntheta++;
			ds->data_ntheta++;
		}
		if (!ds->data_fixed1) {
			mb->theta = Realloc(mb->theta, mb->ntheta + 1, double **);
			mb->theta_hyperid = Realloc(mb->theta_hyperid, mb->ntheta + 1, char *);
			mb->theta_hyperid[mb->ntheta] = ds->data_prior1.hyperid;
			mb->theta_tag = Realloc(mb->theta_tag, mb->ntheta + 1, char *);
			mb->theta_tag_userscale = Realloc(mb->theta_tag_userscale, mb->ntheta + 1, char *);
			mb->theta_dir = Realloc(mb->theta_dir, mb->ntheta + 1, char *);
			mb->theta_tag[mb->ntheta] = inla_make_tag("shape_param_intern for stochvol-nig", mb->ds);
			mb->theta_tag_userscale[mb->ntheta] = inla_make_tag("shape parameter for stochvol-nig", mb->ds);
			GMRFLib_sprintf(&msg, "%s-parameter1", secname);
			mb->theta_dir[mb->ntheta] = msg;

			mb->theta_from = Realloc(mb->theta_from, mb->ntheta + 1, char *);
			mb->theta_to = Realloc(mb->theta_to, mb->ntheta + 1, char *);
			mb->theta_from[mb->ntheta] = Strdup(ds->data_prior1.from_theta);
			mb->theta_to[mb->ntheta] = Strdup(ds->data_prior1.to_theta);

			mb->theta[mb->ntheta] = ds->data_observations.shape_intern_svnig;
			mb->theta_map = Realloc(mb->theta_map, mb->ntheta + 1, map_func_tp *);
			mb->theta_map[mb->ntheta] = map_shape_svnig;
			mb->theta_map_arg = Realloc(mb->theta_map_arg, mb->ntheta + 1, void *);
			mb->theta_map_arg[mb->ntheta] = NULL;
			mb->ntheta++;
			ds->data_ntheta++;
		}
		ds->data_ntheta = (ds->data_fixed0 ? 0 : 1) + (ds->data_fixed1 ? 0 : 1);
	}
		break;

	case L_WEIBULL:
	{
		/*
		 * get options related to the weibull
		 */
		double initial_value = 0.0;

		GMRFLib_ASSERT(ds->variant == 0 || ds->variant == 1, GMRFLib_EPARAMETER);
		tmp = iniparser_getdouble(ini, inla_string_join(secname, "INITIAL"), initial_value);
		ds->data_fixed = iniparser_getboolean(ini, inla_string_join(secname, "FIXED"), 0);
		if (!ds->data_fixed && mb->mode_use_mode) {
			tmp = mb->theta_file[mb->theta_counter_file++];
			if (mb->mode_fixed)
				ds->data_fixed = 1;
		}
		HYPER_NEW(ds->data_observations.alpha_intern, tmp);
		if (mb->verbose) {
			printf("\t\tinitialise alpha_intern[%g]\n", ds->data_observations.alpha_intern[0][0]);
			printf("\t\tfixed=[%1d]\n", ds->data_fixed);
		}
		inla_read_prior(mb, ini, sec, &(ds->data_prior), "LOGGAMMA-ALPHA", NULL);

		/*
		 * add theta 
		 */
		if (!ds->data_fixed) {
			mb->theta = Realloc(mb->theta, mb->ntheta + 1, double **);
			mb->theta_hyperid = Realloc(mb->theta_hyperid, mb->ntheta + 1, char *);
			mb->theta_hyperid[mb->ntheta] = ds->data_prior.hyperid;
			mb->theta_tag = Realloc(mb->theta_tag, mb->ntheta + 1, char *);
			mb->theta_tag_userscale = Realloc(mb->theta_tag_userscale, mb->ntheta + 1, char *);
			mb->theta_dir = Realloc(mb->theta_dir, mb->ntheta + 1, char *);
			mb->theta_tag[mb->ntheta] = inla_make_tag("alpha_intern for weibull", mb->ds);
			mb->theta_tag_userscale[mb->ntheta] = inla_make_tag("alpha parameter for weibull", mb->ds);
			GMRFLib_sprintf(&msg, "%s-parameter", secname);
			mb->theta_dir[mb->ntheta] = msg;

			mb->theta_from = Realloc(mb->theta_from, mb->ntheta + 1, char *);
			mb->theta_to = Realloc(mb->theta_to, mb->ntheta + 1, char *);
			mb->theta_from[mb->ntheta] = Strdup(ds->data_prior.from_theta);
			mb->theta_to[mb->ntheta] = Strdup(ds->data_prior.to_theta);

			mb->theta[mb->ntheta] = ds->data_observations.alpha_intern;
			mb->theta_map = Realloc(mb->theta_map, mb->ntheta + 1, map_func_tp *);
			mb->theta_map[mb->ntheta] = map_alpha_weibull;	/* alpha = exp(alpha.intern) */
			mb->theta_map_arg = Realloc(mb->theta_map_arg, mb->ntheta + 1, void *);
			mb->theta_map_arg[mb->ntheta] = NULL;
			mb->ntheta++;
			ds->data_ntheta++;
		}
	}
		break;

	case L_WEIBULLSURV:
	{
		for (i = 0; i < CURE_MAXTHETA + 1; i++) {
			GMRFLib_sprintf(&ctmp, "FIXED%1d", i);
			iniparser_getstring(ini, inla_string_join(secname, ctmp), NULL);

			GMRFLib_sprintf(&ctmp, "INITIAL%1d", i);
			iniparser_getstring(ini, inla_string_join(secname, ctmp), NULL);

			GMRFLib_sprintf(&ctmp, "PRIOR%1d", i);
			iniparser_getstring(ini, inla_string_join(secname, ctmp), NULL);

			GMRFLib_sprintf(&ctmp, "HYPERID%1d", i);
			iniparser_getstring(ini, inla_string_join(secname, ctmp), NULL);

			GMRFLib_sprintf(&ctmp, "PARAMETERS%1d", i);
			iniparser_getstring(ini, inla_string_join(secname, ctmp), NULL);

			GMRFLib_sprintf(&ctmp, "to.theta%1d", i);
			iniparser_getstring(ini, inla_string_join(secname, ctmp), NULL);

			GMRFLib_sprintf(&ctmp, "from.theta%1d", i);
			iniparser_getstring(ini, inla_string_join(secname, ctmp), NULL);
		}

		ds->data_nfixed = Calloc(CURE_MAXTHETA + 1, int);
		ds->data_nprior = Calloc(CURE_MAXTHETA + 1, Prior_tp);
		ds->data_observations.cure_beta = Calloc(CURE_MAXTHETA, double **);

		double initial_value = 0.0;

		GMRFLib_ASSERT(ds->variant == 0 || ds->variant == 1, GMRFLib_EPARAMETER);
		tmp = iniparser_getdouble(ini, inla_string_join(secname, "INITIAL0"), initial_value);
		ds->data_nfixed[0] = iniparser_getboolean(ini, inla_string_join(secname, "FIXED0"), 0);
		if (!ds->data_nfixed[0] && mb->mode_use_mode) {
			tmp = mb->theta_file[mb->theta_counter_file++];
			if (mb->mode_fixed)
				ds->data_nfixed[0] = 1;
		}
		HYPER_NEW(ds->data_observations.alpha_intern, tmp);
		if (mb->verbose) {
			printf("\t\tinitialise alpha_intern[%g]\n", ds->data_observations.alpha_intern[0][0]);
			printf("\t\tfixed=[%1d]\n", ds->data_fixed);
		}
		inla_read_priorN(mb, ini, sec, &(ds->data_nprior[0]), "LOGGAMMA-ALPHA", 0, NULL);

		if (!ds->data_nfixed[0]) {
			mb->theta = Realloc(mb->theta, mb->ntheta + 1, double **);
			mb->theta_hyperid = Realloc(mb->theta_hyperid, mb->ntheta + 1, char *);
			mb->theta_hyperid[mb->ntheta] = ds->data_nprior[0].hyperid;
			mb->theta_tag = Realloc(mb->theta_tag, mb->ntheta + 1, char *);
			mb->theta_tag_userscale = Realloc(mb->theta_tag_userscale, mb->ntheta + 1, char *);
			mb->theta_dir = Realloc(mb->theta_dir, mb->ntheta + 1, char *);
			mb->theta_tag[mb->ntheta] = inla_make_tag("alpha_intern for weibullsurv", mb->ds);
			mb->theta_tag_userscale[mb->ntheta] = inla_make_tag("alpha parameter for weibullsurv", mb->ds);
			GMRFLib_sprintf(&msg, "%s-parameter", secname);
			mb->theta_dir[mb->ntheta] = msg;

			mb->theta_from = Realloc(mb->theta_from, mb->ntheta + 1, char *);
			mb->theta_to = Realloc(mb->theta_to, mb->ntheta + 1, char *);
			mb->theta_from[mb->ntheta] = Strdup(ds->data_nprior[0].from_theta);
			mb->theta_to[mb->ntheta] = Strdup(ds->data_nprior[0].to_theta);

			mb->theta[mb->ntheta] = ds->data_observations.alpha_intern;
			mb->theta_map = Realloc(mb->theta_map, mb->ntheta + 1, map_func_tp *);
			mb->theta_map[mb->ntheta] = map_alpha_weibull;	/* alpha = exp(alpha.intern) */
			mb->theta_map_arg = Realloc(mb->theta_map_arg, mb->ntheta + 1, void *);
			mb->theta_map_arg[mb->ntheta] = NULL;
			mb->ntheta++;
			ds->data_ntheta++;
		}

		for (i = 1; i < ds->data_observations.cure_ncov + 1; i++) {
			GMRFLib_sprintf(&ctmp, "INITIAL%1d", i);
			tmp = iniparser_getdouble(ini, inla_string_join(secname, ctmp), 0.0);	/* YES! */

			GMRFLib_sprintf(&ctmp, "FIXED%1d", i);
			ds->data_nfixed[i] = iniparser_getboolean(ini, inla_string_join(secname, ctmp), 0);
			if (!ds->data_nfixed[i] && mb->mode_use_mode) {
				tmp = mb->theta_file[mb->theta_counter_file++];
				if (mb->mode_fixed)
					ds->data_nfixed[i] = 1;
			}

			HYPER_NEW(ds->data_observations.cure_beta[i - 1], tmp);
			if (mb->verbose) {
				printf("\t\tbeta[%1d] = %g\n", i, ds->data_observations.cure_beta[i - 1][0][0]);
				printf("\t\tfixed[%1d] = %1d\n", i, ds->data_nfixed[i]);
			}
			inla_read_priorN(mb, ini, sec, &(ds->data_nprior[i]), "GAUSSIAN-std", i, NULL);

			if (!ds->data_nfixed[i]) {
				mb->theta = Realloc(mb->theta, mb->ntheta + 1, double **);
				mb->theta_hyperid = Realloc(mb->theta_hyperid, mb->ntheta + 1, char *);
				mb->theta_hyperid[mb->ntheta] = ds->data_nprior[i].hyperid;
				mb->theta_tag = Realloc(mb->theta_tag, mb->ntheta + 1, char *);
				mb->theta_tag_userscale = Realloc(mb->theta_tag_userscale, mb->ntheta + 1, char *);
				mb->theta_dir = Realloc(mb->theta_dir, mb->ntheta + 1, char *);
				GMRFLib_sprintf(&ctmp, "beta%1d for Weibull-Cure", i);

				mb->theta_tag[mb->ntheta] = inla_make_tag(ctmp, mb->ds);
				mb->theta_tag_userscale[mb->ntheta] = inla_make_tag(ctmp, mb->ds);
				GMRFLib_sprintf(&msg, "%s-parameter%1d", secname, i);
				mb->theta_dir[mb->ntheta] = msg;

				mb->theta_from = Realloc(mb->theta_from, mb->ntheta + 1, char *);
				mb->theta_to = Realloc(mb->theta_to, mb->ntheta + 1, char *);
				mb->theta_from[mb->ntheta] = Strdup(ds->data_nprior[i].from_theta);
				mb->theta_to[mb->ntheta] = Strdup(ds->data_nprior[i].to_theta);

				mb->theta[mb->ntheta] = ds->data_observations.cure_beta[i - 1];
				mb->theta_map = Realloc(mb->theta_map, mb->ntheta + 1, map_func_tp *);

				mb->theta_map[mb->ntheta] = map_identity;
				mb->theta_map_arg = Realloc(mb->theta_map_arg, mb->ntheta + 1, void *);
				mb->theta_map_arg[mb->ntheta] = NULL;

				mb->ntheta++;
				ds->data_ntheta++;
			}
		}

		// all the remaining ones are fixed
		for (i = 1 + ds->data_observations.cure_ncov; i < 1 + CURE_MAXTHETA; i++) {
			ds->data_nfixed[i] = 1;
		}
	}
		break;

	case L_GOMPERTZ:
	{
		/*
		 * get options related to the gompertz
		 */
		double initial_value = 0.0;

		GMRFLib_ASSERT(ds->variant == 0 || ds->variant == 1, GMRFLib_EPARAMETER);
		tmp = iniparser_getdouble(ini, inla_string_join(secname, "INITIAL"), initial_value);
		ds->data_fixed = iniparser_getboolean(ini, inla_string_join(secname, "FIXED"), 0);
		if (!ds->data_fixed && mb->mode_use_mode) {
			tmp = mb->theta_file[mb->theta_counter_file++];
			if (mb->mode_fixed)
				ds->data_fixed = 1;
		}
		HYPER_NEW(ds->data_observations.alpha_intern, tmp);
		if (mb->verbose) {
			printf("\t\tinitialise alpha_intern[%g]\n", ds->data_observations.alpha_intern[0][0]);
			printf("\t\tfixed=[%1d]\n", ds->data_fixed);
		}
		inla_read_prior(mb, ini, sec, &(ds->data_prior), "LOGGAMMA", NULL);

		/*
		 * add theta 
		 */
		if (!ds->data_fixed) {
			mb->theta = Realloc(mb->theta, mb->ntheta + 1, double **);
			mb->theta_hyperid = Realloc(mb->theta_hyperid, mb->ntheta + 1, char *);
			mb->theta_hyperid[mb->ntheta] = ds->data_prior.hyperid;
			mb->theta_tag = Realloc(mb->theta_tag, mb->ntheta + 1, char *);
			mb->theta_tag_userscale = Realloc(mb->theta_tag_userscale, mb->ntheta + 1, char *);
			mb->theta_dir = Realloc(mb->theta_dir, mb->ntheta + 1, char *);
			if (ds->data_id == L_GOMPERTZ) {
				mb->theta_tag[mb->ntheta] = inla_make_tag("alpha_intern for Gompertz", mb->ds);
				mb->theta_tag_userscale[mb->ntheta] = inla_make_tag("alpha parameter for Gompertz", mb->ds);
			} else {
				mb->theta_tag[mb->ntheta] = inla_make_tag("alpha_intern for Gompertz-surv", mb->ds);
				mb->theta_tag_userscale[mb->ntheta] = inla_make_tag("alpha parameter for Gompertz-surv", mb->ds);
			}
			GMRFLib_sprintf(&msg, "%s-parameter", secname);
			mb->theta_dir[mb->ntheta] = msg;

			mb->theta_from = Realloc(mb->theta_from, mb->ntheta + 1, char *);
			mb->theta_to = Realloc(mb->theta_to, mb->ntheta + 1, char *);
			mb->theta_from[mb->ntheta] = Strdup(ds->data_prior.from_theta);
			mb->theta_to[mb->ntheta] = Strdup(ds->data_prior.to_theta);

			mb->theta[mb->ntheta] = ds->data_observations.alpha_intern;
			mb->theta_map = Realloc(mb->theta_map, mb->ntheta + 1, map_func_tp *);
			mb->theta_map[mb->ntheta] = map_alpha_gompertz;	/* alpha = exp(alpha.intern) */
			mb->theta_map_arg = Realloc(mb->theta_map_arg, mb->ntheta + 1, void *);
			mb->theta_map_arg[mb->ntheta] = NULL;
			mb->ntheta++;
			ds->data_ntheta++;
		}
	}
		break;

	case L_GOMPERTZSURV:
	{
		for (i = 0; i < CURE_MAXTHETA + 1; i++) {
			GMRFLib_sprintf(&ctmp, "FIXED%1d", i);
			iniparser_getstring(ini, inla_string_join(secname, ctmp), NULL);

			GMRFLib_sprintf(&ctmp, "INITIAL%1d", i);
			iniparser_getstring(ini, inla_string_join(secname, ctmp), NULL);

			GMRFLib_sprintf(&ctmp, "PRIOR%1d", i);
			iniparser_getstring(ini, inla_string_join(secname, ctmp), NULL);

			GMRFLib_sprintf(&ctmp, "HYPERID%1d", i);
			iniparser_getstring(ini, inla_string_join(secname, ctmp), NULL);

			GMRFLib_sprintf(&ctmp, "PARAMETERS%1d", i);
			iniparser_getstring(ini, inla_string_join(secname, ctmp), NULL);

			GMRFLib_sprintf(&ctmp, "to.theta%1d", i);
			iniparser_getstring(ini, inla_string_join(secname, ctmp), NULL);

			GMRFLib_sprintf(&ctmp, "from.theta%1d", i);
			iniparser_getstring(ini, inla_string_join(secname, ctmp), NULL);
		}

		ds->data_nfixed = Calloc(CURE_MAXTHETA + 1, int);
		ds->data_nprior = Calloc(CURE_MAXTHETA + 1, Prior_tp);
		ds->data_observations.cure_beta = Calloc(CURE_MAXTHETA, double **);

		double initial_value = 0.0;
		GMRFLib_ASSERT(ds->variant == 0 || ds->variant == 1, GMRFLib_EPARAMETER);
		tmp = iniparser_getdouble(ini, inla_string_join(secname, "INITIAL0"), initial_value);
		ds->data_nfixed[0] = iniparser_getboolean(ini, inla_string_join(secname, "FIXED0"), 0);
		if (!ds->data_nfixed[0] && mb->mode_use_mode) {
			tmp = mb->theta_file[mb->theta_counter_file++];
			if (mb->mode_fixed)
				ds->data_nfixed[0] = 1;
		}
		HYPER_NEW(ds->data_observations.alpha_intern, tmp);
		if (mb->verbose) {
			printf("\t\tinitialise alpha_intern[%g]\n", ds->data_observations.alpha_intern[0][0]);
			printf("\t\tfixed=[%1d]\n", ds->data_nfixed[0]);
		}
		inla_read_priorN(mb, ini, sec, &(ds->data_nprior[0]), "LOGGAMMA", 0, NULL);

		/*
		 * add theta 
		 */
		if (!ds->data_nfixed[0]) {
			mb->theta = Realloc(mb->theta, mb->ntheta + 1, double **);
			mb->theta_hyperid = Realloc(mb->theta_hyperid, mb->ntheta + 1, char *);
			mb->theta_hyperid[mb->ntheta] = ds->data_nprior[0].hyperid;
			mb->theta_tag = Realloc(mb->theta_tag, mb->ntheta + 1, char *);
			mb->theta_tag_userscale = Realloc(mb->theta_tag_userscale, mb->ntheta + 1, char *);
			mb->theta_dir = Realloc(mb->theta_dir, mb->ntheta + 1, char *);
			mb->theta_tag[mb->ntheta] = inla_make_tag("alpha_intern for Gompertz-surv", mb->ds);
			mb->theta_tag_userscale[mb->ntheta] = inla_make_tag("alpha parameter for Gompertz-surv", mb->ds);
			GMRFLib_sprintf(&msg, "%s-parameter", secname);
			mb->theta_dir[mb->ntheta] = msg;

			mb->theta_from = Realloc(mb->theta_from, mb->ntheta + 1, char *);
			mb->theta_to = Realloc(mb->theta_to, mb->ntheta + 1, char *);
			mb->theta_from[mb->ntheta] = Strdup(ds->data_nprior[0].from_theta);
			mb->theta_to[mb->ntheta] = Strdup(ds->data_nprior[0].to_theta);

			mb->theta[mb->ntheta] = ds->data_observations.alpha_intern;
			mb->theta_map = Realloc(mb->theta_map, mb->ntheta + 1, map_func_tp *);
			mb->theta_map[mb->ntheta] = map_alpha_gompertz;	/* alpha = exp(alpha.intern) */
			mb->theta_map_arg = Realloc(mb->theta_map_arg, mb->ntheta + 1, void *);
			mb->theta_map_arg[mb->ntheta] = NULL;
			mb->ntheta++;
			ds->data_ntheta++;
		}

		for (i = 1; i < ds->data_observations.cure_ncov + 1; i++) {
			GMRFLib_sprintf(&ctmp, "INITIAL%1d", i);
			tmp = iniparser_getdouble(ini, inla_string_join(secname, ctmp), 0.0);	/* YES! */

			GMRFLib_sprintf(&ctmp, "FIXED%1d", i);
			ds->data_nfixed[i] = iniparser_getboolean(ini, inla_string_join(secname, ctmp), 0);
			if (!ds->data_nfixed[i] && mb->mode_use_mode) {
				tmp = mb->theta_file[mb->theta_counter_file++];
				if (mb->mode_fixed)
					ds->data_nfixed[i] = 1;
			}

			HYPER_NEW(ds->data_observations.cure_beta[i - 1], tmp);
			if (mb->verbose) {
				printf("\t\tbeta[%1d] = %g\n", i, ds->data_observations.cure_beta[i - 1][0][0]);
				printf("\t\tfixed[%1d] = %1d\n", i, ds->data_nfixed[i]);
			}
			inla_read_priorN(mb, ini, sec, &(ds->data_nprior[i]), "GAUSSIAN-std", i, NULL);

			if (!ds->data_nfixed[i]) {
				mb->theta = Realloc(mb->theta, mb->ntheta + 1, double **);
				mb->theta_hyperid = Realloc(mb->theta_hyperid, mb->ntheta + 1, char *);
				mb->theta_hyperid[mb->ntheta] = ds->data_nprior[i].hyperid;
				mb->theta_tag = Realloc(mb->theta_tag, mb->ntheta + 1, char *);
				mb->theta_tag_userscale = Realloc(mb->theta_tag_userscale, mb->ntheta + 1, char *);
				mb->theta_dir = Realloc(mb->theta_dir, mb->ntheta + 1, char *);
				GMRFLib_sprintf(&ctmp, "beta%1d for Gompertz-Cure", i);

				mb->theta_tag[mb->ntheta] = inla_make_tag(ctmp, mb->ds);
				mb->theta_tag_userscale[mb->ntheta] = inla_make_tag(ctmp, mb->ds);
				GMRFLib_sprintf(&msg, "%s-parameter%1d", secname, i);
				mb->theta_dir[mb->ntheta] = msg;

				mb->theta_from = Realloc(mb->theta_from, mb->ntheta + 1, char *);
				mb->theta_to = Realloc(mb->theta_to, mb->ntheta + 1, char *);
				mb->theta_from[mb->ntheta] = Strdup(ds->data_nprior[i].from_theta);
				mb->theta_to[mb->ntheta] = Strdup(ds->data_nprior[i].to_theta);

				mb->theta[mb->ntheta] = ds->data_observations.cure_beta[i - 1];
				mb->theta_map = Realloc(mb->theta_map, mb->ntheta + 1, map_func_tp *);

				mb->theta_map[mb->ntheta] = map_identity;
				mb->theta_map_arg = Realloc(mb->theta_map_arg, mb->ntheta + 1, void *);
				mb->theta_map_arg[mb->ntheta] = NULL;

				mb->ntheta++;
				ds->data_ntheta++;
			}
		}

		// all the remaining ones are fixed
		for (i = 1 + ds->data_observations.cure_ncov; i < 1 + CURE_MAXTHETA; i++) {
			ds->data_nfixed[i] = 1;
		}
	}
		break;

	case L_POISSON_SPECIAL1:
	{
		double initial_value = -1.0;

		tmp = iniparser_getdouble(ini, inla_string_join(secname, "INITIAL"), initial_value);
		ds->data_fixed = iniparser_getboolean(ini, inla_string_join(secname, "FIXED"), 0);
		if (!ds->data_fixed && mb->mode_use_mode) {
			tmp = mb->theta_file[mb->theta_counter_file++];
			if (mb->mode_fixed)
				ds->data_fixed = 1;
		}
		HYPER_NEW(ds->data_observations.prob_intern, tmp);
		if (mb->verbose) {
			printf("\t\tinitialise prob_intern[%g]\n", ds->data_observations.prob_intern[0][0]);
			printf("\t\tfixed=[%1d]\n", ds->data_fixed);
		}
		inla_read_prior(mb, ini, sec, &(ds->data_prior), "GAUSSIAN-std", NULL);

		/*
		 * add theta 
		 */
		if (!ds->data_fixed) {
			mb->theta = Realloc(mb->theta, mb->ntheta + 1, double **);
			mb->theta_hyperid = Realloc(mb->theta_hyperid, mb->ntheta + 1, char *);
			mb->theta_hyperid[mb->ntheta] = ds->data_prior.hyperid;
			mb->theta_tag = Realloc(mb->theta_tag, mb->ntheta + 1, char *);
			mb->theta_tag_userscale = Realloc(mb->theta_tag_userscale, mb->ntheta + 1, char *);
			mb->theta_dir = Realloc(mb->theta_dir, mb->ntheta + 1, char *);
			mb->theta_tag[mb->ntheta] = inla_make_tag("intern one-probability parameter for poisson.special1", mb->ds);
			mb->theta_tag_userscale[mb->ntheta] = inla_make_tag("one-probability parameter for poisson.special1", mb->ds);
			GMRFLib_sprintf(&msg, "%s-parameter", secname);
			mb->theta_dir[mb->ntheta] = msg;

			mb->theta_from = Realloc(mb->theta_from, mb->ntheta + 1, char *);
			mb->theta_to = Realloc(mb->theta_to, mb->ntheta + 1, char *);
			mb->theta_from[mb->ntheta] = Strdup(ds->data_prior.from_theta);
			mb->theta_to[mb->ntheta] = Strdup(ds->data_prior.to_theta);

			mb->theta[mb->ntheta] = ds->data_observations.prob_intern;
			mb->theta_map = Realloc(mb->theta_map, mb->ntheta + 1, map_func_tp *);
			mb->theta_map[mb->ntheta] = map_probability;
			mb->theta_map_arg = Realloc(mb->theta_map_arg, mb->ntheta + 1, void *);
			mb->theta_map_arg[mb->ntheta] = NULL;
			mb->ntheta++;
			ds->data_ntheta++;
		}
	}
		break;

	case L_ZEROINFLATEDPOISSON0:
	case L_ZEROINFLATEDPOISSON1:
	case L_ZEROINFLATEDCENPOISSON0:
	case L_ZEROINFLATEDCENPOISSON1:
	{
		/*
		 * get options related to the zeroinflatedpoisson0/1
		 */
		double initial_value = -1.0;

		tmp = iniparser_getdouble(ini, inla_string_join(secname, "INITIAL"), initial_value);
		ds->data_fixed = iniparser_getboolean(ini, inla_string_join(secname, "FIXED"), 0);
		if (!ds->data_fixed && mb->mode_use_mode) {
			tmp = mb->theta_file[mb->theta_counter_file++];
			if (mb->mode_fixed)
				ds->data_fixed = 1;
		}
		HYPER_NEW(ds->data_observations.prob_intern, tmp);
		if (mb->verbose) {
			printf("\t\tinitialise prob_intern[%g]\n", ds->data_observations.prob_intern[0][0]);
			printf("\t\tfixed=[%1d]\n", ds->data_fixed);
		}
		inla_read_prior(mb, ini, sec, &(ds->data_prior), "GAUSSIAN-std", NULL);

		if (ds->data_id == L_ZEROINFLATEDCENPOISSON0 || ds->data_id == L_ZEROINFLATEDCENPOISSON1) {
			ds->data_observations.cenpoisson_interval = Calloc(2, double);
			ctmp = Strdup(iniparser_getstring(ini, inla_string_join(secname, "CENPOISSON.I"), NULL));
			if (inla_sread_doubles(ds->data_observations.cenpoisson_interval, 2, ctmp) == INLA_FAIL) {
				inla_error_field_is_void(__GMRFLib_FuncName, secname, "CENPOISSON.I", ctmp);
			}

			if ((int) ds->data_observations.cenpoisson_interval[0] == 0) {
				GMRFLib_sprintf(&ctmp,
						"cenpoisson invalid censor-interval [%g, %g]\n",
						ds->data_observations.cenpoisson_interval[0], ds->data_observations.cenpoisson_interval[1]);
				inla_error_general(ctmp);
				exit(1);
			}

			if (mb->verbose) {
				printf("\t\tcenpoisson censor-interval = [%g, %g]\n",
				       ds->data_observations.cenpoisson_interval[0], ds->data_observations.cenpoisson_interval[1]);
			}

			Free(ctmp);
		}

		/*
		 * add theta 
		 */
		if (!ds->data_fixed) {
			mb->theta = Realloc(mb->theta, mb->ntheta + 1, double **);
			mb->theta_hyperid = Realloc(mb->theta_hyperid, mb->ntheta + 1, char *);
			mb->theta_hyperid[mb->ntheta] = ds->data_prior.hyperid;
			mb->theta_tag = Realloc(mb->theta_tag, mb->ntheta + 1, char *);
			mb->theta_tag_userscale = Realloc(mb->theta_tag_userscale, mb->ntheta + 1, char *);
			mb->theta_dir = Realloc(mb->theta_dir, mb->ntheta + 1, char *);
			if (ds->data_id == L_ZEROINFLATEDPOISSON0) {
				mb->theta_tag[mb->ntheta] = inla_make_tag("intern zero-probability parameter for zero-inflated poisson_0", mb->ds);
				mb->theta_tag_userscale[mb->ntheta] =
				    inla_make_tag("zero-probability parameter for zero-inflated poisson_0", mb->ds);
			} else if (ds->data_id == L_ZEROINFLATEDPOISSON1) {
				mb->theta_tag[mb->ntheta] = inla_make_tag("intern zero-probability parameter for zero-inflated poisson_1", mb->ds);
				mb->theta_tag_userscale[mb->ntheta] =
				    inla_make_tag("zero-probability parameter for zero-inflated poisson_1", mb->ds);
			} else if (ds->data_id == L_ZEROINFLATEDCENPOISSON0) {
				mb->theta_tag[mb->ntheta] =
				    inla_make_tag("intern zero-probability parameter for zero-inflated cenpoisson_0", mb->ds);
				mb->theta_tag_userscale[mb->ntheta] =
				    inla_make_tag("zero-probability parameter for zero-inflated cenpoisson_0", mb->ds);
			} else if (ds->data_id == L_ZEROINFLATEDCENPOISSON1) {
				mb->theta_tag[mb->ntheta] =
				    inla_make_tag("intern zero-probability parameter for zero-inflated cenpoisson_1", mb->ds);
				mb->theta_tag_userscale[mb->ntheta] =
				    inla_make_tag("zero-probability parameter for zero-inflated cenpoisson_1", mb->ds);
			} else {
				GMRFLib_ASSERT(0 == 1, GMRFLib_ESNH);
			}
			GMRFLib_sprintf(&msg, "%s-parameter", secname);
			mb->theta_dir[mb->ntheta] = msg;

			mb->theta_from = Realloc(mb->theta_from, mb->ntheta + 1, char *);
			mb->theta_to = Realloc(mb->theta_to, mb->ntheta + 1, char *);
			mb->theta_from[mb->ntheta] = Strdup(ds->data_prior.from_theta);
			mb->theta_to[mb->ntheta] = Strdup(ds->data_prior.to_theta);

			mb->theta[mb->ntheta] = ds->data_observations.prob_intern;
			mb->theta_map = Realloc(mb->theta_map, mb->ntheta + 1, map_func_tp *);
			mb->theta_map[mb->ntheta] = map_probability;
			mb->theta_map_arg = Realloc(mb->theta_map_arg, mb->ntheta + 1, void *);
			mb->theta_map_arg[mb->ntheta] = NULL;
			mb->ntheta++;
			ds->data_ntheta++;
		}
	}
		break;

	case L_ZEROINFLATEDPOISSON2:
	{
		/*
		 * get options related to the zeroinflatedpoisson2
		 */
		double initial_value = log(2.0);

		tmp = iniparser_getdouble(ini, inla_string_join(secname, "INITIAL"), initial_value);
		ds->data_fixed = iniparser_getboolean(ini, inla_string_join(secname, "FIXED"), 0);
		if (!ds->data_fixed && mb->mode_use_mode) {
			tmp = mb->theta_file[mb->theta_counter_file++];
			if (mb->mode_fixed)
				ds->data_fixed = 1;
		}
		HYPER_NEW(ds->data_observations.zeroinflated_alpha_intern, tmp);
		if (mb->verbose) {
			printf("\t\tinitialise alpha_intern[%g]\n", ds->data_observations.zeroinflated_alpha_intern[0][0]);
			printf("\t\tfixed=[%1d]\n", ds->data_fixed);
		}
		inla_read_prior(mb, ini, sec, &(ds->data_prior), "GAUSSIAN-std", NULL);
		/*
		 * add theta 
		 */
		if (!ds->data_fixed) {
			mb->theta = Realloc(mb->theta, mb->ntheta + 1, double **);
			mb->theta_hyperid = Realloc(mb->theta_hyperid, mb->ntheta + 1, char *);
			mb->theta_hyperid[mb->ntheta] = ds->data_prior.hyperid;
			mb->theta_tag = Realloc(mb->theta_tag, mb->ntheta + 1, char *);
			mb->theta_tag_userscale = Realloc(mb->theta_tag_userscale, mb->ntheta + 1, char *);
			mb->theta_dir = Realloc(mb->theta_dir, mb->ntheta + 1, char *);
			mb->theta_tag[mb->ntheta] = inla_make_tag("intern zero-probability parameter for zero-inflated poisson_2", mb->ds);
			mb->theta_tag_userscale[mb->ntheta] = inla_make_tag("zero-probability parameter for zero-inflated poisson_2", mb->ds);
			GMRFLib_sprintf(&msg, "%s-parameter", secname);
			mb->theta_dir[mb->ntheta] = msg;

			mb->theta_from = Realloc(mb->theta_from, mb->ntheta + 1, char *);
			mb->theta_to = Realloc(mb->theta_to, mb->ntheta + 1, char *);
			mb->theta_from[mb->ntheta] = Strdup(ds->data_prior.from_theta);
			mb->theta_to[mb->ntheta] = Strdup(ds->data_prior.to_theta);

			mb->theta[mb->ntheta] = ds->data_observations.zeroinflated_alpha_intern;
			mb->theta_map = Realloc(mb->theta_map, mb->ntheta + 1, map_func_tp *);
			mb->theta_map[mb->ntheta] = map_exp;
			mb->theta_map_arg = Realloc(mb->theta_map_arg, mb->ntheta + 1, void *);
			mb->theta_map_arg[mb->ntheta] = NULL;
			mb->ntheta++;
			ds->data_ntheta++;
		}
	}
		break;

	case L_ZEROINFLATEDBINOMIAL2:
	{
		/*
		 * get options related to the zeroinflatedbinomial2
		 */
		double initial_value = log(2.0);

		tmp = iniparser_getdouble(ini, inla_string_join(secname, "INITIAL"), initial_value);
		ds->data_fixed = iniparser_getboolean(ini, inla_string_join(secname, "FIXED"), 0);
		if (!ds->data_fixed && mb->mode_use_mode) {
			tmp = mb->theta_file[mb->theta_counter_file++];
			if (mb->mode_fixed)
				ds->data_fixed = 1;
		}
		HYPER_NEW(ds->data_observations.zeroinflated_alpha_intern, tmp);
		if (mb->verbose) {
			printf("\t\tinitialise alpha_intern[%g]\n", ds->data_observations.zeroinflated_alpha_intern[0][0]);
			printf("\t\tfixed=[%1d]\n", ds->data_fixed);
		}
		inla_read_prior(mb, ini, sec, &(ds->data_prior), "GAUSSIAN-std", NULL);

		/*
		 * add theta 
		 */
		if (!ds->data_fixed) {
			mb->theta = Realloc(mb->theta, mb->ntheta + 1, double **);
			mb->theta_hyperid = Realloc(mb->theta_hyperid, mb->ntheta + 1, char *);
			mb->theta_hyperid[mb->ntheta] = ds->data_prior.hyperid;
			mb->theta_tag = Realloc(mb->theta_tag, mb->ntheta + 1, char *);
			mb->theta_tag_userscale = Realloc(mb->theta_tag_userscale, mb->ntheta + 1, char *);
			mb->theta_dir = Realloc(mb->theta_dir, mb->ntheta + 1, char *);
			mb->theta_tag[mb->ntheta] = inla_make_tag("intern zero-probability parameter for zero-inflated binomial_2", mb->ds);
			mb->theta_tag_userscale[mb->ntheta] = inla_make_tag("zero-probability parameter for zero-inflated binomial_2", mb->ds);
			GMRFLib_sprintf(&msg, "%s-parameter", secname);
			mb->theta_dir[mb->ntheta] = msg;

			mb->theta_from = Realloc(mb->theta_from, mb->ntheta + 1, char *);
			mb->theta_to = Realloc(mb->theta_to, mb->ntheta + 1, char *);
			mb->theta_from[mb->ntheta] = Strdup(ds->data_prior.from_theta);
			mb->theta_to[mb->ntheta] = Strdup(ds->data_prior.to_theta);

			mb->theta[mb->ntheta] = ds->data_observations.zeroinflated_alpha_intern;
			mb->theta_map = Realloc(mb->theta_map, mb->ntheta + 1, map_func_tp *);
			mb->theta_map[mb->ntheta] = map_exp;
			mb->theta_map_arg = Realloc(mb->theta_map_arg, mb->ntheta + 1, void *);
			mb->theta_map_arg[mb->ntheta] = NULL;
			mb->ntheta++;
			ds->data_ntheta++;
		}
	}
		break;

	case L_ZERO_N_INFLATEDBINOMIAL2:
	{
		/*
		 * get options related to the zero_n_inflatedbinomial2
		 */
		double initial_value = log(2.0);

		tmp = iniparser_getdouble(ini, inla_string_join(secname, "INITIAL0"), initial_value);
		ds->data_fixed0 = iniparser_getboolean(ini, inla_string_join(secname, "FIXED0"), 0);
		if (!ds->data_fixed0 && mb->mode_use_mode) {
			tmp = mb->theta_file[mb->theta_counter_file++];
			if (mb->mode_fixed)
				ds->data_fixed0 = 1;
		}
		HYPER_NEW(ds->data_observations.zero_n_inflated_alpha1_intern, tmp);
		if (mb->verbose) {
			printf("\t\tinitialise alpha1_intern[%g]\n", ds->data_observations.zero_n_inflated_alpha1_intern[0][0]);
			printf("\t\tfixed=[%1d]\n", ds->data_fixed0);
		}
		inla_read_prior0(mb, ini, sec, &(ds->data_prior0), "GAUSSIAN-std", NULL);

		/*
		 * add theta 
		 */
		if (!ds->data_fixed0) {
			mb->theta = Realloc(mb->theta, mb->ntheta + 1, double **);
			mb->theta_hyperid = Realloc(mb->theta_hyperid, mb->ntheta + 1, char *);
			mb->theta_hyperid[mb->ntheta] = ds->data_prior0.hyperid;
			mb->theta_tag = Realloc(mb->theta_tag, mb->ntheta + 1, char *);
			mb->theta_tag_userscale = Realloc(mb->theta_tag_userscale, mb->ntheta + 1, char *);
			mb->theta_dir = Realloc(mb->theta_dir, mb->ntheta + 1, char *);
			mb->theta_tag[mb->ntheta] = inla_make_tag("intern alpha1 parameter for zero-n-inflated binomial_2", mb->ds);
			mb->theta_tag_userscale[mb->ntheta] = inla_make_tag("alpha1 parameter for zero-n-inflated binomial_2", mb->ds);
			GMRFLib_sprintf(&msg, "%s-parameter0", secname);
			mb->theta_dir[mb->ntheta] = msg;

			mb->theta_from = Realloc(mb->theta_from, mb->ntheta + 1, char *);
			mb->theta_to = Realloc(mb->theta_to, mb->ntheta + 1, char *);
			mb->theta_from[mb->ntheta] = Strdup(ds->data_prior0.from_theta);
			mb->theta_to[mb->ntheta] = Strdup(ds->data_prior0.to_theta);

			mb->theta[mb->ntheta] = ds->data_observations.zero_n_inflated_alpha1_intern;
			mb->theta_map = Realloc(mb->theta_map, mb->ntheta + 1, map_func_tp *);
			mb->theta_map[mb->ntheta] = map_exp;
			mb->theta_map_arg = Realloc(mb->theta_map_arg, mb->ntheta + 1, void *);
			mb->theta_map_arg[mb->ntheta] = NULL;
			mb->ntheta++;
			ds->data_ntheta++;
		}

		tmp = iniparser_getdouble(ini, inla_string_join(secname, "INITIAL1"), initial_value);
		ds->data_fixed1 = iniparser_getboolean(ini, inla_string_join(secname, "FIXED1"), 0);
		if (!ds->data_fixed1 && mb->mode_use_mode) {
			tmp = mb->theta_file[mb->theta_counter_file++];
			if (mb->mode_fixed)
				ds->data_fixed1 = 1;
		}
		HYPER_NEW(ds->data_observations.zero_n_inflated_alpha2_intern, tmp);
		if (mb->verbose) {
			printf("\t\tinitialise alpha2_intern[%g]\n", ds->data_observations.zero_n_inflated_alpha2_intern[0][0]);
			printf("\t\tfixed=[%1d]\n", ds->data_fixed1);
		}
		inla_read_prior1(mb, ini, sec, &(ds->data_prior1), "GAUSSIAN-std", NULL);

		/*
		 * add theta 
		 */
		if (!ds->data_fixed1) {
			mb->theta = Realloc(mb->theta, mb->ntheta + 1, double **);
			mb->theta_hyperid = Realloc(mb->theta_hyperid, mb->ntheta + 1, char *);
			mb->theta_hyperid[mb->ntheta] = ds->data_prior1.hyperid;
			mb->theta_tag = Realloc(mb->theta_tag, mb->ntheta + 1, char *);
			mb->theta_tag_userscale = Realloc(mb->theta_tag_userscale, mb->ntheta + 1, char *);
			mb->theta_dir = Realloc(mb->theta_dir, mb->ntheta + 1, char *);
			mb->theta_tag[mb->ntheta] = inla_make_tag("intern alpha2 parameter for zero-n-inflated binomial_2", mb->ds);
			mb->theta_tag_userscale[mb->ntheta] = inla_make_tag("alpha2 parameter for zero-n-inflated binomial_2", mb->ds);
			GMRFLib_sprintf(&msg, "%s-parameter1", secname);
			mb->theta_dir[mb->ntheta] = msg;
			mb->theta_from = Realloc(mb->theta_from, mb->ntheta + 1, char *);
			mb->theta_to = Realloc(mb->theta_to, mb->ntheta + 1, char *);
			mb->theta_from[mb->ntheta] = Strdup(ds->data_prior1.from_theta);
			mb->theta_to[mb->ntheta] = Strdup(ds->data_prior1.to_theta);
			mb->theta[mb->ntheta] = ds->data_observations.zero_n_inflated_alpha2_intern;
			mb->theta_map = Realloc(mb->theta_map, mb->ntheta + 1, map_func_tp *);
			mb->theta_map[mb->ntheta] = map_exp;
			mb->theta_map_arg = Realloc(mb->theta_map_arg, mb->ntheta + 1, void *);
			mb->theta_map_arg[mb->ntheta] = NULL;
			mb->ntheta++;
			ds->data_ntheta++;
		}
	}
		break;

	case L_ZERO_N_INFLATEDBINOMIAL3:
	{
		/*
		 * get options related to the zero_n_inflatedbinomial3
		 */
		double initial_value = log(2.0);

		tmp = iniparser_getdouble(ini, inla_string_join(secname, "INITIAL0"), initial_value);
		ds->data_fixed0 = iniparser_getboolean(ini, inla_string_join(secname, "FIXED0"), 0);
		if (!ds->data_fixed0 && mb->mode_use_mode) {
			tmp = mb->theta_file[mb->theta_counter_file++];
			if (mb->mode_fixed)
				ds->data_fixed0 = 1;
		}
		HYPER_NEW(ds->data_observations.zero_n_inflated_alpha0_intern, tmp);
		if (mb->verbose) {
			printf("\t\tinitialise alpha0_intern[%g]\n", ds->data_observations.zero_n_inflated_alpha0_intern[0][0]);
			printf("\t\tfixed=[%1d]\n", ds->data_fixed0);
		}
		inla_read_prior0(mb, ini, sec, &(ds->data_prior0), "GAUSSIAN-std", NULL);

		/*
		 * add theta 
		 */
		if (!ds->data_fixed0) {
			mb->theta = Realloc(mb->theta, mb->ntheta + 1, double **);
			mb->theta_hyperid = Realloc(mb->theta_hyperid, mb->ntheta + 1, char *);
			mb->theta_hyperid[mb->ntheta] = ds->data_prior0.hyperid;
			mb->theta_tag = Realloc(mb->theta_tag, mb->ntheta + 1, char *);
			mb->theta_tag_userscale = Realloc(mb->theta_tag_userscale, mb->ntheta + 1, char *);
			mb->theta_dir = Realloc(mb->theta_dir, mb->ntheta + 1, char *);
			mb->theta_tag[mb->ntheta] = inla_make_tag("intern alpha0 parameter for zero-n-inflated binomial_3", mb->ds);
			mb->theta_tag_userscale[mb->ntheta] = inla_make_tag("alpha0 parameter for zero-n-inflated binomial_3", mb->ds);
			GMRFLib_sprintf(&msg, "%s-parameter0", secname);
			mb->theta_dir[mb->ntheta] = msg;

			mb->theta_from = Realloc(mb->theta_from, mb->ntheta + 1, char *);
			mb->theta_to = Realloc(mb->theta_to, mb->ntheta + 1, char *);
			mb->theta_from[mb->ntheta] = Strdup(ds->data_prior0.from_theta);
			mb->theta_to[mb->ntheta] = Strdup(ds->data_prior0.to_theta);

			mb->theta[mb->ntheta] = ds->data_observations.zero_n_inflated_alpha0_intern;
			mb->theta_map = Realloc(mb->theta_map, mb->ntheta + 1, map_func_tp *);
			mb->theta_map[mb->ntheta] = map_exp;
			mb->theta_map_arg = Realloc(mb->theta_map_arg, mb->ntheta + 1, void *);
			mb->theta_map_arg[mb->ntheta] = NULL;
			mb->ntheta++;
			ds->data_ntheta++;
		}

		tmp = iniparser_getdouble(ini, inla_string_join(secname, "INITIAL1"), initial_value);
		ds->data_fixed1 = iniparser_getboolean(ini, inla_string_join(secname, "FIXED1"), 0);
		if (!ds->data_fixed1 && mb->mode_use_mode) {
			tmp = mb->theta_file[mb->theta_counter_file++];
			if (mb->mode_fixed)
				ds->data_fixed1 = 1;
		}
		HYPER_NEW(ds->data_observations.zero_n_inflated_alphaN_intern, tmp);
		if (mb->verbose) {
			printf("\t\tinitialise alphaN_intern[%g]\n", ds->data_observations.zero_n_inflated_alphaN_intern[0][0]);
			printf("\t\tfixed=[%1d]\n", ds->data_fixed1);
		}
		inla_read_prior1(mb, ini, sec, &(ds->data_prior1), "GAUSSIAN-std", NULL);

		/*
		 * add theta 
		 */
		if (!ds->data_fixed1) {
			mb->theta = Realloc(mb->theta, mb->ntheta + 1, double **);
			mb->theta_hyperid = Realloc(mb->theta_hyperid, mb->ntheta + 1, char *);
			mb->theta_hyperid[mb->ntheta] = ds->data_prior1.hyperid;
			mb->theta_tag = Realloc(mb->theta_tag, mb->ntheta + 1, char *);
			mb->theta_tag_userscale = Realloc(mb->theta_tag_userscale, mb->ntheta + 1, char *);
			mb->theta_dir = Realloc(mb->theta_dir, mb->ntheta + 1, char *);
			mb->theta_tag[mb->ntheta] = inla_make_tag("intern alphaN parameter for zero-n-inflated binomial_3", mb->ds);
			mb->theta_tag_userscale[mb->ntheta] = inla_make_tag("alphaN parameter for zero-n-inflated binomial_3", mb->ds);
			GMRFLib_sprintf(&msg, "%s-parameter1", secname);
			mb->theta_dir[mb->ntheta] = msg;
			mb->theta_from = Realloc(mb->theta_from, mb->ntheta + 1, char *);
			mb->theta_to = Realloc(mb->theta_to, mb->ntheta + 1, char *);
			mb->theta_from[mb->ntheta] = Strdup(ds->data_prior1.from_theta);
			mb->theta_to[mb->ntheta] = Strdup(ds->data_prior1.to_theta);
			mb->theta[mb->ntheta] = ds->data_observations.zero_n_inflated_alphaN_intern;
			mb->theta_map = Realloc(mb->theta_map, mb->ntheta + 1, map_func_tp *);
			mb->theta_map[mb->ntheta] = map_exp;
			mb->theta_map_arg = Realloc(mb->theta_map_arg, mb->ntheta + 1, void *);
			mb->theta_map_arg[mb->ntheta] = NULL;
			mb->ntheta++;
			ds->data_ntheta++;
		}
	}
		break;

	case L_ZEROINFLATEDBETABINOMIAL2:
	{
		/*
		 * get options related to the zeroinflatedbetabinomial2
		 */
		double initial_value = log(2.0);

		tmp = iniparser_getdouble(ini, inla_string_join(secname, "INITIAL0"), initial_value);
		ds->data_fixed0 = iniparser_getboolean(ini, inla_string_join(secname, "FIXED0"), 0);
		if (!ds->data_fixed0 && mb->mode_use_mode) {
			tmp = mb->theta_file[mb->theta_counter_file++];
			if (mb->mode_fixed)
				ds->data_fixed0 = 1;
		}
		HYPER_NEW(ds->data_observations.zeroinflated_alpha_intern, tmp);
		if (mb->verbose) {
			printf("\t\tinitialise alpha_intern[%g]\n", ds->data_observations.zeroinflated_alpha_intern[0][0]);
			printf("\t\tfixed=[%1d]\n", ds->data_fixed0);
		}
		inla_read_prior0(mb, ini, sec, &(ds->data_prior0), "GAUSSIAN-std", NULL);

		/*
		 * add theta 
		 */
		if (!ds->data_fixed0) {
			mb->theta = Realloc(mb->theta, mb->ntheta + 1, double **);
			mb->theta_hyperid = Realloc(mb->theta_hyperid, mb->ntheta + 1, char *);
			mb->theta_hyperid[mb->ntheta] = ds->data_prior0.hyperid;
			mb->theta_tag = Realloc(mb->theta_tag, mb->ntheta + 1, char *);
			mb->theta_tag_userscale = Realloc(mb->theta_tag_userscale, mb->ntheta + 1, char *);
			mb->theta_dir = Realloc(mb->theta_dir, mb->ntheta + 1, char *);
			mb->theta_tag[mb->ntheta] = inla_make_tag("intern zero-probability parameter for zero-inflated betabinomial_2", mb->ds);
			mb->theta_tag_userscale[mb->ntheta] = inla_make_tag("zero-probability parameter for zero-inflated betabinomial_2", mb->ds);
			GMRFLib_sprintf(&msg, "%s-parameter0", secname);
			mb->theta_dir[mb->ntheta] = msg;

			mb->theta_from = Realloc(mb->theta_from, mb->ntheta + 1, char *);
			mb->theta_to = Realloc(mb->theta_to, mb->ntheta + 1, char *);
			mb->theta_from[mb->ntheta] = Strdup(ds->data_prior0.from_theta);
			mb->theta_to[mb->ntheta] = Strdup(ds->data_prior0.to_theta);

			mb->theta[mb->ntheta] = ds->data_observations.zeroinflated_alpha_intern;
			mb->theta_map = Realloc(mb->theta_map, mb->ntheta + 1, map_func_tp *);
			mb->theta_map[mb->ntheta] = map_exp;
			mb->theta_map_arg = Realloc(mb->theta_map_arg, mb->ntheta + 1, void *);
			mb->theta_map_arg[mb->ntheta] = NULL;
			mb->ntheta++;
			ds->data_ntheta++;
		}

		initial_value = log(1.0);
		tmp = iniparser_getdouble(ini, inla_string_join(secname, "INITIAL1"), initial_value);
		ds->data_fixed1 = iniparser_getboolean(ini, inla_string_join(secname, "FIXED1"), 0);
		if (!ds->data_fixed1 && mb->mode_use_mode) {
			tmp = mb->theta_file[mb->theta_counter_file++];
			if (mb->mode_fixed)
				ds->data_fixed1 = 1;
		}
		HYPER_NEW(ds->data_observations.zeroinflated_delta_intern, tmp);
		if (mb->verbose) {
			printf("\t\tinitialise delta_intern[%g]\n", ds->data_observations.zeroinflated_delta_intern[0][0]);
			printf("\t\tfixed=[%1d]\n", ds->data_fixed1);
		}
		inla_read_prior1(mb, ini, sec, &(ds->data_prior1), "GAUSSIAN-std", NULL);

		/*
		 * add theta 
		 */
		if (!ds->data_fixed1) {
			mb->theta = Realloc(mb->theta, mb->ntheta + 1, double **);
			mb->theta_hyperid = Realloc(mb->theta_hyperid, mb->ntheta + 1, char *);
			mb->theta_hyperid[mb->ntheta] = ds->data_prior1.hyperid;
			mb->theta_tag = Realloc(mb->theta_tag, mb->ntheta + 1, char *);
			mb->theta_tag_userscale = Realloc(mb->theta_tag_userscale, mb->ntheta + 1, char *);
			mb->theta_dir = Realloc(mb->theta_dir, mb->ntheta + 1, char *);
			mb->theta_tag[mb->ntheta] = inla_make_tag("intern overdispersion parameter for zero-inflated betabinomial_2", mb->ds);
			mb->theta_tag_userscale[mb->ntheta] = inla_make_tag("overdispersion parameter for zero-inflated betabinomial_2", mb->ds);
			GMRFLib_sprintf(&msg, "%s-parameter1", secname);
			mb->theta_dir[mb->ntheta] = msg;

			mb->theta_from = Realloc(mb->theta_from, mb->ntheta + 1, char *);
			mb->theta_to = Realloc(mb->theta_to, mb->ntheta + 1, char *);
			mb->theta_from[mb->ntheta] = Strdup(ds->data_prior1.from_theta);
			mb->theta_to[mb->ntheta] = Strdup(ds->data_prior1.to_theta);

			mb->theta[mb->ntheta] = ds->data_observations.zeroinflated_delta_intern;
			mb->theta_map = Realloc(mb->theta_map, mb->ntheta + 1, map_func_tp *);
			mb->theta_map[mb->ntheta] = map_exp;
			mb->theta_map_arg = Realloc(mb->theta_map_arg, mb->ntheta + 1, void *);
			mb->theta_map_arg[mb->ntheta] = NULL;
			mb->ntheta++;
			ds->data_ntheta++;
		}
	}
		break;

	case L_ZEROINFLATEDBINOMIAL0:
	case L_ZEROINFLATEDBINOMIAL1:
	{
		/*
		 * get options related to the zeroinflatedbinomial0/1
		 */
		double initial_value = -1.0;

		tmp = iniparser_getdouble(ini, inla_string_join(secname, "INITIAL"), initial_value);
		ds->data_fixed = iniparser_getboolean(ini, inla_string_join(secname, "FIXED"), 0);
		if (!ds->data_fixed && mb->mode_use_mode) {
			tmp = mb->theta_file[mb->theta_counter_file++];
			if (mb->mode_fixed)
				ds->data_fixed = 1;
		}
		HYPER_NEW(ds->data_observations.prob_intern, tmp);
		if (mb->verbose) {
			printf("\t\tinitialise prob_intern[%g]\n", ds->data_observations.prob_intern[0][0]);
			printf("\t\tfixed=[%1d]\n", ds->data_fixed);
		}
		inla_read_prior(mb, ini, sec, &(ds->data_prior), "GAUSSIAN-std", NULL);

		/*
		 * add theta 
		 */
		if (!ds->data_fixed) {
			mb->theta = Realloc(mb->theta, mb->ntheta + 1, double **);
			mb->theta_hyperid = Realloc(mb->theta_hyperid, mb->ntheta + 1, char *);
			mb->theta_hyperid[mb->ntheta] = ds->data_prior.hyperid;
			mb->theta_tag = Realloc(mb->theta_tag, mb->ntheta + 1, char *);
			mb->theta_tag_userscale = Realloc(mb->theta_tag_userscale, mb->ntheta + 1, char *);
			mb->theta_dir = Realloc(mb->theta_dir, mb->ntheta + 1, char *);
			if (ds->data_id == L_ZEROINFLATEDBINOMIAL0) {
				mb->theta_tag[mb->ntheta] = inla_make_tag("intern zero-probability parameter for zero-inflated binomial_0", mb->ds);
				mb->theta_tag_userscale[mb->ntheta] =
				    inla_make_tag("zero-probability parameter for zero-inflated binomial_0", mb->ds);
			} else {
				mb->theta_tag[mb->ntheta] = inla_make_tag("intern zero-probability parameter for zero-inflated binomial_1", mb->ds);
				mb->theta_tag_userscale[mb->ntheta] =
				    inla_make_tag("zero-probability parameter for zero-inflated binomial_1", mb->ds);
			}
			GMRFLib_sprintf(&msg, "%s-parameter", secname);
			mb->theta_dir[mb->ntheta] = msg;

			mb->theta_from = Realloc(mb->theta_from, mb->ntheta + 1, char *);
			mb->theta_to = Realloc(mb->theta_to, mb->ntheta + 1, char *);
			mb->theta_from[mb->ntheta] = Strdup(ds->data_prior.from_theta);
			mb->theta_to[mb->ntheta] = Strdup(ds->data_prior.to_theta);

			mb->theta[mb->ntheta] = ds->data_observations.prob_intern;
			mb->theta_map = Realloc(mb->theta_map, mb->ntheta + 1, map_func_tp *);
			mb->theta_map[mb->ntheta] = map_probability;
			mb->theta_map_arg = Realloc(mb->theta_map_arg, mb->ntheta + 1, void *);
			mb->theta_map_arg[mb->ntheta] = NULL;
			mb->ntheta++;
			ds->data_ntheta++;
		}
	}
		break;

	case L_NMIX:
	case L_NMIXNB:
	{
		/*
		 * get options related to the nmix and nmixnb
		 */
		char *suff = (ds->data_id == L_NMIX ? Strdup("") : Strdup("nb"));
		if (mb->verbose) {
			printf("\t\tmodel for N in the mixture[%s]\n", (ds->data_id == L_NMIX ? "Poisson" : "NegativeBinomial"));
		}

		// first we need to know 'm'. 
		found = 0;
		ds->data_observations.nmix_m = NMIX_MMAX;
		for (i = 0; i < NMIX_MMAX && !found; i++) {
			for (j = 0; j < mb->predictor_ndata; j++) {
				if (gsl_isnan(ds->data_observations.nmix_x[i][j])) {
					found = 1;
					ds->data_observations.nmix_m = i;
					break;
				}
			}
		}
		if (mb->verbose) {
			printf("\t\tnmix%s.m=[%1d]\n", suff, ds->data_observations.nmix_m);
		}
		assert(ds->data_observations.nmix_m > 0 && ds->data_observations.nmix_m <= NMIX_MMAX);
		ds->data_observations.nmix_beta = Calloc(NMIX_MMAX + 1, double **);	/* yes, its +1 to cover the NB case */
		ds->data_nprior = Calloc(NMIX_MMAX + 1, Prior_tp);
		ds->data_nfixed = Calloc(NMIX_MMAX + 1, int);

		for (int k = 0; k < NMIX_MMAX; k++) {
			ds->data_nfixed[k] = 1;		       /* so that the unused ones are fixed, so we can loop over all in the 'extra'
							        * function */
		}

		// mark all as read
		for (i = 0; i < NMIX_MMAX; i++) {
			for (j = 0; j < keywords_len; j++) {
				GMRFLib_sprintf(&ctmp, "%s%1d", keywords[j], i);
				iniparser_getstring(ini, inla_string_join(secname, ctmp), NULL);
				Free(ctmp);
			}
		}

		for (int k = 0; k < ds->data_observations.nmix_m; k++) {
			GMRFLib_sprintf(&ctmp, "INITIAL%1d", k);
			tmp = iniparser_getdouble(ini, inla_string_join(secname, ctmp), 0.0);

			Free(ctmp);
			GMRFLib_sprintf(&ctmp, "FIXED%1d", k);
			ds->data_nfixed[k] = iniparser_getboolean(ini, inla_string_join(secname, ctmp), 0);

			if (!ds->data_nfixed[k] && mb->mode_use_mode) {
				tmp = mb->theta_file[mb->theta_counter_file++];
				if (mb->mode_fixed)
					ds->data_nfixed[k] = 1;
			}
			HYPER_NEW(ds->data_observations.nmix_beta[k], tmp);
			if (mb->verbose) {
				printf("\t\tinitialise nmix%s.beta[%1d] = %g\n", suff, k, ds->data_observations.nmix_beta[k][0][0]);
				printf("\t\tfixed = %1d\n", ds->data_nfixed[k]);
			}
			inla_read_priorN(mb, ini, sec, &(ds->data_nprior[k]), "GAUSSIAN", k, NULL);

			if (!ds->data_nfixed[k]) {
				mb->theta = Realloc(mb->theta, mb->ntheta + 1, double **);
				mb->theta_hyperid = Realloc(mb->theta_hyperid, mb->ntheta + 1, char *);
				mb->theta_hyperid[mb->ntheta] = ds->data_nprior[k].hyperid;
				mb->theta_tag = Realloc(mb->theta_tag, mb->ntheta + 1, char *);
				mb->theta_tag_userscale = Realloc(mb->theta_tag_userscale, mb->ntheta + 1, char *);
				mb->theta_dir = Realloc(mb->theta_dir, mb->ntheta + 1, char *);

				Free(ctmp);
				if (ds->data_id == L_NMIX) {
					GMRFLib_sprintf(&ctmp, "beta[%1d] for NMix observations", k + 1);
				} else {
					GMRFLib_sprintf(&ctmp, "beta[%1d] for NMixNB observations", k + 1);
				}
				mb->theta_tag[mb->ntheta] = inla_make_tag(ctmp, mb->ds);
				mb->theta_tag_userscale[mb->ntheta] = inla_make_tag(ctmp, mb->ds);
				GMRFLib_sprintf(&msg, "%s-parameter", secname);
				mb->theta_dir[mb->ntheta] = msg;

				mb->theta_from = Realloc(mb->theta_from, mb->ntheta + 1, char *);
				mb->theta_to = Realloc(mb->theta_to, mb->ntheta + 1, char *);
				mb->theta_from[mb->ntheta] = Strdup(ds->data_nprior[k].from_theta);
				mb->theta_to[mb->ntheta] = Strdup(ds->data_nprior[k].to_theta);

				mb->theta[mb->ntheta] = ds->data_observations.nmix_beta[k];
				mb->theta_map = Realloc(mb->theta_map, mb->ntheta + 1, map_func_tp *);
				mb->theta_map[mb->ntheta] = map_identity;
				mb->theta_map_arg = Realloc(mb->theta_map_arg, mb->ntheta + 1, void *);
				mb->theta_map_arg[mb->ntheta] = NULL;
				mb->ntheta++;
				ds->data_ntheta++;
			}
		}

		if (ds->data_id == L_NMIXNB) {
			int k = NMIX_MMAX;		       /* this the overdisperson */

			GMRFLib_sprintf(&ctmp, "INITIAL%1d", k);
			tmp = iniparser_getdouble(ini, inla_string_join(secname, ctmp), 0.0);
			HYPER_NEW(ds->data_observations.nmix_log_overdispersion, tmp);

			Free(ctmp);
			GMRFLib_sprintf(&ctmp, "FIXED%1d", k);
			ds->data_nfixed[k] = iniparser_getboolean(ini, inla_string_join(secname, ctmp), 0);

			if (!ds->data_nfixed[k] && mb->mode_use_mode) {
				tmp = mb->theta_file[mb->theta_counter_file++];
				if (mb->mode_fixed)
					ds->data_nfixed[k] = 1;
			}

			if (mb->verbose) {
				printf("\t\tinitialise nmix%s.log_overdispersion = %g\n", suff,
				       ds->data_observations.nmix_log_overdispersion[0][0]);
				printf("\t\tfixed = %1d\n", ds->data_nfixed[k]);
			}
			inla_read_priorN(mb, ini, sec, &(ds->data_nprior[k]), "LOGGAMMA", k, NULL);

			if (!ds->data_nfixed[k]) {
				mb->theta = Realloc(mb->theta, mb->ntheta + 1, double **);
				mb->theta_hyperid = Realloc(mb->theta_hyperid, mb->ntheta + 1, char *);
				mb->theta_hyperid[mb->ntheta] = ds->data_nprior[k].hyperid;
				mb->theta_tag = Realloc(mb->theta_tag, mb->ntheta + 1, char *);
				mb->theta_tag_userscale = Realloc(mb->theta_tag_userscale, mb->ntheta + 1, char *);
				mb->theta_dir = Realloc(mb->theta_dir, mb->ntheta + 1, char *);

				Free(ctmp);
				GMRFLib_sprintf(&ctmp, "log_overdispersion for NMixNB observations");
				mb->theta_tag[mb->ntheta] = inla_make_tag(ctmp, mb->ds);
				Free(ctmp);
				GMRFLib_sprintf(&ctmp, "overdispersion for NMixNB observations");
				mb->theta_tag_userscale[mb->ntheta] = inla_make_tag(ctmp, mb->ds);
				GMRFLib_sprintf(&msg, "%s-parameter", secname);
				mb->theta_dir[mb->ntheta] = msg;

				mb->theta_from = Realloc(mb->theta_from, mb->ntheta + 1, char *);
				mb->theta_to = Realloc(mb->theta_to, mb->ntheta + 1, char *);
				mb->theta_from[mb->ntheta] = Strdup(ds->data_nprior[k].from_theta);
				mb->theta_to[mb->ntheta] = Strdup(ds->data_nprior[k].to_theta);

				mb->theta[mb->ntheta] = ds->data_observations.nmix_log_overdispersion;
				mb->theta_map = Realloc(mb->theta_map, mb->ntheta + 1, map_func_tp *);
				mb->theta_map[mb->ntheta] = map_exp;
				mb->theta_map_arg = Realloc(mb->theta_map_arg, mb->ntheta + 1, void *);
				mb->theta_map_arg[mb->ntheta] = NULL;
				mb->ntheta++;
				ds->data_ntheta++;
			}
		}
	}
		break;

	case L_OCCUPANCY:
	{
		if (mb->verbose) {
			printf("\t\tny_max=[%1d]\n", ds->data_observations.occ_ny_max);
			printf("\t\tnbeta=[%1d]\n", ds->data_observations.occ_nbeta);
		}
		int nb = ds->data_observations.occ_nbeta;
		ds->data_observations.occ_beta = Calloc(nb, double **);
		ds->data_nprior = Calloc(OCCUPANCY_MAXTHETA, Prior_tp);
		ds->data_nfixed = Calloc(OCCUPANCY_MAXTHETA, int);

		for (int k = 0; k < OCCUPANCY_MAXTHETA; k++) {
			ds->data_nfixed[k] = 1;
		}

		// mark all as read
		for (i = 0; i < OCCUPANCY_MAXTHETA; i++) {
			for (j = 0; j < keywords_len; j++) {
				GMRFLib_sprintf(&ctmp, "%s%1d", keywords[j], i);
				iniparser_getstring(ini, inla_string_join(secname, ctmp), NULL);
				Free(ctmp);
			}
		}

		for (int k = 0; k < nb; k++) {
			GMRFLib_sprintf(&ctmp, "INITIAL%1d", k);
			tmp = iniparser_getdouble(ini, inla_string_join(secname, ctmp), 0.0);
			Free(ctmp);

			GMRFLib_sprintf(&ctmp, "FIXED%1d", k);
			ds->data_nfixed[k] = iniparser_getboolean(ini, inla_string_join(secname, ctmp), 0);

			if (!ds->data_nfixed[k] && mb->mode_use_mode) {
				tmp = mb->theta_file[mb->theta_counter_file++];
				if (mb->mode_fixed)
					ds->data_nfixed[k] = 1;
			}
			HYPER_NEW(ds->data_observations.occ_beta[k], tmp);
			if (mb->verbose) {
				printf("\t\tinitialise occ_beta[%1d] = %g\n", k, ds->data_observations.occ_beta[k][0][0]);
				printf("\t\tfixed = %1d\n", ds->data_nfixed[k]);
			}
			inla_read_priorN(mb, ini, sec, &(ds->data_nprior[k]), "GAUSSIAN", k, NULL);

			if (!ds->data_nfixed[k]) {
				mb->theta = Realloc(mb->theta, mb->ntheta + 1, double **);
				mb->theta_hyperid = Realloc(mb->theta_hyperid, mb->ntheta + 1, char *);
				mb->theta_hyperid[mb->ntheta] = ds->data_nprior[k].hyperid;
				mb->theta_tag = Realloc(mb->theta_tag, mb->ntheta + 1, char *);
				mb->theta_tag_userscale = Realloc(mb->theta_tag_userscale, mb->ntheta + 1, char *);
				mb->theta_dir = Realloc(mb->theta_dir, mb->ntheta + 1, char *);

				Free(ctmp);
				GMRFLib_sprintf(&ctmp, "beta[%1d] for occupancy observations", k);
				mb->theta_tag[mb->ntheta] = inla_make_tag(ctmp, mb->ds);
				mb->theta_tag_userscale[mb->ntheta] = inla_make_tag(ctmp, mb->ds);
				GMRFLib_sprintf(&msg, "%s-parameter", secname);
				mb->theta_dir[mb->ntheta] = msg;

				mb->theta_from = Realloc(mb->theta_from, mb->ntheta + 1, char *);
				mb->theta_to = Realloc(mb->theta_to, mb->ntheta + 1, char *);
				mb->theta_from[mb->ntheta] = Strdup(ds->data_nprior[k].from_theta);
				mb->theta_to[mb->ntheta] = Strdup(ds->data_nprior[k].to_theta);

				mb->theta[mb->ntheta] = ds->data_observations.occ_beta[k];
				mb->theta_map = Realloc(mb->theta_map, mb->ntheta + 1, map_func_tp *);
				mb->theta_map[mb->ntheta] = map_identity;
				mb->theta_map_arg = Realloc(mb->theta_map_arg, mb->ntheta + 1, void *);
				mb->theta_map_arg[mb->ntheta] = NULL;
				mb->ntheta++;
				ds->data_ntheta++;
			}
		}

		char *link_simple = iniparser_getstring(ini, inla_string_join(secname, "LINK.SIMPLE"), NULL);
		ds->data_observations.link_simple_name = link_simple;
		if (!strcasecmp(link_simple, "IDENTITY")) {
			ds->data_observations.link_simple_invlinkfunc = link_identity;
		} else if (!strcasecmp(link_simple, "LOG")) {
			ds->data_observations.link_simple_invlinkfunc = link_log;
		} else if (!strcasecmp(link_simple, "PROBIT")) {
			ds->data_observations.link_simple_invlinkfunc = link_probit;
		} else if (!strcasecmp(link_simple, "CLOGLOG")) {
			ds->data_observations.link_simple_invlinkfunc = link_cloglog;
		} else if (!strcasecmp(link_simple, "LOGIT")) {
			ds->data_observations.link_simple_invlinkfunc = link_logit;
		} else {
			GMRFLib_sprintf(&msg, "%s: 0poisson(S) likelihood: no valid link.simple[%s]", secname, link_simple);
			inla_error_general(msg);
			exit(1);
		}
		if (mb->verbose) {
			printf("\t\tlink.simple[%s]\n", ds->data_observations.link_simple_name);
		}
	}
		break;

	default:
		/*
		 * nothing to do 
		 */
		ds->data_ntheta = 0;
	}

	/*
	 * setup the link-model
	 */
	ds->link_model = Strdup(strupc(iniparser_getstring(ini, inla_string_join(secname, "LINK.MODEL"), Strdup("default"))));
	inla_trim_family(ds->link_model);

	if (!strcasecmp(ds->link_model, "IDENTITY")) {
		ds->link_id = LINK_IDENTITY;
		ds->link_ntheta = 0;
		ds->predictor_invlinkfunc = link_identity;
	} else if (!strcasecmp(ds->link_model, "inverse")) {
		ds->link_id = LINK_INVERSE;
		ds->link_ntheta = 0;
		ds->predictor_invlinkfunc = link_inverse;
	} else if (!strcasecmp(ds->link_model, "LOG")) {
		ds->link_id = LINK_LOG;
		ds->link_ntheta = 0;
		ds->predictor_invlinkfunc = link_log;
	} else if (!strcasecmp(ds->link_model, "LOGa")) {
		ds->link_id = LINK_LOGa;
		ds->link_ntheta = 0;
		ds->predictor_invlinkfunc = link_loga;
	} else if (!strcasecmp(ds->link_model, "NEGLOG")) {
		ds->link_id = LINK_NEGLOG;
		ds->link_ntheta = 0;
		ds->predictor_invlinkfunc = link_neglog;
	} else if (!strcasecmp(ds->link_model, "PROBIT")) {
		ds->link_id = LINK_PROBIT;
		ds->link_ntheta = 0;
		ds->predictor_invlinkfunc = link_probit;
	} else if (!strcasecmp(ds->link_model, "CLOGLOG")) {
		ds->link_id = LINK_CLOGLOG;
		ds->link_ntheta = 0;
		ds->predictor_invlinkfunc = link_cloglog;
	} else if (!strcasecmp(ds->link_model, "CCLOGLOG")) {
		ds->link_id = LINK_CCLOGLOG;
		ds->link_ntheta = 0;
		ds->predictor_invlinkfunc = link_ccloglog;
	} else if (!strcasecmp(ds->link_model, "LOGLOG")) {
		ds->link_id = LINK_LOGLOG;
		ds->link_ntheta = 0;
		ds->predictor_invlinkfunc = link_loglog;
	} else if (!strcasecmp(ds->link_model, "CAUCHIT")) {
		ds->link_id = LINK_CAUCHIT;
		ds->link_ntheta = 0;
		ds->predictor_invlinkfunc = link_cauchit;
	} else if (!strcasecmp(ds->link_model, "LOGIT")) {
		ds->link_id = LINK_LOGIT;
		ds->link_ntheta = 0;
		ds->predictor_invlinkfunc = link_logit;
	} else if (!strcasecmp(ds->link_model, "TAN")) {
		ds->link_id = LINK_TAN;
		ds->link_ntheta = 0;
		ds->predictor_invlinkfunc = link_tan;
	} else if (!strcasecmp(ds->link_model, "LOGOFFSET")) {
		ds->link_id = LINK_LOGOFFSET;
		ds->link_ntheta = 1;
		ds->predictor_invlinkfunc = link_logoffset;
	} else if (!strcasecmp(ds->link_model, "LOGITOFFSET")) {
		ds->link_id = LINK_LOGITOFFSET;
		ds->link_ntheta = 1;
		ds->predictor_invlinkfunc = link_logitoffset;
	} else if (!strcasecmp(ds->link_model, "SSLOGIT")) {
		ds->link_id = LINK_SSLOGIT;
		ds->link_ntheta = 2;
		ds->predictor_invlinkfunc = link_sslogit;
	} else if (!strcasecmp(ds->link_model, "ROBIT")) {
		ds->link_id = LINK_ROBIT;
		ds->link_ntheta = 1;
		ds->predictor_invlinkfunc = link_robit;
	} else if (!strcasecmp(ds->link_model, "SN")) {
		ds->link_id = LINK_SN;
		ds->link_ntheta = 2;
		ds->predictor_invlinkfunc = link_sn;
	} else if (!strcasecmp(ds->link_model, "GEV")) {
		ds->link_id = LINK_GEV;
		ds->link_ntheta = 2;
		ds->predictor_invlinkfunc = link_gev;
	} else if (!strcasecmp(ds->link_model, "CGEV")) {
		ds->link_id = LINK_CGEV;
		ds->link_ntheta = 2;
		ds->predictor_invlinkfunc = link_cgev;
	} else if (!strcasecmp(ds->link_model, "POWERLOGIT")) {
		ds->link_id = LINK_POWER_LOGIT;
		ds->link_ntheta = 2;
		ds->predictor_invlinkfunc = link_power_logit;
	} else if (!strcasecmp(ds->link_model, "TEST1")) {
		ds->link_id = LINK_TEST1;
		ds->link_ntheta = 1;
		ds->predictor_invlinkfunc = link_test1;
	} else if (!strcasecmp(ds->link_model, "SPECIAL1")) {
		ds->link_id = LINK_SPECIAL1;
		ds->link_ntheta = -1;
		ds->predictor_invlinkfunc = link_special1;
	} else if (!strcasecmp(ds->link_model, "SPECIAL2")) {
		ds->link_id = LINK_SPECIAL2;
		ds->link_ntheta = 1;
		ds->predictor_invlinkfunc = link_special2;
	} else if (!strcasecmp(ds->link_model, "QUANTILE")) {
		GMRFLib_ASSERT((ds->data_observations.quantile > 0.0 && ds->data_observations.quantile < 1.0), GMRFLib_EPARAMETER);
		switch (ds->data_id) {
		case L_POISSON:
		case L_XPOISSON:
		{
			ds->link_id = LINK_QPOISSON;
			ds->link_ntheta = 0;
			ds->predictor_invlinkfunc = link_qpoisson;
		}
			break;

		case L_BINOMIAL:
		case L_XBINOMIAL:
		{
			ds->link_id = LINK_QBINOMIAL;
			ds->link_ntheta = 0;
			ds->predictor_invlinkfunc = link_qbinomial;
		}
			break;
		case L_WEIBULL:
		case L_WEIBULLSURV:
		{
			ds->link_id = LINK_QWEIBULL;
			ds->link_ntheta = 0;
			ds->predictor_invlinkfunc = link_qweibull;
		}
			break;

		case L_GAMMA:
		case L_GAMMASURV:
		{
			ds->link_id = LINK_QGAMMA;
			ds->link_ntheta = 0;
			ds->predictor_invlinkfunc = link_qgamma;
			inla_qgamma_cache(0.0, ds->data_observations.quantile);
		}
			break;

		case L_GP:
		{
			ds->link_id = LINK_LOG;
			ds->link_ntheta = 0;
			ds->predictor_invlinkfunc = link_log;
		}
			break;

		case L_EGP:
		{
			ds->link_id = LINK_LOG;
			ds->link_ntheta = 0;
			ds->predictor_invlinkfunc = link_log;
		}
			break;

		case L_DGP:
		{
			ds->link_id = LINK_LOG;
			ds->link_ntheta = 0;
			ds->predictor_invlinkfunc = link_log;
		}
			break;

		case L_EXPPOWER:
		{
			ds->link_id = LINK_QEXPPOWER;
			ds->link_ntheta = 0;
			ds->predictor_invlinkfunc = link_qexppower;
		}
			break;

		default:
			assert(0 == 1);
		}
	} else if (!strcasecmp(ds->link_model, "PQUANTILE")) {
		GMRFLib_ASSERT((ds->data_observations.quantile > 0.0 && ds->data_observations.quantile < 1.0), GMRFLib_EPARAMETER);
		switch (ds->data_id) {
		case L_BINOMIAL:
		case L_XBINOMIAL:
		{
			ds->link_id = LINK_QBINOMIAL;
			ds->link_ntheta = 0;
			ds->predictor_invlinkfunc = link_pqbinomial;
		}
			break;
		default:
			assert(0 == 1);
		}
	} else {
		GMRFLib_sprintf(&msg, "Unknown link-model [%s]\n", ds->link_model);
		inla_error_general(msg);
		exit(1);
	}

	/*
	 * the meaning of 'order' is model dependent
	 */
	ds->link_order = iniparser_getint(ini, inla_string_join(secname, "LINK.ORDER"), -1);
	ds->link_variant = iniparser_getint(ini, inla_string_join(secname, "LINK.VARIANT"), -1);
	ds->link_a = iniparser_getdouble(ini, inla_string_join(secname, "LINK.a"), 1.0);
	if (mb->verbose) {
		printf("\t\tLink model   [%s]\n", ds->link_model);
		printf("\t\tLink order   [%1d]\n", ds->link_order);
		printf("\t\tLink variant [%1d]\n", ds->link_variant);
		printf("\t\tLink a       [%g]\n", ds->link_a);
		printf("\t\tLink ntheta  [%1d]\n", ds->link_ntheta);
	}

	/*
	 * read possible link_covariates 
	 */
	char *link_cov_filename = NULL;
	link_cov_filename = Strdup(iniparser_getstring(ini, inla_string_join(secname, "LINK.COVARIATES"), NULL));

	if (link_cov_filename) {
		ds->link_covariates = GMRFLib_read_fmesher_file(link_cov_filename, (long int) 0, -1);
		if (mb->verbose) {
			int ii, jj;
			printf("\t\tLink_covariates: file[%s] dim=(%1d x %1d)\n", link_cov_filename, ds->link_covariates->nrow,
			       ds->link_covariates->ncol);
			for (ii = 0; ii < IMIN(PREVIEW, ds->link_covariates->nrow); ii++) {
				printf("\t\t\trow=%1d, covariates=[", ii);
				for (jj = 0; jj < ds->link_covariates->ncol; jj++) {
					printf(" %g", GMRFLib_matrix_get(ii, jj, ds->link_covariates));
				}
				printf(" ]\n");
			}
		}

		int len = (mb->predictor_m > 0 ? mb->predictor_m : mb->predictor_n);
		if (len != ds->link_covariates->nrow) {
			char *emsg = NULL;
			GMRFLib_sprintf(&emsg,
					"link.covariates has not the same number of rows as the linear predictor %1d != %1d",
					ds->link_covariates->nrow, len);
			inla_error_general(emsg);
			abort();
		}
	} else {
		ds->link_covariates = NULL;
	}

	/*
	 * read link-parameters, if any
	 */
	ds->predictor_invlinkfunc_arg = Calloc(n_data, void *);

	switch (ds->link_id) {
	case LINK_IDENTITY:
	case LINK_INVERSE:
	case LINK_LOG:
	case LINK_NEGLOG:
	case LINK_PROBIT:
	case LINK_CLOGLOG:
	case LINK_CCLOGLOG:
	case LINK_LOGLOG:
	case LINK_CAUCHIT:
	case LINK_LOGIT:
	case LINK_TAN:
		break;

	case LINK_LOGa:
	{
		Link_param_tp *link_param = Calloc(1, Link_param_tp);
		link_param->a = ds->link_a;
		for (i = 0; i < n_data; i++) {
			ds->predictor_invlinkfunc_arg[i] = (void *) link_param;
		}
	}
		break;

	case LINK_QPOISSON:
	{
		Link_param_tp *link_param = Calloc(1, Link_param_tp);
		link_param->idx = -1;
		link_param->quantile = ds->data_observations.quantile;
		for (i = 0; i < n_data; i++) {
			ds->predictor_invlinkfunc_arg[i] = (void *) link_param;
		}
	}
		break;

	case LINK_QBINOMIAL:
	{
		for (i = 0; i < n_data; i++) {
			Link_param_tp *link_param = Calloc(1, Link_param_tp);
			link_param->idx = i;
			link_param->quantile = ds->data_observations.quantile;
			link_param->Ntrial = ds->data_observations.nb[i];
			ds->predictor_invlinkfunc_arg[i] = (void *) link_param;
		}
	}
		break;

	case LINK_QWEIBULL:
	{
		Link_param_tp *link_param = Calloc(1, Link_param_tp);
		link_param->idx = -1;
		link_param->quantile = ds->data_observations.quantile;
		link_param->alpha_intern = ds->data_observations.alpha_intern;
		link_param->variant = ds->variant;
		for (i = 0; i < n_data; i++) {
			ds->predictor_invlinkfunc_arg[i] = (void *) link_param;
		}
	}
		break;

	case LINK_QGAMMA:
	{
		for (i = 0; i < n_data; i++) {
			Link_param_tp *link_param = Calloc(1, Link_param_tp);
			link_param->idx = i;
			link_param->quantile = ds->data_observations.quantile;
			link_param->scale = ds->data_observations.gamma_scale;
			link_param->log_prec = ds->data_observations.gamma_log_prec;
			ds->predictor_invlinkfunc_arg[i] = (void *) link_param;
		}
	}
		break;

	case LINK_QEXPPOWER:
	{
		for (i = 0; i < n_data; i++) {
			Link_param_tp *link_param = Calloc(1, Link_param_tp);
			link_param->idx = i;
			link_param->quantile = ds->data_observations.quantile;
			link_param->scale = ds->data_observations.weight_gaussian;
			link_param->log_prec = ds->data_observations.log_prec_gaussian;
			link_param->log_power = ds->data_observations.log_power;
			ds->predictor_invlinkfunc_arg[i] = (void *) link_param;
		}
	}
		break;

	case LINK_SSLOGIT:
	{
		/*
		 * logit link with correction for sensitivity and specificity
		 */
		tmp = iniparser_getdouble(ini, inla_string_join(secname, "LINK.INITIAL0"), 0.0);
		ds->link_fixed = Calloc(2, int);
		ds->link_fixed[0] = iniparser_getboolean(ini, inla_string_join(secname, "LINK.FIXED0"), 0);
		if (!ds->link_fixed[0] && mb->mode_use_mode) {
			tmp = mb->theta_file[mb->theta_counter_file++];
			if (mb->mode_fixed)
				ds->link_fixed[0] = 1;
		}
		ds->link_parameters = Calloc(1, Link_param_tp);
		ds->link_parameters->idx = -1;
		ds->link_parameters->order = -1;
		for (i = 0; i < n_data; i++) {
			ds->predictor_invlinkfunc_arg[i] = (void *) (ds->link_parameters);
		}
		HYPER_NEW(ds->link_parameters->sensitivity_intern, tmp);
		if (mb->verbose) {
			printf("\t\tinitialise sslogit link sensitivity_intern[%g]\n", ds->link_parameters->sensitivity_intern[0][0]);
			printf("\t\tfixed=[%1d]\n", ds->link_fixed[0]);
		}
		ds->link_prior = Calloc(2, Prior_tp);
		inla_read_prior_link0(mb, ini, sec, &(ds->link_prior[0]), "LOGITBETA", NULL);	/* read both priors here */
		inla_read_prior_link1(mb, ini, sec, &(ds->link_prior[1]), "LOGITBETA", NULL);	/* ... */

		/*
		 * add theta 
		 */
		if (!ds->link_fixed[0]) {
			mb->theta = Realloc(mb->theta, mb->ntheta + 1, double **);
			mb->theta_hyperid = Realloc(mb->theta_hyperid, mb->ntheta + 1, char *);
			mb->theta_hyperid[mb->ntheta] = ds->link_prior[0].hyperid;
			mb->theta_tag = Realloc(mb->theta_tag, mb->ntheta + 1, char *);
			mb->theta_tag_userscale = Realloc(mb->theta_tag_userscale, mb->ntheta + 1, char *);
			mb->theta_dir = Realloc(mb->theta_dir, mb->ntheta + 1, char *);
			mb->theta_tag[mb->ntheta] = inla_make_tag("Link sslogit sensitivity_intern", mb->ds);
			mb->theta_tag_userscale[mb->ntheta] = inla_make_tag("Link sslogit sensitivity", mb->ds);
			GMRFLib_sprintf(&msg, "%s-parameter0", secname);
			mb->theta_dir[mb->ntheta] = msg;

			mb->theta_from = Realloc(mb->theta_from, mb->ntheta + 1, char *);
			mb->theta_to = Realloc(mb->theta_to, mb->ntheta + 1, char *);
			mb->theta_from[mb->ntheta] = Strdup(ds->link_prior[0].from_theta);
			mb->theta_to[mb->ntheta] = Strdup(ds->link_prior[0].to_theta);
			mb->theta[mb->ntheta] = ds->link_parameters->sensitivity_intern;

			mb->theta_map = Realloc(mb->theta_map, mb->ntheta + 1, map_func_tp *);
			mb->theta_map[mb->ntheta] = map_probability;
			mb->theta_map_arg = Realloc(mb->theta_map_arg, mb->ntheta + 1, void *);
			mb->theta_map_arg[mb->ntheta] = NULL;
			mb->ntheta++;
			ds->link_ntheta++;
		}

		tmp = iniparser_getdouble(ini, inla_string_join(secname, "LINK.INITIAL1"), 0.0);
		ds->link_fixed[1] = iniparser_getboolean(ini, inla_string_join(secname, "LINK.FIXED1"), 0);
		if (!ds->link_fixed[1] && mb->mode_use_mode) {
			tmp = mb->theta_file[mb->theta_counter_file++];
			if (mb->mode_fixed)
				ds->link_fixed[1] = 1;
		}
		HYPER_NEW(ds->link_parameters->specificity_intern, tmp);
		if (mb->verbose) {
			printf("\t\tinitialise sslogit link specificity_intern[%g]\n", ds->link_parameters->specificity_intern[0][0]);
			printf("\t\tfixed=[%1d]\n", ds->link_fixed[1]);
		}
		if (!ds->link_fixed[1]) {
			mb->theta = Realloc(mb->theta, mb->ntheta + 1, double **);
			mb->theta_hyperid = Realloc(mb->theta_hyperid, mb->ntheta + 1, char *);
			mb->theta_hyperid[mb->ntheta] = ds->link_prior[1].hyperid;
			mb->theta_tag = Realloc(mb->theta_tag, mb->ntheta + 1, char *);
			mb->theta_tag_userscale = Realloc(mb->theta_tag_userscale, mb->ntheta + 1, char *);
			mb->theta_dir = Realloc(mb->theta_dir, mb->ntheta + 1, char *);
			mb->theta_tag[mb->ntheta] = inla_make_tag("Link sslogit specificity_intern", mb->ds);
			mb->theta_tag_userscale[mb->ntheta] = inla_make_tag("Link sslogit specificity", mb->ds);
			GMRFLib_sprintf(&msg, "%s-parameter0", secname);
			mb->theta_dir[mb->ntheta] = msg;

			mb->theta_from = Realloc(mb->theta_from, mb->ntheta + 1, char *);
			mb->theta_to = Realloc(mb->theta_to, mb->ntheta + 1, char *);
			mb->theta_from[mb->ntheta] = Strdup(ds->link_prior[1].from_theta);
			mb->theta_to[mb->ntheta] = Strdup(ds->link_prior[1].to_theta);
			mb->theta[mb->ntheta] = ds->link_parameters->specificity_intern;

			mb->theta_map = Realloc(mb->theta_map, mb->ntheta + 1, map_func_tp *);
			mb->theta_map[mb->ntheta] = map_probability;
			mb->theta_map_arg = Realloc(mb->theta_map_arg, mb->ntheta + 1, void *);
			mb->theta_map_arg[mb->ntheta] = NULL;
			mb->ntheta++;
			ds->link_ntheta++;
		}
	}
		break;

	case LINK_ROBIT:
	{
		/*
		 * Robit link with default fixed number of df.
		 */
		tmp = iniparser_getdouble(ini, inla_string_join(secname, "LINK.INITIAL"), 0.0);
		ds->link_fixed = Calloc(2, int);
		ds->link_fixed[0] = iniparser_getboolean(ini, inla_string_join(secname, "LINK.FIXED"), 1);
		if (!ds->link_fixed[0] && mb->mode_use_mode) {
			tmp = mb->theta_file[mb->theta_counter_file++];
			if (mb->mode_fixed)
				ds->link_fixed[0] = 1;
		}
		ds->link_parameters = Calloc(1, Link_param_tp);
		ds->link_parameters->idx = -1;
		ds->link_parameters->order = -1;
		for (i = 0; i < n_data; i++) {
			ds->predictor_invlinkfunc_arg[i] = (void *) (ds->link_parameters);
		}
		HYPER_NEW(ds->link_parameters->dof_intern, tmp);
		if (mb->verbose) {
			printf("\t\tinitialise robit link dof_intern[%g]\n", ds->link_parameters->dof_intern[0][0]);
			printf("\t\tfixed=[%1d]\n", ds->link_fixed[0]);
		}

		ds->link_prior = Calloc(1, Prior_tp);
		inla_read_prior_link(mb, ini, sec, &(ds->link_prior[0]), "PCDOF", NULL);

		/*
		 * add theta 
		 */
		if (!ds->link_fixed[0]) {
			mb->theta = Realloc(mb->theta, mb->ntheta + 1, double **);
			mb->theta_hyperid = Realloc(mb->theta_hyperid, mb->ntheta + 1, char *);
			mb->theta_hyperid[mb->ntheta] = ds->link_prior[0].hyperid;
			mb->theta_tag = Realloc(mb->theta_tag, mb->ntheta + 1, char *);
			mb->theta_tag_userscale = Realloc(mb->theta_tag_userscale, mb->ntheta + 1, char *);
			mb->theta_dir = Realloc(mb->theta_dir, mb->ntheta + 1, char *);
			mb->theta_tag[mb->ntheta] = inla_make_tag("Link robit dof_intern", mb->ds);
			mb->theta_tag_userscale[mb->ntheta] = inla_make_tag("Link robit dof", mb->ds);
			GMRFLib_sprintf(&msg, "%s-parameter", secname);
			mb->theta_dir[mb->ntheta] = msg;

			mb->theta_from = Realloc(mb->theta_from, mb->ntheta + 1, char *);
			mb->theta_to = Realloc(mb->theta_to, mb->ntheta + 1, char *);
			mb->theta_from[mb->ntheta] = Strdup(ds->link_prior[0].from_theta);
			mb->theta_to[mb->ntheta] = Strdup(ds->link_prior[0].to_theta);
			mb->theta[mb->ntheta] = ds->link_parameters->dof_intern;

			mb->theta_map = Realloc(mb->theta_map, mb->ntheta + 1, map_func_tp *);
			mb->theta_map[mb->ntheta] = map_dof;
			mb->theta_map_arg = Realloc(mb->theta_map_arg, mb->ntheta + 1, void *);
			mb->theta_map_arg[mb->ntheta] = NULL;
			mb->ntheta++;
			ds->link_ntheta++;
		}
	}
		break;

	case LINK_SN:
	{
		/*
		 * SN link
		 */

		ds->link_parameters = Calloc(1, Link_param_tp);
		ds->link_parameters->idx = -1;
		ds->link_parameters->order = -1;
		for (i = 0; i < n_data; i++) {
			ds->predictor_invlinkfunc_arg[i] = (void *) (ds->link_parameters);
		}

		ds->link_fixed = Calloc(2, int);
		tmp = iniparser_getdouble(ini, inla_string_join(secname, "LINK.INITIAL0"), 0.0);
		ds->link_fixed[0] = iniparser_getboolean(ini, inla_string_join(secname, "LINK.FIXED0"), 1);
		if (!ds->link_fixed[0] && mb->mode_use_mode) {
			tmp = mb->theta_file[mb->theta_counter_file++];
			if (mb->mode_fixed)
				ds->link_fixed[0] = 1;
		}

		HYPER_NEW(ds->link_parameters->sn_skew, tmp);
		if (mb->verbose) {
			printf("\t\tinitialise link_sn skew[%g]\n", ds->link_parameters->sn_skew[0][0]);
			printf("\t\tfixed=[%1d]\n", ds->link_fixed[0]);
		}

		ds->link_prior = Calloc(2, Prior_tp);
		inla_read_prior_link0(mb, ini, sec, &(ds->link_prior[0]), "PCSN", NULL);

		/*
		 * add theta 
		 */
		if (!ds->link_fixed[0]) {
			mb->theta = Realloc(mb->theta, mb->ntheta + 1, double **);
			mb->theta_hyperid = Realloc(mb->theta_hyperid, mb->ntheta + 1, char *);
			mb->theta_hyperid[mb->ntheta] = ds->link_prior[0].hyperid;
			mb->theta_tag = Realloc(mb->theta_tag, mb->ntheta + 1, char *);
			mb->theta_tag_userscale = Realloc(mb->theta_tag_userscale, mb->ntheta + 1, char *);
			mb->theta_dir = Realloc(mb->theta_dir, mb->ntheta + 1, char *);
			mb->theta_tag[mb->ntheta] = inla_make_tag("Link sn skew", mb->ds);
			mb->theta_tag_userscale[mb->ntheta] = inla_make_tag("Link sn skew", mb->ds);
			GMRFLib_sprintf(&msg, "%s-parameter", secname);
			mb->theta_dir[mb->ntheta] = msg;

			mb->theta_from = Realloc(mb->theta_from, mb->ntheta + 1, char *);
			mb->theta_to = Realloc(mb->theta_to, mb->ntheta + 1, char *);
			mb->theta_from[mb->ntheta] = Strdup(ds->link_prior[0].from_theta);
			mb->theta_to[mb->ntheta] = Strdup(ds->link_prior[0].to_theta);
			mb->theta[mb->ntheta] = ds->link_parameters->sn_skew;

			double *skewmax = Calloc(1, double);
			*skewmax = GMRFLib_SN_SKEWMAX;

			mb->theta_map = Realloc(mb->theta_map, mb->ntheta + 1, map_func_tp *);
			mb->theta_map[mb->ntheta] = map_phi;
			mb->theta_map_arg = Realloc(mb->theta_map_arg, mb->ntheta + 1, void *);
			mb->theta_map_arg[mb->ntheta] = (void *) skewmax;
			mb->ntheta++;
			ds->link_ntheta++;
		}

		tmp = iniparser_getdouble(ini, inla_string_join(secname, "LINK.INITIAL1"), 0);
		ds->link_fixed[1] = iniparser_getboolean(ini, inla_string_join(secname, "LINK.FIXED1"), 1);

		// special option. If 'initial=NA' or 'Inf', then remove intercept from the model. This is done setting fixed=1
		// and then recognising NAN in the map_invsn() function.
		if (ISNAN(tmp) || ISINF(tmp)) {
			tmp = NAN;
			ds->link_fixed[1] = 1;
		}

		if (!ds->link_fixed[1] && mb->mode_use_mode) {
			tmp = mb->theta_file[mb->theta_counter_file++];
			if (mb->mode_fixed)
				ds->link_fixed[1] = 1;
		}
		HYPER_NEW(ds->link_parameters->sn_intercept, tmp);
		if (mb->verbose) {
			printf("\t\tinitialise link_sn intercept[%g]\n", ds->link_parameters->sn_intercept[0][0]);
			if (ISNAN(ds->link_parameters->sn_intercept[0][0])) {
				printf("\t\t *** Intercept is removed from link-model\n");
			}
			printf("\t\tfixed=[%1d]\n", ds->link_fixed[1]);
		}
		inla_read_prior_link1(mb, ini, sec, &(ds->link_prior[1]), "LOGITBETA", NULL);

		/*
		 * add theta 
		 */
		if (!ds->link_fixed[1]) {
			mb->theta = Realloc(mb->theta, mb->ntheta + 1, double **);
			mb->theta_hyperid = Realloc(mb->theta_hyperid, mb->ntheta + 1, char *);
			mb->theta_hyperid[mb->ntheta] = ds->link_prior[1].hyperid;
			mb->theta_tag = Realloc(mb->theta_tag, mb->ntheta + 1, char *);
			mb->theta_tag_userscale = Realloc(mb->theta_tag_userscale, mb->ntheta + 1, char *);
			mb->theta_dir = Realloc(mb->theta_dir, mb->ntheta + 1, char *);
			mb->theta_tag[mb->ntheta] = inla_make_tag("Link sn intercept", mb->ds);
			mb->theta_tag_userscale[mb->ntheta] = inla_make_tag("Link sn intercept", mb->ds);
			GMRFLib_sprintf(&msg, "%s-parameter", secname);
			mb->theta_dir[mb->ntheta] = msg;

			mb->theta_from = Realloc(mb->theta_from, mb->ntheta + 1, char *);
			mb->theta_to = Realloc(mb->theta_to, mb->ntheta + 1, char *);
			mb->theta_from[mb->ntheta] = Strdup(ds->link_prior[1].from_theta);
			mb->theta_to[mb->ntheta] = Strdup(ds->link_prior[1].to_theta);
			mb->theta[mb->ntheta] = ds->link_parameters->sn_intercept;

			mb->theta_map = Realloc(mb->theta_map, mb->ntheta + 1, map_func_tp *);
			mb->theta_map[mb->ntheta] = map_probability;
			mb->theta_map_arg = Realloc(mb->theta_map_arg, mb->ntheta + 1, void *);
			mb->theta_map_arg[mb->ntheta] = NULL;
			mb->ntheta++;
			ds->link_ntheta++;
		}
	}
		break;

	case LINK_GEV:
	case LINK_CGEV:
	{
		char *name = (ds->link_id == LINK_GEV ? Strdup("gev") : Strdup("cgev"));
		ds->link_parameters = Calloc(1, Link_param_tp);
		ds->link_parameters->idx = -1;
		ds->link_parameters->order = -1;
		for (i = 0; i < n_data; i++) {
			ds->predictor_invlinkfunc_arg[i] = (void *) (ds->link_parameters);
		}

		ds->link_fixed = Calloc(2, int);
		tmp = iniparser_getdouble(ini, inla_string_join(secname, "LINK.INITIAL0"), 0.0);
		ds->link_fixed[0] = iniparser_getboolean(ini, inla_string_join(secname, "LINK.FIXED0"), 1);
		if (!ds->link_fixed[0] && mb->mode_use_mode) {
			tmp = mb->theta_file[mb->theta_counter_file++];
			if (mb->mode_fixed)
				ds->link_fixed[0] = 1;
		}

		HYPER_NEW(ds->link_parameters->bgev_tail, tmp);
		if (mb->verbose) {
			printf("\t\tinitialise link_%s_tail[%g]\n", name, ds->link_parameters->bgev_tail[0][0]);
			printf("\t\tfixed=[%1d]\n", ds->link_fixed[0]);
		}

		ds->link_prior = Calloc(2, Prior_tp);
		inla_read_prior_link0(mb, ini, sec, &(ds->link_prior[0]), "PCGEVTAIL", NULL);
		ds->link_parameters->bgev_tail_interval = Calloc(2, double);
		if (ds->link_prior[0].id == P_PC_GEVTAIL) {
			ds->link_parameters->bgev_tail_interval[0] = ds->link_prior[0].parameters[1];
			ds->link_parameters->bgev_tail_interval[1] = ds->link_prior[0].parameters[2];
		} else {
			// used a fixed interval then
			ds->link_parameters->bgev_tail_interval[0] = 0.0;
			ds->link_parameters->bgev_tail_interval[1] = 0.5;
		}

		if (DMIN(ds->link_parameters->bgev_tail_interval[0], ds->link_parameters->bgev_tail_interval[1]) < 0.0 ||
		    DMAX(ds->link_parameters->bgev_tail_interval[0], ds->link_parameters->bgev_tail_interval[1]) > 0.5 ||
		    ds->link_parameters->bgev_tail_interval[0] >= ds->link_parameters->bgev_tail_interval[1]) {
			inla_error_field_is_void(__GMRFLib_FuncName, secname, "BGEV.TAIL.INTERVAL", ctmp);
		}
		if (mb->verbose) {
			printf("\t\t%s.tail.interval [%g %g]\n", name, ds->link_parameters->bgev_tail_interval[0],
			       ds->link_parameters->bgev_tail_interval[1]);
		}

		if (!ds->link_fixed[0]) {
			mb->theta = Realloc(mb->theta, mb->ntheta + 1, double **);
			mb->theta_hyperid = Realloc(mb->theta_hyperid, mb->ntheta + 1, char *);
			mb->theta_hyperid[mb->ntheta] = ds->link_prior[0].hyperid;
			mb->theta_tag = Realloc(mb->theta_tag, mb->ntheta + 1, char *);
			mb->theta_tag_userscale = Realloc(mb->theta_tag_userscale, mb->ntheta + 1, char *);
			mb->theta_dir = Realloc(mb->theta_dir, mb->ntheta + 1, char *);
			if (ds->link_id == LINK_GEV) {
				mb->theta_tag[mb->ntheta] = inla_make_tag("Link gev tail_intern", mb->ds);
				mb->theta_tag_userscale[mb->ntheta] = inla_make_tag("Link gev tail", mb->ds);
			} else {
				mb->theta_tag[mb->ntheta] = inla_make_tag("Link cgev tail_intern", mb->ds);
				mb->theta_tag_userscale[mb->ntheta] = inla_make_tag("Link cgev tail", mb->ds);
			}
			GMRFLib_sprintf(&msg, "%s-parameter", secname);
			mb->theta_dir[mb->ntheta] = msg;

			mb->theta_from = Realloc(mb->theta_from, mb->ntheta + 1, char *);
			mb->theta_to = Realloc(mb->theta_to, mb->ntheta + 1, char *);
			mb->theta_from[mb->ntheta] = Strdup(ds->link_prior[0].from_theta);
			mb->theta_to[mb->ntheta] = Strdup(ds->link_prior[0].to_theta);

			mb->theta[mb->ntheta] = ds->link_parameters->bgev_tail;
			mb->theta_map = Realloc(mb->theta_map, mb->ntheta + 1, map_func_tp *);
			mb->theta_map[mb->ntheta] = map_interval;
			mb->theta_map_arg = Realloc(mb->theta_map_arg, mb->ntheta + 1, void *);
			mb->theta_map_arg[mb->ntheta] = (void *) (ds->link_parameters->bgev_tail_interval);
			mb->ntheta++;
			ds->link_ntheta++;
		}

		tmp = iniparser_getdouble(ini, inla_string_join(secname, "LINK.INITIAL1"), 0);
		ds->link_fixed[1] = iniparser_getboolean(ini, inla_string_join(secname, "LINK.FIXED1"), 1);

		// special option. If 'initial=NA' or 'Inf', then remove intercept from the model. This is done setting fixed=1
		// and then recognising NAN in the map_invsn() function.
		if (ISNAN(tmp) || ISINF(tmp)) {
			tmp = NAN;
			ds->link_fixed[1] = 1;
		}

		if (!ds->link_fixed[1] && mb->mode_use_mode) {
			tmp = mb->theta_file[mb->theta_counter_file++];
			if (mb->mode_fixed)
				ds->link_fixed[1] = 1;
		}
		HYPER_NEW(ds->link_parameters->bgev_intercept, tmp);
		if (mb->verbose) {
			printf("\t\tinitialise link_%s intercept[%g]\n", name, ds->link_parameters->bgev_intercept[0][0]);
			if (ISNAN(ds->link_parameters->bgev_intercept[0][0])) {
				printf("\t\t *** Intercept is removed from link-model\n");
			}
			printf("\t\tfixed=[%1d]\n", ds->link_fixed[1]);
		}
		inla_read_prior_link1(mb, ini, sec, &(ds->link_prior[1]), "NORMAL", NULL);

		/*
		 * add theta 
		 */
		if (!ds->link_fixed[1]) {
			mb->theta = Realloc(mb->theta, mb->ntheta + 1, double **);
			mb->theta_hyperid = Realloc(mb->theta_hyperid, mb->ntheta + 1, char *);
			mb->theta_hyperid[mb->ntheta] = ds->link_prior[1].hyperid;
			mb->theta_tag = Realloc(mb->theta_tag, mb->ntheta + 1, char *);
			mb->theta_tag_userscale = Realloc(mb->theta_tag_userscale, mb->ntheta + 1, char *);
			mb->theta_dir = Realloc(mb->theta_dir, mb->ntheta + 1, char *);

			if (ds->link_id == LINK_GEV) {
				mb->theta_tag[mb->ntheta] = inla_make_tag("Link gev intercept_intern", mb->ds);
				mb->theta_tag_userscale[mb->ntheta] = inla_make_tag("Link gev intercept", mb->ds);
			} else {
				mb->theta_tag[mb->ntheta] = inla_make_tag("Link cgev intercept_intern", mb->ds);
				mb->theta_tag_userscale[mb->ntheta] = inla_make_tag("Link cgev intercept", mb->ds);
			}

			GMRFLib_sprintf(&msg, "%s-parameter", secname);
			mb->theta_dir[mb->ntheta] = msg;

			mb->theta_from = Realloc(mb->theta_from, mb->ntheta + 1, char *);
			mb->theta_to = Realloc(mb->theta_to, mb->ntheta + 1, char *);
			mb->theta_from[mb->ntheta] = Strdup(ds->link_prior[1].from_theta);
			mb->theta_to[mb->ntheta] = Strdup(ds->link_prior[1].to_theta);
			mb->theta[mb->ntheta] = ds->link_parameters->bgev_intercept;

			mb->theta_map = Realloc(mb->theta_map, mb->ntheta + 1, map_func_tp *);
			mb->theta_map[mb->ntheta] = map_probability;
			mb->theta_map_arg = Realloc(mb->theta_map_arg, mb->ntheta + 1, void *);
			mb->theta_map_arg[mb->ntheta] = NULL;
			mb->ntheta++;
			ds->link_ntheta++;
		}
	}
		break;

	case LINK_POWER_LOGIT:
	{
		/*
		 * power logit link
		 */

		ds->link_parameters = Calloc(1, Link_param_tp);
		ds->link_parameters->idx = -1;
		ds->link_parameters->order = -1;
		for (i = 0; i < n_data; i++) {
			ds->predictor_invlinkfunc_arg[i] = (void *) (ds->link_parameters);
		}

		ds->link_fixed = Calloc(2, int);
		tmp = iniparser_getdouble(ini, inla_string_join(secname, "LINK.INITIAL0"), 0.0);
		ds->link_fixed[0] = iniparser_getboolean(ini, inla_string_join(secname, "LINK.FIXED0"), 1);
		if (!ds->link_fixed[0] && mb->mode_use_mode) {
			tmp = mb->theta_file[mb->theta_counter_file++];
			if (mb->mode_fixed)
				ds->link_fixed[0] = 1;
		}

		HYPER_NEW(ds->link_parameters->power_intern, tmp);
		if (mb->verbose) {
			printf("\t\tinitialise link_power_logit power[%g]\n", ds->link_parameters->power_intern[0][0]);
			printf("\t\tfixed=[%1d]\n", ds->link_fixed[0]);
		}

		ds->link_prior = Calloc(2, Prior_tp);
		inla_read_prior_link0(mb, ini, sec, &(ds->link_prior[0]), "NORMAL", NULL);

		/*
		 * add theta 
		 */
		if (!ds->link_fixed[0]) {
			mb->theta = Realloc(mb->theta, mb->ntheta + 1, double **);
			mb->theta_hyperid = Realloc(mb->theta_hyperid, mb->ntheta + 1, char *);
			mb->theta_hyperid[mb->ntheta] = ds->link_prior[0].hyperid;
			mb->theta_tag = Realloc(mb->theta_tag, mb->ntheta + 1, char *);
			mb->theta_tag_userscale = Realloc(mb->theta_tag_userscale, mb->ntheta + 1, char *);
			mb->theta_dir = Realloc(mb->theta_dir, mb->ntheta + 1, char *);
			mb->theta_tag[mb->ntheta] = inla_make_tag("Link power_logit power.intern", mb->ds);
			mb->theta_tag_userscale[mb->ntheta] = inla_make_tag("Link power_logit power", mb->ds);
			GMRFLib_sprintf(&msg, "%s-parameter", secname);
			mb->theta_dir[mb->ntheta] = msg;

			mb->theta_from = Realloc(mb->theta_from, mb->ntheta + 1, char *);
			mb->theta_to = Realloc(mb->theta_to, mb->ntheta + 1, char *);
			mb->theta_from[mb->ntheta] = Strdup(ds->link_prior[0].from_theta);
			mb->theta_to[mb->ntheta] = Strdup(ds->link_prior[0].to_theta);
			mb->theta[mb->ntheta] = ds->link_parameters->power_intern;

			mb->theta_map = Realloc(mb->theta_map, mb->ntheta + 1, map_func_tp *);
			mb->theta_map[mb->ntheta] = map_exp;
			mb->ntheta++;
			ds->link_ntheta++;
		}

		tmp = iniparser_getdouble(ini, inla_string_join(secname, "LINK.INITIAL1"), 0);
		ds->link_fixed[1] = iniparser_getboolean(ini, inla_string_join(secname, "LINK.FIXED1"), 1);

		// special option. If 'initial=NA' or 'Inf', then remove intercept from the model. This is done setting fixed=1
		// and then recognising NAN in the map_invsn() function.
		if (ISNAN(tmp) || ISINF(tmp)) {
			tmp = NAN;
			ds->link_fixed[1] = 1;
		}

		if (!ds->link_fixed[1] && mb->mode_use_mode) {
			tmp = mb->theta_file[mb->theta_counter_file++];
			if (mb->mode_fixed)
				ds->link_fixed[1] = 1;
		}
		HYPER_NEW(ds->link_parameters->intercept_intern, tmp);
		if (mb->verbose) {
			printf("\t\tinitialise link_power_logit intercept.intern[%g]\n", ds->link_parameters->intercept_intern[0][0]);
			if (ISNAN(ds->link_parameters->intercept_intern[0][0])) {
				printf("\t\t *** Intercept is removed from link-model\n");
			}
			printf("\t\tfixed=[%1d]\n", ds->link_fixed[1]);
		}
		inla_read_prior_link1(mb, ini, sec, &(ds->link_prior[1]), "LOGITBETA", NULL);

		/*
		 * add theta 
		 */
		if (!ds->link_fixed[1]) {
			mb->theta = Realloc(mb->theta, mb->ntheta + 1, double **);
			mb->theta_hyperid = Realloc(mb->theta_hyperid, mb->ntheta + 1, char *);
			mb->theta_hyperid[mb->ntheta] = ds->link_prior[1].hyperid;
			mb->theta_tag = Realloc(mb->theta_tag, mb->ntheta + 1, char *);
			mb->theta_tag_userscale = Realloc(mb->theta_tag_userscale, mb->ntheta + 1, char *);
			mb->theta_dir = Realloc(mb->theta_dir, mb->ntheta + 1, char *);
			mb->theta_tag[mb->ntheta] = inla_make_tag("Link power_logit intercept_intern", mb->ds);
			mb->theta_tag_userscale[mb->ntheta] = inla_make_tag("Link power_logit intercept", mb->ds);
			GMRFLib_sprintf(&msg, "%s-parameter", secname);
			mb->theta_dir[mb->ntheta] = msg;

			mb->theta_from = Realloc(mb->theta_from, mb->ntheta + 1, char *);
			mb->theta_to = Realloc(mb->theta_to, mb->ntheta + 1, char *);
			mb->theta_from[mb->ntheta] = Strdup(ds->link_prior[1].from_theta);
			mb->theta_to[mb->ntheta] = Strdup(ds->link_prior[1].to_theta);
			mb->theta[mb->ntheta] = ds->link_parameters->intercept_intern;

			mb->theta_map = Realloc(mb->theta_map, mb->ntheta + 1, map_func_tp *);
			mb->theta_map[mb->ntheta] = map_probability;
			mb->theta_map_arg = Realloc(mb->theta_map_arg, mb->ntheta + 1, void *);
			mb->theta_map_arg[mb->ntheta] = NULL;
			mb->ntheta++;
			ds->link_ntheta++;
		}
	}
		break;

	case LINK_LOGOFFSET:
	{
		/*
		 * exp(beta)*cov + exp(linear.predictor)
		 */
		tmp = iniparser_getdouble(ini, inla_string_join(secname, "LINK.INITIAL"), 0.0);
		ds->link_fixed = Calloc(1, int);
		ds->link_fixed[0] = iniparser_getboolean(ini, inla_string_join(secname, "LINK.FIXED"), 0);
		if (!ds->link_fixed[0] && mb->mode_use_mode) {
			tmp = mb->theta_file[mb->theta_counter_file++];
			if (mb->mode_fixed)
				ds->link_fixed[0] = 1;
		}
		ds->link_parameters = Calloc(1, Link_param_tp);
		ds->link_parameters->idx = -1;
		ds->link_parameters->order = ds->link_order;
		ds->link_parameters->variant = ds->link_variant;
		for (i = 0; i < n_data; i++) {
			ds->predictor_invlinkfunc_arg[i] = (void *) (ds->link_parameters);
		}
		HYPER_NEW(ds->link_parameters->beta_intern, tmp);
		if (mb->verbose) {
			printf("\t\tinitialise link beta[%g]\n", ds->link_parameters->beta_intern[0][0]);
			printf("\t\tfixed=[%1d]\n", ds->link_fixed[0]);
		}
		ds->link_prior = Calloc(1, Prior_tp);
		inla_read_prior_link(mb, ini, sec, ds->link_prior, "GAUSSIAN-std", NULL);

		/*
		 * add theta 
		 */
		if (!ds->link_fixed[0]) {
			mb->theta = Realloc(mb->theta, mb->ntheta + 1, double **);
			mb->theta_hyperid = Realloc(mb->theta_hyperid, mb->ntheta + 1, char *);
			mb->theta_hyperid[mb->ntheta] = ds->link_prior[0].hyperid;
			mb->theta_tag = Realloc(mb->theta_tag, mb->ntheta + 1, char *);
			mb->theta_tag_userscale = Realloc(mb->theta_tag_userscale, mb->ntheta + 1, char *);
			mb->theta_dir = Realloc(mb->theta_dir, mb->ntheta + 1, char *);
			mb->theta_tag[mb->ntheta] = inla_make_tag("Link beta_intern for logoffset", mb->ds);
			mb->theta_tag_userscale[mb->ntheta] = inla_make_tag("Link beta for logoffset", mb->ds);
			GMRFLib_sprintf(&msg, "%s-parameter0", secname);
			mb->theta_dir[mb->ntheta] = msg;

			mb->theta_from = Realloc(mb->theta_from, mb->ntheta + 1, char *);
			mb->theta_to = Realloc(mb->theta_to, mb->ntheta + 1, char *);
			mb->theta_from[mb->ntheta] = Strdup(ds->link_prior[0].from_theta);
			mb->theta_to[mb->ntheta] = Strdup(ds->link_prior[0].to_theta);
			mb->theta[mb->ntheta] = ds->link_parameters->beta_intern;

			mb->theta_map = Realloc(mb->theta_map, mb->ntheta + 1, map_func_tp *);
			mb->theta_map[mb->ntheta] = map_exp;
			mb->theta_map_arg = Realloc(mb->theta_map_arg, mb->ntheta + 1, void *);
			mb->theta_map_arg[mb->ntheta] = NULL;
			mb->ntheta++;
			ds->link_ntheta++;
		}
	}
		break;

	case LINK_LOGITOFFSET:
	{
		/*
		 * p + (1-p) * exp(linear.predictor)/(1 + exp(linear.predictor))
		 */
		tmp = iniparser_getdouble(ini, inla_string_join(secname, "LINK.INITIAL"), 0.0);
		ds->link_fixed = Calloc(1, int);
		ds->link_fixed[0] = iniparser_getboolean(ini, inla_string_join(secname, "LINK.FIXED"), 0);
		if (!ds->link_fixed[0] && mb->mode_use_mode) {
			tmp = mb->theta_file[mb->theta_counter_file++];
			if (mb->mode_fixed)
				ds->link_fixed[0] = 1;
		}
		ds->link_parameters = Calloc(1, Link_param_tp);
		ds->link_parameters->idx = -1;
		ds->link_parameters->order = ds->link_order;
		ds->link_parameters->variant = ds->link_variant;
		for (i = 0; i < n_data; i++) {
			ds->predictor_invlinkfunc_arg[i] = (void *) (ds->link_parameters);
		}
		HYPER_NEW(ds->link_parameters->prob_intern, tmp);
		if (mb->verbose) {
			printf("\t\tinitialise link prob[%g]\n", ds->link_parameters->prob_intern[0][0]);
			printf("\t\tfixed=[%1d]\n", ds->link_fixed[0]);
		}
		ds->link_prior = Calloc(1, Prior_tp);
		inla_read_prior_link(mb, ini, sec, ds->link_prior, "GAUSSIAN-std", NULL);

		/*
		 * add theta 
		 */
		if (!ds->link_fixed[0]) {
			mb->theta = Realloc(mb->theta, mb->ntheta + 1, double **);
			mb->theta_hyperid = Realloc(mb->theta_hyperid, mb->ntheta + 1, char *);
			mb->theta_hyperid[mb->ntheta] = ds->link_prior[0].hyperid;
			mb->theta_tag = Realloc(mb->theta_tag, mb->ntheta + 1, char *);
			mb->theta_tag_userscale = Realloc(mb->theta_tag_userscale, mb->ntheta + 1, char *);
			mb->theta_dir = Realloc(mb->theta_dir, mb->ntheta + 1, char *);
			mb->theta_tag[mb->ntheta] = inla_make_tag("Link prob_intern for logitoffset", mb->ds);
			mb->theta_tag_userscale[mb->ntheta] = inla_make_tag("Link prob for logitoffset", mb->ds);
			GMRFLib_sprintf(&msg, "%s-parameter0", secname);
			mb->theta_dir[mb->ntheta] = msg;

			mb->theta_from = Realloc(mb->theta_from, mb->ntheta + 1, char *);
			mb->theta_to = Realloc(mb->theta_to, mb->ntheta + 1, char *);
			mb->theta_from[mb->ntheta] = Strdup(ds->link_prior[0].from_theta);
			mb->theta_to[mb->ntheta] = Strdup(ds->link_prior[0].to_theta);
			mb->theta[mb->ntheta] = ds->link_parameters->prob_intern;

			mb->theta_map = Realloc(mb->theta_map, mb->ntheta + 1, map_func_tp *);
			mb->theta_map[mb->ntheta] = map_probability;
			mb->theta_map_arg = Realloc(mb->theta_map_arg, mb->ntheta + 1, void *);
			mb->theta_map_arg[mb->ntheta] = NULL;
			mb->ntheta++;
			ds->link_ntheta++;
		}
	}
		break;

	case LINK_TEST1:
	{
		/*
		 * exp(eta - beta*cov
		 */
		tmp = iniparser_getdouble(ini, inla_string_join(secname, "LINK.INITIAL"), 0.0);
		ds->link_fixed = Calloc(1, int);
		ds->link_fixed[0] = iniparser_getboolean(ini, inla_string_join(secname, "LINK.FIXED"), 0);
		if (!ds->link_fixed[0] && mb->mode_use_mode) {
			tmp = mb->theta_file[mb->theta_counter_file++];
			if (mb->mode_fixed)
				ds->link_fixed[0] = 1;
		}
		ds->link_parameters = Calloc(1, Link_param_tp);
		ds->link_parameters->idx = -1;
		ds->link_parameters->order = ds->link_order;
		for (i = 0; i < n_data; i++) {
			ds->predictor_invlinkfunc_arg[i] = (void *) (ds->link_parameters);
		}
		HYPER_NEW(ds->link_parameters->beta, tmp);
		if (mb->verbose) {
			printf("\t\tinitialise link beta[%g]\n", ds->link_parameters->beta[0][0]);
			printf("\t\tfixed=[%1d]\n", ds->link_fixed[0]);
		}
		ds->link_prior = Calloc(1, Prior_tp);
		inla_read_prior_link(mb, ini, sec, ds->link_prior, "GAUSSIAN-std", NULL);

		/*
		 * add theta 
		 */
		if (!ds->link_fixed[0]) {
			mb->theta = Realloc(mb->theta, mb->ntheta + 1, double **);
			mb->theta_hyperid = Realloc(mb->theta_hyperid, mb->ntheta + 1, char *);
			mb->theta_hyperid[mb->ntheta] = ds->link_prior[0].hyperid;
			mb->theta_tag = Realloc(mb->theta_tag, mb->ntheta + 1, char *);
			mb->theta_tag_userscale = Realloc(mb->theta_tag_userscale, mb->ntheta + 1, char *);
			mb->theta_dir = Realloc(mb->theta_dir, mb->ntheta + 1, char *);
			mb->theta_tag[mb->ntheta] = inla_make_tag("Link beta for test1", mb->ds);
			mb->theta_tag_userscale[mb->ntheta] = inla_make_tag("Link beta for test1", mb->ds);
			GMRFLib_sprintf(&msg, "%s-parameter0", secname);
			mb->theta_dir[mb->ntheta] = msg;

			mb->theta_from = Realloc(mb->theta_from, mb->ntheta + 1, char *);
			mb->theta_to = Realloc(mb->theta_to, mb->ntheta + 1, char *);
			mb->theta_from[mb->ntheta] = Strdup(ds->link_prior[0].from_theta);
			mb->theta_to[mb->ntheta] = Strdup(ds->link_prior[0].to_theta);
			mb->theta[mb->ntheta] = ds->link_parameters->beta;

			mb->theta_map = Realloc(mb->theta_map, mb->ntheta + 1, map_func_tp *);
			mb->theta_map[mb->ntheta] = map_identity;
			mb->theta_map_arg = Realloc(mb->theta_map_arg, mb->ntheta + 1, void *);
			mb->theta_map_arg[mb->ntheta] = NULL;
			mb->ntheta++;
			ds->link_ntheta++;
		}
	}
		break;

	case LINK_SPECIAL2:
	{
		/*
		 * exp(eta) * ( 1-x + x*exp(beta) )
		 */
		tmp = iniparser_getdouble(ini, inla_string_join(secname, "LINK.INITIAL"), 0.0);
		ds->link_fixed = Calloc(1, int);
		ds->link_fixed[0] = iniparser_getboolean(ini, inla_string_join(secname, "LINK.FIXED"), 0);
		if (!ds->link_fixed[0] && mb->mode_use_mode) {
			tmp = mb->theta_file[mb->theta_counter_file++];
			if (mb->mode_fixed)
				ds->link_fixed[0] = 1;
		}
		ds->link_parameters = Calloc(1, Link_param_tp);
		ds->link_parameters->idx = -1;
		ds->link_parameters->order = ds->link_order;
		for (i = 0; i < n_data; i++) {
			ds->predictor_invlinkfunc_arg[i] = (void *) (ds->link_parameters);
		}
		HYPER_NEW(ds->link_parameters->beta, tmp);
		if (mb->verbose) {
			printf("\t\tinitialise link beta[%g]\n", ds->link_parameters->beta[0][0]);
			printf("\t\tfixed=[%1d]\n", ds->link_fixed[0]);
		}
		ds->link_prior = Calloc(1, Prior_tp);
		inla_read_prior_link(mb, ini, sec, ds->link_prior, "GAUSSIAN-std", NULL);

		/*
		 * add theta 
		 */
		if (!ds->link_fixed[0]) {
			mb->theta = Realloc(mb->theta, mb->ntheta + 1, double **);
			mb->theta_hyperid = Realloc(mb->theta_hyperid, mb->ntheta + 1, char *);
			mb->theta_hyperid[mb->ntheta] = ds->link_prior[0].hyperid;
			mb->theta_tag = Realloc(mb->theta_tag, mb->ntheta + 1, char *);
			mb->theta_tag_userscale = Realloc(mb->theta_tag_userscale, mb->ntheta + 1, char *);
			mb->theta_dir = Realloc(mb->theta_dir, mb->ntheta + 1, char *);
			mb->theta_tag[mb->ntheta] = inla_make_tag("Link beta for special2", mb->ds);
			mb->theta_tag_userscale[mb->ntheta] = inla_make_tag("Link beta for special2", mb->ds);
			GMRFLib_sprintf(&msg, "%s-parameter0", secname);
			mb->theta_dir[mb->ntheta] = msg;

			mb->theta_from = Realloc(mb->theta_from, mb->ntheta + 1, char *);
			mb->theta_to = Realloc(mb->theta_to, mb->ntheta + 1, char *);
			mb->theta_from[mb->ntheta] = Strdup(ds->link_prior[0].from_theta);
			mb->theta_to[mb->ntheta] = Strdup(ds->link_prior[0].to_theta);
			mb->theta[mb->ntheta] = ds->link_parameters->beta;

			mb->theta_map = Realloc(mb->theta_map, mb->ntheta + 1, map_func_tp *);
			mb->theta_map[mb->ntheta] = map_identity;
			mb->theta_map_arg = Realloc(mb->theta_map_arg, mb->ntheta + 1, void *);
			mb->theta_map_arg[mb->ntheta] = NULL;
			mb->ntheta++;
			ds->link_ntheta++;
		}
	}
		break;

	case LINK_SPECIAL1:
	{
		ds->link_prior = Calloc(2, Prior_tp);
		inla_read_prior_link0(mb, ini, sec, &(ds->link_prior[0]), "LOGGAMMA", NULL);	// log precision
		inla_read_prior_link1(mb, ini, sec, &(ds->link_prior[1]), "MVNORM", NULL);	// the beta's

		if (ds->link_order > 0 && (int) ds->link_prior[1].parameters[0] != ds->link_order) {
			char *ptmp = NULL;
			GMRFLib_sprintf(&ptmp,
					"Dimension of the MVNORM prior is not equal to the order of the link-model: %1d != %1d\n",
					(int) ds->link_prior[1].parameters[0], ds->link_order);
			inla_error_general(ptmp);
			exit(EXIT_FAILURE);
		}
		if (ds->link_order > ds->link_covariates->ncol) {
			char *ptmp = NULL;
			GMRFLib_sprintf(&ptmp, "The link-model %s require more covariates : %1d > %1d\n",
					ds->link_model, ds->link_order, ds->link_covariates->ncol);
			inla_error_general(ptmp);
			exit(EXIT_FAILURE);
		}
		ds->link_ntheta = ds->link_order + 1;
		assert(ds->link_ntheta <= LINK_MAXTHETA + 1);
		if (mb->verbose) {
			printf("\t\tlink ntheta = [%1d]\n", ds->link_ntheta);
		}

		/*
		 * mark all possible as read 
		 */

		// mark all as read
		for (i = 0; i < LINK_MAXTHETA + 1; i++) {
			for (j = 0; j < keywords_len; j++) {
				GMRFLib_sprintf(&ctmp, "LINK.%s%1d", keywords[j], i);
				iniparser_getstring(ini, inla_string_join(secname, ctmp), NULL);
				Free(ctmp);
			}
		}

		ds->link_fixed = Calloc(ds->link_ntheta, int);
		ds->link_initial = Calloc(ds->link_ntheta, double);
		ds->link_parameters = Calloc(1, Link_param_tp);
		ds->link_parameters->idx = -1;
		ds->link_parameters->order = ds->link_order;
		for (i = 0; i < n_data; i++) {
			ds->predictor_invlinkfunc_arg = (void **) ((void *) (ds->link_parameters));
		}
		HYPER_NEW(ds->link_parameters->log_prec, 0.0);
		ds->link_parameters->betas = Calloc(LINK_MAXTHETA, double **);
		for (i = 0; i < LINK_MAXTHETA; i++) {
			HYPER_NEW(ds->link_parameters->betas[i], 0.0);
		}

		/*
		 * then read those we need 
		 */
		int ignore_prior_error = 0;
		for (i = 0; i < ds->link_ntheta; i++) {
			double theta_initial = 0;
			GMRFLib_sprintf(&ctmp, "LINK.FIXED%1d", i);
			ds->link_fixed[i] = iniparser_getboolean(ini, inla_string_join(secname, ctmp), 0);
			GMRFLib_sprintf(&ctmp, "LINK.INITIAL%1d", i);
			theta_initial = iniparser_getdouble(ini, inla_string_join(secname, ctmp), theta_initial);

			if (!ds->link_fixed[i] && mb->mode_use_mode) {
				theta_initial = mb->theta_file[mb->theta_counter_file++];
				if (mb->mode_fixed) {
					ds->link_fixed[i] = 1;
					ignore_prior_error = 1;
				}
			}

			if (mb->verbose) {
				printf("\t\tlink theta fixed  [%1d] = [%1d]\n", i, ds->link_fixed[i]);
				printf("\t\tlink theta initial[%1d] = [%g]\n", i, ds->link_initial[i]);
			}

			if (i == 0) {
				/*
				 * precision 
				 */
				HYPER_INIT(ds->link_parameters->log_prec, theta_initial);
				if (mb->verbose) {
					printf("\t\tlink initialise log_prec = [%g]\n", theta_initial);
					printf("\t\tlink fixed[%1d]=[%1d]\n", i, ds->link_fixed[i]);
				}

				if (!ds->link_fixed[i]) {
					/*
					 * add this \theta 
					 */
					mb->theta = Realloc(mb->theta, mb->ntheta + 1, double **);
					mb->theta_hyperid = Realloc(mb->theta_hyperid, mb->ntheta + 1, char *);
					mb->theta_hyperid[mb->ntheta] = ds->link_prior[0].hyperid;
					mb->theta_tag = Realloc(mb->theta_tag, mb->ntheta + 1, char *);
					mb->theta_tag_userscale = Realloc(mb->theta_tag_userscale, mb->ntheta + 1, char *);
					mb->theta_dir = Realloc(mb->theta_dir, mb->ntheta + 1, char *);
					GMRFLib_sprintf(&msg, "Link log precision for %s", ds->data_likelihood);
					mb->theta_tag[mb->ntheta] = inla_make_tag(msg, mb->ds);
					GMRFLib_sprintf(&msg, "Link precision for %s", ds->data_likelihood);
					mb->theta_tag_userscale[mb->ntheta] = msg;
					GMRFLib_sprintf(&msg, "%s-parameter%1d", msg, i);
					mb->theta_dir[mb->ntheta] = inla_make_tag(msg, mb->ds);

					mb->theta_from = Realloc(mb->theta_from, mb->ntheta + 1, char *);
					mb->theta_to = Realloc(mb->theta_to, mb->ntheta + 1, char *);
					mb->theta_from[mb->ntheta] = Strdup(ds->link_prior[0].from_theta);
					mb->theta_to[mb->ntheta] = Strdup(ds->link_prior[0].to_theta);
					mb->theta[mb->ntheta] = ds->link_parameters->log_prec;
					mb->theta_map = Realloc(mb->theta_map, mb->ntheta + 1, map_func_tp *);
					mb->theta_map[mb->ntheta] = map_precision;
					mb->theta_map_arg = Realloc(mb->theta_map_arg, mb->ntheta + 1, void *);
					mb->theta_map_arg[mb->ntheta] = NULL;
					mb->ntheta++;
				}
			} else {
				/*
				 * beta
				 */
				HYPER_INIT(ds->link_parameters->betas[i - 1], theta_initial);
				if (mb->verbose) {
					printf("\t\tlink initialise beta[%1d]=[%g]\n", i, theta_initial);
					printf("\t\tlink fixed[%1d]=[%1d]\n", i, ds->link_fixed[i]);
				}

				if (!ds->link_fixed[i]) {
					/*
					 * add this \theta 
					 */
					mb->theta = Realloc(mb->theta, mb->ntheta + 1, double **);
					mb->theta_hyperid = Realloc(mb->theta_hyperid, mb->ntheta + 1, char *);
					mb->theta_hyperid[mb->ntheta] = NULL;	/* multivariate normal */
					mb->theta_tag = Realloc(mb->theta_tag, mb->ntheta + 1, char *);
					mb->theta_tag_userscale = Realloc(mb->theta_tag_userscale, mb->ntheta + 1, char *);
					mb->theta_dir = Realloc(mb->theta_dir, mb->ntheta + 1, char *);
					GMRFLib_sprintf(&msg, "Link beta%1d for %s", i, ds->data_likelihood);
					mb->theta_tag[mb->ntheta] = inla_make_tag(msg, mb->ds);
					GMRFLib_sprintf(&msg, "Link beta%1d for %s", i, ds->data_likelihood);
					mb->theta_tag_userscale[mb->ntheta] = msg;
					GMRFLib_sprintf(&msg, "%s-parameter%1d", msg, i);
					mb->theta_dir[mb->ntheta] = inla_make_tag(msg, mb->ds);

					mb->theta_from = Realloc(mb->theta_from, mb->ntheta + 1, char *);
					mb->theta_to = Realloc(mb->theta_to, mb->ntheta + 1, char *);
					mb->theta_from[mb->ntheta] = Strdup(ds->link_prior[1].from_theta);
					mb->theta_to[mb->ntheta] = Strdup(ds->link_prior[1].to_theta);
					mb->theta[mb->ntheta] = ds->link_parameters->betas[i - 1];	/* yes! */
					mb->theta_map = Realloc(mb->theta_map, mb->ntheta + 1, map_func_tp *);
					mb->theta_map[mb->ntheta] = map_identity;
					mb->theta_map_arg = Realloc(mb->theta_map_arg, mb->ntheta + 1, void *);
					mb->theta_map_arg[mb->ntheta] = NULL;
					mb->ntheta++;
				}
			}
		}

		if (ds->link_order > 0 && (int) ds->link_prior[1].parameters[0] != ds->link_order && !ignore_prior_error) {
			char *ptmp = NULL;
			GMRFLib_sprintf(&ptmp,
					"Dimension of the MVNORM prior is not equal to the order of the link-model: %1d != %1d\n",
					(int) ds->link_prior[1].parameters[0], ds->link_order);
			inla_error_general(ptmp);
			exit(EXIT_FAILURE);
		}
	}
		break;

	default:
	{
		GMRFLib_ASSERT(0 == 1, GMRFLib_ESNH);
	}
		break;
	}

	/*
	 * mixing of the likelihood 
	 */
	ds->mix_use = iniparser_getboolean(ini, inla_string_join(secname, "MIX.USE"), 0);
	if (mb->verbose) {
		printf("\t\tmix.use[%1d]\n", ds->mix_use);
	}
	if (ds->mix_use) {
		/*
		 * read mix-parameters
		 */

		char *model = NULL, *integrator = NULL;

		ds->mix_npoints = iniparser_getint(ini, inla_string_join(secname, "MIX.NPOINTS"), 99);
		assert(ds->mix_npoints >= 5);
		model = Strdup(strupc(iniparser_getstring(ini, inla_string_join(secname, "MIX.MODEL"), NULL)));
		integrator = Strdup(strupc(iniparser_getstring(ini, inla_string_join(secname, "MIX.INTEGRATOR"), NULL)));
		if (!strcasecmp(integrator, "default")) {
			ds->mix_integrator = MIX_INT_DEFAULT;
		} else if (!strcasecmp(integrator, "quadrature")) {
			ds->mix_integrator = MIX_INT_QUADRATURE;
		} else if (!strcasecmp(integrator, "simpson")) {
			ds->mix_integrator = MIX_INT_SIMPSON;
		} else {
			assert(0 == 1);
		}

		if (mb->verbose) {
			printf("\t\tmix.npoints [%d]\n", ds->mix_npoints);
			printf("\t\tmix.model [%s]\n", model);
			printf("\t\tmix.integrator [%s]\n", integrator);
		}
		if (!strcasecmp(model, "GAUSSIAN")) {
			ds->mix_id = MIX_GAUSSIAN;
		} else if (!strcasecmp(model, "LOGGAMMA")) {
			ds->mix_id = MIX_LOGGAMMA;
		} else if (!strcasecmp(model, "MLOGGAMMA")) {
			ds->mix_id = MIX_MLOGGAMMA;
		} else {
			inla_error_field_is_void(__GMRFLib_FuncName, secname, "MIX.MODEL", model);
		}

		switch (ds->mix_id) {
		case MIX_GAUSSIAN:
		{
			/*
			 * get options related to the gaussian 
			 */
			tmp = iniparser_getdouble(ini, inla_string_join(secname, "MIX.INITIAL"), G.log_prec_initial);
			ds->mix_fixed = iniparser_getboolean(ini, inla_string_join(secname, "MIX.FIXED"), 0);
			if (!ds->mix_fixed && mb->mode_use_mode) {
				tmp = mb->theta_file[mb->theta_counter_file++];
				if (mb->mode_fixed) {
					ds->mix_fixed = 1;
				}
			}
			HYPER_NEW(ds->data_observations.mix_log_prec_gaussian, tmp);
			if (mb->verbose) {
				printf("\t\tinitialise mix.log_precision[%g]\n", ds->data_observations.mix_log_prec_gaussian[0][0]);
				printf("\t\tmix.fixed=[%1d]\n", ds->mix_fixed);
			}
			inla_read_prior_mix(mb, ini, sec, &(ds->mix_prior), "LOGGAMMA", NULL);

			/*
			 * add theta 
			 */
			if (!ds->mix_fixed) {
				mb->theta = Realloc(mb->theta, mb->ntheta + 1, double **);
				mb->theta_hyperid = Realloc(mb->theta_hyperid, mb->ntheta + 1, char *);
				mb->theta_hyperid[mb->ntheta] = ds->mix_prior.hyperid;
				mb->theta_tag = Realloc(mb->theta_tag, mb->ntheta + 1, char *);
				mb->theta_tag_userscale = Realloc(mb->theta_tag_userscale, mb->ntheta + 1, char *);
				mb->theta_dir = Realloc(mb->theta_dir, mb->ntheta + 1, char *);
				mb->theta_tag[mb->ntheta] = inla_make_tag("Log precision for the Gaussian mix", mb->ds);
				mb->theta_tag_userscale[mb->ntheta] = inla_make_tag("Precision for the Gaussian mix", mb->ds);
				GMRFLib_sprintf(&msg, "%s-parameter", secname);
				mb->theta_dir[mb->ntheta] = msg;

				mb->theta_from = Realloc(mb->theta_from, mb->ntheta + 1, char *);
				mb->theta_to = Realloc(mb->theta_to, mb->ntheta + 1, char *);
				mb->theta_from[mb->ntheta] = Strdup(ds->mix_prior.from_theta);
				mb->theta_to[mb->ntheta] = Strdup(ds->mix_prior.to_theta);

				mb->theta[mb->ntheta] = ds->data_observations.mix_log_prec_gaussian;
				mb->theta_map = Realloc(mb->theta_map, mb->ntheta + 1, map_func_tp *);
				mb->theta_map[mb->ntheta] = map_precision;
				mb->theta_map_arg = Realloc(mb->theta_map_arg, mb->ntheta + 1, void *);
				mb->theta_map_arg[mb->ntheta] = NULL;
				mb->ntheta++;
				ds->mix_ntheta++;
			}

			ds->mix_loglikelihood = ds->loglikelihood;
			ds->loglikelihood = loglikelihood_mix_gaussian;
		}
			break;

		case MIX_LOGGAMMA:
		case MIX_MLOGGAMMA:
		{
			/*
			 * get options related to the loggamma
			 */
			tmp = iniparser_getdouble(ini, inla_string_join(secname, "MIX.INITIAL"), 1.0);
			ds->mix_fixed = iniparser_getboolean(ini, inla_string_join(secname, "MIX.FIXED"), 0);
			if (!ds->mix_fixed && mb->mode_use_mode) {
				tmp = mb->theta_file[mb->theta_counter_file++];
				if (mb->mode_fixed) {
					ds->mix_fixed = 1;
				}
			}
			HYPER_NEW(ds->data_observations.mix_log_prec_loggamma, tmp);
			if (mb->verbose) {
				printf("\t\tinitialise mix.log_precision[%g]\n", ds->data_observations.mix_log_prec_loggamma[0][0]);
				printf("\t\tmix.fixed=[%1d]\n", ds->mix_fixed);
			}
			inla_read_prior_mix(mb, ini, sec, &(ds->mix_prior), "PCMGAMMA", NULL);

			/*
			 * add theta 
			 */
			if (!ds->mix_fixed) {
				mb->theta = Realloc(mb->theta, mb->ntheta + 1, double **);
				mb->theta_hyperid = Realloc(mb->theta_hyperid, mb->ntheta + 1, char *);
				mb->theta_hyperid[mb->ntheta] = ds->mix_prior.hyperid;
				mb->theta_tag = Realloc(mb->theta_tag, mb->ntheta + 1, char *);
				mb->theta_tag_userscale = Realloc(mb->theta_tag_userscale, mb->ntheta + 1, char *);
				mb->theta_dir = Realloc(mb->theta_dir, mb->ntheta + 1, char *);

				if (ds->mix_id == MIX_LOGGAMMA) {
					mb->theta_tag[mb->ntheta] = inla_make_tag("Log precision for the LogGamma mix", mb->ds);
					mb->theta_tag_userscale[mb->ntheta] = inla_make_tag("Precision for the LogGamma mix", mb->ds);
					GMRFLib_sprintf(&msg, "%s-parameter", secname);
					mb->theta_dir[mb->ntheta] = msg;
				} else {
					mb->theta_tag[mb->ntheta] = inla_make_tag("Log precision for the mLogGamma mix", mb->ds);
					mb->theta_tag_userscale[mb->ntheta] = inla_make_tag("Precision for the mLogGamma mix", mb->ds);
					GMRFLib_sprintf(&msg, "%s-parameter", secname);
					mb->theta_dir[mb->ntheta] = msg;
				}
				mb->theta_from = Realloc(mb->theta_from, mb->ntheta + 1, char *);
				mb->theta_to = Realloc(mb->theta_to, mb->ntheta + 1, char *);
				mb->theta_from[mb->ntheta] = Strdup(ds->mix_prior.from_theta);
				mb->theta_to[mb->ntheta] = Strdup(ds->mix_prior.to_theta);

				mb->theta[mb->ntheta] = ds->data_observations.mix_log_prec_loggamma;
				mb->theta_map = Realloc(mb->theta_map, mb->ntheta + 1, map_func_tp *);
				mb->theta_map[mb->ntheta] = map_precision;
				mb->theta_map_arg = Realloc(mb->theta_map_arg, mb->ntheta + 1, void *);
				mb->theta_map_arg[mb->ntheta] = NULL;
				mb->ntheta++;
				ds->mix_ntheta++;
			}

			ds->mix_loglikelihood = ds->loglikelihood;
			ds->loglikelihood = (ds->mix_id == MIX_LOGGAMMA ? loglikelihood_mix_loggamma : loglikelihood_mix_mloggamma);
		}
			break;

		default:
			GMRFLib_ASSERT(0 == 1, GMRFLib_ESNH);
			break;
		}
	}

	if ((ds->data_id != L_GAUSSIAN && ds->data_id != L_AGAUSSIAN && ds->data_id != L_STDGAUSSIAN && ds->data_id != L_GGAUSSIAN
	     && ds->data_id != L_SEM && ds->data_id != L_BC_GAUSSIAN) || ds->predictor_invlinkfunc != link_identity || ds->mix_use
	    || mb->expert_disable_gaussian_check) {
		GMRFLib_gaussian_data = GMRFLib_FALSE;
	}

	if (ds->data_id == L_FL) {
		// mark those indices belonging to 'FL'
		assert(mb->predictor_ndata >= 0);
		for (i = 0; i < mb->predictor_ndata; i++) {
			if (ds->data_observations.d[i]) {
				mb->fl[i] = 1;
			}
		}

		// replace the matrix with its transpose
		double **cc = NULL, **c = ds->data_observations.fl_c;
		int n = mb->predictor_ndata, m = L_FL_NC;
		cc = Calloc(n, double *);
		for (i = 0; i < n; i++) {
			cc[i] = Calloc(m, double);
		}
		for (j = 0; j < m; j++) {
			for (i = 0; i < n; i++) {
				cc[i][j] = c[j][i];
			}
		}
		ds->data_observations.fl_c = cc;
		for (j = 0; j < m; j++) {
			Free(c[j]);
		}
		Free(c);
	}

	return INLA_OK;
}

int inla_parse_ffield(inla_tp *mb, dictionary *ini, int sec)
{
#define _SET(a_, b_) mb->f_ ## a_[mb->nf] = b_
#define _OneOf(a_) (!strcasecmp(model, a_))
#define _OneOf2(a_, b_) (_OneOf(a_) || _OneOf(b_))
#define _OneOf3(a_, b_, c_) (_OneOf(a_) || _OneOf(b_) || _OneOf(c_))
#define _SetInitial(id_, val_) mb->f_initial[mb->nf][id_] = val_
	/*
	 * parse section = ffield 
	 */
	int i, j, k, jj, nlocations = 0, nc = 0, n = 0, zn = 0, zm = 0, s = 0, itmp, id = 0, bvalue = 0, fixed, order, slm_n = -1, slm_m = -1,
	    nstrata = 0, nsubject = 0, cgeneric_n = -1, cgeneric_debug = 0, nbeta = 0, ncov = 0;
	char *filename = NULL, *filenamec = NULL, *secname = NULL, *model = NULL, *ptmp = NULL, *ptmp2 = NULL, *msg =
	    NULL, default_tag[100], *file_loc = NULL, *ctmp = NULL, *rgeneric_filename = NULL, *rgeneric_model = NULL,
	    *cgeneric_shlib = NULL, *cgeneric_model = NULL;
	double **log_prec = NULL, **log_prec0 = NULL, **log_prec1 = NULL, **log_prec2, **phi_intern = NULL, **rho_intern =
	    NULL, **group_rho_intern = NULL, **group_prec_intern = NULL, **rho_intern01 = NULL, **rho_intern02 =
	    NULL, **rho_intern12 = NULL, **range_intern = NULL, tmp, **beta_intern = NULL, **beta = NULL, **h2_intern =
	    NULL, **a_intern = NULL, ***theta_iidwishart = NULL, **log_diag = NULL, rd, **mean_x = NULL, **log_prec_x =
	    NULL, ***pacf_intern = NULL, slm_rho_min = 0.0, slm_rho_max = 0.0, **log_halflife = NULL, **log_shape = NULL, **alpha =
	    NULL, **gama = NULL, **alpha1 = NULL, **alpha2 = NULL, **H_intern = NULL, **nu_intern = NULL, ***intslope_gamma = NULL, *cov = NULL,
	    *loc = NULL, ***betas = NULL;
	GMRFLib_matrix_tp *W = NULL;

	lt_dlhandle handle;
	inla_cgeneric_func_tp *model_func = NULL;
	inla_cgeneric_data_tp *cgeneric_data = NULL;

	GMRFLib_matrix_tp *intslope_def = NULL;
	GMRFLib_crwdef_tp *crwdef = NULL;
	inla_spde_tp *spde_model = NULL;
	inla_spde_tp *spde_model_orig = NULL;
	inla_spde2_tp *spde2_model = NULL;
	inla_spde2_tp *spde2_model_orig = NULL;
	inla_spde3_tp *spde3_model = NULL;
	inla_spde3_tp *spde3_model_orig = NULL;

	int thread_id = 0;
	assert(omp_get_thread_num() == 0);

	if (mb->verbose) {
		printf("\tinla_parse_ffield...\n");
	}
	secname = Strdup(iniparser_getsecname(ini, sec));
	if (mb->verbose) {
		printf("\t\tsection=[%s]\n", secname);
	}
	mb->f_tag = Realloc(mb->f_tag, mb->nf + 1, char *);
	mb->f_dir = Realloc(mb->f_dir, mb->nf + 1, char *);
	mb->f_modelname = Realloc(mb->f_modelname, mb->nf + 1, char *);
	mb->f_c = Realloc(mb->f_c, mb->nf + 1, int *);
	mb->f_n = Realloc(mb->f_n, mb->nf + 1, int);
	mb->f_N = Realloc(mb->f_N, mb->nf + 1, int);
	mb->f_Ntotal = Realloc(mb->f_Ntotal, mb->nf + 1, int);
	mb->f_order = Realloc(mb->f_order, mb->nf + 1, int);
	mb->f_nrep = Realloc(mb->f_nrep, mb->nf + 1, int);
	mb->f_ngroup = Realloc(mb->f_ngroup, mb->nf + 1, int);
	mb->f_group_model = Realloc(mb->f_group_model, mb->nf + 1, int);
	mb->f_group_cyclic = Realloc(mb->f_group_cyclic, mb->nf + 1, int);
	mb->f_group_order = Realloc(mb->f_group_order, mb->nf + 1, int);
	mb->f_group_graph = Realloc(mb->f_group_graph, mb->nf + 1, GMRFLib_graph_tp *);
	mb->f_nrow = Realloc(mb->f_nrow, mb->nf + 1, int);
	mb->f_ncol = Realloc(mb->f_ncol, mb->nf + 1, int);
	mb->f_locations = Realloc(mb->f_locations, mb->nf + 1, double *);
	mb->f_weights = Realloc(mb->f_weights, mb->nf + 1, double *);
	mb->f_scale = Realloc(mb->f_scale, mb->nf + 1, double *);
	mb->f_Qfunc = Realloc(mb->f_Qfunc, mb->nf + 1, GMRFLib_Qfunc_tp *);
	mb->f_Qfunc_orig = Realloc(mb->f_Qfunc_orig, mb->nf + 1, GMRFLib_Qfunc_tp *);
	mb->f_Qfunc_arg = Realloc(mb->f_Qfunc_arg, mb->nf + 1, void *);
	mb->f_Qfunc_arg_orig = Realloc(mb->f_Qfunc_arg_orig, mb->nf + 1, void *);
	mb->f_bfunc2 = Realloc(mb->f_bfunc2, mb->nf + 1, GMRFLib_bfunc2_tp *);
	mb->f_graph = Realloc(mb->f_graph, mb->nf + 1, GMRFLib_graph_tp *);
	mb->f_graph_orig = Realloc(mb->f_graph_orig, mb->nf + 1, GMRFLib_graph_tp *);
	mb->f_prior = Realloc(mb->f_prior, mb->nf + 1, Prior_tp *);
	mb->f_sumzero = Realloc(mb->f_sumzero, mb->nf + 1, char);
	mb->f_constr = Realloc(mb->f_constr, mb->nf + 1, GMRFLib_constr_tp *);
	mb->f_constr_orig = Realloc(mb->f_constr_orig, mb->nf + 1, GMRFLib_constr_tp *);
	mb->f_diag = Realloc(mb->f_diag, mb->nf + 1, double);
	mb->f_compute = Realloc(mb->f_compute, mb->nf + 1, int);
	mb->f_fixed = Realloc(mb->f_fixed, mb->nf + 1, int *);
	mb->f_initial = Realloc(mb->f_initial, mb->nf + 1, double *);
	mb->f_rankdef = Realloc(mb->f_rankdef, mb->nf + 1, double);
	mb->f_id = Realloc(mb->f_id, mb->nf + 1, inla_component_tp);
	mb->f_ntheta = Realloc(mb->f_ntheta, mb->nf + 1, int);
	mb->f_cyclic = Realloc(mb->f_cyclic, mb->nf + 1, int);
	mb->f_nu = Realloc(mb->f_nu, mb->nf + 1, int);
	mb->f_Torder = Realloc(mb->f_Torder, mb->nf + 1, int);
	mb->f_Tmodel = Realloc(mb->f_Tmodel, mb->nf + 1, char *);
	mb->f_Korder = Realloc(mb->f_Korder, mb->nf + 1, int);
	mb->f_Kmodel = Realloc(mb->f_Kmodel, mb->nf + 1, char *);
	mb->f_model = Realloc(mb->f_model, mb->nf + 1, void *);
	mb->f_theta = Realloc(mb->f_theta, mb->nf + 1, double ***);
	mb->f_theta_map = Realloc(mb->f_theta_map, mb->nf + 1, map_func_tp **);
	mb->f_theta_map_arg = Realloc(mb->f_theta_map_arg, mb->nf + 1, void **);
	mb->f_of = Realloc(mb->f_of, mb->nf + 1, char *);
	mb->f_same_as = Realloc(mb->f_same_as, mb->nf + 1, char *);
	mb->f_precision = Realloc(mb->f_precision, mb->nf + 1, double);
	mb->f_output = Realloc(mb->f_output, mb->nf + 1, Output_tp *);
	mb->f_id_names = Realloc(mb->f_id_names, mb->nf + 1, inla_file_contents_tp *);
	mb->f_correct = Realloc(mb->f_correct, mb->nf + 1, int);
	mb->f_vb_correct = Realloc(mb->f_vb_correct, mb->nf + 1, GMRFLib_idx_tp *);
	mb->f_Alocal = Realloc(mb->f_Alocal, mb->nf + 1, GMRFLib_matrix_tp *);

	/*
	 * set everything to `ZERO' initially 
	 */
	_SET(c, NULL);
	_SET(n, 0);
	_SET(N, 0);
	_SET(Ntotal, 0);
	_SET(order, -1);
	_SET(nrow, 0);
	_SET(ncol, 0);
	_SET(locations, NULL);
	_SET(weights, NULL);				       /* this means default weights = 1 */
	_SET(scale, NULL);
	_SET(Qfunc, (GMRFLib_Qfunc_tp *) NULL);
	_SET(Qfunc_orig, (GMRFLib_Qfunc_tp *) NULL);
	_SET(Qfunc_arg, NULL);
	_SET(Qfunc_arg_orig, NULL);
	_SET(bfunc2, NULL);
	_SET(graph, NULL);
	_SET(graph_orig, NULL);
	_SET(prior, NULL);
	_SET(sumzero, 0);
	_SET(constr, NULL);
	_SET(constr_orig, NULL);
	_SET(diag, 0.0);
	_SET(compute, 1);
	_SET(fixed, NULL);
	_SET(initial, NULL);
	_SET(rankdef, 0);
	_SET(output, NULL);
	_SET(id, INVALID_COMPONENT);
	_SET(ntheta, 0);
	_SET(theta, NULL);
	_SET(cyclic, 0);
	_SET(nu, -1);
	_SET(Tmodel, NULL);
	_SET(Torder, -1);
	_SET(Kmodel, NULL);
	_SET(Korder, -1);
	_SET(of, NULL);
	_SET(precision, exp(14.0));
	_SET(nrep, 1);
	_SET(ngroup, 1);
	_SET(group_model, G_EXCHANGEABLE);
	_SET(group_cyclic, 0);
	_SET(group_order, 0);
	_SET(group_graph, NULL);
	_SET(id_names, NULL);
	_SET(correct, -1);
	_SET(vb_correct, NULL);
	_SET(Alocal, NULL);

	sprintf(default_tag, "default tag for ffield %d", mb->nf);
	mb->f_tag[mb->nf] = (secname ? strdup(secname) : strdup(default_tag));
	mb->f_dir[mb->nf] = Strdup(iniparser_getstring(ini, inla_string_join(secname, "DIR"), Strdup(mb->f_tag[mb->nf])));
	if (mb->verbose) {
		printf("\t\tdir=[%s]\n", mb->f_dir[mb->nf]);
	}
	/*
	 * Somewhere, I might need to know this upfront... 
	 */
	HYPER_NEW(log_prec, 0.0);
	HYPER_NEW(log_prec0, 0.0);
	HYPER_NEW(log_prec1, 0.0);
	HYPER_NEW(log_prec2, 0.0);
	HYPER_NEW(H_intern, 0.0);
	HYPER_NEW(phi_intern, 0.0);
	HYPER_NEW(rho_intern, 0.0);
	HYPER_NEW(rho_intern01, 0.0);
	HYPER_NEW(rho_intern02, 0.0);
	HYPER_NEW(rho_intern12, 0.0);
	HYPER_NEW(range_intern, 0.0);
	HYPER_NEW(beta_intern, 0.0);
	HYPER_NEW(beta, 1.0);
	HYPER_NEW(group_rho_intern, 0.0);
	HYPER_NEW(group_prec_intern, 0.0);
	HYPER_NEW(h2_intern, 0.0);
	HYPER_NEW(a_intern, 0.0);
	HYPER_NEW(log_diag, 0.0);
	HYPER_NEW(mean_x, 0.0);
	HYPER_NEW(log_prec_x, 0.0);
	HYPER_NEW(log_halflife, 0.0);
	HYPER_NEW(log_shape, 0.0);
	HYPER_NEW(alpha, 0.0);
	HYPER_NEW(alpha1, 0.0);
	HYPER_NEW(alpha2, 0.0);
	HYPER_NEW(gama, 0.0);
	HYPER_NEW(nu_intern, 0.0);

	/*
	 * start parsing 
	 */
	model = Strdup(iniparser_getstring(ini, inla_string_join(secname, "MODEL"), NULL));
	if (mb->verbose) {
		printf("\t\tmodel=[%s]\n", model);
	}

	if (_OneOf("RW2D")) {
		mb->f_id[mb->nf] = F_RW2D;
		mb->f_ntheta[mb->nf] = 1;
		mb->f_modelname[mb->nf] = Strdup("Random walk 2D");
	} else if (_OneOf("RW2DIID")) {
		mb->f_id[mb->nf] = F_RW2DIID;
		mb->f_ntheta[mb->nf] = 2;
		mb->f_modelname[mb->nf] = Strdup("Random walk 2DIID");
	} else if (_OneOf("BESAG")) {
		mb->f_id[mb->nf] = F_BESAG;
		mb->f_ntheta[mb->nf] = 1;
		mb->f_modelname[mb->nf] = Strdup("Besags ICAR model");
	} else if (_OneOf("BESAG2")) {
		mb->f_id[mb->nf] = F_BESAG2;
		mb->f_ntheta[mb->nf] = 2;
		mb->f_modelname[mb->nf] = Strdup("Besags ICAR model for joint models");
	} else if (_OneOf("BYM")) {
		mb->f_id[mb->nf] = F_BYM;
		mb->f_ntheta[mb->nf] = 2;
		mb->f_modelname[mb->nf] = Strdup("BYM model");
	} else if (_OneOf("BYM2")) {
		mb->f_id[mb->nf] = F_BYM2;
		mb->f_ntheta[mb->nf] = 2;
		mb->f_modelname[mb->nf] = Strdup("BYM2 model");
	} else if (_OneOf("BESAGPROPER")) {
		mb->f_id[mb->nf] = F_BESAGPROPER;
		mb->f_ntheta[mb->nf] = 2;
		mb->f_modelname[mb->nf] = Strdup("Proper version of Besags ICAR model");
	} else if (_OneOf("BESAGPROPER2")) {
		mb->f_id[mb->nf] = F_BESAGPROPER2;
		mb->f_ntheta[mb->nf] = 2;
		mb->f_modelname[mb->nf] = Strdup("Proper version of Besags ICAR model version 2");
	} else if (_OneOf2("GENERIC", "GENERIC0")) {
		mb->f_id[mb->nf] = F_GENERIC0;
		mb->f_ntheta[mb->nf] = 1;
		mb->f_modelname[mb->nf] = Strdup("Generic0 model");
	} else if (_OneOf("GENERIC1")) {
		mb->f_id[mb->nf] = F_GENERIC1;
		mb->f_ntheta[mb->nf] = 2;
		mb->f_modelname[mb->nf] = Strdup("Generic1 model");
	} else if (_OneOf("GENERIC2")) {
		mb->f_id[mb->nf] = F_GENERIC2;
		mb->f_ntheta[mb->nf] = 2;
		mb->f_modelname[mb->nf] = Strdup("Generic2 model");
	} else if (_OneOf("GENERIC3")) {
		mb->f_id[mb->nf] = F_GENERIC3;
		mb->f_ntheta[mb->nf] = GENERIC3_MAXTHETA;
		mb->f_modelname[mb->nf] = Strdup("Generic3 model");
	} else if (_OneOf("SEASONAL")) {
		mb->f_id[mb->nf] = F_SEASONAL;
		mb->f_ntheta[mb->nf] = 1;
		mb->f_modelname[mb->nf] = Strdup("Seasonal model");
	} else if (_OneOf("IID")) {
		mb->f_id[mb->nf] = F_IID;
		mb->f_ntheta[mb->nf] = 1;
		mb->f_modelname[mb->nf] = Strdup("IID model");
	} else if (_OneOf("IID1D")) {
		mb->f_id[mb->nf] = F_IID1D;
		mb->f_ntheta[mb->nf] = inla_iid_wishart_nparam(WISHART_DIM(mb->nf));
		mb->f_modelname[mb->nf] = Strdup("IID1D model");
	} else if (_OneOf("IID2D")) {
		mb->f_id[mb->nf] = F_IID2D;
		mb->f_ntheta[mb->nf] = inla_iid_wishart_nparam(WISHART_DIM(mb->nf));
		mb->f_modelname[mb->nf] = Strdup("IID2D model");
	} else if (_OneOf("IID3D")) {
		mb->f_id[mb->nf] = F_IID3D;
		mb->f_ntheta[mb->nf] = inla_iid_wishart_nparam(WISHART_DIM(mb->nf));
		mb->f_modelname[mb->nf] = Strdup("IID3D model");
	} else if (_OneOf("IID4D")) {
		mb->f_id[mb->nf] = F_IID4D;
		mb->f_ntheta[mb->nf] = inla_iid_wishart_nparam(WISHART_DIM(mb->nf));
		mb->f_modelname[mb->nf] = Strdup("IID4D model");
	} else if (_OneOf("IID5D")) {
		mb->f_id[mb->nf] = F_IID5D;
		mb->f_ntheta[mb->nf] = inla_iid_wishart_nparam(WISHART_DIM(mb->nf));
		mb->f_modelname[mb->nf] = Strdup("IID5D model");
	} else if (_OneOf("IIDKD")) {
		mb->f_id[mb->nf] = F_IIDKD;
		mb->f_order[mb->nf] = iniparser_getint(ini, inla_string_join(secname, "ORDER"), 0);
		mb->f_ntheta[mb->nf] = INLA_WISHARTK_NTHETA(mb->f_order[mb->nf]);
		mb->f_modelname[mb->nf] = Strdup("IIDKD model");
	} else if (_OneOf("2DIID")) {
		mb->f_id[mb->nf] = F_2DIID;
		mb->f_ntheta[mb->nf] = 3;
		mb->f_modelname[mb->nf] = Strdup("2DIID model");
	} else if (_OneOf("INTSLOPE")) {
		mb->f_id[mb->nf] = F_INTSLOPE;
		mb->f_ntheta[mb->nf] = inla_iid_wishart_nparam(2L) + INTSLOPE_MAXTHETA;
		mb->f_modelname[mb->nf] = Strdup("INTSLOPE model");
	} else if (_OneOf("MEC")) {
		mb->f_id[mb->nf] = F_MEC;
		mb->f_ntheta[mb->nf] = 4;
		mb->f_modelname[mb->nf] = Strdup("MEC");
	} else if (_OneOf("MEB")) {
		mb->f_id[mb->nf] = F_MEB;
		mb->f_ntheta[mb->nf] = 2;
		mb->f_modelname[mb->nf] = Strdup("MEB");
	} else if (_OneOf("RGENERIC")) {
		mb->f_id[mb->nf] = F_R_GENERIC;
		mb->f_ntheta[mb->nf] = -1;
		mb->f_modelname[mb->nf] = Strdup("RGeneric2");
	} else if (_OneOf("CGENERIC")) {
		mb->f_id[mb->nf] = F_C_GENERIC;
		mb->f_ntheta[mb->nf] = -1;
		mb->f_modelname[mb->nf] = Strdup("CGeneric");
	} else if (_OneOf("RW1")) {
		mb->f_id[mb->nf] = F_RW1;
		mb->f_ntheta[mb->nf] = 1;
		mb->f_modelname[mb->nf] = Strdup("RW1 model");
	} else if (_OneOf("RW2")) {
		mb->f_id[mb->nf] = F_RW2;
		mb->f_ntheta[mb->nf] = 1;
		mb->f_modelname[mb->nf] = Strdup("RW2 model");
	} else if (_OneOf("CRW2")) {
		mb->f_id[mb->nf] = F_CRW2;
		mb->f_ntheta[mb->nf] = 1;
		mb->f_modelname[mb->nf] = Strdup("CRW2 model");
	} else if (_OneOf("AR1")) {
		mb->f_id[mb->nf] = F_AR1;
		mb->f_ntheta[mb->nf] = 3;
		mb->f_modelname[mb->nf] = Strdup("AR1 model");
	} else if (_OneOf("AR1C")) {
		mb->f_id[mb->nf] = F_AR1C;
		mb->f_ntheta[mb->nf] = 2;
		mb->f_modelname[mb->nf] = Strdup("AR1C model");
	} else if (_OneOf("FGN")) {
		mb->f_id[mb->nf] = F_FGN;
		mb->f_ntheta[mb->nf] = 2;
		mb->f_modelname[mb->nf] = Strdup("FGN model");
	} else if (_OneOf("FGN2")) {
		mb->f_id[mb->nf] = F_FGN2;
		mb->f_ntheta[mb->nf] = 2;
		mb->f_modelname[mb->nf] = Strdup("FGN2 model");
	} else if (_OneOf("OU")) {
		mb->f_id[mb->nf] = F_OU;
		mb->f_ntheta[mb->nf] = 2;
		mb->f_modelname[mb->nf] = Strdup("Ornstein-Uhlenbeck model");
	} else if (_OneOf("MATERN2D")) {
		mb->f_id[mb->nf] = F_MATERN2D;
		mb->f_ntheta[mb->nf] = 2;
		mb->f_modelname[mb->nf] = Strdup("Matern2D model");
	} else if (_OneOf("DMATERN")) {
		mb->f_id[mb->nf] = F_DMATERN;
		mb->f_ntheta[mb->nf] = 3;
		mb->f_modelname[mb->nf] = Strdup("DMatern model");
	} else if (_OneOf("Z")) {
		mb->f_id[mb->nf] = F_Z;
		mb->f_ntheta[mb->nf] = 1;
		mb->f_modelname[mb->nf] = Strdup("Z model");
	} else if (_OneOf("SLM")) {
		mb->f_id[mb->nf] = F_SLM;
		mb->f_ntheta[mb->nf] = 2;
		mb->f_modelname[mb->nf] = Strdup("SLM model");
	} else if (_OneOf("SPDE")) {
		mb->f_id[mb->nf] = F_SPDE;
		mb->f_ntheta[mb->nf] = 4;
		mb->f_modelname[mb->nf] = Strdup("SPDE model");
	} else if (_OneOf("SPDE2")) {
		mb->f_id[mb->nf] = F_SPDE2;
		mb->f_ntheta[mb->nf] = -1;		       /* Not known yet */
		mb->f_modelname[mb->nf] = Strdup("SPDE2 model");
	} else if (_OneOf("SPDE3")) {
		mb->f_id[mb->nf] = F_SPDE3;
		mb->f_ntheta[mb->nf] = -1;		       /* Not known yet */
		mb->f_modelname[mb->nf] = Strdup("SPDE3 model");
	} else if (_OneOf("AR")) {
		mb->f_id[mb->nf] = F_AR;
		mb->f_ntheta[mb->nf] = AR_MAXTHETA + 1;
		mb->f_modelname[mb->nf] = Strdup("AR(p) model");
	} else if (_OneOf("COPY")) {
		mb->f_id[mb->nf] = F_COPY;
		mb->f_ntheta[mb->nf] = 1;
		mb->f_modelname[mb->nf] = Strdup("Copy");
	} else if (_OneOf("SCOPY")) {
		mb->f_id[mb->nf] = F_SCOPY;
		mb->f_ntheta[mb->nf] = SCOPY_MAXTHETA;
		mb->f_modelname[mb->nf] = Strdup("SCopy");
	} else if (_OneOf("CLINEAR")) {
		mb->f_id[mb->nf] = F_CLINEAR;
		mb->f_ntheta[mb->nf] = 1;
		mb->f_modelname[mb->nf] = Strdup("Constrained linear");
	} else if (_OneOf("LOG1EXP")) {
		mb->f_id[mb->nf] = F_LOG1EXP;
		mb->f_ntheta[mb->nf] = 3;
		mb->f_modelname[mb->nf] = Strdup("Log1Exp");
	} else if (_OneOf("LOGDIST")) {
		mb->f_id[mb->nf] = F_LOGDIST;
		mb->f_ntheta[mb->nf] = 3;
		mb->f_modelname[mb->nf] = Strdup("LogDist");
	} else if (_OneOf("SIGM")) {
		mb->f_id[mb->nf] = F_SIGM;
		mb->f_ntheta[mb->nf] = 3;
		mb->f_modelname[mb->nf] = Strdup("Sigmodial");
	} else if (_OneOf("REVSIGM")) {
		mb->f_id[mb->nf] = F_REVSIGM;
		mb->f_ntheta[mb->nf] = 3;
		mb->f_modelname[mb->nf] = Strdup("Reverse sigmodial");
	} else {
		inla_error_field_is_void(__GMRFLib_FuncName, secname, "model", model);
	}
	if (mb->f_ntheta[mb->nf] > 0) {
		mb->f_prior[mb->nf] = Calloc(mb->f_ntheta[mb->nf], Prior_tp);
		mb->f_fixed[mb->nf] = Calloc(mb->f_ntheta[mb->nf], int);
		mb->f_theta[mb->nf] = Calloc(mb->f_ntheta[mb->nf], double **);
	} else {
		mb->f_prior[mb->nf] = NULL;
		mb->f_fixed[mb->nf] = NULL;
		mb->f_theta[mb->nf] = NULL;
	}

	/*
	 * just allocate this here, as its needed all over 
	 */
	if (mb->f_ntheta[mb->nf] > 0) {
		mb->f_initial[mb->nf] = Calloc(mb->f_ntheta[mb->nf], double);
	}

	id = mb->f_id[mb->nf];				       /* shortcut */

	switch (id) {
	case F_RW2D:
	case F_BESAG:
	case F_GENERIC0:
	case F_SEASONAL:
	case F_IID:
	case F_RW1:
	case F_RW2:
	case F_CRW2:
	case F_Z:
	{
		inla_read_prior(mb, ini, sec, &(mb->f_prior[mb->nf][0]), "LOGGAMMA", NULL);
	}
		break;

	case F_SLM:
	{
		inla_read_prior0(mb, ini, sec, &(mb->f_prior[mb->nf][0]), "LOGGAMMA", NULL);	// kappa
		inla_read_prior1(mb, ini, sec, &(mb->f_prior[mb->nf][1]), "NORMAL", NULL);	// rho
	}
		break;

	case F_BESAG2:
	{
		inla_read_prior0(mb, ini, sec, &(mb->f_prior[mb->nf][0]), "LOGGAMMA", NULL);	// kappa
		inla_read_prior1(mb, ini, sec, &(mb->f_prior[mb->nf][1]), "NORMAL-a", NULL);	// a
	}
		break;

	case F_BESAGPROPER:
	{
		inla_read_prior0(mb, ini, sec, &(mb->f_prior[mb->nf][0]), "LOGGAMMA", NULL);	// precision
		inla_read_prior1(mb, ini, sec, &(mb->f_prior[mb->nf][1]), "LOGGAMMA", NULL);	// weight
	}
		break;

	case F_BESAGPROPER2:
	{
		inla_read_prior0(mb, ini, sec, &(mb->f_prior[mb->nf][0]), "LOGGAMMA", NULL);	// precision
		inla_read_prior1(mb, ini, sec, &(mb->f_prior[mb->nf][1]), "GAUSSIAN", NULL);	// lambda
	}
		break;

	case F_SPDE:
	{
		inla_read_prior0(mb, ini, sec, &(mb->f_prior[mb->nf][0]), "NORMAL", NULL);	// T[0]
		inla_read_prior1(mb, ini, sec, &(mb->f_prior[mb->nf][1]), "NORMAL", NULL);	// K[0]
		inla_read_prior2(mb, ini, sec, &(mb->f_prior[mb->nf][2]), "NORMAL", NULL);	// the rest
		inla_read_prior3(mb, ini, sec, &(mb->f_prior[mb->nf][3]), "FLAT", NULL);	// the ocillating cooef
	}
		break;

	case F_SPDE2:
	{
		mb->f_prior[mb->nf] = Calloc(1, Prior_tp);
		inla_read_prior0(mb, ini, sec, &(mb->f_prior[mb->nf][0]), "MVNORM", NULL);	// Just one prior...
	}
		break;

	case F_SPDE3:
	{
		mb->f_prior[mb->nf] = Calloc(1, Prior_tp);
		inla_read_prior0(mb, ini, sec, &(mb->f_prior[mb->nf][0]), "MVNORM", NULL);	// Just one prior...
	}
		break;

	case F_AR:
	{
		assert(11 == AR_MAXTHETA + 1);
		mb->f_prior[mb->nf] = Calloc(11, Prior_tp);
		inla_read_prior0(mb, ini, sec, &(mb->f_prior[mb->nf][0]), "PCPREC", NULL);	// log precision
		inla_read_prior1(mb, ini, sec, &(mb->f_prior[mb->nf][1]), "PCRHO0", NULL);	// the pacf
		inla_read_prior2(mb, ini, sec, &(mb->f_prior[mb->nf][2]), "PCRHO0", NULL);	// the pacf
		inla_read_prior3(mb, ini, sec, &(mb->f_prior[mb->nf][3]), "PCRHO0", NULL);	// the pacf
		inla_read_prior4(mb, ini, sec, &(mb->f_prior[mb->nf][4]), "PCRHO0", NULL);	// the pacf
		inla_read_prior5(mb, ini, sec, &(mb->f_prior[mb->nf][5]), "PCRHO0", NULL);	// the pacf
		inla_read_prior6(mb, ini, sec, &(mb->f_prior[mb->nf][6]), "PCRHO0", NULL);	// the pacf
		inla_read_prior7(mb, ini, sec, &(mb->f_prior[mb->nf][7]), "PCRHO0", NULL);	// the pacf
		inla_read_prior8(mb, ini, sec, &(mb->f_prior[mb->nf][8]), "PCRHO0", NULL);	// the pacf
		inla_read_prior9(mb, ini, sec, &(mb->f_prior[mb->nf][9]), "PCRHO0", NULL);	// the pacf
		inla_read_prior10(mb, ini, sec, &(mb->f_prior[mb->nf][10]), "PCRHO0", NULL);	// the pacf
	}
		break;

	case F_COPY:
	{
		inla_read_prior(mb, ini, sec, &(mb->f_prior[mb->nf][0]), "NORMAL-1", NULL);
	}
		break;

	case F_SCOPY:
	{
		// we do not really use these, as we define priors differently, and we do that in
		// control.scopy
		mb->f_prior[mb->nf] = Calloc(SCOPY_MAXTHETA, Prior_tp);
		for (k = 0; k < SCOPY_MAXTHETA; k++) {
			inla_read_priorN(mb, ini, sec, &(mb->f_prior[mb->nf][k]), "NONE", k, NULL);
		}
	}
		break;

	case F_CLINEAR:
	{
		inla_read_prior(mb, ini, sec, &(mb->f_prior[mb->nf][0]), "NORMAL", NULL);
	}
		break;

	case F_LOG1EXP:
	case F_LOGDIST:
	case F_SIGM:
	case F_REVSIGM:
	{
		inla_read_prior0(mb, ini, sec, &(mb->f_prior[mb->nf][0]), "NORMAL", NULL);	// beta
		inla_read_prior1(mb, ini, sec, &(mb->f_prior[mb->nf][1]), "NORMAL", NULL);	// log_halflife
		inla_read_prior2(mb, ini, sec, &(mb->f_prior[mb->nf][2]), "NORMAL", NULL);	// log_shape
	}
		break;

	case F_IID1D:
	case F_IID2D:
	case F_IID3D:
	case F_IID4D:
	case F_IID5D:
	{
		int dim = WISHART_DIM(mb->nf);
		assert(dim > 0);
		char *pri = NULL, *par = NULL, *to_theta = NULL, *from_theta = NULL, *prifunc = NULL, *hyperid = NULL;
		int nt = inla_iid_wishart_nparam(dim);

		GMRFLib_sprintf(&prifunc, "WISHART%1dD", dim);
		int kk;
		if (dim > 1) {
			for (kk = 0; kk < nt; kk++) {
				GMRFLib_sprintf(&pri, "PRIOR%1d", kk);
				GMRFLib_sprintf(&par, "PARAMETERS%1d", kk);
				GMRFLib_sprintf(&from_theta, "FROM.THETA%1d", kk);
				GMRFLib_sprintf(&to_theta, "TO.THETA%1d", kk);
				GMRFLib_sprintf(&hyperid, "HYPERID%1d", kk);
				inla_read_prior_generic(mb, ini, sec, &(mb->f_prior[mb->nf][kk]), pri, par, from_theta, to_theta,
							hyperid, (kk == 0 ? prifunc : "NONE"), NULL);

				Free(pri);
				Free(par);
				Free(from_theta);
				Free(to_theta);
				Free(hyperid);
			}
		} else {
			kk = 0;
			GMRFLib_sprintf(&pri, "PRIOR");
			GMRFLib_sprintf(&par, "PARAMETERS");
			GMRFLib_sprintf(&from_theta, "FROM.THETA");
			GMRFLib_sprintf(&to_theta, "TO.THETA");
			GMRFLib_sprintf(&hyperid, "HYPERID");
			inla_read_prior_generic(mb, ini, sec, &(mb->f_prior[mb->nf][kk]), pri, par, from_theta, to_theta, hyperid,
						(kk == 0 ? prifunc : "NONE"), NULL);
		}

		Free(pri);
		Free(par);
	}
		break;

	case F_IIDKD:
	{
		int dim = mb->f_order[mb->nf];
		assert(dim >= INLA_WISHARTK_KMIN);
		assert(dim <= INLA_WISHARTK_KMAX);
		char *pri = NULL, *par = NULL, *to_theta = NULL, *from_theta = NULL, *prifunc = NULL, *hyperid = NULL;
		int nt_max = INLA_WISHARTK_NTHETA(INLA_WISHARTK_KMAX);
		int nt = INLA_WISHARTK_NTHETA(dim);

		GMRFLib_sprintf(&prifunc, "WISHARTK%1dD", dim);
		int kk;
		for (kk = 0; kk < nt_max; kk++) {
			GMRFLib_sprintf(&pri, "PRIOR%1d", kk);
			GMRFLib_sprintf(&par, "PARAMETERS%1d", kk);
			GMRFLib_sprintf(&from_theta, "FROM.THETA%1d", kk);
			GMRFLib_sprintf(&to_theta, "TO.THETA%1d", kk);
			GMRFLib_sprintf(&hyperid, "HYPERID%1d", kk);
			if (kk < nt) {
				inla_read_prior_generic(mb, ini, sec, &(mb->f_prior[mb->nf][kk]), pri, par, from_theta, to_theta,
							hyperid, (kk == 0 ? prifunc : "NONE"), NULL);
			} else {
				// mark them as read
				iniparser_getstring(ini, inla_string_join(secname, pri), NULL);
				iniparser_getstring(ini, inla_string_join(secname, par), NULL);
				iniparser_getstring(ini, inla_string_join(secname, from_theta), NULL);
				iniparser_getstring(ini, inla_string_join(secname, to_theta), NULL);
				iniparser_getstring(ini, inla_string_join(secname, hyperid), NULL);
			}

			Free(pri);
			Free(par);
			Free(from_theta);
			Free(to_theta);
			Free(hyperid);
		}
	}
		break;

	case F_INTSLOPE:
	{
		char *pri = NULL, *par = NULL, *to_theta = NULL, *from_theta = NULL, *prifunc = NULL, *hyperid = NULL;
		int dim = 2;
		GMRFLib_sprintf(&prifunc, "WISHART%1dD", dim);
		for (int kk = 0; kk < mb->f_ntheta[mb->nf]; kk++) {
			GMRFLib_sprintf(&pri, "PRIOR%1d", kk);
			GMRFLib_sprintf(&par, "PARAMETERS%1d", kk);
			GMRFLib_sprintf(&from_theta, "FROM.THETA%1d", kk);
			GMRFLib_sprintf(&to_theta, "TO.THETA%1d", kk);
			GMRFLib_sprintf(&hyperid, "HYPERID%1d", kk);
			inla_read_prior_generic(mb, ini, sec, &(mb->f_prior[mb->nf][kk]), pri, par, from_theta, to_theta,
						hyperid, (kk == 0 ? prifunc : "NONE"), NULL);

			Free(pri);
			Free(par);
			Free(from_theta);
			Free(to_theta);
			Free(hyperid);
		}
	}
		break;

	case F_BYM:
	{
		inla_read_prior0(mb, ini, sec, &(mb->f_prior[mb->nf][0]), "LOGGAMMA", NULL);	/* precision0 iid */
		inla_read_prior1(mb, ini, sec, &(mb->f_prior[mb->nf][1]), "LOGGAMMA", NULL);	/* precision1 spatial */
	}
		break;

	case F_BYM2:
	{
		inla_read_prior0(mb, ini, sec, &(mb->f_prior[mb->nf][0]), "LOGGAMMA", NULL);	/* precision */
		inla_read_prior1(mb, ini, sec, &(mb->f_prior[mb->nf][1]), "GAUSSIAN", NULL);	/* phi */
	}
		break;

	case F_RW2DIID:
	{
		inla_read_prior0(mb, ini, sec, &(mb->f_prior[mb->nf][0]), "LOGGAMMA", NULL);	/* precision */
		inla_read_prior1(mb, ini, sec, &(mb->f_prior[mb->nf][1]), "GAUSSIAN", NULL);	/* phi */
	}
		break;

	case F_AR1:
	{
		inla_read_prior0(mb, ini, sec, &(mb->f_prior[mb->nf][0]), "LOGGAMMA", NULL);	/* marginal precision */
		inla_read_prior1(mb, ini, sec, &(mb->f_prior[mb->nf][1]), "GAUSSIAN-rho", NULL);	/* phi (lag-1 correlation) */
		inla_read_prior2(mb, ini, sec, &(mb->f_prior[mb->nf][2]), "GAUSSIAN", NULL);	/* mean */
	}
		break;

	case F_AR1C:
	{
		inla_read_prior0(mb, ini, sec, &(mb->f_prior[mb->nf][0]), "LOGGAMMA", NULL);	/* marginal precision */
		inla_read_prior1(mb, ini, sec, &(mb->f_prior[mb->nf][1]), "GAUSSIAN-rho", NULL);	/* phi (lag-1 correlation) */
	}
		break;

	case F_FGN:
	case F_FGN2:
	{
		inla_read_prior0(mb, ini, sec, &(mb->f_prior[mb->nf][0]), "LOGGAMMA", NULL);	/* marginal precision */
		inla_read_prior1(mb, ini, sec, &(mb->f_prior[mb->nf][1]), "GAUSSIAN", NULL);	/* H */
	}
		break;

	case F_OU:
	{
		inla_read_prior0(mb, ini, sec, &(mb->f_prior[mb->nf][0]), "LOGGAMMA", NULL);	/* marginal precision */
		inla_read_prior1(mb, ini, sec, &(mb->f_prior[mb->nf][1]), "GAUSSIAN", NULL);	/* log(phi) */
	}
		break;

	case F_GENERIC1:
	{
		inla_read_prior0(mb, ini, sec, &(mb->f_prior[mb->nf][0]), "LOGGAMMA", NULL);	/* precision */
		inla_read_prior1(mb, ini, sec, &(mb->f_prior[mb->nf][1]), "GAUSSIAN", NULL);	/* beta */
	}
		break;

	case F_GENERIC2:
	{
		inla_read_prior0(mb, ini, sec, &(mb->f_prior[mb->nf][0]), "LOGGAMMA", NULL);	/* precision Cmatrix */
		inla_read_prior1(mb, ini, sec, &(mb->f_prior[mb->nf][1]), "LOGGAMMA", NULL);	/* the other precision, but theta1 = h^2 */
	}
		break;

	case F_GENERIC3:
	{
		for (i = 0; i < GENERIC3_MAXTHETA; i++) {
			inla_read_priorN(mb, ini, sec, &(mb->f_prior[mb->nf][i]), "LOGGAMMA", i, NULL);
		}
	}
		break;

	case F_2DIID:
	{
		inla_read_prior0(mb, ini, sec, &(mb->f_prior[mb->nf][0]), "LOGGAMMA", NULL);	/* precision0 */
		inla_read_prior1(mb, ini, sec, &(mb->f_prior[mb->nf][1]), "LOGGAMMA", NULL);	/* precision1 */
		inla_read_prior2(mb, ini, sec, &(mb->f_prior[mb->nf][2]), "GAUSSIAN-rho", NULL);	/* rho */
	}
		break;

	case F_MATERN2D:
	{
		inla_read_prior0(mb, ini, sec, &(mb->f_prior[mb->nf][0]), "LOGGAMMA", NULL);	/* precision */
		inla_read_prior1(mb, ini, sec, &(mb->f_prior[mb->nf][1]), "LOGGAMMA", NULL);	/* range */
	}
		break;

	case F_DMATERN:
	{
		inla_read_prior0(mb, ini, sec, &(mb->f_prior[mb->nf][0]), "LOGGAMMA", NULL);	/* precision */
		inla_read_prior1(mb, ini, sec, &(mb->f_prior[mb->nf][1]), "LOGGAMMA", NULL);	/* range */
		inla_read_prior2(mb, ini, sec, &(mb->f_prior[mb->nf][2]), "LOGGAMMA", NULL);	/* nu */
	}
		break;

	case F_MEC:
	{
		inla_read_prior0(mb, ini, sec, &(mb->f_prior[mb->nf][0]), "NORMAL", NULL);	/* beta */
		inla_read_prior1(mb, ini, sec, &(mb->f_prior[mb->nf][1]), "LOGGAMMA", NULL);	/* prec.u */
		inla_read_prior2(mb, ini, sec, &(mb->f_prior[mb->nf][2]), "NORMAL", NULL);	/* mean */
		inla_read_prior3(mb, ini, sec, &(mb->f_prior[mb->nf][3]), "LOGGAMMA", NULL);	/* prec.x */
	}
		break;

	case F_MEB:
	{
		inla_read_prior0(mb, ini, sec, &(mb->f_prior[mb->nf][0]), "NORMAL", NULL);	/* beta */
		inla_read_prior1(mb, ini, sec, &(mb->f_prior[mb->nf][1]), "LOGGAMMA", NULL);	/* prec.u */
	}
		break;

	case F_R_GENERIC:
	case F_C_GENERIC:
		break;

	default:
		abort();
	}

	char *str = iniparser_getstring(ini, inla_string_join(secname, "VB.CORRECT"), NULL);
	int *idxs = NULL;
	int n_idxs = 0;

	inla_sread_ints_q(&idxs, &n_idxs, str);
	if (n_idxs == 0) {
		GMRFLib_idx_add(&(mb->f_vb_correct[mb->nf]), -1);
	} else {
		for (i = 0; i < n_idxs; i++) {
			GMRFLib_idx_add(&(mb->f_vb_correct[mb->nf]), idxs[i]);
		}
	}
	if (mb->verbose) {
		printf("\t\tvb.correct n[%1d] ", mb->f_vb_correct[mb->nf]->n);
		for (i = 0; i < mb->f_vb_correct[mb->nf]->n; i++) {
			printf(" %1d", mb->f_vb_correct[mb->nf]->idx[i]);
		}
		printf("\n");
	}
	mb->f_correct[mb->nf] = iniparser_getint(ini, inla_string_join(secname, "CORRECT"), -1);
	if (mb->verbose) {
		printf("\t\tcorrect=[%1d]\n", mb->f_correct[mb->nf]);
	}
	mb->f_sumzero[mb->nf] = (char) iniparser_getboolean(ini, inla_string_join(secname, "CONSTRAINT"), 0);
	if (mb->verbose) {
		printf("\t\tconstr=[%1d]\n", mb->f_sumzero[mb->nf]);
	}
	mb->f_diag[mb->nf] = iniparser_getdouble(ini, inla_string_join(secname, "DIAGONAL"), 0.0);
	if (mb->verbose) {
		printf("\t\tdiagonal=[%g]\n", mb->f_diag[mb->nf]);
	}
	mb->f_id_names[mb->nf] = inla_read_file_contents(Strdup(iniparser_getstring(ini, inla_string_join(secname, "ID.NAMES"), NULL)));
	if (mb->verbose) {
		printf("\t\tid.names=%s\n", (mb->f_id_names[mb->nf] ? "<read>" : "<not present>"));
	}

	mb->f_compute[mb->nf] = iniparser_getboolean(ini, inla_string_join(secname, "COMPUTE"), 1);
	if (G.mode == INLA_MODE_HYPER) {
		if (mb->f_compute[mb->nf]) {
			fprintf(stderr, "*** Warning: HYPER_MODE require f_compute[%1d] = 0\n", mb->nf);
		}
		mb->f_compute[mb->nf] = 0;
	}
	if (mb->verbose) {
		printf("\t\tcompute=[%1d]\n", mb->f_compute[mb->nf]);
	}
	if (mb->f_ntheta[mb->nf] == 1) {
		mb->f_fixed[mb->nf][0] = iniparser_getboolean(ini, inla_string_join(secname, "FIXED"), 0);
	} else {
		/*
		 * reads FIXED0, FIXED1, FIXED2, etc.... 
		 */
		for (i = 0; i < mb->f_ntheta[mb->nf]; i++) {
			char *fixname = NULL;

			GMRFLib_sprintf(&fixname, "FIXED%1d", i);
			mb->f_fixed[mb->nf][i] = iniparser_getboolean(ini, inla_string_join(secname, fixname), 0);
			Free(fixname);
		}
	}

	mb->f_nrep[mb->nf] = iniparser_getint(ini, inla_string_join(secname, "NREP"), mb->f_nrep[mb->nf]);
	mb->f_ngroup[mb->nf] = iniparser_getint(ini, inla_string_join(secname, "NGROUP"), mb->f_ngroup[mb->nf]);
	if (mb->verbose) {
		printf("\t\tnrep=[%1d]\n", mb->f_nrep[mb->nf]);
		printf("\t\tngroup=[%1d]\n", mb->f_ngroup[mb->nf]);
	}

	char *fnm_Alocal = iniparser_getstring(ini, inla_string_join(secname, "A.local"), NULL);
	if (fnm_Alocal) {
		assert(GMRFLib_is_fmesher_file(fnm_Alocal, (long int) 0, -1) == GMRFLib_SUCCESS);
		mb->f_Alocal[mb->nf] = GMRFLib_read_fmesher_file(fnm_Alocal, (long int) 0, -1);
	}
	if (mb->verbose) {
		printf("\t\tAlocal=[%s]\n", (mb->f_Alocal[mb->nf] ? "yes" : "no"));
		if (mb->f_Alocal[mb->nf]) {
			printf("\t\t\tnrow x ncol =[%1d x %1d]\n", mb->f_Alocal[mb->nf]->nrow, mb->f_Alocal[mb->nf]->ncol);
			printf("\t\t\ttype =[%s]\n", (mb->f_Alocal[mb->nf]->values ? "sparse" : "dense"));
			if (mb->f_Alocal[mb->nf]->values) {
				printf("\t\t\telems = [%1d]\n", mb->f_Alocal[mb->nf]->elems);
			}
		}
	}

	/*
	 * to avoid errors from the R-interface. This just mark them as 'read'. For some reason file_loc is fine here, but later on it just return "", FIXME! 
	 */
	ptmp = Strdup(iniparser_getstring(ini, inla_string_join(secname, "N"), NULL));
	file_loc = Strdup(iniparser_getstring(ini, inla_string_join(secname, "LOCATIONS"), NULL));

	if ((mb->f_id[mb->nf] == F_MATERN2D) || (mb->f_id[mb->nf] == F_RW2D) || (mb->f_id[mb->nf] == F_RW2DIID)) {
		/*
		 * this case is a bit special....
		 *
		 * if predictor->nrow exists this is the default and must be equal if NROW exists. same with NCOL.
		 */
		char *tmp_nrow = Strdup(iniparser_getstring(ini, inla_string_join(secname, "NROW"), NULL));
		char *tmp_ncol = Strdup(iniparser_getstring(ini, inla_string_join(secname, "NCOL"), NULL));

		if (!tmp_nrow) {
			inla_error_missing_required_field(__GMRFLib_FuncName, secname, "nrow");
		} else {
			if (inla_sread_ints(&(mb->f_nrow[mb->nf]), 1, tmp_nrow) == INLA_FAIL) {
				inla_error_field_is_void(__GMRFLib_FuncName, secname, "nrow", tmp_nrow);
			}
			if (mb->f_nrow[mb->nf] <= 0) {
				inla_error_field_is_void(__GMRFLib_FuncName, secname, "nrow", tmp_nrow);
			}
		}
		if (!tmp_ncol) {
			inla_error_missing_required_field(__GMRFLib_FuncName, secname, "ncol");
		} else {
			if (inla_sread_ints(&(mb->f_ncol[mb->nf]), 1, tmp_ncol) == INLA_FAIL) {
				inla_error_field_is_void(__GMRFLib_FuncName, secname, "ncol", tmp_ncol);
			}
			if (mb->f_ncol[mb->nf] <= 0) {
				inla_error_field_is_void(__GMRFLib_FuncName, secname, "ncol", tmp_ncol);
			}
		}
		mb->f_N[mb->nf] = mb->f_n[mb->nf] = mb->f_nrow[mb->nf] * mb->f_ncol[mb->nf];
		if (mb->f_id[mb->nf] == F_RW2DIID) {
			/*
			 * we need this for the RW2DIID model otherwise the locations are not set correctly (the length is wrong...)
			 */
			mb->f_N[mb->nf] = 2 * mb->f_n[mb->nf];
		}

		if (mb->verbose) {
			printf("\t\tnrow[%1d] x ncol[%1d] = n[%1d] \n", mb->f_nrow[mb->nf], mb->f_ncol[mb->nf], mb->f_n[mb->nf]);
		}
		mb->f_cyclic[mb->nf] = iniparser_getboolean(ini, inla_string_join(secname, "CYCLIC"), 0);
		if (mb->verbose) {
			printf("\t\tcyclic=[%1d]\n", mb->f_cyclic[mb->nf]);
		}

		bvalue = iniparser_getint(ini, inla_string_join(secname, "BVALUE"), GMRFLib_BVALUE_DEFAULT);
		if (bvalue != GMRFLib_BVALUE_DEFAULT && bvalue != GMRFLib_BVALUE_ZERO) {
			GMRFLib_sprintf(&msg, "%d", bvalue);
			inla_error_field_is_void(__GMRFLib_FuncName, secname, "BVALUE", (const char *) msg);
		}
		if (mb->verbose) {
			printf("\t\tbvalue=[%s]\n", (bvalue == GMRFLib_BVALUE_DEFAULT ? "Default" : "Zero outside"));
		}
		/*
		 * read the covariates 
		 */
		filenamec = Strdup(iniparser_getstring(ini, inla_string_join(secname, "COVARIATES"), NULL));
		if (!filenamec) {
			inla_error_missing_required_field(__GMRFLib_FuncName, secname, "covariates");
		}
		if (mb->verbose) {
			printf("\t\tread covariates from file=[%s]\n", filenamec);
		}
		inla_read_data_general(NULL, &(mb->f_c[mb->nf]), &nc, filenamec, mb->predictor_n, 0, 1, mb->verbose, -1.0);

		/*
		 * read weights 
		 */
		filenamec = Strdup(iniparser_getstring(ini, inla_string_join(secname, "WEIGHTS"), NULL));
		if (!filenamec) {
			mb->f_weights[mb->nf] = NULL;	       /* default */
		} else {
			if (mb->verbose) {
				printf("\t\tread weights from file=[%s]\n", filenamec);
			}
			inla_read_data_general(&(mb->f_weights[mb->nf]), NULL, NULL, filenamec, mb->predictor_n, 0, 1, mb->verbose, 1.0);
		}
	} else {
		/*
		 * read the covariates 
		 */
		filenamec = Strdup(iniparser_getstring(ini, inla_string_join(secname, "COVARIATES"), NULL));
		if (!filenamec) {
			inla_error_missing_required_field(__GMRFLib_FuncName, secname, "covariates");
		}
		if (mb->verbose) {
			printf("\t\tread covariates from file=[%s]\n", filenamec);
		}
		inla_read_data_general(NULL, &(mb->f_c[mb->nf]), &nc, filenamec, mb->predictor_n, 0, 1, mb->verbose, -1.0);

		/*
		 * read weights 
		 */
		filenamec = Strdup(iniparser_getstring(ini, inla_string_join(secname, "WEIGHTS"), NULL));
		if (!filenamec) {
			mb->f_weights[mb->nf] = NULL;
		} else {
			if (mb->verbose) {
				printf("\t\tread weights from file=[%s]\n", filenamec);
			}
			inla_read_data_general(&(mb->f_weights[mb->nf]), NULL, NULL, filenamec, mb->predictor_n, 0, 1, mb->verbose, 1.0);
		}

		switch (mb->f_id[mb->nf]) {
		case F_GENERIC0:
		{
			/*
			 * use field: QMATRIX to set both graph and n.
			 */
			GMRFLib_tabulate_Qfunc_tp *tab = NULL;

			filename = Strdup(iniparser_getstring(ini, inla_string_join(secname, "CMATRIX"), NULL));
			if (!filename) {
				inla_error_missing_required_field(__GMRFLib_FuncName, secname, "Cmatrix");
			}
			if (mb->verbose) {
				printf("\t\tread Cmatrix from file=[%s]\n", filename);
			}
			GMRFLib_tabulate_Qfunc_from_file(&tab, &(mb->f_graph[mb->nf]), (const char *) filename, -1, log_prec);
			mb->f_Qfunc[mb->nf] = tab->Qfunc;
			mb->f_Qfunc_arg[mb->nf] = tab->Qfunc_arg;
			mb->f_locations[mb->nf] = NULL;
			mb->f_N[mb->nf] = mb->f_n[mb->nf] = mb->f_graph[mb->nf]->n;
			mb->f_rankdef[mb->nf] = 0.0;	       /* default */
			if (mb->verbose) {
				for (i = 0; i < IMIN(PREVIEW, mb->f_graph[mb->nf]->n); i++) {
					printf("\t\t\tQ(%1d,%1d) = %g\n", i, i, tab->Qfunc(thread_id, i, i, NULL, tab->Qfunc_arg));
				}
			}
		}
			break;

		case F_GENERIC1:
		{
			/*
			 * use field: CMATRIX to set both graph and n.
			 */
			inla_generic1_tp *arg = Calloc(1, inla_generic1_tp);
			int nn = -1;
			GMRFLib_graph_tp *g = NULL;

			filename = Strdup(iniparser_getstring(ini, inla_string_join(secname, "CMATRIX"), NULL));
			if (!filename) {
				inla_error_missing_required_field(__GMRFLib_FuncName, secname, "Cmatrix");
			}
			if (mb->verbose) {
				printf("\t\tread Cmatrix from file=[%s]\n", filename);
			}
			GMRFLib_tabulate_Qfunc_from_file(&(arg->tab), &(mb->f_graph[mb->nf]), (const char *) filename, -1, NULL);
			arg->log_prec = log_prec;
			arg->beta = beta_intern;
			arg->n = nn = mb->f_graph[mb->nf]->n;

			/*
			 * we need to compute the eigenvalues 
			 */
			gsl_matrix *C = NULL;
			gsl_vector *evalues = NULL;
			gsl_eigen_symm_workspace *w = NULL;

			g = mb->f_graph[mb->nf];
			C = gsl_matrix_calloc(nn, nn);
			evalues = gsl_vector_calloc(nn);
			w = gsl_eigen_symm_alloc(nn);
			for (i = 0; i < nn; i++) {
				gsl_matrix_set(C, i, i, arg->tab->Qfunc(thread_id, i, i, NULL, arg->tab->Qfunc_arg));
				for (jj = 0; jj < g->lnnbs[i]; jj++) {
					j = g->lnbs[i][jj];
					double val = arg->tab->Qfunc(thread_id, i, j, NULL, arg->tab->Qfunc_arg);
					gsl_matrix_set(C, i, j, val);
					gsl_matrix_set(C, j, i, val);
				}
			}
			gsl_eigen_symm(C, evalues, w);
			arg->eigenvalues = Calloc(nn, double);
			arg->max_eigenvalue = arg->eigenvalues[0];
			arg->min_eigenvalue = arg->eigenvalues[0];
			for (i = 0; i < nn; i++) {
				arg->eigenvalues[i] = gsl_vector_get(evalues, i);
				arg->max_eigenvalue = DMAX(arg->max_eigenvalue, arg->eigenvalues[i]);
				arg->min_eigenvalue = DMIN(arg->min_eigenvalue, arg->eigenvalues[i]);
			}
			assert(arg->max_eigenvalue > 0.0);
			gsl_eigen_symm_free(w);
			gsl_matrix_free(C);
			gsl_vector_free(evalues);

			if (mb->verbose) {
				printf("\t\tMaxmimum eigenvalue = %.12g\n", arg->max_eigenvalue);
				printf("\t\tMinimum  eigenvalue = %.12g\n", arg->min_eigenvalue);
			}

			mb->f_Qfunc[mb->nf] = Qfunc_generic1;
			mb->f_Qfunc_arg[mb->nf] = (void *) arg;
			mb->f_locations[mb->nf] = NULL;
			mb->f_N[mb->nf] = mb->f_n[mb->nf] = nn;
			mb->f_rankdef[mb->nf] = 0.0;	       /* default */
			if (mb->verbose) {
				for (i = 0; i < IMIN(PREVIEW, mb->f_graph[mb->nf]->n); i++) {
					printf("\t\t\tC(%1d,%1d) = %g\n", i, i, arg->tab->Qfunc(thread_id, i, i, NULL, arg->tab->Qfunc_arg));
				}
			}
		}
			break;

		case F_GENERIC2:
		{
			/*
			 * use field: CMATRIX to set both graph and n.
			 */
			inla_generic2_tp *arg = Calloc(1, inla_generic2_tp);
			int nn = -1, ii;
			GMRFLib_graph_tp *g = NULL;

			filename = Strdup(iniparser_getstring(ini, inla_string_join(secname, "CMATRIX"), NULL));
			if (!filename) {
				inla_error_missing_required_field(__GMRFLib_FuncName, secname, "Cmatrix");
			}
			if (mb->verbose) {
				printf("\t\tread Cmatrix from file=[%s]\n", filename);
			}
			GMRFLib_tabulate_Qfunc_from_file(&(arg->tab), &g, (const char *) filename, -1, NULL);
			arg->log_prec = log_prec;
			arg->h2_intern = h2_intern;
			arg->n = nn = g->n;
			arg->N = 2 * arg->n;

			GMRFLib_ged_tp *ged = NULL;
			GMRFLib_ged_init(&ged, NULL);
			for (ii = 0; ii < nn; ii++) {
				GMRFLib_ged_add(ged, ii, ii + nn);
			}
			GMRFLib_ged_insert_graph(ged, g, nn);
			GMRFLib_graph_free(g);
			GMRFLib_ged_build(&g, ged);
			assert(g->n == 2 * nn);
			GMRFLib_ged_free(ged);

			mb->f_graph[mb->nf] = g;
			mb->f_Qfunc[mb->nf] = Qfunc_generic2;
			mb->f_Qfunc_arg[mb->nf] = (void *) arg;
			mb->f_locations[mb->nf] = NULL;
			mb->f_n[mb->nf] = nn;
			mb->f_N[mb->nf] = 2 * nn;
			mb->f_rankdef[mb->nf] = 0.0;	       /* default */
			if (mb->verbose) {
				for (i = 0; i < IMIN(PREVIEW, mb->f_graph[mb->nf]->n); i++) {
					printf("\t\t\tC(%1d,%1d) = %g\n", i, i, arg->tab->Qfunc(thread_id, i, i, NULL, arg->tab->Qfunc_arg));
				}
			}
		}
			break;

		case F_GENERIC3:
		{
			inla_generic3_tp *arg = Calloc(1, inla_generic3_tp);
			inla_generic3_tp *arg_orig = Calloc(1, inla_generic3_tp);

			arg->n = iniparser_getint(ini, inla_string_join(secname, "generic3.n"), -1);
			arg_orig->n = iniparser_getint(ini, inla_string_join(secname, "generic3.n"), -1);
			arg->m = iniparser_getint(ini, inla_string_join(secname, "generic3.m"), -1);
			arg_orig->m = iniparser_getint(ini, inla_string_join(secname, "generic3.m"), -1);
			if (mb->verbose) {
				printf("\t\tread generic3.n[%1d]\n", arg->n);
				printf("\t\tread generic3.m[%1d]\n", arg->m);
			}
			assert(arg->n > 0);
			assert(arg->m > 0 && arg->m < GENERIC3_MAXTHETA);

			arg->g = Calloc(arg->m, GMRFLib_graph_tp *);
			arg_orig->g = Calloc(arg->m, GMRFLib_graph_tp *);
			arg->tab = Calloc(arg->m, GMRFLib_tabulate_Qfunc_tp *);
			arg_orig->tab = Calloc(arg->m, GMRFLib_tabulate_Qfunc_tp *);
			for (k = 0; k < arg->m; k++) {
				ctmp = NULL;
				GMRFLib_sprintf(&ctmp, "generic3.Cmatrix.%1d", k);
				filename = Strdup(iniparser_getstring(ini, inla_string_join(secname, ctmp), NULL));
				assert(filename != NULL);
				if (mb->verbose) {
					printf("\t\tread Cmatrix[[%1d]] from file=[%s]\n", k, filename);
				}
				GMRFLib_tabulate_Qfunc_from_file(&(arg->tab[k]), &(arg->g[k]), filename, arg->n, NULL);
				GMRFLib_tabulate_Qfunc_from_file(&(arg_orig->tab[k]), &(arg_orig->g[k]), filename, arg_orig->n, NULL);
				Free(ctmp);
				Free(filename);
			}
			GMRFLib_graph_union(&(arg->graph), arg->g, arg->m);
			GMRFLib_graph_union(&(arg_orig->graph), arg_orig->g, arg_orig->m);

			mb->f_graph[mb->nf] = arg->graph;
			mb->f_graph_orig[mb->nf] = arg_orig->graph;
			mb->f_Qfunc[mb->nf] = Qfunc_generic3;
			mb->f_Qfunc_orig[mb->nf] = Qfunc_generic3;
			mb->f_Qfunc_arg[mb->nf] = (void *) arg;
			mb->f_Qfunc_arg_orig[mb->nf] = (void *) arg_orig;
			mb->f_locations[mb->nf] = NULL;
			mb->f_n[mb->nf] = arg->n;
			mb->f_N[mb->nf] = arg->n;
			mb->f_rankdef[mb->nf] = 0.0;	       /* default */
		}
			break;

		case F_BESAG:
		{
			/*
			 * use field: GRAPH. use this to set field N 
			 */
			filename = Strdup(iniparser_getstring(ini, inla_string_join(secname, "GRAPH"), NULL));
			if (!filename) {
				inla_error_missing_required_field(__GMRFLib_FuncName, secname, "graph");
			}
			if (mb->verbose) {
				printf("\t\tread graph from file=[%s]\n", filename);
			}
			GMRFLib_graph_read(&(mb->f_graph[mb->nf]), filename);
			if (mb->f_graph[mb->nf]->n <= 0) {
				GMRFLib_sprintf(&msg, "graph=[%s] has zero size", filename);
				inla_error_general(msg);
			}

			mb->f_locations[mb->nf] = NULL;
			mb->f_N[mb->nf] = mb->f_n[mb->nf] = mb->f_graph[mb->nf]->n;
		}
			break;

		case F_BESAG2:
		{
			/*
			 * use field: GRAPH. use this to set field N 
			 */
			filename = Strdup(iniparser_getstring(ini, inla_string_join(secname, "GRAPH"), NULL));
			if (!filename) {
				inla_error_missing_required_field(__GMRFLib_FuncName, secname, "graph");
			}
			if (mb->verbose) {
				printf("\t\tread graph from file=[%s]\n", filename);
			}
			GMRFLib_graph_read(&(mb->f_graph[mb->nf]), filename);
			if (mb->f_graph[mb->nf]->n <= 0) {
				GMRFLib_sprintf(&msg, "graph=[%s] has zero size", filename);
				inla_error_general(msg);
			}
			mb->f_precision[mb->nf] = iniparser_getdouble(ini, inla_string_join(secname, "PRECISION"), mb->f_precision[mb->nf]);
			if (mb->verbose) {
				printf("\t\tprecision=[%f]\n", mb->f_precision[mb->nf]);
			}
			mb->f_locations[mb->nf] = NULL;
			mb->f_N[mb->nf] = mb->f_n[mb->nf] = 2 * mb->f_graph[mb->nf]->n;	/* YES */
		}
			break;

		case F_BYM:
		case F_BYM2:
		{
			/*
			 * use field: GRAPH. use this to set field N 
			 */
			filename = Strdup(iniparser_getstring(ini, inla_string_join(secname, "GRAPH"), NULL));
			if (!filename) {
				inla_error_missing_required_field(__GMRFLib_FuncName, secname, "graph");
			}
			if (mb->verbose) {
				printf("\t\tread graph from file=[%s]\n", filename);
			}
			GMRFLib_graph_read(&(mb->f_graph[mb->nf]), filename);
			if (mb->f_graph[mb->nf]->n <= 0) {
				GMRFLib_sprintf(&msg, "graph=[%s] has zero size", filename);
				inla_error_general(msg);
			}
			mb->f_locations[mb->nf] = NULL;
			mb->f_n[mb->nf] = mb->f_graph[mb->nf]->n;
			mb->f_N[mb->nf] = 2 * mb->f_n[mb->nf]; /* yes */
		}
			break;

		case F_BESAGPROPER:
		case F_BESAGPROPER2:
		{
			/*
			 * use field: GRAPH. use this to set field N 
			 */
			filename = Strdup(iniparser_getstring(ini, inla_string_join(secname, "GRAPH"), NULL));
			if (!filename) {
				inla_error_missing_required_field(__GMRFLib_FuncName, secname, "graph");
			}
			if (mb->verbose) {
				printf("\t\tread graph from file=[%s]\n", filename);
			}
			GMRFLib_graph_read(&(mb->f_graph[mb->nf]), filename);
			if (mb->f_graph[mb->nf]->n <= 0) {
				GMRFLib_sprintf(&msg, "graph=[%s] has zero size", filename);
				inla_error_general(msg);
			}

			mb->f_locations[mb->nf] = NULL;
			mb->f_N[mb->nf] = mb->f_n[mb->nf] = mb->f_graph[mb->nf]->n;
		}
			break;

		case F_SEASONAL:
		{
			/*
			 * seasonal component; need length N, seasonal length SEASON, and a boolean CYCLIC
			 */
			ptmp = Strdup(iniparser_getstring(ini, inla_string_join(secname, "N"), NULL));
			if (!ptmp) {
				inla_error_missing_required_field(__GMRFLib_FuncName, secname, "N");
			}
			n = iniparser_getint(ini, inla_string_join(secname, "N"), 0);
			if (n <= 0) {
				inla_error_field_is_void(__GMRFLib_FuncName, secname, "N", ptmp);
			}
			if (mb->verbose) {
				printf("\t\tn=[%1d]\n", n);
			}
			Free(ptmp);
			ptmp = Strdup(iniparser_getstring(ini, inla_string_join(secname, "SEASON"), NULL));
			if (!ptmp) {
				inla_error_missing_required_field(__GMRFLib_FuncName, secname, "SEASON");
			}
			s = iniparser_getint(ini, inla_string_join(secname, "SEASON"), 0);
			if (s <= 0 || s > n) {
				inla_error_field_is_void(__GMRFLib_FuncName, secname, "SEASON", ptmp);
			}
			if (mb->verbose) {
				printf("\t\tseason=[%1d]\n", s);
			}
			/*
			 * this special option is only valid for model=iid,rw1,rw2 and locations=default. therefore we do not add it
			 * into the mb->f_.... 
			 */
			mb->f_cyclic[mb->nf] = iniparser_getboolean(ini, inla_string_join(secname, "CYCLIC"), 0);
			if (mb->verbose) {
				printf("\t\tcyclic=[%1d]\n", mb->f_cyclic[mb->nf]);
			}
			Free(ptmp);
			mb->f_locations[mb->nf] = NULL;
			mb->f_N[mb->nf] = mb->f_n[mb->nf] = n;
		}
			break;

		case F_FGN:
		case F_FGN2:
		{
			/*
			 * FGN/FGN2-model; need length N 
			 */
			ptmp = Strdup(iniparser_getstring(ini, inla_string_join(secname, "N"), NULL));
			if (!ptmp) {
				inla_error_missing_required_field(__GMRFLib_FuncName, secname, "N");
			}
			n = iniparser_getint(ini, inla_string_join(secname, "N"), 0);
			if (n <= 0) {
				inla_error_field_is_void(__GMRFLib_FuncName, secname, "N", ptmp);
			}
			if (mb->verbose) {
				printf("\t\tn=[%1d]\n", n);
			}
			Free(ptmp);
			// even though we do not need the precision of the FGN2.
			mb->f_precision[mb->nf] = iniparser_getdouble(ini, inla_string_join(secname, "PRECISION"), mb->f_precision[mb->nf]);
			if (mb->verbose) {
				printf("\t\tprecision=[%f]\n", mb->f_precision[mb->nf]);
			}
			mb->f_order[mb->nf] = iniparser_getint(ini, inla_string_join(secname, "ORDER"), 4);
			if (mb->verbose) {
				printf("\t\torder=[%1d]\n", mb->f_order[mb->nf]);
			}

			mb->f_N[mb->nf] = mb->f_n[mb->nf] = n;
		}
			break;

		case F_AR1:
		{
			/*
			 * AR1-model; need length N and a boolean CYCLIC
			 */
			ptmp = Strdup(iniparser_getstring(ini, inla_string_join(secname, "N"), NULL));
			if (!ptmp) {
				inla_error_missing_required_field(__GMRFLib_FuncName, secname, "N");
			}
			n = iniparser_getint(ini, inla_string_join(secname, "N"), 0);
			if (n <= 0) {
				inla_error_field_is_void(__GMRFLib_FuncName, secname, "N", ptmp);
			}
			if (mb->verbose) {
				printf("\t\tn=[%1d]\n", n);
			}
			Free(ptmp);
			mb->f_cyclic[mb->nf] = iniparser_getboolean(ini, inla_string_join(secname, "CYCLIC"), 0);
			if (mb->verbose) {
				printf("\t\tcyclic=[%1d]\n", mb->f_cyclic[mb->nf]);
			}
			mb->f_N[mb->nf] = mb->f_n[mb->nf] = n;
		}
			break;

		case F_AR:
		{
			/*
			 * AR(p)-model; need length N and order P
			 */
			ptmp = Strdup(iniparser_getstring(ini, inla_string_join(secname, "N"), NULL));
			if (!ptmp) {
				inla_error_missing_required_field(__GMRFLib_FuncName, secname, "N");
			}
			n = iniparser_getint(ini, inla_string_join(secname, "N"), 0);
			if (n <= 0) {
				inla_error_field_is_void(__GMRFLib_FuncName, secname, "N", ptmp);
			}
			if (mb->verbose) {
				printf("\t\tn=[%1d]\n", n);
			}
			Free(ptmp);

			order = iniparser_getint(ini, inla_string_join(secname, "ORDER"), 1);
			if (mb->verbose) {
				printf("\t\torder=[%1d]\n", order);
			}
			mb->f_order[mb->nf] = order;
			mb->f_N[mb->nf] = mb->f_n[mb->nf] = n;
		}
			break;

		case F_AR1C:
		{
			/*
			 * AR1-model; need length N and a boolean CYCLIC
			 */
			ptmp = Strdup(iniparser_getstring(ini, inla_string_join(secname, "N"), NULL));
			if (!ptmp) {
				inla_error_missing_required_field(__GMRFLib_FuncName, secname, "N");
			}
			n = iniparser_getint(ini, inla_string_join(secname, "N"), 0);
			if (n <= 0) {
				inla_error_field_is_void(__GMRFLib_FuncName, secname, "N", ptmp);
			}
			if (mb->verbose) {
				printf("\t\tn=[%1d]\n", n);
			}
			Free(ptmp);
			mb->f_cyclic[mb->nf] = iniparser_getboolean(ini, inla_string_join(secname, "CYCLIC"), 0);
			if (mb->verbose) {
				printf("\t\tcyclic=[%1d]\n", mb->f_cyclic[mb->nf]);
			}
			assert(mb->f_cyclic[mb->nf] == 0);     /* not implemented */

			mb->f_N[mb->nf] = mb->f_n[mb->nf] = n;
		}
			break;

		case F_Z:
		{
			/*
			 * Z-model. Here Z is a n x m matrix, and the dimension of the model is (Z*z,z) which is n+m
			 */
			zn = iniparser_getint(ini, inla_string_join(secname, "z.n"), 0);
			zm = iniparser_getint(ini, inla_string_join(secname, "z.m"), 0);
			if (zn == 0) {
				GMRFLib_sprintf(&ctmp, "%1d", zn);
				inla_error_field_is_void(__GMRFLib_FuncName, secname, "z.n", ctmp);
			}
			if (zm == 0) {
				GMRFLib_sprintf(&ctmp, "%1d", zm);
				inla_error_field_is_void(__GMRFLib_FuncName, secname, "z.m", ctmp);
			}
			/*
			 * if given, then argument N must be equal n+m.
			 */
			itmp = iniparser_getint(ini, inla_string_join(secname, "N"), -1);
			if (itmp > -1 && zn + zm != itmp) {
				GMRFLib_sprintf(&ctmp, "Model z: dim(Z)[2] = %1d  !=  argument.N = %1d", zm, itmp);
				inla_error_general(ctmp);
			}
			if (mb->verbose) {
				printf("\t\tz.n=[%1d]\n", zn);
				printf("\t\tz.m=[%1d]\n", zm);
			}
			mb->f_N[mb->nf] = mb->f_n[mb->nf] = zn + zm;	/* Yes, this is correct */
		}
			break;

		case F_SLM:
		{
			/*
			 * SLM-model. 
			 */
			slm_n = iniparser_getint(ini, inla_string_join(secname, "slm.n"), 0);
			slm_m = iniparser_getint(ini, inla_string_join(secname, "slm.m"), 0);
			if (slm_n <= 0) {
				GMRFLib_sprintf(&ctmp, "%1d", slm_n);
				inla_error_field_is_void(__GMRFLib_FuncName, secname, "slm.n", ctmp);
			}
			if (slm_m < 0) {
				GMRFLib_sprintf(&ctmp, "%1d", slm_m);
				inla_error_field_is_void(__GMRFLib_FuncName, secname, "slm.m", ctmp);
			}

			/*
			 * if given, then argument N must be equal n+m.
			 */
			itmp = iniparser_getint(ini, inla_string_join(secname, "N"), -1);
			if (itmp > -1 && slm_n + slm_m != itmp) {
				GMRFLib_sprintf(&ctmp, "Model SLM: N=%1d != slm.n=%1d + slm.m=%1d\n", itmp, slm_n, slm_m);
				inla_error_general(ctmp);
			}
			if (mb->verbose) {
				printf("\t\tslm.n=[%1d]\n", slm_n);
				printf("\t\tslm.m=[%1d]\n", slm_m);
			}

			slm_rho_min = iniparser_getdouble(ini, inla_string_join(secname, "slm.rho.min"), 0);
			slm_rho_max = iniparser_getdouble(ini, inla_string_join(secname, "slm.rho.max"), 0);
			if (mb->verbose) {
				printf("\t\tslm.rho.min=[%g]\n", slm_rho_min);
				printf("\t\tslm.rho.max=[%g]\n", slm_rho_max);
			}

			mb->f_N[mb->nf] = mb->f_n[mb->nf] = slm_n + slm_m;	/* Yes, this is correct */
		}
			break;

		case F_2DIID:
		{
			/*
			 * 2DIID-model; need length N
			 */
			ptmp = Strdup(iniparser_getstring(ini, inla_string_join(secname, "N"), NULL));
			if (!ptmp) {
				inla_error_missing_required_field(__GMRFLib_FuncName, secname, "N");
			}
			n = iniparser_getint(ini, inla_string_join(secname, "N"), 0);
			if (n <= 0) {
				inla_error_field_is_void(__GMRFLib_FuncName, secname, "N", ptmp);
			}
			if (mb->verbose) {
				printf("\t\tn=[%1d]\n", n);
			}
			Free(ptmp);
			mb->f_N[mb->nf] = mb->f_n[mb->nf] = n;
		}
			break;

		case F_IID1D:
		case F_IID2D:
		case F_IID3D:
		case F_IID4D:
		case F_IID5D:
		{
			/*
			 * IID_WISHART-model; need length N
			 */
			int dim = WISHART_DIM(mb->nf);
			assert(dim > 0);

			ptmp = Strdup(iniparser_getstring(ini, inla_string_join(secname, "N"), NULL));
			if (!ptmp) {
				inla_error_missing_required_field(__GMRFLib_FuncName, secname, "N");
			}
			n = iniparser_getint(ini, inla_string_join(secname, "N"), 0);
			if (n <= 0) {
				inla_error_field_is_void(__GMRFLib_FuncName, secname, "N", ptmp);
			}
			if (!inla_divisible(n, dim)) {
				GMRFLib_sprintf(&msg, "%s: N=%1d is not divisible by %1d", secname, n, dim);
				inla_error_general(msg);
				exit(EXIT_FAILURE);
			}
			if (mb->verbose) {
				printf("\t\tdim=[%1d]\n", dim);
				printf("\t\tn=[%1d]\n", n);
			}
			Free(ptmp);
			mb->f_N[mb->nf] = mb->f_n[mb->nf] = n;
		}
			break;

		case F_IIDKD:
		{
			/*
			 * WISHART-model; need length N
			 */
			int dim = mb->f_order[mb->nf];
			ptmp = Strdup(iniparser_getstring(ini, inla_string_join(secname, "N"), NULL));
			if (!ptmp) {
				inla_error_missing_required_field(__GMRFLib_FuncName, secname, "N");
			}
			n = iniparser_getint(ini, inla_string_join(secname, "N"), 0);
			if (n <= 0) {
				inla_error_field_is_void(__GMRFLib_FuncName, secname, "N", ptmp);
			}
			if (!inla_divisible(n, dim)) {
				GMRFLib_sprintf(&msg, "%s: N=%1d is not divisible by %1d", secname, n, dim);
				inla_error_general(msg);
				exit(EXIT_FAILURE);
			}
			if (mb->verbose) {
				printf("\t\tdim=[%1d]\n", dim);
				printf("\t\tn=[%1d]\n", n);
			}
			Free(ptmp);
			mb->f_N[mb->nf] = mb->f_n[mb->nf] = n;
		}
			break;

		case F_INTSLOPE:
		{
			n = iniparser_getint(ini, inla_string_join(secname, "N"), 0);
			if (n <= 0) {
				char *val = NULL;
				GMRFLib_sprintf(&val, "%1d", n);
				inla_error_field_is_void(__GMRFLib_FuncName, secname, "N", val);
			}
			nstrata = iniparser_getint(ini, inla_string_join(secname, "intslope.nstrata"), 1);
			nsubject = iniparser_getint(ini, inla_string_join(secname, "intslope.nsubject"), 1);
			mb->f_n[mb->nf] = n;
			mb->f_N[mb->nf] = n + 2 * nsubject;
			mb->f_precision[mb->nf] = iniparser_getdouble(ini, inla_string_join(secname, "PRECISION"), mb->f_precision[mb->nf]);
			if (mb->verbose) {
				printf("\t\tn=[%1d]\n", mb->f_n[mb->nf]);
				printf("\t\tN=[%1d]\n", mb->f_N[mb->nf]);
				printf("\t\tnstrata=[%1d]\n", nstrata);
				printf("\t\tnsubject=[%1d]\n", nsubject);
				printf("\t\tprecision=[%f]\n", mb->f_precision[mb->nf]);
			}

			filename = Strdup(iniparser_getstring(ini, inla_string_join(secname, "intslope.def"), NULL));
			assert(filename != NULL);
			intslope_def = GMRFLib_read_fmesher_file(filename, (long int) 0, -1);

			if (mb->verbose) {
				printf("\t\tRead intslope.def[%s] nrow[%1d] ncol[%1d]\n", filename, intslope_def->nrow, intslope_def->ncol);
				printf("\t\t    Idx  Subject Strata  Covariate\n");
				assert(intslope_def->ncol == 3);
				for (i = 0; i < IMIN(intslope_def->nrow, PREVIEW); i++) {
					printf("\t\t%6d %6.1g %6.1g  %12.6f\n", i,
					       GMRFLib_matrix_get(i, INTSLOPE_SUBJECT, intslope_def),
					       GMRFLib_matrix_get(i, INTSLOPE_STRATA, intslope_def),
					       GMRFLib_matrix_get(i, INTSLOPE_Z, intslope_def));
				}
			}
		}
			break;

		case F_SPDE:
		case F_SPDE2:
		case F_SPDE3:
			break;

		case F_IID:
		case F_RW1:
		case F_RW2:
		case F_CRW2:
		case F_OU:
		case F_COPY:
		case F_SCOPY:
		case F_MEC:
		case F_MEB:
		case F_CLINEAR:
		case F_SIGM:
		case F_REVSIGM:
		case F_LOG1EXP:
		case F_LOGDIST:
		case F_R_GENERIC:
		case F_C_GENERIC:
		case F_DMATERN:
		{
			/*
			 * RW-models and OU-model and ME: read LOCATIONS, set N from LOCATIONS, else read field N and use LOCATIONS=DEFAULT.
			 */
			filename = Strdup(file_loc);
			if (!filename) {
				if (mb->verbose) {
					printf("\t\tfile for locations=[(NULL)]: read n...\n");
				}
				ptmp = Strdup(iniparser_getstring(ini, inla_string_join(secname, "N"), NULL));
				if (!ptmp) {
					GMRFLib_sprintf(&msg, "%s: section[%s]: LOCATIONS is NULL hence N is required",
							__GMRFLib_FuncName, secname);
					inla_error_general(msg);
				} else {
					mb->f_n[mb->nf] = iniparser_getint(ini, inla_string_join(secname, "N"), -1);
					if (mb->f_n[mb->nf] <= 0) {
						inla_error_field_is_void(__GMRFLib_FuncName, secname, "N", ptmp);
					}
				}
				mb->f_N[mb->nf] = mb->f_n[mb->nf];
				mb->f_locations[mb->nf] = NULL;
				if (mb->verbose) {
					printf("\t\tn=[%1d]: use default locations, if required\n", mb->f_n[mb->nf]);
				}
				/*
				 * this special option is only valid for model=iid,rw1,rw2 and locations=default. therefore we do not add it
				 * into the mb->f_.... 
				 */
				mb->f_cyclic[mb->nf] = iniparser_getboolean(ini, inla_string_join(secname, "CYCLIC"), 0);
				if (mb->verbose) {
					printf("\t\tcyclic=[%1d]\n", mb->f_cyclic[mb->nf]);
				}
			} else {
				if (mb->verbose) {
					printf("\t\tfile for locations=[%s]\n", filename);
				}
				inla_read_data_all(&(mb->f_locations[mb->nf]), &nlocations, filename, NULL);

				/*
				 * if N is set, make sure it match with NLOCATIONS 
				 */
				mb->f_n[mb->nf] = iniparser_getint(ini, inla_string_join(secname, "N"), -99);
				if (mb->f_n[mb->nf] != -99 && nlocations != mb->f_n[mb->nf]) {
					GMRFLib_sprintf(&msg, "Number of locations and N does not match: %d != %d\n", nlocations, mb->f_n[mb->nf]);
					inla_error_general(msg);
					exit(EXIT_FAILURE);
				}
				mb->f_N[mb->nf] = mb->f_n[mb->nf] = nlocations;
				if (mb->verbose) {
					printf("\t\t\tnlocations=[%1d]\n", nlocations);
					for (i = 0; i < IMIN(PREVIEW, nlocations); i++) {
						printf("\t\t\tlocations[%1d]=[%g]\n", i, mb->f_locations[mb->nf][i]);
					}
				}
				/*
				 * the locations must be sorted (for some models only), otherwise, things are messed up!!!!
				 */
				if (mb->f_id[mb->nf] == F_RW1 || mb->f_id[mb->nf] == F_RW2 || mb->f_id[mb->nf] == F_CRW2
				    || mb->f_id[mb->nf] == F_OU) {
					for (i = 0; i < nlocations - 1; i++) {
						if (mb->f_locations[mb->nf][i] >= mb->f_locations[mb->nf][i + 1]) {
							inla_error_file_error_sorted(__GMRFLib_FuncName, filename, nlocations, i,
										     mb->f_locations[mb->nf][i]);
						}
					}
				}
				mb->f_cyclic[mb->nf] = iniparser_getboolean(ini, inla_string_join(secname, "CYCLIC"), 0);
				if (mb->verbose) {
					printf("\t\tcyclic=[%1d]\n", mb->f_cyclic[mb->nf]);
				}
			}
		}
			break;

		default:
			/*
			 * if this happens, its an error, as I have forgot one of the models in the list....
			 */
			GMRFLib_ASSERT(0 == 1, GMRFLib_ESNH);
			break;
		}
	}

	// read locations also here if not read before
	if (!mb->f_locations[mb->nf] && mb->f_n[mb->nf] > 0) {
		filename = Strdup(file_loc);
		if (filename)
			if (mb->verbose) {
				printf("\t\tfile for locations=[%s]\n", filename);
			}
		inla_read_data_all(&(mb->f_locations[mb->nf]), &nlocations, filename, NULL);
		if (mb->f_N[mb->nf] > nlocations) {
			double *t = Calloc(mb->f_N[mb->nf], double);
			Memcpy(t, mb->f_locations[mb->nf], nlocations * sizeof(double));
			Free(mb->f_locations[mb->nf]);
			mb->f_locations[mb->nf] = t;
			for (i = nlocations; i < mb->f_N[mb->nf]; i++) {
				mb->f_locations[mb->nf][i] = i + 1;
			}
		}
		if (nlocations != mb->f_n[mb->nf]) {
			GMRFLib_sprintf(&msg, "Number of locations and N does not match: %d != %d\n", nlocations, mb->f_n[mb->nf]);
			inla_error_general(msg);
			exit(EXIT_FAILURE);
		}
		if (mb->verbose) {
			printf("\t\t\tnlocations=[%1d]\n", nlocations);
			for (i = 0; i < IMIN(PREVIEW, nlocations); i++) {
				printf("\t\t\tlocations[%1d]=[%g]\n", i, mb->f_locations[mb->nf][i]);
			}
		}
	}

	switch (mb->f_id[mb->nf]) {
	case F_RW2D:
	case F_BESAG:
	case F_GENERIC0:
	case F_SEASONAL:
	case F_IID:
	case F_RW1:
	case F_RW2:
	case F_CRW2:
	case F_Z:
	{
		tmp = iniparser_getdouble(ini, inla_string_join(secname, "INITIAL"), G.log_prec_initial);
		if (!mb->f_fixed[mb->nf][0] && mb->mode_use_mode) {
			tmp = mb->theta_file[mb->theta_counter_file++];
			if (mb->mode_fixed) {
				mb->f_fixed[mb->nf][0] = 1;
			}
		}
		_SetInitial(0, tmp);
		HYPER_INIT(log_prec, tmp);
		if (mb->verbose) {
			printf("\t\tinitialise log_precision[%g]\n", tmp);
			printf("\t\tfixed=[%1d]\n", mb->f_fixed[mb->nf][0]);
		}

		mb->f_theta[mb->nf] = Calloc(1, double **);
		mb->f_theta[mb->nf][0] = log_prec;
		if (!mb->f_fixed[mb->nf][0]) {
			/*
			 * add this \theta 
			 */
			mb->theta = Realloc(mb->theta, mb->ntheta + 1, double **);
			mb->theta_hyperid = Realloc(mb->theta_hyperid, mb->ntheta + 1, char *);
			mb->theta_hyperid[mb->ntheta] = mb->f_prior[mb->nf][0].hyperid;
			mb->theta_tag = Realloc(mb->theta_tag, mb->ntheta + 1, char *);
			mb->theta_tag_userscale = Realloc(mb->theta_tag_userscale, mb->ntheta + 1, char *);
			mb->theta_dir = Realloc(mb->theta_dir, mb->ntheta + 1, char *);
			GMRFLib_sprintf(&msg, "Log precision for %s", (secname ? secname : mb->f_tag[mb->nf]));
			mb->theta_tag[mb->ntheta] = msg;
			GMRFLib_sprintf(&msg, "Precision for %s", (secname ? secname : mb->f_tag[mb->nf]));
			mb->theta_tag_userscale[mb->ntheta] = msg;
			GMRFLib_sprintf(&msg, "%s-parameter0", mb->f_dir[mb->nf]);
			mb->theta_dir[mb->ntheta] = msg;

			mb->theta_from = Realloc(mb->theta_from, mb->ntheta + 1, char *);
			mb->theta_to = Realloc(mb->theta_to, mb->ntheta + 1, char *);
			mb->theta_from[mb->ntheta] = Strdup(mb->f_prior[mb->nf][0].from_theta);
			mb->theta_to[mb->ntheta] = Strdup(mb->f_prior[mb->nf][0].to_theta);

			mb->theta[mb->ntheta] = log_prec;
			mb->theta_map = Realloc(mb->theta_map, mb->ntheta + 1, map_func_tp *);
			mb->theta_map[mb->ntheta] = map_precision;
			mb->theta_map_arg = Realloc(mb->theta_map_arg, mb->ntheta + 1, void *);
			mb->theta_map_arg[mb->ntheta] = NULL;
			mb->ntheta++;
		}
	}
		break;

	case F_SPDE:
	{
		char *spde_prefix = NULL;
		int nT, nK;

		spde_prefix = Strdup(".");
		spde_prefix = iniparser_getstring(ini, inla_string_join(secname, "SPDE_PREFIX"), spde_prefix);
		spde_prefix = iniparser_getstring(ini, inla_string_join(secname, "SPDE.PREFIX"), spde_prefix);
		spde_prefix = iniparser_getstring(ini, inla_string_join(secname, "SPDEPREFIX"), spde_prefix);
		if (mb->verbose) {
			printf("\t\tspde.prefix = [%s]\n", spde_prefix);
		}

		/*
		 * 
		 */
		inla_spde_build_model(thread_id, &spde_model_orig, (const char *) spde_prefix);
		mb->f_model[mb->nf] = (void *) spde_model_orig;

		/*
		 * The _userfunc0 must be set directly after the _build_model() call. This is a bit dirty; FIXME later. 
		 */
		inla_spde_build_model(thread_id, &spde_model, (const char *) spde_prefix);
		GMRFLib_ai_INLA_userfunc0 = (GMRFLib_ai_INLA_userfunc0_tp *) inla_spde_userfunc0;
		GMRFLib_ai_INLA_userfunc1 = (GMRFLib_ai_INLA_userfunc1_tp *) inla_spde_userfunc1;
		// GMRFLib_ai_INLA_userfunc0_dim = mb->ntheta; NOT NEEDED; set in userfunc0
		GMRFLib_ai_INLA_userfunc1_dim = mb->ntheta;    /* this is a hack and gives the offset of theta... */

		double initial_t = 0.0, initial_k = 0.0, initial_rest = 0.0, initial_oc = 0.0;

		/*
		 * reread this here, as we need non-std defaults 
		 */
		mb->f_fixed[mb->nf][3] = iniparser_getboolean(ini, inla_string_join(secname, "FIXED3"), 1);	/* default fixed, yes */

		// P(mb->nf);
		// P(mb->f_fixed[mb->nf][3]);

		initial_t = iniparser_getdouble(ini, inla_string_join(secname, "INITIAL0"), 0.0);
		initial_k = iniparser_getdouble(ini, inla_string_join(secname, "INITIAL1"), 0.0);
		initial_rest = iniparser_getdouble(ini, inla_string_join(secname, "INITIAL2"), 0.0);
		initial_oc = iniparser_getdouble(ini, inla_string_join(secname, "INITIAL3"), -20.0);
		_SetInitial(0, initial_t);
		_SetInitial(1, initial_k);
		_SetInitial(2, initial_rest);
		_SetInitial(3, initial_oc);

		nT = spde_model->Tmodel->ntheta;
		nK = spde_model->Kmodel->ntheta;
		mb->f_ntheta[mb->nf] = IMAX(4, nT + nK + 1);
		mb->f_Tmodel[mb->nf] = Strdup("basisT");
		mb->f_Kmodel[mb->nf] = Strdup("basisK");

		if (mb->verbose) {
			printf("\t\tnT=[%d]\n", nT);
			printf("\t\tnK=[%d]\n", nK);
			printf("\t\tinitialise theta_t=[%g]\n", initial_t);
			printf("\t\tinitialise theta_k=[%g]\n", initial_k);
			printf("\t\tinitialise theta_rest=[%g]\n", initial_rest);
			printf("\t\tinitialise theta_oc=[%g]\n", initial_oc);
			printf("\t\tfixed_t=[%1d]\n", mb->f_fixed[mb->nf][0]);
			printf("\t\tfixed_k=[%1d]\n", mb->f_fixed[mb->nf][1]);
			printf("\t\tfixed_rest=[%1d]\n", mb->f_fixed[mb->nf][2]);
			printf("\t\tfixed_oc=[%1d]\n", mb->f_fixed[mb->nf][3]);
		}

		for (k = 0; k < nT; k++) {
			if (k == 0) {
				if (!mb->f_fixed[mb->nf][0] && mb->mode_use_mode) {
					tmp = mb->theta_file[mb->theta_counter_file++];
					if (mb->mode_fixed) {
						mb->f_fixed[mb->nf][0] = 1;
					}
				} else {
					tmp = initial_t;
				}
			} else {
				if (!mb->f_fixed[mb->nf][2] && mb->mode_use_mode) {
					tmp = mb->theta_file[mb->theta_counter_file++];
					if (mb->mode_fixed) {
						mb->f_fixed[mb->nf][2] = 1;
					}
				} else {
					tmp = initial_rest;
				}
			}
			HYPER_INIT(spde_model->Tmodel->theta[k], tmp);
		}
		for (k = 0; k < nK; k++) {
			if (k == 0) {
				if (!mb->f_fixed[mb->nf][1] && mb->mode_use_mode) {
					tmp = mb->theta_file[mb->theta_counter_file++];
					if (mb->mode_fixed) {
						mb->f_fixed[mb->nf][1] = 1;
					}
				} else {
					tmp = initial_k;
				}
			} else {
				if (!mb->f_fixed[mb->nf][2] && mb->mode_use_mode) {
					tmp = mb->theta_file[mb->theta_counter_file++];
					if (mb->mode_fixed) {
						mb->f_fixed[mb->nf][2] = 1;
					}
				} else {
					tmp = initial_rest;
				}
			}
			HYPER_INIT(spde_model->Kmodel->theta[k], tmp);
		}
		if (!mb->f_fixed[mb->nf][3] && mb->mode_use_mode) {
			tmp = mb->theta_file[mb->theta_counter_file++];
			if (mb->mode_fixed) {
				mb->f_fixed[mb->nf][3] = 1;
			}
		} else {
			tmp = initial_oc;
		}
		HYPER_INIT(spde_model->oc, tmp);

		mb->f_theta[mb->nf] = Calloc(nT + nK + 1, double **);
		for (k = 0; k < nT; k++) {
			mb->f_theta[mb->nf][k] = spde_model->Tmodel->theta[k];
		}
		for (k = 0; k < nK; k++) {
			mb->f_theta[mb->nf][k + nT] = spde_model->Kmodel->theta[k];
		}
		mb->f_theta[mb->nf][nK + nT] = spde_model->oc;

		for (k = 0; k < nT + nK + 1; k++) {
			int fx;

			if (k == 0) {
				fx = mb->f_fixed[mb->nf][0];   /* T[0] */
			} else if (k == nT) {
				fx = mb->f_fixed[mb->nf][1];   /* K[0] */
			} else if (k == nT + nK) {
				fx = mb->f_fixed[mb->nf][3];   /* oc */
			} else {
				fx = mb->f_fixed[mb->nf][2];   /* the rest */
			}

			if (!fx) {
				/*
				 * add this \theta 
				 */
				mb->theta = Realloc(mb->theta, mb->ntheta + 1, double **);
				mb->theta_hyperid = Realloc(mb->theta_hyperid, mb->ntheta + 1, char *);
				mb->theta_hyperid[mb->ntheta] = NULL;
				mb->theta_tag = Realloc(mb->theta_tag, mb->ntheta + 1, char *);
				mb->theta_tag_userscale = Realloc(mb->theta_tag_userscale, mb->ntheta + 1, char *);
				mb->theta_dir = Realloc(mb->theta_dir, mb->ntheta + 1, char *);

				if (k == nT + nK) {
					GMRFLib_sprintf(&msg, "%s for %s", "Oc", (secname ? secname : mb->f_tag[mb->nf]));
				} else {
					if (k < nT) {
						GMRFLib_sprintf(&msg, "%s.%1d for %s-%s", "T", k,
								(secname ? secname : mb->f_tag[mb->nf]), mb->f_Tmodel[mb->nf]);
					} else {
						GMRFLib_sprintf(&msg, "%s.%1d for %s-%s", "K", k - nT,
								(secname ? secname : mb->f_tag[mb->nf]), mb->f_Kmodel[mb->nf]);
					}
				}
				mb->theta_tag[mb->ntheta] = msg;
				mb->theta_tag_userscale[mb->ntheta] = msg;

				if (k == nT + nK) {
					GMRFLib_sprintf(&msg, "%s-parameter-Oc", mb->f_dir[mb->nf]);
				} else {
					if (k < nT) {
						GMRFLib_sprintf(&msg, "%s-parameter-T.%1d-%s", mb->f_dir[mb->nf], k, mb->f_Tmodel[mb->nf]);
					} else {
						GMRFLib_sprintf(&msg, "%s-parameter-K.%1d-%s", mb->f_dir[mb->nf], k - nT, mb->f_Kmodel[mb->nf]);
					}
				}
				mb->theta_dir[mb->ntheta] = msg;

				if (k == nT + nK) {
					mb->theta[mb->ntheta] = spde_model->oc;
				} else {
					if (k < nT) {
						mb->theta[mb->ntheta] = spde_model->Tmodel->theta[k];
					} else {
						mb->theta[mb->ntheta] = spde_model->Kmodel->theta[k - nT];
					}
				}

				mb->theta_map = Realloc(mb->theta_map, mb->ntheta + 1, map_func_tp *);
				if (k == nT + nK) {
					mb->theta_map[mb->ntheta] = map_probability;
				} else {
					mb->theta_map[mb->ntheta] = map_identity;
				}
				mb->theta_map_arg = Realloc(mb->theta_map_arg, mb->ntheta + 1, void *);
				mb->theta_map_arg[mb->ntheta] = NULL;

				mb->theta_from = Realloc(mb->theta_from, mb->ntheta + 1, char *);
				mb->theta_to = Realloc(mb->theta_to, mb->ntheta + 1, char *);

				int pri;
				if (k == 0) {
					pri = 0;
				} else if (k == nT) {
					pri = 1;
				} else if (k == nT + nK) {
					pri = 3;
				} else {
					pri = 2;
				}

				mb->theta_from[mb->ntheta] = Strdup(mb->f_prior[mb->nf][pri].from_theta);
				mb->theta_to[mb->ntheta] = Strdup(mb->f_prior[mb->nf][pri].to_theta);

				mb->ntheta++;
			}
		}
	}
		break;

	case F_SPDE2:
	{
		char *spde2_prefix, *transform = NULL;

		spde2_prefix = Strdup(".");
		spde2_prefix = iniparser_getstring(ini, inla_string_join(secname, "SPDE2.PREFIX"), spde2_prefix);
		if (mb->verbose) {
			printf("\t\tspde2.prefix = [%s]\n", spde2_prefix);
		}

		transform = Strdup("logit");
		transform = iniparser_getstring(ini, inla_string_join(secname, "SPDE2.TRANSFORM"), transform);
		if (mb->verbose) {
			printf("\t\tspde2.transform = [%s]\n", transform);
		}

		/*
		 * need to read this twice. can save memory by changing the pointer from spde2_model_orig to spde2_model, like
		 * for B and M matrices and BLC. maybe do later
		 */
		inla_spde2_build_model(thread_id, &spde2_model_orig, (const char *) spde2_prefix, (const char *) transform);
		mb->f_model[mb->nf] = (void *) spde2_model_orig;

		inla_spde2_build_model(thread_id, &spde2_model, (const char *) spde2_prefix, (const char *) transform);

		/*
		 * set up userfunc2 that computes the marginal of BLC*theta.intern 
		 */
		GMRFLib_ai_INLA_userfunc2_n++;

		GMRFLib_ai_INLA_userfunc2_args = Realloc(GMRFLib_ai_INLA_userfunc2_args, GMRFLib_ai_INLA_userfunc2_n, void *);
		GMRFLib_ai_INLA_userfunc2_args[GMRFLib_ai_INLA_userfunc2_n - 1] = (void *) spde2_model;
		GMRFLib_ai_INLA_userfunc2 = Realloc(GMRFLib_ai_INLA_userfunc2, GMRFLib_ai_INLA_userfunc2_n, GMRFLib_ai_INLA_userfunc2_tp *);
		GMRFLib_ai_INLA_userfunc2[GMRFLib_ai_INLA_userfunc2_n - 1] = (GMRFLib_ai_INLA_userfunc2_tp *) inla_spde2_userfunc2;

		char *ltag = NULL;
		GMRFLib_sprintf(&ltag, "%s", secname);
		GMRFLib_ai_INLA_userfunc2_tag = Realloc(GMRFLib_ai_INLA_userfunc2_tag, GMRFLib_ai_INLA_userfunc2_n, char *);
		GMRFLib_ai_INLA_userfunc2_tag[GMRFLib_ai_INLA_userfunc2_n - 1] = ltag;

		/*
		 * we know the number of maximum number of hyperparameters
		 */
		int ntheta = spde2_model->ntheta;
		mb->f_initial[mb->nf] = Calloc(ntheta, double);	/* need to do this here as we do not know n_theta upfront */
		if (mb->verbose) {
			printf("\t\tntheta (max) = [%1d]\n", ntheta);
		}

		mb->f_fixed[mb->nf] = Calloc(ntheta, int);
		mb->f_theta[mb->nf] = Calloc(ntheta, double **);
		mb->f_ntheta[mb->nf] = ntheta;

		/*
		 * mark all possible as read 
		 */
		// mark all as read
		for (i = 0; i < SPDE2_MAXTHETA; i++) {
			for (j = 0; j < keywords_len; j++) {
				GMRFLib_sprintf(&ctmp, "%s%1d", keywords[j], i);
				iniparser_getstring(ini, inla_string_join(secname, ctmp), NULL);
				Free(ctmp);
			}
		}

		/*
		 * need to store where in the theta-list the spde2 parameters are 
		 */
		spde2_model->theta_first_idx = mb->ntheta;

		/*
		 * then read those we need 
		 */
		int ntheta_used = 0;			       /* then we count */
		int ignore_prior_error = 0;

		for (i = 0; i < ntheta; i++) {
			double theta_initial = 0.0;

			GMRFLib_sprintf(&ctmp, "FIXED%1d", i);
			mb->f_fixed[mb->nf][i] = iniparser_getboolean(ini, inla_string_join(secname, ctmp), 0);

			GMRFLib_sprintf(&ctmp, "INITIAL%1d", i);
			theta_initial = iniparser_getdouble(ini, inla_string_join(secname, ctmp), theta_initial);
			if (!mb->f_fixed[mb->nf][i] && mb->mode_use_mode) {
				theta_initial = mb->theta_file[mb->theta_counter_file++];
				if (mb->mode_fixed) {
					mb->f_fixed[mb->nf][i] = 1;
					ignore_prior_error = 1;
				}
			}

			HYPER_INIT(spde2_model->theta[i], theta_initial);
			HYPER_INIT(spde2_model_orig->theta[i], theta_initial);
			if (mb->verbose) {
				printf("\t\tinitialise theta[%1d]=[%g]\n", i, theta_initial);
				printf("\t\tfixed[%1d]=[%1d]\n", i, mb->f_fixed[mb->nf][i]);
			}

			GMRFLib_sprintf(&ctmp, "HYPERID%1d", i);
			char *hid = NULL, *cctmp = iniparser_getstring(ini, inla_string_join(secname, ctmp), NULL);
			GMRFLib_sprintf(&hid, "%s|%s", cctmp, secname);

			/*
			 * add this \theta. make spesific names for the pcmatern-case
			 */
			if (!mb->f_fixed[mb->nf][i]) {
				mb->theta = Realloc(mb->theta, mb->ntheta + 1, double **);
				mb->theta_hyperid = Realloc(mb->theta_hyperid, mb->ntheta + 1, char *);
				mb->theta_hyperid[mb->ntheta] = hid;
				mb->theta_tag = Realloc(mb->theta_tag, mb->ntheta + 1, char *);
				mb->theta_tag_userscale = Realloc(mb->theta_tag_userscale, mb->ntheta + 1, char *);
				mb->theta_dir = Realloc(mb->theta_dir, mb->ntheta + 1, char *);

				if (mb->f_prior[mb->nf][0].id == P_PC_MATERN && (i == 0 || i == 1)) {
					GMRFLib_sprintf(&msg, "%s for %s", (i == 0 ? "log(Range)" : "log(Stdev)"),
							(secname ? secname : mb->f_tag[mb->nf]));
					mb->theta_tag[mb->ntheta] = msg;
					GMRFLib_sprintf(&msg, "%s for %s", (i == 0 ? "Range" : "Stdev"), (secname ? secname : mb->f_tag[mb->nf]));
					mb->theta_tag_userscale[mb->ntheta] = msg;
					GMRFLib_sprintf(&msg, "%s-parameter%1d", mb->f_dir[mb->nf], i + 1);
					mb->theta_dir[mb->ntheta] = msg;
					mb->theta_from = Realloc(mb->theta_from, mb->ntheta + 1, char *);
					mb->theta_to = Realloc(mb->theta_to, mb->ntheta + 1, char *);
					mb->theta_from[mb->ntheta] = Strdup("function (x) <<NEWLINE>>exp(x)");	/* they are not there... */
					mb->theta_to[mb->ntheta] = Strdup("function (x) <<NEWLINE>>log(x)");	/* .... */
					mb->theta[mb->ntheta] = spde2_model->theta[i];
					mb->theta_map = Realloc(mb->theta_map, mb->ntheta + 1, map_func_tp *);
					mb->theta_map[mb->ntheta] = map_exp;
				} else {
					GMRFLib_sprintf(&msg, "Theta%1d for %s", i + 1, (secname ? secname : mb->f_tag[mb->nf]));
					mb->theta_tag[mb->ntheta] = msg;
					GMRFLib_sprintf(&msg, "Theta%1d for %s", i + 1, (secname ? secname : mb->f_tag[mb->nf]));
					mb->theta_tag_userscale[mb->ntheta] = msg;
					GMRFLib_sprintf(&msg, "%s-parameter%1d", mb->f_dir[mb->nf], i + 1);
					mb->theta_dir[mb->ntheta] = msg;
					mb->theta_from = Realloc(mb->theta_from, mb->ntheta + 1, char *);
					mb->theta_to = Realloc(mb->theta_to, mb->ntheta + 1, char *);
					mb->theta_from[mb->ntheta] = Strdup(mb->f_prior[mb->nf][0].from_theta);	/* YES, use prior0, which * is a
														 * joint prior */
					mb->theta_to[mb->ntheta] = Strdup(mb->f_prior[mb->nf][0].to_theta);	/* YES, use prior0, which * is a
														 * joint prior */
					mb->theta[mb->ntheta] = spde2_model->theta[i];
					mb->theta_map = Realloc(mb->theta_map, mb->ntheta + 1, map_func_tp *);
					mb->theta_map[mb->ntheta] = map_identity;
				}
				mb->theta_map_arg = Realloc(mb->theta_map_arg, mb->ntheta + 1, void *);
				mb->theta_map_arg[mb->ntheta] = NULL;
				mb->ntheta++;
				ntheta_used++;
			}
		}

		spde2_model->ntheta_used = ntheta_used;
		if (mb->verbose) {
			printf("\t\tntheta (used) = [%1d]\n", spde2_model->ntheta_used);
		}
		spde2_model->fixed = Calloc(ntheta, int);
		spde2_model->fixed_values = Calloc(ntheta, double);
		for (i = 0; i < ntheta; i++) {
			spde2_model->fixed[i] = mb->f_fixed[mb->nf][i];
			if (spde2_model->fixed[i]) {
				spde2_model->fixed_values[i] = spde2_model->theta[i][0][0];
			} else {
				spde2_model->fixed_values[i] = NAN;	/* so that we get an error if used wrongly */
			}
		}
		spde2_model_orig->ntheta_used = spde2_model->ntheta_used;
		spde2_model_orig->fixed = Calloc(ntheta, int);
		spde2_model_orig->fixed_values = Calloc(ntheta, double);
		Memcpy(spde2_model_orig->fixed, spde2_model->fixed, ntheta * sizeof(int));
		Memcpy(spde2_model_orig->fixed_values, spde2_model->fixed_values, ntheta * sizeof(double));

		if (mb->f_prior[mb->nf][0].id == P_MVNORM) {
			if ((int) mb->f_prior[mb->nf][0].parameters[0] != ntheta_used && !ignore_prior_error) {
				GMRFLib_sprintf(&ptmp,
						"Dimension of the MVNORM prior is not equal to number of used hyperparameters: %1d != %1d\n",
						(int) mb->f_prior[mb->nf][0].parameters[0], ntheta_used);
				inla_error_general(ptmp);
				exit(EXIT_FAILURE);
			}
		} else if (!strcasecmp(mb->f_prior[mb->nf][0].name, "PCSPDEGA")) {
			assert(ntheta_used == 2);	       /* requirement... */
		} else if (mb->f_prior[mb->nf][0].id == P_PC_MATERN) {
			assert(ntheta_used <= 2);	       /* requirement... */
		} else {
			GMRFLib_ASSERT(0 == 1, GMRFLib_ESNH);
		}
	}
		break;

	case F_SPDE3:
	{
		char *spde3_prefix, *transform = NULL;

		spde3_prefix = Strdup(".");
		spde3_prefix = iniparser_getstring(ini, inla_string_join(secname, "SPDE3_PREFIX"), spde3_prefix);
		spde3_prefix = iniparser_getstring(ini, inla_string_join(secname, "SPDE3.PREFIX"), spde3_prefix);
		spde3_prefix = iniparser_getstring(ini, inla_string_join(secname, "SPDE3PREFIX"), spde3_prefix);
		if (mb->verbose) {
			printf("\t\tspde3.prefix = [%s]\n", spde3_prefix);
		}

		transform = Strdup("logit");
		transform = iniparser_getstring(ini, inla_string_join(secname, "SPDE3_TRANSFORM"), transform);
		transform = iniparser_getstring(ini, inla_string_join(secname, "SPDE3.TRANSFORM"), transform);
		transform = iniparser_getstring(ini, inla_string_join(secname, "SPDE3TRANSFORM"), transform);
		if (mb->verbose) {
			printf("\t\tspde3.transform = [%s]\n", transform);
		}

		/*
		 * need to read this twice. can save memory by changing the pointer from spde3_model_orig to spde3_model, like for B and M matrices and
		 * BLC. maybe do later
		 */
		inla_spde3_build_model(thread_id, &spde3_model_orig, (const char *) spde3_prefix, (const char *) transform);
		mb->f_model[mb->nf] = (void *) spde3_model_orig;

		inla_spde3_build_model(thread_id, &spde3_model, (const char *) spde3_prefix, (const char *) transform);

		/*
		 * set up userfunc3 that computes the marginal of BLC*theta.intern 
		 */
		GMRFLib_ai_INLA_userfunc3_n++;

		GMRFLib_ai_INLA_userfunc3_args = Realloc(GMRFLib_ai_INLA_userfunc3_args, GMRFLib_ai_INLA_userfunc3_n, void *);
		GMRFLib_ai_INLA_userfunc3_args[GMRFLib_ai_INLA_userfunc3_n - 1] = (void *) spde3_model;
		GMRFLib_ai_INLA_userfunc3 = Realloc(GMRFLib_ai_INLA_userfunc3, GMRFLib_ai_INLA_userfunc3_n, GMRFLib_ai_INLA_userfunc3_tp *);
		GMRFLib_ai_INLA_userfunc3[GMRFLib_ai_INLA_userfunc3_n - 1] = (GMRFLib_ai_INLA_userfunc3_tp *) inla_spde3_userfunc3;

		char *ltag = NULL;
		GMRFLib_sprintf(&ltag, "%s", secname);
		GMRFLib_ai_INLA_userfunc3_tag = Realloc(GMRFLib_ai_INLA_userfunc3_tag, GMRFLib_ai_INLA_userfunc3_n, char *);
		GMRFLib_ai_INLA_userfunc3_tag[GMRFLib_ai_INLA_userfunc3_n - 1] = ltag;

		/*
		 * now we know the number of hyperparameters ;-) 
		 */
		int ntheta;

		mb->f_ntheta[mb->nf] = ntheta = spde3_model->ntheta;
		mb->f_initial[mb->nf] = Calloc(mb->f_ntheta[mb->nf], double);	/* need to do this here as we do not know n_theta upfront */
		if (mb->verbose) {
			printf("\t\tntheta = [%1d]\n", ntheta);
		}

		mb->f_fixed[mb->nf] = Calloc(ntheta, int);
		mb->f_theta[mb->nf] = Calloc(ntheta, double **);

		/*
		 * mark all possible as read
		 */

		// mark all as read
		for (i = 0; i < SPDE3_MAXTHETA; i++) {
			for (j = 0; j < keywords_len; j++) {
				GMRFLib_sprintf(&ctmp, "%s%1d", keywords[j], i);
				iniparser_getstring(ini, inla_string_join(secname, ctmp), NULL);
				Free(ctmp);
			}
		}

		/*
		 * need to know where in the theta-list the spde3 parameters are 
		 */
		spde3_model->theta_first_idx = mb->ntheta;

		/*
		 * then read those we need 
		 */
		int ignore_prior_error = 0;
		for (i = 0; i < ntheta; i++) {
			double theta_initial = 0.0;

			GMRFLib_sprintf(&ctmp, "FIXED%1d", i);
			mb->f_fixed[mb->nf][i] = iniparser_getboolean(ini, inla_string_join(secname, ctmp), 0);
			if (mb->f_fixed[mb->nf][i]) {
				inla_error_general("Fixed hyperparmaters is not allowed in the SPDE3 model.");
				exit(EXIT_FAILURE);
			}

			GMRFLib_sprintf(&ctmp, "INITIAL%1d", i);
			theta_initial = iniparser_getdouble(ini, inla_string_join(secname, ctmp), theta_initial);
			if (!mb->f_fixed[mb->nf][i] && mb->mode_use_mode) {
				theta_initial = mb->theta_file[mb->theta_counter_file++];
				if (mb->mode_fixed) {
					mb->f_fixed[mb->nf][i] = 1;
					ignore_prior_error = 1;
				}
			}

			HYPER_INIT(spde3_model->theta[i], theta_initial);
			HYPER_INIT(spde3_model_orig->theta[i], theta_initial);

			if (mb->verbose) {
				printf("\t\tinitialise theta[%1d]=[%g]\n", i, theta_initial);
				printf("\t\tfixed[%1d]=[%1d]\n", i, mb->f_fixed[mb->nf][i]);
			}

			/*
			 * add this \theta 
			 */
			mb->theta = Realloc(mb->theta, mb->ntheta + 1, double **);
			mb->theta_hyperid = Realloc(mb->theta_hyperid, mb->ntheta + 1, char *);
			mb->theta_hyperid[mb->ntheta] = NULL;
			mb->theta_tag = Realloc(mb->theta_tag, mb->ntheta + 1, char *);
			mb->theta_tag_userscale = Realloc(mb->theta_tag_userscale, mb->ntheta + 1, char *);
			mb->theta_dir = Realloc(mb->theta_dir, mb->ntheta + 1, char *);
			GMRFLib_sprintf(&msg, "Theta%1d for %s", i + 1, (secname ? secname : mb->f_tag[mb->nf]));
			mb->theta_tag[mb->ntheta] = msg;
			GMRFLib_sprintf(&msg, "Theta%1d for %s", i + 1, (secname ? secname : mb->f_tag[mb->nf]));
			mb->theta_tag_userscale[mb->ntheta] = msg;
			GMRFLib_sprintf(&msg, "%s-parameter%1d", mb->f_dir[mb->nf], i + 1);
			mb->theta_dir[mb->ntheta] = msg;

			mb->theta_from = Realloc(mb->theta_from, mb->ntheta + 1, char *);
			mb->theta_to = Realloc(mb->theta_to, mb->ntheta + 1, char *);
			mb->theta_from[mb->ntheta] = Strdup(mb->f_prior[mb->nf][0].from_theta);	/* YES, use prior0 */
			mb->theta_to[mb->ntheta] = Strdup(mb->f_prior[mb->nf][0].to_theta);	/* YES, use prior0 */

			mb->theta[mb->ntheta] = spde3_model->theta[i];
			mb->theta_map = Realloc(mb->theta_map, mb->ntheta + 1, map_func_tp *);
			mb->theta_map[mb->ntheta] = map_identity;
			mb->theta_map_arg = Realloc(mb->theta_map_arg, mb->ntheta + 1, void *);
			mb->theta_map_arg[mb->ntheta] = NULL;
			mb->ntheta++;
		}

		if ((int) mb->f_prior[mb->nf][0].parameters[0] != ntheta && !ignore_prior_error) {
			GMRFLib_sprintf(&ptmp,
					"Dimension of the MVNORM prior is not equal to number of hyperparameters: %1d != %1d\n",
					(int) mb->f_prior[mb->nf][0].parameters[0], ntheta);
			inla_error_general(ptmp);
			exit(EXIT_FAILURE);
		}
	}
		break;

	case F_AR:
	{
		int ntheta, ntheta_ref = mb->ntheta;

		order = mb->f_order[mb->nf];
		ntheta = mb->f_ntheta[mb->nf] = AR_MAXTHETA + 1;
		mb->f_initial[mb->nf] = Calloc(ntheta, double);
		if (mb->verbose) {
			printf("\t\tntheta.max = [%1d]\n", ntheta);
		}

		/*
		 * add the order as a parameter in the prior
		 */
		switch (mb->f_prior[mb->nf][1].id) {
		case P_REF_AR:
		{
			/*
			 * add the order as the first parameter in the ref-prior! 
			 */
			mb->f_prior[mb->nf][1].parameters[0] = mb->f_order[mb->nf];
		}
			break;
		default:
			break;
		}

		mb->f_fixed[mb->nf] = Calloc(ntheta, int);
		mb->f_theta[mb->nf] = Calloc(ntheta, double **);

		HYPER_NEW(log_prec, 0.0);
		mb->f_theta[mb->nf][0] = log_prec;
		pacf_intern = Calloc(AR_MAXTHETA, double **);
		for (i = 0; i < AR_MAXTHETA; i++) {
			HYPER_NEW(pacf_intern[i], 0.0);
			mb->f_theta[mb->nf][i + 1] = pacf_intern[i];
		}

		// mark all as read
		for (i = 0; i < AR_MAXTHETA + 1; i++) {
			for (j = 0; j < keywords_len; j++) {
				GMRFLib_sprintf(&ctmp, "%s%1d", keywords[j], i);
				iniparser_getstring(ini, inla_string_join(secname, ctmp), NULL);
				Free(ctmp);
			}
		}

		/*
		 * then read those we need 
		 */
		for (i = 0; i < ntheta; i++) {
			double theta_initial = 0.0;

			if (i > order) {
				// by definition, this one must be fixed
				mb->f_fixed[mb->nf][i] = 1;
			} else {
				GMRFLib_sprintf(&ctmp, "FIXED%1d", i);
				mb->f_fixed[mb->nf][i] = iniparser_getboolean(ini, inla_string_join(secname, ctmp), 0);
			}

			GMRFLib_sprintf(&ctmp, "INITIAL%1d", i);
			theta_initial = iniparser_getdouble(ini, inla_string_join(secname, ctmp), theta_initial);
			if (!mb->f_fixed[mb->nf][i] && mb->mode_use_mode) {
				theta_initial = mb->theta_file[mb->theta_counter_file++];
				if (mb->mode_fixed) {
					mb->f_fixed[mb->nf][i] = 1;
				}
			}

			if (i == 0) {
				/*
				 * precision 
				 */
				HYPER_INIT(log_prec, theta_initial);
				if (mb->verbose) {
					printf("\t\tinitialise (log_prec) theta[%1d]=[%g]\n", i, theta_initial);
					printf("\t\tfixed[%1d]=[%1d]\n", i, mb->f_fixed[mb->nf][i]);
				}

				if (!mb->f_fixed[mb->nf][i]) {
					/*
					 * add this \theta 
					 */
					mb->theta = Realloc(mb->theta, mb->ntheta + 1, double **);
					mb->theta_hyperid = Realloc(mb->theta_hyperid, mb->ntheta + 1, char *);
					mb->theta_hyperid[mb->ntheta] = mb->f_prior[mb->nf][i].hyperid;
					mb->theta_tag = Realloc(mb->theta_tag, mb->ntheta + 1, char *);
					mb->theta_tag_userscale = Realloc(mb->theta_tag_userscale, mb->ntheta + 1, char *);
					mb->theta_dir = Realloc(mb->theta_dir, mb->ntheta + 1, char *);
					GMRFLib_sprintf(&msg, "Log precision for %s", (secname ? secname : mb->f_tag[mb->nf]));

					mb->theta_tag[mb->ntheta] = msg;
					GMRFLib_sprintf(&msg, "Precision for %s", (secname ? secname : mb->f_tag[mb->nf]));
					mb->theta_tag_userscale[mb->ntheta] = msg;
					GMRFLib_sprintf(&msg, "%s-parameter%1d", mb->f_dir[mb->nf], i + 1);
					mb->theta_dir[mb->ntheta] = msg;

					mb->theta_from = Realloc(mb->theta_from, mb->ntheta + 1, char *);
					mb->theta_to = Realloc(mb->theta_to, mb->ntheta + 1, char *);
					mb->theta_from[mb->ntheta] = Strdup(mb->f_prior[mb->nf][i].from_theta);
					mb->theta_to[mb->ntheta] = Strdup(mb->f_prior[mb->nf][i].to_theta);
					mb->theta[mb->ntheta] = log_prec;
					mb->theta_map = Realloc(mb->theta_map, mb->ntheta + 1, map_func_tp *);
					mb->theta_map[mb->ntheta] = map_precision;
					mb->theta_map_arg = Realloc(mb->theta_map_arg, mb->ntheta + 1, void *);
					mb->theta_map_arg[mb->ntheta] = NULL;
					mb->ntheta++;
				}
			} else {
				/*
				 * PACF 
				 */
				HYPER_INIT(pacf_intern[i - 1], theta_initial);
				if (mb->verbose) {
					printf("\t\tinitialise (PACF) theta[%1d]=[%g]\n", i, theta_initial);
					printf("\t\tfixed[%1d]=[%1d]\n", i, mb->f_fixed[mb->nf][i]);
				}
				if (!mb->f_fixed[mb->nf][i]) {
					/*
					 * add this \theta 
					 */
					mb->theta = Realloc(mb->theta, mb->ntheta + 1, double **);
					mb->theta_hyperid = Realloc(mb->theta_hyperid, mb->ntheta + 1, char *);
					mb->theta_hyperid[mb->ntheta] = mb->f_prior[mb->nf][i].hyperid;
					mb->theta_tag = Realloc(mb->theta_tag, mb->ntheta + 1, char *);
					mb->theta_tag_userscale = Realloc(mb->theta_tag_userscale, mb->ntheta + 1, char *);
					mb->theta_dir = Realloc(mb->theta_dir, mb->ntheta + 1, char *);
					GMRFLib_sprintf(&msg, "Intern PACF%1d for %s", i, (secname ? secname : mb->f_tag[mb->nf]));

					mb->theta_tag[mb->ntheta] = msg;
					GMRFLib_sprintf(&msg, "PACF%1d for %s", i, (secname ? secname : mb->f_tag[mb->nf]));
					mb->theta_tag_userscale[mb->ntheta] = msg;
					GMRFLib_sprintf(&msg, "%s-parameter%1d", mb->f_dir[mb->nf], i + 1);
					mb->theta_dir[mb->ntheta] = msg;

					mb->theta_from = Realloc(mb->theta_from, mb->ntheta + 1, char *);
					mb->theta_to = Realloc(mb->theta_to, mb->ntheta + 1, char *);
					mb->theta_from[mb->ntheta] = Strdup(mb->f_prior[mb->nf][i].from_theta);
					mb->theta_to[mb->ntheta] = Strdup(mb->f_prior[mb->nf][i].to_theta);
					mb->theta[mb->ntheta] = pacf_intern[i - 1];
					mb->theta_map = Realloc(mb->theta_map, mb->ntheta + 1, map_func_tp *);
					mb->theta_map[mb->ntheta] = ar_map_pacf;
					mb->theta_map_arg = Realloc(mb->theta_map_arg, mb->ntheta + 1, void *);
					mb->theta_map_arg[mb->ntheta] = NULL;
					mb->ntheta++;
				}
			}
		}
		if (mb->verbose) {
			printf("\t\tntheta = [%1d]\n", mb->ntheta - ntheta_ref);
		}
	}
		break;

	case F_MEC:
	{
		tmp = iniparser_getdouble(ini, inla_string_join(secname, "INITIAL0"), 1.0);
		if (!mb->f_fixed[mb->nf][0] && mb->mode_use_mode) {
			tmp = mb->theta_file[mb->theta_counter_file++];
			if (mb->mode_fixed) {
				mb->f_fixed[mb->nf][0] = 1;
			}
		}
		_SetInitial(0, tmp);
		HYPER_INIT(beta, tmp);

		double *range = NULL;
		range = Calloc(2, double);		       /* need this as it will be stored in the map argument */
		range[0] = iniparser_getdouble(ini, inla_string_join(secname, "RANGE.LOW"), 0.0);	/* low = high ==> map = identity */
		range[1] = iniparser_getdouble(ini, inla_string_join(secname, "RANGE.HIGH"), 0.0);

		if (mb->verbose) {
			printf("\t\tinitialise beta[%g]\n", tmp);
			printf("\t\tfixed=[%1d]\n", mb->f_fixed[mb->nf][0]);
			printf("\t\trange[%g, %g]\n", range[0], range[1]);
		}

		mb->f_theta_map[mb->nf] = Calloc(1, map_func_tp *);
		mb->f_theta_map_arg[mb->nf] = Calloc(1, void *);
		mb->f_theta_map[mb->nf][0] = map_beta;	       /* need these */
		mb->f_theta_map_arg[mb->nf][0] = (void *) range;	/* and this one as well */
		mb->f_theta[mb->nf][0] = beta;

		if (!mb->f_fixed[mb->nf][0]) {
			/*
			 * add this \theta 
			 */
			mb->theta = Realloc(mb->theta, mb->ntheta + 1, double **);
			mb->theta_hyperid = Realloc(mb->theta_hyperid, mb->ntheta + 1, char *);
			mb->theta_hyperid[mb->ntheta] = mb->f_prior[mb->nf][0].hyperid;
			mb->theta_tag = Realloc(mb->theta_tag, mb->ntheta + 1, char *);
			mb->theta_tag_userscale = Realloc(mb->theta_tag_userscale, mb->ntheta + 1, char *);
			mb->theta_dir = Realloc(mb->theta_dir, mb->ntheta + 1, char *);
			GMRFLib_sprintf(&msg, "MEC beta for %s", (secname ? secname : mb->f_tag[mb->nf]));
			mb->theta_tag[mb->ntheta] = msg;
			GMRFLib_sprintf(&msg, "MEC beta for %s", (secname ? secname : mb->f_tag[mb->nf]));
			mb->theta_tag_userscale[mb->ntheta] = msg;
			GMRFLib_sprintf(&msg, "%s-parameter0", mb->f_dir[mb->nf]);
			mb->theta_dir[mb->ntheta] = msg;

			mb->theta_from = Realloc(mb->theta_from, mb->ntheta + 1, char *);
			mb->theta_to = Realloc(mb->theta_to, mb->ntheta + 1, char *);
			mb->theta_from[mb->ntheta] = Strdup(mb->f_prior[mb->nf][0].from_theta);
			mb->theta_to[mb->ntheta] = Strdup(mb->f_prior[mb->nf][0].to_theta);

			mb->theta[mb->ntheta] = beta;
			mb->theta_map = Realloc(mb->theta_map, mb->ntheta + 1, map_func_tp *);
			mb->theta_map[mb->ntheta] = map_beta;
			mb->theta_map_arg = Realloc(mb->theta_map_arg, mb->ntheta + 1, void *);
			mb->theta_map_arg[mb->ntheta] = (void *) range;
			mb->ntheta++;
		}

		tmp = iniparser_getdouble(ini, inla_string_join(secname, "INITIAL1"), G.log_prec_initial);
		if (!mb->f_fixed[mb->nf][1] && mb->mode_use_mode) {
			tmp = mb->theta_file[mb->theta_counter_file++];
			if (mb->mode_fixed) {
				mb->f_fixed[mb->nf][1] = 1;
			}
		}
		_SetInitial(1, tmp);
		HYPER_INIT(log_prec, tmp);
		if (mb->verbose) {
			printf("\t\tinitialise log_prec_u[%g]\n", tmp);
			printf("\t\tfixed=[%1d]\n", mb->f_fixed[mb->nf][1]);
		}
		mb->f_theta[mb->nf][1] = log_prec;
		if (!mb->f_fixed[mb->nf][1]) {
			/*
			 * add this \theta 
			 */
			mb->theta = Realloc(mb->theta, mb->ntheta + 1, double **);
			mb->theta_hyperid = Realloc(mb->theta_hyperid, mb->ntheta + 1, char *);
			mb->theta_hyperid[mb->ntheta] = mb->f_prior[mb->nf][1].hyperid;
			mb->theta_tag = Realloc(mb->theta_tag, mb->ntheta + 1, char *);
			mb->theta_tag_userscale = Realloc(mb->theta_tag_userscale, mb->ntheta + 1, char *);
			mb->theta_dir = Realloc(mb->theta_dir, mb->ntheta + 1, char *);
			GMRFLib_sprintf(&msg, "MEC prec_u_intern for %s", (secname ? secname : mb->f_tag[mb->nf]));
			mb->theta_tag[mb->ntheta] = msg;
			GMRFLib_sprintf(&msg, "MEC prec_u for %s", (secname ? secname : mb->f_tag[mb->nf]));
			mb->theta_tag_userscale[mb->ntheta] = msg;
			GMRFLib_sprintf(&msg, "%s-parameter1", mb->f_dir[mb->nf]);
			mb->theta_dir[mb->ntheta] = msg;

			mb->theta_from = Realloc(mb->theta_from, mb->ntheta + 1, char *);
			mb->theta_to = Realloc(mb->theta_to, mb->ntheta + 1, char *);
			mb->theta_from[mb->ntheta] = Strdup(mb->f_prior[mb->nf][1].from_theta);
			mb->theta_to[mb->ntheta] = Strdup(mb->f_prior[mb->nf][1].to_theta);

			mb->theta[mb->ntheta] = log_prec;
			mb->theta_map = Realloc(mb->theta_map, mb->ntheta + 1, map_func_tp *);
			mb->theta_map[mb->ntheta] = map_precision;
			mb->theta_map_arg = Realloc(mb->theta_map_arg, mb->ntheta + 1, void *);
			mb->theta_map_arg[mb->ntheta] = NULL;
			mb->ntheta++;
		}

		tmp = iniparser_getdouble(ini, inla_string_join(secname, "INITIAL2"), 0.0);
		if (!mb->f_fixed[mb->nf][2] && mb->mode_use_mode) {
			tmp = mb->theta_file[mb->theta_counter_file++];
			if (mb->mode_fixed) {
				mb->f_fixed[mb->nf][2] = 1;
			}
		}
		_SetInitial(2, tmp);
		HYPER_INIT(mean_x, tmp);
		if (mb->verbose) {
			printf("\t\tinitialise mean_x[%g]\n", tmp);
			printf("\t\tfixed=[%1d]\n", mb->f_fixed[mb->nf][2]);
		}

		mb->f_theta[mb->nf][2] = mean_x;
		if (!mb->f_fixed[mb->nf][2]) {
			/*
			 * add this \theta 
			 */
			mb->theta = Realloc(mb->theta, mb->ntheta + 1, double **);
			mb->theta_hyperid = Realloc(mb->theta_hyperid, mb->ntheta + 1, char *);
			mb->theta_hyperid[mb->ntheta] = mb->f_prior[mb->nf][2].hyperid;
			mb->theta_tag = Realloc(mb->theta_tag, mb->ntheta + 1, char *);
			mb->theta_tag_userscale = Realloc(mb->theta_tag_userscale, mb->ntheta + 1, char *);
			mb->theta_dir = Realloc(mb->theta_dir, mb->ntheta + 1, char *);
			GMRFLib_sprintf(&msg, "MEC mean_x for %s", (secname ? secname : mb->f_tag[mb->nf]));
			mb->theta_tag[mb->ntheta] = msg;
			GMRFLib_sprintf(&msg, "MEC mean_x for %s", (secname ? secname : mb->f_tag[mb->nf]));
			mb->theta_tag_userscale[mb->ntheta] = msg;
			GMRFLib_sprintf(&msg, "%s-parameter2", mb->f_dir[mb->nf]);
			mb->theta_dir[mb->ntheta] = msg;

			mb->theta_from = Realloc(mb->theta_from, mb->ntheta + 1, char *);
			mb->theta_to = Realloc(mb->theta_to, mb->ntheta + 1, char *);
			mb->theta_from[mb->ntheta] = Strdup(mb->f_prior[mb->nf][2].from_theta);
			mb->theta_to[mb->ntheta] = Strdup(mb->f_prior[mb->nf][2].to_theta);

			mb->theta[mb->ntheta] = mean_x;
			mb->theta_map = Realloc(mb->theta_map, mb->ntheta + 1, map_func_tp *);
			mb->theta_map[mb->ntheta] = map_identity;
			mb->theta_map_arg = Realloc(mb->theta_map_arg, mb->ntheta + 1, void *);
			mb->theta_map_arg[mb->ntheta] = NULL;
			mb->ntheta++;
		}

		tmp = iniparser_getdouble(ini, inla_string_join(secname, "INITIAL3"), G.log_prec_initial);
		if (!mb->f_fixed[mb->nf][3] && mb->mode_use_mode) {
			tmp = mb->theta_file[mb->theta_counter_file++];
			if (mb->mode_fixed) {
				mb->f_fixed[mb->nf][3] = 1;
			}
		}
		_SetInitial(3, tmp);
		HYPER_INIT(log_prec_x, tmp);
		if (mb->verbose) {
			printf("\t\tinitialise log_prec_x[%g]\n", tmp);
			printf("\t\tfixed=[%1d]\n", mb->f_fixed[mb->nf][3]);
		}
		mb->f_theta[mb->nf][3] = log_prec_x;
		if (!mb->f_fixed[mb->nf][3]) {
			/*
			 * add this \theta 
			 */
			mb->theta = Realloc(mb->theta, mb->ntheta + 1, double **);
			mb->theta_hyperid = Realloc(mb->theta_hyperid, mb->ntheta + 1, char *);
			mb->theta_hyperid[mb->ntheta] = mb->f_prior[mb->nf][3].hyperid;
			mb->theta_tag = Realloc(mb->theta_tag, mb->ntheta + 1, char *);
			mb->theta_tag_userscale = Realloc(mb->theta_tag_userscale, mb->ntheta + 1, char *);
			mb->theta_dir = Realloc(mb->theta_dir, mb->ntheta + 1, char *);
			GMRFLib_sprintf(&msg, "MEC prec_x_intern for %s", (secname ? secname : mb->f_tag[mb->nf]));
			mb->theta_tag[mb->ntheta] = msg;
			GMRFLib_sprintf(&msg, "MEC prec_x for %s", (secname ? secname : mb->f_tag[mb->nf]));
			mb->theta_tag_userscale[mb->ntheta] = msg;
			GMRFLib_sprintf(&msg, "%s-parameter3", mb->f_dir[mb->nf]);
			mb->theta_dir[mb->ntheta] = msg;

			mb->theta_from = Realloc(mb->theta_from, mb->ntheta + 1, char *);
			mb->theta_to = Realloc(mb->theta_to, mb->ntheta + 1, char *);
			mb->theta_from[mb->ntheta] = Strdup(mb->f_prior[mb->nf][3].from_theta);
			mb->theta_to[mb->ntheta] = Strdup(mb->f_prior[mb->nf][3].to_theta);

			mb->theta[mb->ntheta] = log_prec_x;
			mb->theta_map = Realloc(mb->theta_map, mb->ntheta + 1, map_func_tp *);
			mb->theta_map[mb->ntheta] = map_precision;
			mb->theta_map_arg = Realloc(mb->theta_map_arg, mb->ntheta + 1, void *);
			mb->theta_map_arg[mb->ntheta] = NULL;
			mb->ntheta++;
		}
	}
		break;

	case F_MEB:
	{
		tmp = iniparser_getdouble(ini, inla_string_join(secname, "INITIAL0"), 1.0);
		if (!mb->f_fixed[mb->nf][0] && mb->mode_use_mode) {
			tmp = mb->theta_file[mb->theta_counter_file++];
			if (mb->mode_fixed) {
				mb->f_fixed[mb->nf][0] = 1;
			}
		}
		_SetInitial(0, tmp);
		HYPER_INIT(beta, tmp);

		double *range = NULL;
		range = Calloc(2, double);		       /* need this as it will be stored in the map argument */
		range[0] = iniparser_getdouble(ini, inla_string_join(secname, "RANGE.LOW"), 0.0);	/* low = high ==> map = identity */
		range[1] = iniparser_getdouble(ini, inla_string_join(secname, "RANGE.HIGH"), 0.0);

		if (mb->verbose) {
			printf("\t\tinitialise beta[%g]\n", tmp);
			printf("\t\tfixed=[%1d]\n", mb->f_fixed[mb->nf][0]);
			printf("\t\trange[%g, %g]\n", range[0], range[1]);
		}

		mb->f_theta_map[mb->nf] = Calloc(1, map_func_tp *);
		mb->f_theta_map_arg[mb->nf] = Calloc(1, void *);
		mb->f_theta_map[mb->nf][0] = map_beta;	       /* need these */
		mb->f_theta_map_arg[mb->nf][0] = (void *) range;	/* and this one as well */
		mb->f_theta[mb->nf][0] = beta;

		if (!mb->f_fixed[mb->nf][0]) {
			/*
			 * add this \theta 
			 */
			mb->theta = Realloc(mb->theta, mb->ntheta + 1, double **);
			mb->theta_hyperid = Realloc(mb->theta_hyperid, mb->ntheta + 1, char *);
			mb->theta_hyperid[mb->ntheta] = mb->f_prior[mb->nf][0].hyperid;
			mb->theta_tag = Realloc(mb->theta_tag, mb->ntheta + 1, char *);
			mb->theta_tag_userscale = Realloc(mb->theta_tag_userscale, mb->ntheta + 1, char *);
			mb->theta_dir = Realloc(mb->theta_dir, mb->ntheta + 1, char *);
			GMRFLib_sprintf(&msg, "MEB beta for %s", (secname ? secname : mb->f_tag[mb->nf]));
			mb->theta_tag[mb->ntheta] = msg;
			GMRFLib_sprintf(&msg, "MEB beta for %s", (secname ? secname : mb->f_tag[mb->nf]));
			mb->theta_tag_userscale[mb->ntheta] = msg;
			GMRFLib_sprintf(&msg, "%s-parameter0", mb->f_dir[mb->nf]);
			mb->theta_dir[mb->ntheta] = msg;

			mb->theta_from = Realloc(mb->theta_from, mb->ntheta + 1, char *);
			mb->theta_to = Realloc(mb->theta_to, mb->ntheta + 1, char *);
			mb->theta_from[mb->ntheta] = Strdup(mb->f_prior[mb->nf][0].from_theta);
			mb->theta_to[mb->ntheta] = Strdup(mb->f_prior[mb->nf][0].to_theta);

			mb->theta[mb->ntheta] = beta;
			mb->theta_map = Realloc(mb->theta_map, mb->ntheta + 1, map_func_tp *);
			mb->theta_map[mb->ntheta] = map_beta;
			mb->theta_map_arg = Realloc(mb->theta_map_arg, mb->ntheta + 1, void *);
			mb->theta_map_arg[mb->ntheta] = (void *) range;
			mb->ntheta++;
		}

		tmp = iniparser_getdouble(ini, inla_string_join(secname, "INITIAL1"), G.log_prec_initial);
		if (!mb->f_fixed[mb->nf][1] && mb->mode_use_mode) {
			tmp = mb->theta_file[mb->theta_counter_file++];
			if (mb->mode_fixed) {
				mb->f_fixed[mb->nf][1] = 1;
			}
		}
		_SetInitial(1, tmp);
		HYPER_INIT(log_prec, tmp);
		if (mb->verbose) {
			printf("\t\tinitialise log_prec_u[%g]\n", tmp);
			printf("\t\tfixed=[%1d]\n", mb->f_fixed[mb->nf][1]);
		}
		mb->f_theta[mb->nf][1] = log_prec;
		if (!mb->f_fixed[mb->nf][1]) {
			/*
			 * add this \theta 
			 */
			mb->theta = Realloc(mb->theta, mb->ntheta + 1, double **);
			mb->theta_hyperid = Realloc(mb->theta_hyperid, mb->ntheta + 1, char *);
			mb->theta_hyperid[mb->ntheta] = mb->f_prior[mb->nf][1].hyperid;
			mb->theta_tag = Realloc(mb->theta_tag, mb->ntheta + 1, char *);
			mb->theta_tag_userscale = Realloc(mb->theta_tag_userscale, mb->ntheta + 1, char *);
			mb->theta_dir = Realloc(mb->theta_dir, mb->ntheta + 1, char *);
			GMRFLib_sprintf(&msg, "MEB prec_u_intern for %s", (secname ? secname : mb->f_tag[mb->nf]));
			mb->theta_tag[mb->ntheta] = msg;
			GMRFLib_sprintf(&msg, "MEB prec_u for %s", (secname ? secname : mb->f_tag[mb->nf]));
			mb->theta_tag_userscale[mb->ntheta] = msg;
			GMRFLib_sprintf(&msg, "%s-parameter1", mb->f_dir[mb->nf]);
			mb->theta_dir[mb->ntheta] = msg;

			mb->theta_from = Realloc(mb->theta_from, mb->ntheta + 1, char *);
			mb->theta_to = Realloc(mb->theta_to, mb->ntheta + 1, char *);
			mb->theta_from[mb->ntheta] = Strdup(mb->f_prior[mb->nf][1].from_theta);
			mb->theta_to[mb->ntheta] = Strdup(mb->f_prior[mb->nf][1].to_theta);

			mb->theta[mb->ntheta] = log_prec;
			mb->theta_map = Realloc(mb->theta_map, mb->ntheta + 1, map_func_tp *);
			mb->theta_map[mb->ntheta] = map_precision;
			mb->theta_map_arg = Realloc(mb->theta_map_arg, mb->ntheta + 1, void *);
			mb->theta_map_arg[mb->ntheta] = NULL;
			mb->ntheta++;
		}
	}
		break;

	case F_R_GENERIC:
	{
		rgeneric_filename = iniparser_getstring(ini, inla_string_join(secname, "RGENERIC.FILE"), NULL);
		rgeneric_model = iniparser_getstring(ini, inla_string_join(secname, "RGENERIC.MODEL"), NULL);

		if (mb->verbose) {
			printf("\t\trgeneric.file   [%s]\n", rgeneric_filename);
			printf("\t\trgeneric.model  [%s]\n", rgeneric_model);
		}

		int n_out, nn_out, nn;
		double *x_out = NULL, *xx_out = NULL;

		/*
		 * we need to know ntheta, therefore we need to initialise and load files etc, here...
		 */

		int zero = 0;
		if (R_load_INLA) {
			inla_R_library("INLA");
			R_load_INLA = 0;
		}
		inla_R_load(rgeneric_filename);
		inla_R_rgeneric(&n_out, &x_out, R_GENERIC_INITIAL, rgeneric_model, &zero, NULL);
		inla_R_rgeneric(&nn_out, &xx_out, R_GENERIC_GRAPH, rgeneric_model, &zero, NULL);	/* need graph->n */

		nn = (int) xx_out[0];
		if (mb->f_n[mb->nf] != nn) {
			int err = 0;
			for (i = 0; i < mb->f_n[mb->nf]; i++) {
				// provide a warning if something could be wrong in the input
				if (mb->f_locations[mb->nf][i] != i + 1)
					err++;
			}
			if (err) {
				GMRFLib_sprintf(&msg, "%s\n\t\t%s, %1d != %1d, \n\t\t%s\n\t\t%s%1d%s",
						"There is a potential issue with the 'rgeneric' model and the indices used:",
						"the dimension of the model is different from expected", mb->f_n[mb->nf], nn,
						"and correctness cannot be verified.",
						"Please use argument n=", nn, ", f.ex, to *define* the dimension of the rgeneric model.");
				inla_error_general(msg);
				assert(0 != 1);
				exit(1);
			}

			Free(mb->f_locations[mb->nf]);
			mb->f_locations[mb->nf] = Calloc(nn, double);
			for (i = 0; i < nn; i++) {
				mb->f_locations[mb->nf][i] = i + 1;	/* set default ones */
			}
			mb->f_n[mb->nf] = mb->f_N[mb->nf] = nn;
		}

		int ntheta = 0;
		double *initial = NULL;

		ntheta = (int) x_out[0];
		if (ntheta) {
			initial = Calloc(ntheta, double);
			Memcpy(initial, &(x_out[1]), ntheta * sizeof(double));
		}

		Free(x_out);
		Free(xx_out);

		mb->f_ntheta[mb->nf] = ntheta;
		mb->f_initial[mb->nf] = initial;

		if (mb->verbose) {
			int ii;

			printf("\t\tntheta = [%1d]\n", ntheta);
			for (ii = 0; ii < ntheta; ii++) {
				printf("\t\tinitial[%1d] = %g\n", ii, initial[ii]);
			}
		}

		if (ntheta) {
			mb->f_fixed[mb->nf] = Calloc(ntheta, int);
			mb->f_theta[mb->nf] = Calloc(ntheta, double **);
		} else {
			mb->f_fixed[mb->nf] = NULL;
			mb->f_theta[mb->nf] = NULL;
		}

		for (i = 0; i < ntheta; i++) {
			double theta_initial;

			mb->f_fixed[mb->nf][i] = 0;
			theta_initial = initial[i];

			if (!mb->f_fixed[mb->nf][i] && mb->mode_use_mode) {
				theta_initial = mb->theta_file[mb->theta_counter_file++];
				if (mb->mode_fixed) {
					mb->f_fixed[mb->nf][i] = 1;
				}
			}
			HYPER_NEW(mb->f_theta[mb->nf][i], theta_initial);
			if (mb->verbose) {
				printf("\t\tinitialise theta[%1d]=[%g]\n", i, theta_initial);
				printf("\t\tfixed[%1d]=[%1d]\n", i, mb->f_fixed[mb->nf][i]);
			}

			/*
			 * add this \theta 
			 */
			mb->theta = Realloc(mb->theta, mb->ntheta + 1, double **);
			mb->theta_hyperid = Realloc(mb->theta_hyperid, mb->ntheta + 1, char *);
			mb->theta_hyperid[mb->ntheta] = NULL;
			mb->theta_tag = Realloc(mb->theta_tag, mb->ntheta + 1, char *);
			mb->theta_tag_userscale = Realloc(mb->theta_tag_userscale, mb->ntheta + 1, char *);
			mb->theta_dir = Realloc(mb->theta_dir, mb->ntheta + 1, char *);
			GMRFLib_sprintf(&msg, "Theta%1d for %s", i + 1, (secname ? secname : mb->f_tag[mb->nf]));
			mb->theta_tag[mb->ntheta] = msg;
			GMRFLib_sprintf(&msg, "Theta%1d for %s", i + 1, (secname ? secname : mb->f_tag[mb->nf]));
			mb->theta_tag_userscale[mb->ntheta] = msg;
			GMRFLib_sprintf(&msg, "%s-parameter%1d", mb->f_dir[mb->nf], i + 1);
			mb->theta_dir[mb->ntheta] = msg;

			mb->theta_from = Realloc(mb->theta_from, mb->ntheta + 1, char *);
			mb->theta_to = Realloc(mb->theta_to, mb->ntheta + 1, char *);
			mb->theta_from[mb->ntheta] = Strdup("function(x) x");
			mb->theta_to[mb->ntheta] = Strdup("function(x) x");

			mb->theta[mb->ntheta] = mb->f_theta[mb->nf][i];
			mb->theta_map = Realloc(mb->theta_map, mb->ntheta + 1, map_func_tp *);
			mb->theta_map[mb->ntheta] = map_identity;
			mb->theta_map_arg = Realloc(mb->theta_map_arg, mb->ntheta + 1, void *);
			mb->theta_map_arg[mb->ntheta] = NULL;
			mb->ntheta++;
		}

	}
		break;

	case F_C_GENERIC:
	{
		const char *emsg = NULL;
		char *cdata_fnm = NULL, *cgeneric_qfnm = NULL;
		int cgeneric_q;
		cgeneric_shlib = iniparser_getstring(ini, inla_string_join(secname, "CGENERIC.SHLIB"), NULL);
		cgeneric_model = iniparser_getstring(ini, inla_string_join(secname, "CGENERIC.MODEL"), NULL);
		cgeneric_n = iniparser_getint(ini, inla_string_join(secname, "CGENERIC.N"), cgeneric_n);
		cgeneric_debug = iniparser_getboolean(ini, inla_string_join(secname, "CGENERIC.DEBUG"), cgeneric_debug);
		cgeneric_q = iniparser_getboolean(ini, inla_string_join(secname, "CGENERIC.Q"), 0);
		cgeneric_qfnm = iniparser_getstring(ini, inla_string_join(secname, "CGENERIC.Q.FILE"), NULL);
		cdata_fnm = iniparser_getstring(ini, inla_string_join(secname, "CGENERIC.DATA"), NULL);

		cgeneric_data = inla_cgeneric_read_data(cdata_fnm, cgeneric_debug);

		if (mb->verbose) {
			printf("\t\tcgeneric.shlib  [%s]\n", cgeneric_shlib);
			printf("\t\tcgeneric.model  [%s]\n", cgeneric_model);
			printf("\t\tcgeneric.data   [%s]\n", cdata_fnm);
			printf("\t\tcgeneric.n      [%1d]\n", cgeneric_n);
			printf("\t\tcgeneric.q      [%1d]\n", cgeneric_q);
			printf("\t\tcgeneric.debug  [%1d]\n", cgeneric_debug);
		}

		int nn;
		double *x_out = NULL, *xx_out = NULL;

		static int ltdl_init = 1;
		if (ltdl_init) {
			lt_dlinit();
			if ((emsg = lt_dlerror())) {
				GMRFLib_sprintf(&msg, "\n *** dlinit error with model[%s] err_msg[%s]\n", cgeneric_model, emsg);
				inla_error_general(msg);
				assert(0 != 1);
				exit(1);
			}
			ltdl_init = 0;
			lt_dlerror();
		}

		handle = lt_dlopen(cgeneric_shlib);
		if (!handle) {
			GMRFLib_sprintf(&msg, "\n *** dlopen error with file[%s] err_msg[%s]\n", cgeneric_shlib, lt_dlerror());
			inla_error_general(msg);
			assert(0 != 1);
			exit(1);
		}
		lt_dlerror();

		model_func = (inla_cgeneric_func_tp *) lt_dlsym(handle, cgeneric_model);
		if ((emsg = lt_dlerror())) {
			GMRFLib_sprintf(&msg, "\n *** dlsym error with model[%s] err_msg[%s]\n", cgeneric_model, emsg);
			inla_error_general(msg);
			assert(0 != 1);
			exit(1);
		}
		lt_dlerror();

		if (cgeneric_q) {
			// this is a hack for `inla.cgeneric.q`. write out info and then exit.

			double *x = NULL, *theta = NULL;
			int nt;

			FILE *fp = fopen(cgeneric_qfnm, "w");

			fprintf(fp, "CGENERIC_BEGIN\n");
			x = model_func(INLA_CGENERIC_INITIAL, NULL, cgeneric_data);
			nt = (int) x[0];
			if (nt > 0) {
				theta = x + 1;
			} else {
				theta = NULL;
			}
			inla_cgeneric_debug(fp, secname, INLA_CGENERIC_INITIAL, x);

			x = model_func(INLA_CGENERIC_GRAPH, NULL, cgeneric_data);
			inla_cgeneric_debug(fp, secname, INLA_CGENERIC_GRAPH, x);

			x = model_func(INLA_CGENERIC_Q, theta, cgeneric_data);
			inla_cgeneric_debug(fp, secname, INLA_CGENERIC_Q, x);

			x = model_func(INLA_CGENERIC_MU, theta, cgeneric_data);
			inla_cgeneric_debug(fp, secname, INLA_CGENERIC_MU, x);

			x = model_func(INLA_CGENERIC_LOG_PRIOR, theta, cgeneric_data);
			inla_cgeneric_debug(fp, secname, INLA_CGENERIC_LOG_PRIOR, x);

			x = model_func(INLA_CGENERIC_LOG_NORM_CONST, theta, cgeneric_data);
			inla_cgeneric_debug(fp, secname, INLA_CGENERIC_LOG_NORM_CONST, x);

			x = model_func(INLA_CGENERIC_QUIT, theta, cgeneric_data);

			fprintf(fp, "CGENERIC_END\n");
			fclose(fp);

			// YES
			exit(0);
		}

		xx_out = model_func(INLA_CGENERIC_GRAPH, NULL, cgeneric_data);
		if (cgeneric_debug) {
			inla_cgeneric_debug(stdout, secname, INLA_CGENERIC_GRAPH, xx_out);
		}
		nn = (int) xx_out[0];

		if (nn != cgeneric_n) {
			GMRFLib_sprintf(&msg, "f(%s): 'n' does not match with the one in cgeneric-function '%s', %1d != %1d\n",
					secname, cgeneric_model, nn, cgeneric_n);
			inla_error_general(msg);
			assert(0 != 1);
			exit(1);
		}

		if (mb->f_n[mb->nf] != nn) {
			int err = 0;
			for (i = 0; i < mb->f_n[mb->nf]; i++) {
				// provide a warning if something could be wrong in the input
				if (mb->f_locations[mb->nf][i] != i + 1)
					err++;
			}
			if (err) {
				GMRFLib_sprintf(&msg, "%s\n\t\t%s, %1d != %1d, \n\t\t%s\n\t\t%s%1d%s",
						"There is a potential issue with the 'cgeneric' model and the indices used:",
						"the dimension of the model is different from expected", mb->f_n[mb->nf], nn,
						"and correctness cannot be verified.",
						"Please use argument n=", nn, ", f.ex, to *define* the dimension of the cgeneric model.");
				inla_error_general(msg);
				assert(0 != 1);
				exit(1);
			}

			Free(mb->f_locations[mb->nf]);
			mb->f_locations[mb->nf] = Calloc(nn, double);
			for (i = 0; i < nn; i++) {
				mb->f_locations[mb->nf][i] = i + 1;	/* set default ones */
			}
			mb->f_n[mb->nf] = mb->f_N[mb->nf] = nn;
		}
		Free(xx_out);

		int ntheta = 0;
		double *initial = NULL;

		x_out = model_func(INLA_CGENERIC_INITIAL, NULL, cgeneric_data);
		if (cgeneric_debug) {
			inla_cgeneric_debug(stdout, secname, INLA_CGENERIC_INITIAL, x_out);
		}
		ntheta = (int) x_out[0];
		if (ntheta) {
			initial = Calloc(ntheta, double);
			Memcpy(initial, &(x_out[1]), ntheta * sizeof(double));
		}
		Free(x_out);

		mb->f_ntheta[mb->nf] = ntheta;
		mb->f_initial[mb->nf] = initial;

		if (mb->verbose) {
			printf("\t\tntheta = [%1d]\n", ntheta);
			for (int ii = 0; ii < ntheta; ii++) {
				printf("\t\tinitial[%1d] = %g\n", ii, initial[ii]);
			}
		}

		if (ntheta) {
			mb->f_fixed[mb->nf] = Calloc(ntheta, int);
			mb->f_theta[mb->nf] = Calloc(ntheta, double **);
		} else {
			mb->f_fixed[mb->nf] = NULL;
			mb->f_theta[mb->nf] = NULL;
		}

		for (i = 0; i < ntheta; i++) {
			double theta_initial;

			mb->f_fixed[mb->nf][i] = 0;
			theta_initial = initial[i];

			if (!mb->f_fixed[mb->nf][i] && mb->mode_use_mode) {
				theta_initial = mb->theta_file[mb->theta_counter_file++];
				if (mb->mode_fixed) {
					mb->f_fixed[mb->nf][i] = 1;
				}
			}
			HYPER_NEW(mb->f_theta[mb->nf][i], theta_initial);
			if (mb->verbose) {
				printf("\t\tinitialise theta[%1d]=[%g]\n", i, theta_initial);
				printf("\t\tfixed[%1d]=[%1d]\n", i, mb->f_fixed[mb->nf][i]);
			}

			/*
			 * add this \theta 
			 */
			mb->theta = Realloc(mb->theta, mb->ntheta + 1, double **);
			mb->theta_hyperid = Realloc(mb->theta_hyperid, mb->ntheta + 1, char *);
			mb->theta_hyperid[mb->ntheta] = NULL;
			mb->theta_tag = Realloc(mb->theta_tag, mb->ntheta + 1, char *);
			mb->theta_tag_userscale = Realloc(mb->theta_tag_userscale, mb->ntheta + 1, char *);
			mb->theta_dir = Realloc(mb->theta_dir, mb->ntheta + 1, char *);
			GMRFLib_sprintf(&msg, "Theta%1d for %s", i + 1, (secname ? secname : mb->f_tag[mb->nf]));
			mb->theta_tag[mb->ntheta] = msg;
			GMRFLib_sprintf(&msg, "Theta%1d for %s", i + 1, (secname ? secname : mb->f_tag[mb->nf]));
			mb->theta_tag_userscale[mb->ntheta] = msg;
			GMRFLib_sprintf(&msg, "%s-parameter%1d", mb->f_dir[mb->nf], i + 1);
			mb->theta_dir[mb->ntheta] = msg;

			mb->theta_from = Realloc(mb->theta_from, mb->ntheta + 1, char *);
			mb->theta_to = Realloc(mb->theta_to, mb->ntheta + 1, char *);
			mb->theta_from[mb->ntheta] = Strdup("function(x) x");
			mb->theta_to[mb->ntheta] = Strdup("function(x) x");

			mb->theta[mb->ntheta] = mb->f_theta[mb->nf][i];
			mb->theta_map = Realloc(mb->theta_map, mb->ntheta + 1, map_func_tp *);
			mb->theta_map[mb->ntheta] = map_identity;
			mb->theta_map_arg = Realloc(mb->theta_map_arg, mb->ntheta + 1, void *);
			mb->theta_map_arg[mb->ntheta] = NULL;
			mb->ntheta++;
		}

	}
		break;

	case F_FGN:
	case F_FGN2:
	{
		tmp = iniparser_getdouble(ini, inla_string_join(secname, "INITIAL0"), G.log_prec_initial);
		if (!mb->f_fixed[mb->nf][0] && mb->mode_use_mode) {
			tmp = mb->theta_file[mb->theta_counter_file++];
			if (mb->mode_fixed) {
				mb->f_fixed[mb->nf][0] = 1;
			}
		}
		_SetInitial(0, tmp);
		HYPER_INIT(log_prec, tmp);
		if (mb->verbose) {
			printf("\t\tinitialise log_precision[%g]\n", tmp);
			printf("\t\tfixed=[%1d]\n", mb->f_fixed[mb->nf][0]);
		}

		mb->f_theta[mb->nf] = Calloc(2, double **);
		mb->f_theta[mb->nf][0] = log_prec;
		if (!mb->f_fixed[mb->nf][0]) {
			/*
			 * add this \theta 
			 */
			mb->theta = Realloc(mb->theta, mb->ntheta + 1, double **);
			mb->theta_hyperid = Realloc(mb->theta_hyperid, mb->ntheta + 1, char *);
			mb->theta_hyperid[mb->ntheta] = mb->f_prior[mb->nf][0].hyperid;
			mb->theta_tag = Realloc(mb->theta_tag, mb->ntheta + 1, char *);
			mb->theta_tag_userscale = Realloc(mb->theta_tag_userscale, mb->ntheta + 1, char *);
			mb->theta_dir = Realloc(mb->theta_dir, mb->ntheta + 1, char *);
			GMRFLib_sprintf(&msg, "Log precision for %s", (secname ? secname : mb->f_tag[mb->nf]));
			mb->theta_tag[mb->ntheta] = msg;
			GMRFLib_sprintf(&msg, "Precision for %s", (secname ? secname : mb->f_tag[mb->nf]));
			mb->theta_tag_userscale[mb->ntheta] = msg;
			GMRFLib_sprintf(&msg, "%s-parameter0", mb->f_dir[mb->nf]);
			mb->theta_dir[mb->ntheta] = msg;

			mb->theta_from = Realloc(mb->theta_from, mb->ntheta + 1, char *);
			mb->theta_to = Realloc(mb->theta_to, mb->ntheta + 1, char *);
			mb->theta_from[mb->ntheta] = Strdup(mb->f_prior[mb->nf][0].from_theta);
			mb->theta_to[mb->ntheta] = Strdup(mb->f_prior[mb->nf][0].to_theta);

			mb->theta[mb->ntheta] = log_prec;
			mb->theta_map = Realloc(mb->theta_map, mb->ntheta + 1, map_func_tp *);
			mb->theta_map[mb->ntheta] = map_precision;
			mb->theta_map_arg = Realloc(mb->theta_map_arg, mb->ntheta + 1, void *);
			mb->theta_map_arg[mb->ntheta] = NULL;
			mb->ntheta++;
		}

		tmp = iniparser_getdouble(ini, inla_string_join(secname, "INITIAL1"), 0.0);
		if (!mb->f_fixed[mb->nf][1] && mb->mode_use_mode) {
			tmp = mb->theta_file[mb->theta_counter_file++];
			if (mb->mode_fixed) {
				mb->f_fixed[mb->nf][1] = 1;
			}
		}
		_SetInitial(1, tmp);
		HYPER_INIT(H_intern, tmp);
		if (mb->verbose) {
			printf("\t\tinitialise H_intern[%g]\n", tmp);
			printf("\t\tfixed=[%1d]\n", mb->f_fixed[mb->nf][1]);
		}
		mb->f_theta[mb->nf][1] = H_intern;
		if (!mb->f_fixed[mb->nf][1]) {
			/*
			 * add this \theta 
			 */
			mb->theta = Realloc(mb->theta, mb->ntheta + 1, double **);
			mb->theta_hyperid = Realloc(mb->theta_hyperid, mb->ntheta + 1, char *);
			mb->theta_hyperid[mb->ntheta] = mb->f_prior[mb->nf][1].hyperid;
			mb->theta_tag = Realloc(mb->theta_tag, mb->ntheta + 1, char *);
			mb->theta_tag_userscale = Realloc(mb->theta_tag_userscale, mb->ntheta + 1, char *);
			mb->theta_dir = Realloc(mb->theta_dir, mb->ntheta + 1, char *);
			GMRFLib_sprintf(&msg, "H_intern for %s", (secname ? secname : mb->f_tag[mb->nf]));
			mb->theta_tag[mb->ntheta] = msg;
			GMRFLib_sprintf(&msg, "H for %s", (secname ? secname : mb->f_tag[mb->nf]));
			mb->theta_tag_userscale[mb->ntheta] = msg;
			GMRFLib_sprintf(&msg, "%s-parameter1", mb->f_dir[mb->nf]);
			mb->theta_dir[mb->ntheta] = msg;

			mb->theta_from = Realloc(mb->theta_from, mb->ntheta + 1, char *);
			mb->theta_to = Realloc(mb->theta_to, mb->ntheta + 1, char *);
			mb->theta_from[mb->ntheta] = Strdup(mb->f_prior[mb->nf][1].from_theta);
			mb->theta_to[mb->ntheta] = Strdup(mb->f_prior[mb->nf][1].to_theta);

			mb->theta[mb->ntheta] = H_intern;
			mb->theta_map = Realloc(mb->theta_map, mb->ntheta + 1, map_func_tp *);
			mb->theta_map[mb->ntheta] = map_H;
			mb->theta_map_arg = Realloc(mb->theta_map_arg, mb->ntheta + 1, void *);
			mb->theta_map_arg[mb->ntheta] = NULL;
			mb->ntheta++;
		}
	}
		break;

	case F_AR1:
	{
		tmp = iniparser_getdouble(ini, inla_string_join(secname, "INITIAL0"), G.log_prec_initial);
		if (!mb->f_fixed[mb->nf][0] && mb->mode_use_mode) {
			tmp = mb->theta_file[mb->theta_counter_file++];
			if (mb->mode_fixed) {
				mb->f_fixed[mb->nf][0] = 1;
			}
		}
		_SetInitial(0, tmp);
		HYPER_INIT(log_prec, tmp);
		if (mb->verbose) {
			printf("\t\tinitialise log_precision[%g]\n", tmp);
			printf("\t\tfixed=[%1d]\n", mb->f_fixed[mb->nf][0]);
		}

		mb->f_theta[mb->nf] = Calloc(3, double **);
		mb->f_theta[mb->nf][0] = log_prec;
		if (!mb->f_fixed[mb->nf][0]) {
			/*
			 * add this \theta 
			 */
			mb->theta = Realloc(mb->theta, mb->ntheta + 1, double **);
			mb->theta_hyperid = Realloc(mb->theta_hyperid, mb->ntheta + 1, char *);
			mb->theta_hyperid[mb->ntheta] = mb->f_prior[mb->nf][0].hyperid;
			mb->theta_tag = Realloc(mb->theta_tag, mb->ntheta + 1, char *);
			mb->theta_tag_userscale = Realloc(mb->theta_tag_userscale, mb->ntheta + 1, char *);
			mb->theta_dir = Realloc(mb->theta_dir, mb->ntheta + 1, char *);
			GMRFLib_sprintf(&msg, "Log precision for %s", (secname ? secname : mb->f_tag[mb->nf]));
			mb->theta_tag[mb->ntheta] = msg;
			GMRFLib_sprintf(&msg, "Precision for %s", (secname ? secname : mb->f_tag[mb->nf]));
			mb->theta_tag_userscale[mb->ntheta] = msg;
			GMRFLib_sprintf(&msg, "%s-parameter0", mb->f_dir[mb->nf]);
			mb->theta_dir[mb->ntheta] = msg;

			mb->theta_from = Realloc(mb->theta_from, mb->ntheta + 1, char *);
			mb->theta_to = Realloc(mb->theta_to, mb->ntheta + 1, char *);
			mb->theta_from[mb->ntheta] = Strdup(mb->f_prior[mb->nf][0].from_theta);
			mb->theta_to[mb->ntheta] = Strdup(mb->f_prior[mb->nf][0].to_theta);

			mb->theta[mb->ntheta] = log_prec;
			mb->theta_map = Realloc(mb->theta_map, mb->ntheta + 1, map_func_tp *);
			mb->theta_map[mb->ntheta] = map_precision;
			mb->theta_map_arg = Realloc(mb->theta_map_arg, mb->ntheta + 1, void *);
			mb->theta_map_arg[mb->ntheta] = NULL;
			mb->ntheta++;
		}

		tmp = iniparser_getdouble(ini, inla_string_join(secname, "INITIAL1"), 0.0);
		if (!mb->f_fixed[mb->nf][1] && mb->mode_use_mode) {
			tmp = mb->theta_file[mb->theta_counter_file++];
			if (mb->mode_fixed) {
				mb->f_fixed[mb->nf][1] = 1;
			}
		}
		_SetInitial(1, tmp);
		HYPER_INIT(phi_intern, tmp);
		if (mb->verbose) {
			printf("\t\tinitialise phi_intern[%g]\n", tmp);
			printf("\t\tfixed=[%1d]\n", mb->f_fixed[mb->nf][1]);
		}
		mb->f_theta[mb->nf][1] = phi_intern;
		if (!mb->f_fixed[mb->nf][1]) {
			/*
			 * add this \theta 
			 */
			mb->theta = Realloc(mb->theta, mb->ntheta + 1, double **);
			mb->theta_hyperid = Realloc(mb->theta_hyperid, mb->ntheta + 1, char *);
			mb->theta_hyperid[mb->ntheta] = mb->f_prior[mb->nf][1].hyperid;
			mb->theta_tag = Realloc(mb->theta_tag, mb->ntheta + 1, char *);
			mb->theta_tag_userscale = Realloc(mb->theta_tag_userscale, mb->ntheta + 1, char *);
			mb->theta_dir = Realloc(mb->theta_dir, mb->ntheta + 1, char *);
			GMRFLib_sprintf(&msg, "Rho_intern for %s", (secname ? secname : mb->f_tag[mb->nf]));
			mb->theta_tag[mb->ntheta] = msg;
			GMRFLib_sprintf(&msg, "Rho for %s", (secname ? secname : mb->f_tag[mb->nf]));
			mb->theta_tag_userscale[mb->ntheta] = msg;
			GMRFLib_sprintf(&msg, "%s-parameter1", mb->f_dir[mb->nf]);
			mb->theta_dir[mb->ntheta] = msg;

			mb->theta_from = Realloc(mb->theta_from, mb->ntheta + 1, char *);
			mb->theta_to = Realloc(mb->theta_to, mb->ntheta + 1, char *);
			mb->theta_from[mb->ntheta] = Strdup(mb->f_prior[mb->nf][1].from_theta);
			mb->theta_to[mb->ntheta] = Strdup(mb->f_prior[mb->nf][1].to_theta);

			mb->theta[mb->ntheta] = phi_intern;
			mb->theta_map = Realloc(mb->theta_map, mb->ntheta + 1, map_func_tp *);
			mb->theta_map[mb->ntheta] = map_rho;
			mb->theta_map_arg = Realloc(mb->theta_map_arg, mb->ntheta + 1, void *);
			mb->theta_map_arg[mb->ntheta] = NULL;
			mb->ntheta++;
		}

		tmp = iniparser_getdouble(ini, inla_string_join(secname, "INITIAL2"), 0.0);
		if (!mb->f_fixed[mb->nf][2] && mb->mode_use_mode) {
			tmp = mb->theta_file[mb->theta_counter_file++];
			if (mb->mode_fixed) {
				mb->f_fixed[mb->nf][2] = 1;
			}
		}
		_SetInitial(2, tmp);
		HYPER_INIT(mean_x, tmp);
		if (mb->verbose) {
			printf("\t\tinitialise mean[%g]\n", tmp);
			printf("\t\tfixed=[%1d]\n", mb->f_fixed[mb->nf][2]);
		}
		mb->f_theta[mb->nf][2] = mean_x;
		if (!mb->f_fixed[mb->nf][2]) {
			/*
			 * add this \theta 
			 */
			mb->theta = Realloc(mb->theta, mb->ntheta + 1, double **);
			mb->theta_hyperid = Realloc(mb->theta_hyperid, mb->ntheta + 1, char *);
			mb->theta_hyperid[mb->ntheta] = mb->f_prior[mb->nf][1].hyperid;
			mb->theta_tag = Realloc(mb->theta_tag, mb->ntheta + 1, char *);
			mb->theta_tag_userscale = Realloc(mb->theta_tag_userscale, mb->ntheta + 1, char *);
			mb->theta_dir = Realloc(mb->theta_dir, mb->ntheta + 1, char *);
			GMRFLib_sprintf(&msg, "Mean for %s", (secname ? secname : mb->f_tag[mb->nf]));
			mb->theta_tag[mb->ntheta] = msg;
			GMRFLib_sprintf(&msg, "Mean for %s", (secname ? secname : mb->f_tag[mb->nf]));
			mb->theta_tag_userscale[mb->ntheta] = msg;
			GMRFLib_sprintf(&msg, "%s-parameter1", mb->f_dir[mb->nf]);
			mb->theta_dir[mb->ntheta] = msg;

			mb->theta_from = Realloc(mb->theta_from, mb->ntheta + 1, char *);
			mb->theta_to = Realloc(mb->theta_to, mb->ntheta + 1, char *);
			mb->theta_from[mb->ntheta] = Strdup(mb->f_prior[mb->nf][1].from_theta);
			mb->theta_to[mb->ntheta] = Strdup(mb->f_prior[mb->nf][1].to_theta);

			mb->theta[mb->ntheta] = mean_x;
			mb->theta_map = Realloc(mb->theta_map, mb->ntheta + 1, map_func_tp *);
			mb->theta_map[mb->ntheta] = map_identity;
			mb->theta_map_arg = Realloc(mb->theta_map_arg, mb->ntheta + 1, void *);
			mb->theta_map_arg[mb->ntheta] = NULL;
			mb->ntheta++;
		}

	}
		break;

	case F_AR1C:
	{
		tmp = iniparser_getdouble(ini, inla_string_join(secname, "INITIAL0"), G.log_prec_initial);
		if (!mb->f_fixed[mb->nf][0] && mb->mode_use_mode) {
			tmp = mb->theta_file[mb->theta_counter_file++];
			if (mb->mode_fixed) {
				mb->f_fixed[mb->nf][0] = 1;
			}
		}
		_SetInitial(0, tmp);
		HYPER_INIT(log_prec, tmp);
		if (mb->verbose) {
			printf("\t\tinitialise log_precision[%g]\n", tmp);
			printf("\t\tfixed=[%1d]\n", mb->f_fixed[mb->nf][0]);
		}

		mb->f_theta[mb->nf] = Calloc(2, double **);
		mb->f_theta[mb->nf][0] = log_prec;
		if (!mb->f_fixed[mb->nf][0]) {
			/*
			 * add this \theta 
			 */
			mb->theta = Realloc(mb->theta, mb->ntheta + 1, double **);
			mb->theta_hyperid = Realloc(mb->theta_hyperid, mb->ntheta + 1, char *);
			mb->theta_hyperid[mb->ntheta] = mb->f_prior[mb->nf][0].hyperid;
			mb->theta_tag = Realloc(mb->theta_tag, mb->ntheta + 1, char *);
			mb->theta_tag_userscale = Realloc(mb->theta_tag_userscale, mb->ntheta + 1, char *);
			mb->theta_dir = Realloc(mb->theta_dir, mb->ntheta + 1, char *);
			GMRFLib_sprintf(&msg, "Log precision for %s", (secname ? secname : mb->f_tag[mb->nf]));
			mb->theta_tag[mb->ntheta] = msg;
			GMRFLib_sprintf(&msg, "Precision for %s", (secname ? secname : mb->f_tag[mb->nf]));
			mb->theta_tag_userscale[mb->ntheta] = msg;
			GMRFLib_sprintf(&msg, "%s-parameter0", mb->f_dir[mb->nf]);
			mb->theta_dir[mb->ntheta] = msg;

			mb->theta_from = Realloc(mb->theta_from, mb->ntheta + 1, char *);
			mb->theta_to = Realloc(mb->theta_to, mb->ntheta + 1, char *);
			mb->theta_from[mb->ntheta] = Strdup(mb->f_prior[mb->nf][0].from_theta);
			mb->theta_to[mb->ntheta] = Strdup(mb->f_prior[mb->nf][0].to_theta);

			mb->theta[mb->ntheta] = log_prec;
			mb->theta_map = Realloc(mb->theta_map, mb->ntheta + 1, map_func_tp *);
			mb->theta_map[mb->ntheta] = map_precision;
			mb->theta_map_arg = Realloc(mb->theta_map_arg, mb->ntheta + 1, void *);
			mb->theta_map_arg[mb->ntheta] = NULL;
			mb->ntheta++;
		}

		tmp = iniparser_getdouble(ini, inla_string_join(secname, "INITIAL1"), 0.0);
		if (!mb->f_fixed[mb->nf][1] && mb->mode_use_mode) {
			tmp = mb->theta_file[mb->theta_counter_file++];
			if (mb->mode_fixed) {
				mb->f_fixed[mb->nf][1] = 1;
			}
		}
		_SetInitial(1, tmp);
		HYPER_INIT(phi_intern, tmp);
		if (mb->verbose) {
			printf("\t\tinitialise phi_intern[%g]\n", tmp);
			printf("\t\tfixed=[%1d]\n", mb->f_fixed[mb->nf][1]);
		}
		mb->f_theta[mb->nf][1] = phi_intern;
		if (!mb->f_fixed[mb->nf][1]) {
			/*
			 * add this \theta 
			 */
			mb->theta = Realloc(mb->theta, mb->ntheta + 1, double **);
			mb->theta_hyperid = Realloc(mb->theta_hyperid, mb->ntheta + 1, char *);
			mb->theta_hyperid[mb->ntheta] = mb->f_prior[mb->nf][1].hyperid;
			mb->theta_tag = Realloc(mb->theta_tag, mb->ntheta + 1, char *);
			mb->theta_tag_userscale = Realloc(mb->theta_tag_userscale, mb->ntheta + 1, char *);
			mb->theta_dir = Realloc(mb->theta_dir, mb->ntheta + 1, char *);
			GMRFLib_sprintf(&msg, "Rho_intern for %s", (secname ? secname : mb->f_tag[mb->nf]));
			mb->theta_tag[mb->ntheta] = msg;
			GMRFLib_sprintf(&msg, "Rho for %s", (secname ? secname : mb->f_tag[mb->nf]));
			mb->theta_tag_userscale[mb->ntheta] = msg;
			GMRFLib_sprintf(&msg, "%s-parameter1", mb->f_dir[mb->nf]);
			mb->theta_dir[mb->ntheta] = msg;

			mb->theta_from = Realloc(mb->theta_from, mb->ntheta + 1, char *);
			mb->theta_to = Realloc(mb->theta_to, mb->ntheta + 1, char *);
			mb->theta_from[mb->ntheta] = Strdup(mb->f_prior[mb->nf][1].from_theta);
			mb->theta_to[mb->ntheta] = Strdup(mb->f_prior[mb->nf][1].to_theta);

			mb->theta[mb->ntheta] = phi_intern;
			mb->theta_map = Realloc(mb->theta_map, mb->ntheta + 1, map_func_tp *);
			mb->theta_map[mb->ntheta] = map_rho;
			mb->theta_map_arg = Realloc(mb->theta_map_arg, mb->ntheta + 1, void *);
			mb->theta_map_arg[mb->ntheta] = NULL;
			mb->ntheta++;
		}
	}
		break;

	case F_SLM:
	{
		tmp = iniparser_getdouble(ini, inla_string_join(secname, "INITIAL0"), G.log_prec_initial);
		if (!mb->f_fixed[mb->nf][0] && mb->mode_use_mode) {
			tmp = mb->theta_file[mb->theta_counter_file++];
			if (mb->mode_fixed) {
				mb->f_fixed[mb->nf][0] = 1;
			}
		}
		_SetInitial(0, tmp);
		HYPER_INIT(log_prec, tmp);
		if (mb->verbose) {
			printf("\t\tinitialise log_precision[%g]\n", tmp);
			printf("\t\tfixed=[%1d]\n", mb->f_fixed[mb->nf][0]);
		}

		mb->f_theta[mb->nf] = Calloc(2, double **);
		mb->f_theta[mb->nf][0] = log_prec;
		if (!mb->f_fixed[mb->nf][0]) {
			/*
			 * add this \theta 
			 */
			mb->theta = Realloc(mb->theta, mb->ntheta + 1, double **);
			mb->theta_hyperid = Realloc(mb->theta_hyperid, mb->ntheta + 1, char *);
			mb->theta_hyperid[mb->ntheta] = mb->f_prior[mb->nf][0].hyperid;
			mb->theta_tag = Realloc(mb->theta_tag, mb->ntheta + 1, char *);
			mb->theta_tag_userscale = Realloc(mb->theta_tag_userscale, mb->ntheta + 1, char *);
			mb->theta_dir = Realloc(mb->theta_dir, mb->ntheta + 1, char *);
			GMRFLib_sprintf(&msg, "Log precision for %s", (secname ? secname : mb->f_tag[mb->nf]));
			mb->theta_tag[mb->ntheta] = msg;
			GMRFLib_sprintf(&msg, "Precision for %s", (secname ? secname : mb->f_tag[mb->nf]));
			mb->theta_tag_userscale[mb->ntheta] = msg;
			GMRFLib_sprintf(&msg, "%s-parameter0", mb->f_dir[mb->nf]);
			mb->theta_dir[mb->ntheta] = msg;

			mb->theta_from = Realloc(mb->theta_from, mb->ntheta + 1, char *);
			mb->theta_to = Realloc(mb->theta_to, mb->ntheta + 1, char *);
			mb->theta_from[mb->ntheta] = Strdup(mb->f_prior[mb->nf][0].from_theta);
			mb->theta_to[mb->ntheta] = Strdup(mb->f_prior[mb->nf][0].to_theta);

			mb->theta[mb->ntheta] = log_prec;
			mb->theta_map = Realloc(mb->theta_map, mb->ntheta + 1, map_func_tp *);
			mb->theta_map[mb->ntheta] = map_precision;
			mb->theta_map_arg = Realloc(mb->theta_map_arg, mb->ntheta + 1, void *);
			mb->theta_map_arg[mb->ntheta] = NULL;
			mb->ntheta++;
		}

		tmp = iniparser_getdouble(ini, inla_string_join(secname, "INITIAL1"), 0.0);
		if (!mb->f_fixed[mb->nf][1] && mb->mode_use_mode) {
			tmp = mb->theta_file[mb->theta_counter_file++];
			if (mb->mode_fixed) {
				mb->f_fixed[mb->nf][1] = 1;
			}
		}
		_SetInitial(1, tmp);
		HYPER_INIT(rho_intern, tmp);
		if (mb->verbose) {
			printf("\t\tinitialise rho_intern[%g]\n", tmp);
			printf("\t\tfixed=[%1d]\n", mb->f_fixed[mb->nf][1]);
		}
		mb->f_theta[mb->nf][1] = rho_intern;
		if (!mb->f_fixed[mb->nf][1]) {
			/*
			 * add this \theta 
			 */
			mb->theta = Realloc(mb->theta, mb->ntheta + 1, double **);
			mb->theta_hyperid = Realloc(mb->theta_hyperid, mb->ntheta + 1, char *);
			mb->theta_hyperid[mb->ntheta] = mb->f_prior[mb->nf][1].hyperid;
			mb->theta_tag = Realloc(mb->theta_tag, mb->ntheta + 1, char *);
			mb->theta_tag_userscale = Realloc(mb->theta_tag_userscale, mb->ntheta + 1, char *);
			mb->theta_dir = Realloc(mb->theta_dir, mb->ntheta + 1, char *);
			GMRFLib_sprintf(&msg, "Rho_intern for %s", (secname ? secname : mb->f_tag[mb->nf]));
			mb->theta_tag[mb->ntheta] = msg;
			GMRFLib_sprintf(&msg, "Rho for %s", (secname ? secname : mb->f_tag[mb->nf]));
			mb->theta_tag_userscale[mb->ntheta] = msg;
			GMRFLib_sprintf(&msg, "%s-parameter1", mb->f_dir[mb->nf]);
			mb->theta_dir[mb->ntheta] = msg;

			mb->theta_from = Realloc(mb->theta_from, mb->ntheta + 1, char *);
			mb->theta_to = Realloc(mb->theta_to, mb->ntheta + 1, char *);
			mb->theta_from[mb->ntheta] = Strdup(mb->f_prior[mb->nf][1].from_theta);
			mb->theta_to[mb->ntheta] = Strdup(mb->f_prior[mb->nf][1].to_theta);

			mb->theta[mb->ntheta] = rho_intern;
			mb->theta_map = Realloc(mb->theta_map, mb->ntheta + 1, map_func_tp *);
			mb->theta_map[mb->ntheta] = map_probability;
			mb->theta_map_arg = Realloc(mb->theta_map_arg, mb->ntheta + 1, void *);
			mb->theta_map_arg[mb->ntheta] = NULL;
			mb->ntheta++;
		}
	}
		break;

	case F_OU:
	{
		tmp = iniparser_getdouble(ini, inla_string_join(secname, "INITIAL0"), G.log_prec_initial);
		if (!mb->f_fixed[mb->nf][0] && mb->mode_use_mode) {
			tmp = mb->theta_file[mb->theta_counter_file++];
			if (mb->mode_fixed) {
				mb->f_fixed[mb->nf][0] = 1;
			}
		}
		_SetInitial(0, tmp);
		HYPER_INIT(log_prec, tmp);
		if (mb->verbose) {
			printf("\t\tinitialise log_precision[%g]\n", tmp);
			printf("\t\tfixed=[%1d]\n", mb->f_fixed[mb->nf][0]);
		}

		mb->f_theta[mb->nf] = Calloc(2, double **);
		mb->f_theta[mb->nf][0] = log_prec;
		if (!mb->f_fixed[mb->nf][0]) {
			/*
			 * add this \theta 
			 */
			mb->theta = Realloc(mb->theta, mb->ntheta + 1, double **);
			mb->theta_hyperid = Realloc(mb->theta_hyperid, mb->ntheta + 1, char *);
			mb->theta_hyperid[mb->ntheta] = mb->f_prior[mb->nf][0].hyperid;
			mb->theta_tag = Realloc(mb->theta_tag, mb->ntheta + 1, char *);
			mb->theta_tag_userscale = Realloc(mb->theta_tag_userscale, mb->ntheta + 1, char *);
			mb->theta_dir = Realloc(mb->theta_dir, mb->ntheta + 1, char *);
			GMRFLib_sprintf(&msg, "Log precision for %s", (secname ? secname : mb->f_tag[mb->nf]));
			mb->theta_tag[mb->ntheta] = msg;
			GMRFLib_sprintf(&msg, "Precision for %s", (secname ? secname : mb->f_tag[mb->nf]));
			mb->theta_tag_userscale[mb->ntheta] = msg;
			GMRFLib_sprintf(&msg, "%s-parameter0", mb->f_dir[mb->nf]);
			mb->theta_dir[mb->ntheta] = msg;

			mb->theta_from = Realloc(mb->theta_from, mb->ntheta + 1, char *);
			mb->theta_to = Realloc(mb->theta_to, mb->ntheta + 1, char *);
			mb->theta_from[mb->ntheta] = Strdup(mb->f_prior[mb->nf][0].from_theta);
			mb->theta_to[mb->ntheta] = Strdup(mb->f_prior[mb->nf][0].to_theta);

			mb->theta[mb->ntheta] = log_prec;
			mb->theta_map = Realloc(mb->theta_map, mb->ntheta + 1, map_func_tp *);
			mb->theta_map[mb->ntheta] = map_precision;
			mb->theta_map_arg = Realloc(mb->theta_map_arg, mb->ntheta + 1, void *);
			mb->theta_map_arg[mb->ntheta] = NULL;
			mb->ntheta++;
		}

		tmp = iniparser_getdouble(ini, inla_string_join(secname, "INITIAL1"), 0.0);
		if (!mb->f_fixed[mb->nf][1] && mb->mode_use_mode) {
			tmp = mb->theta_file[mb->theta_counter_file++];
			if (mb->mode_fixed) {
				mb->f_fixed[mb->nf][1] = 1;
			}
		}
		_SetInitial(1, tmp);
		HYPER_INIT(phi_intern, tmp);
		if (mb->verbose) {
			printf("\t\tinitialise phi_intern[%g]\n", tmp);
			printf("\t\tfixed=[%1d]\n", mb->f_fixed[mb->nf][1]);
		}
		mb->f_theta[mb->nf][1] = phi_intern;
		if (!mb->f_fixed[mb->nf][1]) {
			/*
			 * add this \theta 
			 */
			mb->theta = Realloc(mb->theta, mb->ntheta + 1, double **);
			mb->theta_hyperid = Realloc(mb->theta_hyperid, mb->ntheta + 1, char *);
			mb->theta_hyperid[mb->ntheta] = mb->f_prior[mb->nf][1].hyperid;
			mb->theta_tag = Realloc(mb->theta_tag, mb->ntheta + 1, char *);
			mb->theta_tag_userscale = Realloc(mb->theta_tag_userscale, mb->ntheta + 1, char *);
			mb->theta_dir = Realloc(mb->theta_dir, mb->ntheta + 1, char *);
			GMRFLib_sprintf(&msg, "Phi_intern for %s", (secname ? secname : mb->f_tag[mb->nf]));
			mb->theta_tag[mb->ntheta] = msg;
			GMRFLib_sprintf(&msg, "Phi for %s", (secname ? secname : mb->f_tag[mb->nf]));
			mb->theta_tag_userscale[mb->ntheta] = msg;
			GMRFLib_sprintf(&msg, "%s-parameter1", mb->f_dir[mb->nf]);
			mb->theta_dir[mb->ntheta] = msg;

			mb->theta_from = Realloc(mb->theta_from, mb->ntheta + 1, char *);
			mb->theta_to = Realloc(mb->theta_to, mb->ntheta + 1, char *);
			mb->theta_from[mb->ntheta] = Strdup(mb->f_prior[mb->nf][1].from_theta);
			mb->theta_to[mb->ntheta] = Strdup(mb->f_prior[mb->nf][1].to_theta);

			mb->theta[mb->ntheta] = phi_intern;
			mb->theta_map = Realloc(mb->theta_map, mb->ntheta + 1, map_func_tp *);
			mb->theta_map[mb->ntheta] = map_exp;
			mb->theta_map_arg = Realloc(mb->theta_map_arg, mb->ntheta + 1, void *);
			mb->theta_map_arg[mb->ntheta] = NULL;
			mb->ntheta++;
		}
	}
		break;

	case F_BESAG2:
	{
		tmp = iniparser_getdouble(ini, inla_string_join(secname, "INITIAL0"), G.log_prec_initial);
		if (!mb->f_fixed[mb->nf][0] && mb->mode_use_mode) {
			tmp = mb->theta_file[mb->theta_counter_file++];
			if (mb->mode_fixed) {
				mb->f_fixed[mb->nf][0] = 1;
			}
		}
		_SetInitial(0, tmp);
		HYPER_INIT(log_prec, tmp);
		if (mb->verbose) {
			printf("\t\tinitialise log_precision[%g]\n", tmp);
			printf("\t\tfixed=[%1d]\n", mb->f_fixed[mb->nf][0]);
		}
		mb->f_theta[mb->nf] = Calloc(2, double **);
		mb->f_theta[mb->nf][0] = log_prec;
		if (!mb->f_fixed[mb->nf][0]) {
			/*
			 * add this \theta 
			 */
			mb->theta = Realloc(mb->theta, mb->ntheta + 1, double **);
			mb->theta_hyperid = Realloc(mb->theta_hyperid, mb->ntheta + 1, char *);
			mb->theta_hyperid[mb->ntheta] = mb->f_prior[mb->nf][0].hyperid;
			mb->theta_tag = Realloc(mb->theta_tag, mb->ntheta + 1, char *);
			mb->theta_tag_userscale = Realloc(mb->theta_tag_userscale, mb->ntheta + 1, char *);
			mb->theta_dir = Realloc(mb->theta_dir, mb->ntheta + 1, char *);
			GMRFLib_sprintf(&msg, "Log precision for %s", (secname ? secname : mb->f_tag[mb->nf]));
			mb->theta_tag[mb->ntheta] = msg;
			GMRFLib_sprintf(&msg, "Precision for %s", (secname ? secname : mb->f_tag[mb->nf]));
			mb->theta_tag_userscale[mb->ntheta] = msg;
			GMRFLib_sprintf(&msg, "%s-parameter0", mb->f_dir[mb->nf]);
			mb->theta_dir[mb->ntheta] = msg;

			mb->theta_from = Realloc(mb->theta_from, mb->ntheta + 1, char *);
			mb->theta_to = Realloc(mb->theta_to, mb->ntheta + 1, char *);
			mb->theta_from[mb->ntheta] = Strdup(mb->f_prior[mb->nf][0].from_theta);
			mb->theta_to[mb->ntheta] = Strdup(mb->f_prior[mb->nf][0].to_theta);

			mb->theta[mb->ntheta] = log_prec;
			mb->theta_map = Realloc(mb->theta_map, mb->ntheta + 1, map_func_tp *);
			mb->theta_map[mb->ntheta] = map_precision;
			mb->theta_map_arg = Realloc(mb->theta_map_arg, mb->ntheta + 1, void *);
			mb->theta_map_arg[mb->ntheta] = NULL;
			mb->ntheta++;
		}

		tmp = iniparser_getdouble(ini, inla_string_join(secname, "INITIAL1"), 0.0);
		if (!mb->f_fixed[mb->nf][1] && mb->mode_use_mode) {
			tmp = mb->theta_file[mb->theta_counter_file++];
			if (mb->mode_fixed) {
				mb->f_fixed[mb->nf][1] = 1;
			}
		}
		_SetInitial(1, tmp);
		HYPER_INIT(a_intern, tmp);
		if (mb->verbose) {
			printf("\t\tinitialise a_intern[%g]\n", tmp);
			printf("\t\tfixed=[%1d]\n", mb->f_fixed[mb->nf][1]);
		}
		mb->f_theta[mb->nf][1] = a_intern;
		if (!mb->f_fixed[mb->nf][1]) {
			/*
			 * add this \theta 
			 */
			mb->theta = Realloc(mb->theta, mb->ntheta + 1, double **);
			mb->theta_hyperid = Realloc(mb->theta_hyperid, mb->ntheta + 1, char *);
			mb->theta_hyperid[mb->ntheta] = mb->f_prior[mb->nf][1].hyperid;
			mb->theta_tag = Realloc(mb->theta_tag, mb->ntheta + 1, char *);
			mb->theta_tag_userscale = Realloc(mb->theta_tag_userscale, mb->ntheta + 1, char *);
			mb->theta_dir = Realloc(mb->theta_dir, mb->ntheta + 1, char *);
			GMRFLib_sprintf(&msg, "Scale parameter a_intern for %s", (secname ? secname : mb->f_tag[mb->nf]));
			mb->theta_tag[mb->ntheta] = msg;
			GMRFLib_sprintf(&msg, "Scale paramter a for %s", (secname ? secname : mb->f_tag[mb->nf]));
			mb->theta_tag_userscale[mb->ntheta] = msg;
			GMRFLib_sprintf(&msg, "%s-parameter1", mb->f_dir[mb->nf]);
			mb->theta_dir[mb->ntheta] = msg;

			mb->theta_from = Realloc(mb->theta_from, mb->ntheta + 1, char *);
			mb->theta_to = Realloc(mb->theta_to, mb->ntheta + 1, char *);
			mb->theta_from[mb->ntheta] = Strdup(mb->f_prior[mb->nf][1].from_theta);
			mb->theta_to[mb->ntheta] = Strdup(mb->f_prior[mb->nf][1].to_theta);

			mb->theta[mb->ntheta] = a_intern;
			mb->theta_map = Realloc(mb->theta_map, mb->ntheta + 1, map_func_tp *);
			mb->theta_map[mb->ntheta] = map_exp;
			mb->theta_map_arg = Realloc(mb->theta_map_arg, mb->ntheta + 1, void *);
			mb->theta_map_arg[mb->ntheta] = NULL;
			mb->ntheta++;
		}
	}
		break;

	case F_BESAGPROPER:
	{
		tmp = iniparser_getdouble(ini, inla_string_join(secname, "INITIAL0"), G.log_prec_initial);
		if (!mb->f_fixed[mb->nf][0] && mb->mode_use_mode) {
			tmp = mb->theta_file[mb->theta_counter_file++];
			if (mb->mode_fixed) {
				mb->f_fixed[mb->nf][0] = 1;
			}
		}
		_SetInitial(0, tmp);
		HYPER_INIT(log_prec, tmp);
		if (mb->verbose) {
			printf("\t\tinitialise log_precision[%g]\n", tmp);
			printf("\t\tfixed=[%1d]\n", mb->f_fixed[mb->nf][0]);
		}
		mb->f_theta[mb->nf] = Calloc(2, double **);
		mb->f_theta[mb->nf][0] = log_prec;
		if (!mb->f_fixed[mb->nf][0]) {
			/*
			 * add this \theta 
			 */
			mb->theta = Realloc(mb->theta, mb->ntheta + 1, double **);
			mb->theta_hyperid = Realloc(mb->theta_hyperid, mb->ntheta + 1, char *);
			mb->theta_hyperid[mb->ntheta] = mb->f_prior[mb->nf][0].hyperid;
			mb->theta_tag = Realloc(mb->theta_tag, mb->ntheta + 1, char *);
			mb->theta_tag_userscale = Realloc(mb->theta_tag_userscale, mb->ntheta + 1, char *);
			mb->theta_dir = Realloc(mb->theta_dir, mb->ntheta + 1, char *);
			GMRFLib_sprintf(&msg, "Log precision for %s", (secname ? secname : mb->f_tag[mb->nf]));
			mb->theta_tag[mb->ntheta] = msg;
			GMRFLib_sprintf(&msg, "Precision for %s", (secname ? secname : mb->f_tag[mb->nf]));
			mb->theta_tag_userscale[mb->ntheta] = msg;
			GMRFLib_sprintf(&msg, "%s-parameter0", mb->f_dir[mb->nf]);
			mb->theta_dir[mb->ntheta] = msg;

			mb->theta_from = Realloc(mb->theta_from, mb->ntheta + 1, char *);
			mb->theta_to = Realloc(mb->theta_to, mb->ntheta + 1, char *);
			mb->theta_from[mb->ntheta] = Strdup(mb->f_prior[mb->nf][0].from_theta);
			mb->theta_to[mb->ntheta] = Strdup(mb->f_prior[mb->nf][0].to_theta);

			mb->theta[mb->ntheta] = log_prec;
			mb->theta_map = Realloc(mb->theta_map, mb->ntheta + 1, map_func_tp *);
			mb->theta_map[mb->ntheta] = map_precision;
			mb->theta_map_arg = Realloc(mb->theta_map_arg, mb->ntheta + 1, void *);
			mb->theta_map_arg[mb->ntheta] = NULL;
			mb->ntheta++;
		}

		tmp = iniparser_getdouble(ini, inla_string_join(secname, "INITIAL1"), 0.0);
		if (!mb->f_fixed[mb->nf][1] && mb->mode_use_mode) {
			tmp = mb->theta_file[mb->theta_counter_file++];
			if (mb->mode_fixed) {
				mb->f_fixed[mb->nf][1] = 1;
			}
		}
		_SetInitial(1, tmp);
		HYPER_INIT(log_diag, tmp);
		if (mb->verbose) {
			printf("\t\tinitialise log weight[%g]\n", tmp);
			printf("\t\tfixed=[%1d]\n", mb->f_fixed[mb->nf][1]);
		}
		mb->f_theta[mb->nf][1] = log_diag;
		if (!mb->f_fixed[mb->nf][1]) {
			/*
			 * add this \theta 
			 */
			mb->theta = Realloc(mb->theta, mb->ntheta + 1, double **);
			mb->theta_hyperid = Realloc(mb->theta_hyperid, mb->ntheta + 1, char *);
			mb->theta_hyperid[mb->ntheta] = mb->f_prior[mb->nf][1].hyperid;
			mb->theta_tag = Realloc(mb->theta_tag, mb->ntheta + 1, char *);
			mb->theta_tag_userscale = Realloc(mb->theta_tag_userscale, mb->ntheta + 1, char *);
			mb->theta_dir = Realloc(mb->theta_dir, mb->ntheta + 1, char *);
			GMRFLib_sprintf(&msg, "Log diagonal for %s", (secname ? secname : mb->f_tag[mb->nf]));
			mb->theta_tag[mb->ntheta] = msg;
			GMRFLib_sprintf(&msg, "Diagonal for %s", (secname ? secname : mb->f_tag[mb->nf]));
			mb->theta_tag_userscale[mb->ntheta] = msg;
			GMRFLib_sprintf(&msg, "%s-parameter1", mb->f_dir[mb->nf]);
			mb->theta_dir[mb->ntheta] = msg;

			mb->theta_from = Realloc(mb->theta_from, mb->ntheta + 1, char *);
			mb->theta_to = Realloc(mb->theta_to, mb->ntheta + 1, char *);
			mb->theta_from[mb->ntheta] = Strdup(mb->f_prior[mb->nf][1].from_theta);
			mb->theta_to[mb->ntheta] = Strdup(mb->f_prior[mb->nf][1].to_theta);

			mb->theta[mb->ntheta] = log_diag;
			mb->theta_map = Realloc(mb->theta_map, mb->ntheta + 1, map_func_tp *);
			mb->theta_map[mb->ntheta] = map_exp;
			mb->theta_map_arg = Realloc(mb->theta_map_arg, mb->ntheta + 1, void *);
			mb->theta_map_arg[mb->ntheta] = NULL;
			mb->ntheta++;
		}
	}
		break;

	case F_BESAGPROPER2:
	{
		tmp = iniparser_getdouble(ini, inla_string_join(secname, "INITIAL0"), G.log_prec_initial);
		if (!mb->f_fixed[mb->nf][0] && mb->mode_use_mode) {
			tmp = mb->theta_file[mb->theta_counter_file++];
			if (mb->mode_fixed) {
				mb->f_fixed[mb->nf][0] = 1;
			}
		}
		_SetInitial(0, tmp);
		HYPER_INIT(log_prec, tmp);
		if (mb->verbose) {
			printf("\t\tinitialise log_precision[%g]\n", tmp);
			printf("\t\tfixed=[%1d]\n", mb->f_fixed[mb->nf][0]);
		}
		mb->f_theta[mb->nf] = Calloc(2, double **);
		mb->f_theta[mb->nf][0] = log_prec;
		if (!mb->f_fixed[mb->nf][0]) {
			/*
			 * add this \theta 
			 */
			mb->theta = Realloc(mb->theta, mb->ntheta + 1, double **);
			mb->theta_hyperid = Realloc(mb->theta_hyperid, mb->ntheta + 1, char *);
			mb->theta_hyperid[mb->ntheta] = mb->f_prior[mb->nf][0].hyperid;
			mb->theta_tag = Realloc(mb->theta_tag, mb->ntheta + 1, char *);
			mb->theta_tag_userscale = Realloc(mb->theta_tag_userscale, mb->ntheta + 1, char *);
			mb->theta_dir = Realloc(mb->theta_dir, mb->ntheta + 1, char *);
			GMRFLib_sprintf(&msg, "Log precision for %s", (secname ? secname : mb->f_tag[mb->nf]));
			mb->theta_tag[mb->ntheta] = msg;
			GMRFLib_sprintf(&msg, "Precision for %s", (secname ? secname : mb->f_tag[mb->nf]));
			mb->theta_tag_userscale[mb->ntheta] = msg;
			GMRFLib_sprintf(&msg, "%s-parameter0", mb->f_dir[mb->nf]);
			mb->theta_dir[mb->ntheta] = msg;

			mb->theta_from = Realloc(mb->theta_from, mb->ntheta + 1, char *);
			mb->theta_to = Realloc(mb->theta_to, mb->ntheta + 1, char *);
			mb->theta_from[mb->ntheta] = Strdup(mb->f_prior[mb->nf][0].from_theta);
			mb->theta_to[mb->ntheta] = Strdup(mb->f_prior[mb->nf][0].to_theta);

			mb->theta[mb->ntheta] = log_prec;
			mb->theta_map = Realloc(mb->theta_map, mb->ntheta + 1, map_func_tp *);
			mb->theta_map[mb->ntheta] = map_precision;
			mb->theta_map_arg = Realloc(mb->theta_map_arg, mb->ntheta + 1, void *);
			mb->theta_map_arg[mb->ntheta] = NULL;
			mb->ntheta++;
		}

		tmp = iniparser_getdouble(ini, inla_string_join(secname, "INITIAL1"), 0.0);
		if (!mb->f_fixed[mb->nf][1] && mb->mode_use_mode) {
			tmp = mb->theta_file[mb->theta_counter_file++];
			if (mb->mode_fixed) {
				mb->f_fixed[mb->nf][1] = 1;
			}
		}
		_SetInitial(1, tmp);
		HYPER_INIT(phi_intern, tmp);
		if (mb->verbose) {
			printf("\t\tinitialise logit lambda[%g]\n", tmp);
			printf("\t\tfixed=[%1d]\n", mb->f_fixed[mb->nf][1]);
		}
		mb->f_theta[mb->nf][1] = phi_intern;
		if (!mb->f_fixed[mb->nf][1]) {
			/*
			 * add this \theta 
			 */
			mb->theta = Realloc(mb->theta, mb->ntheta + 1, double **);
			mb->theta_hyperid = Realloc(mb->theta_hyperid, mb->ntheta + 1, char *);
			mb->theta_hyperid[mb->ntheta] = mb->f_prior[mb->nf][1].hyperid;
			mb->theta_tag = Realloc(mb->theta_tag, mb->ntheta + 1, char *);
			mb->theta_tag_userscale = Realloc(mb->theta_tag_userscale, mb->ntheta + 1, char *);
			mb->theta_dir = Realloc(mb->theta_dir, mb->ntheta + 1, char *);
			GMRFLib_sprintf(&msg, "Logit lambda for %s", (secname ? secname : mb->f_tag[mb->nf]));
			mb->theta_tag[mb->ntheta] = msg;
			GMRFLib_sprintf(&msg, "Lambda for %s", (secname ? secname : mb->f_tag[mb->nf]));
			mb->theta_tag_userscale[mb->ntheta] = msg;
			GMRFLib_sprintf(&msg, "%s-parameter1", mb->f_dir[mb->nf]);
			mb->theta_dir[mb->ntheta] = msg;

			mb->theta_from = Realloc(mb->theta_from, mb->ntheta + 1, char *);
			mb->theta_to = Realloc(mb->theta_to, mb->ntheta + 1, char *);
			mb->theta_from[mb->ntheta] = Strdup(mb->f_prior[mb->nf][1].from_theta);
			mb->theta_to[mb->ntheta] = Strdup(mb->f_prior[mb->nf][1].to_theta);

			mb->theta[mb->ntheta] = phi_intern;
			mb->theta_map = Realloc(mb->theta_map, mb->ntheta + 1, map_func_tp *);
			mb->theta_map[mb->ntheta] = map_probability;
			mb->theta_map_arg = Realloc(mb->theta_map_arg, mb->ntheta + 1, void *);
			mb->theta_map_arg[mb->ntheta] = NULL;
			mb->ntheta++;
		}
	}
		break;

	case F_GENERIC1:
	{
		tmp = iniparser_getdouble(ini, inla_string_join(secname, "INITIAL0"), G.log_prec_initial);
		if (!mb->f_fixed[mb->nf][0] && mb->mode_use_mode) {
			tmp = mb->theta_file[mb->theta_counter_file++];
			if (mb->mode_fixed) {
				mb->f_fixed[mb->nf][0] = 1;
			}
		}
		_SetInitial(0, tmp);
		HYPER_INIT(log_prec, tmp);
		if (mb->verbose) {
			printf("\t\tinitialise log_precision[%g]\n", tmp);
			printf("\t\tfixed=[%1d]\n", mb->f_fixed[mb->nf][0]);
		}
		mb->f_theta[mb->nf] = Calloc(2, double **);
		mb->f_theta[mb->nf][0] = log_prec;
		if (!mb->f_fixed[mb->nf][0]) {
			/*
			 * add this \theta 
			 */
			mb->theta = Realloc(mb->theta, mb->ntheta + 1, double **);
			mb->theta_hyperid = Realloc(mb->theta_hyperid, mb->ntheta + 1, char *);
			mb->theta_hyperid[mb->ntheta] = mb->f_prior[mb->nf][0].hyperid;
			mb->theta_tag = Realloc(mb->theta_tag, mb->ntheta + 1, char *);
			mb->theta_tag_userscale = Realloc(mb->theta_tag_userscale, mb->ntheta + 1, char *);
			mb->theta_dir = Realloc(mb->theta_dir, mb->ntheta + 1, char *);
			GMRFLib_sprintf(&msg, "Log precision for %s", (secname ? secname : mb->f_tag[mb->nf]));
			mb->theta_tag[mb->ntheta] = msg;
			GMRFLib_sprintf(&msg, "Precision for %s", (secname ? secname : mb->f_tag[mb->nf]));
			mb->theta_tag_userscale[mb->ntheta] = msg;
			GMRFLib_sprintf(&msg, "%s-parameter0", mb->f_dir[mb->nf]);
			mb->theta_dir[mb->ntheta] = msg;

			mb->theta_from = Realloc(mb->theta_from, mb->ntheta + 1, char *);
			mb->theta_to = Realloc(mb->theta_to, mb->ntheta + 1, char *);
			mb->theta_from[mb->ntheta] = Strdup(mb->f_prior[mb->nf][0].from_theta);
			mb->theta_to[mb->ntheta] = Strdup(mb->f_prior[mb->nf][0].to_theta);

			mb->theta[mb->ntheta] = log_prec;
			mb->theta_map = Realloc(mb->theta_map, mb->ntheta + 1, map_func_tp *);
			mb->theta_map[mb->ntheta] = map_precision;
			mb->theta_map_arg = Realloc(mb->theta_map_arg, mb->ntheta + 1, void *);
			mb->theta_map_arg[mb->ntheta] = NULL;
			mb->ntheta++;
		}

		tmp = iniparser_getdouble(ini, inla_string_join(secname, "INITIAL1"), 0.0);
		if (!mb->f_fixed[mb->nf][1] && mb->mode_use_mode) {
			tmp = mb->theta_file[mb->theta_counter_file++];
			if (mb->mode_fixed) {
				mb->f_fixed[mb->nf][1] = 1;
			}
		}
		_SetInitial(1, tmp);
		HYPER_INIT(beta_intern, tmp);
		if (mb->verbose) {
			printf("\t\tinitialise beta_intern[%g]\n", tmp);
			printf("\t\tfixed=[%1d]\n", mb->f_fixed[mb->nf][1]);
		}
		mb->f_theta[mb->nf][1] = beta_intern;
		if (!mb->f_fixed[mb->nf][1]) {
			/*
			 * add this \theta 
			 */
			mb->theta = Realloc(mb->theta, mb->ntheta + 1, double **);
			mb->theta_hyperid = Realloc(mb->theta_hyperid, mb->ntheta + 1, char *);
			mb->theta_hyperid[mb->ntheta] = mb->f_prior[mb->nf][1].hyperid;
			mb->theta_tag = Realloc(mb->theta_tag, mb->ntheta + 1, char *);
			mb->theta_tag_userscale = Realloc(mb->theta_tag_userscale, mb->ntheta + 1, char *);
			mb->theta_dir = Realloc(mb->theta_dir, mb->ntheta + 1, char *);
			GMRFLib_sprintf(&msg, "Beta_intern for %s", (secname ? secname : mb->f_tag[mb->nf]));
			mb->theta_tag[mb->ntheta] = msg;
			GMRFLib_sprintf(&msg, "Beta for %s", (secname ? secname : mb->f_tag[mb->nf]));
			mb->theta_tag_userscale[mb->ntheta] = msg;
			GMRFLib_sprintf(&msg, "%s-parameter1", mb->f_dir[mb->nf]);
			mb->theta_dir[mb->ntheta] = msg;

			mb->theta_from = Realloc(mb->theta_from, mb->ntheta + 1, char *);
			mb->theta_to = Realloc(mb->theta_to, mb->ntheta + 1, char *);
			mb->theta_from[mb->ntheta] = Strdup(mb->f_prior[mb->nf][1].from_theta);
			mb->theta_to[mb->ntheta] = Strdup(mb->f_prior[mb->nf][1].to_theta);

			mb->theta[mb->ntheta] = beta_intern;
			mb->theta_map = Realloc(mb->theta_map, mb->ntheta + 1, map_func_tp *);
			mb->theta_map[mb->ntheta] = map_probability;
			mb->theta_map_arg = Realloc(mb->theta_map_arg, mb->ntheta + 1, void *);
			mb->theta_map_arg[mb->ntheta] = NULL;
			mb->ntheta++;
		}
	}
		break;

	case F_GENERIC2:
	{
		tmp = iniparser_getdouble(ini, inla_string_join(secname, "INITIAL0"), G.log_prec_initial);
		if (!mb->f_fixed[mb->nf][0] && mb->mode_use_mode) {
			tmp = mb->theta_file[mb->theta_counter_file++];
			if (mb->mode_fixed) {
				mb->f_fixed[mb->nf][0] = 1;
			}
		}
		_SetInitial(0, tmp);
		HYPER_INIT(log_prec, tmp);
		if (mb->verbose) {
			printf("\t\tinitialise log_precision[%g]\n", tmp);
			printf("\t\tfixed=[%1d]\n", mb->f_fixed[mb->nf][0]);
		}
		mb->f_theta[mb->nf] = Calloc(2, double **);
		mb->f_theta[mb->nf][0] = log_prec;
		if (!mb->f_fixed[mb->nf][0]) {
			/*
			 * add this \theta 
			 */
			mb->theta = Realloc(mb->theta, mb->ntheta + 1, double **);
			mb->theta_hyperid = Realloc(mb->theta_hyperid, mb->ntheta + 1, char *);
			mb->theta_hyperid[mb->ntheta] = mb->f_prior[mb->nf][0].hyperid;
			mb->theta_tag = Realloc(mb->theta_tag, mb->ntheta + 1, char *);
			mb->theta_tag_userscale = Realloc(mb->theta_tag_userscale, mb->ntheta + 1, char *);
			mb->theta_dir = Realloc(mb->theta_dir, mb->ntheta + 1, char *);
			GMRFLib_sprintf(&msg, "Log precision-cmatrix for %s", (secname ? secname : mb->f_tag[mb->nf]));
			mb->theta_tag[mb->ntheta] = msg;
			GMRFLib_sprintf(&msg, "Precision-cmatrix for %s", (secname ? secname : mb->f_tag[mb->nf]));
			mb->theta_tag_userscale[mb->ntheta] = msg;
			GMRFLib_sprintf(&msg, "%s-parameter0", mb->f_dir[mb->nf]);
			mb->theta_dir[mb->ntheta] = msg;

			mb->theta_from = Realloc(mb->theta_from, mb->ntheta + 1, char *);
			mb->theta_to = Realloc(mb->theta_to, mb->ntheta + 1, char *);
			mb->theta_from[mb->ntheta] = Strdup(mb->f_prior[mb->nf][0].from_theta);
			mb->theta_to[mb->ntheta] = Strdup(mb->f_prior[mb->nf][0].to_theta);

			mb->theta[mb->ntheta] = log_prec;
			mb->theta_map = Realloc(mb->theta_map, mb->ntheta + 1, map_func_tp *);
			mb->theta_map[mb->ntheta] = map_precision;
			mb->theta_map_arg = Realloc(mb->theta_map_arg, mb->ntheta + 1, void *);
			mb->theta_map_arg[mb->ntheta] = NULL;
			mb->ntheta++;
		}

		tmp = iniparser_getdouble(ini, inla_string_join(secname, "INITIAL1"), 0.0);
		if (!mb->f_fixed[mb->nf][1] && mb->mode_use_mode) {
			tmp = mb->theta_file[mb->theta_counter_file++];
			if (mb->mode_fixed) {
				mb->f_fixed[mb->nf][1] = 1;
			}
		}
		_SetInitial(1, tmp);
		HYPER_INIT(h2_intern, tmp);
		if (mb->verbose) {
			printf("\t\tinitialise h2-intern[%g]\n", tmp);
			printf("\t\tfixed=[%1d]\n", mb->f_fixed[mb->nf][1]);
		}
		mb->f_theta[mb->nf][1] = h2_intern;
		if (!mb->f_fixed[mb->nf][1]) {
			/*
			 * add this \theta 
			 */
			mb->theta = Realloc(mb->theta, mb->ntheta + 1, double **);
			mb->theta_hyperid = Realloc(mb->theta_hyperid, mb->ntheta + 1, char *);
			mb->theta_hyperid[mb->ntheta] = mb->f_prior[mb->nf][1].hyperid;
			mb->theta_tag = Realloc(mb->theta_tag, mb->ntheta + 1, char *);
			mb->theta_tag_userscale = Realloc(mb->theta_tag_userscale, mb->ntheta + 1, char *);
			mb->theta_dir = Realloc(mb->theta_dir, mb->ntheta + 1, char *);
			GMRFLib_sprintf(&msg, "h2-intern for %s", (secname ? secname : mb->f_tag[mb->nf]));
			mb->theta_tag[mb->ntheta] = msg;
			GMRFLib_sprintf(&msg, "h2 for %s", (secname ? secname : mb->f_tag[mb->nf]));
			mb->theta_tag_userscale[mb->ntheta] = msg;
			GMRFLib_sprintf(&msg, "%s-parameter1", mb->f_dir[mb->nf]);
			mb->theta_dir[mb->ntheta] = msg;

			mb->theta_from = Realloc(mb->theta_from, mb->ntheta + 1, char *);
			mb->theta_to = Realloc(mb->theta_to, mb->ntheta + 1, char *);
			mb->theta_from[mb->ntheta] = Strdup(mb->f_prior[mb->nf][1].from_theta);
			mb->theta_to[mb->ntheta] = Strdup(mb->f_prior[mb->nf][1].to_theta);

			mb->theta[mb->ntheta] = h2_intern;
			mb->theta_map = Realloc(mb->theta_map, mb->ntheta + 1, map_func_tp *);
			mb->theta_map[mb->ntheta] = map_probability;
			mb->theta_map_arg = Realloc(mb->theta_map_arg, mb->ntheta + 1, void *);
			mb->theta_map_arg[mb->ntheta] = NULL;
			mb->ntheta++;
		}
	}
		break;

	case F_GENERIC3:
	{
		inla_generic3_tp *a = (inla_generic3_tp *) mb->f_Qfunc_arg[mb->nf];
		a->log_prec = Calloc(GENERIC3_MAXTHETA, double **);
		mb->f_theta[mb->nf] = Calloc(GENERIC3_MAXTHETA, double **);

		for (k = a->m; k < GENERIC3_MAXTHETA - 1; k++) {	/* yes, do not include the common scaling */
			mb->f_fixed[mb->nf][k] = 1;	       /* those not used are set to fixed */
		}
		for (k = 0; k < GENERIC3_MAXTHETA; k++) {
			GMRFLib_sprintf(&ctmp, "INITIAL%1d", k);
			tmp = iniparser_getdouble(ini, inla_string_join(secname, ctmp), G.log_prec_initial);
			if (!mb->f_fixed[mb->nf][k] && mb->mode_use_mode) {
				tmp = mb->theta_file[mb->theta_counter_file++];
				if (mb->mode_fixed) {
					mb->f_fixed[mb->nf][k] = 1;
				}
			}
			_SetInitial(k, tmp);
			HYPER_NEW(a->log_prec[k], tmp);
			if (mb->verbose) {
				printf("\t\tinitialise log_precision[%1d] = %g\n", k, tmp);
				printf("\t\tfixed[%1d] = %1d\n", k, mb->f_fixed[mb->nf][k]);
			}
			mb->f_theta[mb->nf][k] = a->log_prec[k];

			if (!mb->f_fixed[mb->nf][k]) {
				/*
				 * add this \theta 
				 */
				mb->theta = Realloc(mb->theta, mb->ntheta + 1, double **);
				mb->theta_hyperid = Realloc(mb->theta_hyperid, mb->ntheta + 1, char *);
				mb->theta_hyperid[mb->ntheta] = mb->f_prior[mb->nf][k].hyperid;
				mb->theta_tag = Realloc(mb->theta_tag, mb->ntheta + 1, char *);
				mb->theta_tag_userscale = Realloc(mb->theta_tag_userscale, mb->ntheta + 1, char *);
				mb->theta_dir = Realloc(mb->theta_dir, mb->ntheta + 1, char *);
				if (k < GENERIC3_MAXTHETA - 1) {
					GMRFLib_sprintf(&msg, "Log precision for Cmatrix[[%1d]] for %s", k + 1,
							(secname ? secname : mb->f_tag[mb->nf]));
					mb->theta_tag[mb->ntheta] = msg;
					GMRFLib_sprintf(&msg, "Precision for Cmatrix[[%1d]] for %s", k + 1,
							(secname ? secname : mb->f_tag[mb->nf]));
					mb->theta_tag_userscale[mb->ntheta] = msg;
				} else {
					GMRFLib_sprintf(&msg, "Log common precision for %s", (secname ? secname : mb->f_tag[mb->nf]));
					mb->theta_tag[mb->ntheta] = msg;
					GMRFLib_sprintf(&msg, "Common precision for %s", (secname ? secname : mb->f_tag[mb->nf]));
					mb->theta_tag_userscale[mb->ntheta] = msg;
				}
				GMRFLib_sprintf(&msg, "%s-parameter%1d", mb->f_dir[mb->nf], k);
				mb->theta_dir[mb->ntheta] = msg;

				mb->theta_from = Realloc(mb->theta_from, mb->ntheta + 1, char *);
				mb->theta_to = Realloc(mb->theta_to, mb->ntheta + 1, char *);
				mb->theta_from[mb->ntheta] = Strdup(mb->f_prior[mb->nf][k].from_theta);
				mb->theta_to[mb->ntheta] = Strdup(mb->f_prior[mb->nf][k].to_theta);

				mb->theta[mb->ntheta] = a->log_prec[k];
				mb->theta_map = Realloc(mb->theta_map, mb->ntheta + 1, map_func_tp *);
				mb->theta_map[mb->ntheta] = map_precision;
				mb->theta_map_arg = Realloc(mb->theta_map_arg, mb->ntheta + 1, void *);
				mb->theta_map_arg[mb->ntheta] = NULL;
				mb->ntheta++;
			}
		}
	}
		break;

	case F_COPY:
	{
		mb->f_of[mb->nf] = iniparser_getstring(ini, inla_string_join(secname, "OF"), NULL);
		if (mb->verbose && mb->f_of[mb->nf]) {
			printf("\t\tof=[%s]\n", mb->f_of[mb->nf]);
		}

		/*
		 * same_as, says that the beta-parameters is the same as 'same_as', so this is to be determined later on. so we need to add space for it in
		 * f_theta, but not in mb->theta and mb->ntheta. The error-checking is done later. 
		 */
		mb->f_same_as[mb->nf] = iniparser_getstring(ini, inla_string_join(secname, "SAMEAS"), NULL);
		mb->f_same_as[mb->nf] = iniparser_getstring(ini, inla_string_join(secname, "SAME.AS"), mb->f_same_as[mb->nf]);
		if (mb->verbose) {
			printf("\t\tsame.as=[%s]\n", mb->f_same_as[mb->nf]);
		}

		mb->f_precision[mb->nf] = iniparser_getdouble(ini, inla_string_join(secname, "PRECISION"), mb->f_precision[mb->nf]);
		if (mb->verbose) {
			printf("\t\tprecision=[%f]\n", mb->f_precision[mb->nf]);
		}

		int fixed_default = 1;
		fixed_default = iniparser_getint(ini, inla_string_join(secname, "FIXED"), fixed_default);
		if (fixed_default == -1) {
			mb->f_fixed[mb->nf][0] = 1;
		}

		double *range = NULL;
		range = Calloc(2, double);		       /* need this as it will be stored in the map argument */
		range[0] = iniparser_getdouble(ini, inla_string_join(secname, "RANGE.LOW"), 0.0);	/* low = high ==> map = identity */
		range[1] = iniparser_getdouble(ini, inla_string_join(secname, "RANGE.HIGH"), 0.0);

		int aauto = 0;
		tmp = iniparser_getdouble(ini, inla_string_join(secname, "INITIAL"), 1.0);
		if (ISZERO(tmp)) {
			// initial=0.0 means auto-mode: initial=1 if FIXED and 0.1 if not
			tmp = (fixed_default ? 1.0 : 0.1);
			aauto = 1;
		}
		if (!mb->f_fixed[mb->nf][0] && mb->mode_use_mode) {
			tmp = mb->theta_file[mb->theta_counter_file++];
			if (mb->mode_fixed) {
				mb->f_fixed[mb->nf][0] = 1;
			}
		}
		_SetInitial(0, tmp);
		HYPER_INIT(beta, tmp);
		if (mb->verbose) {
			printf("\t\trange[%g, %g]\n", range[0], range[1]);
			if (aauto) {
				printf("\t\tauto-initialise beta[%g]\n", tmp);
			} else {
				printf("\t\tinitialise beta[%g]\n", tmp);
			}
			printf("\t\tfixed=[%1d]\n", mb->f_fixed[mb->nf][0]);
		}
		mb->f_theta[mb->nf] = Calloc(1, double **);
		mb->f_theta[mb->nf][0] = beta;

		mb->f_theta_map[mb->nf] = Calloc(1, map_func_tp *);
		mb->f_theta_map_arg[mb->nf] = Calloc(1, void *);
		mb->f_theta_map[mb->nf][0] = map_beta;	       /* need these */
		mb->f_theta_map_arg[mb->nf][0] = (void *) range;	/* and this one as well */

		if (mb->f_same_as[mb->nf] == NULL && !mb->f_fixed[mb->nf][0]) {
			/*
			 * add this \theta 
			 */
			mb->theta = Realloc(mb->theta, mb->ntheta + 1, double **);
			mb->theta_hyperid = Realloc(mb->theta_hyperid, mb->ntheta + 1, char *);
			mb->theta_hyperid[mb->ntheta] = mb->f_prior[mb->nf][0].hyperid;
			mb->theta_tag = Realloc(mb->theta_tag, mb->ntheta + 1, char *);
			mb->theta_tag_userscale = Realloc(mb->theta_tag_userscale, mb->ntheta + 1, char *);
			mb->theta_dir = Realloc(mb->theta_dir, mb->ntheta + 1, char *);
			GMRFLib_sprintf(&msg, "Beta_intern for %s", (secname ? secname : mb->f_tag[mb->nf]));
			mb->theta_tag[mb->ntheta] = msg;
			GMRFLib_sprintf(&msg, "Beta for %s", (secname ? secname : mb->f_tag[mb->nf]));
			mb->theta_tag_userscale[mb->ntheta] = msg;
			GMRFLib_sprintf(&msg, "%s-parameter", mb->f_dir[mb->nf]);
			mb->theta_dir[mb->ntheta] = msg;

			mb->theta_from = Realloc(mb->theta_from, mb->ntheta + 1, char *);
			mb->theta_to = Realloc(mb->theta_to, mb->ntheta + 1, char *);
			mb->theta_from[mb->ntheta] = Strdup(mb->f_prior[mb->nf][0].from_theta);
			mb->theta_to[mb->ntheta] = Strdup(mb->f_prior[mb->nf][0].to_theta);

			mb->theta[mb->ntheta] = beta;
			mb->theta_map = Realloc(mb->theta_map, mb->ntheta + 1, map_func_tp *);
			mb->theta_map[mb->ntheta] = map_beta;
			mb->theta_map_arg = Realloc(mb->theta_map_arg, mb->ntheta + 1, void *);
			mb->theta_map_arg[mb->ntheta] = (void *) range;
			mb->ntheta++;
		}
	}
		break;


	case F_SCOPY:
	{
		for (i = 0; i < SCOPY_MAXTHETA; i++) {
			GMRFLib_sprintf(&ctmp, "FIXED%1d", i);
			iniparser_getstring(ini, inla_string_join(secname, ctmp), NULL);
			GMRFLib_sprintf(&ctmp, "INITIAL%1d", i);
			iniparser_getstring(ini, inla_string_join(secname, ctmp), NULL);
			GMRFLib_sprintf(&ctmp, "PRIOR%1d", i);
			iniparser_getstring(ini, inla_string_join(secname, ctmp), NULL);
			GMRFLib_sprintf(&ctmp, "HYPERID%1d", i);
			iniparser_getstring(ini, inla_string_join(secname, ctmp), NULL);
			GMRFLib_sprintf(&ctmp, "PARAMETERS%1d", i);
			iniparser_getstring(ini, inla_string_join(secname, ctmp), NULL);
			GMRFLib_sprintf(&ctmp, "to.theta%1d", i);
			iniparser_getstring(ini, inla_string_join(secname, ctmp), NULL);
			GMRFLib_sprintf(&ctmp, "from.theta%1d", i);
			iniparser_getstring(ini, inla_string_join(secname, ctmp), NULL);
		}

		mb->f_of[mb->nf] = iniparser_getstring(ini, inla_string_join(secname, "OF"), NULL);
		if (mb->verbose && mb->f_of[mb->nf]) {
			printf("\t\tof=[%s]\n", mb->f_of[mb->nf]);
		}

		mb->f_precision[mb->nf] = iniparser_getdouble(ini, inla_string_join(secname, "PRECISION"), mb->f_precision[mb->nf]);
		if (mb->verbose) {
			printf("\t\tprecision=[%f]\n", mb->f_precision[mb->nf]);
		}

		nbeta = 0;
		nbeta = iniparser_getint(ini, inla_string_join(secname, "SCOPY.N"), nbeta);
		if (mb->verbose) {
			printf("\t\tnbeta=[%1d]\n", nbeta);
		}
		assert(nbeta <= SCOPY_MAXTHETA);

		char *filenameW = Strdup(iniparser_getstring(ini, inla_string_join(secname, "SCOPY.W"), NULL));
		if (!filenameW) {
			inla_error_missing_required_field(__GMRFLib_FuncName, secname, "W");
		}
		W = GMRFLib_read_fmesher_file(filenameW, (long int) 0, -1);

		filenamec = Strdup(iniparser_getstring(ini, inla_string_join(secname, "SCOPY.COVARIATE"), NULL));
		if (!filenamec) {
			inla_error_missing_required_field(__GMRFLib_FuncName, secname, "covariate");
		}
		GMRFLib_matrix_tp *cov_m = GMRFLib_read_fmesher_file(filenamec, (long int) 0, -1);
		cov = cov_m->A;
		ncov = cov_m->nrow;

		cov_m = NULL;				       /* that is ok */
		assert(cov);
		if (mb->verbose) {
			printf("\t\tread covariates from file=[%s]\n", filenamec);
			for (i = 0; i < IMIN(PREVIEW, ncov); i++) {
				printf("\t\t\tcovariate[%1d] = %g\n", i, cov[i]);
			}
		}
		assert(ncov >= 3);

		loc = Calloc(nbeta, double);
		double cov_min = GMRFLib_min_value(cov, ncov, NULL);
		double cov_max = GMRFLib_max_value(cov, ncov, NULL);
		double beta_step = (cov_max - cov_min) / (nbeta - 1.0);
		for (i = 0; i < nbeta; i++) {
			loc[i] = cov_min + i * beta_step;
		}
		loc[nbeta - 1] = cov_max;		       /* make it exact */
		if (mb->verbose) {
			printf("\t\tUse nbeta = %1d\n", nbeta);
			for (i = 0; i < nbeta; i++) {
				printf("\t\t\tlocation.beta[%1d] = %.3g\n", i, loc[i]);
			}
		}

		betas = Calloc(nbeta, double **);
		for (i = 0; i < nbeta; i++) {
			HYPER_NEW(betas[i], 1.0);
		}

		for (i = 0; i < SCOPY_MAXTHETA; i++) {
			GMRFLib_sprintf(&ctmp, "FIXED%1d", i);
			mb->f_fixed[mb->nf][i] = iniparser_getint(ini, inla_string_join(secname, ctmp), 0);

			GMRFLib_sprintf(&ctmp, "INITIAL%1d", i);
			double init = iniparser_getdouble(ini, inla_string_join(secname, ctmp), 1.0);

			if (i < nbeta) {
				if (!mb->f_fixed[mb->nf][i] && mb->mode_use_mode) {
					init = mb->theta_file[mb->theta_counter_file++];
					if (mb->mode_fixed) {
						mb->f_fixed[mb->nf][i] = 1;
					}
				}
				_SetInitial(i, init);
				HYPER_INIT(betas[i], init);
				if (mb->verbose) {
					printf("\t\tbeta[%1d] init  = %g\n", i, init);
					printf("\t\tbeta[%1d] fixed = %1d\n", i, mb->f_fixed[mb->nf][i]);
				}

				if (!mb->f_fixed[mb->nf][i]) {
					mb->theta = Realloc(mb->theta, mb->ntheta + 1, double **);
					mb->theta_hyperid = Realloc(mb->theta_hyperid, mb->ntheta + 1, char *);
					mb->theta_hyperid[mb->ntheta] = mb->f_prior[mb->nf][i].hyperid;
					mb->theta_tag = Realloc(mb->theta_tag, mb->ntheta + 1, char *);
					mb->theta_tag_userscale = Realloc(mb->theta_tag_userscale, mb->ntheta + 1, char *);
					mb->theta_dir = Realloc(mb->theta_dir, mb->ntheta + 1, char *);

					if (i == 0) {
						GMRFLib_sprintf(&msg, "Beta%1d for %s (scopy mean)", i, (secname ? secname : mb->f_tag[mb->nf]));
					} else if (i == 1) {
						GMRFLib_sprintf(&msg, "Beta%1d for %s (scopy slope)", i, (secname ? secname : mb->f_tag[mb->nf]));
					} else {
						GMRFLib_sprintf(&msg, "Beta%1d for %s (scopy theta)", i, (secname ? secname : mb->f_tag[mb->nf]));
					}
					mb->theta_tag[mb->ntheta] = msg;
					mb->theta_tag_userscale[mb->ntheta] = msg;
					GMRFLib_sprintf(&msg, "%s-parameter", mb->f_dir[mb->nf]);
					mb->theta_dir[mb->ntheta] = msg;

					mb->theta_from = Realloc(mb->theta_from, mb->ntheta + 1, char *);
					mb->theta_to = Realloc(mb->theta_to, mb->ntheta + 1, char *);
					mb->theta_from[mb->ntheta] = Strdup(mb->f_prior[mb->nf][0].from_theta);
					mb->theta_to[mb->ntheta] = Strdup(mb->f_prior[mb->nf][0].to_theta);

					mb->theta[mb->ntheta] = betas[i];
					mb->theta_map = Realloc(mb->theta_map, mb->ntheta + 1, map_func_tp *);
					mb->theta_map[mb->ntheta] = map_identity;
					mb->theta_map_arg = Realloc(mb->theta_map_arg, mb->ntheta + 1, void *);
					mb->ntheta++;
				}
			}
		}
	}
		break;

	case F_CLINEAR:
	{
		mb->f_precision[mb->nf] = iniparser_getdouble(ini, inla_string_join(secname, "PRECISION"), mb->f_precision[mb->nf]);
		if (mb->verbose) {
			printf("\t\tprecision=[%f]\n", mb->f_precision[mb->nf]);
		}

		int fixed_default = 0;
		fixed_default = iniparser_getint(ini, inla_string_join(secname, "FIXED"), fixed_default);
		if (fixed_default == -1) {
			mb->f_fixed[mb->nf][0] = 1;
		}
		if (mb->verbose && fixed_default == -1) {
			printf("\t\tfixed=[%d]\n", mb->f_fixed[mb->nf][0]);
		}

		double *range = NULL;
		range = Calloc(2, double);		       /* need this as it will be stored in the map argument */
		range[0] = iniparser_getdouble(ini, inla_string_join(secname, "RANGE.LOW"), 0.0);	/* low = high ==> map = identity */
		range[1] = iniparser_getdouble(ini, inla_string_join(secname, "RANGE.HIGH"), 0.0);

		tmp = iniparser_getdouble(ini, inla_string_join(secname, "INITIAL"), 0.0);
		if (!mb->f_fixed[mb->nf][0] && mb->mode_use_mode) {
			tmp = mb->theta_file[mb->theta_counter_file++];
			if (mb->mode_fixed) {
				mb->f_fixed[mb->nf][0] = 1;
			}
		}
		_SetInitial(0, tmp);
		HYPER_INIT(beta, tmp);
		if (mb->verbose) {
			printf("\t\trange[%g, %g]\n", range[0], range[1]);
			printf("\t\tinitialise beta[%g]\n", tmp);
			printf("\t\tfixed=[%1d]\n", mb->f_fixed[mb->nf][0]);
		}
		mb->f_theta[mb->nf] = Calloc(1, double **);
		mb->f_theta[mb->nf][0] = beta;

		mb->f_theta_map[mb->nf] = Calloc(1, map_func_tp *);
		mb->f_theta_map_arg[mb->nf] = Calloc(1, void *);
		mb->f_theta_map[mb->nf][0] = map_beta;	       /* need these */
		mb->f_theta_map_arg[mb->nf][0] = (void *) range;	/* and this one as well */

		if (!mb->f_fixed[mb->nf][0]) {
			/*
			 * add this \theta 
			 */
			mb->theta = Realloc(mb->theta, mb->ntheta + 1, double **);
			mb->theta_hyperid = Realloc(mb->theta_hyperid, mb->ntheta + 1, char *);
			mb->theta_hyperid[mb->ntheta] = mb->f_prior[mb->nf][0].hyperid;
			mb->theta_tag = Realloc(mb->theta_tag, mb->ntheta + 1, char *);
			mb->theta_tag_userscale = Realloc(mb->theta_tag_userscale, mb->ntheta + 1, char *);
			mb->theta_dir = Realloc(mb->theta_dir, mb->ntheta + 1, char *);
			GMRFLib_sprintf(&msg, "Beta_intern for %s", (secname ? secname : mb->f_tag[mb->nf]));
			mb->theta_tag[mb->ntheta] = msg;
			GMRFLib_sprintf(&msg, "Beta for %s", (secname ? secname : mb->f_tag[mb->nf]));
			mb->theta_tag_userscale[mb->ntheta] = msg;
			GMRFLib_sprintf(&msg, "%s-parameter", mb->f_dir[mb->nf]);
			mb->theta_dir[mb->ntheta] = msg;

			mb->theta_from = Realloc(mb->theta_from, mb->ntheta + 1, char *);
			mb->theta_to = Realloc(mb->theta_to, mb->ntheta + 1, char *);
			mb->theta_from[mb->ntheta] = Strdup(mb->f_prior[mb->nf][0].from_theta);
			mb->theta_to[mb->ntheta] = Strdup(mb->f_prior[mb->nf][0].to_theta);

			mb->theta[mb->ntheta] = beta;
			mb->theta_map = Realloc(mb->theta_map, mb->ntheta + 1, map_func_tp *);
			mb->theta_map[mb->ntheta] = map_beta;
			mb->theta_map_arg = Realloc(mb->theta_map_arg, mb->ntheta + 1, void *);
			mb->theta_map_arg[mb->ntheta] = (void *) range;
			mb->ntheta++;
		}
	}
		break;

	case F_SIGM:
	case F_REVSIGM:
	{
		char *local_name = (mb->f_id[mb->nf] == F_SIGM ? Strdup("SIGM") : Strdup("REVSIGM"));

		mb->f_precision[mb->nf] = iniparser_getdouble(ini, inla_string_join(secname, "PRECISION"), mb->f_precision[mb->nf]);
		tmp = iniparser_getdouble(ini, inla_string_join(secname, "INITIAL0"), 1.0);
		if (!mb->f_fixed[mb->nf][0] && mb->mode_use_mode) {
			tmp = mb->theta_file[mb->theta_counter_file++];
			if (mb->mode_fixed) {
				mb->f_fixed[mb->nf][0] = 1;
			}
		}
		_SetInitial(0, tmp);
		HYPER_INIT(beta, tmp);
		if (mb->verbose) {
			printf("\t\tinitialise beta[%g]\n", tmp);
			printf("\t\tfixed=[%1d]\n", mb->f_fixed[mb->nf][0]);
		}

		mb->f_theta[mb->nf][0] = beta;
		if (!mb->f_fixed[mb->nf][0]) {
			/*
			 * add this \theta 
			 */
			mb->theta = Realloc(mb->theta, mb->ntheta + 1, double **);
			mb->theta_hyperid = Realloc(mb->theta_hyperid, mb->ntheta + 1, char *);
			mb->theta_hyperid[mb->ntheta] = mb->f_prior[mb->nf][0].hyperid;
			mb->theta_tag = Realloc(mb->theta_tag, mb->ntheta + 1, char *);
			mb->theta_tag_userscale = Realloc(mb->theta_tag_userscale, mb->ntheta + 1, char *);
			mb->theta_dir = Realloc(mb->theta_dir, mb->ntheta + 1, char *);
			GMRFLib_sprintf(&msg, "%s beta for %s", local_name, (secname ? secname : mb->f_tag[mb->nf]));
			mb->theta_tag[mb->ntheta] = msg;
			GMRFLib_sprintf(&msg, "%s beta for %s", local_name, (secname ? secname : mb->f_tag[mb->nf]));
			mb->theta_tag_userscale[mb->ntheta] = msg;
			GMRFLib_sprintf(&msg, "%s-parameter0", mb->f_dir[mb->nf]);
			mb->theta_dir[mb->ntheta] = msg;

			mb->theta_from = Realloc(mb->theta_from, mb->ntheta + 1, char *);
			mb->theta_to = Realloc(mb->theta_to, mb->ntheta + 1, char *);
			mb->theta_from[mb->ntheta] = Strdup(mb->f_prior[mb->nf][0].from_theta);
			mb->theta_to[mb->ntheta] = Strdup(mb->f_prior[mb->nf][0].to_theta);

			mb->theta[mb->ntheta] = beta;
			mb->theta_map = Realloc(mb->theta_map, mb->ntheta + 1, map_func_tp *);
			mb->theta_map[mb->ntheta] = map_identity;
			mb->theta_map_arg = Realloc(mb->theta_map_arg, mb->ntheta + 1, void *);
			mb->theta_map_arg[mb->ntheta] = NULL;
			mb->ntheta++;
		}

		tmp = iniparser_getdouble(ini, inla_string_join(secname, "INITIAL1"), log(20.0));
		if (!mb->f_fixed[mb->nf][1] && mb->mode_use_mode) {
			tmp = mb->theta_file[mb->theta_counter_file++];
			if (mb->mode_fixed) {
				mb->f_fixed[mb->nf][1] = 1;
			}
		}
		_SetInitial(1, tmp);
		HYPER_INIT(log_halflife, tmp);
		if (mb->verbose) {
			printf("\t\tinitialise log_halflife[%g]\n", tmp);
			printf("\t\tfixed=[%1d]\n", mb->f_fixed[mb->nf][1]);
		}
		mb->f_theta[mb->nf][1] = log_halflife;
		if (!mb->f_fixed[mb->nf][1]) {
			/*
			 * add this \theta 
			 */
			mb->theta = Realloc(mb->theta, mb->ntheta + 1, double **);
			mb->theta_hyperid = Realloc(mb->theta_hyperid, mb->ntheta + 1, char *);
			mb->theta_hyperid[mb->ntheta] = mb->f_prior[mb->nf][1].hyperid;
			mb->theta_tag = Realloc(mb->theta_tag, mb->ntheta + 1, char *);
			mb->theta_tag_userscale = Realloc(mb->theta_tag_userscale, mb->ntheta + 1, char *);
			mb->theta_dir = Realloc(mb->theta_dir, mb->ntheta + 1, char *);

			GMRFLib_sprintf(&msg, "%s log_halflife for %s", local_name, (secname ? secname : mb->f_tag[mb->nf]));
			mb->theta_tag[mb->ntheta] = msg;
			GMRFLib_sprintf(&msg, "%s halflife for %s", local_name, (secname ? secname : mb->f_tag[mb->nf]));
			mb->theta_tag_userscale[mb->ntheta] = msg;
			GMRFLib_sprintf(&msg, "%s-parameter1", mb->f_dir[mb->nf]);
			mb->theta_dir[mb->ntheta] = msg;

			mb->theta_from = Realloc(mb->theta_from, mb->ntheta + 1, char *);
			mb->theta_to = Realloc(mb->theta_to, mb->ntheta + 1, char *);
			mb->theta_from[mb->ntheta] = Strdup(mb->f_prior[mb->nf][1].from_theta);
			mb->theta_to[mb->ntheta] = Strdup(mb->f_prior[mb->nf][1].to_theta);

			mb->theta[mb->ntheta] = log_halflife;
			mb->theta_map = Realloc(mb->theta_map, mb->ntheta + 1, map_func_tp *);
			mb->theta_map[mb->ntheta] = map_precision;
			mb->theta_map_arg = Realloc(mb->theta_map_arg, mb->ntheta + 1, void *);
			mb->theta_map_arg[mb->ntheta] = NULL;
			mb->ntheta++;
		}

		tmp = iniparser_getdouble(ini, inla_string_join(secname, "INITIAL2"), log(1.0));
		if (!mb->f_fixed[mb->nf][2] && mb->mode_use_mode) {
			tmp = mb->theta_file[mb->theta_counter_file++];
			if (mb->mode_fixed) {
				mb->f_fixed[mb->nf][2] = 1;
			}
		}
		_SetInitial(2, tmp);
		HYPER_INIT(log_shape, tmp);
		if (mb->verbose) {
			printf("\t\tinitialise log_shape[%g]\n", tmp);
			printf("\t\tfixed=[%1d]\n", mb->f_fixed[mb->nf][2]);
		}
		mb->f_theta[mb->nf][2] = log_shape;
		if (!mb->f_fixed[mb->nf][2]) {
			/*
			 * add this \theta 
			 */
			mb->theta = Realloc(mb->theta, mb->ntheta + 1, double **);
			mb->theta_hyperid = Realloc(mb->theta_hyperid, mb->ntheta + 1, char *);
			mb->theta_hyperid[mb->ntheta] = mb->f_prior[mb->nf][2].hyperid;
			mb->theta_tag = Realloc(mb->theta_tag, mb->ntheta + 1, char *);
			mb->theta_tag_userscale = Realloc(mb->theta_tag_userscale, mb->ntheta + 1, char *);
			mb->theta_dir = Realloc(mb->theta_dir, mb->ntheta + 1, char *);

			GMRFLib_sprintf(&msg, "%s log_shape for %s", local_name, (secname ? secname : mb->f_tag[mb->nf]));
			mb->theta_tag[mb->ntheta] = msg;
			GMRFLib_sprintf(&msg, "%s shape for %s", local_name, (secname ? secname : mb->f_tag[mb->nf]));
			mb->theta_tag_userscale[mb->ntheta] = msg;
			GMRFLib_sprintf(&msg, "%s-parameter1", mb->f_dir[mb->nf]);
			mb->theta_dir[mb->ntheta] = msg;

			mb->theta_from = Realloc(mb->theta_from, mb->ntheta + 1, char *);
			mb->theta_to = Realloc(mb->theta_to, mb->ntheta + 1, char *);
			mb->theta_from[mb->ntheta] = Strdup(mb->f_prior[mb->nf][2].from_theta);
			mb->theta_to[mb->ntheta] = Strdup(mb->f_prior[mb->nf][2].to_theta);

			mb->theta[mb->ntheta] = log_shape;
			mb->theta_map = Realloc(mb->theta_map, mb->ntheta + 1, map_func_tp *);
			mb->theta_map[mb->ntheta] = map_precision;
			mb->theta_map_arg = Realloc(mb->theta_map_arg, mb->ntheta + 1, void *);
			mb->theta_map_arg[mb->ntheta] = NULL;
			mb->ntheta++;
		}

	}
		break;

	case F_LOG1EXP:
	{
		char *local_name = Strdup("LOG1EXP");

		mb->f_precision[mb->nf] = iniparser_getdouble(ini, inla_string_join(secname, "PRECISION"), mb->f_precision[mb->nf]);
		tmp = iniparser_getdouble(ini, inla_string_join(secname, "INITIAL0"), 1.0);
		if (!mb->f_fixed[mb->nf][0] && mb->mode_use_mode) {
			tmp = mb->theta_file[mb->theta_counter_file++];
			if (mb->mode_fixed) {
				mb->f_fixed[mb->nf][0] = 1;
			}
		}
		_SetInitial(0, tmp);
		HYPER_INIT(beta, tmp);
		if (mb->verbose) {
			printf("\t\tinitialise beta[%g]\n", tmp);
			printf("\t\tfixed=[%1d]\n", mb->f_fixed[mb->nf][0]);
		}

		mb->f_theta[mb->nf][0] = beta;
		if (!mb->f_fixed[mb->nf][0]) {
			/*
			 * add this \theta 
			 */
			mb->theta = Realloc(mb->theta, mb->ntheta + 1, double **);
			mb->theta_hyperid = Realloc(mb->theta_hyperid, mb->ntheta + 1, char *);
			mb->theta_hyperid[mb->ntheta] = mb->f_prior[mb->nf][0].hyperid;
			mb->theta_tag = Realloc(mb->theta_tag, mb->ntheta + 1, char *);
			mb->theta_tag_userscale = Realloc(mb->theta_tag_userscale, mb->ntheta + 1, char *);
			mb->theta_dir = Realloc(mb->theta_dir, mb->ntheta + 1, char *);
			GMRFLib_sprintf(&msg, "%s beta for %s", local_name, (secname ? secname : mb->f_tag[mb->nf]));
			mb->theta_tag[mb->ntheta] = msg;
			GMRFLib_sprintf(&msg, "%s beta for %s", local_name, (secname ? secname : mb->f_tag[mb->nf]));
			mb->theta_tag_userscale[mb->ntheta] = msg;
			GMRFLib_sprintf(&msg, "%s-parameter0", mb->f_dir[mb->nf]);
			mb->theta_dir[mb->ntheta] = msg;

			mb->theta_from = Realloc(mb->theta_from, mb->ntheta + 1, char *);
			mb->theta_to = Realloc(mb->theta_to, mb->ntheta + 1, char *);
			mb->theta_from[mb->ntheta] = Strdup(mb->f_prior[mb->nf][0].from_theta);
			mb->theta_to[mb->ntheta] = Strdup(mb->f_prior[mb->nf][0].to_theta);

			mb->theta[mb->ntheta] = beta;
			mb->theta_map = Realloc(mb->theta_map, mb->ntheta + 1, map_func_tp *);
			mb->theta_map[mb->ntheta] = map_identity;
			mb->theta_map_arg = Realloc(mb->theta_map_arg, mb->ntheta + 1, void *);
			mb->theta_map_arg[mb->ntheta] = NULL;
			mb->ntheta++;
		}

		tmp = iniparser_getdouble(ini, inla_string_join(secname, "INITIAL1"), 0.0);
		if (!mb->f_fixed[mb->nf][1] && mb->mode_use_mode) {
			tmp = mb->theta_file[mb->theta_counter_file++];
			if (mb->mode_fixed) {
				mb->f_fixed[mb->nf][1] = 1;
			}
		}
		_SetInitial(1, tmp);
		HYPER_INIT(alpha, tmp);
		if (mb->verbose) {
			printf("\t\tinitialise alpha[%g]\n", tmp);
			printf("\t\tfixed=[%1d]\n", mb->f_fixed[mb->nf][1]);
		}
		mb->f_theta[mb->nf][1] = alpha;
		if (!mb->f_fixed[mb->nf][1]) {
			/*
			 * add this \theta 
			 */
			mb->theta = Realloc(mb->theta, mb->ntheta + 1, double **);
			mb->theta_hyperid = Realloc(mb->theta_hyperid, mb->ntheta + 1, char *);
			mb->theta_hyperid[mb->ntheta] = mb->f_prior[mb->nf][1].hyperid;
			mb->theta_tag = Realloc(mb->theta_tag, mb->ntheta + 1, char *);
			mb->theta_tag_userscale = Realloc(mb->theta_tag_userscale, mb->ntheta + 1, char *);
			mb->theta_dir = Realloc(mb->theta_dir, mb->ntheta + 1, char *);

			GMRFLib_sprintf(&msg, "%s alpha for %s", local_name, (secname ? secname : mb->f_tag[mb->nf]));
			mb->theta_tag[mb->ntheta] = msg;
			GMRFLib_sprintf(&msg, "%s alpha for %s", local_name, (secname ? secname : mb->f_tag[mb->nf]));
			mb->theta_tag_userscale[mb->ntheta] = msg;
			GMRFLib_sprintf(&msg, "%s-parameter1", mb->f_dir[mb->nf]);
			mb->theta_dir[mb->ntheta] = msg;

			mb->theta_from = Realloc(mb->theta_from, mb->ntheta + 1, char *);
			mb->theta_to = Realloc(mb->theta_to, mb->ntheta + 1, char *);
			mb->theta_from[mb->ntheta] = Strdup(mb->f_prior[mb->nf][1].from_theta);
			mb->theta_to[mb->ntheta] = Strdup(mb->f_prior[mb->nf][1].to_theta);

			mb->theta[mb->ntheta] = alpha;
			mb->theta_map = Realloc(mb->theta_map, mb->ntheta + 1, map_func_tp *);
			mb->theta_map[mb->ntheta] = map_identity;
			mb->theta_map_arg = Realloc(mb->theta_map_arg, mb->ntheta + 1, void *);
			mb->theta_map_arg[mb->ntheta] = NULL;
			mb->ntheta++;
		}

		tmp = iniparser_getdouble(ini, inla_string_join(secname, "INITIAL2"), 0.0);
		if (!mb->f_fixed[mb->nf][2] && mb->mode_use_mode) {
			tmp = mb->theta_file[mb->theta_counter_file++];
			if (mb->mode_fixed) {
				mb->f_fixed[mb->nf][2] = 1;
			}
		}
		_SetInitial(2, tmp);
		HYPER_INIT(gama, tmp);
		if (mb->verbose) {
			printf("\t\tinitialise gamma[%g]\n", tmp);
			printf("\t\tfixed=[%1d]\n", mb->f_fixed[mb->nf][2]);
		}
		mb->f_theta[mb->nf][2] = gama;
		if (!mb->f_fixed[mb->nf][2]) {
			/*
			 * add this \theta 
			 */
			mb->theta = Realloc(mb->theta, mb->ntheta + 1, double **);
			mb->theta_hyperid = Realloc(mb->theta_hyperid, mb->ntheta + 1, char *);
			mb->theta_hyperid[mb->ntheta] = mb->f_prior[mb->nf][2].hyperid;
			mb->theta_tag = Realloc(mb->theta_tag, mb->ntheta + 1, char *);
			mb->theta_tag_userscale = Realloc(mb->theta_tag_userscale, mb->ntheta + 1, char *);
			mb->theta_dir = Realloc(mb->theta_dir, mb->ntheta + 1, char *);

			GMRFLib_sprintf(&msg, "%s gamma for %s", local_name, (secname ? secname : mb->f_tag[mb->nf]));
			mb->theta_tag[mb->ntheta] = msg;
			GMRFLib_sprintf(&msg, "%s gamma for %s", local_name, (secname ? secname : mb->f_tag[mb->nf]));
			mb->theta_tag_userscale[mb->ntheta] = msg;
			GMRFLib_sprintf(&msg, "%s-parameter1", mb->f_dir[mb->nf]);
			mb->theta_dir[mb->ntheta] = msg;

			mb->theta_from = Realloc(mb->theta_from, mb->ntheta + 1, char *);
			mb->theta_to = Realloc(mb->theta_to, mb->ntheta + 1, char *);
			mb->theta_from[mb->ntheta] = Strdup(mb->f_prior[mb->nf][2].from_theta);
			mb->theta_to[mb->ntheta] = Strdup(mb->f_prior[mb->nf][2].to_theta);

			mb->theta[mb->ntheta] = gama;
			mb->theta_map = Realloc(mb->theta_map, mb->ntheta + 1, map_func_tp *);
			mb->theta_map[mb->ntheta] = map_identity;
			mb->theta_map_arg = Realloc(mb->theta_map_arg, mb->ntheta + 1, void *);
			mb->theta_map_arg[mb->ntheta] = NULL;
			mb->ntheta++;
		}

	}
		break;

	case F_LOGDIST:
	{
		char *local_name = Strdup("LOGDIST");

		mb->f_precision[mb->nf] = iniparser_getdouble(ini, inla_string_join(secname, "PRECISION"), mb->f_precision[mb->nf]);
		tmp = iniparser_getdouble(ini, inla_string_join(secname, "INITIAL0"), 1.0);
		if (!mb->f_fixed[mb->nf][0] && mb->mode_use_mode) {
			tmp = mb->theta_file[mb->theta_counter_file++];
			if (mb->mode_fixed) {
				mb->f_fixed[mb->nf][0] = 1;
			}
		}
		_SetInitial(0, tmp);
		HYPER_INIT(beta, tmp);
		if (mb->verbose) {
			printf("\t\tinitialise beta[%g]\n", tmp);
			printf("\t\tfixed=[%1d]\n", mb->f_fixed[mb->nf][0]);
		}

		mb->f_theta[mb->nf][0] = beta;
		if (!mb->f_fixed[mb->nf][0]) {
			/*
			 * add this \theta 
			 */
			mb->theta = Realloc(mb->theta, mb->ntheta + 1, double **);
			mb->theta_hyperid = Realloc(mb->theta_hyperid, mb->ntheta + 1, char *);
			mb->theta_hyperid[mb->ntheta] = mb->f_prior[mb->nf][0].hyperid;
			mb->theta_tag = Realloc(mb->theta_tag, mb->ntheta + 1, char *);
			mb->theta_tag_userscale = Realloc(mb->theta_tag_userscale, mb->ntheta + 1, char *);
			mb->theta_dir = Realloc(mb->theta_dir, mb->ntheta + 1, char *);
			GMRFLib_sprintf(&msg, "%s beta for %s", local_name, (secname ? secname : mb->f_tag[mb->nf]));
			mb->theta_tag[mb->ntheta] = msg;
			GMRFLib_sprintf(&msg, "%s beta for %s", local_name, (secname ? secname : mb->f_tag[mb->nf]));
			mb->theta_tag_userscale[mb->ntheta] = msg;
			GMRFLib_sprintf(&msg, "%s-parameter0", mb->f_dir[mb->nf]);
			mb->theta_dir[mb->ntheta] = msg;

			mb->theta_from = Realloc(mb->theta_from, mb->ntheta + 1, char *);
			mb->theta_to = Realloc(mb->theta_to, mb->ntheta + 1, char *);
			mb->theta_from[mb->ntheta] = Strdup(mb->f_prior[mb->nf][0].from_theta);
			mb->theta_to[mb->ntheta] = Strdup(mb->f_prior[mb->nf][0].to_theta);

			mb->theta[mb->ntheta] = beta;
			mb->theta_map = Realloc(mb->theta_map, mb->ntheta + 1, map_func_tp *);
			mb->theta_map[mb->ntheta] = map_identity;
			mb->theta_map_arg = Realloc(mb->theta_map_arg, mb->ntheta + 1, void *);
			mb->theta_map_arg[mb->ntheta] = NULL;
			mb->ntheta++;
		}

		tmp = iniparser_getdouble(ini, inla_string_join(secname, "INITIAL1"), 0.0);
		if (!mb->f_fixed[mb->nf][1] && mb->mode_use_mode) {
			tmp = mb->theta_file[mb->theta_counter_file++];
			if (mb->mode_fixed) {
				mb->f_fixed[mb->nf][1] = 1;
			}
		}
		_SetInitial(1, tmp);
		HYPER_INIT(alpha1, tmp);
		if (mb->verbose) {
			printf("\t\tinitialise alpha1[%g]\n", tmp);
			printf("\t\tfixed=[%1d]\n", mb->f_fixed[mb->nf][1]);
		}
		mb->f_theta[mb->nf][1] = alpha1;
		if (!mb->f_fixed[mb->nf][1]) {
			/*
			 * add this \theta 
			 */
			mb->theta = Realloc(mb->theta, mb->ntheta + 1, double **);
			mb->theta_hyperid = Realloc(mb->theta_hyperid, mb->ntheta + 1, char *);
			mb->theta_hyperid[mb->ntheta] = mb->f_prior[mb->nf][1].hyperid;
			mb->theta_tag = Realloc(mb->theta_tag, mb->ntheta + 1, char *);
			mb->theta_tag_userscale = Realloc(mb->theta_tag_userscale, mb->ntheta + 1, char *);
			mb->theta_dir = Realloc(mb->theta_dir, mb->ntheta + 1, char *);

			GMRFLib_sprintf(&msg, "%s alpha1.intern for %s", local_name, (secname ? secname : mb->f_tag[mb->nf]));
			mb->theta_tag[mb->ntheta] = msg;
			GMRFLib_sprintf(&msg, "%s alpha1 for %s", local_name, (secname ? secname : mb->f_tag[mb->nf]));
			mb->theta_tag_userscale[mb->ntheta] = msg;
			GMRFLib_sprintf(&msg, "%s-parameter1", mb->f_dir[mb->nf]);
			mb->theta_dir[mb->ntheta] = msg;

			mb->theta_from = Realloc(mb->theta_from, mb->ntheta + 1, char *);
			mb->theta_to = Realloc(mb->theta_to, mb->ntheta + 1, char *);
			mb->theta_from[mb->ntheta] = Strdup(mb->f_prior[mb->nf][1].from_theta);
			mb->theta_to[mb->ntheta] = Strdup(mb->f_prior[mb->nf][1].to_theta);

			mb->theta[mb->ntheta] = alpha1;
			mb->theta_map = Realloc(mb->theta_map, mb->ntheta + 1, map_func_tp *);
			mb->theta_map[mb->ntheta] = map_exp;
			mb->theta_map_arg = Realloc(mb->theta_map_arg, mb->ntheta + 1, void *);
			mb->theta_map_arg[mb->ntheta] = NULL;
			mb->ntheta++;
		}

		tmp = iniparser_getdouble(ini, inla_string_join(secname, "INITIAL2"), 0.0);
		if (!mb->f_fixed[mb->nf][2] && mb->mode_use_mode) {
			tmp = mb->theta_file[mb->theta_counter_file++];
			if (mb->mode_fixed) {
				mb->f_fixed[mb->nf][2] = 1;
			}
		}
		_SetInitial(2, tmp);
		HYPER_INIT(alpha2, tmp);
		if (mb->verbose) {
			printf("\t\tinitialise alpha2[%g]\n", tmp);
			printf("\t\tfixed=[%1d]\n", mb->f_fixed[mb->nf][2]);
		}
		mb->f_theta[mb->nf][2] = alpha2;
		if (!mb->f_fixed[mb->nf][2]) {
			/*
			 * add this \theta 
			 */
			mb->theta = Realloc(mb->theta, mb->ntheta + 1, double **);
			mb->theta_hyperid = Realloc(mb->theta_hyperid, mb->ntheta + 1, char *);
			mb->theta_hyperid[mb->ntheta] = mb->f_prior[mb->nf][2].hyperid;
			mb->theta_tag = Realloc(mb->theta_tag, mb->ntheta + 1, char *);
			mb->theta_tag_userscale = Realloc(mb->theta_tag_userscale, mb->ntheta + 1, char *);
			mb->theta_dir = Realloc(mb->theta_dir, mb->ntheta + 1, char *);

			GMRFLib_sprintf(&msg, "%s alpha2.intern for %s", local_name, (secname ? secname : mb->f_tag[mb->nf]));
			mb->theta_tag[mb->ntheta] = msg;
			GMRFLib_sprintf(&msg, "%s alpha2 for %s", local_name, (secname ? secname : mb->f_tag[mb->nf]));
			mb->theta_tag_userscale[mb->ntheta] = msg;
			GMRFLib_sprintf(&msg, "%s-parameter1", mb->f_dir[mb->nf]);
			mb->theta_dir[mb->ntheta] = msg;

			mb->theta_from = Realloc(mb->theta_from, mb->ntheta + 1, char *);
			mb->theta_to = Realloc(mb->theta_to, mb->ntheta + 1, char *);
			mb->theta_from[mb->ntheta] = Strdup(mb->f_prior[mb->nf][2].from_theta);
			mb->theta_to[mb->ntheta] = Strdup(mb->f_prior[mb->nf][2].to_theta);

			mb->theta[mb->ntheta] = alpha2;
			mb->theta_map = Realloc(mb->theta_map, mb->ntheta + 1, map_func_tp *);
			mb->theta_map[mb->ntheta] = map_exp;
			mb->theta_map_arg = Realloc(mb->theta_map_arg, mb->ntheta + 1, void *);
			mb->theta_map_arg[mb->ntheta] = NULL;
			mb->ntheta++;
		}

	}
		break;

	case F_BYM:
	{
		tmp = iniparser_getdouble(ini, inla_string_join(secname, "INITIAL0"), G.log_prec_initial);
		if (!mb->f_fixed[mb->nf][0] && mb->mode_use_mode) {
			tmp = mb->theta_file[mb->theta_counter_file++];
			if (mb->mode_fixed) {
				mb->f_fixed[mb->nf][0] = 1;
			}
		}
		_SetInitial(0, tmp);
		HYPER_INIT(log_prec0, tmp);
		if (mb->verbose) {
			printf("\t\tinitialise log_precision (iid component)[%g]\n", tmp);
			printf("\t\tfixed=[%1d]\n", mb->f_fixed[mb->nf][0]);
		}
		tmp = iniparser_getdouble(ini, inla_string_join(secname, "INITIAL1"), G.log_prec_initial);
		if (!mb->f_fixed[mb->nf][1] && mb->mode_use_mode) {
			tmp = mb->theta_file[mb->theta_counter_file++];
			if (mb->mode_fixed) {
				mb->f_fixed[mb->nf][1] = 1;
			}
		}
		_SetInitial(1, tmp);
		HYPER_INIT(log_prec1, tmp);
		if (mb->verbose) {
			printf("\t\tinitialise log_precision (spatial component)[%g]\n", tmp);
			printf("\t\tfixed=[%1d]\n", mb->f_fixed[mb->nf][1]);
		}
		mb->f_theta[mb->nf] = Calloc(2, double **);
		mb->f_theta[mb->nf][0] = log_prec0;
		if (!mb->f_fixed[mb->nf][0]) {
			mb->theta = Realloc(mb->theta, mb->ntheta + 1, double **);
			mb->theta_hyperid = Realloc(mb->theta_hyperid, mb->ntheta + 1, char *);
			mb->theta_hyperid[mb->ntheta] = mb->f_prior[mb->nf][0].hyperid;
			mb->theta_tag = Realloc(mb->theta_tag, mb->ntheta + 1, char *);
			mb->theta_tag_userscale = Realloc(mb->theta_tag_userscale, mb->ntheta + 1, char *);
			mb->theta_dir = Realloc(mb->theta_dir, mb->ntheta + 1, char *);
			GMRFLib_sprintf(&msg, "Log precision for %s (idd component)", mb->f_tag[mb->nf]);
			mb->theta_tag[mb->ntheta] = msg;
			GMRFLib_sprintf(&msg, "Precision for %s (iid component)", mb->f_tag[mb->nf]);
			mb->theta_tag_userscale[mb->ntheta] = msg;
			GMRFLib_sprintf(&msg, "%s-parameter0", mb->f_dir[mb->nf]);
			mb->theta_dir[mb->ntheta] = msg;

			mb->theta_from = Realloc(mb->theta_from, mb->ntheta + 1, char *);
			mb->theta_to = Realloc(mb->theta_to, mb->ntheta + 1, char *);
			mb->theta_from[mb->ntheta] = Strdup(mb->f_prior[mb->nf][0].from_theta);
			mb->theta_to[mb->ntheta] = Strdup(mb->f_prior[mb->nf][0].to_theta);

			mb->theta[mb->ntheta] = log_prec0;
			mb->theta_map = Realloc(mb->theta_map, mb->ntheta + 1, map_func_tp *);
			mb->theta_map[mb->ntheta] = map_precision;
			mb->theta_map_arg = Realloc(mb->theta_map_arg, mb->ntheta + 1, void *);
			mb->theta_map_arg[mb->ntheta] = NULL;
			mb->ntheta++;
		}
		mb->f_theta[mb->nf][1] = log_prec1;
		if (!mb->f_fixed[mb->nf][1]) {
			mb->theta = Realloc(mb->theta, mb->ntheta + 1, double **);
			mb->theta_hyperid = Realloc(mb->theta_hyperid, mb->ntheta + 1, char *);
			mb->theta_hyperid[mb->ntheta] = mb->f_prior[mb->nf][1].hyperid;
			mb->theta_tag = Realloc(mb->theta_tag, mb->ntheta + 1, char *);
			mb->theta_tag_userscale = Realloc(mb->theta_tag_userscale, mb->ntheta + 1, char *);
			mb->theta_dir = Realloc(mb->theta_dir, mb->ntheta + 1, char *);
			GMRFLib_sprintf(&msg, "Log precision for %s (spatial component)", mb->f_tag[mb->nf]);
			mb->theta_tag[mb->ntheta] = msg;
			GMRFLib_sprintf(&msg, "Precision for %s (spatial component)", mb->f_tag[mb->nf]);
			mb->theta_tag_userscale[mb->ntheta] = msg;
			GMRFLib_sprintf(&msg, "%s-parameter1", mb->f_dir[mb->nf]);
			mb->theta_dir[mb->ntheta] = msg;

			mb->theta_from = Realloc(mb->theta_from, mb->ntheta + 1, char *);
			mb->theta_to = Realloc(mb->theta_to, mb->ntheta + 1, char *);
			mb->theta_from[mb->ntheta] = Strdup(mb->f_prior[mb->nf][1].from_theta);
			mb->theta_to[mb->ntheta] = Strdup(mb->f_prior[mb->nf][1].to_theta);

			mb->theta[mb->ntheta] = log_prec1;
			mb->theta_map = Realloc(mb->theta_map, mb->ntheta + 1, map_func_tp *);
			mb->theta_map[mb->ntheta] = map_precision;
			mb->theta_map_arg = Realloc(mb->theta_map_arg, mb->ntheta + 1, void *);
			mb->theta_map_arg[mb->ntheta] = NULL;
			mb->ntheta++;
		}
	}
		break;

	case F_BYM2:
	{
		tmp = iniparser_getdouble(ini, inla_string_join(secname, "INITIAL0"), G.log_prec_initial);
		if (!mb->f_fixed[mb->nf][0] && mb->mode_use_mode) {
			tmp = mb->theta_file[mb->theta_counter_file++];
			if (mb->mode_fixed) {
				mb->f_fixed[mb->nf][0] = 1;
			}
		}
		_SetInitial(0, tmp);
		HYPER_INIT(log_prec0, tmp);
		if (mb->verbose) {
			printf("\t\tinitialise log_precision [%g]\n", tmp);
			printf("\t\tfixed=[%1d]\n", mb->f_fixed[mb->nf][0]);
		}

		tmp = iniparser_getdouble(ini, inla_string_join(secname, "INITIAL1"), 0.0);
		if (!mb->f_fixed[mb->nf][1] && mb->mode_use_mode) {
			tmp = mb->theta_file[mb->theta_counter_file++];
			if (mb->mode_fixed) {
				mb->f_fixed[mb->nf][1] = 1;
			}
		}
		_SetInitial(1, tmp);
		HYPER_INIT(phi_intern, tmp);
		if (mb->verbose) {
			printf("\t\tinitialise phi_intern [%g]\n", tmp);
			printf("\t\tfixed=[%1d]\n", mb->f_fixed[mb->nf][1]);
		}
		mb->f_theta[mb->nf] = Calloc(2, double **);
		mb->f_theta[mb->nf][0] = log_prec0;
		mb->f_theta[mb->nf][1] = phi_intern;

		if (!mb->f_fixed[mb->nf][0]) {
			mb->theta = Realloc(mb->theta, mb->ntheta + 1, double **);
			mb->theta_hyperid = Realloc(mb->theta_hyperid, mb->ntheta + 1, char *);
			mb->theta_hyperid[mb->ntheta] = mb->f_prior[mb->nf][0].hyperid;
			mb->theta_tag = Realloc(mb->theta_tag, mb->ntheta + 1, char *);
			mb->theta_tag_userscale = Realloc(mb->theta_tag_userscale, mb->ntheta + 1, char *);
			mb->theta_dir = Realloc(mb->theta_dir, mb->ntheta + 1, char *);
			GMRFLib_sprintf(&msg, "Log precision for %s", mb->f_tag[mb->nf]);
			mb->theta_tag[mb->ntheta] = msg;
			GMRFLib_sprintf(&msg, "Precision for %s", mb->f_tag[mb->nf]);
			mb->theta_tag_userscale[mb->ntheta] = msg;
			GMRFLib_sprintf(&msg, "%s-parameter0", mb->f_dir[mb->nf]);
			mb->theta_dir[mb->ntheta] = msg;

			mb->theta_from = Realloc(mb->theta_from, mb->ntheta + 1, char *);
			mb->theta_to = Realloc(mb->theta_to, mb->ntheta + 1, char *);
			mb->theta_from[mb->ntheta] = Strdup(mb->f_prior[mb->nf][0].from_theta);
			mb->theta_to[mb->ntheta] = Strdup(mb->f_prior[mb->nf][0].to_theta);

			mb->theta[mb->ntheta] = log_prec0;
			mb->theta_map = Realloc(mb->theta_map, mb->ntheta + 1, map_func_tp *);
			mb->theta_map[mb->ntheta] = map_precision;
			mb->theta_map_arg = Realloc(mb->theta_map_arg, mb->ntheta + 1, void *);
			mb->theta_map_arg[mb->ntheta] = NULL;
			mb->ntheta++;
		}

		if (!mb->f_fixed[mb->nf][1]) {
			mb->theta = Realloc(mb->theta, mb->ntheta + 1, double **);
			mb->theta_hyperid = Realloc(mb->theta_hyperid, mb->ntheta + 1, char *);
			mb->theta_hyperid[mb->ntheta] = mb->f_prior[mb->nf][1].hyperid;
			mb->theta_tag = Realloc(mb->theta_tag, mb->ntheta + 1, char *);
			mb->theta_tag_userscale = Realloc(mb->theta_tag_userscale, mb->ntheta + 1, char *);
			mb->theta_dir = Realloc(mb->theta_dir, mb->ntheta + 1, char *);
			GMRFLib_sprintf(&msg, "Logit phi for %s", mb->f_tag[mb->nf]);
			mb->theta_tag[mb->ntheta] = msg;
			GMRFLib_sprintf(&msg, "Phi for %s", mb->f_tag[mb->nf]);
			mb->theta_tag_userscale[mb->ntheta] = msg;
			GMRFLib_sprintf(&msg, "%s-parameter1", mb->f_dir[mb->nf]);
			mb->theta_dir[mb->ntheta] = msg;

			mb->theta_from = Realloc(mb->theta_from, mb->ntheta + 1, char *);
			mb->theta_to = Realloc(mb->theta_to, mb->ntheta + 1, char *);
			mb->theta_from[mb->ntheta] = Strdup(mb->f_prior[mb->nf][1].from_theta);
			mb->theta_to[mb->ntheta] = Strdup(mb->f_prior[mb->nf][1].to_theta);

			mb->theta[mb->ntheta] = phi_intern;
			mb->theta_map = Realloc(mb->theta_map, mb->ntheta + 1, map_func_tp *);
			mb->theta_map[mb->ntheta] = map_probability;
			mb->theta_map_arg = Realloc(mb->theta_map_arg, mb->ntheta + 1, void *);
			mb->theta_map_arg[mb->ntheta] = NULL;
			mb->ntheta++;
		}
	}
		break;

	case F_RW2DIID:
	{
		tmp = iniparser_getdouble(ini, inla_string_join(secname, "INITIAL0"), G.log_prec_initial);
		if (!mb->f_fixed[mb->nf][0] && mb->mode_use_mode) {
			tmp = mb->theta_file[mb->theta_counter_file++];
			if (mb->mode_fixed) {
				mb->f_fixed[mb->nf][0] = 1;
			}
		}
		_SetInitial(0, tmp);
		HYPER_INIT(log_prec, tmp);
		if (mb->verbose) {
			printf("\t\tinitialise log_precision [%g]\n", tmp);
			printf("\t\tfixed=[%1d]\n", mb->f_fixed[mb->nf][0]);
		}

		tmp = iniparser_getdouble(ini, inla_string_join(secname, "INITIAL1"), 0.0);
		if (!mb->f_fixed[mb->nf][1] && mb->mode_use_mode) {
			tmp = mb->theta_file[mb->theta_counter_file++];
			if (mb->mode_fixed) {
				mb->f_fixed[mb->nf][1] = 1;
			}
		}
		_SetInitial(1, tmp);
		HYPER_INIT(phi_intern, tmp);
		if (mb->verbose) {
			printf("\t\tinitialise phi_intern [%g]\n", tmp);
			printf("\t\tfixed=[%1d]\n", mb->f_fixed[mb->nf][1]);
		}
		mb->f_theta[mb->nf] = Calloc(2, double **);
		mb->f_theta[mb->nf][0] = log_prec;
		mb->f_theta[mb->nf][1] = phi_intern;

		if (!mb->f_fixed[mb->nf][0]) {
			mb->theta = Realloc(mb->theta, mb->ntheta + 1, double **);
			mb->theta_hyperid = Realloc(mb->theta_hyperid, mb->ntheta + 1, char *);
			mb->theta_hyperid[mb->ntheta] = mb->f_prior[mb->nf][0].hyperid;
			mb->theta_tag = Realloc(mb->theta_tag, mb->ntheta + 1, char *);
			mb->theta_tag_userscale = Realloc(mb->theta_tag_userscale, mb->ntheta + 1, char *);
			mb->theta_dir = Realloc(mb->theta_dir, mb->ntheta + 1, char *);
			GMRFLib_sprintf(&msg, "Log precision for %s", mb->f_tag[mb->nf]);
			mb->theta_tag[mb->ntheta] = msg;
			GMRFLib_sprintf(&msg, "Precision for %s", mb->f_tag[mb->nf]);
			mb->theta_tag_userscale[mb->ntheta] = msg;
			GMRFLib_sprintf(&msg, "%s-parameter0", mb->f_dir[mb->nf]);
			mb->theta_dir[mb->ntheta] = msg;

			mb->theta_from = Realloc(mb->theta_from, mb->ntheta + 1, char *);
			mb->theta_to = Realloc(mb->theta_to, mb->ntheta + 1, char *);
			mb->theta_from[mb->ntheta] = Strdup(mb->f_prior[mb->nf][0].from_theta);
			mb->theta_to[mb->ntheta] = Strdup(mb->f_prior[mb->nf][0].to_theta);

			mb->theta[mb->ntheta] = log_prec;
			mb->theta_map = Realloc(mb->theta_map, mb->ntheta + 1, map_func_tp *);
			mb->theta_map[mb->ntheta] = map_precision;
			mb->theta_map_arg = Realloc(mb->theta_map_arg, mb->ntheta + 1, void *);
			mb->theta_map_arg[mb->ntheta] = NULL;
			mb->ntheta++;
		}

		if (!mb->f_fixed[mb->nf][1]) {
			mb->theta = Realloc(mb->theta, mb->ntheta + 1, double **);
			mb->theta_hyperid = Realloc(mb->theta_hyperid, mb->ntheta + 1, char *);
			mb->theta_hyperid[mb->ntheta] = mb->f_prior[mb->nf][1].hyperid;
			mb->theta_tag = Realloc(mb->theta_tag, mb->ntheta + 1, char *);
			mb->theta_tag_userscale = Realloc(mb->theta_tag_userscale, mb->ntheta + 1, char *);
			mb->theta_dir = Realloc(mb->theta_dir, mb->ntheta + 1, char *);
			GMRFLib_sprintf(&msg, "Logit phi for %s", mb->f_tag[mb->nf]);
			mb->theta_tag[mb->ntheta] = msg;
			GMRFLib_sprintf(&msg, "Phi for %s", mb->f_tag[mb->nf]);
			mb->theta_tag_userscale[mb->ntheta] = msg;
			GMRFLib_sprintf(&msg, "%s-parameter1", mb->f_dir[mb->nf]);
			mb->theta_dir[mb->ntheta] = msg;

			mb->theta_from = Realloc(mb->theta_from, mb->ntheta + 1, char *);
			mb->theta_to = Realloc(mb->theta_to, mb->ntheta + 1, char *);
			mb->theta_from[mb->ntheta] = Strdup(mb->f_prior[mb->nf][1].from_theta);
			mb->theta_to[mb->ntheta] = Strdup(mb->f_prior[mb->nf][1].to_theta);

			mb->theta[mb->ntheta] = phi_intern;
			mb->theta_map = Realloc(mb->theta_map, mb->ntheta + 1, map_func_tp *);
			mb->theta_map[mb->ntheta] = map_probability;
			mb->theta_map_arg = Realloc(mb->theta_map_arg, mb->ntheta + 1, void *);
			mb->theta_map_arg[mb->ntheta] = NULL;
			mb->ntheta++;
		}
	}
		break;

	case F_2DIID:
	{
		tmp = iniparser_getdouble(ini, inla_string_join(secname, "INITIAL0"), G.log_prec_initial);
		if (!mb->f_fixed[mb->nf][0] && mb->mode_use_mode) {
			tmp = mb->theta_file[mb->theta_counter_file++];
			if (mb->mode_fixed) {
				mb->f_fixed[mb->nf][0] = 1;
			}
		}
		_SetInitial(0, tmp);
		HYPER_INIT(log_prec0, tmp);
		if (mb->verbose) {
			printf("\t\tinitialise log_precision (first component)[%g]\n", tmp);
			printf("\t\tfixed=[%1d]\n", mb->f_fixed[mb->nf][0]);
		}
		tmp = iniparser_getdouble(ini, inla_string_join(secname, "INITIAL1"), G.log_prec_initial);
		if (!mb->f_fixed[mb->nf][1] && mb->mode_use_mode) {
			tmp = mb->theta_file[mb->theta_counter_file++];
			if (mb->mode_fixed) {
				mb->f_fixed[mb->nf][1] = 1;
			}
		}
		_SetInitial(1, tmp);
		HYPER_INIT(log_prec1, tmp);
		if (mb->verbose) {
			printf("\t\tinitialise log_precision (second component)[%g]\n", tmp);
			printf("\t\tfixed=[%1d]\n", mb->f_fixed[mb->nf][1]);
		}

		mb->f_theta[mb->nf] = Calloc(3, double **);
		mb->f_theta[mb->nf][0] = log_prec0;
		if (!mb->f_fixed[mb->nf][0]) {
			mb->theta = Realloc(mb->theta, mb->ntheta + 1, double **);
			mb->theta_hyperid = Realloc(mb->theta_hyperid, mb->ntheta + 1, char *);
			mb->theta_hyperid[mb->ntheta] = mb->f_prior[mb->nf][0].hyperid;
			mb->theta_tag = Realloc(mb->theta_tag, mb->ntheta + 1, char *);
			mb->theta_tag_userscale = Realloc(mb->theta_tag_userscale, mb->ntheta + 1, char *);
			mb->theta_dir = Realloc(mb->theta_dir, mb->ntheta + 1, char *);
			GMRFLib_sprintf(&msg, "Log precision for %s (first component)", mb->f_tag[mb->nf]);
			mb->theta_tag[mb->ntheta] = msg;
			GMRFLib_sprintf(&msg, "Precision for %s (first component)", mb->f_tag[mb->nf]);
			mb->theta_tag_userscale[mb->ntheta] = msg;
			GMRFLib_sprintf(&msg, "%s-parameter0", mb->f_dir[mb->nf]);
			mb->theta_dir[mb->ntheta] = msg;

			mb->theta_from = Realloc(mb->theta_from, mb->ntheta + 1, char *);
			mb->theta_to = Realloc(mb->theta_to, mb->ntheta + 1, char *);
			mb->theta_from[mb->ntheta] = Strdup(mb->f_prior[mb->nf][0].from_theta);
			mb->theta_to[mb->ntheta] = Strdup(mb->f_prior[mb->nf][0].to_theta);

			mb->theta[mb->ntheta] = log_prec0;
			mb->theta_map = Realloc(mb->theta_map, mb->ntheta + 1, map_func_tp *);
			mb->theta_map[mb->ntheta] = map_precision;
			mb->theta_map_arg = Realloc(mb->theta_map_arg, mb->ntheta + 1, void *);
			mb->theta_map_arg[mb->ntheta] = NULL;
			mb->ntheta++;
		}
		mb->f_theta[mb->nf][1] = log_prec1;
		if (!mb->f_fixed[mb->nf][1]) {
			mb->theta = Realloc(mb->theta, mb->ntheta + 1, double **);
			mb->theta_hyperid = Realloc(mb->theta_hyperid, mb->ntheta + 1, char *);
			mb->theta_hyperid[mb->ntheta] = mb->f_prior[mb->nf][1].hyperid;
			mb->theta_tag = Realloc(mb->theta_tag, mb->ntheta + 1, char *);
			mb->theta_tag_userscale = Realloc(mb->theta_tag_userscale, mb->ntheta + 1, char *);
			mb->theta_dir = Realloc(mb->theta_dir, mb->ntheta + 1, char *);
			GMRFLib_sprintf(&msg, "Log precision for %s (second component)", mb->f_tag[mb->nf]);
			mb->theta_tag[mb->ntheta] = msg;
			GMRFLib_sprintf(&msg, "Precision for %s (second component)", mb->f_tag[mb->nf]);
			mb->theta_tag_userscale[mb->ntheta] = msg;
			GMRFLib_sprintf(&msg, "%s-parameter1", mb->f_dir[mb->nf]);
			mb->theta_dir[mb->ntheta] = msg;

			mb->theta_from = Realloc(mb->theta_from, mb->ntheta + 1, char *);
			mb->theta_to = Realloc(mb->theta_to, mb->ntheta + 1, char *);
			mb->theta_from[mb->ntheta] = Strdup(mb->f_prior[mb->nf][1].from_theta);
			mb->theta_to[mb->ntheta] = Strdup(mb->f_prior[mb->nf][1].to_theta);

			mb->theta[mb->ntheta] = log_prec1;
			mb->theta_map = Realloc(mb->theta_map, mb->ntheta + 1, map_func_tp *);
			mb->theta_map[mb->ntheta] = map_precision;
			mb->theta_map_arg = Realloc(mb->theta_map_arg, mb->ntheta + 1, void *);
			mb->theta_map_arg[mb->ntheta] = NULL;
			mb->ntheta++;
		}
		tmp = iniparser_getdouble(ini, inla_string_join(secname, "INITIAL2"), 0.0);
		if (!mb->f_fixed[mb->nf][2] && mb->mode_use_mode) {
			tmp = mb->theta_file[mb->theta_counter_file++];
			if (mb->mode_fixed) {
				mb->f_fixed[mb->nf][2] = 1;
			}
		}
		_SetInitial(2, tmp);
		HYPER_INIT(rho_intern, tmp);
		if (mb->verbose) {
			printf("\t\tinitialise rho_intern[%g]\n", tmp);
			printf("\t\tfixed=[%1d]\n", mb->f_fixed[mb->nf][2]);
		}

		mb->f_theta[mb->nf][2] = rho_intern;
		if (!mb->f_fixed[mb->nf][2]) {
			mb->theta = Realloc(mb->theta, mb->ntheta + 1, double **);
			mb->theta_hyperid = Realloc(mb->theta_hyperid, mb->ntheta + 1, char *);
			mb->theta_hyperid[mb->ntheta] = mb->f_prior[mb->nf][2].hyperid;
			mb->theta_tag = Realloc(mb->theta_tag, mb->ntheta + 1, char *);
			mb->theta_tag_userscale = Realloc(mb->theta_tag_userscale, mb->ntheta + 1, char *);
			mb->theta_dir = Realloc(mb->theta_dir, mb->ntheta + 1, char *);
			GMRFLib_sprintf(&msg, "Rho_intern for %s", mb->f_tag[mb->nf]);
			mb->theta_tag[mb->ntheta] = msg;
			GMRFLib_sprintf(&msg, "Rho for %s", mb->f_tag[mb->nf]);
			mb->theta_tag_userscale[mb->ntheta] = msg;
			GMRFLib_sprintf(&msg, "%s-parameter2", mb->f_dir[mb->nf]);
			mb->theta_dir[mb->ntheta] = msg;

			mb->theta_from = Realloc(mb->theta_from, mb->ntheta + 1, char *);
			mb->theta_to = Realloc(mb->theta_to, mb->ntheta + 1, char *);
			mb->theta_from[mb->ntheta] = Strdup(mb->f_prior[mb->nf][2].from_theta);
			mb->theta_to[mb->ntheta] = Strdup(mb->f_prior[mb->nf][2].to_theta);

			mb->theta[mb->ntheta] = rho_intern;
			mb->theta_map = Realloc(mb->theta_map, mb->ntheta + 1, map_func_tp *);
			mb->theta_map[mb->ntheta] = map_rho;
			mb->theta_map_arg = Realloc(mb->theta_map_arg, mb->ntheta + 1, void *);
			mb->theta_map_arg[mb->ntheta] = NULL;
			mb->ntheta++;
		}
	}
		break;

	case F_IID1D:
	case F_IID2D:
	case F_IID3D:
	case F_IID4D:
	case F_IID5D:
	{
		int dim = WISHART_DIM(mb->nf);
		assert(dim > 0);

		int n_theta = mb->f_ntheta[mb->nf];
		theta_iidwishart = Calloc(n_theta, double **);
		for (i = 0; i < n_theta; i++) {
			HYPER_NEW(theta_iidwishart[i], 0.0);
		}

		mb->f_theta[mb->nf] = Calloc(n_theta, double **);
		k = 0;
		for (i = 0; i < dim; i++) {
			/*
			 * first get all the precisions 
			 */
			char *init = NULL;

			if (dim == 1) {
				GMRFLib_sprintf(&init, "INITIAL");
			} else {
				GMRFLib_sprintf(&init, "INITIAL%1d", k);
			}

			tmp = iniparser_getdouble(ini, inla_string_join(secname, init), G.log_prec_initial);

			if (!mb->f_fixed[mb->nf][k] && mb->mode_use_mode) {
				tmp = mb->theta_file[mb->theta_counter_file++];
				if (mb->mode_fixed) {
					mb->f_fixed[mb->nf][k] = 1;
				}
			}
			_SetInitial(k, tmp);
			HYPER_INIT(theta_iidwishart[k], tmp);
			if (mb->verbose) {
				printf("\t\tinitialise log_precision (component %d)[%g]\n", k + 1, tmp);
				printf("\t\tfixed=[%1d]\n", mb->f_fixed[mb->nf][k]);
			}

			mb->f_theta[mb->nf][k] = theta_iidwishart[k];

			if (!mb->f_fixed[mb->nf][k]) {
				mb->theta = Realloc(mb->theta, mb->ntheta + 1, double **);
				mb->theta_hyperid = Realloc(mb->theta_hyperid, mb->ntheta + 1, char *);
				mb->theta_hyperid[mb->ntheta] = mb->f_prior[mb->nf][k].hyperid;
				mb->theta_tag = Realloc(mb->theta_tag, mb->ntheta + 1, char *);
				mb->theta_tag_userscale = Realloc(mb->theta_tag_userscale, mb->ntheta + 1, char *);
				mb->theta_dir = Realloc(mb->theta_dir, mb->ntheta + 1, char *);
				GMRFLib_sprintf(&msg, "Log precision for %s (component %1d)", mb->f_tag[mb->nf], k + 1);
				mb->theta_tag[mb->ntheta] = msg;
				GMRFLib_sprintf(&msg, "Precision for %s (component %1d)", mb->f_tag[mb->nf], k + 1);
				mb->theta_tag_userscale[mb->ntheta] = msg;
				GMRFLib_sprintf(&msg, "%s-parameter%1d", mb->f_dir[mb->nf], k);
				mb->theta_dir[mb->ntheta] = msg;

				mb->theta_from = Realloc(mb->theta_from, mb->ntheta + 1, char *);
				mb->theta_to = Realloc(mb->theta_to, mb->ntheta + 1, char *);
				mb->theta_from[mb->ntheta] = Strdup(mb->f_prior[mb->nf][k].from_theta);
				mb->theta_to[mb->ntheta] = Strdup(mb->f_prior[mb->nf][k].to_theta);

				mb->theta[mb->ntheta] = theta_iidwishart[k];
				mb->theta_map = Realloc(mb->theta_map, mb->ntheta + 1, map_func_tp *);
				mb->theta_map[mb->ntheta] = map_precision;
				mb->theta_map_arg = Realloc(mb->theta_map_arg, mb->ntheta + 1, void *);
				mb->theta_map_arg[mb->ntheta] = NULL;
				mb->ntheta++;
			}

			k++;
		}
		for (i = 0; i < dim; i++) {
			for (j = i + 1; j < dim; j++) {
				/*
				 * all the correlations 
				 */
				char *init = NULL;
				GMRFLib_sprintf(&init, "INITIAL%1d", k);
				tmp = iniparser_getdouble(ini, inla_string_join(secname, init), 0.0);

				if (!mb->f_fixed[mb->nf][k] && mb->mode_use_mode) {
					tmp = mb->theta_file[mb->theta_counter_file++];
					if (mb->mode_fixed) {
						mb->f_fixed[mb->nf][k] = 1;
					}
				}
				_SetInitial(k, tmp);
				HYPER_INIT(theta_iidwishart[k], tmp);
				if (mb->verbose) {
					printf("\t\tinitialise rho_internal%1d:%1d [%g]\n", i + 1, j + 1, tmp);
					printf("\t\tfixed=[%1d]\n", mb->f_fixed[mb->nf][k]);
				}

				mb->f_theta[mb->nf][k] = theta_iidwishart[k];

				if (!mb->f_fixed[mb->nf][k]) {
					mb->theta = Realloc(mb->theta, mb->ntheta + 1, double **);
					mb->theta_hyperid = Realloc(mb->theta_hyperid, mb->ntheta + 1, char *);
					mb->theta_hyperid[mb->ntheta] = mb->f_prior[mb->nf][k].hyperid;
					mb->theta_tag = Realloc(mb->theta_tag, mb->ntheta + 1, char *);
					mb->theta_tag_userscale = Realloc(mb->theta_tag_userscale, mb->ntheta + 1, char *);
					mb->theta_dir = Realloc(mb->theta_dir, mb->ntheta + 1, char *);
					GMRFLib_sprintf(&msg, "Rho_internal%1d:%1d for %s", i + 1, j + 1, mb->f_tag[mb->nf]);
					mb->theta_tag[mb->ntheta] = msg;
					GMRFLib_sprintf(&msg, "Rho%1d:%1d for %s", i + 1, j + 1, mb->f_tag[mb->nf]);
					mb->theta_tag_userscale[mb->ntheta] = msg;
					GMRFLib_sprintf(&msg, "%s-parameter%1d", mb->f_dir[mb->nf], k);
					mb->theta_dir[mb->ntheta] = msg;

					mb->theta_from = Realloc(mb->theta_from, mb->ntheta + 1, char *);
					mb->theta_to = Realloc(mb->theta_to, mb->ntheta + 1, char *);
					mb->theta_from[mb->ntheta] = Strdup(mb->f_prior[mb->nf][k].from_theta);
					mb->theta_to[mb->ntheta] = Strdup(mb->f_prior[mb->nf][k].to_theta);

					mb->theta[mb->ntheta] = theta_iidwishart[k];
					mb->theta_map = Realloc(mb->theta_map, mb->ntheta + 1, map_func_tp *);
					mb->theta_map[mb->ntheta] = map_rho;
					mb->theta_map_arg = Realloc(mb->theta_map_arg, mb->ntheta + 1, void *);
					mb->theta_map_arg[mb->ntheta] = NULL;
					mb->ntheta++;
				}
				k++;
			}
		}
		assert(k == n_theta);
	}
		break;

	case F_IIDKD:
	{
		int dim = mb->f_order[mb->nf];
		int n_theta = mb->f_ntheta[mb->nf];

		theta_iidwishart = Calloc(n_theta, double **);
		for (k = 0; k < n_theta; k++) {
			HYPER_NEW(theta_iidwishart[k], 0.0);
		}

		for (k = 0; k < INLA_WISHARTK_NTHETA(INLA_WISHARTK_KMAX); k++) {
			char *txt = NULL;
			GMRFLib_sprintf(&txt, "INITIAL%1d", k);
			iniparser_getdouble(ini, inla_string_join(secname, txt), 0.0);
			GMRFLib_sprintf(&txt, "FIXED%1d", k);
			iniparser_getint(ini, inla_string_join(secname, txt), 0);
		}

		mb->f_theta[mb->nf] = Calloc(n_theta, double **);
		for (k = 0; k < n_theta; k++) {
			char *init = NULL;
			GMRFLib_sprintf(&init, "INITIAL%1d", k);
			if (k < dim) {
				tmp = iniparser_getdouble(ini, inla_string_join(secname, init), INLA_SPECIAL_NUMBER);
				if (INLA_IS_SPECIAL(tmp)) {
					tmp = G.log_prec_initial / 2.0;
				}
			} else {
				tmp = iniparser_getdouble(ini, inla_string_join(secname, init), INLA_SPECIAL_NUMBER);
				if (INLA_IS_SPECIAL(tmp)) {
					tmp = 0.0;
				}
			}
			if (!mb->f_fixed[mb->nf][k] && mb->mode_use_mode) {
				tmp = mb->theta_file[mb->theta_counter_file++];
				if (mb->mode_fixed) {
					mb->f_fixed[mb->nf][k] = 1;
				}
			}
			_SetInitial(k, tmp);
			HYPER_INIT(theta_iidwishart[k], tmp);
			if (mb->verbose) {
				printf("\t\tinitialise theta%1d=[%g]\n", k + 1, tmp);
				printf("\t\tfixed=[%1d]\n", mb->f_fixed[mb->nf][k]);
			}
			mb->f_theta[mb->nf][k] = theta_iidwishart[k];

			if (!mb->f_fixed[mb->nf][k]) {
				mb->theta = Realloc(mb->theta, mb->ntheta + 1, double **);
				mb->theta_hyperid = Realloc(mb->theta_hyperid, mb->ntheta + 1, char *);
				mb->theta_hyperid[mb->ntheta] = mb->f_prior[mb->nf][k].hyperid;
				mb->theta_tag = Realloc(mb->theta_tag, mb->ntheta + 1, char *);
				mb->theta_tag_userscale = Realloc(mb->theta_tag_userscale, mb->ntheta + 1, char *);
				mb->theta_dir = Realloc(mb->theta_dir, mb->ntheta + 1, char *);
				GMRFLib_sprintf(&msg, "Theta%1d for %s", k + 1, mb->f_tag[mb->nf]);
				mb->theta_tag[mb->ntheta] = msg;
				mb->theta_tag_userscale[mb->ntheta] = msg;
				GMRFLib_sprintf(&msg, "%s-parameter%1d", mb->f_dir[mb->nf], k);
				mb->theta_dir[mb->ntheta] = msg;

				mb->theta_from = Realloc(mb->theta_from, mb->ntheta + 1, char *);
				mb->theta_to = Realloc(mb->theta_to, mb->ntheta + 1, char *);
				mb->theta_from[mb->ntheta] = Strdup(mb->f_prior[mb->nf][k].from_theta);
				mb->theta_to[mb->ntheta] = Strdup(mb->f_prior[mb->nf][k].to_theta);

				mb->theta[mb->ntheta] = theta_iidwishart[k];
				mb->theta_map = Realloc(mb->theta_map, mb->ntheta + 1, map_func_tp *);
				mb->theta_map[mb->ntheta] = map_identity;
				mb->theta_map_arg = Realloc(mb->theta_map_arg, mb->ntheta + 1, void *);
				mb->theta_map_arg[mb->ntheta] = NULL;
				mb->ntheta++;
			}
		}
		assert(k == n_theta);
	}
		break;

	case F_INTSLOPE:
	{
		int dim = 2;
		int n_theta = inla_iid_wishart_nparam(dim);
		theta_iidwishart = Calloc(n_theta, double **);
		for (i = 0; i < n_theta; i++) {
			HYPER_NEW(theta_iidwishart[i], 0.0);
		}

		mb->f_theta[mb->nf] = Calloc(mb->f_ntheta[mb->nf], double **);
		k = 0;
		for (i = 0; i < dim; i++) {
			/*
			 * first get all the precisions 
			 */
			char *init = NULL;

			GMRFLib_sprintf(&init, "INITIAL%1d", k);
			tmp = iniparser_getdouble(ini, inla_string_join(secname, init), G.log_prec_initial);
			if (!mb->f_fixed[mb->nf][k] && mb->mode_use_mode) {
				tmp = mb->theta_file[mb->theta_counter_file++];
				if (mb->mode_fixed) {
					mb->f_fixed[mb->nf][k] = 1;
				}
			}
			_SetInitial(k, tmp);
			HYPER_INIT(theta_iidwishart[k], tmp);
			if (mb->verbose) {
				printf("\t\tinitialise log_precision (component %d)[%g]\n", k + 1, tmp);
				printf("\t\tfixed=[%1d]\n", mb->f_fixed[mb->nf][k]);
			}

			mb->f_theta[mb->nf][k] = theta_iidwishart[k];

			if (!mb->f_fixed[mb->nf][k]) {
				mb->theta = Realloc(mb->theta, mb->ntheta + 1, double **);
				mb->theta_hyperid = Realloc(mb->theta_hyperid, mb->ntheta + 1, char *);
				mb->theta_hyperid[mb->ntheta] = mb->f_prior[mb->nf][k].hyperid;
				mb->theta_tag = Realloc(mb->theta_tag, mb->ntheta + 1, char *);
				mb->theta_tag_userscale = Realloc(mb->theta_tag_userscale, mb->ntheta + 1, char *);
				mb->theta_dir = Realloc(mb->theta_dir, mb->ntheta + 1, char *);
				GMRFLib_sprintf(&msg, "Log precision for %s (component %1d)", mb->f_tag[mb->nf], k + 1);
				mb->theta_tag[mb->ntheta] = msg;
				GMRFLib_sprintf(&msg, "Precision for %s (component %1d)", mb->f_tag[mb->nf], k + 1);
				mb->theta_tag_userscale[mb->ntheta] = msg;
				GMRFLib_sprintf(&msg, "%s-parameter%1d", mb->f_dir[mb->nf], k);
				mb->theta_dir[mb->ntheta] = msg;

				mb->theta_from = Realloc(mb->theta_from, mb->ntheta + 1, char *);
				mb->theta_to = Realloc(mb->theta_to, mb->ntheta + 1, char *);
				mb->theta_from[mb->ntheta] = Strdup(mb->f_prior[mb->nf][k].from_theta);
				mb->theta_to[mb->ntheta] = Strdup(mb->f_prior[mb->nf][k].to_theta);

				mb->theta[mb->ntheta] = theta_iidwishart[k];
				mb->theta_map = Realloc(mb->theta_map, mb->ntheta + 1, map_func_tp *);
				mb->theta_map[mb->ntheta] = map_precision;
				mb->theta_map_arg = Realloc(mb->theta_map_arg, mb->ntheta + 1, void *);
				mb->theta_map_arg[mb->ntheta] = NULL;
				mb->ntheta++;
			}

			k++;
		}
		for (i = 0; i < dim; i++) {
			for (j = i + 1; j < dim; j++) {
				/*
				 * all the correlations 
				 */
				char *init = NULL;
				GMRFLib_sprintf(&init, "INITIAL%1d", k);
				tmp = iniparser_getdouble(ini, inla_string_join(secname, init), 0.0);

				if (!mb->f_fixed[mb->nf][k] && mb->mode_use_mode) {
					tmp = mb->theta_file[mb->theta_counter_file++];
					if (mb->mode_fixed) {
						mb->f_fixed[mb->nf][k] = 1;
					}
				}
				_SetInitial(k, tmp);
				HYPER_INIT(theta_iidwishart[k], tmp);
				if (mb->verbose) {
					printf("\t\tinitialise rho_internal%1d:%1d [%g]\n", i + 1, j + 1, tmp);
					printf("\t\tfixed=[%1d]\n", mb->f_fixed[mb->nf][k]);
				}

				mb->f_theta[mb->nf][k] = theta_iidwishart[k];

				if (!mb->f_fixed[mb->nf][k]) {
					mb->theta = Realloc(mb->theta, mb->ntheta + 1, double **);
					mb->theta_hyperid = Realloc(mb->theta_hyperid, mb->ntheta + 1, char *);
					mb->theta_hyperid[mb->ntheta] = mb->f_prior[mb->nf][k].hyperid;
					mb->theta_tag = Realloc(mb->theta_tag, mb->ntheta + 1, char *);
					mb->theta_tag_userscale = Realloc(mb->theta_tag_userscale, mb->ntheta + 1, char *);
					mb->theta_dir = Realloc(mb->theta_dir, mb->ntheta + 1, char *);
					GMRFLib_sprintf(&msg, "Rho_internal%1d:%1d for %s", i + 1, j + 1, mb->f_tag[mb->nf]);
					mb->theta_tag[mb->ntheta] = msg;
					GMRFLib_sprintf(&msg, "Rho%1d:%1d for %s", i + 1, j + 1, mb->f_tag[mb->nf]);
					mb->theta_tag_userscale[mb->ntheta] = msg;
					GMRFLib_sprintf(&msg, "%s-parameter%1d", mb->f_dir[mb->nf], k);
					mb->theta_dir[mb->ntheta] = msg;

					mb->theta_from = Realloc(mb->theta_from, mb->ntheta + 1, char *);
					mb->theta_to = Realloc(mb->theta_to, mb->ntheta + 1, char *);
					mb->theta_from[mb->ntheta] = Strdup(mb->f_prior[mb->nf][k].from_theta);
					mb->theta_to[mb->ntheta] = Strdup(mb->f_prior[mb->nf][k].to_theta);

					mb->theta[mb->ntheta] = theta_iidwishart[k];
					mb->theta_map = Realloc(mb->theta_map, mb->ntheta + 1, map_func_tp *);
					mb->theta_map[mb->ntheta] = map_rho;
					mb->theta_map_arg = Realloc(mb->theta_map_arg, mb->ntheta + 1, void *);
					mb->theta_map_arg[mb->ntheta] = NULL;
					mb->ntheta++;
				}
				k++;
			}
		}

		intslope_gamma = Calloc(INTSLOPE_MAXTHETA, double **);
		for (i = 0; i < INTSLOPE_MAXTHETA; i++) {
			HYPER_NEW(intslope_gamma[i], 1.0);
		}

		int kk;
		for (i = 0; i < INTSLOPE_MAXTHETA; i++) {
			char *init = NULL;

			kk = k - n_theta;
			if (i >= nstrata) {

				P(nstrata);
				mb->f_fixed[mb->nf][k] = 1;
			}
			GMRFLib_sprintf(&init, "INITIAL%1d", k);
			tmp = iniparser_getdouble(ini, inla_string_join(secname, init), 1.0);

			if (!mb->f_fixed[mb->nf][k] && mb->mode_use_mode) {
				tmp = mb->theta_file[mb->theta_counter_file++];
				if (mb->mode_fixed) {
					mb->f_fixed[mb->nf][k] = 1;
				}
			}
			_SetInitial(k, tmp);
			HYPER_INIT(intslope_gamma[kk], tmp);
			if (mb->verbose) {
				printf("\t\tinitialise gamma[%1d] = [%g]\n", kk, tmp);
				printf("\t\tfixed=[%1d]\n", mb->f_fixed[mb->nf][k]);
			}
			mb->f_theta[mb->nf][k] = intslope_gamma[kk];
			if (!mb->f_fixed[mb->nf][k]) {
				mb->theta = Realloc(mb->theta, mb->ntheta + 1, double **);
				mb->theta_hyperid = Realloc(mb->theta_hyperid, mb->ntheta + 1, char *);
				mb->theta_hyperid[mb->ntheta] = mb->f_prior[mb->nf][k].hyperid;
				mb->theta_tag = Realloc(mb->theta_tag, mb->ntheta + 1, char *);
				mb->theta_tag_userscale = Realloc(mb->theta_tag_userscale, mb->ntheta + 1, char *);
				mb->theta_dir = Realloc(mb->theta_dir, mb->ntheta + 1, char *);
				GMRFLib_sprintf(&msg, "Gamma%1d for %s", kk + 1, mb->f_tag[mb->nf]);
				mb->theta_tag[mb->ntheta] = msg;
				GMRFLib_sprintf(&msg, "Gamma%1d for %s", kk + 1, mb->f_tag[mb->nf]);
				mb->theta_tag_userscale[mb->ntheta] = msg;
				GMRFLib_sprintf(&msg, "%s-parameter%1d", mb->f_dir[mb->nf], k);
				mb->theta_dir[mb->ntheta] = msg;

				mb->theta_from = Realloc(mb->theta_from, mb->ntheta + 1, char *);
				mb->theta_to = Realloc(mb->theta_to, mb->ntheta + 1, char *);
				mb->theta_from[mb->ntheta] = Strdup(mb->f_prior[mb->nf][k].from_theta);
				mb->theta_to[mb->ntheta] = Strdup(mb->f_prior[mb->nf][k].to_theta);

				mb->theta[mb->ntheta] = intslope_gamma[kk];
				mb->theta_map = Realloc(mb->theta_map, mb->ntheta + 1, map_func_tp *);
				mb->theta_map[mb->ntheta] = map_identity;
				mb->theta_map_arg = Realloc(mb->theta_map_arg, mb->ntheta + 1, void *);
				mb->theta_map_arg[mb->ntheta] = NULL;
				mb->ntheta++;
			}
			k++;
		}

		assert(k == n_theta + INTSLOPE_MAXTHETA);
	}
		break;

	case F_MATERN2D:
	{
		itmp = iniparser_getint(ini, inla_string_join(secname, "NU"), 1);
		mb->f_nu[mb->nf] = itmp;
		if (mb->verbose) {
			printf("\t\tnu = [%1d]\n", itmp);
		}
		tmp = iniparser_getdouble(ini, inla_string_join(secname, "INITIAL0"), G.log_prec_initial);
		if (!mb->f_fixed[mb->nf][0] && mb->mode_use_mode) {
			tmp = mb->theta_file[mb->theta_counter_file++];
			if (mb->mode_fixed) {
				mb->f_fixed[mb->nf][0] = 1;
			}
		}
		_SetInitial(0, tmp);
		HYPER_INIT(log_prec, tmp);
		if (mb->verbose) {
			printf("\t\tinitialise log_precision[%g]\n", tmp);
			printf("\t\tfixed=[%1d]\n", mb->f_fixed[mb->nf][0]);
		}

		mb->f_theta[mb->nf] = Calloc(2, double **);
		mb->f_theta[mb->nf][0] = log_prec;
		if (!mb->f_fixed[mb->nf][0]) {
			/*
			 * add this \theta 
			 */
			mb->theta = Realloc(mb->theta, mb->ntheta + 1, double **);
			mb->theta_hyperid = Realloc(mb->theta_hyperid, mb->ntheta + 1, char *);
			mb->theta_hyperid[mb->ntheta] = mb->f_prior[mb->nf][0].hyperid;
			mb->theta_tag = Realloc(mb->theta_tag, mb->ntheta + 1, char *);
			mb->theta_tag_userscale = Realloc(mb->theta_tag_userscale, mb->ntheta + 1, char *);
			mb->theta_dir = Realloc(mb->theta_dir, mb->ntheta + 1, char *);
			GMRFLib_sprintf(&msg, "Log precision for %s", (secname ? secname : mb->f_tag[mb->nf]));
			mb->theta_tag[mb->ntheta] = msg;
			GMRFLib_sprintf(&msg, "Precision for %s", (secname ? secname : mb->f_tag[mb->nf]));
			mb->theta_tag_userscale[mb->ntheta] = msg;
			GMRFLib_sprintf(&msg, "%s-parameter0", mb->f_dir[mb->nf]);
			mb->theta_dir[mb->ntheta] = msg;

			mb->theta_from = Realloc(mb->theta_from, mb->ntheta + 1, char *);
			mb->theta_to = Realloc(mb->theta_to, mb->ntheta + 1, char *);
			mb->theta_from[mb->ntheta] = Strdup(mb->f_prior[mb->nf][0].from_theta);
			mb->theta_to[mb->ntheta] = Strdup(mb->f_prior[mb->nf][0].to_theta);

			mb->theta[mb->ntheta] = log_prec;
			mb->theta_map = Realloc(mb->theta_map, mb->ntheta + 1, map_func_tp *);
			mb->theta_map[mb->ntheta] = map_precision;
			mb->theta_map_arg = Realloc(mb->theta_map_arg, mb->ntheta + 1, void *);
			mb->theta_map_arg[mb->ntheta] = NULL;
			mb->ntheta++;
		}

		tmp = iniparser_getdouble(ini, inla_string_join(secname, "INITIAL1"), 2.0);
		if (!mb->f_fixed[mb->nf][1] && mb->mode_use_mode) {
			tmp = mb->theta_file[mb->theta_counter_file++];
			if (mb->mode_fixed) {
				mb->f_fixed[mb->nf][1] = 1;
			}
		}
		_SetInitial(1, tmp);
		HYPER_INIT(range_intern, tmp);
		if (mb->verbose) {
			printf("\t\tinitialise range_intern[%g]\n", tmp);
			printf("\t\tfixed=[%1d]\n", mb->f_fixed[mb->nf][1]);
		}

		mb->f_theta[mb->nf][1] = range_intern;
		if (!mb->f_fixed[mb->nf][1]) {
			/*
			 * add this \theta 
			 */
			mb->theta = Realloc(mb->theta, mb->ntheta + 1, double **);
			mb->theta_hyperid = Realloc(mb->theta_hyperid, mb->ntheta + 1, char *);
			mb->theta_hyperid[mb->ntheta] = mb->f_prior[mb->nf][1].hyperid;
			mb->theta_tag = Realloc(mb->theta_tag, mb->ntheta + 1, char *);
			mb->theta_tag_userscale = Realloc(mb->theta_tag_userscale, mb->ntheta + 1, char *);
			mb->theta_dir = Realloc(mb->theta_dir, mb->ntheta + 1, char *);
			GMRFLib_sprintf(&msg, "Range_intern for %s", (secname ? secname : mb->f_tag[mb->nf]));
			mb->theta_tag[mb->ntheta] = msg;
			GMRFLib_sprintf(&msg, "Range for %s", (secname ? secname : mb->f_tag[mb->nf]));
			mb->theta_tag_userscale[mb->ntheta] = msg;
			GMRFLib_sprintf(&msg, "%s-parameter1", mb->f_dir[mb->nf]);
			mb->theta_dir[mb->ntheta] = msg;

			mb->theta_from = Realloc(mb->theta_from, mb->ntheta + 1, char *);
			mb->theta_to = Realloc(mb->theta_to, mb->ntheta + 1, char *);
			mb->theta_from[mb->ntheta] = Strdup(mb->f_prior[mb->nf][1].from_theta);
			mb->theta_to[mb->ntheta] = Strdup(mb->f_prior[mb->nf][1].to_theta);

			mb->theta[mb->ntheta] = range_intern;
			mb->theta_map = Realloc(mb->theta_map, mb->ntheta + 1, map_func_tp *);
			mb->theta_map[mb->ntheta] = map_range;
			mb->theta_map_arg = Realloc(mb->theta_map_arg, mb->ntheta + 1, void *);
			mb->theta_map_arg[mb->ntheta] = NULL;
			mb->ntheta++;
		}
	}
		break;

	case F_DMATERN:
	{
		tmp = iniparser_getdouble(ini, inla_string_join(secname, "INITIAL0"), G.log_prec_initial);
		if (!mb->f_fixed[mb->nf][0] && mb->mode_use_mode) {
			tmp = mb->theta_file[mb->theta_counter_file++];
			if (mb->mode_fixed) {
				mb->f_fixed[mb->nf][0] = 1;
			}
		}
		_SetInitial(0, tmp);
		HYPER_INIT(log_prec, tmp);
		if (mb->verbose) {
			printf("\t\tinitialise log_precision[%g]\n", tmp);
			printf("\t\tfixed=[%1d]\n", mb->f_fixed[mb->nf][0]);
		}

		mb->f_theta[mb->nf] = Calloc(3, double **);
		mb->f_theta[mb->nf][0] = log_prec;
		if (!mb->f_fixed[mb->nf][0]) {
			/*
			 * add this \theta 
			 */
			mb->theta = Realloc(mb->theta, mb->ntheta + 1, double **);
			mb->theta_hyperid = Realloc(mb->theta_hyperid, mb->ntheta + 1, char *);
			mb->theta_hyperid[mb->ntheta] = mb->f_prior[mb->nf][0].hyperid;
			mb->theta_tag = Realloc(mb->theta_tag, mb->ntheta + 1, char *);
			mb->theta_tag_userscale = Realloc(mb->theta_tag_userscale, mb->ntheta + 1, char *);
			mb->theta_dir = Realloc(mb->theta_dir, mb->ntheta + 1, char *);
			GMRFLib_sprintf(&msg, "Log precision for %s", (secname ? secname : mb->f_tag[mb->nf]));
			mb->theta_tag[mb->ntheta] = msg;
			GMRFLib_sprintf(&msg, "Precision for %s", (secname ? secname : mb->f_tag[mb->nf]));
			mb->theta_tag_userscale[mb->ntheta] = msg;
			GMRFLib_sprintf(&msg, "%s-parameter0", mb->f_dir[mb->nf]);
			mb->theta_dir[mb->ntheta] = msg;

			mb->theta_from = Realloc(mb->theta_from, mb->ntheta + 1, char *);
			mb->theta_to = Realloc(mb->theta_to, mb->ntheta + 1, char *);
			mb->theta_from[mb->ntheta] = Strdup(mb->f_prior[mb->nf][0].from_theta);
			mb->theta_to[mb->ntheta] = Strdup(mb->f_prior[mb->nf][0].to_theta);

			mb->theta[mb->ntheta] = log_prec;
			mb->theta_map = Realloc(mb->theta_map, mb->ntheta + 1, map_func_tp *);
			mb->theta_map[mb->ntheta] = map_precision;
			mb->theta_map_arg = Realloc(mb->theta_map_arg, mb->ntheta + 1, void *);
			mb->theta_map_arg[mb->ntheta] = NULL;
			mb->ntheta++;
		}

		tmp = iniparser_getdouble(ini, inla_string_join(secname, "INITIAL1"), 2.0);
		if (!mb->f_fixed[mb->nf][1] && mb->mode_use_mode) {
			tmp = mb->theta_file[mb->theta_counter_file++];
			if (mb->mode_fixed) {
				mb->f_fixed[mb->nf][1] = 1;
			}
		}
		_SetInitial(1, tmp);
		HYPER_INIT(range_intern, tmp);
		if (mb->verbose) {
			printf("\t\tinitialise range_intern[%g]\n", tmp);
			printf("\t\tfixed=[%1d]\n", mb->f_fixed[mb->nf][1]);
		}

		mb->f_theta[mb->nf][1] = range_intern;
		if (!mb->f_fixed[mb->nf][1]) {
			/*
			 * add this \theta 
			 */
			mb->theta = Realloc(mb->theta, mb->ntheta + 1, double **);
			mb->theta_hyperid = Realloc(mb->theta_hyperid, mb->ntheta + 1, char *);
			mb->theta_hyperid[mb->ntheta] = mb->f_prior[mb->nf][1].hyperid;
			mb->theta_tag = Realloc(mb->theta_tag, mb->ntheta + 1, char *);
			mb->theta_tag_userscale = Realloc(mb->theta_tag_userscale, mb->ntheta + 1, char *);
			mb->theta_dir = Realloc(mb->theta_dir, mb->ntheta + 1, char *);
			GMRFLib_sprintf(&msg, "log range for %s", (secname ? secname : mb->f_tag[mb->nf]));
			mb->theta_tag[mb->ntheta] = msg;
			GMRFLib_sprintf(&msg, "Range for %s", (secname ? secname : mb->f_tag[mb->nf]));
			mb->theta_tag_userscale[mb->ntheta] = msg;
			GMRFLib_sprintf(&msg, "%s-parameter1", mb->f_dir[mb->nf]);
			mb->theta_dir[mb->ntheta] = msg;

			mb->theta_from = Realloc(mb->theta_from, mb->ntheta + 1, char *);
			mb->theta_to = Realloc(mb->theta_to, mb->ntheta + 1, char *);
			mb->theta_from[mb->ntheta] = Strdup(mb->f_prior[mb->nf][1].from_theta);
			mb->theta_to[mb->ntheta] = Strdup(mb->f_prior[mb->nf][1].to_theta);

			mb->theta[mb->ntheta] = range_intern;
			mb->theta_map = Realloc(mb->theta_map, mb->ntheta + 1, map_func_tp *);
			mb->theta_map[mb->ntheta] = map_range;
			mb->theta_map_arg = Realloc(mb->theta_map_arg, mb->ntheta + 1, void *);
			mb->theta_map_arg[mb->ntheta] = NULL;
			mb->ntheta++;
		}

		tmp = iniparser_getdouble(ini, inla_string_join(secname, "INITIAL2"), log(0.5));
		if (!mb->f_fixed[mb->nf][2] && mb->mode_use_mode) {
			tmp = mb->theta_file[mb->theta_counter_file++];
			if (mb->mode_fixed) {
				mb->f_fixed[mb->nf][2] = 1;
			}
		}
		_SetInitial(2, tmp);
		HYPER_INIT(nu_intern, tmp);
		if (mb->verbose) {
			printf("\t\tinitialise nu_intern[%g]\n", tmp);
			printf("\t\tfixed=[%1d]\n", mb->f_fixed[mb->nf][2]);
		}

		mb->f_theta[mb->nf][2] = nu_intern;
		if (!mb->f_fixed[mb->nf][2]) {
			/*
			 * add this \theta 
			 */
			mb->theta = Realloc(mb->theta, mb->ntheta + 1, double **);
			mb->theta_hyperid = Realloc(mb->theta_hyperid, mb->ntheta + 1, char *);
			mb->theta_hyperid[mb->ntheta] = mb->f_prior[mb->nf][1].hyperid;
			mb->theta_tag = Realloc(mb->theta_tag, mb->ntheta + 1, char *);
			mb->theta_tag_userscale = Realloc(mb->theta_tag_userscale, mb->ntheta + 1, char *);
			mb->theta_dir = Realloc(mb->theta_dir, mb->ntheta + 1, char *);
			GMRFLib_sprintf(&msg, "log nu for %s", (secname ? secname : mb->f_tag[mb->nf]));
			mb->theta_tag[mb->ntheta] = msg;
			GMRFLib_sprintf(&msg, "Nu for %s", (secname ? secname : mb->f_tag[mb->nf]));
			mb->theta_tag_userscale[mb->ntheta] = msg;
			GMRFLib_sprintf(&msg, "%s-parameter1", mb->f_dir[mb->nf]);
			mb->theta_dir[mb->ntheta] = msg;

			mb->theta_from = Realloc(mb->theta_from, mb->ntheta + 1, char *);
			mb->theta_to = Realloc(mb->theta_to, mb->ntheta + 1, char *);
			mb->theta_from[mb->ntheta] = Strdup(mb->f_prior[mb->nf][1].from_theta);
			mb->theta_to[mb->ntheta] = Strdup(mb->f_prior[mb->nf][1].to_theta);

			mb->theta[mb->ntheta] = nu_intern;
			mb->theta_map = Realloc(mb->theta_map, mb->ntheta + 1, map_func_tp *);
			mb->theta_map[mb->ntheta] = map_exp;
			mb->theta_map_arg = Realloc(mb->theta_map_arg, mb->ntheta + 1, void *);
			mb->theta_map_arg[mb->ntheta] = NULL;
			mb->ntheta++;
		}
	}
		break;

	default:
		abort();
	}

	/*
	 ***
	 */

	switch (mb->f_id[mb->nf]) {
	case F_GENERIC0:
	{
		mb->f_N[mb->nf] = mb->f_n[mb->nf];
	}
		break;

	case F_GENERIC1:
	{
		mb->f_N[mb->nf] = mb->f_n[mb->nf];
	}
		break;

	case F_GENERIC2:
	{
		assert(mb->f_N[mb->nf] == 2 * mb->f_n[mb->nf]);
	}
		break;

	case F_GENERIC3:
	{
		assert(mb->f_N[mb->nf] == mb->f_n[mb->nf]);
	}
		break;

	case F_COPY:
	{
		/*
		 * to be filled in later 
		 */
		mb->f_Qfunc[mb->nf] = NULL;
		mb->f_Qfunc_arg[mb->nf] = NULL;
		mb->f_rankdef[mb->nf] = 0.0;
		mb->f_N[mb->nf] = mb->f_n[mb->nf] = -1;
		mb->f_id[mb->nf] = F_COPY;
	}
		break;

	case F_SCOPY:
	{
		/*
		 * these are to be filled in later 
		 */
		mb->f_Qfunc[mb->nf] = NULL;
		mb->f_Qfunc_arg[mb->nf] = NULL;
		mb->f_rankdef[mb->nf] = 0.0;
		mb->f_N[mb->nf] = mb->f_n[mb->nf] = -1;
		mb->f_id[mb->nf] = F_SCOPY;

		// these we can do now
		inla_scopy_arg_tp *def = Calloc(1, inla_scopy_arg_tp);
		def->nbeta = nbeta;
		def->loc_beta = loc;
		def->loc_len = loc[nbeta - 1] - loc[0];
		def->loc_mid = (loc[nbeta - 1] + loc[0]) / 2.0;
		def->cov_beta = cov;
		def->betas = betas;
		def->precision = mb->f_precision[mb->nf];
		def->W = W;
		assert(W);

		def->cache00 = Calloc(GMRFLib_CACHE_LEN(), inla_scopy_cache_tp *);
		for (i = 0; i < GMRFLib_CACHE_LEN(); i++) {
			def->cache00[i] = Calloc(1, inla_scopy_cache_tp);
			def->cache00[i]->betas = Calloc(nbeta, double);
			def->cache00[i]->betas[0] = GMRFLib_uniform();
			def->cache00[i]->betas_tmp = Calloc(nbeta, double);
			def->cache00[i]->splinefun = NULL;
		}

		def->cache01 = Calloc(GMRFLib_CACHE_LEN(), inla_scopy_cache_tp *);
		for (i = 0; i < GMRFLib_CACHE_LEN(); i++) {
			def->cache01[i] = Calloc(1, inla_scopy_cache_tp);
			def->cache01[i]->betas = Calloc(nbeta, double);
			def->cache01[i]->betas[0] = GMRFLib_uniform();
			def->cache01[i]->betas_tmp = Calloc(nbeta, double);
			def->cache01[i]->splinefun = NULL;
		}

		GMRFLib_rwdef_tp *rwdef = Calloc(1, GMRFLib_rwdef_tp);
		rwdef->n = nbeta;
		rwdef->order = 2;

		// we'll add the prior_prec_betas later in the extra() function
		HYPER_NEW(rwdef->log_prec_omp, 1.0);
		GMRFLib_make_rw_graph(&(def->graph_prior), rwdef);
		GMRFLib_rw_scale(thread_id, (void *) rwdef);
		if (mb->verbose) {
			printf("\t\tscale.model: prec_scale[%g]\n", rwdef->prec_scale[0]);
		}
		def->Qfunc_prior = GMRFLib_rw;
		def->Qfunc_arg_prior = (void *) rwdef;
		def->rwdef = rwdef;

		// lets pass 'def' here, maybe chance later
		mb->f_Qfunc_orig[mb->nf] = NULL;
		mb->f_Qfunc_arg_orig[mb->nf] = (void *) def;
	}
		break;

	case F_CLINEAR:
	{
		inla_clinear_tp *def = Calloc(1, inla_clinear_tp);
		def->beta = beta;
		def->beta_arg = mb->f_theta_map_arg[mb->nf][0];
		def->precision = mb->f_precision[mb->nf];
		def->x = mb->f_locations[mb->nf];

		GMRFLib_graph_mk_linear(&(mb->f_graph[mb->nf]), mb->f_n[mb->nf], 0, 0);
		mb->f_Qfunc[mb->nf] = Qfunc_clinear;
		mb->f_Qfunc_arg[mb->nf] = (void *) def;
		mb->f_N[mb->nf] = mb->f_n[mb->nf];
		mb->f_rankdef[mb->nf] = 0.0;
		mb->f_bfunc2[mb->nf] = Calloc(1, GMRFLib_bfunc2_tp);
		mb->f_bfunc2[mb->nf]->graph = mb->f_graph[mb->nf];
		mb->f_bfunc2[mb->nf]->Qfunc = mb->f_Qfunc[mb->nf];
		mb->f_bfunc2[mb->nf]->Qfunc_arg = mb->f_Qfunc_arg[mb->nf];
		mb->f_bfunc2[mb->nf]->diagonal = mb->f_diag[mb->nf];
		mb->f_bfunc2[mb->nf]->mfunc = mfunc_clinear;
		mb->f_bfunc2[mb->nf]->mfunc_arg = mb->f_Qfunc_arg[mb->nf];
		mb->f_bfunc2[mb->nf]->n = mb->f_n[mb->nf];
		mb->f_bfunc2[mb->nf]->nreplicate = 1;
		mb->f_bfunc2[mb->nf]->ngroup = 1;
	}
		break;

	case F_SIGM:
	case F_REVSIGM:
	{
		// mb->f_id[mb->nf]
		inla_sigm_tp *def = Calloc(1, inla_sigm_tp);
		def->beta = beta;
		def->log_halflife = log_halflife;
		def->log_shape = log_shape;
		def->precision = mb->f_precision[mb->nf];
		def->x = mb->f_locations[mb->nf];

		GMRFLib_graph_mk_linear(&(mb->f_graph[mb->nf]), mb->f_n[mb->nf], 0, 0);
		mb->f_Qfunc[mb->nf] = Qfunc_sigm;
		mb->f_Qfunc_arg[mb->nf] = (void *) def;
		mb->f_N[mb->nf] = mb->f_n[mb->nf];
		mb->f_rankdef[mb->nf] = 0.0;
		mb->f_bfunc2[mb->nf] = Calloc(1, GMRFLib_bfunc2_tp);
		mb->f_bfunc2[mb->nf]->graph = mb->f_graph[mb->nf];
		mb->f_bfunc2[mb->nf]->Qfunc = mb->f_Qfunc[mb->nf];
		mb->f_bfunc2[mb->nf]->Qfunc_arg = mb->f_Qfunc_arg[mb->nf];
		mb->f_bfunc2[mb->nf]->diagonal = mb->f_diag[mb->nf];
		mb->f_bfunc2[mb->nf]->mfunc = (mb->f_id[mb->nf] == F_SIGM ? mfunc_sigm : mfunc_revsigm);
		mb->f_bfunc2[mb->nf]->mfunc_arg = mb->f_Qfunc_arg[mb->nf];
		mb->f_bfunc2[mb->nf]->n = mb->f_n[mb->nf];
		mb->f_bfunc2[mb->nf]->nreplicate = 1;
		mb->f_bfunc2[mb->nf]->ngroup = 1;
	}
		break;

	case F_LOG1EXP:
	{
		inla_log1exp_tp *def = Calloc(1, inla_log1exp_tp);
		def->beta = beta;
		def->alpha = alpha;
		def->gamma = gama;
		def->precision = mb->f_precision[mb->nf];
		def->x = mb->f_locations[mb->nf];

		GMRFLib_graph_mk_linear(&(mb->f_graph[mb->nf]), mb->f_n[mb->nf], 0, 0);
		mb->f_Qfunc[mb->nf] = Qfunc_log1exp;
		mb->f_Qfunc_arg[mb->nf] = (void *) def;
		mb->f_N[mb->nf] = mb->f_n[mb->nf];
		mb->f_rankdef[mb->nf] = 0.0;
		mb->f_bfunc2[mb->nf] = Calloc(1, GMRFLib_bfunc2_tp);
		mb->f_bfunc2[mb->nf]->graph = mb->f_graph[mb->nf];
		mb->f_bfunc2[mb->nf]->Qfunc = mb->f_Qfunc[mb->nf];
		mb->f_bfunc2[mb->nf]->Qfunc_arg = mb->f_Qfunc_arg[mb->nf];
		mb->f_bfunc2[mb->nf]->diagonal = mb->f_diag[mb->nf];
		mb->f_bfunc2[mb->nf]->mfunc = mfunc_log1exp;
		mb->f_bfunc2[mb->nf]->mfunc_arg = mb->f_Qfunc_arg[mb->nf];
		mb->f_bfunc2[mb->nf]->n = mb->f_n[mb->nf];
		mb->f_bfunc2[mb->nf]->nreplicate = 1;
		mb->f_bfunc2[mb->nf]->ngroup = 1;
	}
		break;

	case F_LOGDIST:
	{
		inla_logdist_tp *def = Calloc(1, inla_logdist_tp);
		def->beta = beta;
		def->alpha1 = alpha1;
		def->alpha2 = alpha2;
		def->precision = mb->f_precision[mb->nf];
		def->x = mb->f_locations[mb->nf];

		GMRFLib_graph_mk_linear(&(mb->f_graph[mb->nf]), mb->f_n[mb->nf], 0, 0);
		mb->f_Qfunc[mb->nf] = Qfunc_logdist;
		mb->f_Qfunc_arg[mb->nf] = (void *) def;
		mb->f_N[mb->nf] = mb->f_n[mb->nf];
		mb->f_rankdef[mb->nf] = 0.0;
		mb->f_bfunc2[mb->nf] = Calloc(1, GMRFLib_bfunc2_tp);
		mb->f_bfunc2[mb->nf]->graph = mb->f_graph[mb->nf];
		mb->f_bfunc2[mb->nf]->Qfunc = mb->f_Qfunc[mb->nf];
		mb->f_bfunc2[mb->nf]->Qfunc_arg = mb->f_Qfunc_arg[mb->nf];
		mb->f_bfunc2[mb->nf]->diagonal = mb->f_diag[mb->nf];
		mb->f_bfunc2[mb->nf]->mfunc = mfunc_logdist;
		mb->f_bfunc2[mb->nf]->mfunc_arg = mb->f_Qfunc_arg[mb->nf];
		mb->f_bfunc2[mb->nf]->n = mb->f_n[mb->nf];
		mb->f_bfunc2[mb->nf]->nreplicate = 1;
		mb->f_bfunc2[mb->nf]->ngroup = 1;
	}
		break;

	case F_BESAG:
	{
		inla_besag_Qfunc_arg_tp *arg = NULL;

		mb->f_Qfunc[mb->nf] = Qfunc_besag;
		arg = Calloc(1, inla_besag_Qfunc_arg_tp);
		GMRFLib_graph_duplicate(&(arg->graph), mb->f_graph[mb->nf]);
		arg->log_prec = log_prec;

		int adj = iniparser_getint(ini, inla_string_join(secname, "ADJUST.FOR.CON.COMP"), 1);
		int std = iniparser_getint(ini, inla_string_join(secname, "SCALE.MODEL"), 0);
		if (mb->verbose) {
			printf("\t\tadjust.for.con.comp[%1d]\n", adj);
			printf("\t\tscale.model[%1d]\n", std);
		}
		if (std) {
			inla_besag_scale(thread_id, arg, adj, mb->verbose);
			if (mb->verbose) {
				printf("\t\tscale.model: prec_scale[%g]\n", arg->prec_scale[0]);
			}
		}

		mb->f_Qfunc_arg[mb->nf] = (void *) arg;
		mb->f_rankdef[mb->nf] = 1.0;
		mb->f_N[mb->nf] = mb->f_n[mb->nf];
		mb->f_id[mb->nf] = F_BESAG;

		// arg->log_prec[0][0] = 0;
		// GMRFLib_printf_Qfunc(stderr, mb->f_graph[mb->nf], mb->f_Qfunc[mb->nf], mb->f_Qfunc_arg[mb->nf]);
	}
		break;

	case F_BESAG2:
	{
		inla_besag2_Qfunc_arg_tp *arg = NULL;

		mb->f_Qfunc[mb->nf] = Qfunc_besag2;
		arg = Calloc(1, inla_besag2_Qfunc_arg_tp);
		arg->besag_arg = Calloc(1, inla_besag_Qfunc_arg_tp);
		arg->besag_arg->graph = mb->f_graph[mb->nf];

		// do this like this, as only the first half of the diag gets this correction
		if (mb->f_diag[mb->nf]) {
			arg->diag = mb->f_diag[mb->nf];
			mb->f_diag[mb->nf] = 0.0;
		} else {
			arg->diag = 0.0;
		}

		int adj = iniparser_getint(ini, inla_string_join(secname, "ADJUST.FOR.CON.COMP"), 1);
		int std = iniparser_getint(ini, inla_string_join(secname, "SCALE.MODEL"), 0);
		if (mb->verbose) {
			printf("\t\tadjust.for.con.comp[%1d]\n", adj);
			printf("\t\tscale.model[%1d]\n", std);
		}
		if (std) {
			inla_besag_scale(thread_id, arg->besag_arg, adj, mb->verbose);
			if (mb->verbose) {
				printf("\t\tscale.model: prec_scale[%g]\n", arg->besag_arg->prec_scale[0]);
			}
		}

		inla_make_besag2_graph(&(mb->f_graph[mb->nf]), arg->besag_arg->graph);
		arg->precision = mb->f_precision[mb->nf];
		arg->besag_arg->log_prec = log_prec;
		arg->log_a = a_intern;

		mb->f_Qfunc_arg[mb->nf] = (void *) arg;
		mb->f_rankdef[mb->nf] = 1.0;
		mb->f_N[mb->nf] = mb->f_n[mb->nf];
		mb->f_id[mb->nf] = F_BESAG2;
	}
		break;

	case F_BYM:
	{
		inla_bym_Qfunc_arg_tp *arg = NULL;
		GMRFLib_graph_tp *g = NULL;

		arg = Calloc(1, inla_bym_Qfunc_arg_tp);
		arg->besag_arg = Calloc(1, inla_besag_Qfunc_arg_tp);

		/*
		 * make the new augmented graph 
		 */
		g = mb->f_graph[mb->nf];
		inla_make_bym_graph(&(mb->f_graph[mb->nf]), g);

		/*
		 * args to the 'besag' model (spatial) 
		 */
		GMRFLib_graph_duplicate(&(arg->besag_arg->graph), g);
		arg->besag_arg->log_prec = log_prec1;

		int adj = iniparser_getint(ini, inla_string_join(secname, "ADJUST.FOR.CON.COMP"), 1);
		int std = iniparser_getint(ini, inla_string_join(secname, "SCALE.MODEL"), 0);
		if (mb->verbose) {
			printf("\t\tadjust.for.con.comp[%1d]\n", adj);
			printf("\t\tscale.model[%1d]\n", std);
		}
		if (std) {
			inla_besag_scale(thread_id, arg->besag_arg, adj, mb->verbose);
			if (mb->verbose) {
				printf("\t\tscale.model: prec_scale[%g]\n", arg->besag_arg->prec_scale[0]);
			}
		}

		/*
		 * remaing ones 
		 */
		arg->n = mb->f_n[mb->nf];
		arg->N = mb->f_N[mb->nf] = 2 * mb->f_n[mb->nf];
		arg->log_prec_iid = log_prec0;

		/*
		 * general 
		 */
		mb->f_Qfunc[mb->nf] = Qfunc_bym;
		mb->f_Qfunc_arg[mb->nf] = (void *) arg;
		mb->f_rankdef[mb->nf] = 0.0;		       /* since constr=T is converted to extraconstr this will be corrected further below
							        * to 1 */
		break;
	}

	case F_RW2DIID:
	{
		int std;
		inla_rw2diid_Qfunc_arg_tp *arg = NULL;

		arg = Calloc(1, inla_rw2diid_Qfunc_arg_tp);
		arg->rw2ddef = Calloc(1, GMRFLib_rw2ddef_tp);
		arg->rw2ddef->nrow = mb->f_nrow[mb->nf];
		arg->rw2ddef->ncol = mb->f_ncol[mb->nf];
		arg->rw2ddef->order = 2;
		arg->rw2ddef->bvalue = bvalue;
		arg->rw2ddef->cyclic = mb->f_cyclic[mb->nf];

		inla_make_rw2diid_graph(&(mb->f_graph[mb->nf]), arg->rw2ddef);
		std = iniparser_getint(ini, inla_string_join(secname, "SCALE.MODEL"), 1);
		assert(std == 1);			       /* this has to be true for this model */
		GMRFLib_rw2d_scale(thread_id, arg->rw2ddef);
		if (mb->verbose) {
			printf("\t\tscale.model[%1d]\n", std);
			printf("\t\tscale.model: prec_scale[%g]\n", arg->rw2ddef->prec_scale[0]);
		}

		/*
		 * remaing ones 
		 */
		arg->n = mb->f_n[mb->nf];
		arg->N = mb->f_N[mb->nf] = 2 * mb->f_n[mb->nf];
		arg->log_prec = log_prec;
		arg->logit_phi = phi_intern;

		/*
		 * general 
		 */
		mb->f_Qfunc[mb->nf] = Qfunc_rw2diid;
		mb->f_Qfunc_arg[mb->nf] = (void *) arg;
		mb->f_rankdef[mb->nf] = 0.0;		       /* since constr=T is converted to extraconstr this will be corrected further below
							        * to 1 */
	}
		break;

	case F_BYM2:
	{
		inla_bym2_Qfunc_arg_tp *arg = NULL;
		GMRFLib_graph_tp *g = NULL;

		arg = Calloc(1, inla_bym2_Qfunc_arg_tp);
		arg->besag_arg = Calloc(1, inla_besag_Qfunc_arg_tp);

		/*
		 * make the new augmented graph 
		 */
		g = mb->f_graph[mb->nf];
		inla_make_bym_graph(&(mb->f_graph[mb->nf]), g);	/* yes, its the same graph */
		GMRFLib_graph_duplicate(&(arg->besag_arg->graph), g);

		int adj = iniparser_getint(ini, inla_string_join(secname, "ADJUST.FOR.CON.COMP"), 1);
		int std = iniparser_getint(ini, inla_string_join(secname, "SCALE.MODEL"), 1);

		if (mb->verbose) {
			printf("\t\tadjust.for.con.comp[%1d]\n", adj);
			printf("\t\tscale.model[%1d]\n", std);
		}

		if (std) {
			inla_besag_scale(thread_id, arg->besag_arg, adj, mb->verbose);
			if (mb->verbose) {
				printf("\t\tscale.model: prec_scale[%g]\n", arg->besag_arg->prec_scale[0]);
			}
		} else {
			fprintf(stderr,
				"\n\n*** Warning ***\tModel[%s] in Section[%s] use scale.model=FALSE which is NOT recommended!!!\n\n",
				model, secname);
			arg->besag_arg->prec_scale = Calloc(arg->besag_arg->graph->n, double);
			for (k = 0; k < arg->besag_arg->graph->n; k++) {
				arg->besag_arg->prec_scale[k] = 1.0;
			}
		}

		/*
		 * remaing ones 
		 */
		arg->n = mb->f_n[mb->nf];
		arg->N = mb->f_N[mb->nf] = 2 * mb->f_n[mb->nf];
		arg->log_prec = log_prec0;
		arg->logit_phi = phi_intern;

		/*
		 * general 
		 */
		mb->f_Qfunc[mb->nf] = Qfunc_bym2;
		mb->f_Qfunc_arg[mb->nf] = (void *) arg;
		mb->f_rankdef[mb->nf] = 0.0;		       /* since constr=T is converted to extraconstr this will be corrected further below
							        * to 1 */
	}
		break;

	case F_BESAGPROPER:
	{
		inla_besag_proper_Qfunc_arg_tp *arg = NULL, *arg_orig = NULL;
		arg = Calloc(1, inla_besag_proper_Qfunc_arg_tp);
		arg_orig = Calloc(1, inla_besag_proper_Qfunc_arg_tp);

		mb->f_Qfunc[mb->nf] = Qfunc_besagproper;
		mb->f_Qfunc_orig[mb->nf] = Qfunc_besagproper;
		GMRFLib_graph_duplicate(&arg->graph, mb->f_graph[mb->nf]);
		GMRFLib_graph_duplicate(&arg_orig->graph, mb->f_graph[mb->nf]);
		GMRFLib_graph_duplicate(&mb->f_graph_orig[mb->nf], mb->f_graph[mb->nf]);
		arg->log_prec = log_prec;
		arg->log_diag = log_diag;
		arg_orig->log_prec = arg_orig->log_diag = NULL;
		mb->f_Qfunc_arg[mb->nf] = (void *) arg;
		mb->f_Qfunc_arg_orig[mb->nf] = (void *) arg_orig;
		mb->f_rankdef[mb->nf] = 0.0;
		mb->f_N[mb->nf] = mb->f_n[mb->nf];
		mb->f_id[mb->nf] = F_BESAGPROPER;
	}
		break;

	case F_BESAGPROPER2:
	{
		inla_besag_proper2_Qfunc_arg_tp *arg = NULL, *arg_orig = NULL;
		arg = Calloc(1, inla_besag_proper2_Qfunc_arg_tp);
		arg_orig = Calloc(1, inla_besag_proper2_Qfunc_arg_tp);

		mb->f_Qfunc[mb->nf] = Qfunc_besagproper2;
		mb->f_Qfunc_orig[mb->nf] = Qfunc_besagproper2;
		GMRFLib_graph_duplicate(&arg->graph, mb->f_graph[mb->nf]);
		GMRFLib_graph_duplicate(&arg_orig->graph, mb->f_graph[mb->nf]);
		GMRFLib_graph_duplicate(&mb->f_graph_orig[mb->nf], mb->f_graph[mb->nf]);
		arg->log_prec = log_prec;
		arg->logit_lambda = phi_intern;
		arg_orig->log_prec = arg_orig->logit_lambda = NULL;
		mb->f_Qfunc_arg[mb->nf] = (void *) arg;
		mb->f_Qfunc_arg_orig[mb->nf] = (void *) arg_orig;
		mb->f_rankdef[mb->nf] = 0.0;
		mb->f_N[mb->nf] = mb->f_n[mb->nf];
		mb->f_id[mb->nf] = F_BESAGPROPER;
	}
		break;

	case F_SPDE:
	{
		mb->f_Qfunc[mb->nf] = spde_model->Qfunc;
		mb->f_Qfunc_arg[mb->nf] = spde_model->Qfunc_arg;
		mb->f_graph[mb->nf] = spde_model->graph;
		mb->f_rankdef[mb->nf] = 0;
		mb->f_n[mb->nf] = mb->f_N[mb->nf] = spde_model->n;
	}
		break;

	case F_SPDE2:
	{
		mb->f_Qfunc[mb->nf] = spde2_model->Qfunc;
		mb->f_Qfunc_arg[mb->nf] = spde2_model->Qfunc_arg;
		mb->f_graph[mb->nf] = spde2_model->graph;
		mb->f_rankdef[mb->nf] = 0;
		mb->f_n[mb->nf] = mb->f_N[mb->nf] = spde2_model->n;
	}
		break;

	case F_SPDE3:
	{
		mb->f_Qfunc[mb->nf] = spde3_model->Qfunc;
		mb->f_Qfunc_arg[mb->nf] = spde3_model->Qfunc_arg;
		mb->f_graph[mb->nf] = spde3_model->graph;
		mb->f_rankdef[mb->nf] = 0;
		mb->f_n[mb->nf] = mb->f_N[mb->nf] = spde3_model->n;
	}
		break;

	case F_AR:
	{
		GMRFLib_graph_tp *g = NULL;
		ar_def_tp *def = Calloc(1, ar_def_tp);

		def->n = mb->f_n[mb->nf];
		def->p = mb->f_order[mb->nf];
		assert((def->n > def->p) && (def->p > 0));
		def->log_prec = log_prec;
		def->pacf_intern = pacf_intern;

		def->hold_pacf_intern = Calloc(GMRFLib_CACHE_LEN(), double *);
		def->hold_Q = Calloc(GMRFLib_CACHE_LEN(), double *);
		def->hold_Qmarg = Calloc(GMRFLib_CACHE_LEN(), double *);
		for (i = 0; i < GMRFLib_CACHE_LEN(); i++) {
			def->hold_pacf_intern[i] = Calloc(def->p, double);
			for (j = 0; j < def->p; j++) {
				def->hold_pacf_intern[i][j] = GMRFLib_uniform();
			}
		}

		GMRFLib_graph_mk_linear(&g, def->n, def->p, 0);
		mb->f_Qfunc[mb->nf] = Qfunc_ar;
		mb->f_Qfunc_arg[mb->nf] = (void *) def;
		mb->f_rankdef[mb->nf] = 0.0;
		mb->f_graph[mb->nf] = g;
		mb->f_N[mb->nf] = mb->f_n[mb->nf];
	}
		break;

	case F_RW2D:
	{
		GMRFLib_rw2ddef_tp *arg = NULL;

		mb->f_Qfunc[mb->nf] = GMRFLib_rw2d;
		arg = Calloc(1, GMRFLib_rw2ddef_tp);
		arg->nrow = mb->f_nrow[mb->nf];
		arg->ncol = mb->f_ncol[mb->nf];
		arg->order = 2;
		arg->cyclic = mb->f_cyclic[mb->nf];
		arg->bvalue = bvalue;
		arg->log_prec_omp = log_prec;

		int std = iniparser_getint(ini, inla_string_join(secname, "SCALE.MODEL"), 0);
		if (mb->verbose) {
			printf("\t\tscale.model[%1d]\n", std);
		}
		if (std) {
			GMRFLib_rw2d_scale(thread_id, (void *) arg);
			if (mb->verbose) {
				printf("\t\tscale.model: prec_scale[%g]\n", arg->prec_scale[0]);
			}
		}

		mb->f_Qfunc_arg[mb->nf] = (void *) arg;
		mb->f_rankdef[mb->nf] = (bvalue == GMRFLib_BVALUE_ZERO ? 0.0 : (arg->cyclic ? 1.0 : 3.0));
		mb->f_N[mb->nf] = mb->f_n[mb->nf];
		mb->f_id[mb->nf] = F_RW2D;
		GMRFLib_make_rw2d_graph(&(mb->f_graph[mb->nf]), arg);
	}
		break;

	case F_Z:
	{
		char *Am = NULL, *Bm = NULL;
		GMRFLib_tabulate_Qfunc_tp *Qfunc_A = NULL, *Qfunc_B = NULL;
		GMRFLib_graph_tp *graph_A = NULL, *graph_B = NULL, *graph_AB = NULL, *gs[2];
		inla_z_arg_tp *arg = NULL, *arg_orig = NULL;
		double **log_prec_orig = NULL;

		Am = iniparser_getstring(ini, inla_string_join(secname, "z.Amatrix"), NULL);
		Bm = iniparser_getstring(ini, inla_string_join(secname, "z.Bmatrix"), NULL);
		GMRFLib_tabulate_Qfunc_from_file(&Qfunc_A, &graph_A, Am, zn + zm, NULL);
		GMRFLib_tabulate_Qfunc_from_file(&Qfunc_B, &graph_B, Bm, zn + zm, NULL);

		gs[0] = graph_A;
		gs[1] = graph_B;
		GMRFLib_graph_union(&graph_AB, gs, 2);

		arg = Calloc(1, inla_z_arg_tp);
		arg->n = zn;
		arg->m = zm;
		arg->log_prec = log_prec;
		arg->graph_A = graph_A;
		arg->Qfunc_A = Qfunc_A;
		arg->graph_B = graph_B;
		arg->Qfunc_B = Qfunc_B;
		arg->graph_AB = graph_AB;

		HYPER_NEW(log_prec_orig, log_prec[0][0]);
		arg_orig = Calloc(1, inla_z_arg_tp);
		Memcpy(arg_orig, arg, sizeof(inla_z_arg_tp));
		arg_orig->log_prec = log_prec_orig;

		mb->f_Qfunc[mb->nf] = Qfunc_z;
		mb->f_Qfunc_orig[mb->nf] = Qfunc_z;
		mb->f_Qfunc_arg[mb->nf] = (void *) arg;
		mb->f_Qfunc_arg_orig[mb->nf] = (void *) arg_orig;
		mb->f_rankdef[mb->nf] = 0;		       /* default value */
		GMRFLib_graph_duplicate(&(mb->f_graph[mb->nf]), graph_AB);
		GMRFLib_graph_duplicate(&(mb->f_graph_orig[mb->nf]), graph_AB);
	}
		break;

	case F_SLM:
	{
		char *A1m = NULL, *A2m = NULL, *Bm = NULL, *Cm = NULL;
		GMRFLib_tabulate_Qfunc_tp *Qfunc_A1 = NULL, *Qfunc_A2 = NULL, *Qfunc_B = NULL, *Qfunc_C = NULL;
		GMRFLib_graph_tp *graph_A1 = NULL, *graph_A2 = NULL, *graph_B = NULL, *graph_C = NULL, *graph_slm = NULL, *gs[4];
		inla_slm_arg_tp *arg = NULL, *arg_orig = NULL;
		double **log_prec_orig = NULL, **logit_rho_orig = NULL;

		A1m = iniparser_getstring(ini, inla_string_join(secname, "slm.A1matrix"), NULL);
		A2m = iniparser_getstring(ini, inla_string_join(secname, "slm.A2matrix"), NULL);
		Bm = iniparser_getstring(ini, inla_string_join(secname, "slm.Bmatrix"), NULL);
		Cm = iniparser_getstring(ini, inla_string_join(secname, "slm.Cmatrix"), NULL);

		GMRFLib_tabulate_Qfunc_from_file(&Qfunc_A1, &graph_A1, A1m, slm_n + slm_m, NULL);
		GMRFLib_tabulate_Qfunc_from_file(&Qfunc_A2, &graph_A2, A2m, slm_n + slm_m, NULL);
		GMRFLib_tabulate_Qfunc_from_file(&Qfunc_B, &graph_B, Bm, slm_n + slm_m, NULL);
		GMRFLib_tabulate_Qfunc_from_file(&Qfunc_C, &graph_C, Cm, slm_n + slm_m, NULL);

		gs[0] = graph_A1;
		gs[1] = graph_A2;
		gs[2] = graph_B;
		gs[3] = graph_C;
		GMRFLib_graph_union(&graph_slm, gs, 4);

		arg = Calloc(1, inla_slm_arg_tp);
		arg->rho_min = slm_rho_min;
		arg->rho_max = slm_rho_max;
		arg->n = slm_n;
		arg->m = slm_m;
		arg->log_prec = log_prec;
		arg->logit_rho = rho_intern;

		arg->graph_A1 = graph_A1;
		arg->graph_A2 = graph_A2;
		arg->graph_B = graph_B;
		arg->graph_C = graph_C;
		arg->graph_slm = graph_slm;

		arg->Qfunc_A1 = Qfunc_A1;
		arg->Qfunc_A2 = Qfunc_A2;
		arg->Qfunc_B = Qfunc_B;
		arg->Qfunc_C = Qfunc_C;

		HYPER_NEW(log_prec_orig, log_prec[0][0]);
		HYPER_NEW(logit_rho_orig, rho_intern[0][0]);
		arg_orig = Calloc(1, inla_slm_arg_tp);
		Memcpy(arg_orig, arg, sizeof(inla_slm_arg_tp));
		arg_orig->log_prec = log_prec_orig;
		arg_orig->logit_rho = logit_rho_orig;

		mb->f_Qfunc[mb->nf] = Qfunc_slm;
		mb->f_Qfunc_orig[mb->nf] = Qfunc_slm;
		mb->f_Qfunc_arg[mb->nf] = (void *) arg;
		mb->f_Qfunc_arg_orig[mb->nf] = (void *) arg_orig;
		mb->f_rankdef[mb->nf] = 0;		       /* default value */
		GMRFLib_graph_duplicate(&(mb->f_graph[mb->nf]), graph_slm);
		GMRFLib_graph_duplicate(&(mb->f_graph_orig[mb->nf]), graph_slm);
	}
		break;

	case F_2DIID:
	{
		inla_2diid_arg_tp *arg = NULL;

		mb->f_N[mb->nf] = 2 * mb->f_n[mb->nf];
		arg = Calloc(1, inla_2diid_arg_tp);
		arg->n = mb->f_n[mb->nf];
		arg->log_prec0 = log_prec0;
		arg->log_prec1 = log_prec1;
		arg->rho_intern = rho_intern;
		mb->f_Qfunc_arg[mb->nf] = (void *) arg;
		if (mb->f_id[mb->nf] == F_2DIID) {
			mb->f_Qfunc[mb->nf] = Qfunc_2diid;
			inla_make_2diid_graph(&(mb->f_graph[mb->nf]), arg);
		} else {
			mb->f_Qfunc[mb->nf] = Qfunc_2diid_wishart;
			inla_make_2diid_wishart_graph(&(mb->f_graph[mb->nf]), arg);
		}
		mb->f_rankdef[mb->nf] = 0;
	}
		break;

	case F_IID1D:
	case F_IID2D:
	case F_IID3D:
	case F_IID4D:
	case F_IID5D:
	{
		inla_iid_wishart_arg_tp *arg = NULL;
		int dim = WISHART_DIM(mb->nf);
		assert(dim > 0);

		assert(mb->f_N[mb->nf] == mb->f_n[mb->nf]);
		arg = Calloc(1, inla_iid_wishart_arg_tp);
		arg->dim = dim;
		arg->n = mb->f_n[mb->nf] / dim;		       /* yes */
		arg->N = mb->f_N[mb->nf];
		mb->f_Qfunc_arg[mb->nf] = (void *) arg;
		mb->f_rankdef[mb->nf] = 0;
		arg->log_prec = theta_iidwishart;
		arg->rho_intern = theta_iidwishart + dim;
		arg->hold = Calloc(GMRFLib_CACHE_LEN(), inla_wishart_hold_tp *);
		mb->f_Qfunc[mb->nf] = Qfunc_iid_wishart;
		mb->f_Qfunc_arg[mb->nf] = (void *) arg;
		inla_make_iid_wishart_graph(&(mb->f_graph[mb->nf]), arg);
	}
		break;

	case F_IIDKD:
	{
		inla_iid_wishartk_arg_tp *arg = NULL;
		int dim = mb->f_order[mb->nf];

		assert(mb->f_N[mb->nf] == mb->f_n[mb->nf]);
		arg = Calloc(1, inla_iid_wishartk_arg_tp);
		arg->dim = dim;
		arg->n = mb->f_n[mb->nf] / dim;		       /* yes */
		arg->N = mb->f_N[mb->nf];
		arg->ntheta = INLA_WISHARTK_NTHETA(dim);
		mb->f_Qfunc_arg[mb->nf] = (void *) arg;
		mb->f_rankdef[mb->nf] = 0;
		arg->theta = theta_iidwishart;
		arg->vec = Calloc(GMRFLib_CACHE_LEN(), double *);
		for (i = 0; i < GMRFLib_CACHE_LEN(); i++) {
			arg->vec[i] = Calloc(arg->ntheta, double);
		}
		arg->hold = Calloc(GMRFLib_CACHE_LEN(), inla_wishartk_hold_tp *);
		mb->f_Qfunc[mb->nf] = Qfunc_iid_wishartk;
		mb->f_Qfunc_arg[mb->nf] = (void *) arg;
		inla_make_iid_wishartk_graph(&(mb->f_graph[mb->nf]), arg);
	}
		break;

	case F_INTSLOPE:
	{
		int subject;
		inla_intslope_arg_tp *arg = Calloc(1, inla_intslope_arg_tp);

		arg->n = mb->f_n[mb->nf];		       /* n */
		arg->N = mb->f_N[mb->nf];		       /* n + 2 * nsubject */
		arg->precision = mb->f_precision[mb->nf];
		arg->theta_gamma = intslope_gamma;
		arg->def = intslope_def;
		arg->nsubject = nsubject;
		arg->nstrata = nstrata;

		arg->warg = Calloc(1, inla_iid_wishart_arg_tp);
		arg->warg->dim = 2;
		arg->warg->n = nsubject;
		arg->warg->N = arg->warg->dim * arg->warg->n;
		arg->warg->log_prec = theta_iidwishart;
		arg->warg->rho_intern = theta_iidwishart + arg->warg->dim;
		arg->warg->hold = Calloc(GMRFLib_CACHE_LEN(), inla_wishart_hold_tp *);

		// For each subject, we need to know all those using it. Since we have all stored in the 'intslope_def' matrix,
		// we only need to store the references to the rows in that one.
		arg->subject_idx = Calloc(nsubject, GMRFLib_idx_tp *);
		for (subject = 0; subject < nsubject; subject++) {
			GMRFLib_idx_create(&(arg->subject_idx[subject]));
		}
		for (int idx = 0; idx < n; idx++) {
			subject = GMRFLib_matrix_get(idx, INTSLOPE_SUBJECT, arg->def);
			GMRFLib_idx_add(&(arg->subject_idx[subject]), idx);
		}

		mb->f_rankdef[mb->nf] = 0;
		mb->f_Qfunc[mb->nf] = Qfunc_intslope;
		mb->f_Qfunc_arg[mb->nf] = (void *) arg;
		mb->f_Qfunc_arg_orig[mb->nf] = (void *) arg;   /* need access to the internals later */
		inla_make_intslope_graph(&(mb->f_graph[mb->nf]), arg);
	}
		break;

	case F_SEASONAL:
	{
		GMRFLib_seasonaldef_tp *sdef = NULL;
		int std = iniparser_getint(ini, inla_string_join(secname, "SCALE.MODEL"), 0);

		mb->f_Qfunc[mb->nf] = GMRFLib_seasonal;
		sdef = Calloc(1, GMRFLib_seasonaldef_tp);
		sdef->n = n;
		sdef->s = s;
		sdef->cyclic = mb->f_cyclic[mb->nf];
		sdef->log_prec_omp = log_prec;
		GMRFLib_make_seasonal_graph(&(mb->f_graph[mb->nf]), sdef);
		mb->f_Qfunc_arg[mb->nf] = (void *) sdef;

		if (std) {
			thread_id = 0;
			assert(omp_get_thread_num() == 0);
			GMRFLib_seasonal_scale(thread_id, sdef);
		}
		if (mb->verbose) {
			printf("\t\tscale.model[%1d]\n", std);
			if (std) {
				printf("\t\tscale.model: prec_scale[%g]\n", sdef->prec_scale[0]);
			}
		}

		/*
		 * for the rank-deficieny, we know the result for CYCLIC=FALSE, but we need to compute it for CYCLIC=TRUE 
		 */
		if (!mb->f_cyclic[mb->nf]) {
			mb->f_rankdef[mb->nf] = s - 1.0;
		} else {
			double *chol = NULL, eps = 1.0e-8, *Q = NULL;
			int *map = NULL, rank;
			GMRFLib_graph_tp *g = NULL;

			g = mb->f_graph[mb->nf];
			Q = Calloc(ISQR(n), double);

			for (i = 0; i < n; i++) {
				Q[i + i * n] = GMRFLib_seasonal(thread_id, i, i, NULL, (void *) sdef);
				for (jj = 0; jj < g->nnbs[i]; jj++) {
					j = g->nbs[i][jj];
					Q[j + i * n] = Q[i + j * n] = GMRFLib_seasonal(thread_id, i, j, NULL, (void *) sdef);
				}
			}
			GMRFLib_comp_chol_semidef(&chol, &map, &rank, Q, n, NULL, eps);
			if (mb->verbose) {
				printf("\t\tcomputed default rank deficiency [%1d]\n", n - rank);
			}
			mb->f_rankdef[mb->nf] = n - rank;
			Free(Q);
			Free(chol);
			Free(map);
		}
		mb->f_N[mb->nf] = mb->f_n[mb->nf];
		mb->f_id[mb->nf] = F_SEASONAL;
	}
		break;

	case F_MEC:
	{
		/*
		 * MEC
		 */
		char *filename_s = NULL;
		filename_s = Strdup(iniparser_getstring(ini, inla_string_join(secname, "SCALE"), NULL));
		if (filename_s) {
			if (mb->verbose) {
				printf("\t\tread scale from file=[%s]\n", filename_s);
			}
			inla_read_data_general(&(mb->f_scale[mb->nf]), NULL, NULL, filename_s, mb->predictor_n, 0, 1, mb->verbose, 1.0);
		} else {
			if (mb->verbose) {
				printf("\t\tno scale\n");
			}
			mb->f_scale[mb->nf] = Calloc(mb->predictor_n, double);
			int ii;
			for (ii = 0; ii < mb->predictor_n; ii++) {
				mb->f_scale[mb->nf][ii] = 1.0;
			}
		}

		inla_mec_tp *def = Calloc(1, inla_mec_tp);

		def->beta = beta;
		def->log_prec_obs = log_prec;
		def->mean_x = mean_x;
		def->log_prec_x = log_prec_x;
		def->x_obs = mb->f_locations[mb->nf];
		// must make a copy... (realloc)
		def->scale = Calloc(mb->predictor_n, double);
		def->map_beta_arg = mb->f_theta_map_arg[mb->nf][0];
		Memcpy(def->scale, mb->f_scale[mb->nf], mb->predictor_n * sizeof(double));

		GMRFLib_graph_mk_linear(&(mb->f_graph[mb->nf]), mb->f_n[mb->nf], 0, 0);
		mb->f_Qfunc[mb->nf] = Qfunc_mec;
		mb->f_Qfunc_arg[mb->nf] = (void *) def;
		mb->f_N[mb->nf] = mb->f_n[mb->nf];
		mb->f_rankdef[mb->nf] = 0.0;

		mb->f_bfunc2[mb->nf] = Calloc(1, GMRFLib_bfunc2_tp);
		mb->f_bfunc2[mb->nf]->graph = mb->f_graph[mb->nf];
		mb->f_bfunc2[mb->nf]->Qfunc = mb->f_Qfunc[mb->nf];
		mb->f_bfunc2[mb->nf]->Qfunc_arg = mb->f_Qfunc_arg[mb->nf];
		mb->f_bfunc2[mb->nf]->diagonal = mb->f_diag[mb->nf];
		mb->f_bfunc2[mb->nf]->mfunc = mfunc_mec;
		mb->f_bfunc2[mb->nf]->mfunc_arg = mb->f_Qfunc_arg[mb->nf];
		mb->f_bfunc2[mb->nf]->n = mb->f_n[mb->nf];
		mb->f_bfunc2[mb->nf]->nreplicate = 1;
		mb->f_bfunc2[mb->nf]->ngroup = 1;
	}
		break;

	case F_MEB:
	{
		/*
		 * MEB
		 */
		char *filename_s = NULL;
		filename_s = Strdup(iniparser_getstring(ini, inla_string_join(secname, "SCALE"), NULL));
		if (filename_s) {
			if (mb->verbose) {
				printf("\t\tread scale from file=[%s]\n", filename_s);
			}
			inla_read_data_general(&(mb->f_scale[mb->nf]), NULL, NULL, filename_s, mb->predictor_n, 0, 1, mb->verbose, 1.0);
		} else {
			if (mb->verbose) {
				printf("\t\tno scale\n");
			}
			mb->f_scale[mb->nf] = Calloc(mb->predictor_n, double);
			int ii;
			for (ii = 0; ii < mb->predictor_n; ii++) {
				mb->f_scale[mb->nf][ii] = 1.0;
			}
		}

		inla_meb_tp *def = Calloc(1, inla_meb_tp);

		def->beta = beta;
		def->log_prec = log_prec;
		def->x = mb->f_locations[mb->nf];
		// must make a copy... (realloc)
		def->scale = Calloc(mb->predictor_n, double);
		Memcpy(def->scale, mb->f_scale[mb->nf], mb->predictor_n * sizeof(double));
		def->map_beta_arg = mb->f_theta_map_arg[mb->nf][0];

		GMRFLib_graph_mk_linear(&(mb->f_graph[mb->nf]), mb->f_n[mb->nf], 0, 0);
		mb->f_Qfunc[mb->nf] = Qfunc_meb;
		mb->f_Qfunc_arg[mb->nf] = (void *) def;
		mb->f_N[mb->nf] = mb->f_n[mb->nf];
		mb->f_rankdef[mb->nf] = 0.0;

		mb->f_bfunc2[mb->nf] = Calloc(1, GMRFLib_bfunc2_tp);
		mb->f_bfunc2[mb->nf]->graph = mb->f_graph[mb->nf];
		mb->f_bfunc2[mb->nf]->Qfunc = mb->f_Qfunc[mb->nf];
		mb->f_bfunc2[mb->nf]->Qfunc_arg = mb->f_Qfunc_arg[mb->nf];
		mb->f_bfunc2[mb->nf]->diagonal = mb->f_diag[mb->nf];
		mb->f_bfunc2[mb->nf]->mfunc = mfunc_meb;
		mb->f_bfunc2[mb->nf]->mfunc_arg = mb->f_Qfunc_arg[mb->nf];
		mb->f_bfunc2[mb->nf]->n = mb->f_n[mb->nf];
		mb->f_bfunc2[mb->nf]->nreplicate = 1;
		mb->f_bfunc2[mb->nf]->ngroup = 1;
	}
		break;

	case F_R_GENERIC:
	{
		/*
		 * R_GENERIC
		 */
		inla_rgeneric_tp *def = Calloc(1, inla_rgeneric_tp), *def_orig = Calloc(1, inla_rgeneric_tp);
		double ***tptr = NULL;

		def->filename = Strdup(rgeneric_filename);
		def->model = Strdup(rgeneric_model);
		def->mu = Calloc(GMRFLib_CACHE_LEN(), double *);
		def->mu_param = Calloc(GMRFLib_CACHE_LEN(), double *);
		def->ntheta = mb->f_ntheta[mb->nf];
		def->param = Calloc(GMRFLib_CACHE_LEN(), double *);
		def->Q = Calloc(GMRFLib_CACHE_LEN(), GMRFLib_tabulate_Qfunc_tp *);
		def->reset_cache = 0;			       /* only do if = 0 */
		def->graph = NULL;
		if (def->ntheta) {
			tptr = Calloc(def->ntheta, double **);
			for (j = 0; j < def->ntheta; j++)
				tptr[j] = mb->f_theta[mb->nf][j];
			def->theta = tptr;
		} else {
			def->theta = NULL;
		}

		def_orig->filename = Strdup(rgeneric_filename);
		def_orig->model = Strdup(rgeneric_model);
		def_orig->mu = Calloc(GMRFLib_CACHE_LEN(), double *);
		def_orig->mu_param = Calloc(GMRFLib_CACHE_LEN(), double *);
		def_orig->ntheta = mb->f_ntheta[mb->nf];
		def_orig->param = Calloc(GMRFLib_CACHE_LEN(), double *);
		def_orig->Q = Calloc(GMRFLib_CACHE_LEN(), GMRFLib_tabulate_Qfunc_tp *);
		if (def_orig->ntheta) {
			tptr = Calloc(def_orig->ntheta, double **);
			for (j = 0; j < def_orig->ntheta; j++)
				tptr[j] = mb->f_theta[mb->nf][j];
			def_orig->theta = tptr;
		} else {
			def_orig->theta = NULL;
		}

		int n_out, zero = 0;
		double *x_out = NULL;
		inla_R_rgeneric(&n_out, &x_out, R_GENERIC_GRAPH, def->model, &zero, NULL);

		int len, *ilist = NULL, *jlist = NULL;
		k = 0;
		assert(n_out >= 2);
		n = (int) x_out[k++];
		len = (int) x_out[k++];

		ilist = Calloc(len, int);
		for (i = 0; i < len; i++) {
			ilist[i] = (int) x_out[k++];
		}

		jlist = Calloc(len, int);
		for (i = 0; i < len; i++) {
			jlist[i] = (int) x_out[k++];
		}
		assert(k == n_out);

		double *Qijlist = Calloc(len, double);
		for (i = 0; i < len; i++) {
			Qijlist[i] = 1.0;
		}

		GMRFLib_tabulate_Qfunc_tp *tab = NULL;
		GMRFLib_graph_tp *graph = NULL, *ggraph = NULL;

		GMRFLib_tabulate_Qfunc_from_list(&tab, &graph, len, ilist, jlist, Qijlist, n, NULL);
		GMRFLib_free_tabulate_Qfunc(tab);
		Free(ilist);
		Free(jlist);
		Free(Qijlist);
		Free(x_out);

		def->graph = graph;
		mb->f_graph[mb->nf] = graph;
		mb->f_Qfunc[mb->nf] = Qfunc_rgeneric;
		mb->f_Qfunc_arg[mb->nf] = (void *) def;

		// save the indices for the graph, as we need them repeatedly
		def->len_list = graph->nnz / 2 + graph->n;
		def->ilist = Calloc(def->len_list, int);
		def->jlist = Calloc(def->len_list, int);
		for (i = 0, k = 0; i < graph->n; i++) {
			def->ilist[k] = i;
			def->jlist[k] = i;
			k++;
			for (jj = 0; jj < graph->lnnbs[i]; jj++) {
				j = graph->lnbs[i][jj];
				def->ilist[k] = i;
				def->jlist[k] = j;
				k++;
			}
		}

		// we need to revert the order of the list. pretty annoying...
		// GMRFLib_qsort2((void *) def->jlist, (size_t) def->len_list, sizeof(int), (void *) def->ilist, sizeof(int), NULL, 0,
		// GMRFLib_icmp);
		my_sort2_ii(def->jlist, def->ilist, def->len_list);

		// now we need to sort within each value of jlist.
		assert(def->jlist[0] == 0);
		for (j = k = 0; j < graph->n; j++) {
			for (jj = k; jj < def->len_list; jj++) {
				if (def->jlist[jj] > j)
					break;
			}
			QSORT_FUN((void *) &def->ilist[k], (size_t) (jj - k), sizeof(int), GMRFLib_icmp);
			k = jj;
		}
		assert(k == def->len_list);

		def_orig->len_list = def->len_list;
		def_orig->ilist = Calloc(def_orig->len_list, int);
		def_orig->jlist = Calloc(def_orig->len_list, int);
		Memcpy(def_orig->ilist, def->ilist, def->len_list * sizeof(int));
		Memcpy(def_orig->jlist, def->jlist, def->len_list * sizeof(int));

		GMRFLib_graph_duplicate(&ggraph, graph);
		def_orig->graph = graph;
		mb->f_graph_orig[mb->nf] = ggraph;
		mb->f_Qfunc_orig[mb->nf] = Qfunc_rgeneric;
		mb->f_Qfunc_arg_orig[mb->nf] = (void *) def_orig;

		mb->f_N[mb->nf] = mb->f_n[mb->nf] = def->n = graph->n;
		mb->f_rankdef[mb->nf] = 0.0;

		mb->f_bfunc2[mb->nf] = Calloc(1, GMRFLib_bfunc2_tp);
		mb->f_bfunc2[mb->nf]->graph = mb->f_graph[mb->nf];
		mb->f_bfunc2[mb->nf]->Qfunc = mb->f_Qfunc[mb->nf];
		mb->f_bfunc2[mb->nf]->Qfunc_arg = mb->f_Qfunc_arg[mb->nf];
		mb->f_bfunc2[mb->nf]->diagonal = mb->f_diag[mb->nf];
		mb->f_bfunc2[mb->nf]->mfunc = mfunc_rgeneric;
		mb->f_bfunc2[mb->nf]->mfunc_arg = mb->f_Qfunc_arg[mb->nf];
		mb->f_bfunc2[mb->nf]->n = mb->f_n[mb->nf];
		mb->f_bfunc2[mb->nf]->nreplicate = 1;
		mb->f_bfunc2[mb->nf]->ngroup = 1;
	}
		break;

	case F_C_GENERIC:
	{
		/*
		 * C_GENERIC
		 */
		inla_cgeneric_tp *def = Calloc(1, inla_cgeneric_tp), *def_orig = Calloc(1, inla_cgeneric_tp);
		double ***tptr = NULL;

		def->shlib = Strdup(cgeneric_shlib);
		def->model = Strdup(cgeneric_model);
		def->model_func = model_func;
		def->data = cgeneric_data;
		def->secname = Strdup(secname);
		def->debug = cgeneric_debug;
		def->mu = Calloc(GMRFLib_CACHE_LEN(), double *);
		def->mu_param = Calloc(GMRFLib_CACHE_LEN(), double *);
		def->ntheta = mb->f_ntheta[mb->nf];
		def->param = Calloc(GMRFLib_CACHE_LEN(), double *);
		def->Q = Calloc(GMRFLib_CACHE_LEN(), GMRFLib_tabulate_Qfunc_tp *);
		def->reset_cache = 0;			       /* only do if = 0 */
		def->graph = NULL;
		if (def->ntheta) {
			tptr = Calloc(def->ntheta, double **);
			for (j = 0; j < def->ntheta; j++) {
				tptr[j] = mb->f_theta[mb->nf][j];
			}
			def->theta = tptr;
		} else {
			def->theta = NULL;
		}

		def_orig->shlib = Strdup(cgeneric_shlib);
		def_orig->model = Strdup(cgeneric_model);
		def_orig->model_func = model_func;
		def_orig->data = cgeneric_data;
		def_orig->secname = Strdup(secname);
		def_orig->debug = cgeneric_debug;
		def_orig->mu = Calloc(GMRFLib_CACHE_LEN(), double *);
		def_orig->mu_param = Calloc(GMRFLib_CACHE_LEN(), double *);
		def_orig->ntheta = mb->f_ntheta[mb->nf];
		def_orig->param = Calloc(GMRFLib_CACHE_LEN(), double *);
		def_orig->Q = Calloc(GMRFLib_CACHE_LEN(), GMRFLib_tabulate_Qfunc_tp *);
		if (def_orig->ntheta) {
			tptr = Calloc(def_orig->ntheta, double **);
			for (j = 0; j < def_orig->ntheta; j++) {
				tptr[j] = mb->f_theta[mb->nf][j];
			}
			def_orig->theta = tptr;
		} else {
			def_orig->theta = NULL;
		}


		double *x_out = NULL;
		x_out = model_func(INLA_CGENERIC_GRAPH, NULL, cgeneric_data);
		if (cgeneric_debug) {
			inla_cgeneric_debug(stdout, secname, INLA_CGENERIC_GRAPH, x_out);
		}

		int len, *ilist = NULL, *jlist = NULL;
		k = 0;
		n = (int) x_out[k++];
		len = (int) x_out[k++];
		ilist = Calloc(len, int);
		for (i = 0; i < len; i++) {
			ilist[i] = (int) x_out[k++];
		}
		jlist = Calloc(len, int);
		for (i = 0; i < len; i++) {
			jlist[i] = (int) x_out[k++];
		}

		double *Qijlist = Calloc(len, double);
		for (i = 0; i < len; i++) {
			Qijlist[i] = 1.0;
		}

		GMRFLib_tabulate_Qfunc_tp *tab = NULL;
		GMRFLib_graph_tp *graph = NULL, *ggraph = NULL;

		GMRFLib_tabulate_Qfunc_from_list(&tab, &graph, len, ilist, jlist, Qijlist, n, NULL);
		GMRFLib_free_tabulate_Qfunc(tab);
		Free(ilist);
		Free(jlist);
		Free(Qijlist);
		Free(x_out);

		def->graph = graph;
		mb->f_graph[mb->nf] = graph;
		mb->f_Qfunc[mb->nf] = Qfunc_cgeneric;
		mb->f_Qfunc_arg[mb->nf] = (void *) def;

		// save the indices for the graph, as we need them repeatedly
		def->len_list = graph->nnz / 2 + graph->n;
		def->ilist = Calloc(def->len_list, int);
		def->jlist = Calloc(def->len_list, int);
		for (i = 0, k = 0; i < graph->n; i++) {
			def->ilist[k] = i;
			def->jlist[k] = i;
			k++;
			for (jj = 0; jj < graph->lnnbs[i]; jj++) {
				j = graph->lnbs[i][jj];
				def->ilist[k] = i;
				def->jlist[k] = j;
				assert(def->ilist[k] <= def->jlist[k]);
				k++;
			}
		}

		def_orig->len_list = def->len_list;
		def_orig->ilist = Calloc(def_orig->len_list, int);
		def_orig->jlist = Calloc(def_orig->len_list, int);
		Memcpy(def_orig->ilist, def->ilist, def->len_list * sizeof(int));
		Memcpy(def_orig->jlist, def->jlist, def->len_list * sizeof(int));

		GMRFLib_graph_duplicate(&ggraph, graph);
		def_orig->graph = graph;
		mb->f_graph_orig[mb->nf] = ggraph;
		mb->f_Qfunc_orig[mb->nf] = Qfunc_cgeneric;
		mb->f_Qfunc_arg_orig[mb->nf] = (void *) def_orig;

		mb->f_N[mb->nf] = mb->f_n[mb->nf] = def->n = graph->n;
		mb->f_rankdef[mb->nf] = 0.0;

		mb->f_bfunc2[mb->nf] = Calloc(1, GMRFLib_bfunc2_tp);
		mb->f_bfunc2[mb->nf]->graph = mb->f_graph[mb->nf];
		mb->f_bfunc2[mb->nf]->Qfunc = mb->f_Qfunc[mb->nf];
		mb->f_bfunc2[mb->nf]->Qfunc_arg = mb->f_Qfunc_arg[mb->nf];
		mb->f_bfunc2[mb->nf]->diagonal = mb->f_diag[mb->nf];
		mb->f_bfunc2[mb->nf]->mfunc = mfunc_cgeneric;
		mb->f_bfunc2[mb->nf]->mfunc_arg = mb->f_Qfunc_arg[mb->nf];
		mb->f_bfunc2[mb->nf]->n = mb->f_n[mb->nf];
		mb->f_bfunc2[mb->nf]->nreplicate = 1;
		mb->f_bfunc2[mb->nf]->ngroup = 1;
	}
		break;

	case F_FGN:
	{
		inla_fgn_arg_tp *def = NULL, *def_orig = NULL;

		def = Calloc(1, inla_fgn_arg_tp);
		def->n = mb->f_n[mb->nf];
		assert(def->n > 1);
		def->k = mb->f_order[mb->nf];
		def->N = def->n * (def->k + 1);
		def->log_prec = log_prec;
		def->H_intern = H_intern;
		def->prec_eps = mb->f_precision[mb->nf];

		if (mb->f_locations[mb->nf]) {
			double *t = Calloc((def->k + 1) * mb->f_n[mb->nf], double);
			Memcpy(&t[0], mb->f_locations[mb->nf], mb->f_n[mb->nf] * sizeof(double));
			Memcpy(&t[mb->f_n[mb->nf]], mb->f_locations[mb->nf], mb->f_n[mb->nf] * sizeof(double));

			double start = floor(GMRFLib_max_value(t, mb->f_n[mb->nf], NULL) + 1.0);
			for (int ii = 0; ii < def->k * mb->f_n[mb->nf]; ii++) {
				t[ii + mb->f_n[mb->nf]] = start + ii;
			}
			mb->f_locations[mb->nf] = t;
		}

		double **log_prec_orig = NULL, **H_intern_orig = NULL;
		HYPER_NEW(log_prec_orig, log_prec[0][0]);
		HYPER_NEW(H_intern_orig, H_intern[0][0]);

		def_orig = Calloc(1, inla_fgn_arg_tp);
		def_orig->n = def->n;
		def_orig->k = def->k;
		def_orig->N = def->N;
		def_orig->prec_eps = def->prec_eps;
		def_orig->log_prec = log_prec_orig;
		def_orig->H_intern = H_intern_orig;

		inla_make_fgn_graph(&(mb->f_graph[mb->nf]), def);
		mb->f_Qfunc[mb->nf] = Qfunc_fgn;
		mb->f_Qfunc_arg[mb->nf] = (void *) def;

		inla_make_fgn_graph(&(mb->f_graph_orig[mb->nf]), def_orig);
		mb->f_Qfunc_orig[mb->nf] = Qfunc_fgn;
		mb->f_Qfunc_arg_orig[mb->nf] = (void *) def_orig;

		mb->f_N[mb->nf] = mb->f_n[mb->nf] = mb->f_graph[mb->nf]->n;
		assert(mb->f_N[mb->nf] == def->N);
		mb->f_rankdef[mb->nf] = 0.0;
	}
		break;

	case F_FGN2:
	{
		inla_fgn2_arg_tp *def = NULL, *def_orig = NULL;

		def = Calloc(1, inla_fgn2_arg_tp);
		def->n = mb->f_n[mb->nf];
		assert(def->n > 1);
		def->k = mb->f_order[mb->nf];
		def->N = def->n * def->k;
		def->log_prec = log_prec;
		def->H_intern = H_intern;

		if (mb->f_locations[mb->nf]) {
			double *t = Calloc((def->k + 1) * mb->f_n[mb->nf], double);
			Memcpy(&t[0], mb->f_locations[mb->nf], mb->f_n[mb->nf] * sizeof(double));
			Memcpy(&t[mb->f_n[mb->nf]], mb->f_locations[mb->nf], mb->f_n[mb->nf] * sizeof(double));

			double start = floor(GMRFLib_max_value(t, mb->f_n[mb->nf], NULL) + 1.0);
			for (int ii = 0; ii < def->k * mb->f_n[mb->nf]; ii++) {
				t[ii + mb->f_n[mb->nf]] = start + ii;
			}
			mb->f_locations[mb->nf] = t;
		}

		double **log_prec_orig = NULL, **H_intern_orig = NULL;
		HYPER_NEW(log_prec_orig, log_prec[0][0]);
		HYPER_NEW(H_intern_orig, H_intern[0][0]);

		def_orig = Calloc(1, inla_fgn2_arg_tp);
		def_orig->n = def->n;
		def_orig->k = def->k;
		def_orig->N = def->N;
		def_orig->log_prec = log_prec_orig;
		def_orig->H_intern = H_intern_orig;

		inla_make_fgn2_graph(&(mb->f_graph[mb->nf]), def);
		mb->f_Qfunc[mb->nf] = Qfunc_fgn2;
		mb->f_Qfunc_arg[mb->nf] = (void *) def;

		inla_make_fgn2_graph(&(mb->f_graph_orig[mb->nf]), def_orig);
		mb->f_Qfunc_orig[mb->nf] = Qfunc_fgn2;
		mb->f_Qfunc_arg_orig[mb->nf] = (void *) def_orig;

		mb->f_N[mb->nf] = mb->f_n[mb->nf] = mb->f_graph[mb->nf]->n;
		assert(mb->f_N[mb->nf] == def->N);
		mb->f_rankdef[mb->nf] = 0.0;
	}
		break;

	case F_AR1:
	{
		/*
		 * AR1 
		 */
		inla_ar1_arg_tp *def = NULL;

		def = Calloc(1, inla_ar1_arg_tp);
		def->n = mb->f_n[mb->nf];
		assert(def->n > 1);
		def->cyclic = mb->f_cyclic[mb->nf];
		def->log_prec = log_prec;
		def->phi_intern = phi_intern;
		def->mean = mean_x;
		inla_make_ar1_graph(&(mb->f_graph[mb->nf]), def);
		mb->f_Qfunc[mb->nf] = Qfunc_ar1;
		mb->f_Qfunc_arg[mb->nf] = (void *) def;
		mb->f_N[mb->nf] = mb->f_n[mb->nf];
		mb->f_rankdef[mb->nf] = 0.0;

		mb->f_bfunc2[mb->nf] = Calloc(1, GMRFLib_bfunc2_tp);
		mb->f_bfunc2[mb->nf]->graph = mb->f_graph[mb->nf];
		mb->f_bfunc2[mb->nf]->Qfunc = mb->f_Qfunc[mb->nf];
		mb->f_bfunc2[mb->nf]->Qfunc_arg = mb->f_Qfunc_arg[mb->nf];
		mb->f_bfunc2[mb->nf]->diagonal = mb->f_diag[mb->nf];
		mb->f_bfunc2[mb->nf]->mfunc = mfunc_ar1;
		mb->f_bfunc2[mb->nf]->mfunc_arg = mb->f_Qfunc_arg[mb->nf];
		mb->f_bfunc2[mb->nf]->n = mb->f_n[mb->nf];
		mb->f_bfunc2[mb->nf]->nreplicate = 1;
		mb->f_bfunc2[mb->nf]->ngroup = 1;
	}
		break;

	case F_AR1C:
	{
		/*
		 * AR1C
		 */
		inla_ar1c_arg_tp *def = NULL, *def_orig = NULL;

		def = Calloc(1, inla_ar1c_arg_tp);
		def_orig = Calloc(1, inla_ar1c_arg_tp);

		def->n = iniparser_getint(ini, inla_string_join(secname, "ar1c.n"), -1);
		def->m = iniparser_getint(ini, inla_string_join(secname, "ar1c.m"), -1);
		def->N = def->n + def->m;
		assert(def->N == mb->f_n[mb->nf]);
		if (mb->verbose) {
			printf("\t\tn.ar1=[%1d]\n", def->n);
			printf("\t\tm.beta=[%1d]\n", def->m);
			printf("\t\tN=[%1d]\n", def->N);
		}
		if (def->m > 0) {
			filename = Strdup(iniparser_getstring(ini, inla_string_join(secname, "ar1c.Z"), NULL));
			def->Z = GMRFLib_read_fmesher_file(filename, (long int) 0, -1);
			Free(filename);

			filename = Strdup(iniparser_getstring(ini, inla_string_join(secname, "ar1c.ZZ"), NULL));
			def->ZZ = GMRFLib_read_fmesher_file(filename, (long int) 0, -1);
			Free(filename);

			filename = Strdup(iniparser_getstring(ini, inla_string_join(secname, "ar1c.Qbeta"), NULL));
			def->Qbeta = GMRFLib_read_fmesher_file(filename, (long int) 0, -1);
			Free(filename);
			{
				// compute the log|Qbeta| for the normalizing constant
				gsl_matrix *QQ = gsl_matrix_alloc(def->m, def->m);
				for (i = 0; i < def->m; i++) {
					for (j = 0; j < def->m; j++) {
						gsl_matrix_set(QQ, i, j, GMRFLib_matrix_get(i, j, def->Qbeta));
					}
				}
				def->logdet_Qbeta = GMRFLib_gsl_spd_logdet(QQ);
				gsl_matrix_free(QQ);
				if (mb->verbose) {
					printf("\t\tlog(det(Q.beta))=[%.4g]\n", def->logdet_Qbeta);
				}
			}
		} else {
			filename = iniparser_getstring(ini, inla_string_join(secname, "ar1c.Z"), NULL);	/* in case they are there */
			filename = iniparser_getstring(ini, inla_string_join(secname, "ar1c.ZZ"), NULL);
			filename = iniparser_getstring(ini, inla_string_join(secname, "ar1c.Qbeta"), NULL);

			def->Z = NULL;
			def->ZZ = NULL;
			def->Qbeta = NULL;
			def->logdet_Qbeta = 0.0;
		}

		def->log_prec = log_prec;
		def->phi_intern = phi_intern;

		def_orig->n = def->n;			       /* these are the ones I need in extra () */
		def_orig->m = def->m;
		def_orig->logdet_Qbeta = def->logdet_Qbeta;

		inla_make_ar1c_graph(&(mb->f_graph[mb->nf]), def);
		mb->f_Qfunc[mb->nf] = Qfunc_ar1c;
		mb->f_Qfunc_arg[mb->nf] = (void *) def;
		mb->f_Qfunc_arg_orig[mb->nf] = (void *) def_orig;
		mb->f_N[mb->nf] = mb->f_n[mb->nf];
		mb->f_rankdef[mb->nf] = 0.0;
	}
		break;

	case F_OU:
	{
		/*
		 * OU
		 */
		inla_ou_arg_tp *def = NULL;

		def = Calloc(1, inla_ou_arg_tp);
		def->n = mb->f_n[mb->nf];
		def->log_prec = log_prec;
		def->phi_intern = phi_intern;
		assert(mb->f_locations[mb->nf]);
		def->locations = mb->f_locations[mb->nf];
		inla_make_ou_graph(&(mb->f_graph[mb->nf]), def);
		mb->f_Qfunc[mb->nf] = Qfunc_ou;
		mb->f_Qfunc_arg[mb->nf] = (void *) def;

		/*
		 * need this one later to get 'n'. a copy of the original contents is ok. 
		 */
		mb->f_Qfunc_arg_orig[mb->nf] = (void *) Calloc(1, inla_ou_arg_tp);
		Memcpy(mb->f_Qfunc_arg_orig[mb->nf], mb->f_Qfunc_arg[mb->nf], sizeof(inla_ou_arg_tp));

		mb->f_N[mb->nf] = mb->f_n[mb->nf];
		mb->f_rankdef[mb->nf] = 0.0;
	}
		break;

	case F_MATERN2D:
	{
		/*
		 * MATERN2D
		 */
		GMRFLib_matern2ddef_tp *arg = NULL;

		arg = Calloc(1, GMRFLib_matern2ddef_tp);
		arg->nrow = mb->f_nrow[mb->nf];
		arg->ncol = mb->f_ncol[mb->nf];
		arg->cyclic = mb->f_cyclic[mb->nf];
		arg->nu = mb->f_nu[mb->nf];
		arg->log_prec_omp = log_prec;
		arg->log_range_omp = range_intern;

		GMRFLib_matern2ddef_tp *arg_orig = NULL;
		arg_orig = Calloc(1, GMRFLib_matern2ddef_tp);
		Memcpy(arg_orig, arg, sizeof(GMRFLib_matern2ddef_tp));

		mb->f_Qfunc[mb->nf] = GMRFLib_matern2d;
		mb->f_Qfunc_orig[mb->nf] = GMRFLib_matern2d;
		mb->f_Qfunc_arg[mb->nf] = (void *) arg;
		mb->f_Qfunc_arg_orig[mb->nf] = (void *) arg_orig;
		mb->f_rankdef[mb->nf] = 0.0;
		mb->f_N[mb->nf] = mb->f_n[mb->nf];
		mb->f_id[mb->nf] = F_MATERN2D;
		GMRFLib_make_matern2d_graph(&(mb->f_graph[mb->nf]), arg);
		GMRFLib_make_matern2d_graph(&(mb->f_graph_orig[mb->nf]), arg);
	}
		break;

	case F_DMATERN:
	{
		dmatern_arg_tp *arg = Calloc(1, dmatern_arg_tp);
		dmatern_arg_tp *arg_orig = Calloc(1, dmatern_arg_tp);

		filename = Strdup(iniparser_getstring(ini, inla_string_join(secname, "dmatern.locations"), NULL));
		arg->locations = GMRFLib_read_fmesher_file(filename, (long int) 0, -1);
		Free(filename);

		if (!strcasecmp(mb->f_prior[mb->nf][1].name, "PCRANGE")) {
			// In this case, the parameters of the prior depends on the dimension. This is fixed here, where we
			// compute lambda = -U^(dim/2)*log(alpha) and then the parametes are redefined to be (lambda, dim)
			double U = mb->f_prior[mb->nf][1].parameters[0];
			double alpha_local = mb->f_prior[mb->nf][1].parameters[1];
			double dim = arg->locations->ncol;
			mb->f_prior[mb->nf][1].parameters[0] = -pow(U, dim / 2.0) * log(alpha_local);
			mb->f_prior[mb->nf][1].parameters[1] = dim;
		}

		arg->n = arg->locations->nrow;
		arg->dim = arg->locations->ncol;
		arg->log_range = range_intern;
		arg->log_prec = log_prec;
		arg->log_nu = nu_intern;
		Memcpy(arg_orig, arg, sizeof(dmatern_arg_tp));

		mb->f_Qfunc[mb->nf] = Qfunc_dmatern;
		mb->f_Qfunc_orig[mb->nf] = Qfunc_dmatern;
		mb->f_Qfunc_arg[mb->nf] = (void *) arg;
		mb->f_Qfunc_arg_orig[mb->nf] = (void *) arg_orig;
		// dense graph
		GMRFLib_graph_mk_linear(&(mb->f_graph[mb->nf]), arg->n, arg->n, 0);
		GMRFLib_graph_mk_linear(&(mb->f_graph_orig[mb->nf]), arg->n, arg->n, 0);
		mb->f_rankdef[mb->nf] = 0.0;
		assert(mb->f_n[mb->nf] == arg->n);
		mb->f_N[mb->nf] = mb->f_n[mb->nf];
		mb->f_id[mb->nf] = F_DMATERN;

		// setup cache and prefill parameters with random numbers
		arg->param = Calloc(GMRFLib_CACHE_LEN(), double *);
		arg->Q = Calloc(GMRFLib_CACHE_LEN(), gsl_matrix *);
		arg_orig->param = Calloc(GMRFLib_CACHE_LEN(), double *);
		arg_orig->Q = Calloc(GMRFLib_CACHE_LEN(), gsl_matrix *);

		for (i = 0; i < GMRFLib_CACHE_LEN(); i++) {
			int np = 3;
			arg->param[i] = Calloc(np, double);
			arg_orig->param[i] = Calloc(np, double);
			for (j = 0; j < np; j++) {
				arg->param[i][j] = GMRFLib_uniform();
				arg_orig->param[i][j] = GMRFLib_uniform();
			}
		}

		// compute the distance between locations. 
		arg->dist = gsl_matrix_alloc(arg->n, arg->n);
		arg_orig->dist = arg->dist;		       /* read only */
		for (i = 0; i < arg->n; i++) {
			for (j = i; j < arg->n; j++) {
				if (i == j) {
					gsl_matrix_set(arg->dist, i, j, 0.0);
				} else {
					double dist = 0.0;

					for (k = 0; k < arg->dim; k++) {
						dist += SQR(GMRFLib_matrix_get(i, k, arg->locations)
							    - GMRFLib_matrix_get(j, k, arg->locations));
					}
					dist = sqrt(dist);
					gsl_matrix_set(arg->dist, i, j, dist);
					gsl_matrix_set(arg->dist, j, i, dist);
				}
			}
		}
	}
		break;

	default:
	{
		/*
		 * RW-models. do a special test for cyclic, since this require locations = default
		 */
		if ((mb->f_id[mb->nf] == F_IID || mb->f_id[mb->nf] == F_RW1 || mb->f_id[mb->nf] == F_RW2) && mb->f_cyclic[mb->nf]) {
			GMRFLib_rwdef_tp *rwdef = NULL;

			if (mb->f_locations[mb->nf]) {
				int ok = 1;
				double diff = mb->f_locations[mb->nf][1] - mb->f_locations[mb->nf][0];

				for (j = 2; j < mb->f_n[mb->nf]; j++) {
					if (mb->f_locations[mb->nf][j] - mb->f_locations[mb->nf][j - 1] != diff) {
						ok = 0;
						break;
					}
				}
				if (!ok) {
					fprintf(stderr, "\n*** Warning ***\tModel[%s] in Section[%s] has cyclic = TRUE but values != NULL.\n",
						model, secname);
					fprintf(stderr, "*** Warning ***\tCylic = TRUE is not implemented for non-equal spaced values.\n");
					fprintf(stderr, "*** Warning ***\tAssume values are equal spaced.\n\n");
				}
			}

			int std = iniparser_getint(ini, inla_string_join(secname, "SCALE.MODEL"), 0);

			rwdef = Calloc(1, GMRFLib_rwdef_tp);
			rwdef->n = mb->f_n[mb->nf];
			if (mb->f_id[mb->nf] == F_IID) {
				rwdef->order = 0;

				/*
				 * this case has an extra option: scale
				 */
				char *filename_s = NULL;
				filename_s = Strdup(iniparser_getstring(ini, inla_string_join(secname, "SCALE"), NULL));
				if (filename_s) {
					if (mb->verbose) {
						printf("\t\tread scale from file=[%s]\n", filename_s);
					}
					inla_read_data_general(&(rwdef->scale0), NULL, NULL, filename_s, rwdef->n, 0, 1, mb->verbose, 1.0);
					mb->f_scale[mb->nf] = rwdef->scale0;	/* need a copy */
				} else {
					rwdef->scale0 = NULL;
				}
			} else if (mb->f_id[mb->nf] == F_RW1) {
				rwdef->order = 1;
			} else {
				assert(mb->f_id[mb->nf] == F_RW2);
				rwdef->order = 2;
			}
			assert(rwdef->n > rwdef->order);
			rwdef->log_prec_omp = log_prec;
			rwdef->cyclic = mb->f_cyclic[mb->nf];
			if (mb->f_cyclic[mb->nf]) {
				if (rwdef->order == 0) {
					mb->f_rankdef[mb->nf] = 0.0;
				} else {
					mb->f_rankdef[mb->nf] = 1.0;
				}
			} else {
				abort();
			}
			GMRFLib_make_rw_graph(&(mb->f_graph[mb->nf]), rwdef);
			if (mb->verbose) {
				printf("\t\tscale.model[%1d]\n", std);
			}
			if (std) {
				GMRFLib_rw_scale(thread_id, (void *) rwdef);
				if (mb->verbose) {
					printf("\t\tscale.model: prec_scale[%g]\n", rwdef->prec_scale[0]);
				}
			}

			mb->f_Qfunc[mb->nf] = (rwdef->order == 0 ? GMRFLib_rw0 : GMRFLib_rw);
			mb->f_Qfunc_arg[mb->nf] = (void *) rwdef;
			mb->f_N[mb->nf] = mb->f_graph[mb->nf]->n;
		} else if ((mb->f_id[mb->nf] == F_IID || mb->f_id[mb->nf] == F_RW1 ||
			    mb->f_id[mb->nf] == F_RW2 || mb->f_id[mb->nf] == F_CRW2) && !mb->f_cyclic[mb->nf]) {
			crwdef = Calloc(1, GMRFLib_crwdef_tp);
			crwdef->n = mb->f_n[mb->nf];
			crwdef->log_prec_omp = log_prec;
			if (mb->f_id[mb->nf] == F_IID) {
				crwdef->order = 0;
				/*
				 * this case has an extra option: scale
				 */
				char *filename_s = NULL;
				filename_s = Strdup(iniparser_getstring(ini, inla_string_join(secname, "SCALE"), NULL));
				if (filename_s) {
					if (mb->verbose) {
						printf("\t\tread scale from file=[%s]\n", filename_s);
					}
					inla_read_data_general(&(crwdef->scale0), NULL, NULL, filename_s, crwdef->n, 0, 1, mb->verbose, 1.0);
					mb->f_scale[mb->nf] = crwdef->scale0;	/* need a copy */
				} else {
					crwdef->scale0 = NULL;
				}
				crwdef->layout = GMRFLib_CRW_LAYOUT_SIMPLE;
				mb->f_rankdef[mb->nf] = 0.0;
			} else if (mb->f_id[mb->nf] == F_RW1) {
				crwdef->order = 1;
				crwdef->layout = GMRFLib_CRW_LAYOUT_SIMPLE;
				mb->f_rankdef[mb->nf] = 1.0;
			} else if (mb->f_id[mb->nf] == F_RW2) {
				crwdef->order = 2;
				crwdef->layout = GMRFLib_CRW_LAYOUT_SIMPLE;
				mb->f_rankdef[mb->nf] = 2.0;
			} else if (mb->f_id[mb->nf] == F_CRW2) {
				crwdef->order = 2;
				crwdef->layout = GMRFLib_CRW_LAYOUT_BLOCK;
				mb->f_rankdef[mb->nf] = 2.0;

				/*
				 * duplicate the locations and swap the sign, if they are present
				 */
				if (mb->f_locations[mb->nf]) {
					double *t = Calloc(2 * mb->f_n[mb->nf], double);
					Memcpy(&t[0], mb->f_locations[mb->nf], mb->f_n[mb->nf] * sizeof(double));
					Memcpy(&t[mb->f_n[mb->nf]], mb->f_locations[mb->nf], mb->f_n[mb->nf] * sizeof(double));

					int ii;
					for (ii = mb->f_n[mb->nf]; ii < 2 * mb->f_n[mb->nf]; ii++) {
						t[ii] *= -1.0;
					}

					mb->f_locations[mb->nf] = t;
				}
			} else {
				inla_error_field_is_void(__GMRFLib_FuncName, secname, "model", model);
			}
			crwdef->position = mb->f_locations[mb->nf];	/* do this here, as the locations are duplicated for CRW2 */
			assert(crwdef->n > crwdef->order);
			int std = iniparser_getint(ini, inla_string_join(secname, "SCALE.MODEL"), 0);
			if (mb->f_id[mb->nf] == F_RW1 || mb->f_id[mb->nf] == F_RW2 || mb->f_id[mb->nf] == F_CRW2) {
				if (std) {
					GMRFLib_crw_scale(thread_id, (void *) crwdef);
				}
				if (mb->verbose) {
					printf("\t\tscale.model[%1d]\n", std);
					if (std)
						printf("\t\tscale.model: prec_scale[%g]\n", crwdef->prec_scale[0]);
				}
			} else {
				if (std) {
					GMRFLib_sprintf(&msg,
							"model[%s]. scale.model=TRUE but this model cannot be scaled. Contact developers\n", model);
					inla_error_general(msg);
					exit(1);
				}
			}

			GMRFLib_make_crw_graph(&(mb->f_graph[mb->nf]), crwdef);
			mb->f_Qfunc[mb->nf] = GMRFLib_crw;
			mb->f_Qfunc_arg[mb->nf] = (void *) crwdef;
			mb->f_N[mb->nf] = mb->f_graph[mb->nf]->n;
		} else {
			GMRFLib_ASSERT(0 == 1, GMRFLib_ESNH);
		}
	}
	}

	/*
	 * read optional extra constraint 
	 */
	filename = Strdup(iniparser_getstring(ini, inla_string_join(secname, "EXTRACONSTRAINT"), NULL));
	if (filename) {
		if (mb->verbose) {
			printf("\t\tread extra constraint from file=[%s]\n", filename);
		}
		mb->f_constr[mb->nf] = inla_read_constraint(filename, mb->f_N[mb->nf]);
		if (mb->verbose) {
			int nnc = mb->f_constr[mb->nf]->nc;

			for (j = 0; j < nnc; j++) {
				printf("\t\tConstraint[%1d]\n", j);
				k = 0;
				for (i = 0; i < mb->f_N[mb->nf]; i++) {
					double a = mb->f_constr[mb->nf]->a_matrix[i * nnc + j];
					if (!ISZERO(a) || mb->f_N[mb->nf] <= PREVIEW) {
						printf("\t\t\tA[%1d] = %f\n", i, a);
						k++;
					}
					if (k > PREVIEW)
						break;
				}
				printf("\t\t\te[%1d] = %f\n", j, mb->f_constr[mb->nf]->e_vector[j]);
			}
		}
	}

	/*
	 * hold a copy of the original constraints before `group' and `replicate' 
	 */
	mb->f_constr_orig[mb->nf] = inla_make_constraint(mb->f_N[mb->nf], mb->f_sumzero[mb->nf], mb->f_constr[mb->nf]);

	/*
	 * determine the final rankdef. 
	 */
	rd = iniparser_getdouble(ini, inla_string_join(secname, "RANKDEF"), -1.0);
	if (rd >= 0) {
		/*
		 * if RANKDEF is given, then this is used, not matter what! 
		 */
		mb->f_rankdef[mb->nf] = rd;
		if (mb->verbose) {
			printf("\t\trank-deficiency is *defined* [%g]\n", rd);
		}
	} else {
		/*
		 * use the previously set default value for the rankdef.  only in the case of a proper model, correct for sumzero
		 * constraint 
		 */
		if (ISZERO(mb->f_rankdef[mb->nf])) {
			mb->f_rankdef[mb->nf] = (mb->f_sumzero[mb->nf] ? 1.0 : 0.0);
		}
		/*
		 * if extra constraint(s), then correct for this. OOPS: this *can* be wrong, if the extra constraint are in the
		 * NULL-space of Q, but then the RANKDEF *is* required. 
		 */
		mb->f_rankdef[mb->nf] += (mb->f_constr[mb->nf] ? mb->f_constr[mb->nf]->nc : 0.0);
		if (mb->verbose) {
			printf("\t\tcomputed/guessed rank-deficiency = [%g]\n", mb->f_rankdef[mb->nf]);
		}
	}
	inla_parse_output(mb, ini, sec, &(mb->f_output[mb->nf]));

	/*
	 * for all models except the F_COPY/F_SCOPY one, do the group and replicate expansions 
	 */
	if (mb->f_id[mb->nf] != F_COPY && mb->f_id[mb->nf] != F_SCOPY) {

		if (mb->f_ngroup[mb->nf] > 1) {
			/*
			 * add groups! 
			 */
			ptmp = Strdup("EXCHANGEABLE");
			ptmp = Strdup(iniparser_getstring(ini, inla_string_join(secname, "GROUP.MODEL"), ptmp));
			if (!strcasecmp(ptmp, "EXCHANGEABLE")) {
				mb->f_group_model[mb->nf] = G_EXCHANGEABLE;
			} else if (!strcasecmp(ptmp, "EXCHANGEABLEPOS")) {
				mb->f_group_model[mb->nf] = G_EXCHANGEABLE_POS;
			} else if (!strcasecmp(ptmp, "AR1")) {
				mb->f_group_model[mb->nf] = G_AR1;
			} else if (!strcasecmp(ptmp, "AR")) {
				mb->f_group_model[mb->nf] = G_AR;
			} else if (!strcasecmp(ptmp, "RW1")) {
				mb->f_group_model[mb->nf] = G_RW1;
			} else if (!strcasecmp(ptmp, "RW2")) {
				mb->f_group_model[mb->nf] = G_RW2;
			} else if (!strcasecmp(ptmp, "BESAG")) {
				mb->f_group_model[mb->nf] = G_BESAG;
			} else if (!strcasecmp(ptmp, "IID")) {
				mb->f_group_model[mb->nf] = G_IID;
			} else {
				GMRFLib_sprintf(&msg, "%s: Unknown GROUP.MODEL: %s\n", secname, ptmp);
				inla_error_general(msg);
				abort();
			}

			ptmp2 = Strdup(iniparser_getstring(ini, inla_string_join(secname, "GROUP.GRAPH"), NULL));
			if (ptmp2) {
				GMRFLib_graph_read(&(mb->f_group_graph[mb->nf]), ptmp2);
			}

			mb->f_group_cyclic[mb->nf] = iniparser_getint(ini, inla_string_join(secname, "GROUP.CYCLIC"), 0);
			mb->f_group_order[mb->nf] = iniparser_getint(ini, inla_string_join(secname, "GROUP.ORDER"), -1);
			if (mb->verbose) {
				printf("\t\tgroup.model = %s\n", ptmp);
				printf("\t\tgroup.graph = %s\n", (ptmp2 ? ptmp2 : "<NONE>"));
				printf("\t\tgroup.cyclic = %s\n", (mb->f_group_cyclic[mb->nf] ? "True" : "False"));
				printf("\t\tgroup.order = %1d\n", mb->f_group_order[mb->nf]);
			}

			switch (mb->f_group_model[mb->nf]) {
			case G_EXCHANGEABLE:
			case G_EXCHANGEABLE_POS:
			case G_AR1:
			case G_RW1:
			case G_RW2:
			case G_BESAG:
			case G_IID:
				fixed = iniparser_getboolean(ini, inla_string_join(secname, "GROUP.FIXED"), 0);
				tmp = iniparser_getdouble(ini, inla_string_join(secname, "GROUP.INITIAL"), 0.0);
				if (!fixed && mb->mode_use_mode) {
					tmp = mb->theta_file[mb->theta_counter_file++];
					if (mb->mode_fixed) {
						fixed = 1;
					}
				}
				mb->f_initial[mb->nf] = Realloc(mb->f_initial[mb->nf], mb->f_ntheta[mb->nf] + 1, double);
				_SetInitial(mb->f_ntheta[mb->nf], tmp);
				if (mb->f_group_model[mb->nf] == G_AR1 || mb->f_group_model[mb->nf] == G_EXCHANGEABLE ||
				    mb->f_group_model[mb->nf] == G_EXCHANGEABLE_POS) {
					HYPER_INIT(group_rho_intern, tmp);
					if (mb->verbose) {
						printf("\t\tinitialise group_rho_intern[%g]\n", tmp);
						printf("\t\tgroup.fixed=[%1d]\n", fixed);
					}
				} else {
					HYPER_INIT(group_prec_intern, tmp);
					if (mb->verbose) {
						printf("\t\tinitialise group_prec_intern[%g]\n", tmp);
						printf("\t\tgroup.fixed=[%1d]\n", fixed);
					}
				}
				mb->f_theta[mb->nf] = Realloc(mb->f_theta[mb->nf], mb->f_ntheta[mb->nf] + 1, double **);
				mb->f_fixed[mb->nf] = Realloc(mb->f_fixed[mb->nf], mb->f_ntheta[mb->nf] + 1, int);
				mb->f_prior[mb->nf] = Realloc(mb->f_prior[mb->nf], mb->f_ntheta[mb->nf] + 1, Prior_tp);

				if (mb->f_group_model[mb->nf] == G_AR1 || mb->f_group_model[mb->nf] == G_EXCHANGEABLE ||
				    mb->f_group_model[mb->nf] == G_EXCHANGEABLE_POS) {
					mb->f_theta[mb->nf][mb->f_ntheta[mb->nf]] = group_rho_intern;
				} else {
					mb->f_theta[mb->nf][mb->f_ntheta[mb->nf]] = group_prec_intern;
				}
				mb->f_fixed[mb->nf][mb->f_ntheta[mb->nf]] = fixed;

				switch (mb->f_group_model[mb->nf]) {
				case G_EXCHANGEABLE:
				case G_EXCHANGEABLE_POS:
				{
					inla_read_prior_group(mb, ini, sec, &(mb->f_prior[mb->nf][mb->f_ntheta[mb->nf]]), "GAUSSIAN-group", NULL);
					mb->f_ntheta[mb->nf]++;
				}
					break;

				case G_AR1:
				{
					inla_read_prior_group(mb, ini, sec, &(mb->f_prior[mb->nf][mb->f_ntheta[mb->nf]]), "GAUSSIAN-rho", NULL);
					mb->f_ntheta[mb->nf]++;
				}
					break;

				case G_IID:
				case G_RW1:
				case G_RW2:
				case G_BESAG:
				{
					inla_read_prior_group(mb, ini, sec, &(mb->f_prior[mb->nf][mb->f_ntheta[mb->nf]]), "LOGGAMMA", NULL);
					mb->f_ntheta[mb->nf]++;
				}
					break;

				default:
					GMRFLib_ASSERT(0 == 1, GMRFLib_ESNH);
					abort();
				}

				if (!fixed) {
					/*
					 * add this \theta 
					 */
					mb->theta = Realloc(mb->theta, mb->ntheta + 1, double **);
					mb->theta_hyperid = Realloc(mb->theta_hyperid, mb->ntheta + 1, char *);
					mb->theta_hyperid[mb->ntheta] = mb->f_prior[mb->nf][mb->f_ntheta[mb->nf] - 1].hyperid;
					mb->theta_tag = Realloc(mb->theta_tag, mb->ntheta + 1, char *);
					mb->theta_tag_userscale = Realloc(mb->theta_tag_userscale, mb->ntheta + 1, char *);
					mb->theta_dir = Realloc(mb->theta_dir, mb->ntheta + 1, char *);

					if (mb->f_group_model[mb->nf] == G_AR1 || mb->f_group_model[mb->nf] == G_EXCHANGEABLE ||
					    mb->f_group_model[mb->nf] == G_EXCHANGEABLE_POS) {
						GMRFLib_sprintf(&msg, "Group rho_intern for %s", (secname ? secname : mb->f_tag[mb->nf]));
						mb->theta_tag[mb->ntheta] = msg;
						GMRFLib_sprintf(&msg, "GroupRho for %s", (secname ? secname : mb->f_tag[mb->nf]));
						mb->theta_tag_userscale[mb->ntheta] = msg;
					} else {
						GMRFLib_sprintf(&msg, "Group prec_intern for %s", (secname ? secname : mb->f_tag[mb->nf]));
						mb->theta_tag[mb->ntheta] = msg;
						GMRFLib_sprintf(&msg, "GroupPrec for %s", (secname ? secname : mb->f_tag[mb->nf]));
						mb->theta_tag_userscale[mb->ntheta] = msg;
					}

					GMRFLib_sprintf(&msg, "%s-parameter%1d", mb->f_dir[mb->nf], mb->f_ntheta[mb->nf] - 1);
					mb->theta_dir[mb->ntheta] = msg;
					if (mb->f_group_model[mb->nf] == G_AR1 || mb->f_group_model[mb->nf] == G_EXCHANGEABLE ||
					    mb->f_group_model[mb->nf] == G_EXCHANGEABLE_POS) {
						mb->theta[mb->ntheta] = group_rho_intern;
					} else {
						mb->theta[mb->ntheta] = group_prec_intern;
					}

					mb->theta_map = Realloc(mb->theta_map, mb->ntheta + 1, map_func_tp *);
					mb->theta_map_arg = Realloc(mb->theta_map_arg, mb->ntheta + 1, void *);

					int *ngp = NULL;
					switch (mb->f_group_model[mb->nf]) {
					case G_EXCHANGEABLE:
					{
						mb->theta_map[mb->ntheta] = map_group_rho;
						// need to add a pointer that stays fixed, mb->theta_map_arg[mb->nf] does not!
						ngp = Calloc(1, int);
						*ngp = mb->f_ngroup[mb->nf];
						mb->theta_map_arg[mb->ntheta] = (void *) ngp;
					}
						break;

					case G_EXCHANGEABLE_POS:
					{
						mb->theta_map[mb->ntheta] = map_probability;
						// need to add a pointer that stays fixed, mb->theta_map_arg[mb->nf] does not!
						ngp = Calloc(1, int);
						*ngp = mb->f_ngroup[mb->nf];
						mb->theta_map_arg[mb->ntheta] = (void *) ngp;
					}
						break;

					case G_AR1:
					{
						mb->theta_map[mb->ntheta] = map_rho;
						mb->theta_map_arg[mb->ntheta] = NULL;
					}
						break;

					case G_IID:
					case G_RW1:
					case G_RW2:
					case G_BESAG:
					{
						mb->theta_map[mb->ntheta] = map_precision;
						mb->theta_map_arg[mb->ntheta] = NULL;
					}
						break;

					default:
						inla_error_general("this should not happen");
					}

					Prior_tp *pri = &(mb->f_prior[mb->nf][mb->f_ntheta[mb->nf] - 1]);

					mb->theta_from = Realloc(mb->theta_from, mb->ntheta + 1, char *);
					mb->theta_to = Realloc(mb->theta_to, mb->ntheta + 1, char *);
					mb->theta_from[mb->ntheta] = Strdup(pri->from_theta);
					mb->theta_to[mb->ntheta] = Strdup(pri->to_theta);

					mb->ntheta++;
				}
				break;

			case G_AR:
			{
				int ntheta, ntheta_orig;

				ntheta_orig = mb->f_ntheta[mb->nf];
				ntheta = mb->f_group_order[mb->nf] + 1;
				assert(ntheta <= AR_MAXTHETA + 1 && ntheta >= 1);
				assert(11 == AR_MAXTHETA + 1);

				mb->f_prior[mb->nf] = Realloc(mb->f_prior[mb->nf], ntheta_orig + AR_MAXTHETA + 1, Prior_tp);
				inla_read_prior_group0(mb, ini, sec, &(mb->f_prior[mb->nf][ntheta_orig + 0]), "LOGGAMMA", NULL);
				inla_read_prior_group1(mb, ini, sec, &(mb->f_prior[mb->nf][ntheta_orig + 1]), "PCRHO0", NULL);
				inla_read_prior_group2(mb, ini, sec, &(mb->f_prior[mb->nf][ntheta_orig + 2]), "PCRHO0", NULL);
				inla_read_prior_group3(mb, ini, sec, &(mb->f_prior[mb->nf][ntheta_orig + 3]), "PCRHO0", NULL);
				inla_read_prior_group4(mb, ini, sec, &(mb->f_prior[mb->nf][ntheta_orig + 4]), "PCRHO0", NULL);
				inla_read_prior_group5(mb, ini, sec, &(mb->f_prior[mb->nf][ntheta_orig + 5]), "PCRHO0", NULL);
				inla_read_prior_group6(mb, ini, sec, &(mb->f_prior[mb->nf][ntheta_orig + 6]), "PCRHO0", NULL);
				inla_read_prior_group7(mb, ini, sec, &(mb->f_prior[mb->nf][ntheta_orig + 7]), "PCRHO0", NULL);
				inla_read_prior_group8(mb, ini, sec, &(mb->f_prior[mb->nf][ntheta_orig + 8]), "PCRHO0", NULL);
				inla_read_prior_group9(mb, ini, sec, &(mb->f_prior[mb->nf][ntheta_orig + 9]), "PCRHO0", NULL);
				inla_read_prior_group10(mb, ini, sec, &(mb->f_prior[mb->nf][ntheta_orig + 10]), "PCRHO0", NULL);

				mb->f_initial[mb->nf] = Realloc(mb->f_initial[mb->nf], ntheta_orig + AR_MAXTHETA + 1, double);
				if (mb->verbose) {
					printf("\t\tgroup.ntheta = [%1d]\n", ntheta);
				}

				/*
				 * mark all possible as read 
				 */

				// mark all as read
				for (i = 0; i < AR_MAXTHETA + 1; i++) {
					for (j = 0; j < keywords_len; j++) {
						GMRFLib_sprintf(&ctmp, "GROUP.%s%1d", keywords[j], i);
						iniparser_getstring(ini, inla_string_join(secname, ctmp), NULL);
						Free(ctmp);
					}
				}

				mb->f_fixed[mb->nf] = Realloc(mb->f_fixed[mb->nf], ntheta_orig + AR_MAXTHETA + 1, int);
				mb->f_theta[mb->nf] = Realloc(mb->f_theta[mb->nf], ntheta_orig + AR_MAXTHETA + 1, double **);

				HYPER_NEW(log_prec, 0.0);
				mb->f_theta[mb->nf][ntheta_orig] = log_prec;
				pacf_intern = Calloc(AR_MAXTHETA + 1, double **);
				for (i = 0; i < AR_MAXTHETA; i++) {
					HYPER_NEW(pacf_intern[i], 0.0);
					mb->f_theta[mb->nf][ntheta_orig + i + 1] = pacf_intern[i];
				}

				/*
				 * then read those we need 
				 */
				for (i = 0; i < ntheta; i++) {
					double theta_initial = 0;

					GMRFLib_sprintf(&ctmp, "GROUP.FIXED%1d", i);
					mb->f_fixed[mb->nf][ntheta_orig + i] = iniparser_getboolean(ini, inla_string_join(secname, ctmp), 0);

					GMRFLib_sprintf(&ctmp, "GROUP.INITIAL%1d", i);
					theta_initial = iniparser_getdouble(ini, inla_string_join(secname, ctmp), theta_initial);

					if (!mb->f_fixed[mb->nf][ntheta_orig + i] && mb->mode_use_mode) {
						theta_initial = mb->theta_file[mb->theta_counter_file++];
						if (mb->mode_fixed) {
							mb->f_fixed[mb->nf][ntheta_orig + i] = 1;
						}
					}

					if (i == 0) {
						/*
						 * precision 
						 */
						HYPER_INIT(log_prec, theta_initial);
						if (mb->verbose) {
							printf("\t\tinitialise (log_prec) group.theta[%1d]=[%g]\n", i, theta_initial);
							printf("\t\tfixed[%1d]=[%1d]\n", i, mb->f_fixed[mb->nf][ntheta_orig + i]);
						}

						if (!mb->f_fixed[mb->nf][ntheta_orig + i]) {
							/*
							 * add this \theta 
							 */
							mb->theta = Realloc(mb->theta, mb->ntheta + 1, double **);
							mb->theta_hyperid = Realloc(mb->theta_hyperid, mb->ntheta + 1, char *);
							mb->theta_hyperid[mb->ntheta] = mb->f_prior[mb->nf][ntheta_orig + 0].hyperid;

							mb->theta_tag = Realloc(mb->theta_tag, mb->ntheta + 1, char *);
							mb->theta_tag_userscale = Realloc(mb->theta_tag_userscale, mb->ntheta + 1, char *);
							mb->theta_dir = Realloc(mb->theta_dir, mb->ntheta + 1, char *);
							GMRFLib_sprintf(&msg, "Group Log precision for %s",
									(secname ? secname : mb->f_tag[mb->nf]));

							mb->theta_tag[mb->ntheta] = msg;
							GMRFLib_sprintf(&msg, "Group Precision for %s", (secname ? secname : mb->f_tag[mb->nf]));
							mb->theta_tag_userscale[mb->ntheta] = msg;
							GMRFLib_sprintf(&msg, "%s-parameter%1d", mb->f_dir[mb->nf], ntheta_orig + i + 1);
							mb->theta_dir[mb->ntheta] = msg;

							mb->theta_from = Realloc(mb->theta_from, mb->ntheta + 1, char *);
							mb->theta_to = Realloc(mb->theta_to, mb->ntheta + 1, char *);
							mb->theta_from[mb->ntheta] = Strdup(mb->f_prior[mb->nf][ntheta_orig + 0].from_theta);
							mb->theta_to[mb->ntheta] = Strdup(mb->f_prior[mb->nf][ntheta_orig + 0].to_theta);
							mb->theta[mb->ntheta] = log_prec;
							mb->theta_map = Realloc(mb->theta_map, mb->ntheta + 1, map_func_tp *);
							mb->theta_map[mb->ntheta] = map_precision;
							mb->theta_map_arg = Realloc(mb->theta_map_arg, mb->ntheta + 1, void *);
							mb->theta_map_arg[mb->ntheta] = NULL;
							mb->ntheta++;
						}
					} else {
						/*
						 * PACF 
						 */
						HYPER_INIT(pacf_intern[i - 1], theta_initial);
						if (mb->verbose) {
							printf("\t\tinitialise (PACF) theta[%1d]=[%g]\n", i, theta_initial);
							printf("\t\tfixed[%1d]=[%1d]\n", i, mb->f_fixed[mb->nf][ntheta_orig + i]);
						}
						if (!mb->f_fixed[mb->nf][ntheta_orig + i]) {
							/*
							 * add this \theta 
							 */
							mb->theta = Realloc(mb->theta, mb->ntheta + 1, double **);
							mb->theta_hyperid = Realloc(mb->theta_hyperid, mb->ntheta + 1, char *);
							mb->theta_hyperid[mb->ntheta] = mb->f_prior[mb->nf][ntheta_orig + i].hyperid;
							mb->theta_tag = Realloc(mb->theta_tag, mb->ntheta + 1, char *);
							mb->theta_tag_userscale = Realloc(mb->theta_tag_userscale, mb->ntheta + 1, char *);
							mb->theta_dir = Realloc(mb->theta_dir, mb->ntheta + 1, char *);
							GMRFLib_sprintf(&msg, "Group Intern PACF%1d for %s", i,
									(secname ? secname : mb->f_tag[mb->nf]));

							mb->theta_tag[mb->ntheta] = msg;
							GMRFLib_sprintf(&msg, "Group PACF%1d for %s", i, (secname ? secname : mb->f_tag[mb->nf]));
							mb->theta_tag_userscale[mb->ntheta] = msg;
							GMRFLib_sprintf(&msg, "%s-parameter%1d", mb->f_dir[mb->nf], ntheta_orig + i + 1);
							mb->theta_dir[mb->ntheta] = msg;

							mb->theta_from = Realloc(mb->theta_from, mb->ntheta + 1, char *);
							mb->theta_to = Realloc(mb->theta_to, mb->ntheta + 1, char *);
							mb->theta_from[mb->ntheta] = Strdup(mb->f_prior[mb->nf][ntheta_orig + i].from_theta);
							mb->theta_to[mb->ntheta] = Strdup(mb->f_prior[mb->nf][ntheta_orig + i].to_theta);
							mb->theta[mb->ntheta] = pacf_intern[i - 1];
							mb->theta_map = Realloc(mb->theta_map, mb->ntheta + 1, map_func_tp *);
							mb->theta_map[mb->ntheta] = map_rho;
							mb->theta_map_arg = Realloc(mb->theta_map_arg, mb->ntheta + 1, void *);
							mb->theta_map_arg[mb->ntheta] = NULL;
							mb->ntheta++;
						}
					}
				}
			}
				break;

			default:
				GMRFLib_ASSERT(0 == 1, GMRFLib_ESNH);
			}

			/*
			 * make required changes.  oops, the rankdef is for the size-n model, not the size-N one! 
			 */
			int ng = mb->f_ngroup[mb->nf];
			int Norig = mb->f_N[mb->nf];
			GMRFLib_graph_tp *g = NULL;

			inla_make_group_graph(&g, mb->f_graph[mb->nf], ng, mb->f_group_model[mb->nf], mb->f_group_cyclic[mb->nf],
					      mb->f_group_order[mb->nf], mb->f_group_graph[mb->nf]);
			GMRFLib_graph_free(mb->f_graph[mb->nf]);
			mb->f_graph[mb->nf] = g;

			/*
			 * make the constraints 
			 */
			GMRFLib_constr_tp *c = NULL;
			c = inla_make_constraint2(mb->f_N[mb->nf], mb->f_ngroup[mb->nf], mb->f_sumzero[mb->nf], mb->f_constr[mb->nf]);
			if (c) {
				mb->f_sumzero[mb->nf] = 0;
				Free(mb->f_constr[mb->nf]);
				mb->f_constr[mb->nf] = c;
			}

			/*
			 * redefine the N's, also change the rankdef as its defined for `n'. 
			 */
			mb->f_n[mb->nf] *= ng;
			mb->f_N[mb->nf] *= ng;
			mb->f_rankdef[mb->nf] *= ng;

			int adj = iniparser_getint(ini, inla_string_join(secname, "GROUP.ADJUST.FOR.CON.COMP"), 1);
			int std = iniparser_getint(ini, inla_string_join(secname, "GROUP.SCALE.MODEL"), 0);

			/*
			 * setup the new Qfunc++ 
			 */
			inla_group_def_tp *def = Calloc(1, inla_group_def_tp);

			def->N = Norig;
			def->ngroup = ng;
			def->cyclic = mb->f_group_cyclic[mb->nf];
			def->graph = mb->f_group_graph[mb->nf];
			def->type = mb->f_group_model[mb->nf];
			def->Qfunc = mb->f_Qfunc[mb->nf];
			mb->f_Qfunc[mb->nf] = Qfunc_group;
			def->Qfunc_arg = mb->f_Qfunc_arg[mb->nf];
			mb->f_Qfunc_arg[mb->nf] = (void *) def;
			def->group_rho_intern = group_rho_intern;
			def->group_prec_intern = group_prec_intern;
			if (mb->f_group_model[mb->nf] == G_RW1 || mb->f_group_model[mb->nf] == G_RW2) {
				if (def->cyclic) {
					/*
					 * as cyclic is only implemented for GMRFLib_rw()
					 */
					def->rwdef = Calloc(1, GMRFLib_rwdef_tp);
					def->rwdef->n = ng;
					def->rwdef->order = (mb->f_group_model[mb->nf] == G_RW1 ? 1 : 2);
					def->rwdef->cyclic = mb->f_group_cyclic[mb->nf];
					def->rwdef->log_prec_omp = NULL;
					if (std) {
						char *err = NULL;
						GMRFLib_sprintf(&err, "Group: cannot scale.model with option cylic=TRUE. Contact developers.");
						inla_error_general(err);
						exit(1);
					}
				} else {
					/*
					 * otherwise, we use the general function
					 */
					def->crwdef = Calloc(1, GMRFLib_crwdef_tp);
					def->crwdef->n = ng;
					assert(def->crwdef->n > 0);
					def->crwdef->order = (mb->f_group_model[mb->nf] == G_RW1 ? 1 : 2);
					def->crwdef->log_prec_omp = NULL;
					def->crwdef->layout = GMRFLib_CRW_LAYOUT_SIMPLE;
					def->crwdef->position = Calloc(ng, double);
					int kk;
					for (kk = 0; kk < ng; kk++) {
						def->crwdef->position[kk] = (double) kk;
					}
					if (std) {
						GMRFLib_crw_scale(thread_id, (void *) def->crwdef);
					}
					if (mb->verbose) {
						printf("\t\tgroup.scale.model[%1d]\n", std);
						if (std) {
							printf("\t\tgroup.scale.model: prec_scale[%g]\n", def->crwdef->prec_scale[0]);
						}
					}
				}

			} else if (mb->f_group_model[mb->nf] == G_AR) {
				def->ardef = Calloc(1, ar_def_tp);
				def->ardef->n = mb->f_ngroup[mb->nf];
				def->ardef->p = mb->f_group_order[mb->nf];
				def->ardef->log_prec = log_prec;
				def->ardef->pacf_intern = pacf_intern;
				def->ardef->hold_pacf_intern = Calloc(GMRFLib_CACHE_LEN(), double *);
				def->ardef->hold_Q = Calloc(GMRFLib_CACHE_LEN(), double *);
				def->ardef->hold_Qmarg = Calloc(GMRFLib_CACHE_LEN(), double *);
				for (i = 0; i < GMRFLib_CACHE_LEN(); i++) {
					def->ardef->hold_pacf_intern[i] = Calloc(def->ardef->p, double);
					for (j = 0; j < def->ardef->p; j++) {
						def->ardef->hold_pacf_intern[i][j] = GMRFLib_uniform();
					}
				}
			} else if (mb->f_group_model[mb->nf] == G_BESAG) {
				def->besagdef = Calloc(1, inla_besag_Qfunc_arg_tp);
				def->besagdef->graph = mb->f_group_graph[mb->nf];
				if (mb->verbose) {
					printf("\t\tgroup.scale.model[%1d]\n", std);
					printf("\t\tgroup.adjust.for.con.comp[%1d]\n", std);
				}
				if (std) {
					inla_besag_scale(thread_id, (inla_besag_Qfunc_arg_tp *) (def->besagdef), adj, mb->verbose);
				}
			} else {
				def->rwdef = NULL;
				def->crwdef = NULL;
				def->ardef = NULL;
				def->besagdef = NULL;
			}

			if (mb->f_bfunc2[mb->nf]) {
				/*
				 * then revise the contents
				 */
				mb->f_bfunc2[mb->nf]->graph = mb->f_graph[mb->nf];
				mb->f_bfunc2[mb->nf]->Qfunc = mb->f_Qfunc[mb->nf];
				mb->f_bfunc2[mb->nf]->Qfunc_arg = mb->f_Qfunc_arg[mb->nf];
				mb->f_bfunc2[mb->nf]->ngroup = ng;
			}
		}

		/*
		 * Do the replicate stuff; this is nice hack! 
		 */
		int rep = mb->f_nrep[mb->nf];
		if (rep > 1) {
			inla_replicate_tp *rep_arg = Calloc(1, inla_replicate_tp);
			rep_arg->Qfunc = mb->f_Qfunc[mb->nf];
			rep_arg->Qfunc_arg = mb->f_Qfunc_arg[mb->nf];
			rep_arg->n = mb->f_N[mb->nf];
			inla_replicate_graph(&(mb->f_graph[mb->nf]), rep);	/* this also free the old one */
			mb->f_Qfunc[mb->nf] = Qfunc_replicate;
			mb->f_Qfunc_arg[mb->nf] = (void *) rep_arg;

			GMRFLib_constr_tp *c = NULL;
			c = inla_make_constraint2(mb->f_N[mb->nf], mb->f_nrep[mb->nf], mb->f_sumzero[mb->nf], mb->f_constr[mb->nf]);
			if (c) {
				mb->f_sumzero[mb->nf] = 0;
				Free(mb->f_constr[mb->nf]);
				mb->f_constr[mb->nf] = c;
			}
			// GMRFLib_constr_printf(stdout, c, mb->f_graph[mb->nf]);

			if (mb->f_bfunc2[mb->nf]) {
				/*
				 * Then revise the contents
				 */
				mb->f_bfunc2[mb->nf]->graph = mb->f_graph[mb->nf];
				mb->f_bfunc2[mb->nf]->Qfunc = mb->f_Qfunc[mb->nf];
				mb->f_bfunc2[mb->nf]->Qfunc_arg = mb->f_Qfunc_arg[mb->nf];
				mb->f_bfunc2[mb->nf]->nreplicate = rep;
			}
		}
		mb->f_Ntotal[mb->nf] = mb->f_N[mb->nf] * rep;
	} else {
		mb->f_Ntotal[mb->nf] = -1;		       /* yes this is set later */
	}

	mb->nf++;
#undef _SET
#undef _OneOf
#undef _OneOf2
#undef _OneOf3
#undef _SetInitial
	return INLA_OK;
}


int inla_parse_linear(inla_tp *mb, dictionary *ini, int sec)
{
	/*
	 * parse section = LINEAR 
	 */
	int i;
	char *filename = NULL, *secname = NULL, default_tag[100];

	if (mb->verbose) {
		printf("\tinla_parse_linear...\n");
	}
	secname = Strdup(iniparser_getsecname(ini, sec));
	if (mb->verbose) {
		printf("\t\tsection[%s]\n", secname);
	}
	mb->linear_tag = Realloc(mb->linear_tag, mb->nlinear + 1, char *);
	mb->linear_dir = Realloc(mb->linear_dir, mb->nlinear + 1, char *);
	mb->linear_covariate = Realloc(mb->linear_covariate, mb->nlinear + 1, double *);
	mb->linear_precision = Realloc(mb->linear_precision, mb->nlinear + 1, double);
	mb->linear_mean = Realloc(mb->linear_mean, mb->nlinear + 1, double);
	mb->linear_compute = Realloc(mb->linear_compute, mb->nlinear + 1, int);
	mb->linear_output = Realloc(mb->linear_output, mb->nlinear + 1, Output_tp *);
	sprintf(default_tag, "default tag for linear %d", (int) (10000 * GMRFLib_uniform()));
	mb->linear_tag[mb->nlinear] = (secname ? strdup(secname) : strdup(default_tag));
	mb->linear_dir[mb->nlinear] = Strdup(iniparser_getstring(ini, inla_string_join(secname, "DIR"), Strdup(mb->linear_tag[mb->nlinear])));
	if (mb->verbose) {
		printf("\t\tdir=[%s]\n", mb->linear_dir[mb->nlinear]);
	}
	filename = Strdup(iniparser_getstring(ini, inla_string_join(secname, "COVARIATES"), NULL));
	if (!filename) {
		if (mb->verbose) {
			printf("\t\tfile for covariates=[(NULL)]: set all covariates to 1\n");
		}
		mb->linear_covariate[mb->nlinear] = Calloc(mb->predictor_n, double);

		for (i = 0; i < mb->predictor_n; i++) {
			mb->linear_covariate[mb->nlinear][i] = 1.0;
		}
	} else {
		if (mb->verbose) {
			printf("\t\tfile for covariates=[%s]\n", filename);
		}
		inla_read_data_general(&(mb->linear_covariate[mb->nlinear]), NULL, NULL, filename, mb->predictor_n, 0, 1, mb->verbose, -1.0);
	}
	mb->linear_mean[mb->nlinear] = iniparser_getdouble(ini, inla_string_join(secname, "MEAN"), 0.0);
	if (mb->verbose) {
		printf("\t\tprior mean=[%g]\n", mb->linear_mean[mb->nlinear]);
	}
	mb->linear_precision[mb->nlinear] = iniparser_getdouble(ini, inla_string_join(secname, "PRECISION"), DEFAULT_NORMAL_PRIOR_PRECISION);
	if (mb->verbose) {
		printf("\t\tprior precision=[%g]\n", mb->linear_precision[mb->nlinear]);
	}
	mb->linear_compute[mb->nlinear] = iniparser_getboolean(ini, inla_string_join(secname, "COMPUTE"), 1);
	if (G.mode == INLA_MODE_HYPER) {
		if (mb->linear_compute[mb->nlinear]) {
			fprintf(stderr, "*** Warning: HYPER_MODE require linear_compute[%1d] = 0\n", mb->nlinear);
		}
		mb->linear_compute[mb->nlinear] = 0;
	}
	if (mb->verbose) {
		printf("\t\tcompute=[%1d]\n", mb->linear_compute[mb->nlinear]);
	}
	inla_parse_output(mb, ini, sec, &(mb->linear_output[mb->nlinear]));
	mb->nlinear++;
	return INLA_OK;
}


int inla_parse_INLA(inla_tp *mb, dictionary *ini, int sec, int UNUSED(make_dir))
{
	/*
	 * parse section = INLA 
	 */
	char *secname = NULL, *opt = NULL, *msg = NULL, *filename = NULL, *default_int_strategy = NULL, *defname = NULL, *r = NULL, *ctmp = NULL;
	double tmp, tmp_ref;

	if (mb->verbose) {
		printf("\tinla_parse_INLA...\n");
	}
	secname = Strdup(iniparser_getsecname(ini, sec));
	if (mb->verbose) {
		printf("\t\tsection[%s]\n", secname);
	}

	inla_setup_ai_par_default(mb);			       /* most likely already done, but... */
	mb->lc_derived_correlation_matrix = iniparser_getboolean(ini, inla_string_join(secname, "LINCOMB.DERIVED.CORRELATION.MATRIX"), 0);
	if (mb->verbose) {
		printf("\t\t\tlincomb.derived.correlation.matrix = [%s]\n", (mb->lc_derived_correlation_matrix ? "Yes" : "No"));
	}

	opt = Strdup(iniparser_getstring(ini, inla_string_join(secname, "OPTIMISER"), NULL));
	opt = Strdup(iniparser_getstring(ini, inla_string_join(secname, "OPTIMIZER"), opt));
	if (!opt) {
		mb->ai_par->optimiser = GMRFLib_AI_OPTIMISER_DEFAULT;
	} else if (!strcasecmp(opt, "DEFAULT")) {
		mb->ai_par->optimiser = GMRFLib_AI_OPTIMISER_DEFAULT;
	} else if (!strcasecmp(opt, "GSL")) {
		mb->ai_par->optimiser = GMRFLib_AI_OPTIMISER_GSL;
	} else {
		inla_error_field_is_void(__GMRFLib_FuncName, secname, "optimiser", opt);
	}

	/*
	 * if eps. < 0.0 then factory defaults are used. 
	 */
	mb->ai_par->gsl_tol = iniparser_getdouble(ini, inla_string_join(secname, "GSL.TOL"), mb->ai_par->gsl_tol);
	mb->ai_par->gsl_step_size = iniparser_getdouble(ini, inla_string_join(secname, "GSL.STEP.SIZE"), mb->ai_par->gsl_step_size);
	mb->ai_par->gsl_epsg = iniparser_getdouble(ini, inla_string_join(secname, "TOLERANCE.G"), mb->ai_par->gsl_epsg);
	mb->ai_par->gsl_epsf = iniparser_getdouble(ini, inla_string_join(secname, "TOLERANCE.F"), mb->ai_par->gsl_epsf);
	mb->ai_par->gsl_epsx = iniparser_getdouble(ini, inla_string_join(secname, "TOLERANCE.X"), mb->ai_par->gsl_epsx);
	mb->ai_par->optpar_abserr_func = iniparser_getdouble(ini, inla_string_join(secname, "TOLERANCE.F"), mb->ai_par->optpar_abserr_func);
	mb->ai_par->optpar_abserr_step = iniparser_getdouble(ini, inla_string_join(secname, "TOLERANCE.STEP"), mb->ai_par->optpar_abserr_step);
	mb->ai_par->optpar_nr_step_factor =
	    iniparser_getdouble(ini, inla_string_join(secname, "NR.STEP.FACTOR"), mb->ai_par->optpar_nr_step_factor);
	mb->ai_par->restart = iniparser_getint(ini, inla_string_join(secname, "RESTART"), 0);

	if (mb->verbose > 2) {
		ctmp = Strdup("STDOUT");
	} else {
		ctmp = NULL;
	}
	filename = Strdup(iniparser_getstring(ini, inla_string_join(secname, "OPTPAR.FP"), ctmp));

	if (filename) {
		if (!strcasecmp(filename, "STDOUT")) {
			mb->ai_par->optpar_fp = stdout;
		} else if (!strcasecmp(filename, "STDERR")) {
			mb->ai_par->optpar_fp = stderr;
		} else if (!strcasecmp(filename, "NULL")) {
			mb->ai_par->optpar_fp = NULL;
		} else if (!strcasecmp(filename, "/dev/null")) {
			mb->ai_par->optpar_fp = NULL;
		} else {
			static FILE *fp = NULL;

			fp = fopen(filename, "w");
			if (!fp) {
				GMRFLib_sprintf(&msg, "%s: fail to open file[%s]", __GMRFLib_FuncName, filename);
			}
			mb->ai_par->optpar_fp = fp;
		}
	}

	switch (mb->ai_par->int_strategy) {
	case GMRFLib_AI_INT_STRATEGY_AUTO:
	{
		default_int_strategy = Strdup("GMRFLib_AI_INT_STRATEGY_AUTO");
	}
		break;

	case GMRFLib_AI_INT_STRATEGY_GRID:
	{
		default_int_strategy = Strdup("GMRFLib_AI_INT_STRATEGY_GRID");
	}
		break;

	case GMRFLib_AI_INT_STRATEGY_CCD:
	{
		default_int_strategy = Strdup("GMRFLib_AI_INT_STRATEGY_CCD");
	}
		break;

	default:
		GMRFLib_ASSERT(0 == 1, GMRFLib_ESNH);
	}

	opt = Strdup(iniparser_getstring(ini, inla_string_join(secname, "STRATEGY"), Strdup("AUTO")));
	if (!strcasecmp(opt, "AUTO")) {
		mb->ai_par->strategy = (mb->idx_ntot < 5000 ? GMRFLib_AI_STRATEGY_MEANSKEWCORRECTED_GAUSSIAN : GMRFLib_AI_STRATEGY_ADAPTIVE);
	} else if (!strcasecmp(opt, "GMRFLib_AI_STRATEGY_GAUSSIAN") || !strcasecmp(opt, "GAUSSIAN")) {
		mb->ai_par->strategy = GMRFLib_AI_STRATEGY_GAUSSIAN;
	} else if (!strcasecmp(opt, "GMRFLib_AI_STRATEGY_MEANSKEWCORRECTED_GAUSSIAN") ||
		   !strcasecmp(opt, "MEANSKEWCORRECTED_GAUSSIAN") || !strcasecmp(opt, "SLA") || !strcasecmp(opt, "SIMPLIFIED_LAPLACE")
		   || !strcasecmp(opt, "SIMPLIFIED.LAPLACE")) {
		mb->ai_par->strategy = GMRFLib_AI_STRATEGY_MEANSKEWCORRECTED_GAUSSIAN;
	} else if (!strcasecmp(opt, "GMRFLib_AI_STRATEGY_FIT_SCGAUSSIAN") ||
		   !strcasecmp(opt, "FIT_SCGAUSSIAN") ||
		   !strcasecmp(opt, "FIT.SCGAUSSIAN") || !strcasecmp(opt, "SCGAUSSIAN") || !strcasecmp(opt, "LAPLACE")
		   || !strcasecmp(opt, "LA")) {
		mb->ai_par->strategy = GMRFLib_AI_STRATEGY_FIT_SCGAUSSIAN;
	} else if (!strcasecmp(opt, "GMRFLib_AI_STRATEGY_MEANCORRECTED_GAUSSIAN")
		   || !strcasecmp(opt, "MEANCORRECTED_GAUSSIAN")) {
		mb->ai_par->strategy = GMRFLib_AI_STRATEGY_MEANCORRECTED_GAUSSIAN;
	} else if (!strcasecmp(opt, "GMRFLib_AI_STRATEGY_MEANSKEWCORRECTED_GAUSSIAN") ||
		   !strcasecmp(opt, "MEANSKEWCORRECTED_GAUSSIAN") || !strcasecmp(opt, "SLA") || !strcasecmp(opt, "SIMPLIFIED_LAPLACE")
		   || !strcasecmp(opt, "SIMPLIFIED.LAPLACE")) {
		mb->ai_par->strategy = GMRFLib_AI_STRATEGY_MEANSKEWCORRECTED_GAUSSIAN;
	} else if (!strcasecmp(opt, "ADAPTIVE")) {
		mb->ai_par->strategy = GMRFLib_AI_STRATEGY_ADAPTIVE;
	} else {
		inla_error_field_is_void(__GMRFLib_FuncName, secname, "strategy", opt);
	}
	mb->ai_par->adapt_max = iniparser_getint(ini, inla_string_join(secname, "ADAPTIVE.MAX"), 5);

	mb->ai_par->fast = iniparser_getboolean(ini, inla_string_join(secname, "FAST"), mb->ai_par->fast);
	opt = Strdup(iniparser_getstring(ini, inla_string_join(secname, "LINEAR.CORRECTION"), NULL));
	if (opt) {
		if (!strcasecmp(opt, "GMRFLib_AI_LINEAR_CORRECTION_CENTRAL_DIFFERENCE") || !strcasecmp(opt, "CENTRAL_DIFFERENCE")) {
			mb->ai_par->linear_correction = GMRFLib_AI_LINEAR_CORRECTION_CENTRAL_DIFFERENCE;
		} else if (!strcasecmp(opt, "GMRFLib_AI_LINEAR_CORRECTION_FAST") || !strcasecmp(opt, "FAST")
			   || !strcasecmp(opt, "1") || !strcasecmp(opt, "ON") || !strcasecmp(opt, "YES")
			   || !strcasecmp(opt, "TRUE")) {
			mb->ai_par->linear_correction = GMRFLib_AI_LINEAR_CORRECTION_FAST;
		} else if (!strcasecmp(opt, "GMRFLib_AI_LINEAR_CORRECTION_OFF") || !strcasecmp(opt, "OFF") || !strcasecmp(opt, "NO")
			   || !strcasecmp(opt, "0") || !strcasecmp(opt, "FALSE")) {
			mb->ai_par->linear_correction = GMRFLib_AI_LINEAR_CORRECTION_OFF;
		} else {
			inla_error_field_is_void(__GMRFLib_FuncName, secname, "linear_correction", opt);
		}
	}
	mb->ai_par->n_points = iniparser_getint(ini, inla_string_join(secname, "N.POINTS"), mb->ai_par->n_points);
	mb->ai_par->n_points = iniparser_getint(ini, inla_string_join(secname, "NPOINTS"), mb->ai_par->n_points);

	mb->ai_par->stencil = iniparser_getint(ini, inla_string_join(secname, "STENCIL"), mb->ai_par->stencil);
	mb->ai_par->step_len = iniparser_getdouble(ini, inla_string_join(secname, "STEP.LEN"), mb->ai_par->step_len);
	if (ISZERO(mb->ai_par->step_len)) {
		double scale = GSL_DBL_EPSILON / 2.220446049e-16;
		mb->ai_par->step_len = scale * (mb->ai_par->stencil == 5 ? 1.0e-4 : (mb->ai_par->stencil == 7 ? 5.0e-4 : 1.0e-3));
	}

	mb->ai_par->cutoff = iniparser_getdouble(ini, inla_string_join(secname, "CUTOFF"), mb->ai_par->cutoff);
	filename = Strdup(iniparser_getstring(ini, inla_string_join(secname, "FP.LOG"), NULL));
	if (filename) {
		if (!strcasecmp(filename, "STDOUT")) {
			mb->ai_par->fp_log = stdout;
		} else if (!strcasecmp(filename, "STDERR")) {
			mb->ai_par->fp_log = stderr;
		} else if (!strcasecmp(filename, "NULL")) {
			mb->ai_par->fp_log = NULL;
		} else if (!strcasecmp(filename, "/dev/null")) {
			mb->ai_par->fp_log = NULL;
		} else {
			static FILE *fp = NULL;

			fp = fopen(filename, "w");
			if (!fp) {
				GMRFLib_sprintf(&msg, "%s: fail to open file[%s]", __GMRFLib_FuncName, filename);
			}
			mb->ai_par->fp_log = fp;
		}
	}
	GMRFLib_sprintf(&defname, ".inla_hyper");
	filename = Strdup(iniparser_getstring(ini, inla_string_join(secname, "FP.HYPERPARAM"), defname));
	Free(defname);
	if (!filename) {
		filename = Strdup(iniparser_getstring(ini, inla_string_join(secname, "FP.HYPERPARAM"), NULL));
	}
	if (filename) {
		if (!strcasecmp(filename, "STDOUT")) {
			mb->ai_par->fp_hyperparam = stdout;
		} else if (!strcasecmp(filename, "STDERR")) {
			mb->ai_par->fp_hyperparam = stderr;
		} else if (!strcasecmp(filename, "NULL")) {
			mb->ai_par->fp_hyperparam = NULL;
		} else if (!strcasecmp(filename, "/dev/null")) {
			mb->ai_par->fp_hyperparam = NULL;
		} else {
			static FILE *fp = NULL;

			fp = fopen(filename, "w");
			if (!fp) {
				GMRFLib_sprintf(&msg, "%s: fail to open file[%s]", __GMRFLib_FuncName, filename);
			}
			mb->ai_par->fp_hyperparam = fp;
		}
	}
	opt = Strdup(iniparser_getstring(ini, inla_string_join(secname, "INT.STRATEGY"), default_int_strategy));
	if (opt) {
		if (!strcasecmp(opt, "GMRFLib_AI_INT_STRATEGY_AUTO") || !strcasecmp(opt, "AUTO")) {
			mb->ai_par->int_strategy = GMRFLib_AI_INT_STRATEGY_AUTO;
		} else if (!strcasecmp(opt, "GMRFLib_AI_INT_STRATEGY_GRID") || !strcasecmp(opt, "GRID")) {
			mb->ai_par->int_strategy = GMRFLib_AI_INT_STRATEGY_GRID;
		} else if (!strcasecmp(opt, "GMRFLib_AI_INT_STRATEGY_CCD") || !strcasecmp(opt, "CCD")) {
			mb->ai_par->int_strategy = GMRFLib_AI_INT_STRATEGY_CCD;
		} else if (!strcasecmp(opt, "GMRFLib_AI_INT_STRATEGY_USER") || !strcasecmp(opt, "USER")) {
			mb->ai_par->int_strategy = GMRFLib_AI_INT_STRATEGY_USER;
		} else if (!strcasecmp(opt, "GMRFLib_AI_INT_STRATEGY_USER_STD") || !strcasecmp(opt, "USERSTD")
			   || !strcasecmp(opt, "USER.STD")) {
			mb->ai_par->int_strategy = GMRFLib_AI_INT_STRATEGY_USER_STD;
		} else if (!strcasecmp(opt, "GMRFLib_AI_INT_STRATEGY_USER_EXPERT") || !strcasecmp(opt, "USEREXPERT")
			   || !strcasecmp(opt, "USER.EXPERT")) {
			mb->ai_par->int_strategy = GMRFLib_AI_INT_STRATEGY_USER_EXPERT;
		} else if (!strcasecmp(opt, "GMRFLib_AI_INT_STRATEGY_EMPIRICAL_BAYES")
			   || !strcasecmp(opt, "EMPIRICAL_BAYES") || !strcasecmp(opt, "EB")) {
			mb->ai_par->int_strategy = GMRFLib_AI_INT_STRATEGY_EMPIRICAL_BAYES;
		} else {
			inla_error_field_is_void(__GMRFLib_FuncName, secname, "int_strategy", opt);
		}
	}
	if (G.mode == INLA_MODE_HYPER) {
		if (mb->ai_par->int_strategy != GMRFLib_AI_INT_STRATEGY_GRID) {
			fprintf(stderr, "*** Warning: HYPER_MODE require int_strategy = GMRFLib_AI_INT_STRATEGY_GRID\n");
		}
		mb->ai_par->int_strategy = GMRFLib_AI_INT_STRATEGY_GRID;
	}

	if (mb->ai_par->int_strategy == GMRFLib_AI_INT_STRATEGY_USER || mb->ai_par->int_strategy == GMRFLib_AI_INT_STRATEGY_USER_STD ||
	    mb->ai_par->int_strategy == GMRFLib_AI_INT_STRATEGY_USER_EXPERT) {
		GMRFLib_matrix_tp *D = NULL;
		filename = Strdup(iniparser_getstring(ini, inla_string_join(secname, "INT.DESIGN"), NULL));
		if (my_file_exists(filename) != INLA_OK)
			inla_error_field_is_void(__GMRFLib_FuncName, secname, "int.design", filename);
		D = GMRFLib_read_fmesher_file(filename, (long int) 0, -1);
		GMRFLib_design_read(&(mb->ai_par->int_design), D, (mb->ai_par->int_strategy == GMRFLib_AI_INT_STRATEGY_USER_STD ? 1 : 0));
		GMRFLib_matrix_free(D);
	} else {
		// Just mark it as read
		iniparser_getstring(ini, inla_string_join(secname, "INT.DESIGN"), NULL);
	}

	mb->ai_par->f0 = iniparser_getdouble(ini, inla_string_join(secname, "F0"), mb->ai_par->f0);
	tmp = iniparser_getdouble(ini, inla_string_join(secname, "DZ"), mb->ai_par->dz);
	if (G.mode == INLA_MODE_HYPER && tmp > mb->ai_par->dz) {
		/*
		 * cannot set it to a larger value 
		 */
		fprintf(stderr, "*** Warning: HYPER_MODE require dz <= %f\n", mb->ai_par->dz);
	} else {
		mb->ai_par->dz = tmp;
	}
	mb->ai_par->adjust_weights = iniparser_getboolean(ini, inla_string_join(secname, "ADJUST.WEIGHTS"), mb->ai_par->adjust_weights);

	tmp_ref = mb->ai_par->diff_log_dens;
	mb->ai_par->diff_log_dens = iniparser_getdouble(ini, inla_string_join(secname, "DIFF.LOG.DENS"), mb->ai_par->diff_log_dens);
	mb->ai_par->diff_log_dens = iniparser_getdouble(ini, inla_string_join(secname, "DIFF.LOGDENS"), mb->ai_par->diff_log_dens);
	if (G.mode == INLA_MODE_HYPER && mb->ai_par->diff_log_dens < tmp_ref) {
		fprintf(stderr, "*** Warning: HYPER_MODE require diff_log_dens >= %f\n", tmp_ref);
		mb->ai_par->diff_log_dens = tmp_ref;
	}
	mb->ai_par->skip_configurations =
	    iniparser_getboolean(ini, inla_string_join(secname, "SKIP.CONFIGURATIONS"), mb->ai_par->skip_configurations);

	if (G.mode == INLA_MODE_HYPER && mb->ai_par->skip_configurations) {
		fprintf(stderr, "*** Warning: HYPER_MODE require skip_configurations = 0\n");
		mb->ai_par->skip_configurations = 0;
	}

	/*
	 * this is a short version for setting both: grad=H hess=sqrt(H)
	 */
	mb->ai_par->gradient_finite_difference_step_len =
	    iniparser_getdouble(ini, inla_string_join(secname, "H"), mb->ai_par->gradient_finite_difference_step_len);

	/*
	 * if H < 0, use central difference.  FIXME LATER!!! 
	 */
	if (mb->ai_par->gradient_finite_difference_step_len < 0.0) {
		mb->ai_par->gradient_finite_difference_step_len = ABS(mb->ai_par->gradient_finite_difference_step_len);
		mb->ai_par->gradient_forward_finite_difference = GMRFLib_FALSE;
	}

	mb->ai_par->hessian_finite_difference_step_len =
	    sqrt(ABS(iniparser_getdouble(ini, inla_string_join(secname, "H"), SQR(mb->ai_par->hessian_finite_difference_step_len))));

	/*
	 * ...which is overrided by the original names 
	 */
	char *ans = NULL;

	ans = iniparser_getstring(ini, inla_string_join(secname, "NUM.GRADIENT"), Strdup("central"));
	if (!strcasecmp(ans, "central")) {
		mb->ai_par->gradient_forward_finite_difference = GMRFLib_FALSE;
	} else {
		mb->ai_par->gradient_forward_finite_difference = GMRFLib_TRUE;
	}

	ans = iniparser_getstring(ini, inla_string_join(secname, "NUM.HESSIAN"), Strdup("central"));
	if (!strcasecmp(ans, "central")) {
		mb->ai_par->hessian_forward_finite_difference = GMRFLib_FALSE;
	} else {
		mb->ai_par->hessian_forward_finite_difference = GMRFLib_TRUE;
	}

	ans = iniparser_getstring(ini, inla_string_join(secname, "OPTIMISE.STRATEGY"), Strdup("smart"));
	if (!strcasecmp(ans, "smart")) {
		mb->ai_par->optimise_smart = GMRFLib_TRUE;
	} else {
		mb->ai_par->optimise_smart = GMRFLib_FALSE;
	}

	mb->ai_par->optimise_use_directions = iniparser_getboolean(ini, inla_string_join(secname, "USE.DIRECTIONS"),
								   mb->ai_par->optimise_use_directions);
	filename = iniparser_getstring(ini, inla_string_join(secname, "USE.DIRECTIONS.MATRIX"), NULL);
	if (filename) {
		GMRFLib_matrix_tp *mat = GMRFLib_read_fmesher_file(filename, (long int) 0, -1);
		assert(mat->nrow == mat->ncol);

		gsl_matrix *M = gsl_matrix_alloc((size_t) mat->nrow, (size_t) mat->ncol);
		for (int i = 0; i < mat->nrow; i++) {
			for (int j = 0; j < mat->ncol; j++) {
				gsl_matrix_set(M, (size_t) i, (size_t) j, GMRFLib_matrix_get(i, j, mat));
			}
		}
		mb->ai_par->optimise_use_directions_m = M;
		GMRFLib_matrix_free(mat);
	} else {
		mb->ai_par->optimise_use_directions_m = NULL;
	}

	mb->ai_par->gradient_finite_difference_step_len =
	    iniparser_getdouble(ini, inla_string_join(secname, "GRADIENT.FINITE.DIFFERENCE.STEP.LEN"),
				mb->ai_par->gradient_finite_difference_step_len);
	mb->ai_par->hessian_finite_difference_step_len =
	    iniparser_getdouble(ini, inla_string_join(secname, "HESSIAN.FINITE.DIFFERENCE.STEP.LEN"),
				mb->ai_par->hessian_finite_difference_step_len);
	mb->ai_par->hessian_force_diagonal =
	    iniparser_getboolean(ini, inla_string_join(secname, "HESSIAN.FORCE.DIAGONAL"), mb->ai_par->hessian_force_diagonal);

	opt = Strdup(iniparser_getstring(ini, inla_string_join(secname, "INTERPOLATOR"), NULL));
	if (opt) {
		if (!strcasecmp(opt, "GMRFLib_AI_INTERPOLATOR_WEIGHTED_DISTANCE") || !strcasecmp(opt, "WEIGHTED_DISTANCE")
		    || !strcasecmp(opt, "WEIGHTED.DISTANCE")) {
			mb->ai_par->interpolator = GMRFLib_AI_INTERPOLATOR_WEIGHTED_DISTANCE;
		} else if (!strcasecmp(opt, "GMRFLib_AI_INTERPOLATOR_NEAREST") || !strcasecmp(opt, "NEAREST")) {
			mb->ai_par->interpolator = GMRFLib_AI_INTERPOLATOR_NEAREST;
		} else if (!strcasecmp(opt, "GMRFLib_AI_INTERPOLATOR_LINEAR") || !strcasecmp(opt, "LINEAR")) {
			mb->ai_par->interpolator = GMRFLib_AI_INTERPOLATOR_LINEAR;
		} else if (!strcasecmp(opt, "GMRFLib_AI_INTERPOLATOR_QUADRATIC") || !strcasecmp(opt, "QUADRATIC")) {
			mb->ai_par->interpolator = GMRFLib_AI_INTERPOLATOR_QUADRATIC;
		} else if (!strcasecmp(opt, "GMRFLib_AI_INTERPOLATOR_CCD") || !strcasecmp(opt, "CCD")) {
			mb->ai_par->interpolator = GMRFLib_AI_INTERPOLATOR_CCD;
		} else if (!strcasecmp(opt, "GMRFLib_AI_INTERPOLATOR_CCD_INTEGRATE") || !strcasecmp(opt, "CCDINTEGRATE")) {
			mb->ai_par->interpolator = GMRFLib_AI_INTERPOLATOR_CCD_INTEGRATE;
		} else if (!strcasecmp(opt, "GMRFLib_AI_INTERPOLATOR_GRIDSUM") || !strcasecmp(opt, "GRIDSUM")) {
			mb->ai_par->interpolator = GMRFLib_AI_INTERPOLATOR_GRIDSUM;
		} else if (!strcasecmp(opt, "GMRFLib_AI_INTERPOLATOR_AUTO") || !strcasecmp(opt, "AUTO")) {
			mb->ai_par->interpolator = GMRFLib_AI_INTERPOLATOR_AUTO;
		} else if (!strcasecmp(opt, "GMRFLib_AI_INTERPOLATOR_GAUSSIAN") || !strcasecmp(opt, "GAUSSIAN")) {
			mb->ai_par->interpolator = GMRFLib_AI_INTERPOLATOR_GAUSSIAN;
		} else {
			inla_error_field_is_void(__GMRFLib_FuncName, secname, "interpolator", opt);
		}
	}
	if (mb->ai_par->interpolator == GMRFLib_AI_INTERPOLATOR_GRIDSUM) {
		if (!mb->ai_par->hessian_force_diagonal) {
			GMRFLib_sprintf(&msg, "interpolator=GRIDSUM require hessian_force_diagonal=1 (and skip_configurations=0, recommended)");
			inla_error_general(msg);
		}
	}
	mb->ai_par->do_MC_error_check = iniparser_getboolean(ini, inla_string_join(secname, "DO.MC.ERROR.CHECK"), mb->ai_par->do_MC_error_check);
	mb->ai_par->compute_nparam_eff = iniparser_getboolean(ini, inla_string_join(secname, "COMPUTE.NPARAM.EFF"), mb->ai_par->compute_nparam_eff);

	if (G.mode == INLA_MODE_HYPER) {
		if (mb->ai_par->compute_nparam_eff) {
			fprintf(stderr, "*** Warning: HYPER_MODE require compute_nparam_eff = GMRFLib_FALSE\n");
		}
		mb->ai_par->compute_nparam_eff = GMRFLib_FALSE;
	}

	tmp = iniparser_getboolean(ini, inla_string_join(secname, "HUGE"), -1);
	if (tmp != -1) {
		fprintf(stderr, "\n\n*** Warning *** option control.inla(huge=TRUE) is disabled and obsolete.\n");
		fprintf(stderr, "*** Warning *** use control.compute = list(strategy = \"SMALL|MEDIUM|LARGE|HUGE|DEFAULT\") instead.\n\n");
	}

	GMRFLib_global_node.factor = iniparser_getdouble(ini, inla_string_join(secname, "GLOBAL.NODE.FACTOR"), GMRFLib_global_node.factor);
	assert(GMRFLib_global_node.factor >= 0.0);
	if (mb->verbose) {
		printf("\t\tglobal_node.factor = %.3f\n", GMRFLib_global_node.factor);
	}

	GMRFLib_global_node.degree = iniparser_getdouble(ini, inla_string_join(secname, "GLOBAL.NODE.DEGREE"), GMRFLib_global_node.degree);
	assert(GMRFLib_global_node.degree >= 0);
	if (mb->verbose) {
		printf("\t\tglobal_node.degree = %.1d\n", GMRFLib_global_node.degree);
	}

	Memcpy((void *) &(mb->gn), (void *) &GMRFLib_global_node, sizeof(GMRFLib_global_node_tp));

	r = Strdup(iniparser_getstring(ini, inla_string_join(secname, "REORDERING"), NULL));
	if (mb->verbose) {
		printf("\t\treordering = %s\n", (r ? r : "(default)"));
	}

	if (r) {
		int err;

		/*
		 * both these fail if the reordering is void 
		 */
		int itmp;
		err = inla_sread_ints(&itmp, 1, r);
		G.reorder = (GMRFLib_reorder_tp) itmp;
		if (err) {
			itmp = GMRFLib_reorder_id((const char *) r);
			G.reorder = (GMRFLib_reorder_tp) itmp;
		}
		GMRFLib_reorder = G.reorder;		       /* yes! */
	}

	mb->ai_par->cpo_req_diff_logdens =
	    iniparser_getdouble(ini, inla_string_join(secname, "CPO.REQ.DIFF.LOGDENS"), mb->ai_par->cpo_req_diff_logdens);
	mb->ai_par->cpo_req_diff_logdens = iniparser_getdouble(ini, inla_string_join(secname, "CPO.DIFF"), mb->ai_par->cpo_req_diff_logdens);
	mb->ai_par->cpo_req_diff_logdens = DMAX(0.0, mb->ai_par->cpo_req_diff_logdens);

	mb->ai_par->stupid_search_mode = iniparser_getboolean(ini, inla_string_join(secname, "STUPID.SEARCH"), mb->ai_par->stupid_search_mode);
	mb->ai_par->stupid_search_max_iter =
	    iniparser_getint(ini, inla_string_join(secname, "STUPID.SEARCH.MAX.ITER"), mb->ai_par->stupid_search_max_iter);
	mb->ai_par->stupid_search_factor =
	    iniparser_getdouble(ini, inla_string_join(secname, "STUPID.SEARCH.FACTOR"), mb->ai_par->stupid_search_factor);

	mb->expert_diagonal_emergencey = 0.0;
	mb->expert_diagonal_emergencey = iniparser_getdouble(ini, inla_string_join(secname, "DIAGONAL"), mb->expert_diagonal_emergencey);
	mb->expert_diagonal_emergencey = DMAX(0.0, mb->expert_diagonal_emergencey);
	if (mb->expert_diagonal_emergencey && mb->verbose) {
		printf("\tdiagonal (expert emergency) = %g\n", mb->expert_diagonal_emergencey);
	}

	mb->ai_par->numint_max_fn_eval = iniparser_getint(ini, inla_string_join(secname, "NUMINT.MAXFEVAL"), mb->ai_par->numint_max_fn_eval);
	mb->ai_par->numint_rel_err = iniparser_getdouble(ini, inla_string_join(secname, "NUMINT.RELERR"), mb->ai_par->numint_rel_err);
	mb->ai_par->numint_abs_err = iniparser_getdouble(ini, inla_string_join(secname, "NUMINT.ABSERR"), mb->ai_par->numint_abs_err);

	mb->ai_par->cmin = iniparser_getdouble(ini, inla_string_join(secname, "CMIN"), mb->ai_par->cmin);
	mb->ai_par->b_strategy = iniparser_getint(ini, inla_string_join(secname, "B.STRATEGY"), mb->ai_par->b_strategy);

	mb->ai_par->vb_enable = iniparser_getboolean(ini, inla_string_join(secname, "CONTROL.VB.ENABLE"), 0);
	mb->ai_par->vb_verbose = iniparser_getboolean(ini, inla_string_join(secname, "CONTROL.VB.VERBOSE"), 0);
	mb->ai_par->vb_nodes_mean = (mb->ai_par->vb_enable ? Calloc(1, char) : NULL);
	mb->ai_par->vb_nodes_variance = (mb->ai_par->vb_enable ? Calloc(1, char) : NULL);
	mb->ai_par->vb_iter_max = iniparser_getint(ini, inla_string_join(secname, "CONTROL.VB.ITER.MAX"), 5);
	mb->ai_par->vb_emergency = iniparser_getdouble(ini, inla_string_join(secname, "CONTROL.VB.EMERGENCY"), 10.0);
	mb->ai_par->vb_iter_max = IMAX(1, mb->ai_par->vb_iter_max);
	mb->ai_par->vb_f_enable_limit_mean = iniparser_getint(ini, inla_string_join(secname, "CONTROL.VB.F.ENABLE.LIMIT.MEAN"), 20);
	mb->ai_par->vb_f_enable_limit_variance = iniparser_getint(ini, inla_string_join(secname, "CONTROL.VB.F.ENABLE.LIMIT.VARIANCE"), 5);
	mb->ai_par->vb_f_enable_limit_mean_max = iniparser_getint(ini, inla_string_join(secname, "CONTROL.VB.F.ENABLE.LIMIT.MEAN.MAX"), 1024);
	mb->ai_par->vb_f_enable_limit_variance_max =
	    iniparser_getint(ini, inla_string_join(secname, "CONTROL.VB.F.ENABLE.LIMIT.VARIANCE.MAX"), 768);
	mb->ai_par->vb_hessian_update = iniparser_getint(ini, inla_string_join(secname, "CONTROL.VB.HESSIAN.UPDATE"), 1);

	opt = Strdup(iniparser_getstring(ini, inla_string_join(secname, "CONTROL.VB.STRATEGY"), Strdup("MEAN")));
	if (!strcasecmp(opt, "MEAN")) {
		mb->ai_par->vb_strategy = GMRFLib_AI_VB_MEAN;
	} else if (!strcasecmp(opt, "VARIANCE")) {
		mb->ai_par->vb_strategy = GMRFLib_AI_VB_VARIANCE;
	} else {
		inla_error_field_is_void(__GMRFLib_FuncName, secname, "control.vb.strategy", opt);
	}

	opt = Strdup(iniparser_getstring(ini, inla_string_join(secname, "CONTROL.VB.HESSIAN.STRATEGY"), Strdup("FULL")));
	if (!strcasecmp(opt, "FULL")) {
		mb->ai_par->vb_hessian_strategy = GMRFLib_VB_HESSIAN_STRATEGY_FULL;
	} else if (!strcasecmp(opt, "PARTIAL")) {
		mb->ai_par->vb_hessian_strategy = GMRFLib_VB_HESSIAN_STRATEGY_PARTIAL;
	} else if (!strcasecmp(opt, "DIAGONAL")) {
		mb->ai_par->vb_hessian_strategy = GMRFLib_VB_HESSIAN_STRATEGY_DIAGONAL;
	} else {
		inla_error_field_is_void(__GMRFLib_FuncName, secname, "control.vb.hessian.strategy", opt);
	}

	// default value is set in main()
	GMRFLib_aqat_m_diag_add = iniparser_getdouble(ini, inla_string_join(secname, "CONSTR.MARGINAL.DIAGONAL"), GMRFLib_aqat_m_diag_add);
	assert(GMRFLib_aqat_m_diag_add >= 0.0);
	if (mb->verbose) {
		printf("\t\tconstr.marginal.diagonal = %.3g\n", GMRFLib_aqat_m_diag_add);
	}

	mb->ai_par->improved_simplified_laplace = iniparser_getboolean(ini, inla_string_join(secname, "IMPROVED.SIMPLIFIED.LAPLACE"), 0);
	mb->ai_par->parallel_linesearch = iniparser_getboolean(ini, inla_string_join(secname, "PARALLEL.LINESEARCH"), 0);
	mb->ai_par->hessian_correct_skewness_only = iniparser_getboolean(ini, inla_string_join(secname, "HESSIAN.CORRECT.SKEWNESS.ONLY"), 0);
	mb->compute_initial_values = iniparser_getboolean(ini, inla_string_join(secname, "COMPUTE.INITIAL.VALUES"), 1);

	if (mb->verbose) {
		GMRFLib_print_ai_param(stdout, mb->ai_par);
	}

	return INLA_OK;
}

int inla_parse_update(inla_tp *mb, dictionary *ini, int sec, int UNUSED(make_dir))
{
	/*
	 * parse section = UPDATE
	 */
	char *secname = NULL, *filename = NULL;
	GMRFLib_matrix_tp *M = NULL;

	if (mb->verbose) {
		printf("\tinla_parse_update...\n");
	}
	secname = Strdup(iniparser_getsecname(ini, sec));
	if (mb->verbose) {
		printf("\t\tsection[%s]\n", secname);
	}

	filename = Strdup(iniparser_getstring(ini, inla_string_join(secname, "FILENAME"), NULL));
	if (mb->verbose) {
		printf("\t\tfilename[%s]\n", (filename ? filename : "NULL"));
	}

	if (filename) {
		M = GMRFLib_read_fmesher_file(filename, (long int) 0, -1);
		mb->update = Calloc(1, inla_update_tp);
		int i = 0, j = 0, nt, k, kk;

		mb->update->ntheta = nt = (int) GMRFLib_matrix_get(i++, j, M);
		if (mb->verbose) {
			printf("\t\tntheta = %1d\n", nt);
		}

		mb->update->theta_mode = Calloc(nt, double);
		for (k = 0; k < nt; k++) {
			mb->update->theta_mode[k] = GMRFLib_matrix_get(i++, j, M);
			if (mb->verbose) {
				printf("\t\ttheta.mode[%1d] = %.10g\n", k, mb->update->theta_mode[k]);
			}
		}
		mb->update->stdev_corr_pos = Calloc(nt, double);
		for (k = 0; k < nt; k++) {
			mb->update->stdev_corr_pos[k] = GMRFLib_matrix_get(i++, j, M);
			if (mb->verbose) {
				printf("\t\tstdev.corr.pos[%1d] = %.10g\n", k, mb->update->stdev_corr_pos[k]);
			}
		}
		mb->update->stdev_corr_neg = Calloc(nt, double);
		for (k = 0; k < nt; k++) {
			mb->update->stdev_corr_neg[k] = GMRFLib_matrix_get(i++, j, M);
			if (mb->verbose) {
				printf("\t\tstdev.corr.neg[%1d] = %.10g\n", k, mb->update->stdev_corr_neg[k]);
			}
		}
		mb->update->sqrt_eigen_values = gsl_vector_calloc((size_t) nt);
		for (k = 0; k < nt; k++) {
			gsl_vector_set(mb->update->sqrt_eigen_values, k, GMRFLib_matrix_get(i++, j, M));
			if (mb->verbose) {
				printf("\t\tsqrt.eigen.values[%1d] = %.10g\n", k, gsl_vector_get(mb->update->sqrt_eigen_values, k));
			}
		}
		mb->update->eigen_vectors = gsl_matrix_calloc((size_t) nt, (size_t) nt);
		for (kk = 0; kk < nt; kk++) {
			for (k = 0; k < nt; k++) {
				/*
				 * column based storage...
				 */
				gsl_matrix_set(mb->update->eigen_vectors, k, kk, GMRFLib_matrix_get(i++, j, M));
				if (mb->verbose) {
					printf("\t\teigenvectors[%1d,%1d] = %.10g\n", k, kk, gsl_matrix_get(mb->update->eigen_vectors, k, kk));
				}
			}
		}
		assert(1 == M->ncol);
		assert(i == M->nrow);
		GMRFLib_matrix_free(M);
	}

	return INLA_OK;
}

int inla_parse_pardiso(inla_tp *mb, dictionary *ini, int sec, int UNUSED(make_dir))
{
	/*
	 * parse section = PARDISO
	 */
	char *secname = NULL;
	int val;

	if (mb->verbose) {
		printf("\tinla_parse_pardiso...\n");
	}
	secname = Strdup(iniparser_getsecname(ini, sec));
	if (mb->verbose) {
		printf("\t\tsection[%s]\n", secname);
	}

	val = iniparser_getint(ini, inla_string_join(secname, "VERBOSE"), 0);
	if (mb->verbose) {
		printf("\t\tverbose[%1d]\n", val);
	}
	GMRFLib_pardiso_set_verbose(val);

	val = iniparser_getint(ini, inla_string_join(secname, "DEBUG"), 0);
	if (mb->verbose) {
		printf("\t\tdebug[%1d]\n", val);
	}
	GMRFLib_pardiso_set_debug(val);

	val = iniparser_getint(ini, inla_string_join(secname, "PARALLEL.REORDERING"), 0);
	if (mb->verbose) {
		printf("\t\tparallel.reordering[%1d]\n", val);
	}
	GMRFLib_pardiso_set_parallel_reordering(val);

	val = iniparser_getint(ini, inla_string_join(secname, "NRHS"), -1);
	if (mb->verbose) {
		printf("\t\tnrhs[%1d]\n", val);
	}
	GMRFLib_pardiso_set_nrhs(val);

	return INLA_OK;
}

int inla_parse_lp_scale(inla_tp *mb, dictionary *ini, int sec, int UNUSED(make_dir))
{
	/*
	 * parse section = LP.SCALE
	 */
	int i, k;
	char *secname = NULL, *ctmp = NULL, *msg = NULL;
	double tmp = 0.0;
	Data_section_tp *ds = &(mb->data_sections[0]);

	if (mb->verbose) {
		printf("\tinla_parse_lp_scale...\n");
	}
	secname = Strdup(iniparser_getsecname(ini, sec));
	if (mb->verbose) {
		printf("\t\tsection[%s]\n", secname);
	}

	if (!ds->lp_scale) {
		return INLA_OK;
	}

	ds->lp_scale_in_use = Calloc(INLA_LP_SCALE_MAX, int);

	for (i = 0; i < INLA_LP_SCALE_MAX; i++) {
		ds->lp_scale_in_use[i] = 0;
	}
	for (i = 0; i < mb->predictor_ndata; i++) {
		if ((k = (int) ds->lp_scale[i]) >= 0) {
			ds->lp_scale_in_use[k] = 1;
			GMRFLib_ASSERT(k < INLA_LP_SCALE_MAX, GMRFLib_EPARAMETER);
		}
	}
	if (mb->verbose) {
		printf("\t\tlp_scale variables in use = [");
		for (k = 0; k < INLA_LP_SCALE_MAX; k++) {
			if (ds->lp_scale_in_use[k]) {
				printf(" %1d", k);
			}
		}
		printf(" ]\n");
	}

	// mark all hyperpar defs as read
	for (i = 0; i < INLA_LP_SCALE_MAX; i++) {
		for (int j = 0; j < keywords_len; j++) {
			GMRFLib_sprintf(&ctmp, "%s%1d", keywords[j], i);
			iniparser_getstring(ini, inla_string_join(secname, ctmp), NULL);
			Free(ctmp);
		}
	}

	ds->lp_scale_beta = Calloc(INLA_LP_SCALE_MAX, double **);
	ds->lp_scale_nfixed = Calloc(INLA_LP_SCALE_MAX, int);
	ds->lp_scale_nprior = Calloc(INLA_LP_SCALE_MAX, Prior_tp);

	for (k = 0; k < INLA_LP_SCALE_MAX; k++) {
		if (ds->lp_scale_in_use[k]) {
			GMRFLib_sprintf(&ctmp, "INITIAL%1d", k);
			tmp = iniparser_getdouble(ini, inla_string_join(secname, ctmp), 0.0);

			Free(ctmp);
			GMRFLib_sprintf(&ctmp, "FIXED%1d", k);
			ds->lp_scale_nfixed[k] = iniparser_getboolean(ini, inla_string_join(secname, ctmp), 0);

			if (!ds->lp_scale_nfixed[k] && mb->mode_use_mode) {
				tmp = mb->theta_file[mb->theta_counter_file++];
				if (mb->mode_fixed) {
					ds->lp_scale_nfixed[k] = 1;
				}
			}
			HYPER_NEW(ds->lp_scale_beta[k], tmp);
			if (mb->verbose) {
				printf("\t\tinitialise lp_scale_beta[%1d] = %g\n", k, ds->lp_scale_beta[k][0][0]);
				printf("\t\tfixed = %1d\n", ds->lp_scale_nfixed[k]);
			}
			inla_read_priorN(mb, ini, sec, &(ds->lp_scale_nprior[k]), "GAUSSIAN", k, NULL);

			if (!ds->lp_scale_nfixed[k]) {
				mb->theta = Realloc(mb->theta, mb->ntheta + 1, double **);
				mb->theta_hyperid = Realloc(mb->theta_hyperid, mb->ntheta + 1, char *);
				mb->theta_hyperid[mb->ntheta] = ds->lp_scale_nprior[k].hyperid;
				mb->theta_tag = Realloc(mb->theta_tag, mb->ntheta + 1, char *);
				mb->theta_tag_userscale = Realloc(mb->theta_tag_userscale, mb->ntheta + 1, char *);
				mb->theta_dir = Realloc(mb->theta_dir, mb->ntheta + 1, char *);

				Free(ctmp);
				GMRFLib_sprintf(&ctmp, "beta[%1d] for lp_scale", k + 1);
				mb->theta_tag[mb->ntheta] = ctmp;
				mb->theta_tag_userscale[mb->ntheta] = ctmp;
				GMRFLib_sprintf(&msg, "%s-parameter", secname);
				mb->theta_dir[mb->ntheta] = msg;

				mb->theta_from = Realloc(mb->theta_from, mb->ntheta + 1, char *);
				mb->theta_to = Realloc(mb->theta_to, mb->ntheta + 1, char *);
				mb->theta_from[mb->ntheta] = Strdup(ds->lp_scale_nprior[k].from_theta);
				mb->theta_to[mb->ntheta] = Strdup(ds->lp_scale_nprior[k].to_theta);

				mb->theta[mb->ntheta] = ds->lp_scale_beta[k];
				mb->theta_map = Realloc(mb->theta_map, mb->ntheta + 1, map_func_tp *);
				mb->theta_map[mb->ntheta] = map_identity;
				mb->theta_map_arg = Realloc(mb->theta_map_arg, mb->ntheta + 1, void *);
				mb->theta_map_arg[mb->ntheta] = NULL;
				mb->ntheta++;
			}
		} else {
			ds->lp_scale_nfixed[k] = 1;
		}
	}

	for (i = 1; i < mb->nds; i++) {
		mb->data_sections[i].lp_scale = mb->data_sections[0].lp_scale;
		mb->data_sections[i].lp_scale_in_use = mb->data_sections[0].lp_scale_in_use;
		mb->data_sections[i].lp_scale_beta = mb->data_sections[0].lp_scale_beta;
		mb->data_sections[i].lp_scale_nfixed = mb->data_sections[0].lp_scale_nfixed;
		mb->data_sections[i].lp_scale_nprior = mb->data_sections[0].lp_scale_nprior;
	}

	return INLA_OK;
}

int inla_parse_expert(inla_tp *mb, dictionary *ini, int sec)
{
	/*
	 * parse section = expert
	 */
	char *secname = NULL;

	if (mb->verbose) {
		printf("\tinla_parse_expert...\n");
	}
	secname = Strdup(iniparser_getsecname(ini, sec));
	if (mb->verbose) {
		printf("\t\tsection[%s]\n", secname);
	}

	mb->expert_disable_gaussian_check = iniparser_getint(ini, inla_string_join(secname, "DISABLE.GAUSSIAN.CHECK"), 0);
	if (mb->verbose) {
		printf("\t\t\tdisable.gaussian.check=[%1d]\n", mb->expert_disable_gaussian_check);
	}

	int dot_product_gain = iniparser_getboolean(ini, inla_string_join(secname, "DOT.PRODUCT.GAIN"), 0);
	// >=0 will measure, <0 will not (default off)
	GMRFLib_dot_product_gain = (dot_product_gain ? 0.0 : -1.0);
	if (mb->verbose) {
		printf("\t\t\tMeasure dot.product.gain=[%s]\n", (dot_product_gain ? "Yes" : "No"));
	}

	GMRFLib_opt_solve = iniparser_getboolean(ini, inla_string_join(secname, "OPT.SOLVE"), 0);
	if (mb->verbose) {
		printf("\t\t\tOptimise linear solve=[%s]\n", (GMRFLib_opt_solve ? "Yes" : "No"));
	}

	/*
	 * do error-checking later on 
	 */
	mb->expert_cpo_manual = iniparser_getint(ini, inla_string_join(secname, "CPO.MANUAL"), 0);

	char *str = NULL;
	str = iniparser_getstring(ini, inla_string_join(secname, "CPO.IDX"), str);

	int n = 0;
	int *idx = NULL;
	inla_sread_ints_q(&idx, &n, (const char *) str);

	mb->expert_n_cpo_idx = n;
	mb->expert_cpo_idx = idx;

	if (mb->verbose) {
		int i;

		printf("\t\t\tcpo.manual=[%1d]\n", mb->expert_cpo_manual);
		for (i = 0; i < mb->expert_n_cpo_idx; i++) {
			printf("\t\t\tcpo.idx=[%1d]\n", mb->expert_cpo_idx[i]);
		}
	}

	/*
	 * joint prior?
	 */
	char *file = NULL, *model = NULL;

	file = iniparser_getstring(ini, inla_string_join(secname, "JP.FILE"), file);
	model = iniparser_getstring(ini, inla_string_join(secname, "JP.MODEL"), model);
	if (mb->verbose) {
		printf("\t\t\tjp.file=[%s]\n", file);
		printf("\t\t\tjp.model=[%s]\n", model);
	}
	if (model) {
		mb->jp = Calloc(1, inla_jp_tp);
		mb->jp->file = Strdup(file);
		mb->jp->model = Strdup(model);
		if (R_load_INLA) {
			inla_R_library("INLA");		       /* initialize here */
			R_load_INLA = 0;
		}
	} else {
		mb->jp = NULL;
	}

	Free(file);
	file = iniparser_getstring(ini, inla_string_join(secname, "GLOBALCONSTR.A.FILE"), NULL);
	if (file) {
		mb->global_constr = Calloc(2, GMRFLib_matrix_tp *);
		mb->global_constr[0] = GMRFLib_read_fmesher_file(file, 0, -1);
		assert(mb->global_constr[0]);

		file = iniparser_getstring(ini, inla_string_join(secname, "GLOBALCONSTR.E.FILE"), NULL);
		mb->global_constr[1] = GMRFLib_read_fmesher_file(file, 0, -1);
		assert(mb->global_constr[1]);
		printf("\t\t\tnumber of global.constr=[%1d]\n", mb->global_constr[1]->nrow);
		assert(GMRFLib_inla_mode == GMRFLib_MODE_COMPACT);
	}

	return INLA_OK;
}
