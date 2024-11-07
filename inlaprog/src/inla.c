
/* inla.c
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

#ifndef GITCOMMIT
#define GITCOMMIT "devel"
#endif

#if defined(__sun__)
#include <stdlib.h>
#endif
#if defined(__linux__)
#include <getopt.h>
#endif
#include <float.h>
#include <assert.h>
#include <ctype.h>
#include <errno.h>
#include <math.h>
#include <signal.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <time.h>
#include <ltdl.h>

#if !defined(WINDOWS)
#include <sys/resource.h>
#endif

#if defined(__APPLE__)
#include <sys/types.h>
#include <sys/sysctl.h>
#endif

#include "rmath.h"

#include "GMRFLib/GMRFLib.h"
#include "GMRFLib/GMRFLibP.h"
#if !defined(ISNAN)
#define ISNAN(x) (isnan(x)!=0)
#endif

#if !defined(INLA_TAG)
#define INLA_TAG "devel"
#endif

#include <unistd.h>
#include <stdlib.h>
#if defined(WIN32) || defined(WINDOWS)
#include <windows.h>
#include <direct.h>
#endif

#include "inla.h"
#include "my.h"
#include "spde.h"
#include "spde2.h"
#include "spde3.h"
#include "eval.h"
#include "ar.h"
#include "pc-priors.h"
#include "R-interface.h"
#include "fgn.h"
#include "tweedie.h"
#include "pc-powerlink.h"
#include "link-gev.h"
#include "cgeneric.h"

#define PREVIEW (10)
#define MODEFILENAME ".inla-mode"
#define MODEFILENAME_FMT "%02x"


// as given in models.R 
#define TSTRATA_MAXTHETA (11L)
#define SPDE2_MAXTHETA (100L)
#define SPDE3_MAXTHETA (100L)
#define GENERIC3_MAXTHETA (11L)
#define AR_MAXTHETA (10L)
#define LINK_MAXTHETA (10L)
#define STRATA_MAXTHETA (10L)
#define NMIX_MMAX (15L)
#define POM_MAXTHETA (10L)
#define INTSLOPE_MAXTHETA (50L)
#define BGEV_MAXTHETA (10L)
#define POISSON0_MAXTHETA (10L)
#define BINOMIAL0_MAXTHETA (10L)
#define GGAUSSIAN_MAXTHETA (10L)
#define CURE_MAXTHETA (10L)
#define SCOPY_MAXTHETA (15L)
#define RCPOISSON_MAXTHETA (10L)
#define TPOISSON_MAXTHETA (10L)
#define OCCUPANCY_MAXTHETA (10L)
#define BINOMIALMIX_NBETA (9L)
#define L_FL_NC (9L)

G_tp G = { 1, INLA_MODE_DEFAULT, 4.0, 0.5, 2, 0, GMRFLib_REORDER_DEFAULT, 0, 0 };

const int keywords_len = 7;
const char *keywords[] = {
	"FIXED", "INITIAL", "PRIOR", "HYPERID", "PARAMETERS", "TO.THETA", "FROM.THETA", NULL
};

// defined in R-interface.c
extern double R_rgeneric_cputime;

int G_norm_const_len = 0;
double *G_norm_const = NULL;				       /* store static normalization constants for likelihoods */
void **G_norm_const_v = NULL;
char *G_norm_const_compute = NULL;			       /* to be computed */

int R_load_INLA = 0;

/* 
   default values for priors
 */
#define DEFAULT_GAMMA_PRIOR_A          1.0
#define DEFAULT_GAMMA_PRIOR_B          0.00005
#define DEFAULT_NORMAL_PRIOR_PRECISION 0.001

/* 
   these are the same, just that the interface is cleaner with the `ds'
 */
#define OFFSET(idx_) ds->offset[idx_]
#define OFFSET3(idx_) mb->offset[idx_]

#define LINK_INIT							\
	double off = OFFSET(idx);					\
	double *_link_covariates = NULL;				\
	Link_param_tp *predictor_invlinkfunc_arg = (Link_param_tp *) (ds->predictor_invlinkfunc_arg[idx]); \
	if (ds->link_covariates) {					\
		_link_covariates = Calloc(ds->link_covariates->ncol, double); \
		GMRFLib_matrix_get_row(_link_covariates, idx, ds->link_covariates); \
	}								\
	double _lp_scale = 1.0;						\
	if (ds->lp_scale && ds->lp_scale[idx] >= 0) {					\
		_lp_scale = ds->lp_scale_beta[(int)ds->lp_scale[idx]][thread_id][0]; \
	}


#define LINK_END Free(_link_covariates)
#define PREDICTOR_SCALE _lp_scale
#define PREDICTOR_LINK_EQ(_fun) (ds->predictor_invlinkfunc == (_fun))
#define PREDICTOR_SIMPLE_LINK_EQ(_fun) (ds->data_observations.link_simple_invlinkfunc ==  (_fun))
#define PREDICTOR_INVERSE_LINK(xx_)					\
	ds->predictor_invlinkfunc(thread_id, _lp_scale * (xx_), MAP_FORWARD, (void *)predictor_invlinkfunc_arg, _link_covariates)
#define PREDICTOR_INVERSE_LINK_NO_SCALE(xx_)				\
	ds->predictor_invlinkfunc(thread_id, (xx_), MAP_FORWARD, (void *)predictor_invlinkfunc_arg, _link_covariates)
#define PREDICTOR_INVERSE_IDENTITY_LINK(xx_) (_lp_scale * (xx_))
#define PREDICTOR_LINK(xx_)						\
	(ds->predictor_invlinkfunc(thread_id, (xx_), MAP_BACKWARD, (void *)predictor_invlinkfunc_arg, _link_covariates) / _lp_scale)
#define PREDICTOR_INVERSE_LINK_LOGJACOBIAN(xx_)  \
	log(ABS(_lp_scale * ds->predictor_invlinkfunc(thread_id, _lp_scale * (xx_), MAP_DFORWARD, (void *)predictor_invlinkfunc_arg, _link_covariates)))

#define PENALTY (-100.0)				       /* wishart3d: going over limit... */

#define READ(fd, buf, num, type)					\
	if (1) {							\
		if (num > 0) {						\
			int bytes_read;					\
			bytes_read = read(fd, (void *) buf, (size_t) ((num)*sizeof(type))); \
			assert(bytes_read == (int) ((num)*sizeof(type))); \
		}							\
	}

#define WRITE(fd, buf, num, type)					\
	if (1) {							\
		if (num > 0) {						\
			int bytes_written;				\
			bytes_written = write(fd, (void *) buf, (size_t) ((num)*sizeof(type))); \
			assert(bytes_written == (int) ((num)*sizeof(type))); \
		}							\
	}

#define WRITE_MSG(fd, Id, msg)					\
	if (1) {						\
		int _id = (Id);					\
		WRITE(fd, &_id, 1, int);			\
		char *dup_msg = Strdup(msg);		\
		WRITE(fd, dup_msg, strlen(dup_msg)+1, char);	\
		Free(dup_msg);					\
	}

#define CODE_NEEDED \
	if (1) {							\
		fprintf(stderr, "\n\n%s: %d  CODE NEEDED HERE!\n\n", __FILE__, __LINE__); \
		exit(1);						\
		assert(0==1);						\
	}

// these versions (in my.c) cache the result when the argument is integer
#define gsl_sf_lngamma(_x) lgamma(_x)
#define gsl_sf_lnchoose_e(_a, _b, _c) my_gsl_sf_lnchoose_e(_a, _b, _c)

#include "inla-sys.c"
#include "inla-map-and-link.c"
#include "inla-testit.c"
#include "inla-graph.c"
#include "inla-Qfunc.c"
#include "inla-likelihood.c"
#include "inla-priors.c"
#include "inla-special-functions.c"
#include "inla-output.c"
#include "inla-modes.c"
#include "inla-classic.c"
#include "inla-read.c"
#include "inla-parse.c"
#include "param-constr.c"

double inla_interpolate_mode(double *x, double *y)
{
	// give 3 values (x,y), return the mode. truncate at the boundary

	double xm = (y[0] * x[1] * x[1] - y[0] * x[2] * x[2] - y[1] * x[0] * x[0] + y[1] * x[2] * x[2] + y[2] * x[0] * x[0] - y[2] * x[1] * x[1]) /
	    (y[0] * x[1] - y[0] * x[2] - y[1] * x[0] + y[1] * x[2] + y[2] * x[0] - x[1] * y[2]) / 0.2e1;

	double xmin = DMIN(x[0], DMIN(x[1], x[2]));
	double xmax = DMAX(x[0], DMAX(x[1], x[2]));
	xm = TRUNCATE(xm, xmin, xmax);

	return (xm);
}


inla_tp *inla_build(const char *dict_filename, int verbose, int make_dir)
{
	/*
	 * This function builds the model from the contents in INI 
	 */
	int found, sec, nsec, count, len, i, idx, j, k = -1;
	char *secname = NULL, *sectype = NULL, *sec_read = NULL, *msg = NULL;
	dictionary *ini = NULL;
	inla_tp *mb = NULL;

	if (verbose) {
		printf("%s...\n", __GMRFLib_FuncName);
	}
	mb = Calloc(1, inla_tp);
	mb->verbose = verbose;
	mb->mode_restart = 1;
	mb->mode_fixed = mb->mode_use_mode = 0;

	ini = iniparser_load(dict_filename);
	if (!ini) {
		GMRFLib_sprintf(&msg, "Fail to parse ini-file[%s]....", dict_filename);
		inla_error_general(msg);
	}
	nsec = iniparser_getnsec(ini);
	if (mb->verbose) {
		printf("\tnumber of sections=[%1d]\n", nsec);
	}
	/*
	 * first check that "type" is present in each section 
	 */
	for (sec = 0; sec < nsec; sec++) {
		secname = Strdup(iniparser_getsecname(ini, sec));
		sectype = Strdup(strupc(iniparser_getstring(ini, inla_string_join((const char *) secname, "type"), NULL)));
		if (!sectype) {
			inla_error_missing_required_field(__GMRFLib_FuncName, secname, "type");
		}
		Free(secname);
		Free(sectype);
	}
	sec_read = Calloc(nsec, char);

	/*
	 * default: gaussian data is on, then its turned off... unless we chose to disable the check
	 */
	if (mb->expert_disable_gaussian_check) {
		GMRFLib_gaussian_data = GMRFLib_FALSE;
	} else {
		GMRFLib_gaussian_data = GMRFLib_TRUE;
	}

	/*
	 * ...then parse the sections in this order: RLIB, EXPERT, MODE, PROBLEM, PREDICTOR, DATA, FFIELD, LINEAR, INLA, UPDATE, LINCOMB, OUTPUT
	 * 
	 * it is easier to do it like this, instead of insisting the user to write the section in a spesific order.
	 * 
	 */
	for (sec = found = 0; sec < nsec; sec++) {
		secname = Strdup(iniparser_getsecname(ini, sec));
		sectype = Strdup(strupc(iniparser_getstring(ini, inla_string_join((const char *) secname, "TYPE"), NULL)));
		if (!strcmp(sectype, "LIBR")) {
			if (mb->verbose) {
				printf("\tparse section=[%1d] name=[%s] type=[LIBR]\n", sec, iniparser_getsecname(ini, sec));
			}
			if (found++) {
				GMRFLib_sprintf(&msg, "%s: two or more sections of type = [LIBR]. Exit.\n", __GMRFLib_FuncName);
				inla_error_general(msg);
			}
			sec_read[sec] = 1;
			inla_parse_libR(mb, ini, sec);
		}
		Free(secname);
		Free(sectype);
	}

	for (sec = found = 0; sec < nsec; sec++) {
		secname = Strdup(iniparser_getsecname(ini, sec));
		sectype = Strdup(strupc(iniparser_getstring(ini, inla_string_join((const char *) secname, "TYPE"), NULL)));
		if (!strcmp(sectype, "EXPERT")) {
			if (mb->verbose) {
				printf("\tparse section=[%1d] name=[%s] type=[EXPERT]\n", sec, iniparser_getsecname(ini, sec));
			}
			if (found++) {
				GMRFLib_sprintf(&msg, "%s: two or more sections of type = [EXPERT]. Exit.\n", __GMRFLib_FuncName);
				inla_error_general(msg);
			}
			sec_read[sec] = 1;
			inla_parse_expert(mb, ini, sec);
		}
		Free(secname);
		Free(sectype);
	}

	for (sec = found = 0; sec < nsec; sec++) {
		secname = Strdup(iniparser_getsecname(ini, sec));
		sectype = Strdup(strupc(iniparser_getstring(ini, inla_string_join((const char *) secname, "TYPE"), NULL)));
		if (!strcmp(sectype, "MODE")) {
			if (mb->verbose) {
				printf("\tparse section=[%1d] name=[%s] type=[MODE]\n", sec, iniparser_getsecname(ini, sec));
			}
			if (found++) {
				GMRFLib_sprintf(&msg, "%s: two or more sections of type = [MODE]. Exit.\n", __GMRFLib_FuncName);
				inla_error_general(msg);
			}
			sec_read[sec] = 1;
			inla_parse_mode(mb, ini, sec);
		}
		Free(secname);
		Free(sectype);
	}

	for (sec = found = 0; sec < nsec; sec++) {
		secname = Strdup(iniparser_getsecname(ini, sec));
		sectype = Strdup(strupc(iniparser_getstring(ini, inla_string_join((const char *) secname, "TYPE"), NULL)));
		if (!strcmp(sectype, "PROBLEM")) {
			if (mb->verbose) {
				printf("\tparse section=[%1d] name=[%s] type=[PROBLEM]\n", sec, iniparser_getsecname(ini, sec));
			}
			if (found++) {
				GMRFLib_sprintf(&msg, "%s: two or more sections of type = [PROBLEM]. Exit.\n", __GMRFLib_FuncName);
				inla_error_general(msg);
			}
			sec_read[sec] = 1;
			inla_parse_problem(mb, ini, sec, make_dir);
		}
		Free(secname);
		Free(sectype);
	}
	if (!found) {
		GMRFLib_sprintf(&msg, "%s: no section of type = [PROBLEM]", __GMRFLib_FuncName);
		inla_error_general(msg);
	}

	/*
	 * type = predictor 
	 */
	for (sec = found = 0; sec < nsec; sec++) {
		secname = Strdup(iniparser_getsecname(ini, sec));
		sectype = Strdup(strupc(iniparser_getstring(ini, inla_string_join((const char *) secname, "TYPE"), NULL)));
		if (!strcmp(sectype, "PREDICTOR")) {
			if (mb->verbose) {
				printf("\tparse section=[%1d] name=[%s] type=[PREDICTOR]\n", sec, iniparser_getsecname(ini, sec));
			}
			if (found++) {
				GMRFLib_sprintf(&msg, "%s: two or more sections of type = [PREDICTOR]. Exit.\n", __GMRFLib_FuncName);
				inla_error_general(msg);
			}
			sec_read[sec] = 1;
			inla_parse_predictor(mb, ini, sec);
		}
		Free(secname);
		Free(sectype);
	}
	if (!found) {
		GMRFLib_sprintf(&msg, "%s: no section of type = [PREDICTOR]", __GMRFLib_FuncName);
		inla_error_general(msg);
	}

	/*
	 * type = DATA 
	 */
	mb->ds = 0;
	for (sec = 0; sec < nsec; sec++) {
		secname = Strdup(iniparser_getsecname(ini, sec));
		sectype = Strdup(strupc(iniparser_getstring(ini, inla_string_join((const char *) secname, "TYPE"), NULL)));
		if (!strcmp(sectype, "DATA")) {
			if (mb->verbose) {
				printf("\tparse section=[%1d] name=[%s] type=[DATA]\n", sec, iniparser_getsecname(ini, sec));
			}
			sec_read[sec] = 1;
			inla_parse_data(mb, ini, sec);
			mb->ds++;
		}
		Free(secname);
		Free(sectype);
	}
	if (!found) {
		GMRFLib_sprintf(&msg, "%s: no section of type [DATA] found", __GMRFLib_FuncName);
		inla_error_general(msg);
	}

	found = 0;
	/*
	 * type = ffield 
	 */
	for (sec = 0; sec < nsec; sec++) {
		secname = Strdup(iniparser_getsecname(ini, sec));
		sectype = Strdup(strupc(iniparser_getstring(ini, inla_string_join((const char *) secname, "TYPE"), NULL)));
		if (!strcmp(sectype, "FFIELD")) {
			if (mb->verbose) {
				printf("\tparse section=[%1d] name=[%s] type=[FFIELD]\n", sec, iniparser_getsecname(ini, sec));
			}
			found++;
			sec_read[sec] = 1;
			inla_parse_ffield(mb, ini, sec);
		}
		Free(secname);
		Free(sectype);
	}

	inla_add_copyof(mb);
	inla_add_scopyof(mb);

	/*
	 * type = linear 
	 */
	for (sec = 0; sec < nsec; sec++) {
		secname = Strdup(iniparser_getsecname(ini, sec));
		sectype = Strdup(strupc(iniparser_getstring(ini, inla_string_join((const char *) secname, "TYPE"), NULL)));
		if (!strcmp(sectype, "LINEAR")) {
			if (mb->verbose) {
				printf("\tsection=[%1d] name=[%s] type=[LINEAR]\n", sec, iniparser_getsecname(ini, sec));
			}
			found++;
			sec_read[sec] = 1;
			inla_parse_linear(mb, ini, sec);
		}
		Free(secname);
		Free(sectype);
	}

	/*
	 * type = UPDATE
	 */
	inla_setup_ai_par_default(mb);			       /* need this if there is no INLA section */
	mb->ai_par->mode_fixed = mb->mode_fixed;
	mb->ai_par->mode_restart = mb->mode_restart;
	mb->ai_par->mode_use_mode = mb->mode_use_mode;

	for (sec = 0; sec < nsec; sec++) {
		secname = Strdup(iniparser_getsecname(ini, sec));
		sectype = Strdup(strupc(iniparser_getstring(ini, inla_string_join((const char *) secname, "TYPE"), NULL)));
		if (!strcmp(sectype, "UPDATE")) {
			if (mb->verbose) {
				printf("\tparse section=[%1d] name=[%s] type=[UPDATE]\n", sec, iniparser_getsecname(ini, sec));
			}
			sec_read[sec] = 1;
			inla_parse_update(mb, ini, sec, make_dir);
		}
		Free(secname);
		Free(sectype);
	}

	/*
	 * type = PARDISO
	 */
	for (sec = 0; sec < nsec; sec++) {
		secname = Strdup(iniparser_getsecname(ini, sec));
		sectype = Strdup(strupc(iniparser_getstring(ini, inla_string_join((const char *) secname, "TYPE"), NULL)));
		if (!strcmp(sectype, "PARDISO")) {
			if (mb->verbose) {
				printf("\tparse section=[%1d] name=[%s] type=[PARDISO]\n", sec, iniparser_getsecname(ini, sec));
			}
			sec_read[sec] = 1;
			inla_parse_pardiso(mb, ini, sec, make_dir);
		}
		Free(secname);
		Free(sectype);
	}

	/*
	 * type = LPSCALE
	 */
	for (sec = 0; sec < nsec; sec++) {
		secname = Strdup(iniparser_getsecname(ini, sec));
		sectype = Strdup(strupc(iniparser_getstring(ini, inla_string_join((const char *) secname, "TYPE"), NULL)));
		if (!strcmp(sectype, "LP.SCALE")) {
			if (mb->verbose) {
				printf("\tparse section=[%1d] name=[%s] type=[LP.SCALE]\n", sec, iniparser_getsecname(ini, sec));
			}
			sec_read[sec] = 1;
			inla_parse_lp_scale(mb, ini, sec, make_dir);
		}
		Free(secname);
		Free(sectype);
	}

	inla_parse_param_constraints(mb);

	/*
	 * build the index table and the hash; need this before reading the lincomb sections
	 */
	len = 1 + (mb->predictor_m > 0 ? 1 : 0) + mb->nf + mb->nlinear;
	mb->idx_tag = Calloc(len, char *);
	mb->idx_start = Calloc(len, int);
	mb->idx_n = Calloc(len, int);

	j = idx = 0;
	if (mb->predictor_m > 0) {
		mb->idx_tag[j] = Strdup(mb->Apredictor_tag);
		mb->idx_start[j] = idx;
		mb->idx_n[j] = mb->predictor_m;
		idx += mb->idx_n[j++];
	}
	mb->idx_tag[j] = Strdup(mb->predictor_tag);
	mb->idx_start[j] = idx;
	mb->idx_n[j] = mb->predictor_n;
	idx += mb->idx_n[j++];
	for (i = 0; i < mb->nf; i++) {
		mb->idx_tag[j] = Strdup(mb->f_tag[i]);
		mb->idx_start[j] = idx;
		mb->idx_n[j] = mb->f_Ntotal[i];
		idx += mb->idx_n[j++];
	}
	for (i = 0; i < mb->nlinear; i++) {
		mb->idx_tag[j] = Strdup(mb->linear_tag[i]);
		mb->idx_start[j] = idx;
		mb->idx_n[j] = 1;
		idx += mb->idx_n[j++];
	}
	mb->idx_tot = j;
	mb->idx_ntot = idx;

	map_stri_init_hint(&(mb->idx_hash), mb->idx_tot);
	for (i = 0; i < mb->idx_tot; i++) {
		map_stri_set(&(mb->idx_hash), Strdup(mb->idx_tag[i]), i);
	}

	if (mb->verbose) {
		printf("\n\tIndex table: number of entries[%1d], total length[%1d]\n", mb->idx_tot, mb->idx_ntot);
		printf("\t\t%-30s %10s %10s\n", "tag", "start-index", "length");
		for (i = 0; i < mb->idx_tot; i++) {
			printf("\t\t%-30s %10d %10d\n", mb->idx_tag[i], mb->idx_start[i], mb->idx_n[i]);
		}
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

	// copy ptr to these, in case we need to correct in strategy="prior" later with gcpo
	if (mb->gcpo_param) {
		mb->gcpo_param->idx_tot = mb->idx_tot;
		mb->gcpo_param->idx_tag = mb->idx_tag;
		mb->gcpo_param->idx_start = mb->idx_start;
		mb->gcpo_param->idx_n = mb->idx_n;
	}

	/*
	 * type = INLA 
	 */
	inla_setup_ai_par_default(mb);			       /* need this if there is no INLA section */
	for (sec = 0; sec < nsec; sec++) {
		secname = Strdup(iniparser_getsecname(ini, sec));
		sectype = Strdup(strupc(iniparser_getstring(ini, inla_string_join((const char *) secname, "TYPE"), NULL)));
		if (!strcmp(sectype, "INLA")) {
			if (mb->verbose) {
				printf("\tparse section=[%1d] name=[%s] type=[INLA]\n", sec, iniparser_getsecname(ini, sec));
			}
			sec_read[sec] = 1;
			inla_parse_INLA(mb, ini, sec, make_dir);
		}
		Free(secname);
		Free(sectype);
	}

	/*
	 * type = lincomb
	 */
	int numsec = 0, *secmap = NULL, isec = -1;

	/*
	 * first we count them
	 */
	for (sec = 0; sec < nsec; sec++) {
		secname = Strdup(iniparser_getsecname(ini, sec));
		sectype = Strdup(strupc(iniparser_getstring(ini, inla_string_join((const char *) secname, "TYPE"), NULL)));
		if (!strcmp(sectype, "LINCOMB")) {
			numsec++;
		}
	}
	if (numsec > 0) {
		/*
		 * then we find out which order to read them using 'LINCOMB.ORDER'
		 */
		secmap = Calloc(numsec, int);
		for (sec = 0; sec < nsec; sec++) {
			secname = Strdup(iniparser_getsecname(ini, sec));
			sectype = Strdup(strupc(iniparser_getstring(ini, inla_string_join((const char *) secname, "TYPE"), NULL)));
			if (!strcmp(sectype, "LINCOMB")) {
				int ordering;
				ordering = (int) iniparser_getdouble(ini, inla_string_join((const char *) secname, "LINCOMB.ORDER"), -1);
				GMRFLib_ASSERT_RETVAL(ordering > 0, GMRFLib_ESNH, (inla_tp *) NULL);
				secmap[ordering - 1] = sec;    /* ordering in the Model.ini is from 1...n */
			}
		}

		/*
		 * then we read them
		 */
		for (isec = 0; isec < numsec; isec++) {
			sec = secmap[isec];
			secname = Strdup(iniparser_getsecname(ini, sec));
			sectype = Strdup(strupc(iniparser_getstring(ini, inla_string_join((const char *) secname, "TYPE"), NULL)));
			if (!strcmp(sectype, "LINCOMB")) {
				/*
				 * we need to implement this here, as the number of linear combinations can get really huge and we need to surpress the verbose mode just
				 * for these sections. 
				 */
				int verbose_save = mb->verbose;

				// This option can surpress mb->verbose locally, but not the other way around.
				mb->verbose = iniparser_getint(ini, inla_string_join(secname, "VERBOSE"), mb->verbose)
				    && mb->verbose;

				if (mb->verbose) {
					printf("\tsection=[%1d] name=[%s] type=[LINCOMB]\n", sec, iniparser_getsecname(ini, sec));
				}
				found++;
				sec_read[sec] = 1;
				inla_parse_lincomb(mb, ini, sec);

				mb->verbose = verbose_save;    /* set it back */
			} else {
				/*
				 * should not happen
				 */
				GMRFLib_ASSERT_RETVAL(0 == 1, GMRFLib_ESNH, NULL);
			}
			Free(secname);
			Free(sectype);
		}
	}
	if (mb->verbose) {
		if (numsec) {
			printf("\tRead [%1d] sections with mode=[LINCOMB]\n", numsec);
		}
	}
	Free(secmap);

	/*
	 * check that all sections are read 
	 */
	for (sec = 0; sec < nsec; sec++) {
		if (!sec_read[sec]) {
			secname = Strdup(iniparser_getsecname(ini, sec));
			sectype = Strdup(strupc(iniparser_getstring(ini, inla_string_join((const char *) secname, "TYPE"), NULL)));
			GMRFLib_sprintf(&msg, "%s: section=[%s] is not used; please check its type=[%s]", __GMRFLib_FuncName, secname, sectype);
			inla_error_general(msg);
		}
	}
	if (mb->verbose) {
		printf("%s: check for unused entries in[%s]\n", __GMRFLib_FuncName, dict_filename);
	}
	if ((count = dictionary_dump_unused(ini, stderr))) {
		fprintf(stderr, "\n\ninla_build: [%s] contain[%1d] unused entries. PLEASE CHECK\n", dict_filename, count);
		exit(EXIT_FAILURE);
	}

	if (0) {
		dictionary_dump(ini, stdout);
		exit(0);
	}

	if (mb->mode_use_mode) {
		/*
		 * if the test fail, its a good idea to provide some debug information which might be helpful to help what is wrong in the spesification. 
		 */
		if (mb->theta_counter_file != mb->ntheta_file) {
			char *ctmp = NULL;
			GMRFLib_sprintf(&ctmp,
					"Your model has %1d hyperparameter(s) which is different from the %1d hyperparameter(s) given in 'control.mode'",
					mb->theta_counter_file, mb->ntheta_file);
			inla_error_general(ctmp);
			assert(mb->theta_counter_file == mb->ntheta_file);
		}
	}

	/*
	 * make the final likelihood from all the data-sections 
	 */
	mb->loglikelihood = Calloc(mb->predictor_ndata, GMRFLib_logl_tp *);
	mb->loglikelihood_arg = Calloc(mb->predictor_ndata, void *);
	mb->d = Calloc(mb->predictor_ndata, double);
	mb->family_idx = Calloc(mb->predictor_ndata, double);
	mb->len_family_idx = mb->predictor_ndata;

	for (i = 0; i < mb->predictor_ndata; i++) {
		for (j = found = 0; j < mb->nds; j++) {
			if (mb->data_sections[j].data_observations.d[i]) {
				k = j;
				found++;
			}
		}
		if (found > 1) {
			GMRFLib_sprintf(&msg, "Observation %d occurs in more than one data-section\n", i);
			inla_error_general(msg);
			exit(EXIT_FAILURE);
		}

		if (found) {
			mb->family_idx[i] = k;
			mb->loglikelihood[i] = mb->data_sections[k].loglikelihood;
			mb->loglikelihood_arg[i] = (void *) &(mb->data_sections[k]);
			mb->d[i] = mb->data_sections[k].data_observations.d[i];
		} else {
			mb->family_idx[i] = NAN;
			mb->loglikelihood[i] = NULL;
			mb->loglikelihood_arg[i] = NULL;
			mb->d[i] = 0.0;
		}
	}
	mb->data_ntheta_all = 0;
	for (j = 0; j < mb->nds; j++) {
		mb->data_ntheta_all += mb->data_sections[j].data_ntheta + mb->data_sections[j].mix_ntheta + mb->data_sections[j].link_ntheta;
		mb->data_sections[j].offset = mb->offset;      /* just a copy */
		mb->data_sections[j].mb = mb;		       /* just a copy */
	}

	/*
	 * make the final predictor_... from all the data-sections 
	 */
	int need_link = 0;
	mb->predictor_invlinkfunc = Calloc(mb->predictor_n + mb->predictor_m, link_func_tp *);
	mb->predictor_invlinkfunc_arg = Calloc(mb->predictor_n + mb->predictor_m, void *);
	mb->predictor_invlinkfunc_covariates = Calloc(mb->predictor_n + mb->predictor_m, GMRFLib_matrix_tp *);
	mb->predictor_family = Calloc(mb->predictor_n + mb->predictor_m, double);	/* as we use NAN */
	for (i = 0; i < mb->predictor_ndata; i++) {
		for (j = found = 0; j < mb->nds; j++) {
			if (mb->data_sections[j].data_observations.d[i]) {
				k = j;
				found++;
			}
		}
		if (found > 1) {
			GMRFLib_sprintf(&msg, "Observation %d occurs in more than one data-section\n", i);
			inla_error_general(msg);
			exit(EXIT_FAILURE);
		}
		mb->predictor_family[i] = (found == 1 ? (double) k :
					   ((mb->link_fitted_values && !gsl_isnan(mb->link_fitted_values[i])) ? mb->link_fitted_values[i] : NAN));
		mb->predictor_invlinkfunc[i] = (found == 1 ? mb->data_sections[k].predictor_invlinkfunc :
						((mb->link_fitted_values && !gsl_isnan(mb->link_fitted_values[i])) ?
						 mb->data_sections[(int) (mb->link_fitted_values[i])].predictor_invlinkfunc : NULL));
		mb->predictor_invlinkfunc_arg[i] = (found == 1 ? mb->data_sections[k].predictor_invlinkfunc_arg[i] :
						    ((mb->link_fitted_values && !gsl_isnan(mb->link_fitted_values[i])) ?
						     mb->data_sections[(int) (mb->link_fitted_values[i])].predictor_invlinkfunc_arg[i] : NULL));
		mb->predictor_invlinkfunc_covariates[i] = (found == 1 ? mb->data_sections[k].link_covariates :
							   ((mb->link_fitted_values && !gsl_isnan(mb->link_fitted_values[i])) ?
							    mb->data_sections[(int) (mb->link_fitted_values[i])].link_covariates : NULL));
		if (found == 0 && mb->predictor_invlinkfunc[i] == NULL)
			need_link++;
	}

	// fprintf(stderr, "%d\n", GMRFLib_gaussian_data);
	if (need_link && !GMRFLib_gaussian_data) {
		fprintf(stderr, "\n\n*** Warning *** You might want to consider to setting ``control.predictor=list(link=...)''\n");
		fprintf(stderr, "*** Warning *** otherwise the identity link will be used to compute the fitted values for NA data\n\n\n");
	}

	iniparser_freedict(ini);
	return mb;
}

GMRFLib_constr_tp *inla_make_constraint(int n, int sumzero, GMRFLib_constr_tp *constr)
{
	/*
	 * merge the constraints, if any. yes, i have to do this manually. 
	 */

	int i, j, nc;
	GMRFLib_constr_tp *c = NULL;

	if (!sumzero && !constr) {
		return NULL;
	}

	GMRFLib_make_empty_constr(&c);
	if (sumzero && !constr) {
		/*
		 * just the sumzero constraint 
		 */
		nc = 1;
		c->nc = nc;
		c->a_matrix = Calloc(nc * n, double);
		c->e_vector = Calloc(nc, double);
		for (i = 0; i < n; i++) {
			c->a_matrix[i] = 1.0;
		}
	} else if (!sumzero && constr) {
		nc = constr->nc;
		c->nc = nc;
		c->a_matrix = Calloc(nc * n, double);
		c->e_vector = Calloc(nc, double);
		Memcpy(c->a_matrix, constr->a_matrix, n * nc * sizeof(double));
		Memcpy(c->e_vector, constr->e_vector, nc * sizeof(double));
	} else {
		assert(sumzero && constr);

		nc = constr->nc + 1;
		c->nc = nc;
		c->a_matrix = Calloc(nc * n, double);
		c->e_vector = Calloc(nc, double);

		for (j = 0; j < nc - 1; j++) {
			for (i = 0; i < n; i++) {
				c->a_matrix[i * nc + j] = constr->a_matrix[i * constr->nc + j];
			}
			c->e_vector[j] = constr->e_vector[j];
		}
		j = nc - 1;
		for (i = 0; i < n; i++) {
			c->a_matrix[i * nc + j] = 1.0;
		}
		c->e_vector[j] = 0.0;
	}

	/*
	 * this is a stupid construction (i blame myself...): i have to pass the graph in order to pass `n'... 
	 */
	GMRFLib_graph_tp *g = NULL;

	GMRFLib_graph_mk_linear(&g, n, 0, 0);
	GMRFLib_prepare_constr(c, g, 0);
	GMRFLib_graph_free(g);

	return c;
}

GMRFLib_constr_tp *inla_make_constraint2(int n, int replicate, int sumzero, GMRFLib_constr_tp *constr)
{
	/*
	 * merge the constraints, if any. yes, i have to do this manually.
	 * this is the second version for which n is the basic size which has to be replicated.
	 */

	int i, ii, j, k, ccount, Ntotal;
	GMRFLib_constr_tp *c = NULL;

	if (!sumzero && !constr) {
		return NULL;
	}

	Ntotal = n * replicate;
	GMRFLib_make_empty_constr(&c);
	c->nc = ((constr ? constr->nc : 0) + (sumzero ? 1 : 0)) * replicate;
	c->a_matrix = Calloc(c->nc * Ntotal, double);
	c->e_vector = Calloc(c->nc, double);

	ccount = 0;
	if (sumzero) {
		for (j = 0; j < replicate; j++) {
			for (i = 0; i < n; i++) {
				ii = i + j * n;
				c->a_matrix[ii * c->nc + ccount] = 1.0;
			}
			ccount++;
		}
	}
	if (constr) {
		for (k = 0; k < constr->nc; k++) {
			for (j = 0; j < replicate; j++) {
				for (i = 0; i < n; i++) {
					ii = i + j * n;
					c->a_matrix[ii * c->nc + ccount] = constr->a_matrix[i * constr->nc + k];
				}
				c->e_vector[ccount] = constr->e_vector[k];
				ccount++;
			}
		}
	}
	assert(ccount == c->nc);

	/*
	 * this is a stupid construction (i blame myself...): i have to pass the graph in order to pass `n'... 
	 */
	GMRFLib_graph_tp *g = NULL;

	GMRFLib_graph_mk_linear(&g, Ntotal, 0, 0);
	GMRFLib_prepare_constr(c, g, 0);
	GMRFLib_graph_free(g);

	return c;
}


int inla_cgeneric_debug(FILE *fp, char *secname, inla_cgeneric_cmd_tp cmd, double *out)
{
	int i, n, m;

#pragma omp critical (Name_3bf40dfd2018962c499124a04e65c58fff8542a3)
	{
		fprintf(fp, "cgeneric[ %s ]  cmd = %s\n", secname, INLA_CGENERIC_CMD_NAME(cmd));
		switch (cmd) {
		case INLA_CGENERIC_GRAPH:
		{
			int ii = -1, jj, ii_prev;
			n = (int) out[0];
			m = (int) out[1];
			fprintf(fp, "\tdimension = %1d   number.of.elements = %1d\n", n, m);
			for (i = 0; i < m; i++) {
				ii_prev = ii;
				ii = (int) out[2 + i];
				jj = (int) out[2 + m + i];
				assert(ii >= 0 && jj >= 0);
				assert(ii >= ii_prev);
				assert(jj >= ii);
				fprintf(fp, "\tidx = %1d i = %1d j = %1d\n", i, (int) out[2 + i], (int) out[2 + m + i]);
			}
		}
			break;

		case INLA_CGENERIC_Q:
		{
			n = (int) out[0];
			m = (int) out[1];
			fprintf(fp, "\tn = %d  m = %d\n", n, m);
			if (n == -1) {
				fprintf(fp, "\toptimized format (same as the graph)\n");
				for (i = 0; i < m; i++) {
					fprintf(fp, "\tidx = %1d Q = %.8f\n", i, out[2 + i]);
				}
			} else {
				for (i = 0; i < m; i++) {
					fprintf(fp, "\tidx = %1d i = %1d, j = %1d, Qij = %.8f\n", i, (int) out[2 + i], (int) out[2 + n + i],
						out[2 + 2 * n + i]);
				}
			}
		}
			break;

		case INLA_CGENERIC_MU:
		{
			n = (int) out[0];
			fprintf(fp, "n = %1d\n", n);
			for (i = 0; i < n; i++) {
				fprintf(fp, "\ti = %1d mu_i = %.8f\n", i, out[1 + i]);
			}
		}
			break;

		case INLA_CGENERIC_INITIAL:
		{
			n = (int) out[0];
			fprintf(fp, "n = %1d\n", n);
			for (i = 0; i < n; i++) {
				fprintf(fp, "\tidx = %1d initial = %.8f\n", i, out[1 + i]);
			}
		}
			break;

		case INLA_CGENERIC_LOG_NORM_CONST:
		{
			fprintf(fp, "\tlog.norm.const = %.8f\n", (out ? out[0] : NAN));
		}
			break;

		case INLA_CGENERIC_LOG_PRIOR:
		{
			if (out) {
				fprintf(fp, "\tlog.prior = %.8f\n", out[0]);
			}
		}
			break;

		default:
			break;
		}
	}

	return GMRFLib_SUCCESS;
}

int inla_add_copyof(inla_tp *mb)
{
	int i, k, kk, kkk, nf = mb->nf;
	const int debug = 0;
	char *msg = NULL;

	if (debug) {
		for (k = 0; k < nf; k++) {
			printf("k= %1d tag= %s of= %s\n", k, mb->f_tag[k], mb->f_of[k]);
		}
	}

	for (k = 0; k < nf; k++) {
		if (mb->f_id[k] == F_COPY) {
			if (debug) {
				printf("ffield %d is F_COPY\n", k);
			}

			kk = find_tag(mb, mb->f_of[k]);
			if (kk < 0 || k == kk) {
				GMRFLib_sprintf(&msg, "ffield %1d is F_COPY and a copy of %s which is not found", k, mb->f_of[k]);
				inla_error_general(msg);
				exit(EXIT_FAILURE);
			}
			if (mb->f_id[kk] == F_COPY && kk > k) {
				GMRFLib_sprintf(&msg, "ffield [%s] is a copy of a (later defined) F_COPY field [%s]; please swap",
						mb->f_tag[k], mb->f_tag[kk]);
				inla_error_general(msg);
				exit(EXIT_FAILURE);
			}

			if (mb->f_same_as[k]) {
				kkk = find_tag(mb, mb->f_same_as[k]);
				if (kkk < 0) {
					GMRFLib_sprintf(&msg, "ffield %1d is F_COPY but same.as=[%s] is not found", k, mb->f_same_as[k]);
					inla_error_general(msg);
					exit(EXIT_FAILURE);
				}
				if (mb->f_id[kkk] != F_COPY) {
					GMRFLib_sprintf(&msg,
							"ffield [%s] is a copy of [%s], but same.as=[%s] which is not F_COPY\n",
							mb->f_tag[k], mb->f_tag[kk], mb->f_same_as[k]);
					inla_error_general(msg);
					exit(EXIT_FAILURE);
				}
				if (kkk == k) {
					GMRFLib_sprintf(&msg,
							"ffield [%s] is a copy of [%s], but same.as=[%s] which is not allowed.\n",
							mb->f_tag[k], mb->f_tag[kk], mb->f_same_as[k]);
					inla_error_general(msg);
					exit(EXIT_FAILURE);
				}
				if (kkk > k) {
					GMRFLib_sprintf(&msg,
							"ffield [%s] is a copy of [%s], but same.as=[%s] which is after; please swap. %1d > %1d\n",
							mb->f_tag[k], mb->f_tag[kk], mb->f_same_as[k], kkk, k);
					inla_error_general(msg);
					exit(EXIT_FAILURE);
				}
			} else {
				kkk = k;
			}

			if (debug) {
				if (mb->f_same_as[k]) {
					printf("found name %s at ffield %1d [same.as %s = ffield %1d]\n", mb->f_of[k], kk, mb->f_same_as[k], kkk);
				} else {
					printf("found name %s at ffield %1d\n", mb->f_of[k], kk);
				}
			}

			/*
			 * this is required! 
			 */
			if (!mb->ff_Qfunc) {
				mb->ff_Qfunc = Calloc(nf, GMRFLib_Qfunc_tp **);
				mb->ff_Qfunc_arg = Calloc(nf, void **);
				for (i = 0; i < nf; i++) {
					mb->ff_Qfunc[i] = Calloc(nf, GMRFLib_Qfunc_tp *);
					mb->ff_Qfunc_arg[i] = Calloc(nf, void *);
				}
			}

			/*
			 * yes, just use that size 
			 */
			GMRFLib_graph_mk_linear(&(mb->f_graph[k]), mb->f_Ntotal[kk], 0, 0);
			GMRFLib_free_constr(mb->f_constr[k]);  /* if its any */
			mb->f_constr[k] = NULL;
			mb->f_sumzero[k] = 0;
			mb->f_rankdef[k] = 0;

			inla_copy_arg_tp *arg = Calloc(1, inla_copy_arg_tp);

			arg->Qfunc = mb->f_Qfunc[kk];
			arg->Qfunc_arg = mb->f_Qfunc_arg[kk];
			arg->precision = mb->f_precision[k];
			arg->beta = mb->f_theta[kkk][0];

			arg->map_beta = mb->f_theta_map[kkk][0];
			arg->map_beta_arg = mb->f_theta_map_arg[kkk][0];

			if (0) {
				if (arg->map_beta_arg) {
					printf("range %g %g\n", ((double *) (arg->map_beta_arg))[0], ((double *) (arg->map_beta_arg))[1]);
				}
			}

			/*
			 * zero this out if its not needed anymore 
			 */
			if (k != kkk) {
				mb->f_theta[k] = NULL;
			}

			mb->f_Qfunc[kk] = Qfunc_copy_part00;
			mb->f_Qfunc_arg[kk] = (void *) arg;

			mb->f_Qfunc[k] = Qfunc_copy_part11;
			mb->f_Qfunc_arg[k] = (void *) arg;

			mb->ff_Qfunc[k][kk] = mb->ff_Qfunc[kk][k] = Qfunc_copy_part01;
			mb->ff_Qfunc_arg[k][kk] = mb->ff_Qfunc_arg[kk][k] = (void *) arg;

			mb->f_n[k] = mb->f_n[kk];
			mb->f_N[k] = mb->f_N[kk];
			mb->f_Ntotal[k] = mb->f_Ntotal[kk];
		}
	}
	return 0;
}

int inla_add_scopyof(inla_tp *mb)
{
	int i, k, kk, nf = mb->nf;
	const int debug = 0;
	char *msg = NULL;

	if (debug) {
		for (k = 0; k < nf; k++) {
			printf("k= %1d tag= %s of= %s\n", k, mb->f_tag[k], mb->f_of[k]);
		}
	}

	for (k = 0; k < nf; k++) {
		if (mb->f_id[k] == F_SCOPY) {
			if (debug) {
				printf("ffield %d is F_SCOPY\n", k);
			}

			kk = find_tag(mb, mb->f_of[k]);
			if (kk < 0 || k == kk) {
				GMRFLib_sprintf(&msg, "ffield %1d is F_SCOPY and a copy of %s which is not found", k, mb->f_of[k]);
				inla_error_general(msg);
				exit(EXIT_FAILURE);
			}
			if (mb->f_id[kk] == F_SCOPY && kk > k) {
				GMRFLib_sprintf(&msg, "ffield [%s] is a copy of a (later defined) F_SCOPY field [%s]; please swap",
						mb->f_tag[k], mb->f_tag[kk]);
				inla_error_general(msg);
				exit(EXIT_FAILURE);
			}

			if (debug) {
				printf("found name %s at ffield %1d\n", mb->f_of[k], kk);
			}

			/*
			 * this is required! 
			 */
			if (!mb->ff_Qfunc) {
				mb->ff_Qfunc = Calloc(nf, GMRFLib_Qfunc_tp **);
				mb->ff_Qfunc_arg = Calloc(nf, void **);
				for (i = 0; i < nf; i++) {
					mb->ff_Qfunc[i] = Calloc(nf, GMRFLib_Qfunc_tp *);
					mb->ff_Qfunc_arg[i] = Calloc(nf, void *);
				}
			}

			/*
			 * yes, just use that size 
			 */
			GMRFLib_graph_mk_linear(&(mb->f_graph[k]), mb->f_Ntotal[kk], 0, 0);
			GMRFLib_free_constr(mb->f_constr[k]);  /* if its any */
			mb->f_constr[k] = NULL;
			mb->f_sumzero[k] = 0;
			mb->f_rankdef[k] = 0;

			// yes, this is where its stored
			inla_scopy_arg_tp *arg = (inla_scopy_arg_tp *) mb->f_Qfunc_arg_orig[k];

			arg->Qfunc = mb->f_Qfunc[kk];
			arg->Qfunc_arg = mb->f_Qfunc_arg[kk];
			arg->precision = mb->f_precision[k];

			mb->f_Qfunc[kk] = Qfunc_scopy_part00;
			mb->f_Qfunc_arg[kk] = (void *) arg;
			mb->f_Qfunc[k] = Qfunc_scopy_part11;
			mb->f_Qfunc_arg[k] = (void *) arg;
			mb->ff_Qfunc[k][kk] = mb->ff_Qfunc[kk][k] = Qfunc_scopy_part01;
			mb->ff_Qfunc_arg[k][kk] = mb->ff_Qfunc_arg[kk][k] = (void *) arg;

			mb->f_n[k] = mb->f_n[kk];
			mb->f_N[k] = mb->f_N[kk];
			mb->f_Ntotal[k] = mb->f_Ntotal[kk];
		}
	}
	return 0;
}

inla_iarray_tp *find_all_f(inla_tp *mb, inla_component_tp id)
{
	inla_iarray_tp *ia = Calloc(1, inla_iarray_tp);

	ia->n = count_f(mb, id);
	if (ia->n) {
		int i, j;

		ia->array = Calloc(ia->n, int);
		for (i = j = 0; i < mb->nf; i++) {
			if (mb->f_id[i] == id) {
				ia->array[j++] = i;
			}
		}
		assert(j == ia->n);
	}

	return ia;
}

int find_f(inla_tp *mb, inla_component_tp id)
{
	for (int i = 0; i < mb->nf; i++) {
		if (mb->f_id[i] == id) {
			return i;
		}
	}
	return -1;
}

int find_tag(inla_tp *mb, const char *name)
{
	for (int i = 0; i < mb->nf; i++) {
		if (!strcasecmp((const char *) mb->f_tag[i], name))
			return i;
	}
	return -1;
}

int count_f(inla_tp *mb, inla_component_tp id)
{
	int n = 0;
	for (int i = 0; i < mb->nf; i++) {
		if (mb->f_id[i] == id) {
			n++;
		}
	}
	return n;
}

int inla_setup_ai_par_default(inla_tp *mb)
{
	/*
	 * change some these values to provide a inla-spesific defaults:
	 * 
	 * - the verbose controls the output - use CCD as default integrator, except for ntheta=1, where the GRID is used. 
	 */
	int i;

	if (!mb->ai_par) {
		GMRFLib_default_ai_param(&(mb->ai_par));

		if (!(G.mode == INLA_MODE_HYPER)) {
			/*
			 * default mode 
			 */

			if (mb->verbose) {
				mb->ai_par->fp_log = stdout;
			} else {
				mb->ai_par->fp_log = NULL;
			}
			if (mb->ntheta == 1) {
				mb->ai_par->int_strategy = GMRFLib_AI_INT_STRATEGY_GRID;
			} else {
				mb->ai_par->int_strategy = GMRFLib_AI_INT_STRATEGY_CCD;
			}

			for (i = 0; i < mb->nds; i++) {
				if (mb->data_sections[i].data_id == L_T) {
					/*
					 * use special options for the additive student-t 
					 */
					mb->ai_par->strategy = GMRFLib_AI_STRATEGY_FIT_SCGAUSSIAN;
					mb->ai_par->linear_correction = GMRFLib_AI_LINEAR_CORRECTION_FAST;
				}
			}
		} else {
			/*
			 * hyperparameter mode: special options 
			 */

			if (mb->verbose) {
				mb->ai_par->fp_log = stdout;
			} else {
				mb->ai_par->fp_log = NULL;
			}
			for (i = 0; i < mb->nds; i++) {
				if (mb->data_sections[i].data_id == L_T) {
					/*
					 * use special options for the additive student-t 
					 */
					mb->ai_par->strategy = GMRFLib_AI_STRATEGY_FIT_SCGAUSSIAN;
					mb->ai_par->linear_correction = GMRFLib_AI_LINEAR_CORRECTION_FAST;
				}
			}
			mb->ai_par->int_strategy = GMRFLib_AI_INT_STRATEGY_GRID;
			// mb->ai_par->skip_configurations = GMRFLib_FALSE;
			mb->ai_par->hessian_force_diagonal = GMRFLib_TRUE;
			switch (mb->ntheta) {
			case 0:
			case 1:
			{
				mb->ai_par->dz = 0.75;
				mb->ai_par->diff_log_dens = 6;
			}
				break;

			default:
				mb->ai_par->dz = 1.0;
				mb->ai_par->diff_log_dens = 5;
			}
			mb->ai_par->compute_nparam_eff = GMRFLib_FALSE;
		}
	}

	return INLA_OK;
}

double inla_ar1_cyclic_logdet(int N_orig, double phi)
{
	double logdet = 0.0;
	double tpon = 2.0 * M_PI / N_orig;

	for (int jj = 0; jj < N_orig; jj++) {
		logdet += log1p(SQR(phi) - phi * (cos(tpon * jj) + cos(tpon * (N_orig - 1.0) * jj)));
	}
	return (logdet);
}

double extra(int thread_id, double *theta, int ntheta, void *argument)
{
	int i, j, count = 0, nfixed = 0, fail, fixed0, fixed1, fixed2, fixed3, evaluate_hyper_prior = 1;
	const int debug = 0;

	double val = 0.0, log_precision, log_precision0, log_precision1, rho, rho_intern, beta, beta_intern, logit_rho,
	    group_rho = NAN, group_rho_intern = NAN, ngroup = NAN, normc_g = 0.0, n_orig = NAN, N_orig = NAN, rankdef_orig = NAN,
	    h2_intern, phi, phi_intern, a_intern, dof_intern, logdet, group_prec = NAN, group_prec_intern = NAN, grankdef = 0.0,
	    gcorr = 1.0, log_halflife, log_shape, alpha, gama, alpha1, alpha2;

	inla_tp *mb = NULL;
	gsl_matrix *Q = NULL;
	gsl_matrix *L = NULL;

#define _NOT_FIXED(_fx) (!mb->mode_fixed && !mb->_fx)
#define _SET_GROUP_RHO(_nt_)						\
	if (mb->f_group_model[i] != G_AR) {				\
		if (mb->f_ngroup[i] == 1) {				\
			assert(mb->f_ntheta[i] == (_nt_));		\
		} else {						\
			assert(mb->f_ntheta[i] == (_nt_)+1);		\
		}							\
	}								\
	ngroup = mb->f_ngroup[i];					\
	n_orig = mb->f_n[i]/ngroup;					\
	N_orig = mb->f_N[i]/ngroup;					\
	rankdef_orig = mb->f_rankdef[i]/ngroup;				\
	grankdef = 0.0;							\
	if (mb->f_ngroup[i] > 1) {					\
		if (mb->f_group_model[i] == G_AR) {			\
			double _log_precision, *_pacf = NULL,  *_pacf_intern = NULL; \
			int _p;						\
									\
			_p = mb->f_group_order[i];			\
			if (_NOT_FIXED(f_fixed[i][(_nt_) + 0])) {	\
				_log_precision = theta[count];		\
				count++;				\
			} else {					\
				_log_precision = mb->f_theta[i][(_nt_) + 0][thread_id][0]; \
			}						\
			_pacf = Calloc(_p, double);			\
			_pacf_intern = Calloc(_p, double);		\
			for (j = 0; j < _p; j++) {			\
				if (_NOT_FIXED(f_fixed[i][(_nt_) + j + 1])) { \
					_pacf_intern[j] = theta[count];	\
					count++;			\
				} else {				\
					_pacf_intern[j] = mb->f_theta[i][(_nt_) + j + 1][thread_id][0]; \
				}					\
				_pacf[j] = ar_map_pacf(_pacf_intern[j], MAP_FORWARD, NULL); \
			}						\
			double _marginal_prec, _conditional_prec, *_marginal_Q = NULL,  *_param = NULL,  *_zero = NULL;  \
			_marginal_Q = Calloc(ISQR(_p), double);		\
			ar_marginal_distribution(_p, _pacf, &_marginal_prec, _marginal_Q); \
			_conditional_prec = exp(_log_precision) / _marginal_prec; \
			_param = Calloc(1 + _p + ISQR(_p), double);	\
			_zero = Calloc(_p, double);			\
			_param[0] = _p;					\
			for (j = 0; j < ISQR(_p); j++) {		\
				_param[1 + _p + j] = _marginal_Q[j] * exp(_log_precision); \
			}						\
			normc_g = priorfunc_mvnorm(_zero, _param) + (ngroup - _p) * (-0.5 * log(2 * M_PI) + 0.5 * log(_conditional_prec)); \
			normc_g -= LOG_NORMC_GAUSSIAN * ngroup; /* This term goes into the main code therefore its removed here	*/ \
			normc_g *= (N_orig - rankdef_orig);		\
			if (_NOT_FIXED(f_fixed[i][(_nt_) + 0])) {	\
				val += PRIOR_EVAL(mb->f_prior[i][(_nt_) + 0], &_log_precision); \
			}						\
			for (j = 0; j < _p; j++) {			\
				if (_NOT_FIXED(f_fixed[i][(_nt_) + j + 1])) { \
					val += PRIOR_EVAL(mb->f_prior[i][(_nt_) + j + 1], &(_pacf_intern[j])); \
				}					\
			}						\
			Free(_param);					\
			Free(_zero);					\
			Free(_pacf);					\
			Free(_pacf_intern);				\
			Free(_marginal_Q);				\
		} else {						\
			if (_NOT_FIXED(f_fixed[i][(_nt_)])){		\
				group_rho_intern = group_prec_intern = theta[count]; \
				count++;				\
			} else {					\
				group_rho_intern = group_prec_intern= mb->f_theta[i][(_nt_)][thread_id][0]; \
			}						\
			if (mb->f_group_model[i] == G_EXCHANGEABLE) {	\
				int ingroup = (int) ngroup;		\
				group_rho = map_group_rho(group_rho_intern, MAP_FORWARD, (void *) &ingroup); \
				normc_g = - 0.5 * (log1p((ngroup - 1.0) * group_rho) + (ngroup-1)*LOG_1mp(group_rho)); \
				if (_NOT_FIXED(f_fixed[i][(_nt_)])) {	\
					val += PRIOR_EVAL(mb->f_prior[i][(_nt_)], &group_rho_intern); \
				}					\
			}else if (mb->f_group_model[i] == G_EXCHANGEABLE_POS) {	\
				int ingroup = (int) ngroup;		\
				group_rho = map_probability(group_rho_intern, MAP_FORWARD, (void *) &ingroup); \
				normc_g = - 0.5 * (log1p((ngroup - 1.0) * group_rho) + (ngroup-1)*LOG_1mp(group_rho)); \
				if (_NOT_FIXED(f_fixed[i][(_nt_)])) {	\
					val += PRIOR_EVAL(mb->f_prior[i][(_nt_)], &group_rho_intern); \
				}					\
			} else if (mb->f_group_model[i] == G_AR1) {	\
				group_rho = map_rho(group_rho_intern, MAP_FORWARD, NULL); \
				if (mb->f_group_cyclic[i]) {		\
					normc_g = - (ngroup - 0.0) * 0.5 * LOG_1mp(SQR(group_rho)) + 0.5 * inla_ar1_cyclic_logdet(ngroup, group_rho); \
				} else {				\
					normc_g = - (ngroup - 1.0) * 0.5 * LOG_1mp(SQR(group_rho)); \
				}					\
				if (_NOT_FIXED(f_fixed[i][(_nt_)])) {	\
					val += PRIOR_EVAL(mb->f_prior[i][(_nt_)], &group_rho_intern); \
				}					\
			} else if (mb->f_group_model[i] == G_RW1 || mb->f_group_model[i] == G_RW2 || mb->f_group_model[i] == G_BESAG) { \
				grankdef = (mb->f_group_model[i] == G_RW1 || mb->f_group_model[i] == G_BESAG ? 1.0 : 2.0); \
				group_prec = map_precision(group_prec_intern, MAP_FORWARD, NULL); \
				normc_g = 0.5 * (ngroup - grankdef) * log(group_prec); \
				if (_NOT_FIXED(f_fixed[i][(_nt_)])) {	\
					val += PRIOR_EVAL(mb->f_prior[i][(_nt_)], &group_prec_intern); \
				}					\
			} else if (mb->f_group_model[i] == G_IID) {	\
				grankdef = 0.0;				\
				group_prec = map_precision(group_prec_intern, MAP_FORWARD, NULL); \
				normc_g = 0.5 * (ngroup - grankdef) * log(group_prec); \
				if (_NOT_FIXED(f_fixed[i][(_nt_)])) {	\
					val += PRIOR_EVAL(mb->f_prior[i][(_nt_)], &group_prec_intern); \
				}					\
			} else {					\
				abort();				\
			}						\
			normc_g *= (N_orig - rankdef_orig);		\
		}							\
		gcorr = 1.0 - grankdef / mb->f_ngroup[i];		\
	} else {							\
		group_rho = group_rho_intern = 0.0;			\
		normc_g = 0.0;						\
		ngroup = 1.0;						\
		n_orig = mb->f_n[i];					\
		N_orig = mb->f_N[i];					\
		rankdef_orig = mb->f_rankdef[i];			\
		grankdef = 0.0;						\
		gcorr = 1.0;						\
	}

	mb = (inla_tp *) argument;

	/*
	 * this will evaluate all the hyperparameters and disable EVAL_PRIOR...
	 */
	if (mb->update) {
		val += inla_update_density(theta, mb->update);
		evaluate_hyper_prior = 0;
	}
	// joint prior evaluated in R
	static int jp_first_time = 1;
	static void *jp_vec_sexp = NULL;

	if (mb->jp) {
		{
			if (jp_first_time) {
#pragma omp critical (Name_886e56613636db7114635f9aacba5e4fc01cca40)
				if (jp_first_time) {
					char **vec_str = NULL;

					inla_R_load(mb->jp->file);
					if (mb->ntheta > 0) {
						vec_str = Calloc(mb->ntheta, char *);
						for (i = 0; i < mb->ntheta; i++) {
							vec_str[i] = Strdup(mb->theta_tag[i]);
						}
					}
					jp_vec_sexp = inla_R_vector_of_strings(mb->ntheta, vec_str);
					if (vec_str) {
						for (i = 0; i < mb->ntheta; i++) {
							Free(vec_str[i]);
						}
						Free(vec_str);
					}
				}
				jp_first_time = 0;
			}
			assert(!(mb->update));		       /* only one at the time... */
			evaluate_hyper_prior = 0;

			int verbose = 0;
			if (ntheta > 0) {
				int n_out = 0;
				double *lprior = NULL;

				inla_R_funcall_jp(&n_out, &lprior, (const char *) mb->jp->model, &ntheta, theta, jp_vec_sexp);
				assert(n_out == 1);
				assert(lprior);
				val += *lprior;
				if (verbose) {
					printf("got lprior = %g\n", *lprior);
				}
				Free(lprior);
			} else {
				if (verbose) {
					printf("lprior: ntheta= 0, no prior needed\n");
				}
			}
		}
	}

	/*
	 * this is for the linear predictor
	 */

	if (GMRFLib_inla_mode == GMRFLib_MODE_CLASSIC) {
		if (!mb->predictor_fixed) {
			log_precision = theta[count];
			count++;
		} else {
			log_precision = mb->predictor_log_prec[thread_id][0];
		}
		val += mb->predictor_n * (LOG_NORMC_GAUSSIAN + 1.0 / 2.0 * log_precision);
		if (!mb->predictor_fixed) {
			val += PRIOR_EVAL(mb->predictor_prior, &log_precision);
		}

		/*
		 * this is for the A-matrix
		 */
		if (mb->predictor_m > 0) {
			log_precision = log(mb->predictor_Aext_precision);
			val += mb->predictor_m * (LOG_NORMC_GAUSSIAN + 1.0 / 2.0 * log_precision);
		}
	}

	if (mb->data_ntheta_all) {
		int check = 0;

		for (j = 0; j < mb->nds; j++) {
			Data_section_tp *ds = &(mb->data_sections[j]);

			check += ds->data_ntheta;
			switch (ds->data_id) {

			case L_SEM:
				break;

			case L_GAUSSIAN:
			{
				if (!ds->data_fixed0) {
					/*
					 * we only need to add the prior, since the normalisation constant due to the likelihood, is included in the likelihood
					 * function.
					 */
					log_precision = theta[count];
					val += PRIOR_EVAL(ds->data_prior0, &log_precision);
					count++;
				}
				if (!ds->data_fixed1) {
					// this should not happen, as data_fixed1 must be TRUE. We prepare the code for it, in
					// any case. 
					log_precision = theta[count];
					val += PRIOR_EVAL(ds->data_prior1, &log_precision);
					count++;
				}
			}
				break;

			case L_GAUSSIANJW:
			{
				for (int k = 0; k < 3; k++) {
					if (!ds->data_nfixed[k]) {
						beta = theta[count];
						val += PRIOR_EVAL(ds->data_nprior[k], &beta);
						count++;
					}
				}
			}
				break;

			case L_AGAUSSIAN:
			{
				if (!ds->data_fixed0) {
					log_precision = theta[count];
					val += PRIOR_EVAL(ds->data_prior, &log_precision);
					count++;
				}
			}
				break;

			case L_BC_GAUSSIAN:
			{
				if (!ds->data_fixed0) {
					log_precision = theta[count];
					val += PRIOR_EVAL(ds->data_prior0, &log_precision);
					count++;
				}
				if (!ds->data_fixed1) {
					double lambda = theta[count];
					val += PRIOR_EVAL(ds->data_prior1, &lambda);
					count++;
				}
			}
				break;

			case L_LOGNORMAL:
			{
				if (!ds->data_fixed) {
					log_precision = theta[count];
					val += PRIOR_EVAL(ds->data_prior, &log_precision);
					count++;
				}
			}
				break;

			case L_LOGNORMALSURV:
			{
				if (!ds->data_nfixed[0]) {
					log_precision = theta[count];
					val += PRIOR_EVAL(ds->data_nprior[0], &log_precision);
					count++;
				}
				for (int k = 0; k < ds->data_observations.cure_ncov; k++) {
					if (!ds->data_nfixed[1 + k]) {
						beta = theta[count];
						val += PRIOR_EVAL(ds->data_nprior[1 + k], &beta);
						count++;
					}
				}
			}
				break;

			case L_EXPPOWER:
			{
				if (!ds->data_fixed0) {
					log_precision = theta[count];
					val += PRIOR_EVAL(ds->data_prior0, &log_precision);
					count++;
				}
				if (!ds->data_fixed1) {
					double log_power = theta[count];
					val += PRIOR_EVAL(ds->data_prior1, &log_power);
					count++;
				}
			}
				break;

			case L_SIMPLEX:
			{
				if (!ds->data_fixed) {
					log_precision = theta[count];
					val += PRIOR_EVAL(ds->data_prior, &log_precision);
					count++;
				}
			}
				break;

			case L_GPOISSON:
			{
				if (!ds->data_fixed0) {
					double log_overdispersion = theta[count];
					val += PRIOR_EVAL(ds->data_prior0, &log_overdispersion);
					count++;
				}
				if (!ds->data_fixed1) {
					double p = theta[count];
					val += PRIOR_EVAL(ds->data_prior1, &p);
					count++;
				}
			}
				break;

			case L_CIRCULAR_NORMAL:
			{
				if (!ds->data_fixed) {
					log_precision = theta[count];
					val += PRIOR_EVAL(ds->data_prior, &log_precision);
					count++;
				}
			}
				break;

			case L_WRAPPED_CAUCHY:
			{
				if (!ds->data_fixed) {
					log_precision = theta[count];
					val += PRIOR_EVAL(ds->data_prior, &log_precision);
					count++;
				}
			}
				break;

			case L_TWEEDIE:
			{
				if (!ds->data_fixed0) {
					double p_intern = theta[count];
					val += PRIOR_EVAL(ds->data_prior0, &p_intern);
					count++;
				}
				if (!ds->data_fixed1) {
					double log_phi = theta[count];
					val += PRIOR_EVAL(ds->data_prior1, &log_phi);
					count++;
				}
			}
				break;

			case L_FMRI:
			case L_FMRISURV:
			{
				if (!ds->data_fixed0) {
					double lprec = theta[count];
					val += PRIOR_EVAL(ds->data_prior0, &lprec);
					count++;
				}
				if (!ds->data_fixed1) {
					double ldof = theta[count];
					val += PRIOR_EVAL(ds->data_prior1, &ldof);
					count++;
				}
			}
				break;

			case L_GP:
			case L_DGP:
			{
				if (!ds->data_fixed) {
					double log_tail = theta[count];
					val += PRIOR_EVAL(ds->data_prior, &log_tail);
					count++;
				}
			}
				break;

			case L_EGP:
			{
				if (!ds->data_fixed0) {
					double internal_tail = theta[count];
					val += PRIOR_EVAL(ds->data_prior0, &internal_tail);
					count++;
				}
				if (!ds->data_fixed1) {
					double internal_shape = theta[count];
					val += PRIOR_EVAL(ds->data_prior1, &internal_shape);
					count++;
				}
			}
				break;

			case L_IID_GAMMA:
			{
				if (!ds->data_fixed0) {
					log_shape = theta[count];
					val += PRIOR_EVAL(ds->data_prior0, &log_shape);
					count++;
				}
				if (!ds->data_fixed1) {
					double log_rate = theta[count];
					val += PRIOR_EVAL(ds->data_prior1, &log_rate);
					count++;
				}
			}
				break;

			case L_IID_LOGITBETA:
			{
				if (!ds->data_fixed0) {
					double log_a = theta[count];
					val += PRIOR_EVAL(ds->data_prior0, &log_a);
					count++;
				}
				if (!ds->data_fixed1) {
					double log_b = theta[count];
					val += PRIOR_EVAL(ds->data_prior1, &log_b);
					count++;
				}
			}
				break;

			case L_LOGGAMMA_FRAILTY:
			{
				if (!ds->data_fixed) {
					log_precision = theta[count];
					val += PRIOR_EVAL(ds->data_prior, &log_precision);
					count++;
				}
			}
				break;

			case L_LOGISTIC:
			{
				if (!ds->data_fixed) {
					log_precision = theta[count];
					val += PRIOR_EVAL(ds->data_prior, &log_precision);
					count++;
				}
			}
				break;

			case L_SKEWNORMAL:
			{
				if (!ds->data_fixed0) {
					log_precision = theta[count];

					val += PRIOR_EVAL(ds->data_prior0, &log_precision);
					count++;
				}
				if (!ds->data_fixed1) {
					double skewness = theta[count];

					val += PRIOR_EVAL(ds->data_prior1, &skewness);
					count++;
				}
			}
				break;

			case L_GEV:
			{
				if (!ds->data_fixed0) {
					log_precision = theta[count];
					val += PRIOR_EVAL(ds->data_prior0, &log_precision);
					count++;
				}
				if (!ds->data_fixed1) {
					/*
					 * this is the gev-parameter. Note that we need to scale it back to the real scale, as
					 * 'scale_xi' is there to help the numerics only
					 */
					double xi = theta[count] * ds->data_observations.gev_scale_xi;
					val += PRIOR_EVAL(ds->data_prior1, &xi) + log(ds->data_observations.gev_scale_xi);
					count++;
				}
			}
				break;

			case L_BGEV:
			{
				{
					if (!ds->data_nfixed[0]) {
						double spread = theta[count];
						val += PRIOR_EVAL(ds->data_nprior[0], &spread);
						count++;
					}
					if (!ds->data_nfixed[1]) {
						double intern_tail = theta[count];
						val += PRIOR_EVAL(ds->data_nprior[1], &intern_tail);
						count++;
					}

					int nbetas = ds->data_observations.bgev_nbetas[0] + ds->data_observations.bgev_nbetas[1];
					int off = 2;
					for (int k = off; k < off + nbetas; k++) {
						if (!ds->data_nfixed[k]) {
							double b = theta[count];
							if (k < ds->data_observations.bgev_nbetas[0]) {
								val += PRIOR_EVAL(ds->data_nprior[k], &b);
							} else {
								val += PRIOR_EVAL(ds->data_nprior[k], &b);
							}
							count++;
						}
					}
				}
			}
				break;

			case L_GGAUSSIAN:
			case L_GGAUSSIANS:
			{
				int nbeta = ds->data_observations.ggaussian_nbeta;
				for (int k = 0; k < nbeta; k++) {
					if (!ds->data_nfixed[k]) {
						double b = theta[count];
						val += PRIOR_EVAL(ds->data_nprior[k], &b);
						count++;
					}
				}
			}
				break;

			case L_RCPOISSON:
			{
				int nbeta = ds->data_observations.rcp_nbeta;
				for (int k = 0; k < nbeta; k++) {
					if (!ds->data_nfixed[k]) {
						double b = theta[count];
						val += PRIOR_EVAL(ds->data_nprior[k], &b);
						count++;
					}
				}
			}
				break;

			case L_TPOISSON:
			{
				int nbeta = ds->data_observations.tp_nbeta;
				for (int k = 0; k < nbeta; k++) {
					if (!ds->data_nfixed[k]) {
						double b = theta[count];
						val += PRIOR_EVAL(ds->data_nprior[k], &b);
						count++;
					}
				}
			}
				break;

			case L_0POISSON:
			case L_0POISSONS:
			{
				int nbeta = ds->data_observations.poisson0_nbeta;
				for (int k = 0; k < nbeta; k++) {
					if (!ds->data_nfixed[k]) {
						double b = theta[count];
						val += PRIOR_EVAL(ds->data_nprior[k], &b);
						count++;
					}
				}
			}
				break;

			case L_0BINOMIAL:
			case L_0BINOMIALS:
			{
				int nbeta = ds->data_observations.binomial0_nbeta;
				for (int k = 0; k < nbeta; k++) {
					if (!ds->data_nfixed[k]) {
						double b = theta[count];
						val += PRIOR_EVAL(ds->data_nprior[k], &b);
						count++;
					}
				}
			}
				break;

			case L_BINOMIALMIX:
			{
				for (int k = 0; k < BINOMIALMIX_NBETA; k++) {
					if (!ds->data_nfixed[k]) {
						double b = theta[count];
						val += PRIOR_EVAL(ds->data_nprior[k], &b);
						count++;
					}
				}
			}
				break;

			case L_GAMMA:
			case L_MGAMMA:
			{
				if (!ds->data_fixed) {
					double precision_intern = theta[count];

					val += PRIOR_EVAL(ds->data_prior, &precision_intern);
					count++;
				}
			}
				break;

			case L_GAMMASURV:
			case L_MGAMMASURV:
			{
				if (!ds->data_nfixed[0]) {
					double precision_intern = theta[count];
					val += PRIOR_EVAL(ds->data_nprior[0], &precision_intern);
					count++;
				}
				for (int k = 0; k < ds->data_observations.cure_ncov; k++) {
					if (!ds->data_nfixed[1 + k]) {
						beta = theta[count];
						val += PRIOR_EVAL(ds->data_nprior[1 + k], &beta);
						count++;
					}
				}
			}
				break;

			case L_GAMMAJW:
				break;

			case L_GAMMACOUNT:
			{
				if (!ds->data_fixed) {
					double log_alpha = theta[count];
					val += PRIOR_EVAL(ds->data_prior, &log_alpha);
					count++;
				}
			}
				break;

			case L_QKUMAR:
			{
				if (!ds->data_fixed) {
					double precision_intern = theta[count];
					val += PRIOR_EVAL(ds->data_prior, &precision_intern);
					count++;
				}
			}
				break;

			case L_LOGLOGISTIC:
			case L_QLOGLOGISTIC:
			{
				if (!ds->data_fixed) {
					double alpha_intern = theta[count];
					val += PRIOR_EVAL(ds->data_prior, &alpha_intern);
					count++;
				}
			}
				break;

			case L_LOGLOGISTICSURV:
			case L_QLOGLOGISTICSURV:
			{
				if (!ds->data_nfixed[0]) {
					double alpha_intern = theta[count];
					val += PRIOR_EVAL(ds->data_nprior[0], &alpha_intern);
					count++;
				}
				for (int k = 0; k < ds->data_observations.cure_ncov; k++) {
					if (!ds->data_nfixed[1 + k]) {
						beta = theta[count];
						val += PRIOR_EVAL(ds->data_nprior[1 + k], &beta);
						count++;
					}
				}
			}
				break;

			case L_BETA:
			{
				if (!ds->data_fixed) {
					double precision_intern = theta[count];
					val += PRIOR_EVAL(ds->data_prior, &precision_intern);
					count++;
				}
			}
				break;

			case L_BETABINOMIAL:
			case L_BETABINOMIALNA:
			{
				if (!ds->data_fixed) {
					double intern_overdispersion = theta[count];
					val += PRIOR_EVAL(ds->data_prior, &intern_overdispersion);
					count++;
				}
			}
				break;

			case L_NBINOMIAL:
			case L_CENNBINOMIAL2:
			{
				if (!ds->data_fixed) {
					double log_size = theta[count];
					val += PRIOR_EVAL(ds->data_prior, &log_size);
					count++;
				}
			}
				break;

			case L_ZEROINFLATEDNBINOMIAL0:
			case L_ZEROINFLATEDNBINOMIAL1:
			{
				if (!ds->data_fixed0) {
					double log_size = theta[count];
					val += PRIOR_EVAL(ds->data_prior0, &log_size);
					count++;
				}
				if (!ds->data_fixed1) {
					/*
					 * this is the probability-parameter in the zero-inflated nbinomial_0/1
					 */
					double prob_intern = theta[count];
					val += PRIOR_EVAL(ds->data_prior1, &prob_intern);
					count++;
				}
			}
				break;

			case L_ZEROINFLATEDNBINOMIAL2:
			{
				if (!ds->data_fixed0) {
					double log_size = theta[count];
					val += PRIOR_EVAL(ds->data_prior0, &log_size);
					count++;
				}
				if (!ds->data_fixed1) {
					/*
					 * this is the alpha-parameter in the zero-inflated nbinomial_0/1
					 */
					double alpha_intern = theta[count];
					val += PRIOR_EVAL(ds->data_prior1, &alpha_intern);
					count++;
				}
			}
				break;

			case L_ZEROINFLATEDBETABINOMIAL0:
			case L_ZEROINFLATEDBETABINOMIAL1:
			{
				if (!ds->data_fixed0) {
					rho_intern = theta[count];
					val += PRIOR_EVAL(ds->data_prior0, &rho_intern);
					count++;
				}
				if (!ds->data_fixed1) {
					/*
					 * this is the probability-parameter in the zero-inflated betabinomial_0/1
					 */
					double prob_intern = theta[count];
					val += PRIOR_EVAL(ds->data_prior1, &prob_intern);
					count++;
				}
			}
				break;

			case L_ZEROINFLATEDNBINOMIAL1STRATA2:
			{
				if (!ds->data_fixed0) {
					double log_size = theta[count];
					val += PRIOR_EVAL(ds->data_prior0, &log_size);
					count++;
				}
				for (int icount = 0; icount < STRATA_MAXTHETA; icount++) {
					if (!ds->data_nfixed[icount]) {
						/*
						 * this is the probability-parameter in the zero-inflated nbinomial_strata2
						 */
						double prob_intern = theta[count];
						val += PRIOR_EVAL(ds->data_nprior[icount], &prob_intern);
						count++;
					}
				}
			}
				break;

			case L_ZEROINFLATEDNBINOMIAL1STRATA3:
			{
				if (!ds->data_fixed0) {
					double prob_intern = theta[count];
					val += PRIOR_EVAL(ds->data_prior0, &prob_intern);
					count++;
				}
				for (int icount = 0; icount < STRATA_MAXTHETA; icount++) {
					if (!ds->data_nfixed[icount]) {
						/*
						 * this is the size-parameter in the zero-inflated nbinomial_strata2
						 */
						double log_size = theta[count];
						val += PRIOR_EVAL(ds->data_nprior[icount], &log_size);
						count++;
					}
				}
			}
				break;

			case L_ZERO_N_INFLATEDBINOMIAL2:
			{
				if (!ds->data_fixed0) {
					double log_alpha1 = theta[count];
					val += PRIOR_EVAL(ds->data_prior0, &log_alpha1);
					count++;
				}
				if (!ds->data_fixed1) {
					double log_alpha2 = theta[count];
					val += PRIOR_EVAL(ds->data_prior1, &log_alpha2);
					count++;
				}
			}
				break;

			case L_ZERO_N_INFLATEDBINOMIAL3:
			{
				if (!ds->data_fixed0) {
					double log_alpha0 = theta[count];
					val += PRIOR_EVAL(ds->data_prior0, &log_alpha0);
					count++;
				}
				if (!ds->data_fixed1) {
					double log_alphaN = theta[count];
					val += PRIOR_EVAL(ds->data_prior1, &log_alphaN);
					count++;
				}
			}
				break;

			case L_T:
			{
				if (!ds->data_fixed0) {
					log_precision = theta[count];
					val += PRIOR_EVAL(ds->data_prior0, &log_precision);
					count++;
				}
				if (!ds->data_fixed1) {
					dof_intern = theta[count];
					val += PRIOR_EVAL(ds->data_prior1, &dof_intern);
					count++;
				}
			}
				break;

			case L_TSTRATA:
			{
				int k;
				for (k = 0; k < TSTRATA_MAXTHETA; k++) {
					if (!ds->data_nfixed[k]) {
						double th = theta[count];
						val += PRIOR_EVAL(ds->data_nprior[k], &th);
						count++;
					}
				}
			}
				break;

			case L_STOCHVOL:
			{
				if (!ds->data_fixed) {
					log_precision = theta[count];
					val += PRIOR_EVAL(ds->data_prior, &log_precision);
					count++;
				}
			}
				break;

			case L_STOCHVOL_LN:
			{
				if (!ds->data_fixed) {
					double off = theta[count];
					val += PRIOR_EVAL(ds->data_prior, &off);
					count++;
				}
			}
				break;

			case L_STOCHVOL_SN:
			{
				if (!ds->data_fixed0) {
					double skew_intern = theta[count];
					val += PRIOR_EVAL(ds->data_prior0, &skew_intern);
					count++;
				}
				if (!ds->data_fixed1) {
					double log_prec_offset = theta[count];
					val += PRIOR_EVAL(ds->data_prior1, &log_prec_offset);
					count++;
				}
			}
				break;

			case L_STOCHVOL_T:
			{
				if (!ds->data_fixed) {
					dof_intern = theta[count];
					val += PRIOR_EVAL(ds->data_prior, &dof_intern);
					count++;
				}
			}
				break;

			case L_STOCHVOL_NIG:
			{
				if (!ds->data_fixed0) {
					/*
					 * this is the skewness 
					 */
					double skew = theta[count];
					val += PRIOR_EVAL(ds->data_prior0, &skew);
					count++;
				}
				if (!ds->data_fixed1) {
					/*
					 * this is the shape 
					 */
					double shape_intern = theta[count];
					val += PRIOR_EVAL(ds->data_prior1, &shape_intern);
					count++;
				}
			}
				break;

			case L_WEIBULL:
			{
				if (!ds->data_fixed) {
					double alpha_intern = theta[count];
					val += PRIOR_EVAL(ds->data_prior, &alpha_intern);
					count++;
				}
			}
				break;

			case L_WEIBULLSURV:
			{
				if (!ds->data_nfixed[0]) {
					double alpha_intern = theta[count];
					val += PRIOR_EVAL(ds->data_nprior[0], &alpha_intern);
					count++;
				}
				for (int k = 0; k < ds->data_observations.cure_ncov; k++) {
					if (!ds->data_nfixed[1 + k]) {
						beta = theta[count];
						val += PRIOR_EVAL(ds->data_nprior[1 + k], &beta);
						count++;
					}
				}
			}
				break;

			case L_GOMPERTZ:
			{
				if (!ds->data_fixed) {
					double alpha_intern = theta[count];
					val += PRIOR_EVAL(ds->data_prior, &alpha_intern);
					count++;
				}
			}
				break;

			case L_GOMPERTZSURV:
			{
				if (!ds->data_nfixed[0]) {
					double alpha_intern = theta[count];
					val += PRIOR_EVAL(ds->data_nprior[0], &alpha_intern);
					count++;
				}
				for (int k = 0; k < ds->data_observations.cure_ncov; k++) {
					if (!ds->data_nfixed[1 + k]) {
						beta = theta[count];
						val += PRIOR_EVAL(ds->data_nprior[1 + k], &beta);
						count++;
					}
				}
			}
				break;

			case L_ZEROINFLATEDPOISSON0:
			case L_ZEROINFLATEDPOISSON1:
			case L_ZEROINFLATEDCENPOISSON0:
			case L_ZEROINFLATEDCENPOISSON1:
			case L_POISSON_SPECIAL1:
			{
				if (!ds->data_fixed) {
					/*
					 * this is the probability-parameter in the zero-inflated Poisson_0/1 or special1
					 */
					double prob_intern = theta[count];
					val += PRIOR_EVAL(ds->data_prior, &prob_intern);
					count++;
				}
			}
				break;

			case L_ZEROINFLATEDPOISSON2:
			{
				if (!ds->data_fixed) {
					/*
					 * this is the probability-parameter in the zero-inflated Poisson_2
					 */
					double alpha_intern = theta[count];
					val += PRIOR_EVAL(ds->data_prior, &alpha_intern);
					count++;
				}
			}
				break;

			case L_ZEROINFLATEDBINOMIAL2:
				if (!ds->data_fixed) {
					/*
					 * this is the probability-parameter in the zero-inflated Binomial_2
					 */
					double alpha_intern = theta[count];
					val += PRIOR_EVAL(ds->data_prior, &alpha_intern);
					count++;
				}
				break;

			case L_ZEROINFLATEDBINOMIAL0:
			case L_ZEROINFLATEDBINOMIAL1:
			{
				if (!ds->data_fixed) {
					/*
					 * this is the probability-parameter in the zero-inflated Binomial_0/1
					 */
					double prob_intern = theta[count];
					val += PRIOR_EVAL(ds->data_prior, &prob_intern);
					count++;
				}
			}
				break;

			case L_ZEROINFLATEDBETABINOMIAL2:
			{
				if (!ds->data_fixed0) {
					/*
					 * this is the probability-related-parameter 
					 */
					double alpha_intern = theta[count];
					val += PRIOR_EVAL(ds->data_prior0, &alpha_intern);
					count++;
				}
				if (!ds->data_fixed1) {
					/*
					 * this is the delta-parameter 
					 */
					double delta_intern = theta[count];
					val += PRIOR_EVAL(ds->data_prior1, &delta_intern);
					count++;
				}
			}
				break;

			case L_OCCUPANCY:
			{
				for (int k = 0; k < OCCUPANCY_MAXTHETA; k++) {
					if (!ds->data_nfixed[k]) {
						beta = theta[count];
						val += PRIOR_EVAL(ds->data_nprior[k], &beta);
						count++;
					}
				}
			}
				break;

			case L_NMIX:
			{
				for (int k = 0; k < NMIX_MMAX; k++) {
					if (!ds->data_nfixed[k]) {
						beta = theta[count];
						val += PRIOR_EVAL(ds->data_nprior[k], &beta);
						count++;
					}
				}
			}
				break;

			case L_NMIXNB:
			{
				/*
				 *  the last one here is the log_overdispersion, which I do not rename to, for simplicity
				 */
				for (int k = 0; k < NMIX_MMAX + 1; k++) {
					if (!ds->data_nfixed[k]) {
						beta = theta[count];
						val += PRIOR_EVAL(ds->data_nprior[k], &beta);
						count++;
					}
				}
			}
				break;

			case L_BELL:
				break;

			case L_STDGAUSSIAN:
			case L_BINOMIAL:
			case L_XBINOMIAL:
			case L_EXPONENTIAL:
			case L_POISSON:
			case L_NPOISSON:
			case L_NZPOISSON:
			case L_XPOISSON:
			case L_CONTPOISSON:
			case L_QCONTPOISSON:
				break;

			case L_EXPONENTIALSURV:
			case L_GAMMAJWSURV:
			{
				for (int k = 0; k < ds->data_observations.cure_ncov; k++) {
					if (!ds->data_nfixed[k]) {
						beta = theta[count];
						val += PRIOR_EVAL(ds->data_nprior[k], &beta);
						count++;
					}
				}
			}
				break;

			case L_POM:
			{
				double *v = Calloc(POM_MAXTHETA, double);
				int v_count = 0;
				for (int k = 0; k < POM_MAXTHETA; k++) {
					if (!ds->data_nfixed[k]) {
						v[v_count] = theta[count];
						v_count++;
						count++;
					}
				}
				if (v_count > 0) {
					val += PRIOR_EVAL(ds->data_nprior[0], v);
				}
				Free(v);
			}
				break;

			case L_FL:
				break;

			default:
				break;
			}

			/*
			 * link-models 
			 */
			check += ds->link_ntheta;

			if (debug) {
				P(ds->link_id);
				P(ds->link_ntheta);
			}

			switch (ds->link_id) {
			case LINK_IDENTITY:
			case LINK_INVERSE:
			case LINK_LOG:
			case LINK_LOGa:
			case LINK_NEGLOG:
			case LINK_PROBIT:
			case LINK_CLOGLOG:
			case LINK_CCLOGLOG:
			case LINK_LOGLOG:
			case LINK_CAUCHIT:
			case LINK_LOGIT:
			case LINK_TAN:
			case LINK_QPOISSON:
			case LINK_QBINOMIAL:
			case LINK_QWEIBULL:
			case LINK_QGAMMA:
			case LINK_QEXPPOWER:
				break;

			case LINK_LOGOFFSET:
			{
				if (!ds->link_fixed[0]) {
					beta_intern = theta[count];
					val += PRIOR_EVAL(ds->link_prior[0], &beta_intern);
					count++;
				}
			}
				break;

			case LINK_LOGITOFFSET:
			{
				if (!ds->link_fixed[0]) {
					double prob_intern = theta[count];
					val += PRIOR_EVAL(ds->link_prior[0], &prob_intern);
					count++;
				}
			}
				break;

			case LINK_SSLOGIT:
			{
				if (!ds->link_fixed[0]) {
					double sensitivity_intern = theta[count];
					val += PRIOR_EVAL(ds->link_prior[0], &sensitivity_intern);
					count++;
				}
				if (!ds->link_fixed[1]) {
					double specificity_intern = theta[count];
					val += PRIOR_EVAL(ds->link_prior[1], &specificity_intern);
					count++;
				}
			}
				break;

			case LINK_ROBIT:
			{
				if (!ds->link_fixed[0]) {
					dof_intern = theta[count];
					val += PRIOR_EVAL(ds->link_prior[0], &dof_intern);
					count++;
				}
			}
				break;

			case LINK_SN:
			{
				if (!ds->link_fixed[0]) {
					double skew = theta[count];
					val += PRIOR_EVAL(ds->link_prior[0], &skew);
					count++;
				}
				if (!ds->link_fixed[1]) {
					double intercept = theta[count];
					val += PRIOR_EVAL(ds->link_prior[1], &intercept);
					count++;
				}
			}
				break;

			case LINK_GEV:
			case LINK_CGEV:
			{
				if (!ds->link_fixed[0]) {
					double tail_intern = theta[count];
					val += PRIOR_EVAL(ds->link_prior[0], &tail_intern);
					count++;
				}
				if (!ds->link_fixed[1]) {
					double intercept_intern = theta[count];
					val += PRIOR_EVAL(ds->link_prior[1], &intercept_intern);
					count++;
				}
			}
				break;

			case LINK_POWER_LOGIT:
			{
				if (!ds->link_fixed[0]) {
					double power = theta[count];
					val += PRIOR_EVAL(ds->link_prior[0], &power);
					count++;
				}
				if (!ds->link_fixed[1]) {
					double intercept = theta[count];
					val += PRIOR_EVAL(ds->link_prior[1], &intercept);
					count++;
				}
			}
				break;

			case LINK_TEST1:
			{
				if (!ds->link_fixed[0]) {
					beta = theta[count];
					val += PRIOR_EVAL(ds->link_prior[0], &beta);
					count++;
				}
			}
				break;

			case LINK_SPECIAL2:
			{
				if (!ds->link_fixed[0]) {
					beta = theta[count];
					val += PRIOR_EVAL(ds->link_prior[0], &beta);
					count++;
				}
			}
				break;

			case LINK_SPECIAL1:
			{
				if (!ds->link_fixed[0]) {
					log_precision = theta[count];
					val += PRIOR_EVAL(ds->link_prior[0], &log_precision);
					count++;
				} else {
					// log_precision = ds->link_parameters->log_prec[thread_id][0]; 
				}

				if (ds->link_order > 0) {
					double *bbeta = Calloc(ds->link_order, double);
					for (j = 0; j < ds->link_order; j++) {
						if (!ds->link_fixed[j + 1]) {
							bbeta[j] = theta[count];
							count++;
						} else {
							bbeta[j] = ds->link_parameters->betas[j][thread_id][0];
						}
					}
					val += PRIOR_EVAL(ds->link_prior[1], bbeta);
					Free(bbeta);
				}
			}
				break;

			default:
				GMRFLib_ASSERT(0 == 1, GMRFLib_ESNH);
			}

			/*
			 * re-models 
			 */
			if (ds->mix_use) {
				check += ds->mix_ntheta;

				switch (ds->mix_id) {
				case MIX_GAUSSIAN:
				case MIX_LOGGAMMA:
				case MIX_MLOGGAMMA:
					if (!ds->mix_fixed) {
						log_precision = theta[count];
						val += PRIOR_EVAL(ds->mix_prior, &log_precision);
						count++;
					}
					break;
				default:
					GMRFLib_ASSERT(0 == 1, GMRFLib_ESNH);
				}
			}
		}
		assert(mb->data_ntheta_all == check);
	}

	for (i = 0; i < mb->nlinear; i++) {
		if (mb->linear_precision[i] > 0.0) {
			val += LOG_NORMC_GAUSSIAN + 1.0 / 2.0 * log(mb->linear_precision[i])
			    - 1.0 / 2.0 * mb->linear_precision[i] * SQR(mb->linear_mean[i]);
		}
	}


	typedef struct {
		GMRFLib_store_tp *store;
	} Store_tp;
	static Store_tp ***sstore = NULL;

	if (!sstore) {
#pragma omp critical (Name_87d8c02a8a06b017c5015b7132be14e8b5996507)
		if (!sstore) {
			sstore = Calloc(GMRFLib_CACHE_LEN(), Store_tp **);
		}
	}
	int cidx = 0;
	GMRFLib_CACHE_SET_ID(cidx);

	if (!sstore[cidx]) {
		sstore[cidx] = Calloc(mb->nf, Store_tp *);
	}

	for (i = 0; i < mb->nf; i++) {

		GMRFLib_store_tp *store = NULL;
		if (!sstore[cidx][i]) {
			sstore[cidx][i] = Calloc(1, Store_tp);
			sstore[cidx][i]->store = Calloc(1, GMRFLib_store_tp);
		}
		store = sstore[cidx][i]->store;

		switch (mb->f_id[i]) {
		case F_RW2D:
		case F_BESAG:
		case F_GENERIC0:
		case F_SEASONAL:
		case F_IID:
		case F_RW1:
		case F_RW2:
		case F_CRW2:
		{
			if (_NOT_FIXED(f_fixed[i][0])) {
				log_precision = theta[count];
				count++;
			} else {
				log_precision = mb->f_theta[i][0][thread_id][0];
			}
			_SET_GROUP_RHO(1);

			double scale_correction = 0.0;
			if (mb->f_id[i] == F_IID && mb->f_scale[i]) {
				int ii, nii = mb->f_N[i] / mb->f_ngroup[i];

				for (ii = 0; ii < nii; ii++) {
					scale_correction += log(mb->f_scale[i][ii]);
				}
				scale_correction /= nii;
			}

			val += mb->f_nrep[i] * (normc_g + gcorr * (LOG_NORMC_GAUSSIAN * (mb->f_N[i] - mb->f_rankdef[i]) +
								   (mb->f_N[i] - mb->f_rankdef[i]) / 2.0 * (log_precision + scale_correction)));
			if (_NOT_FIXED(f_fixed[i][0])) {
				val += PRIOR_EVAL(mb->f_prior[i][0], &log_precision);
			}
		}
			break;

		case F_SPDE:
		{
			int k, nT, nK, nt, Toffset = 0, Koffset = 0;
			double t;
			inla_spde_tp *spde = NULL;
			double *Tpar = NULL, *Kpar = NULL, init0, init1, init2, init3;

			spde = (inla_spde_tp *) mb->f_model[i];
			assert(spde->Qfunc_arg == spde);

			nT = spde->Tmodel->ntheta;
			nK = spde->Kmodel->ntheta;
			nt = mb->f_ntheta[i];

			fixed0 = !_NOT_FIXED(f_fixed[i][0]);
			fixed1 = !_NOT_FIXED(f_fixed[i][1]);
			fixed2 = !_NOT_FIXED(f_fixed[i][2]);
			fixed3 = !_NOT_FIXED(f_fixed[i][3]);
			init0 = mb->f_initial[i][0];
			init1 = mb->f_initial[i][1];
			init2 = mb->f_initial[i][2];
			init3 = mb->f_initial[i][3];

			Tpar = Calloc(nT, double);
			Kpar = Calloc(nK, double);

			if (debug) {
				P(i);
				P(mb->f_fixed[i][0]);
				P(mb->f_fixed[i][1]);
				P(mb->f_fixed[i][2]);
				P(mb->f_fixed[i][3]);
				P(nT);
				P(nK);
				P(nt);
				P(fixed0);
				P(fixed1);
				P(fixed2);
				P(fixed3);
				P(init0);
				P(init1);
				P(init2);
				P(init3);
			}

			if (nT) {
				if (fixed0) {
					Tpar[0] = init0;
				} else {
					Tpar[0] = theta[count];
					Toffset++;
				}
				if (nT > 1) {
					if (fixed2) {
						for (k = 1; k < nT; k++)
							Tpar[k] = init2;
					} else {
						for (k = 1; k < nT; k++) {
							Tpar[k] = theta[count + Toffset];
							Toffset++;
						}
					}
				}
				spde->Tmodel->theta_extra[thread_id] = Tpar;
			}

			if (debug) {
				P(count);
				P(Toffset);
				for (k = 0; k < nT; k++) {
					printf("Tpar[%d] = %g\n", k, Tpar[k]);
				}
			}

			if (nK) {
				if (fixed1) {
					Kpar[0] = init1;
				} else {
					Kpar[0] = theta[count + Toffset];
					Koffset++;
				}
				if (nK > 1) {
					if (fixed2) {
						for (k = 1; k < nK; k++)
							Kpar[k] = init2;
					} else {
						for (k = 1; k < nK; k++) {
							Kpar[k] = theta[count + Koffset + Toffset];
							Koffset++;
						}
					}
				}
				spde->Kmodel->theta_extra[thread_id] = Kpar;
			}

			if (debug) {
				P(count);
				P(Toffset);
				P(Koffset);
				for (k = 0; k < nT; k++) {
					printf("Tpar[%d] = %g\n", k, Tpar[k]);
				}
				for (k = 0; k < nK; k++) {
					printf("Kpar[%d] = %g\n", k, Kpar[k]);
				}
			}

			if (fixed3) {
				spde->oc[thread_id][0] = init3;
			} else {
				spde->oc[thread_id][0] = theta[count + Koffset + Toffset];
			}

			if (debug) {
				printf("call extra() with\n");
				for (k = 0; k < nT; k++) {
					printf("Tmodel %d %g\n", k, spde->Tmodel->theta_extra[thread_id][k]);
				}
				for (k = 0; k < nK; k++) {
					printf("Kmodel %d %g\n", k, spde->Kmodel->theta_extra[thread_id][k]);
				}
				printf("Oc %g\n", spde->oc[thread_id][0]);
			}

			/*
			 * T 
			 */
			if (nT) {
				if (_NOT_FIXED(f_fixed[i][0])) {
					t = theta[count];
					val += PRIOR_EVAL(mb->f_prior[i][0], &t);
					count++;
				}
				for (k = 1; k < nT; k++) {
					if (_NOT_FIXED(f_fixed[i][2])) {
						t = theta[count];
						val += PRIOR_EVAL(mb->f_prior[i][2], &t);
						count++;
					}
				}
			}

			/*
			 * K 
			 */
			if (nK) {
				if (_NOT_FIXED(f_fixed[i][1])) {
					t = theta[count];
					val += PRIOR_EVAL(mb->f_prior[i][1], &t);
					count++;
				}
				for (k = 1; k < nK; k++) {
					if (_NOT_FIXED(f_fixed[i][2])) {
						t = theta[count];
						val += PRIOR_EVAL(mb->f_prior[i][2], &t);
						count++;
					}
				}
			}
			/*
			 * Ocillating coeff 
			 */
			if (!fixed3) {
				t = theta[count];
				val += PRIOR_EVAL(mb->f_prior[i][3], &t);
				count++;
			}

			if (debug) {
				P(nT);
				P(nK);
				P(mb->f_ntheta[i]);
			}
			assert(IMAX(4, nT + nK + 1) + (mb->f_ngroup[i] > 1 ? 1 : 0) == mb->f_ntheta[i]);

			_SET_GROUP_RHO(IMAX(4, nT + nK + 1));

			GMRFLib_problem_tp *problem = NULL;
			/*
			 * do a check for numerical not pos def matrix here, as its so close to being singular 
			 */
			int retval = GMRFLib_SUCCESS, ok = 0, num_try = 0, num_try_max = 100;;
			GMRFLib_error_handler_tp *old_handler = GMRFLib_set_error_handler_off();
			double *cc_add = Calloc(spde->graph->n, double);

			if (mb->f_diag[i]) {
				int ii;
				for (ii = 0; ii < spde->graph->n; ii++) {
					cc_add[ii] = mb->f_diag[i];
				}
			}

			while (!ok) {
				retval = GMRFLib_init_problem_store(thread_id, &problem, NULL, NULL, cc_add, NULL,
								    spde->graph, spde->Qfunc, spde->Qfunc_arg, mb->f_constr_orig[i], store);
				switch (retval) {
				case GMRFLib_EPOSDEF:
				{
					int ii;
					double eps = GSL_SQRT_DBL_EPSILON;

					for (ii = 0; ii < spde->graph->n; ii++) {
						cc_add[ii] = (cc_add[ii] == 0.0 ? eps : cc_add[ii] * 10.0);
					}
					break;
				}

				case GMRFLib_SUCCESS:
					ok = 1;
					break;

				default:
					/*
					 * some other error 
					 */
					GMRFLib_set_error_handler(old_handler);
					abort();
					break;
				}

				if (++num_try >= num_try_max) {
					FIXME("This should not happen. Contact developers...");
					abort();
				}
			}
			Free(cc_add);
			GMRFLib_set_error_handler(old_handler);

			GMRFLib_evaluate(problem);
			val += mb->f_nrep[i] * (problem->sub_logdens * (ngroup - grankdef) + normc_g);

			if (nT) {
				spde->Tmodel->theta_extra[thread_id] = NULL;
			}
			if (nK) {
				spde->Kmodel->theta_extra[thread_id] = NULL;
			}

			GMRFLib_free_problem(problem);
			Free(Tpar);
			Free(Kpar);
		}
			break;

		case F_SPDE2:
		{
			int k, kk;
			inla_spde2_tp *spde2 = NULL;

			spde2 = (inla_spde2_tp *) mb->f_model[i];
			assert(spde2->Qfunc_arg == spde2);

			if (debug) {
				static int first = 1;
				if (first) {
					P(spde2->ntheta_used);
					P(spde2->ntheta);
					P(mb->f_ntheta[i]);
					first = 0;
				}
			}

			spde2->debug = 0;
			if (!mb->mode_fixed) {
				for (k = kk = 0; k < spde2->ntheta_used; k++, kk++) {
					while (mb->f_fixed[i][kk])
						kk++;
					spde2->theta[kk][thread_id][0] = theta[count + k];
				}
			}
			int count_ref = count;

			if (!mb->mode_fixed) {
				count += spde2->ntheta_used;   /* as _SET_GROUP_RHO need 'count' */
			}
			_SET_GROUP_RHO(spde2->ntheta);

			/*
			 * do a check for numerical not pos def matrix here, as its so close to being singular 
			 */
			GMRFLib_problem_tp *problem = NULL;
			int retval = GMRFLib_SUCCESS, ok = 0, num_try = 0, num_try_max = 100;
			GMRFLib_error_handler_tp *old_handler = GMRFLib_set_error_handler_off();
			double *cc_add = Calloc(spde2->graph->n, double);

			if (mb->f_diag[i]) {
				int ii;
				for (ii = 0; ii < spde2->graph->n; ii++) {
					cc_add[ii] = mb->f_diag[i];
				}
			}

			while (!ok) {
				retval = GMRFLib_init_problem_store(thread_id, &problem, NULL, NULL, cc_add, NULL,
								    spde2->graph, spde2->Qfunc, spde2->Qfunc_arg, mb->f_constr_orig[i], store);
				switch (retval) {
				case GMRFLib_EPOSDEF:
				{
					int ii;
					double eps = GSL_SQRT_DBL_EPSILON;

					for (ii = 0; ii < spde2->graph->n; ii++) {
						cc_add[ii] = (cc_add[ii] == 0.0 ? eps : cc_add[ii] * 10.0);
					}
					break;
				}

				case GMRFLib_SUCCESS:
					ok = 1;
					break;

				default:
					/*
					 * some other error 
					 */
					GMRFLib_set_error_handler(old_handler);
					abort();
					break;
				}

				if (++num_try >= num_try_max) {
					FIXME("This should not happen. Contact developers...");
					abort();
				}
			}
			Free(cc_add);
			GMRFLib_set_error_handler(old_handler);

			GMRFLib_evaluate(problem);
			val += mb->f_nrep[i] * (problem->sub_logdens * (ngroup - grankdef) + normc_g);

			/*
			 * this is a multivariate prior...  'count_ref' is the 'first theta'
			 */
			if (!mb->mode_fixed) {
				if (mb->f_prior[i][0].id == P_PC_MATERN) {
					/*
					 *  This is a special case: the pc_matern prior. pass NAN for fixed values and the prior will do the correct thing.
					 */
					double local_theta[2];
					int local_count = 0;

					local_theta[0] = (_NOT_FIXED(f_fixed[i][0]) ? theta[count_ref + local_count++] : NAN);
					local_theta[1] = (_NOT_FIXED(f_fixed[i][1]) ? theta[count_ref + local_count++] : NAN);
					assert(local_count == spde2->ntheta_used);
					val += PRIOR_EVAL(mb->f_prior[i][0], local_theta);
				} else {
					// the mvnorm prior, defined on the _USED_ thetas!
					val += PRIOR_EVAL(mb->f_prior[i][0], &theta[count_ref]);
				}
			}

			GMRFLib_free_problem(problem);
		}
			break;

		case F_SPDE3:
		{
			int k, spde3_ntheta;
			inla_spde3_tp *spde3 = NULL;

			spde3 = (inla_spde3_tp *) mb->f_model[i];
			assert(spde3->Qfunc_arg == spde3);

			spde3->debug = 0;
			spde3_ntheta = spde3->ntheta;
			if (!mb->mode_fixed) {
				for (k = 0; k < spde3_ntheta; k++) {
					spde3->theta[k][thread_id][0] = theta[count + k];
				}
			}
			int count_ref = count;

			if (!mb->mode_fixed) {
				count += spde3_ntheta;	       /* as _SET_GROUP_RHO need 'count' */
			}
			_SET_GROUP_RHO(spde3_ntheta);

			/*
			 * do a check for numerical not pos def matrix here, as its so close to being singular 
			 */
			GMRFLib_problem_tp *problem = NULL;
			int retval = GMRFLib_SUCCESS, ok = 0, num_try = 0, num_try_max = 100;
			GMRFLib_error_handler_tp *old_handler = GMRFLib_set_error_handler_off();
			double *cc_add = Calloc(spde3->graph->n, double);

			if (mb->f_diag[i]) {
				int ii;
				for (ii = 0; ii < spde3->graph->n; ii++) {
					cc_add[ii] = mb->f_diag[i];
				}
			}

			while (!ok) {
				retval = GMRFLib_init_problem_store(thread_id, &problem, NULL, NULL, cc_add, NULL,
								    spde3->graph, spde3->Qfunc, spde3->Qfunc_arg, mb->f_constr_orig[i], store);
				switch (retval) {
				case GMRFLib_EPOSDEF:
				{
					int ii;
					double eps = GSL_SQRT_DBL_EPSILON;

					for (ii = 0; ii < spde3->graph->n; ii++) {
						cc_add[ii] = (cc_add[ii] == 0.0 ? eps : cc_add[ii] * 10.0);
					}
					break;
				}

				case GMRFLib_SUCCESS:
					ok = 1;
					break;

				default:
					/*
					 * some other error 
					 */
					GMRFLib_set_error_handler(old_handler);
					abort();
					break;
				}

				if (++num_try >= num_try_max) {
					FIXME("This should not happen. Contact developers...");
					abort();
				}
			}
			Free(cc_add);
			GMRFLib_set_error_handler(old_handler);

			GMRFLib_evaluate(problem);
			val += mb->f_nrep[i] * (problem->sub_logdens * (ngroup - grankdef) + normc_g);

			/*
			 * this is the mvnormal prior...  'count_ref' is the 'first theta as this is a mutivariate prior.
			 */
			if (!mb->mode_fixed) {
				val += PRIOR_EVAL(mb->f_prior[i][0], &theta[count_ref]);
			}

			GMRFLib_free_problem(problem);
		}
			break;

		case F_AR:
		{
			double *pacf = NULL, *pacf_intern = NULL;
			int p;

			p = mb->f_order[i];
			if (_NOT_FIXED(f_fixed[i][0])) {
				log_precision = theta[count];
				count++;
			} else {
				log_precision = mb->f_theta[i][0][thread_id][0];
			}
			pacf = Calloc(p, double);
			pacf_intern = Calloc(p, double);
			for (j = 0; j < p; j++) {
				if (_NOT_FIXED(f_fixed[i][j + 1])) {
					pacf_intern[j] = theta[count];
					count++;
				} else {
					pacf_intern[j] = mb->f_theta[i][j + 1][thread_id][0];
				}
				pacf[j] = ar_map_pacf(pacf_intern[j], MAP_FORWARD, NULL);
			}
			_SET_GROUP_RHO(AR_MAXTHETA + 1);

			int n_ar;
			double marginal_prec, conditional_prec, *marginal_Q = NULL, *param = NULL, *zero = NULL, ldens;

			marginal_Q = Calloc(ISQR(p), double);
			ar_marginal_distribution(p, pacf, &marginal_prec, marginal_Q);
			conditional_prec = exp(log_precision) / marginal_prec;

			param = Calloc(1 + p + ISQR(p), double);
			zero = Calloc(p, double);
			param[0] = p;
			for (j = 0; j < ISQR(p); j++) {
				param[1 + p + j] = marginal_Q[j] * exp(log_precision);
			}

			n_ar = mb->f_n[i] / mb->f_ngroup[i];
			ldens = priorfunc_mvnorm(zero, param) +
			    (n_ar - p - mb->f_rankdef[i]) * (-0.5 * log(2 * M_PI) + 0.5 * log(conditional_prec));
			val += mb->f_nrep[i] * (ldens * (ngroup - grankdef) + normc_g);

			if (_NOT_FIXED(f_fixed[i][0])) {
				val += PRIOR_EVAL(mb->f_prior[i][0], &log_precision);
			}

			if (mb->f_prior[i][1].id == P_REF_AR) {
				// joint reference prior
				val += PRIOR_EVAL(mb->f_prior[i][1], pacf_intern);
			} else {
				// sequential pc prior. No need to use AR_MAXTHETA there
				for (j = 0; j < p; j++) {
					if (_NOT_FIXED(f_fixed[i][j + 1])) {
						val += PRIOR_EVAL(mb->f_prior[i][j + 1], &(pacf_intern[j]));
					}
				}
			}
			Free(pacf);
			Free(pacf_intern);
			Free(marginal_Q);
			Free(param);
			Free(zero);
		}
			break;

		case F_GENERIC1:
		{
			if (_NOT_FIXED(f_fixed[i][0])) {
				log_precision = theta[count];
				count++;
			} else {
				log_precision = mb->f_theta[i][0][thread_id][0];
			}
			if (_NOT_FIXED(f_fixed[i][1])) {
				beta_intern = theta[count];
				count++;
			} else {
				beta_intern = mb->f_theta[i][1][thread_id][0];
			}
			beta = map_probability(beta_intern, MAP_FORWARD, NULL);
			_SET_GROUP_RHO(2);

			double logdet_Q = 0.0;
			inla_generic1_tp *a = (inla_generic1_tp *) mb->f_Qfunc_arg[i];
			for (j = 0; j < n_orig; j++) {
				logdet_Q += LOG_1mp(beta * a->eigenvalues[j] / a->max_eigenvalue);
			}

			val += mb->f_nrep[i] * (normc_g + gcorr * (LOG_NORMC_GAUSSIAN * (mb->f_n[i] - mb->f_rankdef[i])
								   + (mb->f_n[i] - mb->f_rankdef[i]) / 2.0 * log_precision +
								   ngroup * 0.5 * logdet_Q));
			if (_NOT_FIXED(f_fixed[i][0])) {
				val += PRIOR_EVAL(mb->f_prior[i][0], &log_precision);
			}
			if (_NOT_FIXED(f_fixed[i][1])) {
				val += PRIOR_EVAL(mb->f_prior[i][1], &beta_intern);
			}
		}
			break;

		case F_GENERIC2:
		{
			/*
			 * OOPS: even though the parameters are (log_prec, h2_inter), the prior is defined on (log_prec, log_prec_unstruct), with the proper
			 * Jacobian added. 
			 */
			double h2;

			if (_NOT_FIXED(f_fixed[i][0])) {
				log_precision = theta[count];
				count++;
			} else {
				log_precision = mb->f_theta[i][0][thread_id][0];
			}
			if (_NOT_FIXED(f_fixed[i][1])) {
				h2_intern = theta[count];
				count++;
			} else {
				h2_intern = mb->f_theta[i][1][thread_id][0];
			}
			h2 = map_probability(h2_intern, MAP_FORWARD, NULL);
			_SET_GROUP_RHO(2);

			double log_prec_unstruct = log(h2 / (1.0 - h2)) + log_precision;
			double n = (double) mb->f_n[i];

			val += mb->f_nrep[i] * (normc_g + gcorr * (LOG_NORMC_GAUSSIAN * (n / 2.0 + (n - mb->f_rankdef[i]) / 2.0) +
								   +(n - mb->f_rankdef[i]) / 2.0 * log_precision + n / 2.0 * log_prec_unstruct));

			if (_NOT_FIXED(f_fixed[i][0])) {
				val += PRIOR_EVAL(mb->f_prior[i][0], &log_precision);
			}
			if (_NOT_FIXED(f_fixed[i][1])) {
				val += PRIOR_EVAL(mb->f_prior[i][1], &log_prec_unstruct);
			}
			/*
			 * The Jacobian for the change of variables, is
			 * 
			 * | d log_prec_unstruct / d h2_intern | = 1, so no need to correct for the Jacobian from the change of variables. 
			 */
		}
			break;

		case F_FGN:
		{
			double H_intern;

			if (_NOT_FIXED(f_fixed[i][0])) {
				log_precision = theta[count];
				count++;
			} else {
				log_precision = mb->f_theta[i][0][thread_id][0];
			}
			if (_NOT_FIXED(f_fixed[i][1])) {
				H_intern = theta[count];
				count++;
			} else {
				H_intern = mb->f_theta[i][1][thread_id][0];
			}
			_SET_GROUP_RHO(2);

			inla_fgn_arg_tp *arg = (inla_fgn_arg_tp *) mb->f_Qfunc_arg_orig[i];

			arg->log_prec[thread_id][0] = log_precision;
			arg->H_intern[thread_id][0] = H_intern;

			int n = mb->f_graph_orig[i]->n;

			GMRFLib_problem_tp *problem = NULL;
			int retval = GMRFLib_SUCCESS, ok = 0, num_try = 0, num_try_max = 100;
			GMRFLib_error_handler_tp *old_handler = GMRFLib_set_error_handler_off();
			double *cc_add = Calloc(n, double);

			if (mb->f_diag[i]) {
				int ii;
				for (ii = 0; ii < n; ii++) {
					cc_add[ii] = mb->f_diag[i];
				}
			}

			while (!ok) {
				retval =
				    GMRFLib_init_problem_store(thread_id, &problem, NULL, NULL, cc_add, NULL, mb->f_graph_orig[i],
							       mb->f_Qfunc_orig[i], mb->f_Qfunc_arg_orig[i], mb->f_constr_orig[i], store);
				switch (retval) {
				case GMRFLib_EPOSDEF:
				{
					int ii;
					double eps = GSL_SQRT_DBL_EPSILON;

					for (ii = 0; ii < arg->n; ii++) {
						cc_add[ii] = (cc_add[ii] == 0.0 ? eps : cc_add[ii] * 10.0);
					}
				}
					break;

				case GMRFLib_SUCCESS:
				{
					ok = 1;
				}
					break;

				default:
					/*
					 * some other error 
					 */
					GMRFLib_set_error_handler(old_handler);
					abort();
					break;
				}

				if (++num_try >= num_try_max) {
					FIXME("This should not happen. Contact developers...");
					abort();
				}
			}
			Free(cc_add);
			GMRFLib_set_error_handler(old_handler);

			GMRFLib_evaluate(problem);
			val += mb->f_nrep[i] * (problem->sub_logdens * (ngroup - grankdef) + normc_g);

			if (_NOT_FIXED(f_fixed[i][0])) {
				val += PRIOR_EVAL(mb->f_prior[i][0], &log_precision);
			}
			if (_NOT_FIXED(f_fixed[i][1])) {
				val += PRIOR_EVAL(mb->f_prior[i][1], &H_intern);
			}
			GMRFLib_free_problem(problem);
		}
			break;

		case F_FGN2:
		{
			double H_intern;

			if (_NOT_FIXED(f_fixed[i][0])) {
				log_precision = theta[count];
				count++;
			} else {
				log_precision = mb->f_theta[i][0][thread_id][0];
			}
			if (_NOT_FIXED(f_fixed[i][1])) {
				H_intern = theta[count];
				count++;
			} else {
				H_intern = mb->f_theta[i][1][thread_id][0];
			}
			_SET_GROUP_RHO(2);

			inla_fgn2_arg_tp *arg = (inla_fgn2_arg_tp *) mb->f_Qfunc_arg_orig[i];

			arg->log_prec[thread_id][0] = log_precision;
			arg->H_intern[thread_id][0] = H_intern;

			int n = mb->f_graph_orig[i]->n;

			/*
			 * do a check for numerical not pos def matrix here, as it may be close to being singular 
			 */
			GMRFLib_problem_tp *problem = NULL;
			int retval = GMRFLib_SUCCESS, ok = 0, num_try = 0, num_try_max = 100;
			GMRFLib_error_handler_tp *old_handler = GMRFLib_set_error_handler_off();
			double *cc_add = Calloc(n, double);

			if (mb->f_diag[i]) {
				int ii;
				for (ii = 0; ii < n; ii++) {
					cc_add[ii] = mb->f_diag[i];
				}
			}

			while (!ok) {
				retval =
				    GMRFLib_init_problem_store(thread_id, &problem, NULL, NULL, cc_add, NULL, mb->f_graph_orig[i],
							       mb->f_Qfunc_orig[i], mb->f_Qfunc_arg_orig[i], mb->f_constr_orig[i], store);
				switch (retval) {
				case GMRFLib_EPOSDEF:
				{
					int ii;
					double eps = GSL_SQRT_DBL_EPSILON;

					for (ii = 0; ii < arg->n; ii++) {
						cc_add[ii] = (cc_add[ii] == 0.0 ? eps : cc_add[ii] * 10.0);
					}
					break;
				}

				case GMRFLib_SUCCESS:
					ok = 1;
					break;

				default:
					/*
					 * some other error 
					 */
					GMRFLib_set_error_handler(old_handler);
					abort();
					break;
				}

				if (++num_try >= num_try_max) {
					FIXME("This should not happen. Contact developers...");
					abort();
				}
			}
			Free(cc_add);
			GMRFLib_set_error_handler(old_handler);

			GMRFLib_evaluate(problem);
			val += mb->f_nrep[i] * (problem->sub_logdens * (ngroup - grankdef) + normc_g);

			if (_NOT_FIXED(f_fixed[i][0])) {
				val += PRIOR_EVAL(mb->f_prior[i][0], &log_precision);
			}
			if (_NOT_FIXED(f_fixed[i][1])) {
				val += PRIOR_EVAL(mb->f_prior[i][1], &H_intern);
			}
			GMRFLib_free_problem(problem);
		}
			break;

		case F_Z:
		{
			if (_NOT_FIXED(f_fixed[i][0])) {
				log_precision = theta[count];
				count++;
			} else {
				log_precision = mb->f_theta[i][0][thread_id][0];
			}
			_SET_GROUP_RHO(1);

			inla_z_arg_tp *arg = (inla_z_arg_tp *) mb->f_Qfunc_arg_orig[i];
			arg->log_prec[thread_id][0] = log_precision;

			GMRFLib_problem_tp *problem = NULL;
			int retval = GMRFLib_SUCCESS, ok = 0, num_try = 0, num_try_max = 100;
			GMRFLib_error_handler_tp *old_handler = GMRFLib_set_error_handler_off();
			double *cc_add = Calloc(arg->n + arg->m, double);

			assert(mb->f_graph_orig[i]->n == arg->n + arg->m);

			if (mb->f_diag[i]) {
				int ii;
				for (ii = 0; ii < arg->n + arg->m; ii++) {
					cc_add[ii] = mb->f_diag[i];
				}
			}

			while (!ok) {
				retval =
				    GMRFLib_init_problem_store(thread_id, &problem, NULL, NULL, cc_add, NULL, mb->f_graph_orig[i],
							       mb->f_Qfunc_orig[i], mb->f_Qfunc_arg_orig[i], mb->f_constr_orig[i], store);
				switch (retval) {
				case GMRFLib_EPOSDEF:
				{
					int ii;
					double eps = GSL_SQRT_DBL_EPSILON;

					/*
					 * only need to add for the z-part; the last m components.
					 */
					for (ii = arg->n; ii < arg->n + arg->m; ii++) {
						cc_add[ii] = (cc_add[ii] == 0.0 ? eps : cc_add[ii] * 10.0);
					}
					break;
				}

				case GMRFLib_SUCCESS:
					ok = 1;
					break;

				default:
					/*
					 * some other error 
					 */
					GMRFLib_set_error_handler(old_handler);
					abort();
					break;
				}

				if (++num_try >= num_try_max) {
					FIXME("This should not happen. Contact developers...");
					abort();
				}
			}
			Free(cc_add);
			GMRFLib_set_error_handler(old_handler);

			GMRFLib_evaluate(problem);
			val += mb->f_nrep[i] * (problem->sub_logdens * (ngroup - grankdef) + normc_g);

			if (_NOT_FIXED(f_fixed[i][0])) {
				val += PRIOR_EVAL(mb->f_prior[i][0], &log_precision);
			}

			GMRFLib_free_problem(problem);
		}
			break;

		case F_GENERIC3:
		{
			int k, kk;
			inla_generic3_tp *a = NULL;

			a = Calloc(1, inla_generic3_tp);
			Memcpy((void *) a, (void *) mb->f_Qfunc_arg_orig[i], sizeof(inla_generic3_tp));
			a->log_prec = Calloc(GENERIC3_MAXTHETA, double **);
			for (k = 0; k < GENERIC3_MAXTHETA; k++) {
				if (_NOT_FIXED(f_fixed[i][k])) {
					HYPER_NEW(a->log_prec[k], theta[count]);
					count++;
				} else {
					HYPER_NEW(a->log_prec[k], mb->f_theta[i][k][thread_id][0]);
				}
			}
			_SET_GROUP_RHO(GENERIC3_MAXTHETA);

			GMRFLib_problem_tp *problem = NULL;
			int retval = GMRFLib_SUCCESS, ok = 0, num_try = 0, num_try_max = 100;
			GMRFLib_error_handler_tp *old_handler = GMRFLib_set_error_handler_off();
			double *cc_add = Calloc(a->n, double);

			if (mb->f_diag[i]) {
				int ii;
				for (ii = 0; ii < a->n; ii++) {
					cc_add[ii] = mb->f_diag[i];
				}
			}

			while (!ok) {
				retval = GMRFLib_init_problem_store(thread_id, &problem, NULL, NULL, cc_add, NULL,
								    mb->f_graph_orig[i], mb->f_Qfunc_orig[i], (void *) a, mb->f_constr_orig[i],
								    store);
				switch (retval) {
				case GMRFLib_EPOSDEF:
				{
					int ii;
					double eps = GSL_SQRT_DBL_EPSILON;

					/*
					 * only need to add for the z-part; the last m components.
					 */
					for (ii = 0; ii < a->n; ii++) {
						cc_add[ii] = (cc_add[ii] == 0.0 ? eps : cc_add[ii] * 10.0);
					}
					break;
				}

				case GMRFLib_SUCCESS:
					ok = 1;
					break;

				default:
					/*
					 * some other error 
					 */
					GMRFLib_set_error_handler(old_handler);
					abort();
					break;
				}

				if (++num_try >= num_try_max) {
					FIXME("This should not happen. Contact developers...");
					abort();
				}
			}
			Free(cc_add);
			GMRFLib_set_error_handler(old_handler);

			GMRFLib_evaluate(problem);
			val += mb->f_nrep[i] * (problem->sub_logdens * (ngroup - grankdef) + normc_g);

			for (k = 0; k < GENERIC3_MAXTHETA; k++) {
				if (_NOT_FIXED(f_fixed[i][k])) {
					log_precision = a->log_prec[k][0][0];
					val += PRIOR_EVAL(mb->f_prior[i][k], &log_precision);
				}
			}

			/*
			 * Free 'a' 
			 */
			for (k = 0; k < GENERIC3_MAXTHETA; k++) {
				for (kk = 0; kk < GMRFLib_MAX_THREADS(); kk++) {
					Free(a->log_prec[k][kk]);
				}
				Free(a->log_prec[k]);
			}
			Free(a->log_prec);
			Free(a);
			GMRFLib_free_problem(problem);
		}
			break;

		case F_SLM:
		{
			if (_NOT_FIXED(f_fixed[i][0])) {
				log_precision = theta[count];
				count++;
			} else {
				log_precision = mb->f_theta[i][0][thread_id][0];
			}
			if (_NOT_FIXED(f_fixed[i][1])) {
				logit_rho = theta[count];
				count++;
			} else {
				logit_rho = mb->f_theta[i][1][thread_id][0];
			}
			_SET_GROUP_RHO(2);

			inla_slm_arg_tp *arg = (inla_slm_arg_tp *) mb->f_Qfunc_arg_orig[i];
			arg->log_prec[thread_id][0] = log_precision;
			arg->logit_rho[thread_id][0] = logit_rho;

			GMRFLib_problem_tp *problem = NULL;
			int retval = GMRFLib_SUCCESS, ok = 0, num_try = 0, num_try_max = 100;
			GMRFLib_error_handler_tp *old_handler = GMRFLib_set_error_handler_off();
			double *cc_add = Calloc(arg->n + arg->m, double);

			assert(mb->f_graph_orig[i]->n == arg->n + arg->m);

			if (mb->f_diag[i]) {
				int ii;
				for (ii = 0; ii < arg->n + arg->m; ii++) {
					cc_add[ii] = mb->f_diag[i];
				}
			}

			while (!ok) {
				retval = GMRFLib_init_problem_store(thread_id, &problem, NULL, NULL, cc_add, NULL,
								    mb->f_graph_orig[i],
								    mb->f_Qfunc_orig[i], mb->f_Qfunc_arg_orig[i], mb->f_constr_orig[i], store);
				switch (retval) {
				case GMRFLib_EPOSDEF:
				{
					int ii;
					double eps = GSL_SQRT_DBL_EPSILON;
					for (ii = 0; ii < mb->f_graph_orig[i]->n; ii++) {
						cc_add[ii] = (cc_add[ii] == 0.0 ? eps : cc_add[ii] * 10.0);
					}
					break;
				}

				case GMRFLib_SUCCESS:
					ok = 1;
					break;

				default:
					/*
					 * some other error 
					 */
					GMRFLib_set_error_handler(old_handler);
					abort();
					break;
				}

				if (++num_try >= num_try_max) {
					FIXME("This should not happen. Contact developers...");
					abort();
				}
			}
			Free(cc_add);
			GMRFLib_set_error_handler(old_handler);

			GMRFLib_evaluate(problem);
			val += mb->f_nrep[i] * (problem->sub_logdens * (ngroup - grankdef) + normc_g);

			if (_NOT_FIXED(f_fixed[i][0])) {
				val += PRIOR_EVAL(mb->f_prior[i][0], &log_precision);
			}
			if (_NOT_FIXED(f_fixed[i][1])) {
				val += PRIOR_EVAL(mb->f_prior[i][1], &logit_rho);
			}
			GMRFLib_free_problem(problem);
		}
			break;

		case F_MEC:
		{
			double mean_x, log_precision_x, log_precision_obs;

			if (_NOT_FIXED(f_fixed[i][0])) {
				beta_intern = theta[count];
				val += PRIOR_EVAL(mb->f_prior[i][0], &beta_intern);
				count++;
			} else {
				beta_intern = mb->f_theta[i][0][thread_id][0];
			}
			beta = mb->f_theta_map[i][0] (beta_intern, MAP_FORWARD, mb->f_theta_map_arg[i][0]);

			if (_NOT_FIXED(f_fixed[i][1])) {
				log_precision_obs = theta[count];
				val += PRIOR_EVAL(mb->f_prior[i][1], &log_precision_obs);
				count++;
			} else {
				log_precision_obs = mb->f_theta[i][1][thread_id][0];
			}

			if (_NOT_FIXED(f_fixed[i][2])) {
				mean_x = theta[count];
				val += PRIOR_EVAL(mb->f_prior[i][2], &mean_x);
				count++;
			} else {
				mean_x = mb->f_theta[i][2][thread_id][0];
			}
			if (_NOT_FIXED(f_fixed[i][3])) {
				log_precision_x = theta[count];
				val += PRIOR_EVAL(mb->f_prior[i][3], &log_precision_x);
				count++;
			} else {
				log_precision_x = mb->f_theta[i][3][thread_id][0];
			}

			_SET_GROUP_RHO(4);

			double AA = 0.0, BB = 0.0, CC = 0.0;
			int ii, nii = mb->f_N[i] / mb->f_ngroup[i];

			/*
			 * we need to do it like this as the scale[i] varies with 'i'.
			 */
			assert(mb->f_scale[i]);
			for (ii = 0; ii < nii; ii++) {
				AA += log((exp(log_precision_obs) * mb->f_scale[i][ii] + exp(log_precision_x)) / SQR(beta));
				BB += log(1.0 / (1.0 / (mb->f_scale[i][ii] * exp(log_precision_obs)) + 1.0 / exp(log_precision_x)));
				CC += SQR(mb->f_locations[i][ii] - mean_x) *
				    (1.0 / (1.0 / (mb->f_scale[i][ii] * exp(log_precision_obs)) + 1.0 / exp(log_precision_x)));
			}
			val += mb->f_nrep[i] * (normc_g + gcorr * (2.0 * LOG_NORMC_GAUSSIAN * (mb->f_N[i] - mb->f_rankdef[i])
								   + mb->f_ngroup[i] * 0.5 * (AA + BB)
								   /*
								    * and the exponent which depends on the precisions. we need to correct with ngroup here as its
								    * not scaled with f_N that already is corrected for ngroup.
								    */
								   - mb->f_ngroup[i] * 0.5 * CC));
		}
			break;

		case F_MEB:
		{
			if (_NOT_FIXED(f_fixed[i][0])) {
				beta_intern = theta[count];
				val += PRIOR_EVAL(mb->f_prior[i][0], &beta_intern);
				count++;
			} else {
				beta_intern = mb->f_theta[i][0][thread_id][0];
			}
			beta = mb->f_theta_map[i][0] (beta_intern, MAP_FORWARD, mb->f_theta_map_arg[i][0]);

			if (_NOT_FIXED(f_fixed[i][1])) {
				log_precision = theta[count];
				val += PRIOR_EVAL(mb->f_prior[i][1], &log_precision);
				count++;
			} else {
				log_precision = mb->f_theta[i][1][thread_id][0];
			}

			_SET_GROUP_RHO(2);

			double scale_correction = 0.0;
			int ii, nii = mb->f_N[i] / mb->f_ngroup[i];

			/*
			 * we need to do it like this as the scale[i] varies with 'i'.
			 */
			if (mb->f_scale[i]) {
				for (ii = 0; ii < nii; ii++) {
					scale_correction += log(mb->f_scale[i][ii]);
				}
				scale_correction /= nii;
			}
			val += mb->f_nrep[i] * (normc_g + gcorr * (LOG_NORMC_GAUSSIAN * (mb->f_N[i] - mb->f_rankdef[i]) +
								   (mb->f_N[i] - mb->f_rankdef[i]) / 2.0 * (log_precision +
													    scale_correction - log(SQR(beta)))));
		}
			break;

		case F_R_GENERIC:
		{
			int n_out, nn_out, ii, nt = 0, all_fixed = 0;

			double *x_out = NULL, *xx_out = NULL, *param = NULL, log_norm_const = 0.0, log_prior = 0.0;
			inla_rgeneric_tp *def = NULL;
			def = (inla_rgeneric_tp *) mb->f_Qfunc_arg_orig[i];

			nt = def->ntheta;
			param = Calloc(nt, double);
			if (nt) {
				all_fixed = 1;
				for (ii = 0; ii < nt; ii++) {
					if (_NOT_FIXED(f_fixed[i][ii])) {
						param[ii] = theta[count];
						all_fixed = 0;
						count++;
					} else {
						param[ii] = mb->f_theta[i][ii][thread_id][0];
					}
				}
			}

			inla_R_rgeneric(&n_out, &x_out, R_GENERIC_LOG_NORM_CONST, def->model, &nt, param);
			inla_R_rgeneric(&nn_out, &xx_out, R_GENERIC_LOG_PRIOR, def->model, &nt, param);

			switch (nn_out) {
			case 0:
				log_prior = 0.0;
				break;
			case 1:
				// we need to add a check for 'all_fixed' here, as with control.mode=list(...,fixed=TRUE) will
				// trigger all_fixed=1.
				log_prior = (evaluate_hyper_prior && !all_fixed ? xx_out[0] : 0.0);
				break;
			default:
				assert(0 == 1);
			}
			if (nn_out) {
				Free(xx_out);
			}

			switch (n_out) {
			case 0:
			{
				/*
				 * if it is the standard norm.const, the user can request us to compute it here if numeric(0) is
				 * returned from R_rgeneric.
				 */
				int *ilist = NULL, *jlist = NULL, n, len, k = 0, jj;
				double *Qijlist = NULL;
				GMRFLib_tabulate_Qfunc_tp *Qf = NULL;

				inla_R_rgeneric(&nn_out, &xx_out, R_GENERIC_Q, def->model, &nt, param);
				assert(nn_out >= 2);

				if ((int) xx_out[0] == -1) {
					// optimized output
					k = 1;
					len = (int) xx_out[k++];
					assert(len == def->len_list);
					assert(def->graph);
					assert(def->ilist);
					assert(def->jlist);
					n = def->graph->n;
					GMRFLib_tabulate_Qfunc_from_list2(&Qf, def->graph, def->len_list, def->ilist, def->jlist, &(xx_out[k]),
									  def->graph->n, NULL);
				} else {
					n = (int) xx_out[k++];
					len = (int) xx_out[k++];
					ilist = (int *) &xx_out[k];
					jlist = (int *) &xx_out[k + len];
					Qijlist = &xx_out[k + 2 * len];
					for (jj = 0; jj < len; jj++) {
						ilist[jj] = (int) xx_out[k + jj];
						jlist[jj] = (int) xx_out[k + len + jj];
					}
					GMRFLib_tabulate_Qfunc_from_list2(&Qf, def->graph, len, ilist, jlist, Qijlist, n, NULL);
				}

				int retval = GMRFLib_SUCCESS, ok = 0, num_try = 0, num_try_max = 100;
				GMRFLib_problem_tp *problem = NULL;
				GMRFLib_error_handler_tp *old_handler = GMRFLib_set_error_handler_off();
				double *cc_add = Calloc(n, double);

				if (mb->f_diag[i]) {
					for (jj = 0; jj < n; jj++) {
						cc_add[jj] = mb->f_diag[i];
					}
				}

				while (!ok) {
					retval = GMRFLib_init_problem_store(thread_id, &problem, NULL, NULL, cc_add, NULL,
									    def->graph, Qf->Qfunc, Qf->Qfunc_arg, mb->f_constr_orig[i], store);
					switch (retval) {
					case GMRFLib_EPOSDEF:
					{
						double eps = GSL_SQRT_DBL_EPSILON;
						for (jj = 0; jj < n; jj++) {
							cc_add[jj] = (cc_add[jj] == 0.0 ? eps : cc_add[jj] * 10.0);
						}
						break;
					}

					case GMRFLib_SUCCESS:
						ok = 1;
						break;

					default:
						/*
						 * some other error 
						 */
						GMRFLib_set_error_handler(old_handler);
						assert(0 == 1);
						abort();
					}

					if (++num_try >= num_try_max) {
						FIXME("This should not happen. Contact developers...");
						abort();
					}
				}
				Free(cc_add);
				GMRFLib_set_error_handler(old_handler);
				GMRFLib_evaluate(problem);
				log_norm_const = problem->sub_logdens;

				GMRFLib_free_problem(problem);
				GMRFLib_free_tabulate_Qfunc(Qf);
				Free(xx_out);
			}
				break;

			case 1:
			{
				log_norm_const = x_out[0];
				break;
			}

			default:
				assert(0 == 1);
			}

			if (debug) {
				for (ii = 0; ii < nt; ii++) {
					printf("p %.12g ", param[ii]);
				}
				printf(" %.12g  prior %.12g\n", log_norm_const, log_prior);
			}

			_SET_GROUP_RHO(nt);
			val += mb->f_nrep[i] * (normc_g + log_norm_const * (mb->f_ngroup[i] - grankdef)) + log_prior;

			Free(param);
			if (n_out) {
				Free(x_out);
			}
		}
			break;

		case F_C_GENERIC:
		{
			int n_out, nn_out, ii, nt = 0, all_fixed = 0;

			double *x_out = NULL, *xx_out = NULL, *param = NULL, log_norm_const = 0.0, log_prior = 0.0;
			inla_cgeneric_tp *def = NULL;
			def = (inla_cgeneric_tp *) mb->f_Qfunc_arg_orig[i];

			nt = def->ntheta;
			param = Calloc(nt, double);
			if (nt) {
				all_fixed = 1;
				for (ii = 0; ii < nt; ii++) {
					if (_NOT_FIXED(f_fixed[i][ii])) {
						param[ii] = theta[count];
						all_fixed = 0;
						count++;
					} else {
						param[ii] = mb->f_theta[i][ii][thread_id][0];
					}
				}
			}

			x_out = def->model_func(INLA_CGENERIC_LOG_NORM_CONST, param, def->data);
			if (def->debug) {
				inla_cgeneric_debug(stdout, def->secname, INLA_CGENERIC_LOG_NORM_CONST, x_out);
			}

			xx_out = def->model_func(INLA_CGENERIC_LOG_PRIOR, param, def->data);
			if (def->debug) {
				inla_cgeneric_debug(stdout, def->secname, INLA_CGENERIC_LOG_PRIOR, xx_out);
			}

			nn_out = (xx_out ? 1 : 0);
			switch (nn_out) {
			case 0:
			{
				log_prior = 0.0;
			}
				break;
			case 1:
			{
				// we need to add a check for 'all_fixed' here, as with control.mode=list(...,fixed=TRUE) will
				// trigger all_fixed=1.
				log_prior = (evaluate_hyper_prior && !all_fixed ? xx_out[0] : 0.0);
			}
				break;
			default:
				assert(0 == 1);
			}
			Free(xx_out);

			n_out = (x_out ? 1 : 0);
			switch (n_out) {
			case 0:
			{
				/*
				 * if it is the standard norm.const, the user can request us to compute it here if NULL is returned
				 */
				int *ilist = NULL, *jlist = NULL, n, len, k = 0, jj;
				double *Qijlist = NULL;
				GMRFLib_tabulate_Qfunc_tp *Qf = NULL;
				xx_out = def->model_func(INLA_CGENERIC_Q, param, def->data);
				if (def->debug) {
					inla_cgeneric_debug(stdout, def->secname, INLA_CGENERIC_Q, xx_out);
				}

				if ((int) xx_out[0] == -1) {
					// optimized output
					k = 1;
					len = (int) xx_out[k++];
					assert(len == def->len_list);
					assert(def->graph);
					assert(def->ilist);
					assert(def->jlist);
					n = def->graph->n;
					GMRFLib_tabulate_Qfunc_from_list2(&Qf, def->graph, def->len_list, def->ilist, def->jlist, &(xx_out[k]),
									  def->graph->n, NULL);
				} else {
					n = (int) xx_out[k++];
					len = (int) xx_out[k++];
					ilist = (int *) &xx_out[k];
					jlist = (int *) &xx_out[k + len];
					Qijlist = &xx_out[k + 2 * len];
					for (jj = 0; jj < len; jj++) {
						ilist[jj] = (int) xx_out[k + jj];
						jlist[jj] = (int) xx_out[k + len + jj];
					}
					GMRFLib_tabulate_Qfunc_from_list2(&Qf, def->graph, len, ilist, jlist, Qijlist, n, NULL);
				}

				int retval = GMRFLib_SUCCESS, ok = 0, num_try = 0, num_try_max = 100;
				GMRFLib_problem_tp *problem = NULL;
				GMRFLib_error_handler_tp *old_handler = GMRFLib_set_error_handler_off();
				double *cc_add = Calloc(n, double);

				if (mb->f_diag[i]) {
					for (jj = 0; jj < n; jj++) {
						cc_add[jj] = mb->f_diag[i];
					}
				}

				while (!ok) {
					retval = GMRFLib_init_problem_store(thread_id, &problem, NULL, NULL, cc_add, NULL,
									    def->graph, Qf->Qfunc, Qf->Qfunc_arg, mb->f_constr_orig[i], store);
					switch (retval) {
					case GMRFLib_EPOSDEF:
					{
						double eps = GSL_SQRT_DBL_EPSILON;
						for (jj = 0; jj < n; jj++) {
							cc_add[jj] = (cc_add[jj] == 0.0 ? eps : cc_add[jj] * 10.0);
						}
						break;
					}

					case GMRFLib_SUCCESS:
						ok = 1;
						break;

					default:
						/*
						 * some other error 
						 */
						GMRFLib_set_error_handler(old_handler);
						assert(0 == 1);
						abort();
					}

					if (++num_try >= num_try_max) {
						FIXME("This should not happen. Contact developers...");
						abort();
					}
				}
				Free(cc_add);
				GMRFLib_set_error_handler(old_handler);
				GMRFLib_evaluate(problem);
				log_norm_const = problem->sub_logdens;

				GMRFLib_free_problem(problem);
				GMRFLib_free_tabulate_Qfunc(Qf);
				Free(xx_out);
			}
				break;

			case 1:
			{
				log_norm_const = x_out[0];
			}
				break;

			default:
				assert(0 == 1);
			}

			if (debug) {
				for (ii = 0; ii < nt; ii++) {
					printf("p %.12g ", param[ii]);
				}
				printf(" %.12g  prior %.12g\n", log_norm_const, log_prior);
			}

			_SET_GROUP_RHO(nt);
			val += mb->f_nrep[i] * (normc_g + log_norm_const * (mb->f_ngroup[i] - grankdef)) + log_prior;

			Free(param);
			Free(x_out);
		}
			break;

		case F_AR1:
		{
			double mean_x;

			if (_NOT_FIXED(f_fixed[i][0])) {
				log_precision = theta[count];
				count++;
			} else {
				log_precision = mb->f_theta[i][0][thread_id][0];
			}

			if (_NOT_FIXED(f_fixed[i][1])) {
				phi_intern = theta[count];
				count++;
			} else {
				phi_intern = mb->f_theta[i][1][thread_id][0];
			}

			if (_NOT_FIXED(f_fixed[i][2])) {
				mean_x = theta[count];
				count++;
			} else {
				mean_x = mb->f_theta[i][2][thread_id][0];
			}

			phi = map_phi(phi_intern, MAP_FORWARD, NULL);
			_SET_GROUP_RHO(3);

			double log_precision_noise = log_precision - LOG_1mp(SQR(phi));

			if (mb->f_cyclic[i]) {
				logdet = inla_ar1_cyclic_logdet(N_orig, phi);
				val += mb->f_nrep[i] * (normc_g + gcorr * (LOG_NORMC_GAUSSIAN * (mb->f_N[i] - mb->f_rankdef[i])
									   + (mb->f_N[i] -
									      mb->f_rankdef[i]) / 2.0 * log_precision_noise +
									   ngroup * 0.5 * logdet));
			} else {
				val += mb->f_nrep[i] * (normc_g + gcorr * (LOG_NORMC_GAUSSIAN * (mb->f_N[i] - mb->f_rankdef[i])
									   + (mb->f_N[i] -
									      mb->f_rankdef[i]) / 2.0 * log_precision_noise +
									   ngroup * 0.5 * LOG_1mp(SQR(phi))));
			}
			if (_NOT_FIXED(f_fixed[i][0])) {
				val += PRIOR_EVAL(mb->f_prior[i][0], &log_precision);
			}
			if (_NOT_FIXED(f_fixed[i][1])) {
				val += PRIOR_EVAL(mb->f_prior[i][1], &phi_intern);
			}
			if (_NOT_FIXED(f_fixed[i][2])) {
				val += PRIOR_EVAL(mb->f_prior[i][2], &mean_x);
			}
		}
			break;

		case F_AR1C:
		{
			inla_ar1c_arg_tp *a = (inla_ar1c_arg_tp *) mb->f_Qfunc_arg_orig[i];

			if (_NOT_FIXED(f_fixed[i][0])) {
				log_precision = theta[count];
				count++;
			} else {
				log_precision = mb->f_theta[i][0][thread_id][0];
			}

			if (_NOT_FIXED(f_fixed[i][1])) {
				phi_intern = theta[count];
				count++;
			} else {
				phi_intern = mb->f_theta[i][1][thread_id][0];
			}

			phi = map_phi(phi_intern, MAP_FORWARD, NULL);
			_SET_GROUP_RHO(2);

			double log_precision_noise = log_precision - LOG_1mp(SQR(phi));
			val += mb->f_nrep[i] * (normc_g + gcorr * (LOG_NORMC_GAUSSIAN * (mb->f_N[i] - mb->f_rankdef[i])
								   + ((mb->f_N[i] - a->m * ngroup) -
								      mb->f_rankdef[i]) / 2.0 * log_precision_noise +
								   ngroup * 0.5 * LOG_1mp(SQR(phi))
								   + ngroup * 0.5 * a->logdet_Qbeta));
			if (_NOT_FIXED(f_fixed[i][0])) {
				val += PRIOR_EVAL(mb->f_prior[i][0], &log_precision);
			}
			if (_NOT_FIXED(f_fixed[i][1])) {
				val += PRIOR_EVAL(mb->f_prior[i][1], &phi_intern);
			}
		}
			break;

		case F_OU:
		{
			if (_NOT_FIXED(f_fixed[i][0])) {
				log_precision = theta[count];
				count++;
			} else {
				log_precision = mb->f_theta[i][0][thread_id][0];
			}
			if (_NOT_FIXED(f_fixed[i][1])) {
				phi_intern = theta[count];
				count++;
			} else {
				phi_intern = mb->f_theta[i][1][thread_id][0];
			}
			phi = map_exp(phi_intern, MAP_FORWARD, NULL);
			_SET_GROUP_RHO(2);

			int ii;
			double ou_nc = 0.0;
			int nn = ((inla_ou_arg_tp *) mb->f_Qfunc_arg_orig[i])->n;
			for (ii = 1; ii < nn; ii++) {
				ou_nc -= LOG_1mp(exp(-2.0 * phi * (mb->f_locations[i][ii] - mb->f_locations[i][ii - 1])));
			}
			val += mb->f_nrep[i] * (normc_g + gcorr * (LOG_NORMC_GAUSSIAN * (mb->f_N[i] - mb->f_rankdef[i])
								   + (mb->f_N[i] - mb->f_rankdef[i]) / 2.0 * log_precision + ngroup * ou_nc / 2.0));

			if (_NOT_FIXED(f_fixed[i][0])) {
				val += PRIOR_EVAL(mb->f_prior[i][0], &log_precision);
			}
			if (_NOT_FIXED(f_fixed[i][1])) {
				val += PRIOR_EVAL(mb->f_prior[i][1], &phi_intern);
			}
		}
			break;

		case F_BESAG2:
		{
			if (_NOT_FIXED(f_fixed[i][0])) {
				log_precision = theta[count];
				count++;
			} else {
				log_precision = mb->f_theta[i][0][thread_id][0];
			}
			if (_NOT_FIXED(f_fixed[i][1])) {
				// a_intern= log(a)
				a_intern = theta[count];
				count++;
			} else {
				a_intern = mb->f_theta[i][1][thread_id][0];
			}
			_SET_GROUP_RHO(2);
			// N is 2*graph->n here. 
			// Add the high precision contribution? 
			double high_precision = mb->f_precision[i];
			double n = mb->f_N[i] / 2.0;

			val += mb->f_nrep[i] * (normc_g + gcorr * (LOG_NORMC_GAUSSIAN * (2.0 * (n - mb->f_rankdef[i]))
								   + (n - mb->f_rankdef[i]) / 2.0 * log(high_precision)
								   + (n - mb->f_rankdef[i]) / 2.0 * (log_precision - 2.0 * a_intern)));
			if (_NOT_FIXED(f_fixed[i][0])) {
				val += PRIOR_EVAL(mb->f_prior[i][0], &log_precision);
			}
			if (_NOT_FIXED(f_fixed[i][1])) {
				val += PRIOR_EVAL(mb->f_prior[i][1], &a_intern);
			}
		}
			break;

		case F_BYM:
		{
			if (_NOT_FIXED(f_fixed[i][0])) {       /* iid */
				log_precision0 = theta[count];
				count++;
			} else {
				log_precision0 = mb->f_theta[i][0][thread_id][0];
			}
			if (_NOT_FIXED(f_fixed[i][1])) {       /* spatial */
				log_precision1 = theta[count];
				count++;
			} else {
				log_precision1 = mb->f_theta[i][1][thread_id][0];
			}
			_SET_GROUP_RHO(2);

			double n = (double) mb->f_n[i];
			val += mb->f_nrep[i] * (normc_g + gcorr * (LOG_NORMC_GAUSSIAN * (n / 2.0 + (n - mb->f_rankdef[i]) / 2.0)
								   + n / 2.0 * log_precision0	/* iid */
								   + (n - mb->f_rankdef[i]) / 2.0 * log_precision1));	/* spatial */
			if (_NOT_FIXED(f_fixed[i][0])) {
				val += PRIOR_EVAL(mb->f_prior[i][0], &log_precision0);
			}
			if (_NOT_FIXED(f_fixed[i][1])) {
				val += PRIOR_EVAL(mb->f_prior[i][1], &log_precision1);
			}
		}
			break;

		case F_RW2DIID:
		{
			if (_NOT_FIXED(f_fixed[i][0])) {
				log_precision = theta[count];
				count++;
			} else {
				log_precision = mb->f_theta[i][0][thread_id][0];
			}
			if (_NOT_FIXED(f_fixed[i][1])) {
				phi_intern = theta[count];
				count++;
			} else {
				phi_intern = mb->f_theta[i][1][thread_id][0];
			}
			_SET_GROUP_RHO(2);

			double n = (double) mb->f_n[i];
			phi = map_probability(phi_intern, MAP_FORWARD, NULL);

			val += mb->f_nrep[i] * (normc_g + gcorr * (LOG_NORMC_GAUSSIAN * n + n / 2.0 * (log_precision - LOG_1mp(phi))));

			if (_NOT_FIXED(f_fixed[i][0])) {
				val += PRIOR_EVAL(mb->f_prior[i][0], &log_precision);
			}
			if (_NOT_FIXED(f_fixed[i][1])) {
				val += PRIOR_EVAL(mb->f_prior[i][1], &phi_intern);
			}
		}
			break;

		case F_BYM2:
		{
			if (_NOT_FIXED(f_fixed[i][0])) {
				log_precision = theta[count];
				count++;
			} else {
				log_precision = mb->f_theta[i][0][thread_id][0];
			}
			if (_NOT_FIXED(f_fixed[i][1])) {
				phi_intern = theta[count];
				count++;
			} else {
				phi_intern = mb->f_theta[i][1][thread_id][0];
			}
			_SET_GROUP_RHO(2);

			double n = (double) mb->f_n[i];
			phi = map_probability(phi_intern, MAP_FORWARD, NULL);

			val += mb->f_nrep[i] * (normc_g + gcorr * (LOG_NORMC_GAUSSIAN * n + n / 2.0 * (log_precision - LOG_1mp(phi))));

			if (_NOT_FIXED(f_fixed[i][0])) {
				val += PRIOR_EVAL(mb->f_prior[i][0], &log_precision);
			}
			if (_NOT_FIXED(f_fixed[i][1])) {
				val += PRIOR_EVAL(mb->f_prior[i][1], &phi_intern);
			}
		}
			break;

		case F_2DIID:
		{
			if (mb->f_ngroup[i] > 1) {
				fprintf(stderr, "\n\n F_2DIID is not yet prepared for ngroup > 1\n");
				exit(EXIT_FAILURE);
			}

			assert(mb->f_ntheta[i] == 3);	       /* yes */
			if (_NOT_FIXED(f_fixed[i][0])) {
				log_precision0 = theta[count];
				count++;
			} else {
				log_precision0 = mb->f_theta[i][0][thread_id][0];
			}
			if (_NOT_FIXED(f_fixed[i][1])) {
				log_precision1 = theta[count];
				count++;
			} else {
				log_precision1 = mb->f_theta[i][1][thread_id][0];
			}
			if (_NOT_FIXED(f_fixed[i][2])) {
				rho_intern = theta[count];
				count++;
			} else {
				rho_intern = mb->f_theta[i][2][thread_id][0];
			}
			rho = map_rho(rho_intern, MAP_FORWARD, NULL);
			double n = (double) mb->f_n[i];
			assert(mb->f_ngroup[i] == 1);
			val += mb->f_nrep[i] * (LOG_NORMC_GAUSSIAN * 2.0 * (n - mb->f_rankdef[i])	/* yes, the total length is * N=2n */
						+(n - mb->f_rankdef[i]) / 2.0 * log_precision0	/* and there is n-pairs... */
						+ (n - mb->f_rankdef[i]) / 2.0 * log_precision1 - (n - mb->f_rankdef[i]) / 2.0 * LOG_1mp(SQR(rho)));
			if (_NOT_FIXED(f_fixed[i][0])) {
				val += PRIOR_EVAL(mb->f_prior[i][0], &log_precision0);
			}
			if (_NOT_FIXED(f_fixed[i][1])) {
				val += PRIOR_EVAL(mb->f_prior[i][1], &log_precision1);
			}
			if (_NOT_FIXED(f_fixed[i][2])) {
				val += PRIOR_EVAL(mb->f_prior[i][2], &rho_intern);
			}
		}
			break;

		case F_IID1D:
		case F_IID2D:
		case F_IID3D:
		case F_IID4D:
		case F_IID5D:
		{
			int jj, count_ref = count;
			int dim = WISHART_DIM(i);
			assert(dim > 0);

			int nt = inla_iid_wishart_nparam(dim);
			double log_jacobian = 0.0;
			double *theta_vec = Calloc(nt, double);
			int k = 0;
			nfixed = 0;
			for (j = 0; j < dim; j++) {
				if (_NOT_FIXED(f_fixed[i][k])) {
					theta_vec[k] = theta[count];
					count++;
				} else {
					nfixed++;
					theta_vec[k] = mb->f_theta[i][k][thread_id][0];
				}
				log_jacobian += log(map_precision(theta_vec[k], MAP_DFORWARD, NULL));
				theta_vec[k] = map_precision(theta_vec[k], MAP_FORWARD, NULL);
				k++;
			}

			for (j = 0; j < dim; j++) {
				for (jj = j + 1; jj < dim; jj++) {
					if (_NOT_FIXED(f_fixed[i][k])) {
						theta_vec[k] = theta[count];
						count++;
					} else {
						nfixed++;
						theta_vec[k] = mb->f_theta[i][k][thread_id][0];
					}
					log_jacobian += log(map_rho(theta_vec[k], MAP_DFORWARD, NULL));
					theta_vec[k] = map_rho(theta_vec[k], MAP_FORWARD, NULL);
					k++;
				}
			}
			assert(k == nt);

			fail = inla_iid_wishart_adjust(dim, theta_vec);
			Q = gsl_matrix_calloc(dim, dim);
			k = 0;
			for (j = 0; j < dim; j++) {
				gsl_matrix_set(Q, j, j, 1.0 / theta_vec[k]);
				k++;
			}
			for (j = 0; j < dim; j++) {
				for (jj = j + 1; jj < dim; jj++) {
					double value = theta_vec[k] / sqrt(theta_vec[j] * theta_vec[jj]);
					gsl_matrix_set(Q, j, jj, value);
					gsl_matrix_set(Q, jj, j, value);
					k++;
				}
			}
			assert(k == nt);

			GMRFLib_gsl_spd_inverse(Q);
			logdet = GMRFLib_gsl_spd_logdet(Q);
			gsl_matrix_free(Q);

			_SET_GROUP_RHO(nt);

			/*
			 * n is the small length. yes, the total length is N=dim*n
			 */
			double n = (double) (mb->f_n[i] / dim);	/* YES! */
			val += mb->f_nrep[i] * (normc_g + gcorr * (LOG_NORMC_GAUSSIAN * dim * (n - mb->f_rankdef[i])
								   + (n - mb->f_rankdef[i]) / 2.0 * logdet));
			if (fail) {
				val += PENALTY;
			}

			/*
			 * if all parameters are fixed, there is no prior to add
			 */
			if (count - count_ref > 0) {
				if (nfixed) {
					static char first = 1;
					if (first) {
						fprintf(stderr,
							"\n\n\nWARNING: Wishart prior is not corrected to account for %d fixed hyperparameters.\n\n",
							nfixed);
						first = 0;
					}
				}
				/*
				 * prior density wrt theta. Include here the Jacobian from going from (precision0, precision1, rho), to theta = (log_precision0,
				 * log_precision1, rho_intern). 
				 */
				val += PRIOR_EVAL(mb->f_prior[i][0], theta_vec) + log_jacobian;
			}
		}
			break;

		case F_IIDKD:
		{
			int dim = mb->f_order[i];
			assert(dim > 1);
			int nt = INLA_WISHARTK_NTHETA(dim);
			double *theta_vec = Calloc(nt, double);

			nfixed = 0;
			for (j = 0; j < nt; j++) {
				if (_NOT_FIXED(f_fixed[i][j])) {
					theta_vec[j] = theta[count];
					count++;
				} else {
					nfixed++;
					theta_vec[j] = mb->f_theta[i][j][thread_id][0];
				}
			}
			L = gsl_matrix_calloc(dim, dim);
			Q = gsl_matrix_calloc(dim, dim);
			inla_wishartk_build_Q(dim, theta_vec, Q, L);
			logdet = GMRFLib_gsl_spd_logdet(Q);
			gsl_matrix_free(Q);
			gsl_matrix_free(L);

			_SET_GROUP_RHO(nt);

			/*
			 * n is the small length. yes, the total length is N=dim*n
			 */
			double n = (double) (mb->f_n[i] / dim);	/* YES! */
			val += mb->f_nrep[i] * (normc_g + gcorr * (LOG_NORMC_GAUSSIAN * dim * (n - mb->f_rankdef[i])
								   + (n - mb->f_rankdef[i]) / 2.0 * logdet));
			val += PRIOR_EVAL(mb->f_prior[i][0], theta_vec);

			Free(theta_vec);
		}
			break;

		case F_INTSLOPE:
		{
			// this first part is just a copy from F_IID2D
			int jj, count_ref = count;
			int dim = 2;
			assert(dim > 0);

			int nt = inla_iid_wishart_nparam(dim);
			double log_jacobian = 0.0;
			double *theta_vec = Calloc(nt, double);
			int k = 0;
			nfixed = 0;
			for (j = 0; j < dim; j++) {
				if (_NOT_FIXED(f_fixed[i][k])) {
					theta_vec[k] = theta[count];
					count++;
				} else {
					nfixed++;
					theta_vec[k] = mb->f_theta[i][k][thread_id][0];
				}
				log_jacobian += log(map_precision(theta_vec[k], MAP_DFORWARD, NULL));
				theta_vec[k] = map_precision(theta_vec[k], MAP_FORWARD, NULL);
				k++;
			}

			for (j = 0; j < dim; j++) {
				for (jj = j + 1; jj < dim; jj++) {
					if (_NOT_FIXED(f_fixed[i][k])) {
						theta_vec[k] = theta[count];
						count++;
					} else {
						nfixed++;
						theta_vec[k] = mb->f_theta[i][k][thread_id][0];
					}
					log_jacobian += log(map_rho(theta_vec[k], MAP_DFORWARD, NULL));
					theta_vec[k] = map_rho(theta_vec[k], MAP_FORWARD, NULL);
					k++;
				}
			}
			assert(k == nt);

			fail = inla_iid_wishart_adjust(dim, theta_vec);	/* should not be needed */
			Q = gsl_matrix_calloc(dim, dim);
			k = 0;
			for (j = 0; j < dim; j++) {
				gsl_matrix_set(Q, j, j, 1.0 / theta_vec[k]);
				k++;
			}
			for (j = 0; j < dim; j++) {
				for (jj = j + 1; jj < dim; jj++) {
					double value = theta_vec[k] / sqrt(theta_vec[j] * theta_vec[jj]);
					gsl_matrix_set(Q, j, jj, value);
					gsl_matrix_set(Q, jj, j, value);
					k++;
				}
			}
			GMRFLib_gsl_spd_inverse(Q);
			assert(k == nt);

			if (count - count_ref > 0) {
				if (nfixed) {
					static char first = 1;
					if (first) {
						fprintf(stderr,
							"\n\n\nWARNING: Wishart prior is not corrected to account for %d fixed hyperparameters.\n\n",
							nfixed);
						first = 0;
					}
				}
				/*
				 * prior density wrt theta. Include here the Jacobian from going from (precision0, precision1,
				 * rho), to theta = (log_precision0, log_precision1, rho_intern).
				 */
				val += PRIOR_EVAL(mb->f_prior[i][0], theta_vec) + log_jacobian;
			}

			Free(theta_vec);

			for (j = 0; j < INTSLOPE_MAXTHETA; j++) {
				if (_NOT_FIXED(f_fixed[i][k + j])) {
					double gam = theta[count];
					count++;
					nt++;
					val += PRIOR_EVAL(mb->f_prior[i][k + j], &gam);
				}
			}
			_SET_GROUP_RHO(mb->f_ntheta[i]);

			inla_intslope_arg_tp *arg = (inla_intslope_arg_tp *) mb->f_Qfunc_arg_orig[i];
			ngroup = mb->f_ngroup[i];

			assert(mb->f_rankdef[i] == 0);	       /* as this does not make sense if its not 0 */

			logdet = GMRFLib_gsl_spd_logdet(Q) * arg->nsubject + log(arg->precision) * arg->n;
			val += mb->f_nrep[i] * (normc_g + gcorr * ngroup * (LOG_NORMC_GAUSSIAN * arg->N + 0.5 * logdet));
			gsl_matrix_free(Q);
			if (fail) {
				val += PENALTY;
			}
		}
			break;

		case F_MATERN2D:
		{
			/*
			 * this is the safe version, which also works correctly for diagonal > 0. It is possible to reuse the calculations for the same
			 * range and different precision, provided diagonal = 0, but care must be taken about the constraints.
			 */

			typedef struct {
				int n;
				int N;
				int ngroup;
				int nrep;
				double precision;
				double range;
				double *c;
				double rankdef1;
				GMRFLib_matern2ddef_tp *matern2ddef;
				GMRFLib_problem_tp *problem;
			} Hold_tp;

			int jj;
			Hold_tp *h = NULL;
			GMRFLib_matern2ddef_tp *q = NULL;

			h = Calloc(1, Hold_tp);

			h->nrep = mb->f_nrep[i];
			h->ngroup = mb->f_ngroup[i];
			h->n = mb->f_n[i] / h->ngroup;
			h->N = mb->f_N[i] / h->ngroup;

			assert(h->N == mb->f_graph_orig[i]->n);

			if (debug) {
				P(h->n);
				P(h->N);
				P(h->nrep);
				P(h->ngroup);
			}

			if (mb->f_diag[i]) {
				h->c = Calloc(h->N, double);
				for (jj = 0; jj < h->N; jj++) {
					h->c[jj] = mb->f_diag[i];
				}
			} else {
				h->c = NULL;
			}

			q = (GMRFLib_matern2ddef_tp *) mb->f_Qfunc_arg_orig[i];
			h->matern2ddef = Calloc(1, GMRFLib_matern2ddef_tp);
			Memcpy(h->matern2ddef, q, sizeof(GMRFLib_matern2ddef_tp));
			if (_NOT_FIXED(f_fixed[i][0])) {
				h->precision = map_precision(theta[count], MAP_FORWARD, NULL);
				val += PRIOR_EVAL(mb->f_prior[i][0], &theta[count]);
				count++;
			} else {
				h->precision = map_precision(mb->f_theta[i][0][thread_id][0], MAP_FORWARD, NULL);
			}
			if (_NOT_FIXED(f_fixed[i][1])) {
				h->range = map_range(theta[count], MAP_FORWARD, NULL);
				val += PRIOR_EVAL(mb->f_prior[i][1], &theta[count]);
				count++;
			} else {
				h->range = map_range(mb->f_theta[i][1][thread_id][0], MAP_FORWARD, NULL);
			}
			HYPER_NEW(h->matern2ddef->log_prec_omp, map_precision(h->precision, MAP_BACKWARD, NULL));
			HYPER_NEW(h->matern2ddef->log_range_omp, map_precision(h->range, MAP_BACKWARD, NULL));

			_SET_GROUP_RHO(2);

			GMRFLib_init_problem_store(thread_id, &(h->problem), NULL, NULL, h->c, NULL, mb->f_graph_orig[i], mb->f_Qfunc_orig[i],
						   (void *) h->matern2ddef, mb->f_constr_orig[i], store);
			if (debug) {
				P(h->precision);
				P(h->range);
				P(h->problem->sub_logdens);
			}
			val += h->nrep * (h->problem->sub_logdens * (ngroup - grankdef) + normc_g);

			Free(h->c);
			HYPER_FREE(h->matern2ddef->log_prec_omp);
			HYPER_FREE(h->matern2ddef->log_range_omp);
			Free(h->matern2ddef);
			GMRFLib_free_problem(h->problem);
			Free(h);
		}
			break;

		case F_DMATERN:
		{
			double log_range, log_nu, range, var, nu, prec;

			if (_NOT_FIXED(f_fixed[i][0])) {
				log_precision = theta[count];
				count++;
			} else {
				log_precision = mb->f_theta[i][0][thread_id][0];
			}
			if (_NOT_FIXED(f_fixed[i][1])) {
				log_range = theta[count];
				count++;
			} else {
				log_range = mb->f_theta[i][1][thread_id][0];
			}
			if (_NOT_FIXED(f_fixed[i][2])) {
				log_nu = theta[count];
				count++;
			} else {
				log_nu = mb->f_theta[i][2][thread_id][0];
			}

			dmatern_arg_tp *a = (dmatern_arg_tp *) (mb->f_Qfunc_arg_orig[i]);
			gsl_matrix *S = gsl_matrix_calloc(a->n, a->n);
			prec = map_exp(log_precision, MAP_FORWARD, NULL);
			range = map_range(log_range, MAP_FORWARD, NULL);
			nu = map_exp(log_nu, MAP_FORWARD, NULL);
			var = 1.0 / prec;

			for (int ii = 0; ii < a->n; ii++) {
				for (int jj = ii; jj < a->n; jj++) {
					double dist;
					dist = gsl_matrix_get(a->dist, ii, jj);
					val = var * inla_dmatern_cf(dist, range, nu);
					gsl_matrix_set(S, ii, jj, val);
					gsl_matrix_set(S, jj, ii, val);
				}
			}
			// need a '-' as we're computing the |S| instead of |Q|.
			logdet = -GMRFLib_gsl_spd_logdet(S) / a->n;	/* logdet(Q), as if a->n=1. makes its easier below */
			gsl_matrix_free(S);

			_SET_GROUP_RHO(3);
			val += mb->f_nrep[i] * (normc_g + gcorr * (LOG_NORMC_GAUSSIAN * (mb->f_N[i] - mb->f_rankdef[i]) +
								   (mb->f_N[i] - mb->f_rankdef[i]) / 2.0 * logdet));

			if (_NOT_FIXED(f_fixed[i][0])) {
				val += PRIOR_EVAL(mb->f_prior[i][0], &log_precision);
			}
			if (_NOT_FIXED(f_fixed[i][1])) {
				val += PRIOR_EVAL(mb->f_prior[i][1], &log_range);
			}
			if (_NOT_FIXED(f_fixed[i][2])) {
				val += PRIOR_EVAL(mb->f_prior[i][2], &log_nu);
			}
		}
			break;

		case F_BESAGPROPER:
		{
			typedef struct {
				int n;
				int N;
				int ngroup;
				int nrep;
				double **log_prec;
				double **log_diag;
				double *c;
				double rankdef1;
				inla_besag_proper_Qfunc_arg_tp *def;
				GMRFLib_problem_tp *problem;
			} Hold_tp;
			static Hold_tp ***hhold = NULL;

			if (!hhold) {
#pragma omp critical (Name_35784cb53aa98d636cf2d0897410586e2705f61e)
				{
					if (!hhold) {
						hhold = Calloc(GMRFLib_CACHE_LEN(), Hold_tp **);
					}
				}
			}
			int idx = 0;
			GMRFLib_CACHE_SET_ID(idx);

			int jj;
			Hold_tp *h = NULL, **hold = NULL;

			hold = hhold[idx];
			if (!hold) {
				hold = Calloc(mb->nf, Hold_tp *);
			}

			if (!hold[i]) {
				h = hold[i] = Calloc(1, Hold_tp);

				h->nrep = mb->f_nrep[i];
				h->ngroup = mb->f_ngroup[i];
				h->n = mb->f_n[i] / h->ngroup;
				h->N = mb->f_N[i] / h->ngroup;

				assert(h->N == mb->f_graph_orig[i]->n);

				if (debug) {
					P(h->n);
					P(h->N);
					P(h->nrep);
					P(h->ngroup);
				}

				HYPER_NEW(h->log_prec, 0.0);
				HYPER_NEW(h->log_diag, 0.0);

				if (mb->f_diag[i]) {
					h->c = Calloc(h->N, double);
					for (jj = 0; jj < h->N; jj++) {
						h->c[jj] = mb->f_diag[i];
					}
				}

				h->def = Calloc(1, inla_besag_proper_Qfunc_arg_tp);
				Memcpy(h->def, mb->f_Qfunc_arg_orig[i], sizeof(inla_besag_proper_Qfunc_arg_tp));
				h->def->log_prec = h->log_prec;
				h->def->log_diag = h->log_diag;
			} else {
				h = hold[i];
			}

			if (_NOT_FIXED(f_fixed[i][0])) {
				h->log_prec[thread_id][0] = theta[count];
				val += PRIOR_EVAL(mb->f_prior[i][0], &theta[count]);
				count++;
			} else {
				h->log_prec[thread_id][0] = mb->f_theta[i][0][thread_id][0];
			}
			if (_NOT_FIXED(f_fixed[i][1])) {
				h->log_diag[thread_id][0] = theta[count];
				val += PRIOR_EVAL(mb->f_prior[i][1], &theta[count]);
				count++;
			} else {
				h->log_diag[thread_id][0] = mb->f_theta[i][1][thread_id][0];
			}

			_SET_GROUP_RHO(2);

			GMRFLib_init_problem_store(thread_id, &(h->problem), NULL, NULL, h->c, NULL, mb->f_graph_orig[i], mb->f_Qfunc_orig[i],
						   (void *) h->def, mb->f_constr_orig[i], store);
			if (debug) {
				P(h->log_prec[thread_id][0]);
				P(h->log_diag[thread_id][0]);
				P(h->problem->sub_logdens);
			}
			val += h->nrep * (h->problem->sub_logdens * (ngroup - grankdef) + normc_g);

			GMRFLib_free_problem(h->problem);
			h->problem = NULL;
		}
			break;

		case F_BESAGPROPER2:
		{
			typedef struct {
				int n;
				int N;
				int ngroup;
				int nrep;
				double **log_prec;
				double **logit_lambda;
				double *c;
				double rankdef1;
				inla_besag_proper2_Qfunc_arg_tp *def;
				GMRFLib_problem_tp *problem;
			} Hold_tp;
			static Hold_tp ***hhold = NULL;

			if (!hhold) {
#pragma omp critical (Name_7acab2f371bbea723e9820a667f70647967dbd17)
				if (!hhold) {
					hhold = Calloc(GMRFLib_CACHE_LEN(), Hold_tp **);
				}
			}
			int idx = 0;
			GMRFLib_CACHE_SET_ID(idx);

			int jj;
			Hold_tp *h = NULL, **hold = NULL;

			hold = hhold[idx];
			if (!hold) {
				hold = Calloc(mb->nf, Hold_tp *);
			}

			if (!hold[i]) {
				h = hold[i] = Calloc(1, Hold_tp);
				h->nrep = mb->f_nrep[i];
				h->ngroup = mb->f_ngroup[i];
				h->n = mb->f_n[i] / h->ngroup;
				h->N = mb->f_N[i] / h->ngroup;

				assert(h->N == mb->f_graph_orig[i]->n);

				if (debug) {
					P(h->n);
					P(h->N);
					P(h->nrep);
					P(h->ngroup);
				}

				HYPER_NEW(h->log_prec, 0.0);
				HYPER_NEW(h->logit_lambda, 0.0);

				if (mb->f_diag[i]) {
					h->c = Calloc(h->N, double);
					for (jj = 0; jj < h->N; jj++) {
						h->c[jj] = mb->f_diag[i];
					}
				}

				h->def = Calloc(1, inla_besag_proper2_Qfunc_arg_tp);
				Memcpy(h->def, mb->f_Qfunc_arg_orig[i], sizeof(inla_besag_proper2_Qfunc_arg_tp));
				h->def->log_prec = h->log_prec;
				h->def->logit_lambda = h->logit_lambda;
			} else {
				h = hold[i];
			}

			if (_NOT_FIXED(f_fixed[i][0])) {
				h->log_prec[thread_id][0] = theta[count];
				val += PRIOR_EVAL(mb->f_prior[i][0], &theta[count]);
				count++;
			} else {
				h->log_prec[thread_id][0] = mb->f_theta[i][0][thread_id][0];
			}
			if (_NOT_FIXED(f_fixed[i][1])) {
				h->logit_lambda[thread_id][0] = theta[count];
				val += PRIOR_EVAL(mb->f_prior[i][1], &theta[count]);
				count++;
			} else {
				h->logit_lambda[thread_id][0] = mb->f_theta[i][1][thread_id][0];
			}

			_SET_GROUP_RHO(2);

			GMRFLib_init_problem_store(thread_id, &(h->problem), NULL, NULL, h->c, NULL, mb->f_graph_orig[i], mb->f_Qfunc_orig[i],
						   (void *) h->def, mb->f_constr_orig[i], store);

			val += h->nrep * (h->problem->sub_logdens * (ngroup - grankdef) + normc_g);

			GMRFLib_free_problem(h->problem);
			h->problem = NULL;
		}
			break;

		case F_COPY:
		{
			if (_NOT_FIXED(f_fixed[i][0]) && !mb->f_same_as[i]) {
				beta = theta[count];
				count++;
				val += PRIOR_EVAL(mb->f_prior[i][0], &beta);
			}
			val += mb->f_Ntotal[i] * (LOG_NORMC_GAUSSIAN + 0.5 * log(mb->f_precision[i]));
		}
			break;

		case F_SCOPY:
		{
			inla_scopy_arg_tp *a = (inla_scopy_arg_tp *) mb->f_Qfunc_arg_orig[i];
			for (int k = 0; k < 2; k++) {	       /* mean and slope */
				if (_NOT_FIXED(f_fixed[i][k])) {
					double b = theta[count];
					count++;
					val += PRIOR_EVAL(mb->f_prior[i][k], &b);
				}
			}
			for (int k = 2; k < a->nbeta; k++) {   /* mean and slope */
				if (_NOT_FIXED(f_fixed[i][k])) {
					double b = theta[count];
					count++;
					// yes, we're using the prior for beta[2]
					val += PRIOR_EVAL(mb->f_prior[i][2], &b);
				}
			}
		}
			break;

		case F_CLINEAR:
		{
			if (_NOT_FIXED(f_fixed[i][0])) {
				beta = theta[count];
				count++;
				val += PRIOR_EVAL(mb->f_prior[i][0], &beta);
			}
			val += mb->f_Ntotal[i] * (LOG_NORMC_GAUSSIAN + 0.5 * log(mb->f_precision[i]));
		}
			break;

		case F_SIGM:
		case F_REVSIGM:
		{
			if (_NOT_FIXED(f_fixed[i][0])) {
				beta = theta[count];
				count++;
				val += PRIOR_EVAL(mb->f_prior[i][0], &beta);
			}
			if (_NOT_FIXED(f_fixed[i][1])) {
				log_halflife = theta[count];
				count++;
				val += PRIOR_EVAL(mb->f_prior[i][1], &log_halflife);
			}
			if (_NOT_FIXED(f_fixed[i][2])) {
				log_shape = theta[count];
				count++;
				val += PRIOR_EVAL(mb->f_prior[i][2], &log_shape);
			}
			val += mb->f_Ntotal[i] * (LOG_NORMC_GAUSSIAN + 0.5 * log(mb->f_precision[i]));
		}
			break;

		case F_LOG1EXP:
		{
			if (_NOT_FIXED(f_fixed[i][0])) {
				beta = theta[count];
				count++;
				val += PRIOR_EVAL(mb->f_prior[i][0], &beta);
			}
			if (_NOT_FIXED(f_fixed[i][1])) {
				alpha = theta[count];
				count++;
				val += PRIOR_EVAL(mb->f_prior[i][1], &alpha);
			}
			if (_NOT_FIXED(f_fixed[i][2])) {
				gama = theta[count];
				count++;
				val += PRIOR_EVAL(mb->f_prior[i][2], &gama);
			}
			val += mb->f_Ntotal[i] * (LOG_NORMC_GAUSSIAN + 0.5 * log(mb->f_precision[i]));
		}
			break;

		case F_LOGDIST:
		{
			if (_NOT_FIXED(f_fixed[i][0])) {
				beta = theta[count];
				count++;
				val += PRIOR_EVAL(mb->f_prior[i][0], &beta);
			}
			if (_NOT_FIXED(f_fixed[i][1])) {
				alpha1 = theta[count];
				count++;
				val += PRIOR_EVAL(mb->f_prior[i][1], &alpha1);
			}
			if (_NOT_FIXED(f_fixed[i][2])) {
				alpha2 = theta[count];
				count++;
				val += PRIOR_EVAL(mb->f_prior[i][2], &alpha2);
			}
			val += mb->f_Ntotal[i] * (LOG_NORMC_GAUSSIAN + 0.5 * log(mb->f_precision[i]));
		}
			break;

		default:
		{
			P(mb->f_id[i]);
			abort();
			GMRFLib_ASSERT(0 == 1, GMRFLib_ESNH);
		}
			break;
		}
	}

	if (mb->data_sections[0].lp_scale) {
		for (int k = 0; k < INLA_LP_SCALE_MAX; k++) {
			if (mb->data_sections[0].lp_scale_in_use[k]) {
				if (_NOT_FIXED(data_sections[0].lp_scale_nfixed[k])) {
					beta = theta[count];
					val = PRIOR_EVAL(mb->data_sections[0].lp_scale_nprior[k], &beta);
					count++;
				}
			}
		}
	}

	if (debug) {
		P(val);
		P(count);
		P(mb->ntheta);
		P(ntheta);
	}

	if (!(mb->mode_fixed)) {
		assert((count == mb->ntheta) && (count == ntheta));	/* check... */
	}
#undef _SET_GROUP_RHO
#undef _NOT_FIXED
	return val;
}

double inla_compute_initial_value(int idx, GMRFLib_logl_tp *loglfunc, double *x_vec, void *arg)
{
	/*
	 * solve arg min logl(x[i]) - prec * 0.5*(x[i]-mean)^2. But we have no option of what PREC is, so I set it to 10.0
	 */
	Data_section_tp *ds = (Data_section_tp *) arg;
	double prec, prec_max = 1.0E6, prec_min = 10.0, w, x, xnew, f, deriv, dderiv, arr[3], eps = 1.0E-4, steplen = 1.0E-4, mean = -OFFSET(idx);
	int niter = 0, niter_min = 25, niter_max = 100, stencil = 3;
	const int debug = 0;

	int thread_id = 0;
	x = xnew = mean;

	while (1) {
		w = (double) DMIN(niter_min, niter) / (double) niter_min;
		prec = exp(w * log(prec_min) + (1.0 - w) * log(prec_max));
		GMRFLib_2order_taylor(thread_id, &arr[0], &arr[1], &arr[2], NULL, 1.0, x, idx, x_vec, loglfunc, arg, &steplen, &stencil);
		f = arr[0] - 0.5 * prec * SQR((x - mean));
		deriv = arr[1] - prec * (x - mean);
		dderiv = DMIN(0.0, arr[2]) - prec;
		xnew = x - DMIN(0.25 + niter * 0.25, 1.0) * deriv / dderiv;
		if (debug) {
			printf("idx %d x %.4g xnew %.4g f %.4g deriv %.4g dderiv %.4g mean %.6g prec %.2g\n",
			       idx, x, xnew, f, deriv, dderiv, mean, prec);
		}
		x = xnew;

		if (niter > niter_min && ABS(deriv / dderiv) < eps) {
			break;
		}
		if (++niter > niter_max) {
			x = mean;
			break;
		}
	}

	return x;
}

int inla_INLA_preopt_experimental(inla_tp *mb)
{
	double *c = NULL, *x = NULL, *b = NULL;
	int N, count, i, j;
	GMRFLib_bfunc_tp **bfunc = NULL;
	GMRFLib_preopt_tp *preopt = NULL;

	if (mb->verbose) {
		printf("%s...\n", __GMRFLib_FuncName);
	}

	/*
	 * We need to determine the strategy if strategy is default 
	 */
	int ntot = 0;

	ntot = mb->nlinear;
	for (i = 0; i < mb->nf; i++) {
		ntot += mb->f_graph[i]->n;
	}
	N = ntot;
	Calloc_init(3 * N, 3);

	ntot += mb->predictor_m + mb->predictor_n;
	if (mb->strategy == GMRFLib_OPENMP_STRATEGY_DEFAULT) {
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

	GMRFLib_density_storage_strategy = (mb->strategy == GMRFLib_OPENMP_STRATEGY_HUGE || mb->strategy == GMRFLib_OPENMP_STRATEGY_LARGE ?
					    GMRFLib_DENSITY_STORAGE_STRATEGY_LOW : GMRFLib_DENSITY_STORAGE_STRATEGY_HIGH);
	GMRFLib_openmp->strategy = mb->strategy;

	b = Calloc_get(N);
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
	for (i = 0; i < mb->nlinear; i++) {
		b[count] = mb->linear_precision[i] * mb->linear_mean[i];
		count++;
	}
	assert(count == N);

	GMRFLib_prior_mean_tp **prior_mean = Calloc(N, GMRFLib_prior_mean_tp *);
	for (count = 0, i = 0; i < mb->nf; i++) {
		if (mb->f_bfunc2[i]) {
			for (j = 0; j < mb->f_Ntotal[i]; j++) {
				prior_mean[count + j] = Calloc(1, GMRFLib_prior_mean_tp);
				prior_mean[count + j]->bdef = mb->f_bfunc2[i];
				prior_mean[count + j]->idx = j;
				prior_mean[count + j]->fixed_mean = 0.0;
			}
		}
		count += mb->f_Ntotal[i];
	}
	for (i = 0; i < mb->nlinear; i++) {
		prior_mean[count] = Calloc(1, GMRFLib_prior_mean_tp);
		prior_mean[count]->bdef = NULL;
		prior_mean[count]->idx = -1;
		prior_mean[count]->fixed_mean = mb->linear_mean[i];
		count++;
	}
	assert(count == N);

	// VB corrections
	if (mb->ai_par->vb_enable) {
		// tp = 0 is mean, tp = 1 is variance
		for (int tp = 0; tp < 2; tp++) {
			char *vb_nodes = Calloc(N, char);
			int debug = 0;
			int local_count = 0;
			count = 0;

			for (i = 0; i < mb->nf; i++) {

				int n = mb->f_N[i] / mb->f_ngroup[i];
				int ngroup = mb->f_ngroup[i];
				int nrep = mb->f_Ntotal[i] / mb->f_N[i];
				ntot = mb->f_Ntotal[i];
				int lim = (tp == 0 ? mb->ai_par->vb_f_enable_limit_mean : mb->ai_par->vb_f_enable_limit_variance);
				int nngroup = n * ngroup;
				assert(ntot == n * ngroup * nrep);

				if (lim > 0) {
					// yes, integer division here && we correct also for mb->nf
					if (tp == 0) {
						lim = IMAX(0, IMIN(lim, mb->ai_par->vb_f_enable_limit_mean_max / (mb->nf * ngroup * nrep)));
					} else {
						lim = IMAX(0, IMIN(lim, mb->ai_par->vb_f_enable_limit_variance_max / (mb->nf * ngroup * nrep)));
					}
				}

				if (debug) {
					P(n);
					P(ngroup);
					P(nrep);
					P(ntot);
					P(lim);
				}

				if (lim > 0) {
					GMRFLib_idx_tp *vb = mb->f_vb_correct[i];

					if (debug) {
						P(vb->idx[0]);
					}

					if ((vb->idx[0] == -1L && n <= lim)) {
						for (j = 0; j < ntot; j++) {
							vb_nodes[count + j] = (char) 1;
							local_count++;
						}
					} else if (vb->idx[0] == -1L) {
						int len = IMAX(1, n / lim);
						int k = IMAX(1, len / 2);
						for (int r = 0; r < nrep; r++) {
							for (int g = 0; g < ngroup; g++) {
								for (j = 0; j < lim; j++) {
									int jj = (j * len + k) % n + g * n + r * nngroup;
									if (debug)
										printf("%d %d %d %d\n", g, r, j, jj);
									vb_nodes[count + jj] = (char) 1;
									local_count++;
								}
							}
						}
					} else if (vb->idx[0] >= 0) {
						for (int r = 0; r < nrep; r++) {
							for (int g = 0; g < ngroup; g++) {
								for (j = 0; j < vb->n; j++) {
									if (LEGAL(vb->idx[j], n)) {
										int jj = vb->idx[j] + g * n + r * nngroup;
										if (debug)
											printf("%d %d %d %d\n", g, r, j, jj);
										vb_nodes[count + jj] = (char) 1;
										local_count++;
									}
								}
							}
						}
					}
				}
				count += mb->f_Ntotal[i];
			}
			for (i = 0; i < mb->nlinear; i++) {
				vb_nodes[count++] = (char) 1;
				local_count++;
			}
			if (local_count == 0) {		       /* then there is nothting to correct for */
				Free(vb_nodes);
				vb_nodes = NULL;
			}

			if (tp == 0) {
				mb->ai_par->vb_nodes_mean = vb_nodes;
			} else if (tp == 1) {
				mb->ai_par->vb_nodes_variance = vb_nodes;
			} else {
				assert(0 == 1);
			}
		}
	}

	double tref = GMRFLib_timer();
	GMRFLib_openmp_implement_strategy(GMRFLib_OPENMP_PLACES_GCPO_BUILD, NULL, NULL);
	GMRFLib_preopt_init(&preopt,
			    mb->predictor_n, mb->nf, mb->f_c, mb->f_weights,
			    mb->f_graph, mb->f_Qfunc, mb->f_Qfunc_arg, mb->f_sumzero, mb->f_constr,
			    mb->f_diag,
			    mb->ff_Qfunc, mb->ff_Qfunc_arg, mb->f_Alocal,
			    mb->nlinear, mb->linear_covariate, mb->linear_precision, bfunc, mb->ai_par, mb->predictor_A_fnm, mb->global_constr);
	mb->preopt = preopt;
	assert(preopt->latent_graph->n == N);

	// time the two versions of Qfunc_like
	double time_used_Qx[2] = { 0.0, 0.0 };
	double time_used_pred[2] = { 0.0, 0.0 };

	if (GMRFLib_internal_opt) {
		// cannot run this in parallel as we're changing global variables
		GMRFLib_openmp_implement_strategy(GMRFLib_OPENMP_PLACES_TIMING, NULL, NULL);
		int thread_id = 0;
		int nn = preopt->preopt_graph->n;
		assert(omp_get_thread_num() == 0);
		double res[4] = { 0, 0, 0, 0 };
		double *test_vector = Calloc(nn, double);
		for (i = 0; i < nn; i++) {
			test_vector[i] = GMRFLib_uniform();
		}
		for (int time = -2; time < 4; time++) {
			for (int mett = 0; mett < 2; mett++) {
				GMRFLib_Qx_strategy = mett;
				double *cpu = GMRFLib_preopt_measure_time(thread_id, preopt, res + mett * 2, test_vector);
				if (time > 0) {
					time_used_Qx[mett] += cpu[1];
				}
				// printf("%d %d %f %f\n", time, mett, cpu[0], cpu[1]);
				Free(cpu);
			}
		}
		Free(test_vector);
		double eps = sqrt(nn) * 1.0e-6;
		if (!(ABS(res[0] - res[2]) < DMAX(1.0, ABS(res[0])) * eps) || !(ABS(res[1] - res[3]) < DMAX(1.0, ABS(res[1])) * eps)) {
			P(res[0]);
			P(res[2]);
			P(res[1]);
			P(res[3]);
			assert(ABS(res[0] - res[2]) < DMAX(1, ABS(res[0])) * eps);
			assert(ABS(res[1] - res[3]) < DMAX(1, ABS(res[1])) * eps);
		}

		// we have a slight preference for the simpler/serial ones
		GMRFLib_Qx_strategy = (time_used_Qx[0] / time_used_Qx[1] < 1.1 ? 0 : 1);

		// do this alone as this strategy depends on the previous choices
		for (int time = -2; time < 4; time++) {
			for (int mettt = 0; mettt < 2; mettt++) {
				GMRFLib_preopt_predictor_strategy = mettt;
				double *cpu = GMRFLib_preopt_measure_time2(preopt);
				if (time > 0) {
					time_used_pred[mettt] += cpu[0];
				}
				// printf("%d %f\n", mettt, cpu[0]);
				Free(cpu);
			}
		}
		// we have a slight preference for the simpler/serial ones
		GMRFLib_preopt_predictor_strategy = (time_used_pred[0] / time_used_pred[1] < 1.1 ? 0 : 1);
	} else {
		GMRFLib_Qx_strategy = 0;
		GMRFLib_preopt_predictor_strategy = 0;
		GMRFLib_sort2_id_cut_off = 128;		       // override value found 
		GMRFLib_sort2_dd_cut_off = 128;		       // override value found 
	}

#if 0
	// report timings
	double time_loop[13] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
	if (GMRFLib_internal_opt && GMRFLib_dot_product_optim_report) {
		for (i = 0; i < GMRFLib_CACHE_LEN(); i++) {
			for (j = 0; j < 13; j++) {
				time_loop[j] += GMRFLib_dot_product_optim_report[i][j];
			}
		}
		double time_sum = GMRFLib_dsum(6, time_loop);
		if (time_sum > 0.0) {
			time_sum = 1.0 / time_sum;
			GMRFLib_dscale(6, time_sum, time_loop);
			time_loop[6] *= time_sum;
		}
		time_sum = GMRFLib_dsum(6, time_loop + 7);
		if (time_sum > 0.0) {
			time_sum = 1.0 / time_sum;
			GMRFLib_dscale(6, time_sum, time_loop + 7);
		}
	}
#endif
#if !defined(INLA_WITH_MKL) && !defined(INLA_WITH_ARMPL)
	// report timings
	double time_loop[5] = { 0.0, 0.0, 0.0, 0.0, 0.0 };
	if (GMRFLib_internal_opt && GMRFLib_dot_product_optim_report) {
		for (i = 0; i < GMRFLib_CACHE_LEN(); i++) {
			for (j = 0; j < 5; j++) {
				time_loop[j] += GMRFLib_dot_product_optim_report[i][j];
			}
		}
		double time_sum = GMRFLib_dsum(2, time_loop);
		if (time_sum > 0.0) {
			time_sum = 1.0 / time_sum;
			GMRFLib_dscale(2, time_sum, time_loop);
			time_loop[2] *= time_sum;
		}
		time_sum = GMRFLib_dsum(2, time_loop + 3);
		if (time_sum > 0.0) {
			time_sum = 1.0 / time_sum;
			GMRFLib_dscale(2, time_sum, time_loop + 3);
		}
	}
#endif

	GMRFLib_openmp_implement_strategy(GMRFLib_OPENMP_PLACES_DEFAULT, NULL, NULL);
	if (mb->verbose) {
		printf("\tMode....................... [%s]\n", GMRFLib_MODE_NAME());
		printf("\tSetup...................... [%.2fs]\n", GMRFLib_timer() - tref);
		printf("\tSparse-matrix library...... [%s]\n", mb->smtp);
		if (!strcasecmp("taucs", mb->smtp)) {
			printf("\tsort L..................... [%s]\n", (GMRFLib_taucs_sort_L ? "yes" : "no"));
		}
		printf("\tOpenMP strategy............ [%s]\n", GMRFLib_OPENMP_STRATEGY_NAME(GMRFLib_openmp->strategy));
		printf("\tnum.threads................ [%1d:%1d]\n", GMRFLib_openmp->max_threads_nested[0], GMRFLib_openmp->max_threads_nested[1]);
		if (GMRFLib_openmp->adaptive) {
			printf("\tnum.threads (adaptive)..... [%1d]\n", GMRFLib_PARDISO_MAX_NUM_THREADS());
		}
		if (GMRFLib_openmp->blas_num_threads_force) {
			printf("\tblas.num.threads........... [%1d]\n", GMRFLib_openmp->blas_num_threads_force);
		} else {
			printf("\tblas.num.threads........... [%s]\n", "adaptive");
		}
		printf("\tDensity-strategy........... [%s]\n",
		       (GMRFLib_density_storage_strategy == GMRFLib_DENSITY_STORAGE_STRATEGY_LOW ? "Low" : "High"));
		printf("\tSize of graph.............. [%d]\n", N);
		printf("\tNumber of constraints...... [%d]\n", (preopt->latent_constr ? preopt->latent_constr->nc : 0));
		if (GMRFLib_internal_opt) {
			printf("\tOptimizing sort2_id........ [%1d]\n", GMRFLib_sort2_id_cut_off);
			printf("\tOptimizing sort2_dd........ [%1d]\n", GMRFLib_sort2_dd_cut_off);
			printf("\tOptimizing Qx-strategy..... serial[%.3f] parallel [%.3f] choose[%s]\n",
			       time_used_Qx[0] / (time_used_Qx[0] + time_used_Qx[1]), time_used_Qx[1] / (time_used_Qx[0] + time_used_Qx[1]),
			       (GMRFLib_Qx_strategy == 0 ? "serial" : "parallel"));
			printf("\tOptimizing pred-strategy... plain [%.3f] data-rich[%.3f] choose[%s]\n",
			       time_used_pred[0] / (time_used_pred[0] + time_used_pred[1]),
			       time_used_pred[1] / (time_used_pred[0] + time_used_pred[1]),
			       (GMRFLib_preopt_predictor_strategy == 0 ? "plain" : "data-rich"));
#if 0
			printf("\tOptimizing dot-products.... serial[%.3f] serial.mkl[%.3f] serial.mkl.alt[%.3f]\n", time_loop[0], time_loop[1],
			       time_loop[2]);
			printf("\t                            group [%.3f] group.mkl [%.3f] group.mkl.alt [%.3f]\n", time_loop[3], time_loop[4],
			       time_loop[5]);
			printf("\t                            ==> optimal.mix.strategy[%.3f]\n", time_loop[6]);
			printf("\t                                serial[%4.1f] serial.mkl[%4.1f] serial.mkl.alt[%4.1f]\n",
			       100 * time_loop[7], 100 * time_loop[8], 100 * time_loop[9]);
			printf("\t                                group [%4.1f] group.mkl [%4.1f] group.mkl.alt [%4.1f]\n",
			       100 * time_loop[10], 100 * time_loop[11], 100 * time_loop[12]);
#endif
#if !defined(INLA_WITH_MKL) && !defined(INLA_WITH_ARMPL)
			printf("\tOptimizing dot-products.... serial[%.3f] group[%.3f]\n", time_loop[0], time_loop[1]);
			printf("\t                            ==> optimal.mix.strategy[%.3f]\n", time_loop[2]);
			printf("\t                                serial[%4.1f] group[%4.1f]\n", 100 * time_loop[3], 100 * time_loop[4]);
#endif
		}
	}
	GMRFLib_openmp_implement_strategy(GMRFLib_OPENMP_PLACES_OPTIMIZE, NULL, NULL);

	c = Calloc_get(N);
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
			printf("\n\tList of hyperparameters: \n");
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

	mb->misc_output = Calloc(1, GMRFLib_ai_misc_output_tp);
	if (mb->output->config) {
		mb->misc_output->configs_preopt = Calloc(GMRFLib_MAX_THREADS(), GMRFLib_store_configs_preopt_tp *);
		mb->misc_output->config_lite = mb->output->config_lite;
	} else {
		mb->misc_output->configs_preopt = NULL;
		mb->misc_output->config_lite = 0;
	}

	mb->misc_output->likelihood_info = mb->output->likelihood_info;

	if (mb->lc_derived_correlation_matrix) {
		mb->misc_output->compute_corr_lin = mb->nlc;   /* yes, pass the dimension */
	} else {
		mb->misc_output->compute_corr_lin = 0;
	}

	x = Calloc_get(N);
	if (mb->mode_use_mode && mb->x_file) {
		Memcpy(x, mb->x_file + preopt->mnpred, N * sizeof(double));
	}

	mb->ai_par->compute_nparam_eff = 1;
	mb->predictor_compute = GMRFLib_TRUE;
	mb->ai_par->strategy = GMRFLib_AI_STRATEGY_GAUSSIAN;
	if (GMRFLib_gaussian_data) {
		mb->ai_par->vb_enable = GMRFLib_FALSE;
		mb->ai_par->step_len = 1.0;		       /* override the default in this particular case */
	}

	mb->transform_funcs = Calloc(N + preopt->mnpred, GMRFLib_transform_array_func_tp *);
	for (i = 0; i < preopt->mnpred; i++) {
		/*
		 * only where we have data (ie n or m), can be a invlinkfunc different from the identity. 
		 */
		if (i < mb->predictor_ndata && mb->predictor_invlinkfunc[i]) {
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

	if (mb->nlc > 0) {
		// postprocess the lincombs to convert APredictor and Predictor into sums of the latent
		int debug = 0;
		int mpred = preopt->mpred;		       // length(pApredictor)
		int mnpred = preopt->mnpred;		       // length(c(pApredictor, Apredictor))

		for (int k = 0; k < mb->nlc; k++) {
			GMRFLib_idxval_tp *idx = NULL;
			GMRFLib_lc_tp *lc = mb->lc_lc[k];
			for (int ii = 0; ii < lc->n; ii++) {
				i = lc->idx[ii];
				double w = lc->weight[ii];
				if (debug) {
					printf("lc[%1d] decode [idx= %1d, weight= %.8f]\n", k, i, w);
				}
				if (lc->idx[ii] < mnpred) {
					// replace this statement with a row of either pAA or A
					GMRFLib_idxval_tp *AA = NULL;
					if (lc->idx[ii] < mpred) {
						AA = preopt->pAA_idxval[lc->idx[ii]];
					} else {
						AA = preopt->A_idxval[lc->idx[ii] - mpred];
					}

					for (j = 0; j < AA->n; j++) {
						GMRFLib_idxval_addto(&idx, AA->idx[j] + mnpred, w * AA->val[j]);
					}
				} else {
					GMRFLib_idxval_addto(&idx, i, w);
				}
			}
			GMRFLib_idxval_sort(idx);
			if (debug) {
				GMRFLib_idxval_printf(stdout, idx, "");
			}

			Free(lc->idx);
			Free(lc->weight);
			lc->n = idx->n;
			lc->idx = idx->idx;
			lc->weight = idx->val;
		}
	}

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

	if (!mb->x_file && mb->compute_initial_values) {
		tref = -GMRFLib_timer();
		double *eta_pseudo = Calloc(preopt->Npred, double);

#pragma omp parallel for private(i) num_threads(GMRFLib_openmp->max_threads_outer)
		for (i = 0; i < preopt->Npred; i++) {
			if (mb->d[i]) {
				eta_pseudo[i] = inla_compute_initial_value(i, mb->loglikelihood[i], x, (void *) mb->loglikelihood_arg[i]);
			} else {
				eta_pseudo[i] = 0.0;
			}
		}

		double *eta = Calloc(preopt->Npred, double);
		double *Ad = Calloc(preopt->Npred, double);
		double *e = Calloc(preopt->Npred, double);
		double *bb = Calloc(preopt->n, double);
		double *scale = Calloc(preopt->n, double);
		double s0 = 1.0;

		if (preopt->pAAt_idxval) {
#pragma omp parallel for private(i) num_threads(GMRFLib_openmp->max_threads_outer)
			for (i = 0; i < preopt->n; i++) {
				double s = GMRFLib_dssqr(preopt->pAAt_idxval[i]->n, preopt->pAAt_idxval[i]->val);
				scale[i] = 1.0 / (s0 + DMAX(0.0, s));
			}
		} else {
#pragma omp parallel for private(i) num_threads(GMRFLib_openmp->max_threads_outer)
			for (i = 0; i < preopt->n; i++) {
				double s = GMRFLib_dsum(preopt->AtA_idxval[i][0]->n, preopt->AtA_idxval[i][0]->val);
				scale[i] = 1.0 / (s0 + DMAX(0.0, s));
			}
		}

		int iter_max = 10;
		double norm_initial = 0.0;
		double *d = bb;				       /* use the same space for both */
		for (int iter = 0; iter < iter_max; iter++) {
			double norm = 0.0, sum1 = 0.0, sum2 = 0.0, gamma;

			GMRFLib_preopt_predictor(eta, x, preopt);
			GMRFLib_daxpbyz(preopt->Npred, 1.0, eta_pseudo, -1.0, eta, e);
			norm = sqrt(GMRFLib_ddot(preopt->Npred, e, e) / preopt->Npred);
			GMRFLib_preopt_bnew_like(bb, e, preopt);
			GMRFLib_mul(preopt->n, d, scale, d);
			GMRFLib_preopt_predictor(Ad, d, preopt);
			sum1 = GMRFLib_ddot(preopt->Npred, Ad, e);
			sum2 = GMRFLib_ddot(preopt->Npred, Ad, Ad);
			gamma = DMAX(0.0, DMIN(2.0, sum1 / (FLT_EPSILON + sum2)));

			if (iter == 0) {
				norm_initial = norm;
			}
			if (mb->verbose) {
				if (iter == 0)
					printf("\nCompute initial values...\n");
				printf("\tIter[%1d] RMS(err) = %.3f, update with step-size = %.3f\n", iter, norm / norm_initial, gamma);
			}

			GMRFLib_daxpy(preopt->n, gamma, d, x);
			if (iter > 0) {
				if ((norm_initial - norm) / norm_initial < 0.05)
					break;
			}
			norm_initial = norm;
		}

		tref += GMRFLib_timer();
		if (mb->verbose) {
			printf("\tInitial values computed in %.4f seconds\n", tref);
			for (i = 0; i < IMIN(preopt->n, PREVIEW / 2L); i++) {
				printf("\t\tx[%1d] = %.4f\n", i, x[i]);
			}
			for (i = IMAX(0, preopt->n - PREVIEW / 2L); i < preopt->n; i++) {
				printf("\t\tx[%1d] = %.4f\n", i, x[i]);
			}
			printf("\n");
		}

		Free(eta_pseudo);
		Free(eta);
		Free(Ad);
		Free(e);
		Free(bb);
		Free(scale);
	}

	GMRFLib_ai_INLA_experimental(&(mb->density),
				     NULL, NULL,
				     (mb->output->hyperparameters ? &(mb->density_hyper) : NULL),
				     (mb->output->gcpo ? &(mb->gcpo) : NULL), mb->gcpo_param,
				     (mb->output->cpo || mb->expert_cpo_manual ? &(mb->cpo) : NULL),
				     (mb->output->po ? &(mb->po) : NULL),
				     mb->dic,
				     (mb->output->mlik ? &(mb->mlik) : NULL),
				     mb->theta, mb->ntheta,
				     extra, (void *) mb,
				     x, b, c, NULL, bfunc, prior_mean, mb->d, mb->fl,
				     loglikelihood_inla, (void *) mb,
				     preopt->preopt_graph, preopt->preopt_Qfunc, preopt->preopt_Qfunc_arg, preopt->latent_constr,
				     mb->ai_par, ai_store, mb->nlc, mb->lc_lc, &(mb->density_lin), mb->misc_output, preopt);


	/*
	 * add the offsets to the linear predictor. Add the offsets to the 'configs' (if any), at a later stage. 
	 */
#pragma omp parallel for private(i) num_threads(GMRFLib_openmp->max_threads_outer)
	for (i = 0; i < mb->predictor_n + mb->predictor_m; i++) {
		GMRFLib_density_tp *d = NULL;
		if (mb->density[i] && !ISZERO(OFFSET3(i))) {
			d = mb->density[i];
			GMRFLib_density_new_mean(&(mb->density[i]), d, d->std_mean + OFFSET3(i));
			GMRFLib_free_density(d);
		}
	}

	GMRFLib_free_ai_store(ai_store);
	Free(bfunc);
	Calloc_free();

	return INLA_OK;
}

int inla_computed(GMRFLib_density_tp **d, int n)
{
	/*
	 * return 0 if all d[i]'s are NULL, and 1 otherwise. 
	 */
	int i;

	if (!d) {
		return INLA_OK;
	}
	for (i = 0; i < n; i++) {
		if (d[i]) {
			return INLA_FAIL;
		}
	}
	return INLA_OK;
}

int inla_integrate_func(double *d_mean, double *d_stdev, double *d_mode, GMRFLib_density_tp *density, map_func_tp *func,
			void *func_arg, GMRFLib_transform_array_func_tp *tfunc)
{
	// this require 'i_max', 'np', 'z' and 'ldz'
#define COMPUTE_MODE()							\
	if (d_mode) {							\
		int ii = i_max;						\
		if (ii > 0 && ii < np - 1) {				\
			ii--;						\
		} else if (ii == np - 1) {				\
			ii -= 2;					\
		}							\
		double zm = inla_interpolate_mode(z + ii, ldz + ii);	\
		double zm_orig = zm;					\
		int m = 5;						\
		int low_ = IMAX(0, i_max - m);				\
		int high_ = IMIN(np-1, i_max + m);			\
		int len = high_ - low_ + 1;				\
		GMRFLib_spline_tp *lds = GMRFLib_spline_create_x(z + low_, ldz + low_, len, GMRFLib_INTPOL_TRANS_NONE, GMRFLib_INTPOL_CACHE_LEVEL1); \
		if (lds == NULL) {					\
			*d_mode = zm_orig;				\
		} else {						\
			double step_size[] = {0.0, 0.0};		\
			for(int iter = 0; iter < 2; iter++) {		\
				step_size[iter] = GMRFLib_spline_eval_deriv(zm, lds) / GMRFLib_spline_eval_deriv2(zm, lds); \
				zm -= step_size[iter];			\
			}						\
			if ((ABS(step_size[0]) >= ABS(step_size[1]))) {	\
				*d_mode = zm;				\
			} else {					\
				/* emergency option */			\
				*d_mode = zm_orig;			\
			}						\
			GMRFLib_spline_free(lds);			\
			}						\
	}

#define _MAP_X(_x_user) (plain_case ? (_x_user) : \
			 (func ? func(_x_user, MAP_FORWARD, func_arg) : \
			  (tfunc ? tfunc->func(thread_id, _x_user, GMRFLib_TRANSFORM_FORWARD, tfunc->arg, tfunc->cov) : \
			   (_x_user))))
#define _MAP_X_plain(_x_user) (_x_user)
#define _MAP_X_func(_x_user) func(_x_user, MAP_FORWARD, func_arg)
#define _MAP_X_tfunc(_x_user) tfunc->func(thread_id, _x_user, GMRFLib_TRANSFORM_FORWARD, tfunc->arg, tfunc->cov)

#define _MAP_DX(_x_user) (plain_case ? DSIGN(_x_user) : \
			  (func ? func(_x_user, MAP_DFORWARD, func_arg) : \
			   (tfunc ? tfunc->func(thread_id, _x_user, GMRFLib_TRANSFORM_DFORWARD, tfunc->arg, tfunc->cov) : \
			    DSIGN(_x_user))))
#define _MAP_DX_plain(_x_user) DSIGN(_x_user)
#define _MAP_DX_func(_x_user) func(_x_user, MAP_DFORWARD, func_arg)
#define _MAP_DX_tfunc(_x_user) tfunc->func(thread_id, _x_user, GMRFLib_TRANSFORM_DFORWARD, tfunc->arg, tfunc->cov)

	int thread_id = 0;
	const int plain = (!func && !tfunc);

	/*
	 * We need to integrate to get the transformed mean and variance. Use a simple Simpsons-rule.  The simple mapping we did before was not good enough,
	 * obviously... 
	 */
	if (!func && !tfunc) {
		*d_mean = density->user_mean;
		*d_stdev = density->user_stdev;
		if (d_mode) {
			*d_mode = density->user_mode;
		}
		return GMRFLib_SUCCESS;
	}

	if (!density) {
		*d_mean = NAN;
		*d_stdev = NAN;
		if (d_mode) {
			*d_mode = NAN;
		}
		return GMRFLib_SUCCESS;
	}

	// GMRFLib_ENTER_ROUTINE;

	int np = GMRFLib_INT_NUM_POINTS;
	int npm = GMRFLib_INT_NUM_INTERPOL * np - (GMRFLib_INT_NUM_INTERPOL - 1);
	double low = 0.0, high = 0.0, *xpm = NULL, *ld = NULL, *ldm = NULL, *xp = NULL, *xx = NULL, dx = 0.0, m0 = 0.0, m1 = 0.0, m2 = 0.0;

	static double *w = NULL;
	if (!w) {
#pragma omp critical (Name_1eca1953ff3f841bda5736d32498e385b608fa2c)
		if (!w) {
			double wref[] = { 4.0, 2.0 };
			double *ww = Calloc(npm, double);
			ww[0] = ww[npm - 1] = 1.0;
			for (int i = 1, j = 0; i < npm - 1; i++, j = (j + 1) % 2L) {
				ww[i] = wref[j];
			}
			w = ww;
		}
	}

	if (density->type == GMRFLib_DENSITY_TYPE_GAUSSIAN) {
		// then we can do better
		np = GMRFLib_INT_GHQ_POINTS;
		double *wp = NULL;
		double mean = density->user_mean;
		double stdev = density->user_stdev;

		GMRFLib_ghq(&xp, &wp, np);

		Calloc_init(2 * np, 2);
		double *ldz = Calloc_get(np);
		double *z = Calloc_get(np);

		int i_max = 0;
		m1 = 0.0;
		m2 = 0.0;

		if (d_mode) {
			if (plain) {
				for (int i = 0; i < np; i++) {
					double x = xp[i] * stdev + mean;
					double f = _MAP_X_plain(x);
					double df = _MAP_DX_plain(x);
					m1 += wp[i] * f;
					m2 += wp[i] * SQR(f);

					z[i] = f;
					ldz[i] = -0.5 * SQR(xp[i]) - log(ABS(df));
					if ((i == 0) || ldz[i] > ldz[i_max]) {
						i_max = i;
					}
				}
			} else if (func) {
				for (int i = 0; i < np; i++) {
					double x = xp[i] * stdev + mean;
					double f = _MAP_X_func(x);
					double df = _MAP_DX_func(x);
					m1 += wp[i] * f;
					m2 += wp[i] * SQR(f);

					z[i] = f;
					ldz[i] = -0.5 * SQR(xp[i]) - log(ABS(df));
					if ((i == 0) || ldz[i] > ldz[i_max]) {
						i_max = i;
					}
				}
			} else if (tfunc) {
				for (int i = 0; i < np; i++) {
					double x = xp[i] * stdev + mean;
					double f = _MAP_X_tfunc(x);
					double df = _MAP_DX_tfunc(x);
					m1 += wp[i] * f;
					m2 += wp[i] * SQR(f);

					z[i] = f;
					ldz[i] = -0.5 * SQR(xp[i]) - log(ABS(df));
					if ((i == 0) || ldz[i] > ldz[i_max]) {
						i_max = i;
					}
				}
			} else {
				assert(0 == 1);
			}
		} else {
			if (plain) {
				for (int i = 0; i < np; i++) {
					double x = xp[i] * stdev + mean;
					double f = _MAP_X_plain(x);
					m1 += wp[i] * f;
					m2 += wp[i] * SQR(f);
				}
			} else if (func) {
				for (int i = 0; i < np; i++) {
					double x = xp[i] * stdev + mean;
					double f = _MAP_X_func(x);
					m1 += wp[i] * f;
					m2 += wp[i] * SQR(f);
				}
			} else if (tfunc) {
				for (int i = 0; i < np; i++) {
					double x = xp[i] * stdev + mean;
					double f = _MAP_X_tfunc(x);
					m1 += wp[i] * f;
					m2 += wp[i] * SQR(f);
				}
			} else {
				assert(0 == 1);
			}
		}

		*d_mean = m1;
		*d_stdev = sqrt(DMAX(0.0, m2 - SQR(m1)));

		COMPUTE_MODE();
		Calloc_free();
	} else {
		Calloc_init(3 * npm + 4 * np, 7);
		low = density->x_min;
		high = density->x_max;
		dx = (high - low) / (np - 1.0);
		xp = Calloc_get(np);
		ld = Calloc_get(np);

#pragma omp simd
		for (int i = 0; i < np; i++) {
			xp[i] = low + i * dx;
		}
		GMRFLib_evaluate_nlogdensity(ld, xp, np, density);

		int i_max = 0;
		double *z = Calloc_get(np);
		double *ldz = Calloc_get(np);

		if (d_mode) {
			// reusing 'z' for x_user here
			GMRFLib_density_std2user_n(z, xp, np, density);

			if (plain) {
				Memcpy(ldz, ld, np * sizeof(double));
			} else {
				if (plain) {
					for (int i = 0; i < np; i++) {
						ldz[i] = ld[i] - log(ABS(_MAP_DX_plain(z[i])));
						z[i] = _MAP_X_plain(z[i]);
					}
				} else if (func) {
					for (int i = 0; i < np; i++) {
						ldz[i] = ld[i] - log(ABS(_MAP_DX_func(z[i])));
						z[i] = _MAP_X_func(z[i]);
					}
				} else if (tfunc) {
					for (int i = 0; i < np; i++) {
						ldz[i] = ld[i] - log(ABS(_MAP_DX_tfunc(z[i])));
						z[i] = _MAP_X_tfunc(z[i]);
					}
				} else {
					assert(0 == 1);
				}
			}
			GMRFLib_max_value(ldz, np, &i_max);
		}

		// interpolate
		xpm = Calloc_get(npm);
		ldm = Calloc_get(npm);

		if (GMRFLib_INT_NUM_INTERPOL == 3) {
			const double div3 = 1.0 / 3.0;
#pragma omp simd
			for (int i = 0; i < np - 1; i++) {
				xpm[3 * i] = xp[i];
				xpm[3 * i + 1] = (2.0 * xp[i] + xp[i + 1]) * div3;
				xpm[3 * i + 2] = (xp[i] + 2.0 * xp[i + 1]) * div3;
			}
#pragma omp simd
			for (int i = 0; i < np - 1; i++) {
				ldm[3 * i] = ld[i];
				ldm[3 * i + 1] = (2.0 * ld[i] + ld[i + 1]) * div3;
				ldm[3 * i + 2] = (ld[i] + 2.0 * ld[i + 1]) * div3;
			}
			xpm[3 * (np - 2) + 3] = xp[np - 1];
			ldm[3 * (np - 2) + 3] = ld[np - 1];
			assert(3 * (np - 2) + 3 == npm - 1);
		} else if (GMRFLib_INT_NUM_INTERPOL == 2) {
			const double div2 = 0.5;
#pragma omp simd
			for (int i = 0; i < np - 1; i++) {
				xpm[2 * i] = xp[i];
				xpm[2 * i + 1] = (xp[i] + xp[i + 1]) * div2;
			}
#pragma omp simd
			for (int i = 0; i < np - 1; i++) {
				ldm[2 * i] = ld[i];
				ldm[2 * i + 1] = (ld[i] + ld[i + 1]) * div2;
			}
			xpm[2 * (np - 2) + 2] = xp[np - 1];
			ldm[2 * (np - 2) + 2] = ld[np - 1];
			assert(2 * (np - 2) + 2 == npm - 1);
		} else {
			assert(GMRFLib_INT_NUM_INTERPOL == 2 || GMRFLib_INT_NUM_INTERPOL == 3);
		}

		GMRFLib_exp(npm, ldm, ldm);
		xx = Calloc_get(npm);
		GMRFLib_density_std2user_n(xx, xpm, npm, density);
		if (plain) {
			// for (int i = 0; i < npm; i++) {
			// xx[i] = xx[i];
			// }
		} else if (func) {
			for (int i = 0; i < npm; i++) {
				xx[i] = _MAP_X_func(xx[i]);
			}
		} else if (tfunc) {
			for (int i = 0; i < npm; i++) {
				xx[i] = _MAP_X_tfunc(xx[i]);
			}
		} else {
			assert(0 == 1);
		}

		GMRFLib_mul(npm, ldm, w, ldm);
		m0 = GMRFLib_dsum(npm, ldm);
		m1 = GMRFLib_ddot(npm, ldm, xx);
		GMRFLib_sqr(npm, xx, xx);
		m2 = GMRFLib_ddot(npm, ldm, xx);

		m1 /= m0;
		m2 /= m0;

		*d_mean = m1;
		*d_stdev = sqrt(DMAX(0.0, m2 - SQR(m1)));

		COMPUTE_MODE();
		Calloc_free();
	}

#undef COMPUTE_MODE
#undef _MAP_X
#undef _MAP_X_plain
#undef _MAP_X_func
#undef _MAP_X_tfunc
#undef _MAP_DX
#undef _MAP_DX_plain
#undef _MAP_DX_func
#undef _MAP_DX_tfunc

//      GMRFLib_LEAVE_ROUTINE;
	return GMRFLib_SUCCESS;
}

int inla_divisible(int n, int by)
{
	/*
	 * same function as inla.divisi 
	 */

	if (by == 0)
		return GMRFLib_TRUE;

	if (by > 0)
		return (by * (n / by) == n ? GMRFLib_TRUE : GMRFLib_FALSE);
	else
		return ((-by) * (n / (-by)) == n ? GMRFLib_FALSE : GMRFLib_TRUE);
}

int inla_besag_scale(int thread_id, inla_besag_Qfunc_arg_tp *arg, int adj, int verbose)
{
	// if VERBOSE, write out the scalings.
	inla_besag_Qfunc_arg_tp *def = Calloc(1, inla_besag_Qfunc_arg_tp);
	int i, k, *cc = NULL, n = arg->graph->n;
	const int debug = 0;
	arg->prec_scale = Calloc(arg->graph->n, double);

	if (debug)
		P(adj);

	if (adj) {
		// use the cc in the graph
		cc = GMRFLib_graph_cc(arg->graph);
	} else {
		// treat the whole graph as one cc, with index 0
		cc = Calloc(n, int);
	}

	if (!adj) {
		// special for nnbs=0
		for (i = 0; i < n; i++) {
			if (arg->graph->nnbs[i] == 0) {
				cc[i] = 1;		       // so it will be disabled in the computations below
			}
		}
	}
	// if !adj, then we will not use the nodes where nnbs=0, and we do this by forcing ncc=0,
	// since cc[i]=1 for those nodes as set above.
	int ncc;
	ncc = (adj ? 1 + GMRFLib_imax_value(cc, arg->graph->n, NULL) : 1);
	if (debug)
		P(ncc);

	// work with each cc at the time
	for (k = 0; k < ncc; k++) {
		if (debug) {
			printf("Working with cc %d\n", k);
		}

		GMRFLib_constr_tp *constr = NULL;
		GMRFLib_make_empty_constr(&constr);

		char *remove = Calloc(n, char);
		int num = 0;
		for (i = num = 0; i < n; i++) {
			if (cc[i] == k) {
				remove[i] = 0;
				num++;
			} else {
				remove[i] = 1;
			}
		}
		if (debug) {
			printf("\tsize is %d\n", num);
		}

		if (num == 1) {
			// only one node in this component. then we know that the precision i 1.
			for (i = 0; i < n; i++) {
				if (cc[i] == k) {
					arg->prec_scale[i] = -1.0;	/* this is code for treating this case specially */
				}
			}
			if (verbose)
				printf("\t\tconnected component[%1d] size[%1d] scale[%.6g]\n", k, num, -1.0);
		} else {
			// compute the subgraph and find the scaling for this connected component
			GMRFLib_graph_comp_subgraph(&(def->graph), arg->graph, remove, NULL);

			constr->nc = 1;
			constr->a_matrix = Calloc(def->graph->n, double);
			for (i = 0; i < def->graph->n; i++) {
				constr->a_matrix[i] = 1.0;
			}
			constr->e_vector = Calloc(1, double);
			GMRFLib_prepare_constr(constr, def->graph, GMRFLib_TRUE);

			GMRFLib_problem_tp *problem = NULL;
			int retval = GMRFLib_SUCCESS, ok = 0, num_try = 0, num_try_max = 100;
			GMRFLib_error_handler_tp *old_handler = GMRFLib_set_error_handler_off();

			double *c = Calloc(def->graph->n, double), eps = GSL_SQRT_DBL_EPSILON;
			for (i = 0; i < def->graph->n; i++) {
				c[i] = eps;
			}

			while (!ok) {
				retval =
				    GMRFLib_init_problem(thread_id, &problem, NULL, NULL, c, NULL, def->graph, Qfunc_besag, (void *) def, constr);
				switch (retval) {
				case GMRFLib_EPOSDEF:
				{
					for (i = 0; i < def->graph->n; i++) {
						c[i] *= 10.0;
					}
					problem = NULL;
				}
					break;

				case GMRFLib_SUCCESS:
				{
					ok = 1;
				}
					break;
				default:
				{
					GMRFLib_set_error_handler(old_handler);
					GMRFLib_ERROR(retval);
					abort();
				}
					break;
				}

				if (++num_try >= num_try_max) {
					FIXME("This should not happen. Contact developers...");
					abort();
				}
			}
			GMRFLib_set_error_handler(old_handler);
			GMRFLib_Qinv(problem);

			if (debug)
				P(def->graph->n);
			double sum = 0.0, value;

			for (i = 0; i < def->graph->n; i++) {
				sum += log(*(GMRFLib_Qinv_get(problem, i, i)));
			}
			value = exp(sum / def->graph->n);
			if (debug) {
				printf("\tprec_scale is %f\n", value);
			}
			for (i = 0; i < n; i++) {
				if (cc[i] == k) {
					arg->prec_scale[i] = value;
				}
			}
			if (verbose)
				printf("\t\tconnected component[%1d] size[%1d] scale[%.6g]\n", k, def->graph->n, value);

			GMRFLib_graph_free(def->graph);
			GMRFLib_free_problem(problem);
			GMRFLib_free_constr(constr);
			Free(c);
		}
	}

	if (!adj) {
		for (i = 0; i < n; i++) {
			if (cc[i] > 0) {
				assert(arg->graph->nnbs[i] == 0);
				arg->prec_scale[i] = -1.0;     /* this is code for treating this case specially */
			}
		}
	}

	Free(def);
	Free(cc);

	return GMRFLib_SUCCESS;
}

double inla_update_density(double *theta, inla_update_tp *arg)
{
	/*
	 * joint posterior for theta
	 */

	int i, corr = (arg->stdev_corr_pos && arg->stdev_corr_neg);
	double value = 0.0, sd, log_nc, update_dens, *z = NULL;

	z = Calloc(arg->ntheta, double);
	GMRFLib_ai_theta2z(z, arg->ntheta, arg->theta_mode, theta, arg->sqrt_eigen_values, arg->eigen_vectors);

	for (i = 0; i < arg->ntheta; i++) {
		if (corr) {
			sd = (z[i] > 0 ? arg->stdev_corr_pos[i] : arg->stdev_corr_neg[i]);
		} else {
			sd = 1.0;
		}
		value += -0.5 * SQR(z[i] / sd);
	}

	/*
	 * this is the normalizing constant
	 */
	log_nc = 0.5 * arg->ntheta * log(2.0 * M_PI);
	for (i = 0; i < arg->ntheta; i++) {
		if (corr) {
			log_nc += -0.5 * (2.0 * log(gsl_vector_get(arg->sqrt_eigen_values, (unsigned int) i)) +
					  0.5 * (log(SQR(arg->stdev_corr_pos[i])) + log(SQR(arg->stdev_corr_neg[i]))));
		} else {
			log_nc += -0.5 * (2.0 * log(gsl_vector_get(arg->sqrt_eigen_values, (unsigned int) i)));
		}
	}

	/*
	 * and then the log-joint-posterior-density
	 */
	update_dens = value - log_nc;

	Free(z);
	return update_dens;
}

double inla_dmatern_cf(double dist, double range, double nu)
{
	double kappa = sqrt(8.0 * nu) / range;
	double dd = kappa * dist;
	double cf;

	cf = (dist <= 1e-12 ? 1.0 : 1.0 / pow(2.0, nu - 1.0) / MATHLIB_FUN(gammafn) (nu) * pow(dd, nu) * MATHLIB_FUN(bessel_k) (dd, nu, 1.0));

	return (cf);
}

double inla_sn_intercept(double intern_quantile, double skew)
{
	// testing only
	double a3[2] = { 0.0, 0.0 }, val;
	a3[0] = POW3(inla_pc_sn_skew2alpha(skew));
	val = map_invsn(intern_quantile, MAP_BACKWARD, (void *) a3);
	P(intern_quantile);
	P(skew);
	P(inla_pc_sn_skew2alpha(skew));
	P(a3[0]);
	P(val);

	return (0);

	a3[0] = POW3(inla_pc_sn_skew2alpha(skew));
	return (map_invsn(intern_quantile, MAP_FORWARD, (void *) a3));
}

int inla_reset(void)
{
	// reset static variables various places as need need to call _ai_INLA() twice in preopt_mode

	GMRFLib_opt_exit();
	GMRFLib_pardiso_exit();
	R_rgeneric_cputime = 0.0;

	return GMRFLib_SUCCESS;
}

int main(int argc, char **argv)
{
#define _USAGE_intern(fp)  fprintf(fp, "\nUsage: %s [-v] [-V] [-h] [-f] [-eVAR=VALUE] [-tA:B] [-mMODE] FILE.INI\n", program)
#define _USAGE _USAGE_intern(stderr)
#define _HELP _USAGE_intern(stdout);					\
	printf("\t\t-v\t: Verbose output.\n");				\
	printf("\t\t-V\t: Print version and exit.\n");			\
	printf("\t\t-b\t: Use binary output-files.\n");			\
	printf("\t\t-s\t: Be silent.\n");				\
	printf("\t\t-c\t: Create core-file if needed (and allowed). (Linux/MacOSX only.)\n"); \
	printf("\t\t-R\t: Restart using previous mode.\n");		\
	printf("\t\t-e VAR=VALUE\t: Set variable VAR to VALUE.\n");	\
	printf("\t\t-t A:B\t: the number of threads (A=outer,B=inner), 0 means auto\n"); \
	printf("\t\t-m MODE\t: Enable special mode:\n");		\
	printf("\t\t\tHYPER :  Enable HYPERPARAMETER mode\n");		\
	printf("\t\t-h\t: Print (this) help.\n")

#define _BUGS_intern(fp) fprintf(fp, "Report bugs to <help@r-inla.org>\n")
#define _BUGS _BUGS_intern(stdout)
	int i, verbose = 0, silent = 0, opt, arg, ntt[2] = { 0, 0 }, err;
#if !defined(WINDOWS)
	int enable_core_file = 0;			       /* allow for core files */
#endif
	char *program = argv[0];
	double time_used[4] = { -1, -1, -1, -1 };
#if !defined(WINDOWS)
	double eff_nt = 0.0;
#endif
	clock_t atime_used[4] = { 0, 0, 0, 0 };
	inla_tp *mb = NULL;

	int host_max_threads = IMAX(omp_get_max_threads(), omp_get_num_procs());

	GMRFLib_malloc_debug_check();

	GMRFLib_openmp = Calloc(1, GMRFLib_openmp_tp);
	GMRFLib_openmp->max_threads = host_max_threads;
	GMRFLib_openmp->blas_num_threads_force = 0;
	GMRFLib_openmp->max_threads_nested = Calloc(2, int);
	GMRFLib_openmp->max_threads_nested[0] = GMRFLib_openmp->max_threads;
	GMRFLib_openmp->max_threads_nested[1] = 1;
	GMRFLib_openmp->adaptive = GMRFLib_FALSE;
	GMRFLib_openmp->strategy = GMRFLib_OPENMP_STRATEGY_DEFAULT;
	GMRFLib_openmp_implement_strategy(GMRFLib_OPENMP_PLACES_DEFAULT, NULL, NULL);

	GMRFLib_bitmap_max_dimension = 512;
	GMRFLib_bitmap_swap = GMRFLib_TRUE;
	GMRFLib_aqat_m_diag_add = GSL_SQRT_DBL_EPSILON;
	GMRFLib_gaussian_data = 1;
	GMRFLib_taucs_sort_L = 0;

	GMRFLib_init_constr_store();
	GMRFLib_init_constr_store_logdet();		       /* no need to reset this with preopt */
	GMRFLib_graph_init_store();			       /* no need to reset this with pretop */
	GMRFLib_csr_init_store();
	GMRFLib_trace_functions(NULL);
	GMRFLib_debug_functions(NULL);
	GMRFLib_reorder = G.reorder;
	GMRFLib_inla_mode = GMRFLib_MODE_COMPACT;
	my_sort2_id_test_cutoff(0);
	my_sort2_dd_test_cutoff(0);

	GMRFLib_cachelinesize = GMRFLib_get_cachelinesize();

	/*
	 * special option: if one of the arguments is `--ping', then just return INLA[<VERSION>] IS ALIVE 
	 */
	for (i = 1; i < argc; i++) {
		if (!strcasecmp(argv[i], "-ping") || !strcasecmp(argv[i], "--ping")) {
			printf("INLA[%s] IS ALIVE\n", GITCOMMIT);
			exit(EXIT_SUCCESS);
		}
	}

#if !defined(WINDOWS)
	signal(SIGUSR1, inla_signal);
	signal(SIGUSR2, inla_signal);
	signal(SIGINT, inla_signal);
#endif
	while ((opt = getopt(argc, argv, "vVe:t:B:m:S:z:hsr:R:cpLP:")) != -1) {
		switch (opt) {
		case 'P':
		{
			if (!strcasecmp(optarg, "CLASSIC") || !strcasecmp(optarg, "CLASSICAL")) {
				GMRFLib_inla_mode = GMRFLib_MODE_CLASSIC;
			} else if (!strcasecmp(optarg, "EXPERIMENTAL") || !strcasecmp(optarg, "COMPACT")) {
				GMRFLib_inla_mode = GMRFLib_MODE_COMPACT;
			} else {
				assert(0 == 1);
			}
		}
			break;

		case 'v':
		{
			silent = 1;
			verbose++;
		}
			break;

		case 'V':
		{
			printf("This program has version:\n\t%s\nand is linked with ", GITCOMMIT);
			GMRFLib_version(stdout);
			_BUGS;
			exit(EXIT_SUCCESS);
		}
			break;
		case 'e':
		{
			my_setenv(optarg, 1);
		}
			break;

		case 'B':
		{
			int bnt;
			if (inla_sread_ints(&bnt, 1, optarg) == INLA_OK) {
				bnt = IMAX(bnt, 0);
				GMRFLib_openmp->blas_num_threads_force = bnt;
			} else {
				fprintf(stderr, "Fail to read BLAS_NUM_THREADS from %s\n", optarg);
				exit(EXIT_SUCCESS);
			}
		}
			break;

		case 'm':
		{
			if (!strncasecmp(optarg, "HYPER", 5)) {
				G.mode = INLA_MODE_HYPER;
			} else if (!strncasecmp(optarg, "QINV", 4)) {
				G.mode = INLA_MODE_QINV;
			} else if (!strncasecmp(optarg, "QSOLVE", 6)) {
				G.mode = INLA_MODE_QSOLVE;
			} else if (!strncasecmp(optarg, "QREORDERING", 11)) {
				G.mode = INLA_MODE_QREORDERING;
			} else if (!strncasecmp(optarg, "QSAMPLE", 7)) {
				G.mode = INLA_MODE_QSAMPLE;
			} else if (!strncasecmp(optarg, "FINN", 4)) {
				G.mode = INLA_MODE_FINN;
			} else if (!strncasecmp(optarg, "GRAPH", 5)) {
				G.mode = INLA_MODE_GRAPH;
			} else if (!strncasecmp(optarg, "R", 1)) {
				G.mode = INLA_MODE_R;
			} else if (!strncasecmp(optarg, "FGN", 3)) {
				G.mode = INLA_MODE_FGN;
			} else if (!strncasecmp(optarg, "PARDISO", 7)) {
				G.mode = INLA_MODE_PARDISO;
			} else if (!strncasecmp(optarg, "OPENMP", 6)) {
				G.mode = INLA_MODE_OPENMP;
			} else if (!strncasecmp(optarg, "DRYRUN", 6)) {
				G.mode = INLA_MODE_DRYRUN;
			} else if (!strncasecmp(optarg, "TESTIT", 6)) {
				G.mode = INLA_MODE_TESTIT;
			} else {
				fprintf(stderr, "\n*** Error: Unknown mode (argument to '-m') : %s\n", optarg);
				exit(EXIT_FAILURE);
			}
		}
			break;

		case 'S':
		{
			// this option is only used for other MODES than INLA, like qsample
			inla_tolower(optarg);
			if (!strcasecmp(optarg, "default")) {
				if (GMRFLib_pardiso_check_install(1, 1) == GMRFLib_SUCCESS ? 1 : 0) {
					GMRFLib_smtp = GMRFLib_SMTP_PARDISO;
				} else {
					GMRFLib_smtp = GMRFLib_SMTP_TAUCS;
				}
			} else if (!strcasecmp(optarg, "taucs")) {
				GMRFLib_smtp = GMRFLib_SMTP_TAUCS;
			} else if (!strcasecmp(optarg, "band")) {
				GMRFLib_smtp = GMRFLib_SMTP_BAND;
			} else if (!strcasecmp(optarg, "pardiso")) {
				GMRFLib_smtp = GMRFLib_SMTP_PARDISO;
				GMRFLib_openmp->strategy = GMRFLib_OPENMP_STRATEGY_PARDISO;
			}
			GMRFLib_openmp_implement_strategy(GMRFLib_OPENMP_PLACES_EXTERNAL, NULL, &GMRFLib_smtp);
		}
			break;

		case 't':
		{
			if (inla_sread_colon_ints(&ntt[0], &ntt[1], optarg) == INLA_OK || inla_sread(ntt, 1, optarg, 0) == INLA_OK) {

				if (verbose > 0) {
					printf("\tRead ntt %d %d with max.threads %d\n", ntt[0], ntt[1], GMRFLib_openmp->max_threads);
				}

				// a hidden option...  enable also if ntt[1] > 1, not only if < 0.
				if (ntt[1] > 1 || ntt[1] < 0) {
					ntt[1] = IABS(ntt[1]);
					GMRFLib_openmp->adaptive = GMRFLib_TRUE;
				}

				for (i = 0; i < 2; i++) {
					ntt[i] = IMAX(0, ntt[i]);
				}

				// replace 0 with auto-values
				if (ntt[0] == 0 && ntt[1] == 0) {
					ntt[0] = GMRFLib_openmp->max_threads;
					ntt[1] = 1;
				} else if (ntt[0] == 0 && ntt[1] > 0) {
					ntt[0] = GMRFLib_openmp->max_threads / ntt[1];
				} else if (ntt[1] == 0) {
					// let 0 means 1 for the moment. only larger problems gives a speedup,
					// much likely it will slow things down.
					// ntt[1] = GMRFLib_openmp->max_threads / ntt[0];
					ntt[1] = 1;
				}

				if (ntt[0] * ntt[1] > GMRFLib_MAX_THREADS()) {
					fprintf(stderr, "\n\n\tYou ask for %1d x %1d = %1d number of threads,\n", ntt[0], ntt[1], ntt[0] * ntt[1]);
					fprintf(stderr, "\twhich is more that I got from the system: %1d\n", GMRFLib_MAX_THREADS());

					if (ntt[0] > GMRFLib_MAX_THREADS()) {
						ntt[0] = GMRFLib_MAX_THREADS();
						ntt[1] = 1;
					} else if (ntt[1] > GMRFLib_MAX_THREADS()) {
						ntt[1] = GMRFLib_MAX_THREADS();
						ntt[0] = 1;
					} else {
						// something gotta give
						while (ntt[0] * ntt[1] > GMRFLib_MAX_THREADS()) {
							ntt[1]--;
						}
					}
					fprintf(stderr, "\tNumber of threads is reduced to %1d:%1d\n\n", ntt[0], ntt[1]);
				}

				for (i = 0; i < 2; i++) {
					ntt[i] = IMAX(1, ntt[i]);
					GMRFLib_openmp->max_threads_nested[i] = ntt[i];
				}
				GMRFLib_openmp->max_threads = IMIN(GMRFLib_MAX_THREADS(), ntt[0] * ntt[1]);
			} else {
				fprintf(stderr, "Fail to read A:B from [%s]\n", optarg);
				fprintf(stderr, "Will continue with '4:1'\n");
				ntt[0] = 4;
				ntt[1] = 1;
				for (i = 0; i < 2; i++) {
					ntt[i] = IMIN(GMRFLib_openmp->max_threads, IMAX(1, ntt[i]));
					GMRFLib_openmp->max_threads_nested[i] = ntt[i];
				}
				GMRFLib_openmp->max_threads = ntt[0] * ntt[1];
			}
			if (verbose > 0) {
				printf("\tFound num.threads = %1d:%1d max_threads = %1d\n", GMRFLib_openmp->max_threads_nested[0],
				       GMRFLib_openmp->max_threads_nested[1], GMRFLib_openmp->max_threads);
			}
			omp_set_num_threads(GMRFLib_MAX_THREADS());
			GMRFLib_openmp_implement_strategy(GMRFLib_OPENMP_PLACES_DEFAULT, NULL, NULL);
		}
			break;

		case 'z':
		{
			if (!(G.mode == INLA_MODE_FINN || G.mode == INLA_MODE_QSAMPLE || G.mode == INLA_MODE_TESTIT)) {
				fprintf(stderr, "\n *** ERROR *** Option `-z seed' only available in selected modes\n");
				exit(EXIT_FAILURE);
			} else {
				int int_seed = 0;
				if (inla_sread_ints(&int_seed, 1, optarg) == INLA_OK) {
					;
				} else {
					fprintf(stderr, "Fail to read FINN_SEED from %s\n", optarg);
					exit(EXIT_SUCCESS);
				}
				if (int_seed != 0) {
					/*
					 * this is only for the main thread...
					 */
					GMRFLib_rng_init((unsigned long int) int_seed);
				}
			}
		}
			break;

		case 'h':
		{
			_HELP;
			_BUGS;
			exit(EXIT_SUCCESS);
		}
			break;

		case 's':
		{
			verbose = 0;
			silent = 1;
		}
			break;

		case 'r':
		{
			int itmp = 0;
			err = inla_sread_ints(&itmp, 1, optarg);
			G.reorder = (GMRFLib_reorder_tp) itmp;
			if (err) {
				itmp = GMRFLib_reorder_id((const char *) optarg);
				G.reorder = (GMRFLib_reorder_tp) itmp;
			}
			GMRFLib_reorder = G.reorder;	       /* yes! */
		}
			break;

		case 'R':
		{
			int nrhs = 0;
			err = inla_sread_ints(&nrhs, 1, optarg);
			if (err || nrhs < 0) {
				GMRFLib_ASSERT(err, GMRFLib_EPARAMETER);
				exit(1);
			}
			GMRFLib_pardiso_set_nrhs(nrhs);
		}
			break;

		case 'c':
		{
#if !defined(WINDOWS)
			enable_core_file = 1;		       /* allow for core files */
#endif
		}
			break;

		case 'p':
		{
#if !defined(WINDOWS)
			long int pid = (long int) getpid();
			FILE *fp_pid = fopen(".inla.pid", "w");
			if (fp_pid) {
				fprintf(fp_pid, "%ld\n", pid);
				fclose(fp_pid);
			}
#endif
		}
			break;

		case 'L':
		{
			// link with vecLib on Mac. this is just dummy option
		}
			break;

		default:
			_USAGE;
			exit(EXIT_FAILURE);
		}
	}

#if !defined(INLA_WITH_MKL) && !defined(INLA_WITH_ARMPL)
	// I need to set it here as it depends on MAX_THREADS
	GMRFLib_dot_product_optim_report = Calloc(GMRFLib_CACHE_LEN(), double *);
	for (i = 0; i < GMRFLib_CACHE_LEN(); i++) {
		// GMRFLib_dot_product_optim_report[i] = Calloc(13, double);
		GMRFLib_dot_product_optim_report[i] = Calloc(5, double);
	}
#endif

#if !defined(WINDOWS)
	/*
	 * disable the creation of core-file, unless explicite asked for by the argument '-c'.
	 */
	struct rlimit rlim;
	getrlimit(RLIMIT_CORE, &rlim);
	rlim.rlim_cur = (enable_core_file ? rlim.rlim_max : (rlim_t) 0L);
	setrlimit(RLIMIT_CORE, (const struct rlimit *) &rlim);
	if (0) {
		getrlimit(RLIMIT_CORE, &rlim);
		printf("NEW cur %lld max %lld\n", (long long) rlim.rlim_cur, (long long) rlim.rlim_max);
	}
#endif

	/*
	 * these options does not belong here in this program, but it makes all easier... and its undocumented.
	 */
	switch (G.mode) {
	case INLA_MODE_OPENMP:
	{
		printf("export OMP_NUM_THREADS=%1d,%1d,1,1; ", GMRFLib_openmp->max_threads_nested[0], GMRFLib_openmp->max_threads_nested[1]);
		printf("export OMP_MAX_ACTIVE_LEVELS=%1d; ", (GMRFLib_openmp->max_threads_nested[1] <= 1 ? 1 : 2));
		exit(EXIT_SUCCESS);
	}
		break;

	case INLA_MODE_QINV:
	{
		inla_qinv(argv[optind], argv[optind + 1], argv[optind + 2]);
		exit(EXIT_SUCCESS);
	}
		break;

	case INLA_MODE_QSOLVE:
	{
		inla_qsolve(argv[optind], argv[optind + 1], argv[optind + 2], argv[optind + 3]);
		exit(EXIT_SUCCESS);
	}
		break;

	case INLA_MODE_QREORDERING:
	{
		inla_qreordering(argv[optind]);
		exit(EXIT_SUCCESS);
	}
		break;

	case INLA_MODE_QSAMPLE:
	{
		inla_qsample(argv[optind], argv[optind + 1], argv[optind + 2], argv[optind + 3], argv[optind + 4], argv[optind + 5],
			     argv[optind + 6], argv[optind + 7], argv[optind + 8], argv[optind + 9], verbose);
		exit(EXIT_SUCCESS);
	}
		break;

	case INLA_MODE_FINN:
	{
		inla_finn(argv[optind]);
		exit(EXIT_SUCCESS);
	}
		break;

	case INLA_MODE_GRAPH:
	{
		inla_read_graph(argv[optind]);
		exit(EXIT_SUCCESS);
	}
		break;

	case INLA_MODE_R:
	{
		inla_R(&(argv[optind]));
		exit(EXIT_SUCCESS);
	}
		break;

	case INLA_MODE_FGN:
	{
		inla_fgn(argv[optind], argv[optind + 1]);
		exit(EXIT_SUCCESS);
	}
		break;

	case INLA_MODE_PARDISO:
	{
		inla_check_pardiso();
		exit(EXIT_SUCCESS);
	}
		break;

	case INLA_MODE_TESTIT:
	{
		testit(argc - optind, &(argv[optind]));
		exit(EXIT_SUCCESS);
	}
		break;

	case INLA_MODE_HYPER:
	case INLA_MODE_DEFAULT:
		break;

	default:
		break;
	}

	if (!silent || verbose) {
		fprintf(stdout, "\n\t%s\n", GITCOMMIT);
	}
	if (verbose) {
		_BUGS_intern(stdout);
	}
	if (verbose && G.reorder >= 0) {
		printf("Set reordering to id=[%d] and name=[%s]\n", G.reorder, GMRFLib_reorder_name(G.reorder));
	}
	if (optind >= argc) {
		fprintf(stderr, "\n*** Error: Expected argument after options.\n");
		_USAGE;
		exit(EXIT_FAILURE);
	}
	if (optind < argc - 1) {
		fprintf(stderr, "\n");
		for (i = 0; i < argc; i++) {
			fprintf(stderr, "\targv[%1d] = [%s]\n", i, argv[i]);
		}
		fprintf(stderr, "\targc=[%1d] optind=[%1d]\n\n", argc, optind);
		fprintf(stderr, "\n*** Error: Can only process one .INI-file at the time.\n");
		exit(EXIT_SUCCESS);
	}
	if (verbose) {
		if (G.mode == INLA_MODE_HYPER) {
			fprintf(stderr, "\nRun in mode=[%s]\n", "HYPER");
		}
	}

	if (G.mode == INLA_MODE_DRYRUN) {
		for (arg = optind; arg < argc; arg++) {
			mb = inla_build(argv[arg], verbose, 1);

			char *nndir = NULL;
			FILE *fp = NULL;
			GMRFLib_sprintf(&nndir, "%s/%s", mb->dir, "dryrun");
			fp = fopen(nndir, "w");
			if (!fp) {
				inla_error_open_file(nndir);
			}
			fprintf(fp, "%1d %1d\n", mb->idx_tot, mb->idx_ntot);
			for (i = 0; i < mb->idx_tot; i++) {
				fprintf(fp, "%s %1d %1d\n", mb->idx_tag[i], mb->idx_start[i], mb->idx_n[i]);
			}
			fprintf(fp, "%1d\n", mb->ntheta);
			if (mb->ntheta) {
				for (i = 0; i < mb->ntheta; i++) {
					fprintf(fp, "%s\n", mb->theta_tag[i]);
				}
			}
			fclose(fp);
		}
	}

	if (G.mode == INLA_MODE_DEFAULT || G.mode == INLA_MODE_HYPER) {
		for (arg = optind; arg < argc; arg++) {
			if (verbose) {
				printf("Process file[%s] threads[%1d] max.threads[%1d] blas_threads_force[%1d]",
				       argv[arg], GMRFLib_MAX_THREADS(), host_max_threads, GMRFLib_openmp->blas_num_threads_force);
				if (GMRFLib_openmp->max_threads_nested) {
					printf(" nested[%1d:%1d]\n", GMRFLib_openmp->max_threads_nested[0], GMRFLib_openmp->max_threads_nested[1]);
				} else {
					printf("\n");
				}
			}
			if (!silent) {
				printf("\nWall-clock time used on [%s] max_threads=[%1d]\n", argv[arg], GMRFLib_MAX_THREADS());
			}
			time_used[0] = GMRFLib_timer();
			atime_used[0] = clock();

			GMRFLib_openmp_implement_strategy(GMRFLib_OPENMP_PLACES_PARSE_MODEL, NULL, NULL);
			mb = inla_build(argv[arg], verbose, 1);
			GMRFLib_openmp_implement_strategy(GMRFLib_OPENMP_PLACES_BUILD_MODEL, NULL, NULL);
			time_used[0] = GMRFLib_timer() - time_used[0];
			atime_used[0] = clock() - atime_used[0];
			if (!silent) {
				printf("\tPreparations             : %7.3f seconds\n", time_used[0]);
				fflush(stdout);
			}
			time_used[1] = GMRFLib_timer();
			atime_used[1] = clock();

			int nfunc[2] = { 0, 0 };
			double rgeneric_cpu[2] = { 0.0, 0.0 };

			if (GMRFLib_inla_mode == GMRFLib_MODE_COMPACT) {
				time_used[3] = GMRFLib_timer();
				inla_INLA_preopt_experimental(mb);
				time_used[3] = GMRFLib_timer() - time_used[1];
				atime_used[3] = clock() - atime_used[1];
				nfunc[0] = mb->misc_output->nfunc;
				rgeneric_cpu[0] = R_rgeneric_cputime;

				Free(mb->theta_file);
				Free(mb->x_file);
				if (mb->preopt->mode_theta) {
					mb->theta_file = Calloc(mb->ntheta, double);
					Memcpy(mb->theta_file, mb->preopt->mode_theta, mb->ntheta * sizeof(double));
				} else {
					mb->theta_file = NULL;
				}
				mb->x_file = Calloc(mb->preopt->n + mb->preopt->mnpred, double);
				Memcpy(mb->x_file, mb->preopt->mode_x, (mb->preopt->n + mb->preopt->mnpred) * sizeof(double));
				for (i = 0; i < mb->preopt->mnpred; i++) {
					mb->x_file[i] += OFFSET3(i);
				}

				mb->ntheta_file = mb->ntheta;
				mb->nx_file = mb->preopt->n + mb->preopt->mnpred;
				mb->mode_restart = mb->ai_par->mode_restart = GMRFLib_TRUE;
				inla_reset();
			} else if (GMRFLib_inla_mode == GMRFLib_MODE_CLASSIC) {
				GMRFLib_inla_mode = GMRFLib_MODE_CLASSIC;
				inla_INLA(mb);
				nfunc[0] = mb->misc_output->nfunc;
			} else {
				assert(0 == 1);
			}

			time_used[1] = GMRFLib_timer() - time_used[1];
			atime_used[1] = clock() - atime_used[1];
			if (!silent) {
				if (GMRFLib_inla_mode == GMRFLib_MODE_CLASSIC) {
					printf("\tApprox inference         : %7.3f seconds\n", time_used[1]);
				} else {
					printf("\tApprox inference (stage1): %7.3f seconds\n", time_used[3]);
					printf("\tApprox inference (stage2): %7.3f seconds\n", time_used[1] - time_used[3]);
					printf("\tApprox inference (total) : %7.3f seconds\n", time_used[1]);
				}
				fflush(stdout);
			}
			time_used[2] = GMRFLib_timer();
			atime_used[2] = clock();
			GMRFLib_openmp_implement_strategy(GMRFLib_OPENMP_PLACES_DEFAULT, NULL, NULL);
			inla_output(mb);
			time_used[2] = GMRFLib_timer() - time_used[2];
			atime_used[2] = clock() - atime_used[2];

#define PEFF_OUTPUT(fp_)						\
			if (1) {					\
				char *tab = Strdup("");			\
				if (fp_ == stdout) tab = Strdup("\t");	\
				fprintf(fp_, "Total:");			\
				eff_nt = ((double)(atime_used[0] + atime_used[1]))/CLOCKS_PER_SEC/(time_used[0] + time_used[1]); \
				fprintf(fp_, "%sAccumulated CPU-time is equivalent to %.2f threads running at 100%%\n", tab, eff_nt); \
				fprintf(fp_, "%sEfficiency using %1d threads = %.2f%%\n", tab, GMRFLib_MAX_THREADS(), \
					100.0 * eff_nt/GMRFLib_MAX_THREADS()); \
			}
#define PEFF_PREOPT_OUTPUT(fp_)						\
			if (1) {					\
				char *tab = Strdup("");			\
				if (fp_ == stdout) tab = Strdup("\t");	\
				fprintf(fp_, "Stage1:");		\
				eff_nt = ((double)(atime_used[0] + atime_used[3]))/CLOCKS_PER_SEC/(time_used[0] + time_used[3]); \
				fprintf(fp_,"%sAccumulated CPU-time is equivalent to %.2f threads running at 100%%\n", tab, eff_nt); \
				fprintf(fp_,"%sEfficiency using %1d threads = %.2f%%\n", tab, GMRFLib_MAX_THREADS(), \
					100.0 * eff_nt/GMRFLib_MAX_THREADS()); \
				fprintf(fp_,"Stage2:");			\
				eff_nt = ((double)(atime_used[0] + atime_used[1] - atime_used[3]))/CLOCKS_PER_SEC/(time_used[0] + time_used[1] - time_used[3]); \
				fprintf(fp_,"%sAccumulated CPU-time is equivalent to %.2f threads running at 100%%\n", tab, eff_nt); \
				fprintf(fp_,"%sEfficiency using %1d threads = %.2f%%\n", tab, GMRFLib_MAX_THREADS(), \
					100.0 * eff_nt/GMRFLib_MAX_THREADS()); \
			}

			if (!silent) {
				printf("\tOutput                   : %7.3f seconds\n", time_used[2]);
				printf("\t------------------------------------------\n");
				printf("\tTotal                    : %7.3f seconds\n\n", time_used[0] + time_used[1] + time_used[2]);
#if !defined(WINDOWS)
				if (GMRFLib_inla_mode != GMRFLib_MODE_CLASSIC) {
					PEFF_PREOPT_OUTPUT(stdout);
				}
				PEFF_OUTPUT(stdout);
#endif
				fflush(stdout);
			}
			if (verbose) {
				printf("\nWall-clock time used on [%s]\n", argv[arg]);
				printf("\tPreparations             : %7.3f seconds\n", time_used[0]);
				if (GMRFLib_inla_mode == GMRFLib_MODE_CLASSIC) {
					printf("\tApprox inference         : %7.3f seconds\n", time_used[1]);
				} else {
					printf("\tApprox inference (stage1): %7.3f seconds\n", time_used[3]);
					printf("\tApprox inference (stage2): %7.3f seconds\n", time_used[1] - time_used[3]);
					printf("\tApprox inference (total) : %7.3f seconds\n", time_used[1]);
				}
				printf("\tOutput                   : %7.3f seconds\n", time_used[2]);
				printf("\t------------------------------------------\n");
				printf("\tTotal                    : %7.3f seconds\n\n", time_used[0] + time_used[1] + time_used[2]);

				if (GMRFLib_inla_mode == GMRFLib_MODE_CLASSIC) {
					printf("\nNumber of fn-calls= %1d with %.4f sec/fn-call\n",
					       mb->misc_output->nfunc, time_used[1] / IMAX(1, mb->misc_output->nfunc));
					if (R_rgeneric_cputime > 0.0) {
						printf("rgeneric-time= %.4f seconds, with %.4f sec/fn-call and %.4f%% of the total time\n",
						       R_rgeneric_cputime,
						       R_rgeneric_cputime / IMAX(1, mb->misc_output->nfunc),
						       R_rgeneric_cputime / time_used[1] * 100.0);
					}
				} else {
					printf("Stage1:");
					printf("\tNumber of fn-calls= %1d with %.4f sec/fn-call\n", nfunc[0], time_used[3] / IMAX(1, nfunc[0]));
					if (rgeneric_cpu[0] > 0.0) {
						printf("\trgeneric-time= %.4f seconds, with %.4f sec/fn-call and %.4f%% of the total time\n",
						       rgeneric_cpu[0],
						       rgeneric_cpu[0] / IMAX(1, nfunc[0]), rgeneric_cpu[0] / time_used[3] * 100.0);
					}
					printf("Stage2:");
					printf("\tNumber of fn-calls= %1d with %.4f sec/fn-call\n", nfunc[1], time_used[1] / IMAX(1, nfunc[1]));
					if (rgeneric_cpu[1] > 0.0) {
						printf("\trgeneric-time= %.4f seconds, with %.4f sec/fn-call and %.4f%% of the total time\n",
						       rgeneric_cpu[1],
						       rgeneric_cpu[1] / IMAX(1, nfunc[1]),
						       rgeneric_cpu[1] / (time_used[1] - time_used[3]) * 100.0);
					}
				}
				if (GMRFLib_dot_product_gain > 0.0) {
					printf("\nDot-product gain: %.3f seconds, %.6f seconds/fn-call\n\n", GMRFLib_dot_product_gain,
					       GMRFLib_dot_product_gain / nfunc[0]);
				}
#if !defined(WINDOWS)
				if (GMRFLib_inla_mode != GMRFLib_MODE_CLASSIC) {
					PEFF_PREOPT_OUTPUT(stdout);
				}
				PEFF_OUTPUT(stdout);
#endif
				printf("\n");
			}

			if (mb->dir) {
				// just a copy of what is above
				char *nfile = NULL;
				GMRFLib_sprintf(&nfile, "%s/cpu-intern", mb->dir);
				FILE *fp = fopen(nfile, "w");
				if (fp) {
					fprintf(fp, "Wall-clock time used on [%s]\n", argv[arg]);
					fprintf(fp, "Preparations             : %7.3f seconds\n", time_used[0]);
					if (GMRFLib_inla_mode == GMRFLib_MODE_CLASSIC) {
						fprintf(fp, "Approx inference         : %7.3f seconds\n", time_used[1]);
					} else {
						fprintf(fp, "Approx inference (stage1): %7.3f seconds\n", time_used[3]);
						fprintf(fp, "Approx inference (stage2): %7.3f seconds\n", time_used[1] - time_used[3]);
						fprintf(fp, "Approx inference (total) : %7.3f seconds\n", time_used[1]);
					}
					fprintf(fp, "Output                   : %7.3f seconds\n", time_used[2]);
					fprintf(fp, "------------------------------------------\n");
					fprintf(fp, "Total                    : %7.3f seconds\n", time_used[0] + time_used[1] + time_used[2]);
				}

				if (GMRFLib_inla_mode == GMRFLib_MODE_CLASSIC) {
					fprintf(fp, "Number of fn-calls= %1d with %.4f sec/fn-call\n",
						mb->misc_output->nfunc, time_used[1] / IMAX(1, mb->misc_output->nfunc));
					if (R_rgeneric_cputime > 0.0) {
						fprintf(fp, "rgeneric-time= %.4f seconds, with %.4f sec/fn-call and %.4f%% of the total time\n",
							R_rgeneric_cputime,
							R_rgeneric_cputime / IMAX(1, mb->misc_output->nfunc),
							R_rgeneric_cputime / time_used[1] * 100.0);
					}
				} else {
					fprintf(fp, "Stage1:");
					fprintf(fp, "Number of fn-calls= %1d with %.4f sec/fn-call\n", nfunc[0], time_used[3] / IMAX(1, nfunc[0]));
					if (rgeneric_cpu[0] > 0.0) {
						fprintf(fp, "rgeneric-time= %.4f seconds, with %.4f sec/fn-call and %.4f%% of the total time\n",
							rgeneric_cpu[0],
							rgeneric_cpu[0] / IMAX(1, nfunc[0]), rgeneric_cpu[0] / time_used[3] * 100.0);
					}
					fprintf(fp, "Stage2:");
					fprintf(fp, "Number of fn-calls= %1d with %.4f sec/fn-call\n", nfunc[1], time_used[1] / IMAX(1, nfunc[1]));
					if (rgeneric_cpu[1] > 0.0) {
						fprintf(fp, "rgeneric-time= %.4f seconds, with %.4f sec/fn-call and %.4f%% of the total time\n",
							rgeneric_cpu[1],
							rgeneric_cpu[1] / IMAX(1, nfunc[1]),
							rgeneric_cpu[1] / (time_used[1] - time_used[3]) * 100.0);
					}
				}
#if !defined(WINDOWS)
				if (GMRFLib_inla_mode != GMRFLib_MODE_CLASSIC) {
					PEFF_PREOPT_OUTPUT(fp);
				}
				PEFF_OUTPUT(fp);
#endif
				fclose(fp);
				Free(nfile);
			}
		}
	}

	if (mb) {
		inla_output_ok(mb->dir);
	}

	return EXIT_SUCCESS;
#undef _USAGE_intern
#undef _USAGE
#undef _HELP
#undef _BUGS_intern
#undef _BUGS
}
