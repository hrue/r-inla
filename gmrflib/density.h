
/* density.h
 * 
 * Copyright (C) 2006 Havard Rue
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
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
 *
 */

/*!
  \file density.h
  \brief Typedefs used to represent densities
*/

#ifndef __GMRFLib_DENSITY_H__
#define __GMRFLib_DENSITY_H__

#include <stdlib.h>
#include <stddef.h>
#include <math.h>
#if !defined(__FreeBSD__)
#include <malloc.h>
#endif

#undef __BEGIN_DECLS
#undef __END_DECLS
#ifdef __cplusplus
#define __BEGIN_DECLS extern "C" {
#define __END_DECLS }
#else
#define __BEGIN_DECLS					       /* empty */
#define __END_DECLS					       /* empty */
#endif

__BEGIN_DECLS

/*
 * do not integrate beyond +- GMRFLib_DENSITY_INTEGRATION_LIMIT stdev's, unless requested to do so.
 */
#define GMRFLib_DENSITY_INTEGRATION_LIMIT (10.0)

/* 
 *  length of work the GSL-integration routine
 */
#define GMRFLib_DENSITY_LENGTH_WORK (2048)

/* 
 *  this is a wrapper around the gsl-integration routine, so that if the relative high accurancy FAILS due to roundoff errors,
 *  or similar, the increase the 'eps' with a factor of 10 and try again until success. we need in general to use the _qags
 *  routine which can cope with 'nice' singularities, and further, we split the region into subregions as well; sometimes this
 *  is needed, epspecially if the density is peaked and low and high are far away...
 */
#define GMRFLib_gsl_integration_wrapper_qag(a0, a1, a2, a3, a4, a5, a6, a7, a8) \
	if (1) {							\
		int ssstatus;						\
		gsl_error_handler_t *ehandler = gsl_set_error_handler_off(); /* turn off the error-handler */ \
		do							\
		{							\
			ssstatus = gsl_integration_qag(a0, a1, a2, a3, a4, a5, GSL_INTEG_GAUSS41, a6, a7, a8); \
			if (ssstatus) eps *= 10.0;			\
		} while(ssstatus == GSL_ETOL || ssstatus == GSL_ELOSS || ssstatus == GSL_EROUND); \
		gsl_set_error_handler(ehandler);	  /*  turn it on again */ \
	}
#define GMRFLib_gsl_integration_wrapper_qags(a0, a1, a2, a3, a4, a5, a6, a7, a8) \
	if (1) {							\
		int ssstatus;						\
		gsl_error_handler_t *ehandler = gsl_set_error_handler_off(); /* turn off the error-handler */ \
		do							\
		{							\
			ssstatus = gsl_integration_qags(a0, a1, a2, a3, a4, a5, a6, a7, a8); \
			if (ssstatus) eps *= 10.0;			\
		} while(ssstatus == GSL_ETOL || ssstatus == GSL_ELOSS || ssstatus == GSL_EROUND); \
		gsl_set_error_handler(ehandler);	  /*  turn it on again */ \
	}
#define GMRFLib_gsl_integration_wrapper(My_F, My_lower, My_upper, My_epsrel, My_epsabs, My_result, My_error) \
	if (1) {							\
		int My_i, My_ndiv = 3;					\
		double My_newlower, My_newupper, My_res, My_err, My_step; \
									\
		GMRFLib_gsl_integration_fix_limits(&My_newlower, &My_newupper, My_F, My_lower, My_upper); \
		My_step = (My_newupper-My_newlower)/(double)My_ndiv;	\
		*(My_result) = *(My_error) = 0.0;			\
									\
		for(My_i=0; My_i < My_ndiv; My_i++){			\
			gsl_integration_workspace *My_work = NULL;	\
			My_work = gsl_integration_workspace_alloc(GMRFLib_DENSITY_LENGTH_WORK); \
			GMRFLib_gsl_integration_wrapper_qag(My_F, My_newlower + My_i * My_step, \
							    My_newlower + (My_i + 1.0) * My_step, My_epsrel, My_epsabs, \
							    GMRFLib_DENSITY_LENGTH_WORK, My_work, &My_res, &My_err); \
			gsl_integration_workspace_free(My_work);	\
			*(My_result) += My_res;				\
			*(My_error) += My_err;				\
		}							\
	}

/**
 *   \brief Available densities.
 */
    typedef enum {

	/**
	 * \brief Density is Gaussian
	 */
	GMRFLib_DENSITY_TYPE_GAUSSIAN = 1,

	/**
	 * \brief Density is Skew-Normal
	 */
	GMRFLib_DENSITY_TYPE_SKEWNORMAL = 2,

	/**
	 * \brief Density is Spline-corrected Normal
	 */
	GMRFLib_DENSITY_TYPE_SCGAUSSIAN = 3
} GMRFLib_density_type_tp;

/* 
   parameters for the skew-normal, where

   Y ~ SN(xi, omega^2, alpha)

   and (Y-xi)/omega is standard skew-normal with density

   2 \phi(x) \Phi(alpha*x)
*/
typedef struct {
	double xi;
	double omega;
	double alpha;
} GMRFLib_sn_param_tp;

typedef struct {
	size_t n;					       /* number of data */
	size_t m;					       /* number of parameters */
	double *x;					       /* x-values */
	double *log_density;				       /* log_densities */
} GMRFLib_sn_fit_data_tp;


/* 
   the scaling of the weight for SN-fit, stdev = exp(- log_density[i] * GMRFLib_SN_WEIGHT_SCALING )
*/
#define GMRFLib_SN_WEIGHT_SCALING (0.5)

typedef struct {
	double xmin;
	double xmax;
	gsl_interp_accel *accel;
	gsl_spline *spline;
} GMRFLib_spline_tp;

/**
 * \brief The density-object
 */
typedef struct {

	/**
	 * \brief The mean in standarised scale.
	 */
	double mean;

	/**
	 * \brief The standard deviation in standardised scale.
	 */
	double stdev;

	/**
	 * \brief The standardised skewness
	 */
	double skewness;

	/**
	 * \brief The mean in users own scale (not standardised).
	 */
	double user_mean;

	/**
	 * \brief The standard deviation in users own scale (not standardised).
	 */
	double user_stdev;

	/** 
	 * \brief The mode in user scale (set if !NAN)
	 */
	double user_mode;

	/**
	 * \brief The offset that is used for standarisation.
	 *
	 * The density is for the random variable (x-m)/sd where \c m is \c GMRFLib_density_tp::std_mean and \c sd is
	 * GMRFLib_density_tp::std_stdev.
	 */
	double std_mean;

	/** 
	 * \brief The scale that is used for standarisation.
	 *
	 * The density is for the random variable (x-m)/sd where \c m is \c GMRFLib_density_tp::std_mean and \c sd is
	 * GMRFLib_density_tp::std_stdev.
	 */
	double std_stdev;				       /* this is the the standarisation */

	/*
	 * The rest of the paramers is for internal use only. 
	 */

	/*
	 * the type of the density:
	 */
	GMRFLib_density_type_tp type;

	double x_min, x_max;				       /* range for the log_correction, ALSO used by the others */

	/*
	 * params for the GMRFLib_DENSITY_TYPE_GAUSSIAN 
	 */
	double mean_gaussian;
	double stdev_gaussian;

	/*
	 * params for the GMRFLib_DENSITY_TYPE_SKEWNORMAL 
	 */
	GMRFLib_sn_param_tp *sn_param;

	/*
	 * params for the GMRFLib_DENSITY_TYPE_SCGAUSSIAN 
	 */
	double log_norm_const;				       /* log(norm_const), divide by norm_const to get the normalised
							        * density.  */
	GMRFLib_spline_tp *log_correction;

	GMRFLib_spline_tp *P;
	GMRFLib_spline_tp *Pinv;

	GMRFLib_uchar flags;				       /* can set various flags */

} GMRFLib_density_tp;

typedef enum {
	DENSITY_FLAGS_FAILURE = 0
} GMRFLib_density_flag_tp;

typedef enum {

	/**
	 * \brief Default storage strategy
	 */
	GMRFLib_DENSITY_STORAGE_STRATEGY_DEFAULT = 0,

	/**
	 * \brief Low storage strategy
	 */
	GMRFLib_DENSITY_STORAGE_STRATEGY_LOW,

	/**
	 * \brief High storage strategy
	 */
	GMRFLib_DENSITY_STORAGE_STRATEGY_HIGH
} GMRFLib_density_storage_strategy_tp;


/* 
   this object is used for computing properties of the density
*/
typedef struct {
	GMRFLib_density_tp *density;			       /* the density itself */
	int power;					       /* power; for moments calculations */
	double prob;					       /* probability; for solving P(X<x)=prob */
} GMRFLib_density_properties_tp;

const gsl_interp_type *GMRFLib_density_interp_type(int n);
double GMRFLib_density_Pinv_df(double x, void *param);
double GMRFLib_density_Pinv_f(double x, void *param);
double GMRFLib_density_std2user(double x, GMRFLib_density_tp * density);
double GMRFLib_density_std2user_n(double *x_user, double *x, int n, GMRFLib_density_tp * density);
double GMRFLib_density_user2std(double x, GMRFLib_density_tp * density);
double GMRFLib_evaluate_density__intern(double x, void *param);
double GMRFLib_evaluate_density_kld2__intern(double x, void *param);
double GMRFLib_evaluate_density_kld__intern(double x, void *param);
double GMRFLib_evaluate_density_power__intern(double x, void *param);
double GMRFLib_log_gsl_cdf_ugaussian_P(double z);
double GMRFLib_sn_logdensity_diff_alpha(double x, void *param);
double GMRFLib_sn_logdensity_diff_omega(double x, void *param);
double GMRFLib_sn_logdensity_diff_xi(double x, void *param);
int GMRFLib_density_P(double *px, double x, GMRFLib_density_tp * density);
int GMRFLib_density_Pinv(double *xp, double p, GMRFLib_density_tp * density);
int GMRFLib_density_adjust_vector(double *ldens, int n);
int GMRFLib_density_combine(GMRFLib_density_tp ** density, GMRFLib_density_tp ** gdensity, int n, GMRFLib_density_tp ** densities,
			    double *weights);
int GMRFLib_density_create(GMRFLib_density_tp ** density, int type, int n, double *x, double *logdens, double std_mean,
			   double std_stdev, int lookup_tables);
int GMRFLib_density_create_normal(GMRFLib_density_tp ** density, double mean, double stdev, double std_mean, double std_stdev);
int GMRFLib_density_create_sn(GMRFLib_density_tp ** density, GMRFLib_sn_param_tp sn_param, double std_mean, double std_stdev,
			      int lookup_tables);
int GMRFLib_density_duplicate(GMRFLib_density_tp ** density_to, GMRFLib_density_tp * density_from);
int GMRFLib_density_layout_x(double **x_vec, int *len_x, GMRFLib_density_tp * density);
int GMRFLib_density_new_mean(GMRFLib_density_tp ** new_density, GMRFLib_density_tp * density, double new_mean);
int GMRFLib_density_printf(FILE * fp, GMRFLib_density_tp * density);
int GMRFLib_density_prune_weights(int *n_idx, int *idx, double *weights, int n);
int GMRFLib_density_user2std_n(double *x_std, double *x, GMRFLib_density_tp * density, int n);
int GMRFLib_evaluate_densities(double *dens, double x_user, int n, GMRFLib_density_tp ** densities, double *weights);
int GMRFLib_evaluate_density(double *dens, double x, GMRFLib_density_tp * density);
int GMRFLib_evaluate_gdensities(double *dens, double x_user, int n, GMRFLib_density_tp ** densities, double *weights);
int GMRFLib_evaluate_logdensity(double *logdens, double x, GMRFLib_density_tp * density);
int GMRFLib_evaluate_ndensities(double *dens, int nd, double *x_user, int nx, GMRFLib_density_tp ** densities, double *weights);
int GMRFLib_evaluate_ndensity(double *dens, double *x, int n, GMRFLib_density_tp * density);
int GMRFLib_evaluate_nlogdensity(double *logdens, double *x, int n, GMRFLib_density_tp * density);
int GMRFLib_free_density(GMRFLib_density_tp * density);
int GMRFLib_gsl_integration_fix_limits(double *new_lower, double *new_upper, gsl_function * F, double lower, double upper);
int GMRFLib_init_density(GMRFLib_density_tp * density, int lookup_tables);
int GMRFLib_kld(double *kld, GMRFLib_density_tp * density, GMRFLib_density_tp * ddensity);
int GMRFLib_kld_sym(double *kld_sym, GMRFLib_density_tp * density, GMRFLib_density_tp * ddensity);
int GMRFLib_mkld(double *mkld, GMRFLib_density_tp * density, GMRFLib_density_tp * ddensity);
int GMRFLib_mkld_sym(double *mkld_sym, GMRFLib_density_tp * density, GMRFLib_density_tp * ddensity);
int GMRFLib_normal_fit(double *mean, double *variance, double *fval, double *x, double *log_density, int n);
int GMRFLib_sn_density(double *dens, double x, void *param);
int GMRFLib_sn_fit(GMRFLib_sn_param_tp * param, double *fval, double *x, double *log_density, int n);
int GMRFLib_sn_fit__intern(void *param, double *fval, double *x, double *log_density, size_t n, size_t m);
int GMRFLib_sn_fit_df(const gsl_vector * param, void *data, gsl_matrix * J);
int GMRFLib_sn_fit_f(const gsl_vector * param, void *data, gsl_vector * f);
int GMRFLib_sn_fit_fdf(const gsl_vector * param, void *data, gsl_vector * f, gsl_matrix * J);
int GMRFLib_sn_logdensity(double *ldens, double x, void *param);
int GMRFLib_sn_moments(double *mean, double *stdev, double *skewness, GMRFLib_sn_param_tp * p);
void GMRFLib_density_Pinv_fdf(double x, void *param, double *f, double *df);

__END_DECLS
#endif
