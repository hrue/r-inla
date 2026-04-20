#ifndef __GMRFLib_AUX_H__
#define __GMRFLib_AUX_H__

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
  To represent a log_gamma as a mixture of normals
*/
    typedef struct {
	double weights[10];
	double means[10];
	double variances[10];
	int ncomp;
} GMRFLib_lgamma_mixture_tp;

/* 
   model y ~ Po( E * exp(x) )
*/
typedef struct {
	double *y;					       /* observed data, y */
	double *E;					       /* scaling, y ~ Po(E*exp(x)) */
	double *x;					       /* x as above */
	double tau1, tau2;				       /* inter-arrival times */
	GMRFLib_lgamma_mixture_tp *mix1, *mix2;		       /* the normal mixtures */
	int r1, r2;					       /* mixture indicators */
} GMRFLib_poisson_aux_tp;

/* 
   model y ~ Binomial(nb, p), where logit(p) = x
 */
typedef struct {
	double *y;					       /* y ~ Bin(nb, ...) */
	double *nb;					       /* nb in the Bin(...) */
	double y_star;					       /* the aggregated utility */
	double *x;					       /* x as above */
	int r;						       /* the component indicator */
	GMRFLib_lgamma_mixture_tp *mix;			       /* the normal mixtures */
} GMRFLib_binomial_aux_tp;

/* 
   more aux-schemes can be added here
*/
typedef struct {
	GMRFLib_poisson_aux_tp *poisson;
	GMRFLib_binomial_aux_tp *binomial;
} GMRFLib_aux_tp;

/* 
   this hold them all
 */
typedef struct {
	int n;
	double *d;
	GMRFLib_aux_tp **aux;
} GMRFLib_auxs_tp;

int GMRFLib_aux_free(GMRFLib_aux_tp * aux);
int GMRFLib_aux_free_all(GMRFLib_auxs_tp * auxs);
int GMRFLib_aux_gauss_approx(double *b, double *c, GMRFLib_aux_tp * aux);
int GMRFLib_aux_gauss_approx_all(double *b, double *c, GMRFLib_auxs_tp * auxs);
int GMRFLib_aux_gauss_approx_binomial(double *b, double *c, GMRFLib_binomial_aux_tp * binomial);
int GMRFLib_aux_gauss_approx_poisson(double *b, double *c, GMRFLib_poisson_aux_tp * poisson);
int GMRFLib_aux_init(GMRFLib_aux_tp * aux);
int GMRFLib_aux_init_all(GMRFLib_auxs_tp * auxs);
int GMRFLib_aux_init_binomial(GMRFLib_binomial_aux_tp * binomial);
int GMRFLib_aux_init_poisson(GMRFLib_poisson_aux_tp * poisson);
int GMRFLib_aux_setup_binomial(GMRFLib_aux_tp ** aux, double *y, double *nb, double *x);
int GMRFLib_aux_setup_binomial_all(GMRFLib_auxs_tp ** auxs, int n, double *d, double *y, double *nb, double *x);
int GMRFLib_aux_setup_poisson(GMRFLib_aux_tp ** aux, double *y, double *E, double *x);
int GMRFLib_aux_setup_poisson_all(GMRFLib_auxs_tp ** auxs, int n, double *d, double *y, double *E, double *x);
int GMRFLib_aux_update(GMRFLib_aux_tp * aux);
int GMRFLib_aux_update_all(GMRFLib_auxs_tp * auxs);
int GMRFLib_aux_update_binomial(GMRFLib_binomial_aux_tp * binomial);
int GMRFLib_aux_update_poisson(GMRFLib_poisson_aux_tp * poisson);
int GMRFLib_mixture_lgamma(GMRFLib_lgamma_mixture_tp ** mixture, double n);
int GMRFLib_mixture_lgamma_testing__intern(void);
int GMRFLib_mixture_update_indicator(int *r, double tau, double lambda, GMRFLib_lgamma_mixture_tp * mix);

__END_DECLS
#endif
