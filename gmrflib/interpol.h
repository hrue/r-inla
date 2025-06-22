#ifndef __GMRFLib_INTERPOL_H__
#define __GMRFLib_INTERPOL_H__
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
 *
 */
    typedef enum {
	GMRFLib_INTPOL_TRANS_NONE = 0,
	GMRFLib_INTPOL_TRANS_P,
	GMRFLib_INTPOL_TRANS_Pinv
} GMRFLib_intpol_transform_tp;

typedef enum {
	GMRFLib_INTPOL_CACHE_LEVEL12 = 0,		       /* one for both levels */
	GMRFLib_INTPOL_CACHE_LEVEL1 = 1,		       /* level 1 only */
	GMRFLib_INTPOL_CACHE_SIMPLE = 2,		       /* serial, not thread-safe */
	GMRFLib_INTPOL_CACHE_NONE = 3			       /* none */
} GMRFLib_intpol_cache_tp;


#define GMRFLib_SN_SKEWMAX (0.988)
typedef struct {
	GMRFLib_intpol_transform_tp trans;
	GMRFLib_intpol_cache_tp cache;
	int cache_len;
	double xmin;
	double xmax;
	gsl_interp_accel **accel;
	gsl_spline *spline;
} GMRFLib_spline_tp;

GMRFLib_spline_tp *GMRFLib_spline_create(double *x, double *y, int n);
GMRFLib_spline_tp *GMRFLib_spline_create_x(double *x, double *y, int n, GMRFLib_intpol_transform_tp trans, GMRFLib_intpol_cache_tp cache,
					   int skip_checks);
GMRFLib_spline_tp *GMRFLib_spline_create_from_matrix(GMRFLib_matrix_tp * M);
double GMRFLib_spline_eval(double x, GMRFLib_spline_tp * s);
double GMRFLib_spline_eval_deriv(double x, GMRFLib_spline_tp * s);
double GMRFLib_spline_eval_deriv_x(double x, GMRFLib_spline_tp * s);
double GMRFLib_spline_eval_deriv2(double x, GMRFLib_spline_tp * s);
double GMRFLib_spline_eval_deriv2_x(double x, GMRFLib_spline_tp * s);
int GMRFLib_spline_free(GMRFLib_spline_tp * s);
int GMRFLib_spline_eval_x(int n, double *x, GMRFLib_spline_tp * s, double *values);


__END_DECLS
#endif
