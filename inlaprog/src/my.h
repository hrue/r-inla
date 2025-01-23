#ifndef __INLA_MY_H__
#define __INLA_MY_H__
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
double *my_compute_lbell(int nmax);
double my_betabinomial(int y, int n, double a, double b, double *work);
double my_betabinomial2(int y, int n, double a, double b);
double my_betabinomial_helper(int n, double a, double *work);
double my_betabinomial_helper4(int n, double a, double *work);
double my_betabinomial_helper8(int n, double a, double *work);
double my_betabinomial_helper16(int n, double a, double *work);
double my_betabinomial_helper_core(int n, double a, double *work, int roll);
double my_gsl_sf_lnbeta(double a, double b);
double my_gsl_sf_lnchoose(unsigned int n, unsigned int m);
double my_gsl_sf_lnfact(int x);
double my_gsl_sf_lngamma(double x);
double my_lambert_W0(double y);
double my_lbell(int y);
int my_dir_exists(const char *dirname);
int my_file_exists(const char *filename);
int my_gsl_sf_lnchoose_e(unsigned int n, unsigned int m, gsl_sf_result * result);
int my_gsl_sf_lnfact_e(const unsigned int n, gsl_sf_result * result);
int my_setenv(char *str, int prefix);
void my_lambert_W0s(int m, double *y, double *res);

__END_DECLS
#endif
