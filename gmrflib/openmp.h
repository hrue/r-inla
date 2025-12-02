#ifndef __GMRFLib_OPENMP_H__
#define __GMRFLib_OPENMP_H__
#include <assert.h>

#undef __BEGIN_DECLS
#undef __END_DECLS
#ifdef __cplusplus
#define __BEGIN_DECLS  extern "C" {
#define __END_DECLS }
#else
#define __BEGIN_DECLS					       /* empty */
#define __END_DECLS					       /* empty */
#endif


// workaround for the moment (gcc-15.2.1)

#if defined(__cplusplus)
# define my_tmp__ __cplusplus
# undef __cplusplus
# include <omp.h> 
# define __cplusplus my_tmp__
#else
# include <omp.h> 
#endif

__BEGIN_DECLS

typedef struct {
	char *tag;
	int max_nt;
	int best_nt;
	int min_num_try;
	int try_nt;
	int done;
	int step;
	double tot_times;
	double *ntimes;
	double *acc_wtime;
} GMRFLib_adapt_nt_tp;

typedef enum {
	GMRFLib_OPENMP_STRATEGY_SMALL = 1,
	GMRFLib_OPENMP_STRATEGY_MEDIUM,
	GMRFLib_OPENMP_STRATEGY_LARGE,
	GMRFLib_OPENMP_STRATEGY_HUGE,
	GMRFLib_OPENMP_STRATEGY_PARDISO,
	GMRFLib_OPENMP_STRATEGY_STILES,
	GMRFLib_OPENMP_STRATEGY_DEFAULT,
	GMRFLib_OPENMP_STRATEGY_NONE
} GMRFLib_openmp_strategy_tp;

#define GMRFLib_OPENMP_STRATEGY_NAME(num)				\
	((num) == GMRFLib_OPENMP_STRATEGY_SMALL ? "small" :		\
	 ((num) == GMRFLib_OPENMP_STRATEGY_MEDIUM ? "medium" :		\
	  ((num) == GMRFLib_OPENMP_STRATEGY_LARGE ? "large" :		\
	   ((num) == GMRFLib_OPENMP_STRATEGY_HUGE ? "huge" :		\
	    ((num) == GMRFLib_OPENMP_STRATEGY_PARDISO ? "pardiso" :	\
	     ((num) == GMRFLib_OPENMP_STRATEGY_STILES ? "sTiles" :	\
	      ((num) == GMRFLib_OPENMP_STRATEGY_DEFAULT ? "default" :	\
	       ((num) == GMRFLib_OPENMP_STRATEGY_NONE ? "none" : "THIS SHOULD NOT HAPPEN"))))))))

typedef enum {
	GMRFLib_OPENMP_PLACES_PARSE_MODEL = 1,
	GMRFLib_OPENMP_PLACES_BUILD_MODEL,
	GMRFLib_OPENMP_PLACES_BUILD_MODEL2,
	GMRFLib_OPENMP_PLACES_OPTIMIZE,
	GMRFLib_OPENMP_PLACES_HESSIAN,
	GMRFLib_OPENMP_PLACES_HESSIAN_SCALE,
	GMRFLib_OPENMP_PLACES_INTEGRATE_HYPERPAR,
	GMRFLib_OPENMP_PLACES_COMBINE,
	GMRFLib_OPENMP_PLACES_EXTERNAL,
	GMRFLib_OPENMP_PLACES_TIMING,
	GMRFLib_OPENMP_PLACES_GCPO_BUILD,
	GMRFLib_OPENMP_PLACES_SPECIAL,
	GMRFLib_OPENMP_PLACES_DEFAULT,
	GMRFLib_OPENMP_PLACES_NONE
} GMRFLib_openmp_place_tp;

#define GMRFLib_OPENMP_PLACE_NAME(num)					\
	(num == GMRFLib_OPENMP_PLACES_PARSE_MODEL ? "parse.model" :	\
		(num == GMRFLib_OPENMP_PLACES_BUILD_MODEL ? "build.model" : \
		 (num == GMRFLib_OPENMP_PLACES_BUILD_MODEL2 ? "build.model2" : \
		  (num == GMRFLib_OPENMP_PLACES_OPTIMIZE ? "optimize" :	\
		   (num == GMRFLib_OPENMP_PLACES_HESSIAN ? "hessian" :	\
		    (num == GMRFLib_OPENMP_PLACES_HESSIAN_SCALE ? "hessian.scale" : \
		     (num == GMRFLib_OPENMP_PLACES_INTEGRATE_HYPERPAR ? "integrate.hyperpar" : \
		      (num == GMRFLib_OPENMP_PLACES_COMBINE ? "combine" : \
		       (num == GMRFLib_OPENMP_PLACES_EXTERNAL ? "external" : \
			(num == GMRFLib_OPENMP_PLACES_TIMING ? "timing" : \
			 (num == GMRFLib_OPENMP_PLACES_GCPO_BUILD ? "gcpo.build" : \
			  (num == GMRFLib_OPENMP_PLACES_SPECIAL ? "special" : \
			   (num == GMRFLib_OPENMP_PLACES_DEFAULT ? "default" : \
			    (num == GMRFLib_OPENMP_PLACES_NONE ? "none" : \
			     "This should not happen"))))))))))))))

typedef struct {
	GMRFLib_openmp_place_tp place;
	GMRFLib_openmp_strategy_tp strategy;
	int max_threads;
	int max_threads2;
	int *max_threads_nested;
	int blas_num_threads_force;
	// for PARDISO, like _outer is the number of threads in the outer loop, while _inner is the number of threads for
	// pardiso. the _inner is only relevant if nested=1.
	int max_threads_outer;
	int max_threads_inner;
	// when this is TRUE, then do PARDISO is parallel if the function call is serial
	int adaptive;
	// default schedule
	omp_sched_t schedule;
	int chunk_size;
	// optimal number of threads for likelihood computations
	int likelihood_nt;
} GMRFLib_openmp_tp;

#define GMRFLib_MAX_THREADS() (GMRFLib_openmp->max_threads)
#define GMRFLib_MAX_THREADS2() (GMRFLib_openmp->max_threads2)

#define GMRFLib_ADAPTIVE_NUM_THREADS() (GMRFLib_openmp->adaptive ? GMRFLib_openmp->adaptive : GMRFLib_openmp->max_threads_nested[1])

#define GMRFLib_OPENMP_IN_SERIAL()                  ((omp_get_num_threads() == 1) && (omp_get_level() == 0))
#define GMRFLib_OPENMP_IN_PARALLEL()                (!GMRFLib_OPENMP_IN_SERIAL())
#define GMRFLib_OPENMP_IN_PARALLEL_ONE_THREAD()     ((omp_get_num_threads() == 1) && (omp_get_level() == 1))
#define GMRFLib_OPENMP_IN_PARALLEL_ONEPLUS_THREAD() (omp_in_parallel() == 1)

#define GMRFLib_OPENMP_IN_OUTER() (omp_get_level() == 0)
#define GMRFLib_OPENMP_IN_INNER() (omp_get_level() == 1)

#define GMRFLib_OPENMP_NUM_THREADS_LEVEL() (GMRFLib_OPENMP_IN_OUTER() ? GMRFLib_openmp->max_threads_outer : \
					    GMRFLib_openmp->max_threads_inner)

#define GMRFLib_STOP_IF_NOT_SERIAL() assert(GMRFLib_OPENMP_IN_SERIAL() || GMRFLib_OPENMP_IN_PARALLEL_ONE_THREAD())

int GMRFLib_set_blas_num_threads(int threads);
int GMRFLib_openmp_nested_fix(void);
int GMRFLib_openmp_implement_strategy(GMRFLib_openmp_place_tp place, void *arg, GMRFLib_smtp_tp * smtp);
int GMRFLib_openmp_implement_strategy_special(int outer, int inner);

#if defined(INLA_WITH_MKL)
void MKL_Set_Num_Threads(int);
#endif

void GMRFLib_openmp_chunk(int n, double *A, double *b);
void GMRFLib_openmp_timing(void);

void GMRFLib_adapt_nt_init(int max_levels);
int GMRFLib_adapt_nt_get(char *tag, int thread_num, int level, int default_num_threads);
void GMRFLib_adapt_nt_update(char *tag, int thread_num, int level, double wtime);
void GMRFLib_adapt_nt_print(FILE * fp);

__END_DECLS
#endif
