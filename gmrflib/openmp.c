#include <assert.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>

#include "GMRFLib/GMRFLib.h"

#define omp_max_max_nested() 3

// deprecated functions
#define omp_set_nested(v_) omp_set_max_active_levels(((v_) ? 3 : 1))
#define omp_get_nested() (omp_get_max_active_levels() > 1 ? 1 : 0)

static int blas_num_threads = 1;
int GMRFLib_get_blas_num_threads(void)
{
	return blas_num_threads;
}

int GMRFLib_set_blas_num_threads(int threads)
{
	if (threads < 1) {
		return GMRFLib_SUCCESS;
	}
	blas_num_threads = threads;
#if defined(INLA_WITH_MKL)
	void MKL_Set_Num_Threads(int);
	MKL_Set_Num_Threads(threads);
#endif
#if defined(INLA_WITH_OPENBLAS)
	void openblas_set_num_threads(int);
	openblas_set_num_threads(threads);
#endif
#if defined(INLA_WITH_ARMPL)
	void armpl_omp_set_num_threads(int);
	armpl_omp_set_num_threads(threads);
#endif
	return GMRFLib_SUCCESS;
}

int GMRFLib_openmp_implement_strategy_special(int outer, int inner)
{
	GMRFLib_openmp->place = GMRFLib_OPENMP_PLACES_SPECIAL;
	GMRFLib_openmp->max_threads_outer = IMAX(1, outer);
	GMRFLib_openmp->max_threads_inner = IMAX(1, inner);
	omp_set_nested((inner > 1));

	return GMRFLib_SUCCESS;
}

int GMRFLib_openmp_implement_strategy(GMRFLib_openmp_place_tp place, void *arg, GMRFLib_smtp_tp *smtp)
{
	GMRFLib_ENTER_FUNCTION;

	int debug = GMRFLib_DEBUG_IF_TRUE();
	int nt;
	int ntmax = GMRFLib_MAX_THREADS();
	int strategy = GMRFLib_openmp->strategy;
	int *nhyper = (int *) arg;
	int nhyper_def = 5;

	if (nhyper == NULL) {
		nhyper = &nhyper_def;
	}
	// this check is done once only
	if (GMRFLib_pardiso_ok < 0) {
		GMRFLib_pardiso_ok = (GMRFLib_pardiso_check_install(0, 1) == GMRFLib_SUCCESS ? 1 : 0);
		if (debug) {
			printf("%s:%1d: pardiso-library installed and working? [%s]\n", __FILE__, __LINE__, (GMRFLib_pardiso_ok ? "YES" : "NO"));
		}
	}

	static GMRFLib_smtp_tp smtp_store = GMRFLib_SMTP_DEFAULT;
	if (smtp) {
		smtp_store = *smtp;
	}

	if (GMRFLib_pardiso_ok && (smtp_store == GMRFLib_SMTP_PARDISO || smtp_store == GMRFLib_SMTP_DEFAULT)) {
		strategy = GMRFLib_OPENMP_STRATEGY_PARDISO;
		if (debug) {
			printf("%s:%1d: Switch to strategy [%s]\n", __FILE__, __LINE__, GMRFLib_OPENMP_STRATEGY_NAME(strategy));
		}
	}
	// default values
	nt = ntmax;

	GMRFLib_openmp->place = place;
	switch (place) {
	case GMRFLib_OPENMP_PLACES_PARSE_MODEL:
	{
		// this is serial section, except for _scale_model computations which
		switch (strategy) {
		case GMRFLib_OPENMP_STRATEGY_SMALL:
		case GMRFLib_OPENMP_STRATEGY_MEDIUM:
		case GMRFLib_OPENMP_STRATEGY_LARGE:
		case GMRFLib_OPENMP_STRATEGY_DEFAULT:
		case GMRFLib_OPENMP_STRATEGY_HUGE:
		{
			GMRFLib_openmp->max_threads_outer = 1;
			GMRFLib_openmp->max_threads_inner = 1;
		}
			break;

		case GMRFLib_OPENMP_STRATEGY_PARDISO:
		{
			GMRFLib_openmp->max_threads_outer = 1;
			GMRFLib_openmp->max_threads_inner = GMRFLib_openmp->max_threads_nested[1];

		}
			break;

		case GMRFLib_OPENMP_STRATEGY_STILES:
		{
			GMRFLib_openmp->max_threads_outer = GMRFLib_openmp->max_threads_nested[0];
			GMRFLib_openmp->max_threads_inner = GMRFLib_openmp->max_threads_nested[1];
		}
			break;

		case GMRFLib_OPENMP_STRATEGY_NONE:
		default:
			assert(0 == 1);
		}
	}
		break;

	case GMRFLib_OPENMP_PLACES_GCPO_BUILD:
	{
		switch (strategy) {
		case GMRFLib_OPENMP_STRATEGY_SMALL:
		case GMRFLib_OPENMP_STRATEGY_MEDIUM:
		case GMRFLib_OPENMP_STRATEGY_LARGE:
		case GMRFLib_OPENMP_STRATEGY_DEFAULT:
		case GMRFLib_OPENMP_STRATEGY_HUGE:
		case GMRFLib_OPENMP_STRATEGY_PARDISO:
		{
			GMRFLib_openmp->max_threads_outer = nt;
			GMRFLib_openmp->max_threads_inner = 1;
		}
			break;

		case GMRFLib_OPENMP_STRATEGY_STILES:
		{
			GMRFLib_openmp->max_threads_outer = GMRFLib_openmp->max_threads_nested[0];
			GMRFLib_openmp->max_threads_inner = GMRFLib_openmp->max_threads_nested[1];
		}
			break;

		case GMRFLib_OPENMP_STRATEGY_NONE:
		default:
			assert(0 == 1);
		}
	}
		break;

	case GMRFLib_OPENMP_PLACES_BUILD_MODEL:
	{
		switch (strategy) {
		case GMRFLib_OPENMP_STRATEGY_SMALL:
		case GMRFLib_OPENMP_STRATEGY_MEDIUM:
		case GMRFLib_OPENMP_STRATEGY_LARGE:
		case GMRFLib_OPENMP_STRATEGY_DEFAULT:
		case GMRFLib_OPENMP_STRATEGY_HUGE:
		case GMRFLib_OPENMP_STRATEGY_PARDISO:
		{
			GMRFLib_openmp->max_threads_outer = nt;
			GMRFLib_openmp->max_threads_inner = 1;
		}
			break;

		case GMRFLib_OPENMP_STRATEGY_STILES:
		{
			GMRFLib_openmp->max_threads_outer = GMRFLib_openmp->max_threads_nested[0];
			GMRFLib_openmp->max_threads_inner = GMRFLib_openmp->max_threads_nested[1];
		}
			break;

		case GMRFLib_OPENMP_STRATEGY_NONE:
		default:
			assert(0 == 1);
		}
	}
		break;

	case GMRFLib_OPENMP_PLACES_BUILD_MODEL2:
	{
		switch (strategy) {
		case GMRFLib_OPENMP_STRATEGY_SMALL:
		case GMRFLib_OPENMP_STRATEGY_MEDIUM:
		case GMRFLib_OPENMP_STRATEGY_LARGE:
		case GMRFLib_OPENMP_STRATEGY_DEFAULT:
		case GMRFLib_OPENMP_STRATEGY_HUGE:
		case GMRFLib_OPENMP_STRATEGY_PARDISO:
		{
			GMRFLib_openmp->max_threads_outer = (nt > 1 ? 2 : 1);
			GMRFLib_openmp->max_threads_inner = IMAX(1, nt / 2L);
		}
			break;

		case GMRFLib_OPENMP_STRATEGY_STILES:
		{
			GMRFLib_openmp->max_threads_outer = GMRFLib_openmp->max_threads_nested[0];
			GMRFLib_openmp->max_threads_inner = GMRFLib_openmp->max_threads_nested[1];
		}
			break;

		case GMRFLib_OPENMP_STRATEGY_NONE:
		default:
			assert(0 == 1);
		}
	}
		break;

	case GMRFLib_OPENMP_PLACES_OPTIMIZE:
	{
		switch (strategy) {
		case GMRFLib_OPENMP_STRATEGY_SMALL:
		case GMRFLib_OPENMP_STRATEGY_MEDIUM:
		{
			nt = IMIN(*nhyper + 1, ntmax);
			GMRFLib_openmp->max_threads_outer = nt;
			GMRFLib_openmp->max_threads_inner = 1;
		}
			break;

		case GMRFLib_OPENMP_STRATEGY_LARGE:
		case GMRFLib_OPENMP_STRATEGY_DEFAULT:
		{
			nt = IMIN(2 * *nhyper + 1, ntmax);
			GMRFLib_openmp->max_threads_outer = nt;
			GMRFLib_openmp->max_threads_inner = 1;
		}
			break;

		case GMRFLib_OPENMP_STRATEGY_HUGE:
		case GMRFLib_OPENMP_STRATEGY_PARDISO:
		{
			GMRFLib_openmp->max_threads_outer = GMRFLib_openmp->max_threads_nested[0];
			GMRFLib_openmp->max_threads_inner = GMRFLib_openmp->max_threads_nested[1];
		}
			break;

		case GMRFLib_OPENMP_STRATEGY_STILES:
		{
			GMRFLib_openmp->max_threads_outer = GMRFLib_openmp->max_threads_nested[0];
			GMRFLib_openmp->max_threads_inner = GMRFLib_openmp->max_threads_nested[1];
		}
			break;

		case GMRFLib_OPENMP_STRATEGY_NONE:
		default:
			assert(0 == 1);
		}
	}
		break;

	case GMRFLib_OPENMP_PLACES_EXTERNAL:
	{
		switch (strategy) {
		case GMRFLib_OPENMP_STRATEGY_SMALL:
		case GMRFLib_OPENMP_STRATEGY_MEDIUM:
		case GMRFLib_OPENMP_STRATEGY_LARGE:
		case GMRFLib_OPENMP_STRATEGY_DEFAULT:
		{
			GMRFLib_openmp->max_threads_outer = nt;
			GMRFLib_openmp->max_threads_inner = 1;
		}
			break;

		case GMRFLib_OPENMP_STRATEGY_HUGE:
		case GMRFLib_OPENMP_STRATEGY_PARDISO:
		{
			GMRFLib_openmp->max_threads_outer = GMRFLib_openmp->max_threads_nested[0];
			GMRFLib_openmp->max_threads_inner = GMRFLib_openmp->max_threads_nested[1];
		}
			break;

		case GMRFLib_OPENMP_STRATEGY_STILES:
		{
			GMRFLib_openmp->max_threads_outer = GMRFLib_openmp->max_threads_nested[0];
			GMRFLib_openmp->max_threads_inner = GMRFLib_openmp->max_threads_nested[1];
		}
			break;

		case GMRFLib_OPENMP_STRATEGY_NONE:
		default:
			assert(0 == 1);
		}
	}
		break;

	case GMRFLib_OPENMP_PLACES_HESSIAN:
	{
		switch (strategy) {
		case GMRFLib_OPENMP_STRATEGY_SMALL:
		{
			nt = 1;
			GMRFLib_openmp->max_threads_outer = 1;
			GMRFLib_openmp->max_threads_inner = 1;
		}
			break;

		case GMRFLib_OPENMP_STRATEGY_MEDIUM:
		{
			GMRFLib_openmp->max_threads_outer = nt;
			GMRFLib_openmp->max_threads_inner = 1;
		}
			break;

		case GMRFLib_OPENMP_STRATEGY_LARGE:
		case GMRFLib_OPENMP_STRATEGY_DEFAULT:
		{
			GMRFLib_openmp->max_threads_outer = nt;
			GMRFLib_openmp->max_threads_inner = 1;
		}
			break;

		case GMRFLib_OPENMP_STRATEGY_HUGE:
		case GMRFLib_OPENMP_STRATEGY_PARDISO:
		{
			GMRFLib_openmp->max_threads_outer = GMRFLib_openmp->max_threads_nested[0];
			GMRFLib_openmp->max_threads_inner = GMRFLib_openmp->max_threads_nested[1];
		}
			break;

		case GMRFLib_OPENMP_STRATEGY_STILES:
		{
			GMRFLib_openmp->max_threads_outer = GMRFLib_openmp->max_threads_nested[0];
			GMRFLib_openmp->max_threads_inner = GMRFLib_openmp->max_threads_nested[1];
		}
			break;

		case GMRFLib_OPENMP_STRATEGY_NONE:
		default:
			assert(0 == 1);
		}
	}
		break;

	case GMRFLib_OPENMP_PLACES_HESSIAN_SCALE:
	case GMRFLib_OPENMP_PLACES_INTEGRATE_HYPERPAR:
	{
		switch (strategy) {
		case GMRFLib_OPENMP_STRATEGY_SMALL:
		{
			GMRFLib_openmp->max_threads_outer = nt;
			GMRFLib_openmp->max_threads_inner = 1;
		}
			break;

		case GMRFLib_OPENMP_STRATEGY_MEDIUM:
		{
			GMRFLib_openmp->max_threads_outer = nt;
			GMRFLib_openmp->max_threads_inner = 1;
		}
			break;

		case GMRFLib_OPENMP_STRATEGY_LARGE:
		case GMRFLib_OPENMP_STRATEGY_DEFAULT:
		case GMRFLib_OPENMP_STRATEGY_HUGE:
		case GMRFLib_OPENMP_STRATEGY_PARDISO:
		{
			GMRFLib_openmp->max_threads_outer = GMRFLib_openmp->max_threads_nested[0];
			GMRFLib_openmp->max_threads_inner = GMRFLib_openmp->max_threads_nested[1];
		}
			break;

		case GMRFLib_OPENMP_STRATEGY_STILES:
		{
			GMRFLib_openmp->max_threads_outer = GMRFLib_openmp->max_threads_nested[0];
			GMRFLib_openmp->max_threads_inner = GMRFLib_openmp->max_threads_nested[1];
		}
			break;

		case GMRFLib_OPENMP_STRATEGY_NONE:
		default:
			assert(0 == 1);
		}
	}
		break;

	case GMRFLib_OPENMP_PLACES_COMBINE:
	{
		switch (strategy) {
		case GMRFLib_OPENMP_STRATEGY_SMALL:
		case GMRFLib_OPENMP_STRATEGY_MEDIUM:
		case GMRFLib_OPENMP_STRATEGY_LARGE:
		case GMRFLib_OPENMP_STRATEGY_DEFAULT:
		case GMRFLib_OPENMP_STRATEGY_HUGE:
		case GMRFLib_OPENMP_STRATEGY_PARDISO:
		{
			GMRFLib_openmp->max_threads_outer = nt;
			GMRFLib_openmp->max_threads_inner = 1;
		}
			break;

		case GMRFLib_OPENMP_STRATEGY_STILES:
		{
			GMRFLib_openmp->max_threads_outer = nt;
			GMRFLib_openmp->max_threads_inner = 1;
		}
			break;

		case GMRFLib_OPENMP_STRATEGY_NONE:
		default:
			assert(0 == 1);
		}
	}
		break;

	case GMRFLib_OPENMP_PLACES_TIMING:
	{
		GMRFLib_openmp->max_threads_outer = GMRFLib_openmp->max_threads_nested[1];	/* YES! */
		GMRFLib_openmp->max_threads_inner = GMRFLib_openmp->max_threads_nested[1];
	}
		break;

	case GMRFLib_OPENMP_PLACES_DEFAULT:
	{
		switch (strategy) {
		case GMRFLib_OPENMP_STRATEGY_SMALL:
		{
			nt = 1;
			GMRFLib_openmp->max_threads_outer = nt;
			GMRFLib_openmp->max_threads_inner = 1;
		}
			break;

		case GMRFLib_OPENMP_STRATEGY_MEDIUM:
		case GMRFLib_OPENMP_STRATEGY_LARGE:
		case GMRFLib_OPENMP_STRATEGY_DEFAULT:
		{
			GMRFLib_openmp->max_threads_outer = nt;
			GMRFLib_openmp->max_threads_inner = 1;
		}
			break;

		case GMRFLib_OPENMP_STRATEGY_HUGE:
		case GMRFLib_OPENMP_STRATEGY_PARDISO:
		{
			GMRFLib_openmp->max_threads_outer = GMRFLib_openmp->max_threads_nested[0];
			GMRFLib_openmp->max_threads_inner = GMRFLib_openmp->max_threads_nested[1];
		}
			break;

		case GMRFLib_OPENMP_STRATEGY_STILES:
		{
			GMRFLib_openmp->max_threads_outer = GMRFLib_openmp->max_threads_nested[0];
			GMRFLib_openmp->max_threads_inner = GMRFLib_openmp->max_threads_nested[1];
		}
			break;

		case GMRFLib_OPENMP_STRATEGY_NONE:
		default:
			assert(0 == 1);
		}
	}
		break;

	case GMRFLib_OPENMP_PLACES_NONE:
	case GMRFLib_OPENMP_PLACES_SPECIAL:
	{
		switch (strategy) {
		case GMRFLib_OPENMP_STRATEGY_NONE:	       /* only one option allowed */
			break;
		case GMRFLib_OPENMP_STRATEGY_SMALL:
		case GMRFLib_OPENMP_STRATEGY_MEDIUM:
		case GMRFLib_OPENMP_STRATEGY_LARGE:
		case GMRFLib_OPENMP_STRATEGY_DEFAULT:
		case GMRFLib_OPENMP_STRATEGY_HUGE:
		case GMRFLib_OPENMP_STRATEGY_PARDISO:
		case GMRFLib_OPENMP_STRATEGY_STILES:
		default:
			assert(0 == 1);
		}
	}
	}


	if (0) {
		P(place);
		if (place > GMRFLib_OPENMP_PLACES_NONE)
			abort();
		printf("place = %s\n", GMRFLib_OPENMP_PLACE_NAME(place));
		P(GMRFLib_openmp->max_threads_outer);
		P(GMRFLib_openmp->max_threads_inner);
	}

	int nested = (GMRFLib_openmp->max_threads_outer > 1 && GMRFLib_openmp->max_threads_inner > 1);
	// P(nested);

	// only set if changed
	if ((nested && !omp_get_nested()) || (!nested && omp_get_nested())) {
		omp_set_nested(nested);
	}

	omp_sched_t kind;
	int chunk_size;
	omp_get_schedule(&kind, &chunk_size);
	if (kind != GMRFLib_openmp->schedule || chunk_size != GMRFLib_openmp->chunk_size) {
		omp_set_schedule(GMRFLib_openmp->schedule, GMRFLib_openmp->chunk_size);
	}

	omp_set_num_threads(GMRFLib_openmp->max_threads_outer);
	if (GMRFLib_openmp->blas_num_threads_force) {
		GMRFLib_set_blas_num_threads(GMRFLib_openmp->blas_num_threads_force);
	} else {
		GMRFLib_set_blas_num_threads(GMRFLib_openmp->max_threads_inner);
	}

	static char init_dynamic = 1;
	if (init_dynamic) {
		GMRFLib_openmp_dynamic_init(omp_max_max_nested());
		init_dynamic = 0;
	}

	if (debug) {
		printf("%s:%1d: smtp[%s] strategy[%s] place[%s] nested[%1d]\n", __FILE__, __LINE__,
		       GMRFLib_SMTP_NAME(smtp_store), GMRFLib_OPENMP_STRATEGY_NAME(strategy), GMRFLib_OPENMP_PLACE_NAME(place), omp_get_nested());
		printf("%s:%1d: max.threads[%1d] num.threads[%1d] blas.num.threads[%1d] max.inner[%1d] max.outer[%1d]\n", __FILE__, __LINE__,
		       GMRFLib_MAX_THREADS(), nt, blas_num_threads, GMRFLib_openmp->max_threads_inner, GMRFLib_openmp->max_threads_outer);
	}

	GMRFLib_LEAVE_FUNCTION;
	return GMRFLib_SUCCESS;
}

void GMRFLib_openmp_chunk(int n, double *A_vec, double *b_vec)
{
	double *AA = Malloc(ISQR(n), double);
	Memcpy(AA, A_vec, ISQR(n) * sizeof(double));
	gsl_matrix_view m = gsl_matrix_view_array(AA, n, n);
	gsl_vector_view b = gsl_vector_view_array(b_vec, n);
	gsl_vector *x = gsl_vector_alloc(n);
	int s;
	for (int i = 0; i < n; i++) {
		gsl_matrix_set(&m.matrix, i, i, n);
	}
	gsl_permutation *p = gsl_permutation_alloc(n);
	gsl_linalg_LU_decomp(&m.matrix, p, &s);
	gsl_linalg_LU_solve(&m.matrix, p, &b.vector, x);
	gsl_permutation_free(p);
	gsl_vector_free(x);
	Free(AA);
}

void GMRFLib_openmp_timing(void)
{
	int nmax = 18;
	int m = 50;
	// double s = 1.0E6;

	double *A = Malloc(ISQR(nmax), double);
	double *b = Malloc(nmax, double);
	for (int i = 0; i < ISQR(nmax); i++) {
		A[i] = GMRFLib_stdnormal();
	}
	for (int i = 0; i < nmax; i++) {
		b[i] = GMRFLib_stdnormal();
	}

	int nt_def = GMRFLib_MAX_THREADS();
	char *tag = Strdup("8117db4e4a6ae84f37bb33bd6760734bcf122e0b");
	while (1) {
		int nt = GMRFLib_openmp_dynamic_get_nt(tag, 0, 0, nt_def);
		double tref = -GMRFLib_timer();
#pragma omp parallel for num_threads(nt) schedule(static)
		for (int k = 0; k < m; k++) {
			GMRFLib_openmp_chunk(nmax, A, b);
		}
		tref += GMRFLib_timer();
		GMRFLib_openmp_dynamic_update(tag, 0, 0, tref);
	}
}

static map_strvp ***dyn_nt = NULL;
static int dyn_nt_max_levels = -1;

void GMRFLib_openmp_dynamic_init(int max_levels)
{
	dyn_nt_max_levels = max_levels;
	dyn_nt = Calloc(max_levels, map_strvp **);
	for (int i = 0; i < max_levels; i++) {
		dyn_nt[i] = Calloc(1 + GMRFLib_MAX_THREADS(), map_strvp *);
		for (int j = 0; j <= GMRFLib_MAX_THREADS(); j++) {
			dyn_nt[i][j] = Calloc(1, map_strvp);
			map_strvp_init(dyn_nt[i][j]);
		}
	}
}

int GMRFLib_openmp_dynamic_get_nt(char *tag, int thread_num, int level, int default_num_threads)
{
	void **p = map_strvp_ptr(dyn_nt[level][thread_num], tag);
	GMRFLib_openmp_dynamic_num_threads_tp *obj = NULL;
	if (!p) {
		obj = Calloc(1, GMRFLib_openmp_dynamic_num_threads_tp);
		obj->tag = Strdup(tag);
		obj->ntimes = Calloc(1 + GMRFLib_MAX_THREADS(), double);
		obj->min_num_try = 2;
		obj->done = 0;
		// code below assumes we need to start with default_num_threads
		obj->max_nt = obj->best_nt = obj->try_next_nt = default_num_threads;
		obj->acc_wtime = Calloc(1 + GMRFLib_MAX_THREADS(), double);
		map_strvp_set(dyn_nt[level][thread_num], obj->tag, (void *) obj);
	} else {
		obj = *((GMRFLib_openmp_dynamic_num_threads_tp **) p);
	}
	return obj->try_next_nt;
}

void GMRFLib_openmp_dynamic_update(char *tag, int thread_num, int level, double wtime)
{
	if (GMRFLib_opt_num_threads) {
		int debug = 0;
		void **p = map_strvp_ptr(dyn_nt[level][thread_num], tag);
		assert(p);
		GMRFLib_openmp_dynamic_num_threads_tp *obj = *((GMRFLib_openmp_dynamic_num_threads_tp **) p);
		obj->tot_times++;
		obj->ntimes[obj->try_next_nt]++;
		obj->acc_wtime[obj->try_next_nt] += wtime;

		if (!(obj->done) && obj->max_nt == 1) {
			obj->done = 1;
		}

		if (!(obj->done)) {
			double time_try_next = obj->acc_wtime[obj->try_next_nt] / obj->ntimes[obj->try_next_nt];
			double time_best = obj->acc_wtime[obj->best_nt] / obj->ntimes[obj->best_nt];

			if (obj->ntimes[obj->try_next_nt] < obj->min_num_try && time_try_next > time_best) {
				// no point of trying more, abort early
				obj->try_next_nt = obj->best_nt;
				obj->done = 1;
			} else if (obj->ntimes[obj->try_next_nt] >= obj->min_num_try) {
				int step = 2;

				if (obj->best_nt == obj->max_nt && obj->ntimes[obj->try_next_nt] == obj->min_num_try) {
					// when its all started, this happens, as we always start with _nt = max_nt
					obj->try_next_nt = IMAX(1, obj->best_nt - step);
				} else {
					if (time_try_next < time_best) {
						if (debug) {
							printf("\nFound new best (%.1g, %1d) < (%.1g, %1d)\n\n", time_try_next * 1E6, obj->try_next_nt,
							       time_best * 1E6, obj->best_nt);
						}
						int itmp = obj->try_next_nt;
						obj->try_next_nt = IMAX(1, IMIN(obj->max_nt, itmp - step));
						obj->best_nt = itmp;
					} else {
						// its over...
						obj->try_next_nt = obj->best_nt;
						obj->done = 1;
					}
				}
			}
			if (debug && !(obj->done)) {
				printf("[%s][%1d][%1d] UPDATE: next.nt %1d best.nt %1d time*1E6 %1g\n", obj->tag, level, thread_num,
				       obj->try_next_nt, obj->best_nt, 1.0E6 * obj->acc_wtime[obj->best_nt] / obj->ntimes[obj->best_nt]);
			}
		}
	}
}

void GMRFLib_openmp_dynamic_print(FILE *fp)
{
	if (GMRFLib_opt_num_threads) {
		fp = (fp ? fp : stdout);
		fprintf(fp, "\nDump of dyn_nt\n");
		double tot_save = 0.0;
		for (int i = 0; i < dyn_nt_max_levels; i++) {
			for (int j = 0; j <= GMRFLib_MAX_THREADS(); j++) {
				if (dyn_nt[i][j]) {
					map_strvp_storage *ptr = NULL;
					for (ptr = NULL; (ptr = map_strvp_nextptr(dyn_nt[i][j], ptr)) != NULL;) {
						GMRFLib_openmp_dynamic_num_threads_tp *r = ((GMRFLib_openmp_dynamic_num_threads_tp *) ptr->value);
						if (r) {
							double wtime = 0.0;
							int ntimes = 0.0;
							for (int k = 0; k <= r->max_nt; k++) {
								wtime += r->acc_wtime[k];
								ntimes += r->ntimes[k];
							}
							double tsave = (r->acc_wtime[r->max_nt] / r->ntimes[r->max_nt]) * ntimes - wtime;
							fprintf(fp, "\t[%s][lev=%1d][th=%1d] best=%1d max=%1d wtime=%.3fs n=%1d (save=%.3fs)\n",
								r->tag, i, j, r->best_nt, r->max_nt, wtime, ntimes, tsave);

							tot_save += tsave;
						}
					}
				}
			}
		}
		fprintf(fp, "\tTotal wtime saved = %.3fs, estimated wtime = %.3fs saved\n", tot_save, tot_save / GMRFLib_openmp->max_threads_nested[0]);
		fprintf(fp, "\n");
	}
}
