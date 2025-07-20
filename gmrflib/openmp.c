#include <assert.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>

#include "GMRFLib/GMRFLib.h"

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

	if (debug) {
		printf("%s:%1d: smtp[%s] strategy[%s] place[%s] nested[%1d]\n", __FILE__, __LINE__,
		       GMRFLib_SMTP_NAME(smtp_store), GMRFLib_OPENMP_STRATEGY_NAME(strategy), GMRFLib_OPENMP_PLACE_NAME(place), omp_get_nested());
		printf("%s:%1d: max.threads[%1d] num.threads[%1d] blas.num.threads[%1d] max.inner[%1d] max.outer[%1d]\n", __FILE__, __LINE__,
		       GMRFLib_MAX_THREADS(), nt, blas_num_threads, GMRFLib_openmp->max_threads_inner, GMRFLib_openmp->max_threads_outer);
	}

	GMRFLib_LEAVE_FUNCTION;
	return GMRFLib_SUCCESS;
}
