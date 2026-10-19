#include "GMRFLib/GMRFLib.h"

#if defined(INLA_WITH_LIBR) && !defined(INLA_WITH_LIBR_DLOPEN)
#       include <R.h>
#       include <Rembedded.h>
#       include <Rinternals.h>
#       include <Rdefines.h>

// Rinterface.h is Unix-only; R.dll still exports this symbol
extern uintptr_t R_CStackLimit;
#endif

#include "GMRFLib/timer.h"
extern char *GMRFLib_tmpdir;

// two copies...
#define R_GENERIC_WRAPPER "inla.rgeneric.wrapper"
#define INLA_OK (0)

#undef __BEGIN_DECLS
#undef __END_DECLS
#ifdef __cplusplus
#       define __BEGIN_DECLS extern "C" {
#       define __END_DECLS }
#else
#       define __BEGIN_DECLS				       /* empty */
#       define __END_DECLS				       /* empty */
#endif

__BEGIN_DECLS
//
#define MATHLIB_STANDALONE
#if defined(ISNAN)
#       undef ISNAN
#endif
#include <Rmath.h>					       // define gammafn/bessel_k
char *Strdup(const char *s);
int my_file_exists(const char *filename);
int my_dir_exists(const char *filename);
int my_setenv(char *str, int prefix);
int GMRFLib_sprintf(char **ptr, const char *fmt, ...);

__END_DECLS
#include "R-interface.h"
#if defined(INLA_WITH_LIBR) && defined(INLA_WITH_LIBR_DLOPEN)

/*
 * Load libR at runtime instead of linking it. the binary then starts on machines with any R or with none, and rgeneric resolves the
 * R that actually runs it (through R_HOME) on first use. Only the small embedding API below is needed, so building this mode
 * requires no R at all.
 */
#       include <ltdl.h>
typedef void *SEXP;
typedef ptrdiff_t R_xlen_t;

#       define REALSXP 14				       /* SEXPTYPE values are fixed in R's ABI */
#       define STRSXP  16
#       ifndef TRUE
#              define TRUE  1
#              define FALSE 0
#       endif

static lt_dlhandle R_dlhandle = NULL;
static SEXP(*p_Rf_protect) (SEXP);
static void (*p_Rf_unprotect)(int);
static SEXP(*p_Rf_mkString) (const char *);
static SEXP(*p_Rf_mkChar) (const char *);
static SEXP(*p_Rf_install) (const char *);
static SEXP(*p_Rf_lang2) (SEXP, SEXP);
static SEXP(*p_Rf_lang3) (SEXP, SEXP, SEXP);
static SEXP(*p_Rf_lang4) (SEXP, SEXP, SEXP, SEXP);
static SEXP(*p_Rf_ScalarLogical) (int);
static SEXP(*p_Rf_allocVector) (unsigned int, R_xlen_t);
static SEXP(*p_R_tryEval) (SEXP, SEXP, int *);
static double *(*p_REAL)(SEXP);
static R_xlen_t(*p_XLENGTH) (SEXP);
static void (*p_SET_STRING_ELT)(SEXP, R_xlen_t, SEXP);
static int (*p_Rf_initEmbeddedR)(int, char **);
static void (*p_Rf_endEmbeddedR)(int);
static SEXP *p_R_GlobalEnv;
static uintptr_t *p_R_CStackLimit;

#       define PROTECT(x_)      p_Rf_protect(x_)
#       define UNPROTECT(n_)    p_Rf_unprotect(n_)
#       define mkString         p_Rf_mkString
#       define mkChar           p_Rf_mkChar
#       define install          p_Rf_install
#       define lang2            p_Rf_lang2
#       define lang3            p_Rf_lang3
#       define lang4            p_Rf_lang4
#       define ScalarLogical    p_Rf_ScalarLogical
#       define allocVector      p_Rf_allocVector
#       define R_tryEval        p_R_tryEval
#       define REAL             p_REAL
#       define XLENGTH          p_XLENGTH
#       define SET_STRING_ELT   p_SET_STRING_ELT
#       define Rf_initEmbeddedR p_Rf_initEmbeddedR
#       define Rf_endEmbeddedR  p_Rf_endEmbeddedR
#       define R_GlobalEnv      (*p_R_GlobalEnv)
#       define R_CStackLimit    (*p_R_CStackLimit)

static void *inla_R_dlsym_(const char *name, const char *alt)
{
	void *p = (void *) lt_dlsym(R_dlhandle, name);

	if (!p && alt) {
		p = (void *) lt_dlsym(R_dlhandle, alt);
	}
	if (!p) {
		LOG_ERROR("libR has no symbol [%s]\n", name);
		exit(1);
	}
	return p;
}

static void inla_R_dlopen_(void)
{
	if (R_dlhandle) {
		return;
	}
	if (lt_dlinit() != 0) {
		LOG_ERROR("lt_dlinit failed: %s\n", lt_dlerror());
		exit(1);
	}
	// R_HOME is validated (or guessed) by the caller before this runs.
	// Global loading: shared objects that R loads later (packages)
	// resolve their R symbols from the global namespace.
	lt_dladvise advise;

	lt_dladvise_init(&advise);
	lt_dladvise_global(&advise);
	// let libltdl append the platform's own extension: .so, .dylib, .dll
	lt_dladvise_ext(&advise);

	// where R keeps its library, per platform layout
	static const char *rel[] = { "lib/libR", "bin/x64/R", "bin/R", NULL };
	char *rhome = getenv((const char *) "R_HOME");

	for (int i = 0; rhome && rel[i] && !R_dlhandle; i++) {
		char *path = NULL;

		GMRFLib_sprintf(&path, "%s/%s", rhome, rel[i]);
		R_dlhandle = lt_dlopenadvise(path, advise);
	}
	if (!R_dlhandle) {
		R_dlhandle = lt_dlopenadvise("libR", advise);
	}
	lt_dladvise_destroy(&advise);
	if (!R_dlhandle) {
		LOG_ERROR("rgeneric needs R with a shared libR: %s\n", lt_dlerror());
		exit(1);
	}
	p_Rf_protect = (SEXP(*)(SEXP)) (intptr_t) inla_R_dlsym_("Rf_protect", NULL);
	p_Rf_unprotect = (void (*)(int)) (intptr_t) inla_R_dlsym_("Rf_unprotect", NULL);
	p_Rf_mkString = (SEXP(*)(const char *)) (intptr_t) inla_R_dlsym_("Rf_mkString", NULL);
	p_Rf_mkChar = (SEXP(*)(const char *)) (intptr_t) inla_R_dlsym_("Rf_mkChar", NULL);
	p_Rf_install = (SEXP(*)(const char *)) (intptr_t) inla_R_dlsym_("Rf_install", NULL);
	p_Rf_lang2 = (SEXP(*)(SEXP, SEXP)) (intptr_t) inla_R_dlsym_("Rf_lang2", NULL);
	p_Rf_lang3 = (SEXP(*)(SEXP, SEXP, SEXP)) (intptr_t) inla_R_dlsym_("Rf_lang3", NULL);
	p_Rf_lang4 = (SEXP(*)(SEXP, SEXP, SEXP, SEXP)) (intptr_t) inla_R_dlsym_("Rf_lang4", NULL);
	p_Rf_ScalarLogical = (SEXP(*)(int)) (intptr_t) inla_R_dlsym_("Rf_ScalarLogical", NULL);
	p_Rf_allocVector = (SEXP(*)(unsigned int, R_xlen_t)) (intptr_t) inla_R_dlsym_("Rf_allocVector", NULL);
	p_R_tryEval = (SEXP(*)(SEXP, SEXP, int *)) (intptr_t) inla_R_dlsym_("R_tryEval", NULL);
	p_REAL = (double *(*)(SEXP))(intptr_t) inla_R_dlsym_("REAL", NULL);
	p_XLENGTH = (R_xlen_t(*)(SEXP)) (intptr_t) inla_R_dlsym_("XLENGTH", "Rf_xlength");
	p_SET_STRING_ELT = (void (*)(SEXP, R_xlen_t, SEXP)) (intptr_t) inla_R_dlsym_("SET_STRING_ELT", NULL);
	p_Rf_initEmbeddedR = (int (*)(int, char **)) (intptr_t) inla_R_dlsym_("Rf_initEmbeddedR", NULL);
	p_Rf_endEmbeddedR = (void (*)(int)) (intptr_t) inla_R_dlsym_("Rf_endEmbeddedR", NULL);
	p_R_GlobalEnv = (SEXP *) (intptr_t) inla_R_dlsym_("R_GlobalEnv", NULL);
	p_R_CStackLimit = (uintptr_t *) (intptr_t) inla_R_dlsym_("R_CStackLimit", NULL);
}

/*
 * The external model packages compile against R's headers and so reference
 * the Rf_-prefixed names, which normally come from libR. Here libR is not
 * linked: forward the math ones to the standalone Rmath library (which
 * exports the unprefixed names), and make Rf_error fail loudly -- R's
 * error unwinding does not exist when these models run outside R anyway.
 */
#       include <stdarg.h>
//extern double gammafn(double);
//extern double bessel_k(double, double, double);

double Rf_gammafn(double x)
{
	return gammafn(x);
	// return exp(lgamma(x));
}

double Rf_bessel_k(double x, double alpha, double expo)
{
	return bessel_k(x, alpha, expo);
}

void Rf_error(const char *fmt, ...)
{
	va_list ap;

	va_start(ap, fmt);
	vfprintf(stderr, fmt, ap);
	va_end(ap);
	fprintf(stderr, "\n");
	exit(1);
}
#endif							       /* INLA_WITH_LIBR_DLOPEN */

static int POSSIBLY_UNUSED(R_init) = 1;
static int POSSIBLY_UNUSED(R_debug) = 0;
static char *R_home = NULL;
double R_rgeneric_cputime = 0.0;

#if !defined(_OPENMP)
#       error "OpenMP must be enabled."
#endif

void inla_set_R_home(char *home)
{
	R_home = (home ? Strdup(home) : NULL);
}

#if defined(INLA_WITH_LIBR)

int inla_R_do_(inla_R_cmd_tp cmd, void *a1, void *a2, void *a3, void *a4, void *a5, void *a6)
{
	if (R_init) {
		inla_R_init_();
	}

	int ret = 0;

#       pragma omp critical (Name_95227b3fc78ae25be9b977b6385cae68f179f781)
	{
		R_rgeneric_cputime -= GMRFLib_timer();

		switch (cmd) {
		case INLA_R_INIT:
			break;

		case INLA_R_ASSIGN:
		{
			ret = inla_R_assign_((const char *) a1, (int *) a2, (double *) a3);
		}
			break;

		case INLA_R_FUNCALL1:
		{
			ret = inla_R_funcall1_((int *) a1, (double **) a2, (const char *) a3, (int *) a4, (double *) a5);
		}
			break;

		case INLA_R_FUNCALL2:
		{
			ret = inla_R_funcall2_((int *) a1, (double **) a2, (const char *) a3, (const char *) a4, (int *) a5, (double *) a6);
		}
			break;

		case INLA_R_FUNCALL_JP:
		{
			ret = inla_R_funcall_jp_((int *) a1, (double **) a2, (const char *) a3, (int *) a4, (double *) a5, (void *) a6);
		}
			break;

		case INLA_R_GET:
		{
			ret = inla_R_get_((int *) a1, (double **) a2, (const char *) a3);
		}
			break;

		case INLA_R_INLALOAD:
		{
			ret = inla_R_inlaload_((const char *) a1);
		}
			break;

		case INLA_R_LIBRARY:
		{
			ret = inla_R_library_((const char *) a1);
		}
			break;

		case INLA_R_LOAD:
		{
			ret = inla_R_load_((const char *) a1);
		}
			break;

		case INLA_R_RGENERIC:
		{
			ret = inla_R_rgeneric_((int *) a1, (double **) a2, (const char *) a3, (const char *) a4, (int *) a5, (double *) a6);
		}
			break;

		case INLA_R_SOURCE:
		{
			ret = inla_R_source_((const char *) a1);
		}
			break;

		case INLA_R_EXIT:
		{
			ret = inla_R_exit_();
		}
			break;

		default:
			Assert(0 == 1);
		}

		R_rgeneric_cputime += GMRFLib_timer();
	}

	return ret;
}

int inla_R_exit_(void)
{
	if (!R_init) {
		if (R_debug) {
			LOG_INFO("exit");
		}

		Rf_endEmbeddedR(0);
		R_rgeneric_cputime = 0.0;
		R_init = 1;
		R_debug = 0;
		R_home = NULL;
	}

	return INLA_OK;
}

int inla_R_init_(void)
{
	if (R_init) {
#       pragma omp critical (Name_aac7e80b592e4a6319788827c44116e831460cd6)
		if (R_init) {
			R_debug = (getenv((const char *) "INLA_DEBUG_R") ? 1 : 0);

			// Check if R_HOME is set. If not, try to guess it, otherwise fail.
			char *rhome = (R_home ? R_home : getenv((const char *) "R_HOME"));

			if (!rhome || (rhome && (my_dir_exists(rhome) != INLA_OK))) {
				if (my_dir_exists("/Library/Frameworks/R.framework/Resources") == INLA_OK) {
					GMRFLib_sprintf(&rhome, "R_HOME=/Library/Frameworks/R.framework/Resources");
				} else if (my_dir_exists("/usr/lib64/R") == INLA_OK) {
					GMRFLib_sprintf(&rhome, "R_HOME=/usr/lib64/R");
				} else if (my_dir_exists("/usr/lib/R") == INLA_OK) {
					GMRFLib_sprintf(&rhome, "R_HOME=/usr/lib/R");
				} else if (my_dir_exists("/usr/local/lib64/R") == INLA_OK) {
					GMRFLib_sprintf(&rhome, "R_HOME=/usr/local/lib64/R");
				} else if (my_dir_exists("/usr/local/lib/R") == INLA_OK) {
					GMRFLib_sprintf(&rhome, "R_HOME=/usr/local/lib/R");
				} else {
					LOG_ERROR_WRAPPER(
						fprintf(stderr, "\n\n");
						fprintf(stderr, "*** ERROR: Environment variable R_HOME is not set or invalid.\n");
						fprintf(stderr, "*** ERROR: Evaluate this in R:  Sys.getenv(\"R_HOME\")\n");
						fprintf(stderr, "\n\n");
						);
					exit(1);
				}
				LOG_ERROR_WRAPPER(
					fprintf(stderr, "\n\n");
					fprintf(stderr, "*** WARNING: Environment variable R_HOME is not set or invalid.\n");
					fprintf(stderr, "*** WARNING: Set it to a _GUESSED_ value [%s]\n\n", rhome);
					);
				my_setenv(rhome, 0);
			} else {
				char *rrhome = NULL;

				GMRFLib_sprintf(&rrhome, "R_HOME=%s", rhome);
				my_setenv(rrhome, 0);
			}

#       if defined(INLA_WITH_LIBR_DLOPEN)
			// R_HOME is set by now; load libR and resolve the API.
			inla_R_dlopen_();
#       endif
			char *Rargv[4];

			Rargv[0] = Strdup("REmbeddedPostgres");
			Rargv[1] = Strdup("--gui=none");
			Rargv[2] = Strdup("--vanilla");
			Rargv[3] = Strdup("--quiet");
			Rf_initEmbeddedR((R_debug ? 3 : 4), Rargv);

			if (R_debug) {
				LOG_INFO("init");
			}

			// Disable C stack limit check
			R_CStackLimit = (uintptr_t) (-1);

			char *filename = NULL;

			GMRFLib_sprintf(&filename, "%s/inla_rgeneric_wrapper_XXXXXX", GMRFLib_tmpdir);
			int fd = mkstemp(filename);

			close(fd);
			FILE *fp = fopen(filename, "w");

			if (R_debug) {
				LOG_INFO("base::searchpaths()");
				LOG_INFO("utils::sessionInfo()");
			}
			fprintf(fp, "%s <- function(cmd, model, theta = NULL) INLA::%s(cmd, model, theta)\n", R_GENERIC_WRAPPER, R_GENERIC_WRAPPER);
			fclose(fp);
			inla_R_source_quiet_(filename);

			R_init = 0;
		}
	}

	return (INLA_OK);
}

int inla_R_library_(const char *library)
{
	if (!library) {
		return (INLA_OK);
	}

	if (R_debug) {
		LOG_INFO("load library [%s]", library);
	}

	SEXP e, result, yy;
	int error;

	yy = PROTECT(mkString(library));
	e = PROTECT(lang2(install("library"), yy));
	result = PROTECT(R_tryEval(e, R_GlobalEnv, &error));
	if (result == NULL || error) {
		LOG_ERROR("load library [%s] failed", library);
		exit(1);
	}
	UNPROTECT(3);

	return (INLA_OK);
}

int inla_R_source_(const char *filename)
{
	if (!filename)
		return (INLA_OK);

	if (R_debug) {
		LOG_INFO("source file [%s]", filename);
	}

	SEXP e, result, yy, ffalse, ttrue;
	int error;

	ffalse = PROTECT(ScalarLogical(FALSE));
	ttrue = PROTECT(ScalarLogical(TRUE));
	yy = PROTECT(mkString(filename));
	e = PROTECT(lang4(install("source"), yy, ffalse, ttrue));
	result = PROTECT(R_tryEval(e, R_GlobalEnv, &error));
	if (result == NULL || error) {
		LOG_ERROR("source R-file [%s] failed", filename);
		exit(1);
	}
	UNPROTECT(5);

	return (INLA_OK);
}

int inla_R_source_quiet_(const char *filename)
{
	if (!filename)
		return (INLA_OK);

	if (R_debug) {
		LOG_INFO("source file [%s]", filename);
	}

	SEXP e, result, yy, logical;
	int error;

	logical = PROTECT(ScalarLogical((R_debug ? TRUE : FALSE)));
	yy = PROTECT(mkString(filename));
	e = PROTECT(lang4(install("source"), yy, logical, logical));
	result = PROTECT(R_tryEval(e, R_GlobalEnv, &error));
	if (result == NULL || error) {
		LOG_ERROR("source R-file [%s] failed", filename);
		exit(1);
	}
	UNPROTECT(4);

	return (INLA_OK);
}

int inla_R_load_(const char *filename)
{
	if (!filename)
		return (INLA_OK);

	if (R_debug) {
		LOG_INFO("load file [%s]", filename);
	}

	SEXP e, result, yy;
	int error;

	yy = PROTECT(mkString(filename));
	e = PROTECT(lang2(install("load"), yy));
	result = PROTECT(R_tryEval(e, R_GlobalEnv, &error));
	if (result == NULL || error) {
		LOG_ERROR("load RData-file [%s] failed", filename);
		exit(1);
	}
	UNPROTECT(3);

	return (INLA_OK);
}

int inla_R_inlaload_(const char *filename)
{
	if (!filename)
		return (INLA_OK);

	if (R_debug) {
		LOG_INFO("load file [%s]", filename);
	}

	SEXP e, result, yy;
	int error;

	yy = PROTECT(mkString(filename));
	e = PROTECT(lang2(install("inla.load"), yy));
	result = PROTECT(R_tryEval(e, R_GlobalEnv, &error));
	if (result == NULL || error) {
		LOG_ERROR("inla.load file [%s] failed", filename);
		exit(1);
	}
	UNPROTECT(3);

	return (INLA_OK);
}

int inla_R_funcall2_(int *n_out, double **x_out, const char *function, const char *tag, int *n, double *x)
{
	/*
	 * Call function(tag,x), where x is a double vector of length n. output is 'x_out' with length 'n_out'
	 */

	if (R_debug) {
		LOG_INFO("function [%s] tag [%s] n [%1d]", function, tag, *n);
	}

	int error, i;
	SEXP xx, result, e;

	xx = PROTECT(allocVector(REALSXP, *n));
	for (i = 0; i < *n; i++) {
		REAL(xx)[i] = x[i];
	}
	if (tag) {
		SEXP yy = PROTECT(mkString(tag));

		e = PROTECT(lang3(install(function), yy, xx));
	} else {
		e = PROTECT(lang2(install(function), xx));
	}
	result = PROTECT(R_tryEval(e, R_GlobalEnv, &error));
	if (result == NULL || error) {
		LOG_ERROR("Calling R-function [%s] with tag [%s] and [%1d] arguments", function, tag, *n);
		exit(1);
	}
	*n_out = (int) XLENGTH(result);
	Assert(*n_out >= 0);
	if (*n_out > 0) {
		*x_out = (double *) calloc((size_t) (*n_out), sizeof(double));	/* otherwise I'' use the R-version... */
		for (i = 0; i < *n_out; i++) {
			(*x_out)[i] = REAL(result)[i];
		}
	} else {
		*x_out = NULL;
	}

	UNPROTECT((tag ? 4 : 3));

	return (INLA_OK);
}

int inla_R_funcall1_(int *n_out, double **x_out, const char *function, int *n, double *x)
{
	return inla_R_funcall2_(n_out, x_out, function, NULL, n, x);
}

int inla_R_funcall_jp_(int *n_out, double **x_out, const char *function, int *n, double *x, void *sexp)
{
	/*
	 * Call function(x, sexp), where x is a double vector of length n, and 'sexp' is some SEXP ptr.
	 * Output is 'x_out' with length 'n_out'
	 */

	if (R_debug) {
		LOG_INFO("function[%s] n[%1d]", function, *n);
	}

	int error, i;
	SEXP xx, result, e;

	xx = PROTECT(allocVector(REALSXP, *n));
	for (i = 0; i < *n; i++) {
		REAL(xx)[i] = x[i];
	}

	e = PROTECT(lang3(install(function), xx, (SEXP) sexp));
	result = PROTECT(R_tryEval(e, R_GlobalEnv, &error));
	if (result == NULL || error) {
		LOG_ERROR("Calling R-function [%s] with [%1d] arguments", function, *n);
		exit(1);
	}
	*n_out = (int) XLENGTH(result);
	Assert(*n_out >= 0);
	if (*n_out > 0) {
		*x_out = (double *) calloc((size_t) (*n_out), sizeof(double));	/* otherwise I'' use the R-version... */
		for (i = 0; i < *n_out; i++) {
			(*x_out)[i] = REAL(result)[i];
		}
	} else {
		*x_out = NULL;
	}

	UNPROTECT(3);

	return (INLA_OK);
}

int inla_R_assign_(const char *variable, int *n, double *x)
{
	/*
	 * variable = x
	 */

	if (R_debug) {
		LOG_INFO("assign [%s] n [%1d]", variable, *n);
	}

	int error, i;
	SEXP xx, result, e, yy;

	yy = PROTECT(mkString(variable));
	xx = PROTECT(allocVector(REALSXP, *n));
	for (i = 0; i < *n; i++) {
		REAL(xx)[i] = x[i];
	}
	e = PROTECT(lang3(install("assign"), yy, xx));
	result = PROTECT(R_tryEval(e, R_GlobalEnv, &error));
	if (result == NULL || error) {
		LOG_ERROR("assign [%s] with n [%1d] failed\n", variable, *n);
		exit(1);
	}
	UNPROTECT(4);

	return (INLA_OK);
}

int inla_R_get_(int *n_out, double **x_out, const char *variable)
{
	/*
	 * return variable
	 */

	if (R_debug) {
		LOG_INFO("get: [%s]\n", variable);
	}

	int error, i;
	SEXP result, e, yy;

	yy = PROTECT(mkString(variable));
	e = PROTECT(lang2(install("get"), yy));
	result = PROTECT(R_tryEval(e, R_GlobalEnv, &error));
	if (result == NULL || error) {
		LOG_ERROR("get [%s]", variable);
		exit(1);
	}
	*n_out = (int) XLENGTH(result);
	Assert(*n_out >= 0);
	if (*n_out > 0) {
		*x_out = (double *) calloc((size_t) (*n_out), sizeof(double));	/* otherwise I'' use the R-version... */
		for (i = 0; i < *n_out; i++) {
			(*x_out)[i] = REAL(result)[i];
		}
	} else {
		*x_out = NULL;
	}

	UNPROTECT(3);

	return (INLA_OK);
}

int inla_R_rgeneric_(int *n_out, double **x_out, const char *cmd, const char *model, int *n, double *theta)
{
	/*
	 * do the rgeneric call with CMD and THETA and given MODEL name for the model definition
	 */

	if (R_debug) {
		LOG_INFO("cmd [%s] model [%s]", cmd, model);
	}

	int error, i;
	SEXP xx_theta, result, e, yy, yyy;

	xx_theta = PROTECT(allocVector(REALSXP, *n));
	for (i = 0; i < *n; i++) {
		REAL(xx_theta)[i] = theta[i];
	}

	yy = PROTECT(mkString(cmd));
	yyy = PROTECT(mkString(model));
	e = PROTECT(lang4(install(R_GENERIC_WRAPPER), yy, yyy, xx_theta));
	result = PROTECT(R_tryEval(e, R_GlobalEnv, &error));
	if (result == NULL || error) {
		LOG_ERROR("rgeneric [%s] with model [%s] failed\n", cmd, model);
		UNPROTECT(5);
		exit(1);
	}
	*n_out = (int) XLENGTH(result);
	Assert(*n_out >= 0);
	if (*n_out > 0) {
		*x_out = (double *) calloc((size_t) (*n_out), sizeof(double));	/* otherwise I'' use the R-version... */
		for (i = 0; i < *n_out; i++) {
			(*x_out)[i] = REAL(result)[i];
		}
	} else {
		*x_out = NULL;
	}

	UNPROTECT(5);

	return (INLA_OK);
}

void *inla_R_vector_of_strings(int n, char **s)
{
	// make and return an SEXP object which is a vector of given strings

	if (n == 0) {
		return ((void *) NULL);
	}
	SEXP sexp = PROTECT(allocVector(STRSXP, n));
	for (int i = 0; i < n; i++) {
		SEXP str = PROTECT(mkChar(s[i]));
		SET_STRING_ELT(sexp, i, str);
	}
	return ((void *) sexp);
}

#else							       /* if defined(INLA_WITH_LIBR) */

void inla_R_no_lib(void)
{
	LOG_ERROR("libR is not supported in this build");
	exit(1);
}

int inla_R_do_(inla_R_cmd_tp POSSIBLY_UNUSED(cmd), void POSSIBLY_UNUSED(*a1), void POSSIBLY_UNUSED(*a2),
	       void POSSIBLY_UNUSED(*a3), void POSSIBLY_UNUSED(*a4), void POSSIBLY_UNUSED(*a5), void POSSIBLY_UNUSED(*a6))
{
	inla_R_no_lib();
	return 1;
}

int inla_R_exit_(void)
{
	inla_R_no_lib(); 
	return 1;
}

int inla_R_init_(void)
{
	inla_R_no_lib(); 
	return 1;
}

int inla_R_library_(const char POSSIBLY_UNUSED(*library))
{
	inla_R_no_lib(); 
	return 1;
}

int inla_R_source_(const char POSSIBLY_UNUSED(*filename))
{
	inla_R_no_lib(); 
	return 1;
}

int inla_R_load_(const char POSSIBLY_UNUSED(*filename))
{
	inla_R_no_lib(); 
	return 1;
}

int inla_R_inlaload_(const char POSSIBLY_UNUSED(*filename))
{
	inla_R_no_lib(); 
	return 1;
}

int inla_R_funcall2_(int POSSIBLY_UNUSED(*n_out), double POSSIBLY_UNUSED(**x_out), const char POSSIBLY_UNUSED(*function),
		     const char POSSIBLY_UNUSED(*tag), int POSSIBLY_UNUSED(*n), double POSSIBLY_UNUSED(*x))
{
	inla_R_no_lib(); 
	return 1;
}

int inla_R_funcall1_(int POSSIBLY_UNUSED(*n_out), double POSSIBLY_UNUSED(**x_out), const char POSSIBLY_UNUSED(*function),
		     int POSSIBLY_UNUSED(*n), double POSSIBLY_UNUSED(*x))
{
	inla_R_no_lib(); 
	return 1;
}

int inla_R_funcall_jp_(int POSSIBLY_UNUSED(*n_out), double POSSIBLY_UNUSED(**x_out), const char POSSIBLY_UNUSED(*function),
		       int POSSIBLY_UNUSED(*n), double POSSIBLY_UNUSED(*x), void POSSIBLY_UNUSED(*sexp))
{
	inla_R_no_lib(); 
	return 1;
}

int inla_R_assign_(const char POSSIBLY_UNUSED(*variable), int POSSIBLY_UNUSED(*n), double POSSIBLY_UNUSED(*x))
{
	inla_R_no_lib(); 
	return 1;
}

int inla_R_get_(int POSSIBLY_UNUSED(*n_out), double POSSIBLY_UNUSED(**x_out), const char POSSIBLY_UNUSED(*variable))
{
	inla_R_no_lib(); 
	return 1;
}

int inla_R_rgeneric_(int POSSIBLY_UNUSED(*n_out), double POSSIBLY_UNUSED(**x_out), const char POSSIBLY_UNUSED(*cmd),
		     const char POSSIBLY_UNUSED(*model), int POSSIBLY_UNUSED(*n), double POSSIBLY_UNUSED(*theta))
{
	inla_R_no_lib(); 
	return 1;
}

void *inla_R_vector_of_strings(int POSSIBLY_UNUSED(n), char POSSIBLY_UNUSED(**s))
{
	inla_R_no_lib(); 
	return (void *)NULL;
}

#endif
