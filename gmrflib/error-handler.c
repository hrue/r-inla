#include <stdarg.h>
#include <strings.h>
#include <stdio.h>
#include <stdlib.h>
#include <omp.h>

#include "GMRFLib/GMRFLib.h"

#ifndef GITCOMMIT
#define GITCOMMIT "devel"
#endif

static GMRFLib_error_handler_tp *handler = NULL;

#pragma omp threadprivate(handler)

const char *GMRFLib_error_reason(int errorno)
{
	/*
	 * return a pointer to the reason for error=errorno 
	 */
#define NMSG 33
	static const char *reasons[NMSG] = {
		"Warning",
		"Alloc failed",
		"Matrix is not (numerical) positive definite",
		"Matrix is (numerical) singular",
		"Invalid argument or argument combination",
		"Graph is invalid",
		"Error reordering the graph",
		"Error while reading file",
		"Error opening file, file not found or readable",
		"Invalid value of parameter",
		"No index found",
		"Pointer is required to be non-NULL",
		"The Newton-Raphson optimizer did not converge",
		"The Conjugate-Gradient optimizer did not converge",
		"The line-optimizer in the Conjugate-Gradient optimizer did not converge",
		"Error occured in the mapkit-library (hash)",
		"Sparse-matrix-type (or function) not implemented yet",
		"Covariance matrix for Ax is singular",
		"Fitting of the Skew-Normal distribution failed to converge",
		"Geo-coefficients are not available",
		"Error occured in the GSL-Library",
		"This should not happen",
		"Error writing file",
		"Misc error",
		"License file for PARDISO: Not found",
		"License file for PARDISO: Expired",
		"License file for PARDISO: Wrong username",
		"PARDISO: Internal error",
		"PARDISO: Library not available",
		"Singular constraints: Matrix AA' is singular (input error)",
		"dlopen error",
		"dlsym error",
		"(((this is an unknown errorcode)))"
	};
	if (errorno < 0 || errorno >= NMSG - 1)
		return (const char *) reasons[NMSG - 1];       /* return the last message in this case */
	else
		return (const char *) reasons[errorno];
#undef NMSG
}

int GMRFLib_error_handler(const char *reason, const char *UNUSED(file), const char *function, int line, int errorno, const char *msg)
{
	/*
	 * this is the default error-handler 
	 */
	fprintf(stderr, "\n\n\tGitId: %s\n", __GMRFLib_symbol_to_string(GITCOMMIT));
	if (reason) {
		fprintf(stderr, "\tError:%1d Reason: %s\n", errorno, reason);
	}
	if (msg) {
		fprintf(stderr, "\tMessage: %s\n", msg);
	}
	if (function) {
		fprintf(stderr, "\tLine:%1d Function: %s\n", line, function);
	}
	fprintf(stderr, "\n");

	if (errorno != GMRFLib_SUCCESS) {
		abort();				       /* no reason to abort if ok */
	}

	return GMRFLib_SUCCESS;
}

int GMRFLib_error_handler_null(const char *reason, const char *file, const char *function, int line, int errorno, const char *msg)
{
	/*
	 * gcc produce warnings here, so i have to fake some code to make them go away 
	 */

	int dummy = (reason || file || function || line == 0 || msg);

	return (dummy ? errorno : errorno);
}

GMRFLib_error_handler_tp *GMRFLib_set_error_handler_off(void)
{
	GMRFLib_error_handler_tp *h = handler;

	GMRFLib_set_error_handler(GMRFLib_error_handler_null);
	return h;
}

int GMRFLib_handle_error(const char *file, const char *function, int line, int errorno, const char *msg)
{
	/*
	 * this function handle an error, and is the one called from the library through the macro 'GMRFLib_ERROR' 
	 */
	if (handler) {
		(*handler) (GMRFLib_error_reason(errorno), file, function, line, errorno, msg);
	} else {
		GMRFLib_error_handler(GMRFLib_error_reason(errorno), file, function, line, errorno, msg);
	}
	return GMRFLib_SUCCESS;
}

int GMRFLib_set_error_handler(GMRFLib_error_handler_tp *new_error_handler)
{
	/*
	 * set new error-handler. if the new error-handler is NULL, then the default error-handler is used. 
	 */
	handler = new_error_handler;
	return GMRFLib_SUCCESS;
}
