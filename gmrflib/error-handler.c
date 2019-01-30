
/* error-handler.c
 * 
 * Copyright (C) 2001-2006 Havard Rue
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
  \file error-handler.c
  \brief Functions to handle error messages.

  If an error occurs running the routines from \em GMRFLib, an error message will be issued stating
  the type of error, and the function, the file and the line at which the error has occurred.  By
  default, the program is core dumped and aborted at the occurrence of an error, using the function
  \c GMRFLib_error_handler().  This side effect can be overridden by the user, re-defining the error
  handling function as described in \c GMRFLib_set_error_handler(), or simply turn-off the
  error-handling system using \c GMRFLib_set_error_handler_off().

  \note The behaviour when the error-handling system is turned off and the code is \a returnd \a in \a the \a call.  Temporary
   allocated memory may not be free'd if the error-handling system is turned off or altered such that the program continues. We
   hope to fix this in later releases.

  Example
  \verbinclude example-doxygen-error-handler.txt
  
*/

#include <stdarg.h>
#include <strings.h>
#include <stdio.h>
#if !defined(__FreeBSD__)
#include <malloc.h>
#endif
#include <stdlib.h>
#include "GMRFLib/GMRFLib.h"
#include "GMRFLib/GMRFLibP.h"

#ifndef HGVERSION
#define HGVERSION
#endif
static const char RCSId[] = "file: " __FILE__ "  " HGVERSION;

/* Pre-hg-Id: $Id: error-handler.c,v 1.49 2008/08/26 07:07:11 hrue Exp $ */
static GMRFLib_error_handler_tp *handler = NULL;

#pragma omp threadprivate(handler)

/*!
  \brief Returns the type of error, given the error code \a errorno.
  \param errorno The error id number (see below).
  \return A character string, giving the error message corresponding to error number
  \a errorno.

  \par Error numbers and error messages
  0 :    No error, please ignore \n
  1 :    Alloc failed \n
  2 :    Matrix is not positive definite \n
  3 :    Matrix is singular \n
  4 :    Invalid argument \n
  5 :    Graph is invalid \n
  6 :    Error reordering the graph \n
  7 :    Error while reading file \n
  8 :    Error opening file, file not found or readable \n
  9 :    Invalid value of parameter \n
  10 :   No index found \n
  11 :   Pointer is required to be non-NULL \n
  12 :   The Newton-Raphson optimizer did not converge  \n
  13 :   The Conjugate-Gradient optimizer did not converge \n
  14 :   The line-optimizer in the Conjugate-Gradient optimizer did not converge \n
  15 :   Error occured in the mapkit-library (hash)\n
  16 :   Sparse-matrix-type not implemented yet\n
  17 :   Constraints or its covariance matrix is singular\n
  18 :   Fitting of the Skew-Normal distribution fail to converge\n
  19 :   Geo-coefficients are not available\n
  20 :   Error occured in the GSL-Library\n
  21 :   This should not happen\n
  22 :   Error writing file
  23 :   Misc error\n
  24 :   PARDISO License file: not found\n
  25 :   PARDISO License file: expired\n
  26 :   PARDISO License file: wrong username\
  27 :   PARDISO: Internal error\
  28 :   PARDISO: Library not loaded or available\
  29 :   (this is an unknown error code) \n

  \remarks This function is used within the library generating the \a reason argument 
  in the default error handling function \c GMRFLib_error_handler().
  \sa GMRFLib_error_handler.
 */
const char *GMRFLib_error_reason(int errorno)
{
	/*
	 * return a pointer to the reason for error=errorno 
	 */
#define NMSG 30
	static const char *reasons[NMSG] = {
		"No error, please ignore",
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
		"Constraints or its covariance matrix is singular",
		"Fitting of the Skew-Normal distribution fail to converge",
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
		"(((this is an unknown errorcode)))"
	};
	if (errorno < 0 || errorno >= NMSG - 1)
		return (const char *) reasons[NMSG - 1];       /* return the last message in this case */
	else
		return (const char *) reasons[errorno];
#undef NMSG
}

/*!
  \brief The default error handling function, of type \c GMRFLib_error_handler_tp().

  \param reason A character string describing the reason for the error 
  (see \c GMRFLib_error_reason()).
  \param function The name of the function in which the error occurs.
  \param file The file name.
  \param line The number of the file.
  \param id The RCS Id (version control Id) of the file.
  \param errorno The error id number (see <tt> \ref GMRFLib_error_reason()) </tt>.
  \param msg Error message.

  \return The error number. The function is aborted if reporting an error, so a value
  is returned only if \a errorno = 0.

  \remarks This is the default error handling function integrated in the library routines.  It is
  called by the function \c GMRFLib_handle_error(). When an error occurs, the function prints an
  error message of the form \verbinclude doxygen_error_1.txt and the program is aborted,
  with a core dump. In the example above, the error message states that the command line arguments
  of the main program on the file \c test-graph.c is missing. To avoid aborting the program, a new
  error handling function can be defined, returning the error code without aborting. To turn off or
  change to another error handling function, use \c GMRFLib_set_error_handler_off() or \c
  GMRFLib_set_error_handler().

  \sa GMRFLib_error_handler_tp, GMRFLib_handle_error, GMRFLib_set_error_handler.
 */
int GMRFLib_error_handler(const char *reason, const char *file, const char *function, int line, const char *id, int errorno, const char *msg)
{
	/*
	 * this is the default error-handler 
	 */
	fprintf(stderr, "\n\n\tGMRFLib version %1s, has recived error no [%1d]\n", GMRFLib_VERSION, errorno);
	if (reason) {
		fprintf(stderr, "\tReason    : %s\n", reason);
	}
	if (msg) {
		fprintf(stderr, "\tMessage   : %s\n", msg);
	}
	if (function) {
		fprintf(stderr, "\tFunction  : %s\n", function);
	}
	if (file) {
		fprintf(stderr, "\tFile      : %s\n", file);
	}
	if (1) {
		fprintf(stderr, "\tLine      : %1d\n", line);
	}
	if (id) {
		fprintf(stderr, "\tRCSId     : %s\n", id);
	}
	fprintf(stderr, "\n");

	if (errorno != GMRFLib_SUCCESS) {
		abort();				       /* no reason to abort if ok */
	}

	return GMRFLib_SUCCESS;
}
int GMRFLib_error_handler_null(const char *reason, const char *file, const char *function, int line, const char *id, int errorno, const char *msg)
{
	/*
	 * gcc produce warnings here, so i have to fake some code to make them go away 
	 */

	int dummy = (reason || file || function || line == 0 || id || msg);

	return (dummy ? errorno : errorno);
}

/*!
  \brief Turn off the error handler and return the current one.

  This function turn off the error handling system and return the current error handler. This
  pointer can be passed to \c GMRFLib_set_error_handler() to restore the original beahviour.

  \sa GMRFLib_set_error_handler(), GMRFLib_error_handler()
*/
GMRFLib_error_handler_tp *GMRFLib_set_error_handler_off(void)
{
	GMRFLib_error_handler_tp *h = handler;

	GMRFLib_set_error_handler(GMRFLib_error_handler_null);
	return h;
}

/*!
  \brief Handles the error by calling an internal macro, issuing an error message.

  This is the function that is actually called at the occurrence of an error.
  \param file The file name.
  \param function The name of the function in which the error occurs.
  \param line The number of the file.
  \param id The RCS Id (version control Id) of the file.
  \param errorno The error id number (see <tt> \ref GMRFLib_error_reason()) </tt>.
  \param msg Error message.

  \sa GMRFLib_error_handler_tp, GMRFLib_error_handler, GMRFLib_set_error_handler.
 */
int GMRFLib_handle_error(const char *file, const char *function, int line, const char *id, int errorno, const char *msg)
{
	/*
	 * this function handle an error, and is the one called from the library through the macro 'GMRFLib_ERROR' 
	 */
	if (handler) {
		(*handler) (GMRFLib_error_reason(errorno), file, function, line, id, errorno, msg);
	} else {
		GMRFLib_error_handler(GMRFLib_error_reason(errorno), file, function, line, id, errorno, msg);
	}
	return GMRFLib_SUCCESS;
}

/*!
  \brief Sets the error handler function to be the default function or a user specified function.

  The default error messages can be overridden by the user by defining a new 
  error handling function of a specified format, \c GMRFLib_error_handler_tp(). 
  Within the library, a static function pointer \em handler is defined, pointing 
  to a function of type \c GMRFLib_error_handler_tp(). By default, this is \c NULL, 
  indicating that the default error handling function, \c GMRFLib_error_handler() 
  is to be used. To override the default function, \em handler should be set to
  point at the user-specified \c GMRFLib_error_handler_tp()-function, which can 
  be done using \c GMRFLib_set_error_handler().

  \param new_error_handler A pointer to a user-specified \c GMRFLib_error_handler_tp()
  -function. If \a new_error_handler = \c NULL, the default error handler function is chosen.

  \remarks To override the error handling function, the user should write a function with
  return value and argument list as defined by the \c GMRFLib_error_handler_tp() template,
  and then call \c GMRFLib_set_error_handler(). As a result, the user-specified function
  will be called at the occurrence of an error.

  \sa GMRFLib_error_handler_tp, GMRFLib_error_handler(), GMRFLib_set_error_handler_off()
 */
int GMRFLib_set_error_handler(GMRFLib_error_handler_tp * new_error_handler)
{
	/*
	 * set new error-handler. if the new error-handler is NULL, then the default error-handler is used. 
	 */
	handler = new_error_handler;
	return GMRFLib_SUCCESS;
}

/* 
   

void Error(char *fmt, ...) {
  va_list args;
  va_start(args, fmt);
  vfprintf(stderr, fmt, args);
  va_end(args);
  fprintf(stderr, "\nExiting ...\n");
  exit(1);
}

*/
