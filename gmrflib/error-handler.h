
/* error-handler.h
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

/**
  \file error-handler.h
  \brief Typedefs and defines for \ref error-handler.c
*/

#ifndef __GMRFLib_ERROR_HANDLER_H__
#define __GMRFLib_ERROR_HANDLER_H__

#include <strings.h>
#include <stdio.h>
#if !defined(__FreeBSD__)
#include <malloc.h>
#endif
#include <stdlib.h>

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
  These are documented in error-handler.c
*/
#define  GMRFLib_SUCCESS     (0)			       /* OK !!!MUST BE ZERO!!!! */
#define  GMRFLib_EMEMORY     (1)			       /* alloc failed */
#define  GMRFLib_EPOSDEF     (2)			       /* matrix not positive definite */
#define  GMRFLib_ESINGMAT    (3)			       /* singular matrix */
#define  GMRFLib_EINVARG     (4)			       /* invalid argument */
#define  GMRFLib_EGRAPH      (5)			       /* graph is invalid */
#define  GMRFLib_EREORDER    (6)			       /* error reordering the graph */
#define  GMRFLib_EREADFILE   (7)			       /* error reading file */
#define  GMRFLib_EOPENFILE   (8)			       /* error opening file, file not found or readable */
#define  GMRFLib_EPARAMETER  (9)			       /* invalid parameter */
#define  GMRFLib_EINDEX      (10)			       /* no index found */
#define  GMRFLib_EPTR        (11)			       /* pointer is required to be non-NULL */
#define  GMRFLib_EOPTNR      (12)			       /* The Newton-Raphson optimizer did not converge */
#define  GMRFLib_EOPTCG      (13)			       /* The Conjugate-Gradient optimizer did not converge */
#define  GMRFLib_EOPTCGLINE  (14)			       /* The Conjugate-Gradient line-optimizer did not converge */
#define  GMRFLib_EMAPKIT     (15)			       /* Error occured in the mapkit-library (hash) */
#define  GMRFLib_ESMTP       (16)			       /* Sparse-matrix type (or function) not implemented */
#define  GMRFLib_ESINGCONSTR (17)			       /* Constraints (A) or its covariance matrix is singular */
#define  GMRFLib_ESN         (18)			       /* Fitting of the Skew-Normal distribution fail to converge */
#define  GMRFLib_EGEOCOOF    (19)			       /* Geo-coefficients are not available */
#define  GMRFLib_EGSL        (20)			       /* Error occured in the GSL-Library */
#define  GMRFLib_ESNH        (21)			       /* This should not happen */
#define  GMRFLib_EWRITE      (22)			       /* Error writing to file */
#define  GMRFLib_EMISC       (23)			       /* Misc error */
#define  GMRFLib_EPARDISO_LICENSE_NOTFOUND (24)		       /*  */
#define  GMRFLib_EPARDISO_LICENSE_EXPIRED  (25)		       /*  */
#define  GMRFLib_EPARDISO_LICENSE_ERR_USERNAME  (26)	       /*  */
#define  GMRFLib_EPARDISO_INTERNAL_ERROR  (27)		       /*  */
#define  GMRFLib_EPARDISO_NO_LIBRARY  (28)		       /*  */


/**
  \brief A template function declaration for specifying error-handling functions.

  \param reason A character string describing the reason for the error 
  (see \c GMRFLib_error_reason()).
  \param function The name of the function in which the error occurs.
  \param file The file name.
  \param line The number of the file.
  \param id The RCS Id (version control Id) of the file.
  \param errorno The error id number (see <tt> \ref GMRFLib_error_reason()) </tt>.
  \param msg Error message. 

  \return The function is to return the error number \a errorno.
  \sa GMRFLib_error_handler, GMRFLib_set_error_handler.
 */
typedef int GMRFLib_error_handler_tp(const char *reason, const char *file, const char *function, int line, const char *id, int errorno, const char *msg);

const char *GMRFLib_error_reason(int errorno);

int GMRFLib_error_handler(const char *reason, const char *file, const char *function, int line, const char *id, int errorno, const char *msg);
int GMRFLib_handle_error(const char *file, const char *function, int line, const char *id, int errorno, const char *msg);
int GMRFLib_set_error_handler(GMRFLib_error_handler_tp * new_error_handler);
int GMRFLib_error_handler_null(const char *reason, const char *file, const char *function, int line, const char *id, int errorno, const char *msg);
GMRFLib_error_handler_tp *GMRFLib_set_error_handler_off(void);

/* 
   long general versions
*/
#define GMRFLib_ERROR_MSG(errorno,msg) \
       if (1) { GMRFLib_handle_error(__FILE__, __GMRFLib_FuncName, __LINE__, \
                                     (const char *)RCSId, errorno,msg); return errorno; }
#define GMRFLib_ERROR_MSG_NO_RETURN(errorno,msg) \
       if (1) { GMRFLib_handle_error(__FILE__, __GMRFLib_FuncName, __LINE__, \
                                     (const char *)RCSId, errorno,msg); }

/* 
   short versions, no `msg'
*/
#define GMRFLib_ERROR(errorno) GMRFLib_ERROR_MSG(errorno,NULL)
#define GMRFLib_ERROR_NO_RETURN(errorno) GMRFLib_ERROR_MSG_NO_RETURN(errorno,NULL)

/* 
   easy use macros for everyday use!
*/
#define GMRFLib_ASSERT_RETVAL(condition,errorno,retval) \
            if (!(condition)) { \
                GMRFLib_handle_error(__FILE__, __GMRFLib_FuncName, __LINE__,\
		 		     (const char *)RCSId, errorno,  "Condition `" __GMRFLib_STRINGIFY(condition) "' is not TRUE"); \
                return retval;}
#define GMRFLib_ASSERT(condition,errorno) GMRFLib_ASSERT_RETVAL(condition,errorno,errorno)

/* 
   macro for writing a message
*/
#define GMRFLib_msg(fp,msg) if (1) {fprintf(fp,"\n%s: %s: %1d: %s\n",__FILE__,__GMRFLib_FuncName, __LINE__,msg);}

__END_DECLS
#endif
