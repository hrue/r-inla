
/* blockupdate.h
 * 
 * Copyright (C) 2001-2023 Havard Rue
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
  \file blockupdate.h
  \brief Typedefs and defines for \ref blockupdate.c
*/

#ifndef __GMRFLib_BLOCKUPDATE_H__
#define __GMRFLib_BLOCKUPDATE_H__

#include <math.h>
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

/*!
  \brief Expand around the mode
*/
#define GMRFLib_MODEOPTION_MODE    0

/*!
  \brief Expand around the current configuration
*/
#define GMRFLib_MODEOPTION_CURRENT 1

/*!
  \struct GMRFLib_blockupdate_param_tp blockupdate.h
  \brief Holding information on how to set up the block-sampling algorithm and the optimizer used to
  compute the current mode of the posterior distribution. 

  Within the block-sampling algorithm, a quadratic approximation to the posterior distribution,
  <b>(GMRF-30)</b> in \ref block, is computed such that the sampling algorithms for GMRF's declared
  in \ref problem-setup.c can be used. The quadratic approximation is computed by a Taylor expansion
  around the mode of the posterior distribution, found by an internal optimizer, or around the
  current value of <em>\b x</em>. Which method to choose is defined by the value of the data
  structure \c GMRFLib_blockupdate_param_tp. \n Note, if there is no data present with
  <em>d_new</em> and/or <em>d_old = NULL</em>, then no optimization is done as it not needed for any
  purpose.  This data structure is also used to choose whether or not to display the calculated
  values of likelihoods and acceptance rates for the updates, and the name of the file on which to
  print the calculations. \n The simplest way to create a \c GMRFLib_blockupdate_param_tp -object
  holding the default values, is by calling the function \c GMRFLib_default_blockupdate_param().

  \note If your log-likelihood is e.g. picewise linear, then the second derivative will be zero,
  hence the Taylor expantion will be a bad approximation to the log-likelihood curve! You can avoid
  (partially) this problem if you increase \a step_len.
*/
    typedef struct {

	/**
	 *  \brief Configuration to Taylor-expand around.
	 * 
	 * If equal to #GMRFLib_MODEOPTION_MODE, the mode of the posterior distribution is located and the Taylor-expansion is
	 * done around the mode. If equal to #GMRFLib_MODEOPTION_CURRENT, the expansion is done around the current value of
	 * <em>\b x</em>. 
	 */
	int modeoption;

	/**
	 *  \brief Display intermediate calculations on \c fp, if \c fp != \c NULL 
	 */
	FILE *fp;

	/**
	 *  \brief Step-length to compute the Taylor-expantions
	 * 
	 * Step length in the computation of a Taylor expansion or second order approximation of the log-likelihood. 
	 */
	double step_len;
	int stencil;

} GMRFLib_blockupdate_param_tp;

int GMRFLib_default_blockupdate_param(GMRFLib_blockupdate_param_tp ** blockupdate_par);
int GMRFLib_2order_approx(int thread_id, double *a, double *b, double *c, double *dd, double d, double x0, int indx,
			  double *x_vec, GMRFLib_logl_tp * loglFunc, void *loglFunc_arg, double *step_len, int *stencil, double *cmin);
int GMRFLib_2order_taylor(int thread_id, double *a, double *b, double *c, double *dd, double d, double x0, int indx,
			  double *x_vec, GMRFLib_logl_tp * loglFunc, void *loglFunc_arg, double *step_len, int *stencil);
int GMRFLib_2order_approx_core(int thread_id, double *a, double *b, double *c, double *dd, double x0, int indx,
			       double *x_vec, GMRFLib_logl_tp * loglFunc, void *loglFunc_arg, double *step_len, int *stencil);

__END_DECLS
#endif
