
/* GMRFLib.h
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
  \file GMRFLib.h
  \brief Include all include-files needed for using the GMRFLib in the correct order.
*/

#ifndef __GMRFLib_H__
#define __GMRFLib_H__

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
#define GMRFLib_VERSION_MAJOR    "3"
#define GMRFLib_VERSION_MINOR    "0"
#define GMRFLib_VERSION_REVISION "0-snapshot"
#define GMRFLib_VERSION          "3.0-0-snapshot"
#if defined(WINDOWS)
#define GMRFLib_NEED_DRAND48  1				       /* include implementation of drand48() */
#define GMRFLib_NEED_SRAND48  1				       /* include implementation of srand48() */
#endif
void daxpby_(int *n, double *a, double *x, int *incx, double *b, double *y, int *incy);


#include <stdlib.h>

/* 
 *  include files we need from GSL
 */
#include <gsl/gsl_types.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_multifit_nlin.h>
#include <gsl/gsl_pow_int.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_sys.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_sf_psi.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_sort_vector.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_multifit.h>

/* 
 * include all the include-files in GMRFLib
 */
#include "GMRFLib/init.h"
#include "GMRFLib/error-handler.h"
#include "GMRFLib/utils.h"
#include "GMRFLib/idxval.h"
#include "GMRFLib/lapack-interface.h"
#include "GMRFLib/dot.h"
#include "GMRFLib/timer.h"
#include "GMRFLib/io.h"
#include "GMRFLib/taucs.h"
#include "GMRFLib/random.h"
#include "GMRFLib/graph.h"
#include "GMRFLib/seasonal.h"
#include "GMRFLib/rw.h"
#include "GMRFLib/tabulate-Qfunc.h"
#include "GMRFLib/smtp-pardiso.h"
#include "GMRFLib/sparse-interface.h"
#include "GMRFLib/problem-setup.h"
#include "GMRFLib/openmp.h"
#include "GMRFLib/fmesher-io.h"
#include "GMRFLib/interpol.h"
#include "GMRFLib/density.h"
#include "GMRFLib/globals.h"
#include "GMRFLib/hash.h"
#include "GMRFLib/optimize.h"
#include "GMRFLib/blockupdate.h"
#include "GMRFLib/distributions.h"
#include "GMRFLib/smtp-band.h"
#include "GMRFLib/smtp-taucs.h"
#include "GMRFLib/bitmap.h"				       /* needs both graph and problem and sparse */
#include "GMRFLib/ghq.h"
#include "GMRFLib/design.h"
#include "GMRFLib/approx-inference.h"
#include "GMRFLib/pre-opt.h"
#include "GMRFLib/graph-edit.h"
#include "GMRFLib/domin-interface.h"
#include "GMRFLib/integrator.h"
#include "GMRFLib/version.h"
#include "GMRFLib/hgmrfm.h"
#include "GMRFLib/matern.h"

#if defined(INLA_LINK_WITH_MKL)
void vdExp(int, const double *, double *);
void vdLog1p(int, const double *, double *);
void vsAdd(int n, const double *a, const double *b, double *y);
void vdLn(int, const double *, double *);
void vdLog1p(int, const double *, double *);
void vdLinearFrac(int n, double *a, double *b, double scalea, double shifta, double scaleb, double shiftb, double *y);
void vdSqr(int, double *, double *);
void vdMul(int, double *, double *, double *);
void vdAdd(int, double *, double *, double *);

void mkl_dcsrmv(const char *transa, const int *m, const int *k, const double *alpha,
		const char *matdescra, const double *val, const int *indx,
		const int *pntrb, const int *pntre, const double *x, const double *beta, double *y);
double cblas_ddoti(const int nz, const double *x, const int *indx, const double *y);
#endif



__END_DECLS
#endif
