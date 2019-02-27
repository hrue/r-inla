
/* optimize.h
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
  \file optimize.h
  \brief Typedefs and defines for \ref optimize.c 
*/

#ifndef __GMRFLib_OPTIMIZE_H__
#define __GMRFLib_OPTIMIZE_H__

#include <math.h>
#include <strings.h>
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

/*!
  \brief Use a conjugate gradiant method for optimisation
*/
#define GMRFLib_OPTTYPE_CG     0

/*!
  \brief Use a Newton-Raphson method for optimisation
*/
#define GMRFLib_OPTTYPE_NR     1

/*!
  \brief Use a  two-step  conjugate gradiant method for optimisation
*/
#define GMRFLib_OPTTYPE_SAFECG 2

/*!
  \brief Use a two-step Newton-Raphson method for optimisation
*/
#define GMRFLib_OPTTYPE_SAFENR 3

/*!

  \brief A template for routines computing the functions \f$ f(x_i) \f$ (and optionally the
  derivaties) of <b>(GMRF-30)</b> in \ref block.

  In most applications, these functions are equal to the elements of the log-likelihood of the GMRF
  <em>\b x</em>, given the data.

  The common usage is to let these functions are equal to the elements of the log-likelihood of the GMRF
  <em>\b x</em>, given the data, in which case the arguments are as follows.
  
  \param[out] logll At output, \a logll is a length \a m array holding 
  the values of \f$ f(x_i) \f$, computed for \a m different values of 
  element \a idx of <em>\b x</em>.
  \param[in] x_i A length \a m array of values of 
  \f$ x_{\mbox{\small\tt idx}} \f$, for which to compute
  \f$ f(x_{\mbox{\small\tt idx}}) \f$.
  \param[in] m The number of values of \f$ x_{\mbox{\small\tt idx}} \f$
  for which to compute the function values (either 0, 1 or 3).
  \param[in] idx The index of <em>\b x</em> corresponding to the values
  in \a x_i.
  \param[in] x_vec A length \em n array holding the current values 
  of all elements of the GMRF <em>\b x</em>. The value of
  <tt>x_vec[idx]</tt> is substituted by each element of
  \a x_i in turn.
  \param[in] logl_arg Additional arguments to the function. The
  \a logl_arg \em void -pointer is to hold the address of a variable 
  or data structure holding the additional arguments to the function, 
  including the data.

  \return The function is to return an error code (or 0).

  \par Example:
  The function \c loglik() computes the Poisson log-likelihood 
  \f$ y_{\mbox{\small\tt idx}} \sim \mbox{Pois}(E_{\mbox{\small\tt idx}} 
  \exp(x_{\mbox{\small\tt idx}})) \f$.
  \verbinclude doxygen_optimize_1.txt

  It is \b optional to let these functions also compute the exact first and second order derivaties
  (wrt \b x_i), in which case the following procedure must be adopted.

  - When this function is called with \f$m=0\f$, then the function \b MUST return \c GMRFLib_LOGL_COMPUTE_DERIVATIES, which
    indicate that for for this \b idx, the exact derivatie and second order derivaties are available. When this function is
    called with \f$m=3\f$, then it must return the function value in \f$logl[0]\f$, the derivative in \f$logl[1]\f$ and the
    second order derivative in \f$logl[2]\f$, all those evaluated at \b x_i[0], \b x_i[1] and \b x_i[2], respectively. (However,
    <em> x_i[0] = x_i[1] = x_i[2]</em> in this case.) If called with \f$m=2\f$, the it must return the function value in \f$logl[0]\f$
    and the derivative in \f$logl[1]\f$.
    OPS: if called with \f$m>3\f$ or $\f$m=1\f$, then it must return the log-likelihood for all <em>x_i[j]</em> for
    \f$j=0,\ldots, m-1\f$.

    \note Example
    \verbatim 
    int loglik_poisson(double *logll, double *x, int m, int idx, double *x_vec, char *arg)
    {
        //implement the log-likelihood for a poisson, y_idx|... ~ Po(\exp(x_idx)) with exact derivaties
        int i;
        double *y;
        y = (double *)arg;				 // for example

        if (m == 0) return GMRFLib_LOGL_COMPUTE_DERIVATIES;	// tell that exact derivatives are provided

	if (m > 3) {                                     // evaluate log-likelihoods
	    for(i=0;i<m;i++) {
	         logll[i] = y[idx]*x[i]-exp(x[i]);
	    }
	} else {
	                                                 // special section
            logll[0] = y[idx]*x[0]-exp(x[0]);            // f
            if (m > 1) {
	        logll[1] = y[idx] - exp(x[1]);		 // df
	    }
	    if (m > 2) {
	        logll[2] =        - exp(x[2]);		 // ddf 
	    }
	}
        return GMRFLib_SUCCESS;
    }
    \endverbatim

  \sa GMRFLib_blockupdate
*/
typedef int GMRFLib_logl_tp(double *logll, double *x_i, int m, int idx, double *x_vec, double *y_cdf, void *logl_arg);

/*!
  \brief The value to return to tell that exact derivaties are provided.

  \sa GMRFLib_logl_tp
*/
#define GMRFLib_LOGL_COMPUTE_DERIVATIES (135792467)
#define GMRFLib_LOGL_COMPUTE_DERIVATIES_AND_CDF (135792468)
#define GMRFLib_LOGL_COMPUTE_CDF (135792469)

/*!
  \struct GMRFLib_optimize_param_tp optimize.h
  \brief To specify the options for the optimizer.

  If no object of type \c GMRFLib_optimize_param_tp is specified and supplied to the block-updating-routine \c
  GMRFLib_blockupdate(), default values are used (see \c GMRFLib_default_optimize_param()). \n To change the default settings, a
  \c GMRFLib_optimize_param_tp -object should be dynamically allocated and the members should be set to the preferred values.  A
  pointer to this object should then be supplied as an argument to \c GMRFLib_blockupdate().

  \note If your log-likelihood is f.ex picewise linear, then the second derivative will be zero, hence the Taylor expantion will
  be a bad approximation to the log-likelihood curve! You can avoid (partially) this problem if you increase \a step_len.
 */
typedef struct {

	/**
	 *  \brief If \c != \c NULL, then printing output from the optimizer to \c fp 
	 */
	FILE *fp;

	/**
	 *  \brief The optmizing method to be used
	 * 
	 * There are four optimisation methods are available: - \c #GMRFLib_OPTTYPE_CG The conjugate gradient method - \c
	 * #GMRFLib_OPTTYPE_NR The Newton-Raphson method - \c #GMRFLib_OPTTYPE_SAFECG A two-step conjugate gradient method - \c 
	 * #GMRFLib_OPTTYPE_SAFENR A two-step Newton-Raphson method 
	 */
	int opt_type;

	/**
	 *  \brief Number of search directions
	 * 
	 * Indicates the number of previously search gradient directions on which the current search direction should be
	 * orthogonal (using the conjugate gradient (CG) method). 
	 */
	int nsearch_dir;

	/**
	 *  \brief Restart interval
	 * 
	 * If <em>restart_interval = r </em>, the CG search will be restarted every <em>r</em>'th iteration. 
	 */
	int restart_interval;

	/**
	 *  \brief Maximum number of iterations in each linesearch in the conjugate gradient method 
	 */
	int max_linesearch_iter;

	/**
	 *  \brief The abs(error) stopping criteria
	 * 
	 * The absolute error tolerance, relative to the number of nodes, for the size of one step of the optimizer. \n\n 
	 */
	double abserr_step;

	/**
	 *  \brief Maximum number of iterations 
	 */
	int max_iter;

	/**
	 *  \brief Run for a fixed number of iterations
	 * 
	 * If \c fixed_iter > 0, then this fix the number of iterations, whatever all other stopping-options. 
	 */
	int fixed_iter;

	/**
	 *  \brief abs(error) in function
	 * 
	 * The absolute error tolerance for the value of the function to be optimized. 
	 */
	double abserr_func;

	/**
	 *  \brief The step-length \f$h\f$ in discrete Taylor-approximation
	 * 
	 * Step length in the (discrete) computation of a Taylor expansion or second order approximation of the log-likelihood
	 * around a point <em> \b x_0</em> (CG and NR). \n\n 
	 */
	double step_len;
	int stencil;

	/**
	 * \brief Use reduced step-len in the Newton-Raphson routine
	 *
	 * Use reduced step-len in the Newton-Raphson routine, where the step-len for iteration \c iter, is minimum of \c (iter+1)
	 * \c times \c nr_step_factor and \c 1.
	 */
	double nr_step_factor;
} GMRFLib_optimize_param_tp;

typedef struct {
	int *map;

	double *mode;
	double *b;
	double *d;
	double *x_vec;
	GMRFLib_graph_tp *sub_graph;			       /* the subgraph */
	GMRFLib_Qfunc_tp *sub_Qfunc;			       /* the Qfunc */
	GMRFLib_Qfunc_arg_tp *sub_Qfunc_arg;		       /* ant its arguments */

	GMRFLib_constr_tp *sub_constr;			       /* stochastic constraint */

	GMRFLib_logl_tp *loglFunc;
	void *loglFunc_arg;

	GMRFLib_optimize_param_tp *optpar;
} GMRFLib_optimize_problem_tp;

int GMRFLib_default_optimize_param(GMRFLib_optimize_param_tp ** optpar);
int GMRFLib_optimize_set_store_flags(GMRFLib_store_tp * store);
int GMRFLib_optimize(double *mode, double *b, double *c, double *mean,
		     GMRFLib_graph_tp * graph, GMRFLib_Qfunc_tp * Qfunc, void *Qfunc_args,
		     char *fixed_value, GMRFLib_constr_tp * constr, double *d, GMRFLib_logl_tp * loglFunc, void *loglFunc_arg, GMRFLib_optimize_param_tp * optpar);
int GMRFLib_optimize_store(double *mode, double *b, double *c, double *mean,
			   GMRFLib_graph_tp * graph, GMRFLib_Qfunc_tp * Qfunc, void *Qfunc_args,
			   char *fixed_value, GMRFLib_constr_tp * constr,
			   double *d, GMRFLib_logl_tp * loglFunc, void *loglFunc_arg, GMRFLib_optimize_param_tp * optpar, GMRFLib_store_tp * store);
int GMRFLib_optimize2(GMRFLib_optimize_problem_tp * opt_problem, GMRFLib_store_tp * store);
int GMRFLib_Qadjust(double *dir, double *odir, GMRFLib_graph_tp * graph, GMRFLib_Qfunc_tp * Qfunc, void *Qfunc_arg);

int GMRFLib_linesearch(GMRFLib_optimize_problem_tp * opt_problem, double *dir);
double GMRFLib_linesearch_func(double length, double *dir, GMRFLib_optimize_problem_tp * opt_problem);
int GMRFLib_optimize3(GMRFLib_optimize_problem_tp * opt_problem, GMRFLib_store_tp * store);

__END_DECLS
#endif
