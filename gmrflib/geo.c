
/* geo.c
 * 
 * Copyright (C) 2005-2006 Havard Rue
 * 
 * This program is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation; either version 2 of the License, or (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
 * of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA
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
 */

/* This code was initially contributed by Hanne T. Wist, later rewritten by H. Rue */

/*!
  \file geo.c
  \brief Setup for the Gaussian field approximation problems.
 

  \note To use the routines in \ref geo.c, you need also to link with the "geo" module, with default name \c libGMRFLib-geo.a,
  see the \c Makefile in the \c $PREFIX/doc/GMRFLib/examples/ directory.  The \ref geo.c routines is collected in a separate
  module due to the huge tables of coefficients.

  Commonly used Gaussian fields on regular lattices can be well approximated using GMRF models with a small neighbourhood.\n In
  geostatistics commonly used isotropic covariance function (CF) are
  - Matern CF  \f[ C(h)=\frac{1}{\Gamma(\nu)2^{\nu-1}}(s_{\nu}h)^{\nu}K_{\nu}(s_{\nu}h)\f]
  - Exponential  \f[ C(h) = \exp\{-3h\} \f]
  - Gaussian  \f[ C(h) = \exp\left\{-3h^2\right\} \f]
  
  \f$ K_{\nu} \f$ is the modified Bessel function of the second kind and order \f$ \nu>0 \f$, and \f$ s_{\nu}\f$ is a function
  of \f$ \nu \f$ such that the covariance function is scaled to \f$ C(1)=0.05 \f$.\n The Matern CF also includes the exponential
  CF, with \f$ \nu=1/2 \f$, and the Gaussian CF, with \f$ \nu\rightarrow\infty \f$.\n A further parameter <em>r</em> (range) is
  often introduced to scale the Euclidean distance. \n\n

  Assume a Gaussian field that is isotropic with zero mean and restricted to the torus \f$\mathcal{T}_{n}\f$ (or to a
  lattice). Let <em>\b z</em> denote a zero mean Gaussian field on \f$\mathcal{T}_{n}\f$. Then \f[
  \mbox{Cov}(z_{ij},z_{i'j'})=C\left(\frac{d((i,j),(i',j'))}{r}\right) \f] where <em>d((i,j),(i',j'))</em> is the Euclidian
  distance between \f$ z_{ij}\f$ and \f$ z_{i'j'} \f$ on the torus.\n\n
  
  Let <em>\b x</em> be a GMRF defined on \f$\mathcal{T}_{n}\f$ with precision matrix \f$ \mbox{\boldmath $Q(\theta)$} \f$
  depending on some parameter vector.  Regarding the neighbourhood and parameterisation, we will use a square window of size \f$
  (2m+1)\times(2m+1) \f$ centered at each <em>(i,j)</em>, where <em> m=</em> 2 or 3. The coefficients are numbered as \f[
  \left[\begin{array}{ccc} & & \theta_6\\ & \theta_3 & \theta_5 \\ \theta_1 & \theta_2 & \theta_4 \end{array}\right] \mbox{ for
  }m=2\mbox{, and} \left[\begin{array}{cccc} & & & \theta_{10}\\ & & \theta_6 & \theta_9 \\ & \theta_3 & \theta_5 & \theta_8 \\
  \theta_1 & \theta_2 & \theta_4 & \theta_7 \end{array}\right] \mbox{ for }m=3 \f] representing the lower part of the upper
  right quadrant.\n\n

  The coefficients giving the best fit are found as \f[ \mbox{\boldmath $\theta^*$}=\arg \min_{\theta\in\Theta_{\infty}^+}
  \parallel \mbox{\boldmath $\rho-\rho(\theta)$}\parallel^2_w =\arg \min_{\theta\in\Theta_{\infty}^+}
  \sum_{ij}(\rho_{ij}-\rho_{ij}(\mbox{\boldmath $\theta$}))^2w_{ij} \f] where \f$ \mbox{\boldmath $\rho$} \f$ is the base of the
  correlation matrix of the Gaussian field and \f$ \mbox{\boldmath $\rho(\theta)$} \f$ is the base of the correlation matrix of
  the GMRF, and \f$w_{ij}>0\f$ are weights.\n\n

  See Rue and Held (2005), Ch. 5, for further details.

  \par Example:
  \verbinclude example-doxygen-geo.txt

*/

#include <math.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#if !defined(__FreeBSD__)
#include <malloc.h>
#endif

#include "GMRFLib/GMRFLib.h"
#include "GMRFLib/GMRFLibP.h"
#include "GMRFLib/hashP.h"

#include "GMRFLib/geo-coefs2.h"				       /* these contains the BIG tables */
#include "GMRFLib/geo-coefs3.h"

#ifndef HGVERSION
#define HGVERSION
#endif
static const char RCSId[] = "file: " __FILE__ "  " HGVERSION;

/* Pre-hg-Id: $Id: geo.c,v 1.42 2008/08/26 07:07:12 hrue Exp $ */

static GMRFLib_Global_geo_tp Geo_coef;			       /* hash-table for the coefficient lookup */

#pragma omp threadprivate(Geo_coef)

/*!
  \brief Initializes and specifies a \c GMRFLib_geo_problem_tp -object.
 
  \param[in,out] geo_problem At output, <em>(*geo_problem)</em> is a pointer to
  a \c GMRFLib_geo_problem_tp -object, initialised and defined according to the
  problem specification.
  \param[in] name Name of the covariance function in the Gaussian field
  - #GMRFLib_CORTP_MATERN : Matern CF
  - #GMRFLib_CORTP_EXP : Exponential CF
  - #GMRFLib_CORTP_GAUSS : Gaussian CF
  \param[in] neigh Neighbourhood parameter (<em>m</em>)
  - \a neigh = 2 corresponds to \f$5\times5\f$ neighbourhood
  - \a neigh = 3 corresponds to \f$7\times7\f$ neighbourhood.
  \param[in] param Parameter in covariance function in the Gaussian field (\f$ =\nu \f$
             for Matern CF). For the Exponential and Gaussian, the value of \c param is not used.
  \param[in] range Range (for scaling of the Euclidean distance) (\a range >= 0)
  \param[in] nrow Number of rows in the lattice (or torus) (\a nrow > 0)
  \param[in] ncol Number of columns in the lattice (or torus) (\a ncol > 0)
  \param[in] prec Precision. If \c NULL, then \a prec = 1.0.
  \param[in] cyclic_flag If TRUE, then the graph is made cyclic (using a torus)
                                                                                                    
  \par Example:
  See \ref ex_geo.
  \sa GMRFLib_free_geo_problem(), GMRFLib_revise_geo_problem()
*/
int GMRFLib_init_geo_problem(GMRFLib_geo_problem_tp ** geo_problem, int name, int neigh,
			     double param, double range, int nrow, int ncol, double *prec, int cyclic_flag)
{
	GMRFLib_geoQfunc_arg_tp *geo_arg;
	GMRFLib_graph_tp *graph;
	map_ii *hash_table;
	int i, j, jj, row_i, col_i, row_j, col_j, d_row, d_col, d, *nelm, n;
	int map[] = { 0, 1, 2, 0, 3, 4, 0, 0, 5, 6, 7, 0, 0, 8, 0, 0, 0, 0, 9 };
	static double one = 1.0;			       /* yes, need a static one */

	n = nrow * ncol;

	/*
	 * Make the argument list to the geo-function 
	 */
	geo_arg = Calloc(1, GMRFLib_geoQfunc_arg_tp);
	geo_arg->name = name;
	geo_arg->neigh = neigh;
	geo_arg->param = param;
	geo_arg->range = range;
	geo_arg->nrow = nrow;
	geo_arg->ncol = ncol;
	geo_arg->cyclic_flag = cyclic_flag;
	geo_arg->prec = (prec ? prec : &one);

	switch (neigh) {
	case 2:
		GMRFLib_EWRAP0(GMRFLib_get_geo_coefs2(&(geo_arg->coef), name, param, range));
		break;
	case 3:
		GMRFLib_EWRAP0(GMRFLib_get_geo_coefs3(&(geo_arg->coef), name, param, range));
		break;
	default:
		GMRFLib_ASSERT(0, GMRFLib_ESNH);
	}

	/*
	 * make the graph 
	 */
	GMRFLib_make_lattice_graph(&graph, nrow, ncol, neigh, neigh, cyclic_flag);

	/*
	 * allocate the hash-table. allocate it with the correct number of elements 
	 */
	nelm = Calloc(n, int);

	for (i = 0; i < n; i++) {
		nelm[i] = 1;				       /* the diagonal */
		for (jj = 0; jj < graph->nnbs[i]; jj++)
			if (graph->nbs[i][jj] > i)
				nelm[i]++;
	}
	hash_table = Calloc(n, map_ii);
	for (i = 0; i < n; i++)
		map_ii_init_hint(&hash_table[i], nelm[i]);
	Free(nelm);

	for (i = 0; i < n; i++) {
		map_ii_set(&hash_table[i], i, map[0]);
		GMRFLib_node2lattice(i, &row_i, &col_i, nrow, ncol);
		for (j = 0; j < graph->nnbs[i]; j++) {
			jj = graph->nbs[i][j];
			if (jj > i) {			       /* since symmetric */
				GMRFLib_node2lattice(jj, &row_j, &col_j, nrow, ncol);
				d_row = IABS(row_i - row_j);
				d_col = IABS(col_i - col_j);
				if (cyclic_flag) {
					d_row = IMIN(d_row, nrow - d_row);
					d_col = IMIN(d_col, ncol - d_col);
				}
				d = ISQR(d_row) + ISQR(d_col);
				map_ii_set(&hash_table[i], jj, map[d]);
			}
		}
	}

	geo_arg->hash_table = hash_table;
	*geo_problem = Calloc(1, GMRFLib_geo_problem_tp);
	(*geo_problem)->graph = graph;
	(*geo_problem)->Qfunc = GMRFLib_geoQfunc;
	(*geo_problem)->Qfunc_arg = (void *) geo_arg;

	return GMRFLib_SUCCESS;
}

/*!
  \brief Revise a \c GMRFLib_geo_problem_tp -object.
 
  This function revise a \c GMRFLib_geo_problem_tp -object created by \c GMRFLib_init_geo_problem() without the need to
  rebuilding the graph and internal hash-tables. The parameters that can be revised are \c name, \c param, \c range and \c prec,
  which are the same as in \c GMRFLib_init_geo_problem(), i.e. the dimension, cyclic_flag and the size of the neighbourhood must
  be the same.

  \param[in,out] geo_problem At output, a revised \c GMRFLib_geo_problem_tp -object.
  \param[in] name Name of the covariance function in the Gaussian field
  - #GMRFLib_CORTP_MATERN : Matern CF
  - #GMRFLib_CORTP_EXP : Exponential CF
  - #GMRFLib_CORTP_GAUSS : Gaussian CF
  \param[in] param Parameter in covariance function in the Gaussian field (\f$ =\nu \f$
  for Matern CF).
  \param[in] range Range (for scaling of the Euclidean distance) (\a range >= 0)
  \param[in] prec Precision. If \c NULL, then \a prec = 1.0.
                                                                                                    
  \par Example:
  See \ref ex_geo.

  \sa GMRFLib_free_geo_problem(), GMRFLib_init_geo_problem()
*/
int GMRFLib_revise_geo_problem(GMRFLib_geo_problem_tp * geo_problem, int name, double param, double range, double *prec)
{
	static double one = 1.0;
	GMRFLib_geoQfunc_arg_tp *geo_arg;

	if (!geo_problem)
		return GMRFLib_SUCCESS;

	geo_arg = (GMRFLib_geoQfunc_arg_tp *) geo_problem->Qfunc_arg;

	geo_arg->name = name;
	geo_arg->param = param;
	geo_arg->range = range;
	geo_arg->prec = (prec ? prec : &one);

	switch (geo_arg->neigh) {
	case 2:
		GMRFLib_EWRAP0(GMRFLib_get_geo_coefs2(&(geo_arg->coef), name, param, range));
		break;
	case 3:
		GMRFLib_EWRAP0(GMRFLib_get_geo_coefs3(&(geo_arg->coef), name, param, range));
		break;
	default:
		GMRFLib_ASSERT(0, GMRFLib_ESNH);
	}

	return GMRFLib_SUCCESS;
}
double GMRFLib_geoQfunc(int node, int nnode, void *arg)
{
	int idx;
	GMRFLib_geoQfunc_arg_tp *geo_arg = (GMRFLib_geoQfunc_arg_tp *) arg;

	map_ii_get(&(geo_arg->hash_table)[IMIN(node, nnode)], IMAX(node, nnode), &idx);
	return *(geo_arg->prec) * geo_arg->coef[idx];
}

/*!
  \brief Deallocates all allocated arrays in a \c GMRFLib_geo_problem_tp -object.
*/
int GMRFLib_free_geo_problem(GMRFLib_geo_problem_tp * geo_problem)
{
	int i;
	GMRFLib_geoQfunc_arg_tp *arg;

	if (!geo_problem)
		return GMRFLib_SUCCESS;

	GMRFLib_free_graph(geo_problem->graph);
	arg = (GMRFLib_geoQfunc_arg_tp *) (geo_problem->Qfunc_arg);
	for (i = 0; i < arg->nrow * arg->ncol; i++)
		map_ii_free(&(arg->hash_table[i]));
	Free(arg->hash_table);
	Free(geo_problem->Qfunc_arg);
	geo_problem->Qfunc_arg = NULL;

	return GMRFLib_SUCCESS;
}
const char *GMRFLib_geo_translate_cortp(int name)
{
	/*
	 * return the cortp as text 
	 */
	const char *map[] = { "Matern", "Exponential", "Gaussian" };

	switch (name) {
	case GMRFLib_CORTP_MATERN:
		return map[0];
	case GMRFLib_CORTP_EXP:
		return map[1];
	case GMRFLib_CORTP_GAUSS:
		return map[2];
	default:
		GMRFLib_ASSERT_RETVAL(0, GMRFLib_ESNH, NULL);
	}
	return NULL;
}
const char *GMRFLib_geo_translate_neigh(int neigh)
{
	/*
	 * return the neigh as text 
	 */
	const char *map[] = { "5x5", "7x7" };

	switch (neigh) {
	case 2:
		return map[0];
	case 3:
		return map[1];
	default:
		GMRFLib_ASSERT_RETVAL(0, GMRFLib_ESNH, NULL);
	}
	return NULL;
}
int GMRFLib_get_geo_coefs2(double **coef, int name, double param, double range)
{
	char *key, *ind_s;
	mapkit_error err;
	int i, n;
	static int first = 1;

#pragma omp threadprivate(first)

	/*
	 * param not in use for these, but must be zero 
	 */
	if (name == GMRFLib_CORTP_GAUSS || name == GMRFLib_CORTP_EXP)
		param = 0.0;

	if (first) {
		/*
		 * build hash-table 
		 */

		n = sizeof(Geo_coefs2) / sizeof(Geo_coefs2[0]);
		map_stri_init_hint(&(Geo_coef.coef2), n);
		for (i = 0; i < n; i++) {
			GMRFLib_EWRAP0(GMRFLib_sprintf(&key, GMRFLib_GEO_FMT, Geo_coefs2[i].name, Geo_coefs2[i].param, Geo_coefs2[i].range));
			map_stri_set(&(Geo_coef.coef2), key, i);
		}
		first = 0;
	}

	GMRFLib_sprintf(&ind_s, GMRFLib_GEO_FMT, name, param, range);
	err = map_stri_get(&(Geo_coef.coef2), ind_s, &i);
	Free(ind_s);

	if (err) {
		*coef = NULL;
		GMRFLib_EWRAP0(GMRFLib_sprintf(&ind_s, "%s %s%s %5.3f%s %5.3f%s", "Coefficients not available (name =",
					       GMRFLib_geo_translate_cortp(name), ", param =", param, ", range =", range, ")"));
		GMRFLib_ERROR_MSG(GMRFLib_EGEOCOOF, ind_s);
		Free(ind_s);
		return GMRFLib_EGEOCOOF;
	} else {
		*coef = &(Geo_coefs2[i].coefs[0]);
		return GMRFLib_SUCCESS;
	}

	return GMRFLib_SUCCESS;
}
int GMRFLib_get_geo_coefs3(double **coef, int name, double param, double range)
{
	char *key, *ind_s;
	mapkit_error err;
	int i, n;
	static int first = 1;

#pragma omp threadprivate(first)

	/*
	 * param not in use for these, but must be zero 
	 */
	if (name == GMRFLib_CORTP_GAUSS || name == GMRFLib_CORTP_EXP)
		param = 0.0;

	if (first) {
		/*
		 * build hash-table 
		 */

		n = sizeof(Geo_coefs3) / sizeof(Geo_coefs3[0]);
		map_stri_init_hint(&(Geo_coef.coef3), n);
		for (i = 0; i < n; i++) {
			GMRFLib_EWRAP0(GMRFLib_sprintf(&key, GMRFLib_GEO_FMT, Geo_coefs3[i].name, Geo_coefs3[i].param, Geo_coefs3[i].range));
			map_stri_set(&(Geo_coef.coef3), key, i);
		}
		first = 0;
	}

	GMRFLib_sprintf(&ind_s, GMRFLib_GEO_FMT, name, param, range);
	err = map_stri_get(&(Geo_coef.coef3), ind_s, &i);
	Free(ind_s);

	if (err != MAPKIT_OK) {
		*coef = NULL;
		GMRFLib_EWRAP0(GMRFLib_sprintf(&ind_s, "%s %s%s %5.3f%s %5.3f%s", "Coefficients not available (name =",
					       GMRFLib_geo_translate_cortp(name), ", param =", param, ", range =", range, ")"));
		GMRFLib_ERROR_MSG(GMRFLib_EGEOCOOF, ind_s);
		Free(ind_s);
		return GMRFLib_EGEOCOOF;
	} else {
		*coef = &(Geo_coefs3[i].coefs[0]);
		return GMRFLib_SUCCESS;
	}

	return GMRFLib_SUCCESS;
}

/*!
  \brief Return the coefficients the given parameters.

  \param[out] coef The coefficients are return in \a *coef.
  \param[in] name Name of the covariance function in the Gaussian field
  - #GMRFLib_CORTP_MATERN : Matern CF
  - #GMRFLib_CORTP_EXP : Exponential CF
  - #GMRFLib_CORTP_GAUSS : Gaussian CF
  \param[in] neigh Neighbourhood parameter (<em>m</em>)
  - \a neigh = 2 corresponds to \f$5\times5\f$ neighbourhood
  - \a neigh = 3 corresponds to \f$7\times7\f$ neighbourhood.
  \param[in] param Parameter in covariance function in the Gaussian field (\f$ =\nu \f$
  for Matern CF).
  \param[in] range Range (for scaling of the Euclidean distance) (\a range >= 0)

  \return Returns in \a *coef a pointer to the values of the coefficients if they exist in the library for the given parameters,
  otherwise \a *coef is \c NULL. The error-code is either #GMRFLib_SUCCESS or an error.

  \sa GMRFLib_is_geo_coefs.
*/
int GMRFLib_get_geo_coefs(double **coef, int name, int neigh, double param, double range)
{
	*coef = NULL;					       /* default */
	switch (neigh) {
	case 2:
		GMRFLib_EWRAP0(GMRFLib_get_geo_coefs2(coef, name, param, range));
		break;
	case 3:
		GMRFLib_EWRAP0(GMRFLib_get_geo_coefs3(coef, name, param, range));
		break;
	default:
		GMRFLib_ASSERT(0, GMRFLib_ESNH);
		return GMRFLib_ESNH;
	}

	return (*coef ? GMRFLib_SUCCESS : GMRFLib_EGEOCOOF);
}

/*!
  \brief Check if coefficients are available for the given parameters.

  \param[in] name Name of the covariance function in the Gaussian field
  - #GMRFLib_CORTP_MATERN : Matern CF
  - #GMRFLib_CORTP_EXP : Exponential CF
  - #GMRFLib_CORTP_GAUSS : Gaussian CF
  \param[in] neigh Neighbourhood parameter (<em>m</em>)
  - \a neigh = 2 corresponds to \f$5\times5\f$ neighbourhood
  - \a neigh = 3 corresponds to \f$7\times7\f$ neighbourhood.
  \param[in] param Parameter in covariance function in the Gaussian field (\f$ =\nu \f$
  for Matern CF).
  \param[in] range Range (for scaling of the Euclidean distance) (\a range >= 0)

  \return This function returns \c GMRFLib_TRUE if the coefficients exist in the library for the given parameters, otherwise it
  returns \c GMRFLib_EGEOCOOF .\n

  \sa GMRFLib_get_geo_coefs.
*/
int GMRFLib_is_geo_coefs(int name, int neigh, double param, double range)
{
	/*
	 * turn off the error-system and check the error-code. restore it afterwards 
	 */
	int retval;
	double *coef;
	GMRFLib_error_handler_tp *ehandler = GMRFLib_set_error_handler_off();

	retval = GMRFLib_get_geo_coefs(&coef, name, neigh, param, range);
	GMRFLib_set_error_handler(ehandler);		       /* restore it */

	return (retval == GMRFLib_SUCCESS ? GMRFLib_SUCCESS : GMRFLib_EGEOCOOF);
}

/*!
  \brief Prints the parameters for which coefficients are available.

  \param[in] fp The FILE* on which to print the parameters. If \c fp is \c NULL, then use \c stdout.
*/
int GMRFLib_print_geo_coefs(FILE * fp)
{
	int i;

	if (!fp)
		fp = stdout;

	fprintf(fp, "\nCoefficients available:\n\n");
	i = 0;
	while (strlen(Geo_desc2[i].desc))
		fprintf(fp, "%s\n", Geo_desc2[i++].desc);
	i = 0;
	while (strlen(Geo_desc3[i].desc))
		fprintf(fp, "%s\n", Geo_desc3[i++].desc);

	return GMRFLib_SUCCESS;
}

/*
  
  Example for manual
 */

/*! \page ex_geo An example using GMRFs as approximation to Gaussian fields
  
  Consider a Gaussian field on a \f$ 40\times40 \f$ torus. We use a Matern CF with
  parameter \f$ \nu=0.5 \f$, with range = 15, and precision = 1.0.
  Further, we use a GMRF with a \f$ 5\times5 \f$ neighbourhood that is fitted to the
  Gaussian field.
  Figure 1 (a) shows the approximative precision matrix
  \f$ \mbox{\boldmath $Q(\theta)$} \f$, where the coefficients are given by
  \f[ \mbox{\boldmath $\theta$} = \left[\begin{array}{ccc} & & 0.049183 \\
  & -0.005295 & -0.059234 \\ 1.0 & -0.289775 & 0.114575 \end{array}\right] \f]  
  Figure 1 (b) shows the precision matrix for the Gaussian field. Note that the
  precision matrix for the Gaussian field is a full matrix, while the precision matrix
  for the GMRF is quite sparse. \n
  This demonstrates that using GMRFs as approximation to Gaussian fields can be quite
  effective. The approximated CF is also accurate, which can be seen in Figure 2.\n
  Figure 3 shows a sample from the GMRF (a) and a sample from the Gaussian
  field (b). \n\n
  \htmlonly
  <table>
  <tr><td align="center">
    <table>
    <tr><td width="350" align="center"><img src="figs/geo-Q.gif" width="300" height="300"></td>  
    <td width="350" align="center"><img src="figs/geo-Q-matern.gif" width="300" height="300"></td></tr>
    </table></td></tr>
  <tr><td align="center">
    <table>
    <tr><td width="350" align="center">(a)</td>  
    <td width="350" align="center">(b)</td></tr>
    </table></td></tr>
  </table>
  \endhtmlonly
  <b>Figure 1</b> (a) Approximative precision matrix \f$ \mbox{\boldmath $Q(\theta)$} \f$
  with 40000 nonzero elements.
  (b) Precision matrix with Matern CF with 2560000 nonzero elements.\n\n

  \htmlonly
  <table>
  <tr><td><img src="figs/CFMatern15.gif" width="500"> </td></tr>
  </table> 
  \endhtmlonly
  <b>Figure 2</b> Correlation function for the fitted GMRF (red line) and the Matern CF
  (blue line) with range 15 and \f$ 5\times5 \f$ neighbourhood.\n\n
  
  \htmlonly
  <table>
  <tr><td align="center">
    <table>
    <tr><td width="350" align="center"><img src="figs/geo-sample.gif" width="300" height="300"></td>  
    <td width="350" align="center"><img src="figs/geo-sample-matern.gif" width="300" height="300"></td></tr>
    </table></td></tr>
  <tr><td align="center">
    <table>
    <tr><td width="350" align="center">(a)</td>  
    <td width="350" align="center">(b)</td></tr>
    </table></td></tr>
  </table>
  \endhtmlonly
  <b>Figure 3</b> (a) Sample from GMRF using <em>\b Q</em> in Figure 1 (a).
  (b) Sample from Gaussian field using <em>\b Q</em> in Figure 1 (b). \n\n
  
  \par Program code:
 
  \verbinclude example-doxygen-geo.txt
 
*/
