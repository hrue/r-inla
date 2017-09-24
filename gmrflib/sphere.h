
/* sphere.h
 * 
 * Copyright (C) 2005-2006 Havard Rue
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
  \file sphere.h
  \brief Typedefs used to define spherical IGMRFs
*/

#ifndef __GMRFLib_SPHERE_H__
#define __GMRFLib_SPHERE_H__

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
  \struct GMRFLib_coord_tp sphere.h
  \brief  The typedef used to defined coordinates of each node of the spherical iGMRFs
*/

/*
 */
    typedef struct {

	/**
	 *  \brief The \c x-coordinate 
	 */
	double x;

	/**
	 *  \brief The \c y-coordinate 
	 */
	double y;

	/**
	 *  \brief The \c z-coordinate 
	 */
	double z;
} GMRFLib_sphere_coord_tp;

/*!
  \struct GMRFLib_sphere_tp sphere.h
  \brief The typedef used to define spherical iGMRFs on a unit sphere
*/
typedef struct {

	/**
	 *  \brief The graph of the spherical IGMRF 
	 */
	GMRFLib_graph_tp *graph;

	/**
	 *  \brief The Qfunction of the spherical IGMRF 
	 */
	GMRFLib_Qfunc_tp *Qfunc;

	/**
	 *  \brief The arguments to GMRFLib_sphere_tp::Qfunc 
	 */
	void *Qfunc_arg;

	/**
	 *  \brief The areas of the corrsponding spherical triangles divided by 3 
	 */
	double *area;

	/**
	 *  \brief The coordinates of each node.
	 * 
	 * For each node\c i, then \c coord[i].x, \c coord[i].y and \c coord[i].z, gives its \c x, \c y and \c z coordinate. 
	 */
	GMRFLib_sphere_coord_tp *coord;
} GMRFLib_sphere_tp;

int GMRFLib_make_spherical_igmrf(GMRFLib_sphere_tp ** sphere, int resolution, double *scale, double *log_scale, double **log_scale_omp);
int GMRFLib_free_spherical_igmrf(GMRFLib_sphere_tp * sphere);

__END_DECLS
#endif
