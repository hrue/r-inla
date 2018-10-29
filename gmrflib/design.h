
/* design.h
 * 
 * Copyright (C) 2006 Havard Rue
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
  \file design.h
  \brief Typedefs and defs for design
*/

#ifndef __GMRFLib_DESIGN_H__
#define __GMRFLib_DESIGN_H__

#include <stdlib.h>
#include <stddef.h>
#include <math.h>
#if !defined(__FreeBSD__)
#include <malloc.h>
#endif

#include "GMRFLib/fmesher-io.h"

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

//

typedef struct {
	double **experiment;
	double *int_weight;
	int nexperiments;
	int nfactors;
	int std_scale;					       /* if true, then the weights are on a standardized scale. this
								* is the normal case */
} GMRFLib_design_tp;

int GMRFLib_get_design(GMRFLib_design_tp ** design, int nfactors);
int GMRFLib_read_design(GMRFLib_design_tp ** design, GMRFLib_matrix_tp *D, int std_scale);
int GMRFLib_free_design(GMRFLib_design_tp * design);
int GMRFLib_print_design(FILE * fp, GMRFLib_design_tp * design);

__END_DECLS
#endif
