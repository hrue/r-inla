
/* simplemvspde.h
 * WRITTEN BY DAN - DO NOT EXPECT THIS TO WORK!
 * Copyright (C) 2011 Havard Rue
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
 *       H{\aa}vard Rue
 *       Department of Mathematical Sciences
 *       The Norwegian University of Science and Technology
 *       N-7491 Trondheim, Norway
 *       Voice: +47-7359-3533    URL  : http://www.math.ntnu.no/~hrue  
 *       Fax  : +47-7359-3524    Email: havard.rue@math.ntnu.no
 *
 *
 */
#ifndef __INLA_SIMPLEMVSPDE_H__
#define __INLA_SIMPLEMVSPDE_H__
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
#include "GMRFLib/GMRFLib.h"


typedef struct {
        int n;
  int n_block;
  int ntheta;  // kappa11, kappa21, kappa22, b11, b21, b22
	int debug;
	int theta_first_idx;

	GMRFLib_matrix_tp **M;

  double  ***theta;  //theta[i][thread_id][0]
	GMRFLib_Qfunc_tp *Qfunc;
	void *Qfunc_arg;
	GMRFLib_graph_tp *graph;
} inla_simplemvspde_tp;


double inla_simplemvspde_Qfunction(int node, int nnode, void *arg);
int inla_simplemvspde_build_model(inla_simplemvspde_tp ** smodel, const char *prefix);

__END_DECLS
#endif
