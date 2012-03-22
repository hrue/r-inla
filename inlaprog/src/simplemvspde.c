
/* simplemvspde.c
 * WRITTEN BY DAN.  DON'T EXPECT THIS TO WORK! 
 * Copyright (C) 2011  Havard Rue
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
 */
#ifndef HGVERSION
#define HGVERSION
#endif
static const char RCSId[] = "file: " __FILE__ "  " HGVERSION;


#if !defined(__FreeBSD__)
#include <malloc.h>
#endif
#include <stdlib.h>

#include "GMRFLib/GMRFLib.h"
#include "GMRFLib/GMRFLibP.h"
#include "inla.h"
#include "simplemvspde.h"


extern G_tp G;						       /* import some global parametes from inla */


double inla_simplemvspde_Qfunction(int node, int nnode, void *arg) {

  inla_simplemvspde_tp *model = (inla_simplemvspde_tp *) arg;
  double value =0.0;
  double kappa2[3], b[3];
  int i_;
  for (i_ =0; i_<3; i_++) {
    kappa2[i] = model->theta[i][GMRFLib_thread_id][0];
    b[i] = model->theta[i+3][GMRFLib_thread_id][0];
  }

  int i,j;
  
  
  int n = model->n_block;

  if ( (node <n ) && (nnode < n)) {
    i = node;
    j = nnode;  
    
    // modify for more complex models.  Simplifying assumption: alpha === 2 and mesh is the same everywhere. That means that each block consists of linear combinations of the same 3 matrices!!
    value = b[0]* (kappa2[0]*kappa2[0] * GMRFLib_matrix_get(i, j, model->M[0]) +
			  kappa2[0]* GMRFLib_matrix_get(i,j, model->M[1]) +
			  kappa2[0]* GMRFLib_matrix_get(j,i, model->M[1]) +
			  GMRFLib_matrix_get(i,j, model->M[2])  ) + 
            b[1]* (kappa2[1]*kappa2[1]* GMRFLib_matrix_get(i, j, model->M[0]) +
		   kappa2[1]* GMRFLib_matrix_get(i,j, model->M[1]) +
		   kappa2[1]*GMRFLib_matrix_get(j,i,model->M[1])+ 
		   GMRFLib_matrix_get(i,j, model->M[2])  );
    

  } else if ( (node >= n) && (nnode >= n) ) {

    i = node-n;  //transform coordinates to i,j
    j = nnode - n;  
    value = b[2]* (kappa2[2]*kappa2[2]* GMRFLib_matrix_get(i, j, model->M[0]) +
		   kappa2[2]* GMRFLib_matrix_get(i,j, model->M[1]) + 
		   kappa2[2]* GMRFLib_matrix_get(j,i, model->M[1]) +
		   GMRFLib_matrix_get(i,j, model->M[2])  );

  } else {
    // what the fuck is mod in c?  God i hate this language!
    i = node;
    j = nnode;
    if (i  >= n) {
      i = i-n;
    } else {
      j = j-n;
    }

  
    
  

	value =  b[1]*b[2]*(kappa2[1]*kappa2[2] * GMRFLib_matrix_get(i, j, model->M[0]) +
				   kappa2[1] * GMRFLib_matrix_get(i, j, model->M[1]) +
				   kappa2[2] * GMRFLib_matrix_get(j, i, model->M[1]) + 
				   GMRFLib_matrix_get(i, j, model->M[2]));

	
  }

  return value;
}


int inla_simplemvspde_build_model(inla_simplemvspde_tp ** smodel, const char *prefix)
{

  inla_simplemvspde_tp *model = NULL;
  char *fnm = NULL;

  model = Calloc(1, inla_simplemvspde_tp);

  model->M = Calloc(3, GMRFLib_matrix_tp *);
       for (i = 0; i < 3; i++) {
		GMRFLib_sprintf(&fnm, "%s%s%1d", prefix, "M", i);
		model->M[i] = GMRFLib_read_fmesher_file((const char *) fnm, 0, -1);
	}

	for (i = 1; i < 3; i++) {
		/*
		 * all are square with the same dimension n x n 
		 */
		assert(model->M[i]->nrow == model->M[i]->ncol);
		assert(model->M[0]->nrow == model->M[i]->nrow);

	}

	model->n = model->M[0]->nrow*2;
	model->n_block = model->M[0]->nrow;
	model->ntheta = 6; //Hard coded!!



	/*
	 * I need to build the graph. Need to add both M_ij and M_ji as M1 can be non-symmetric. 
	 */
	GMRFLib_ged_tp *ged = NULL;
	GMRFLib_ged_init(&ged, NULL);

#define ADD_GRAPH(_G)							\
	{								\
		int i_, j_, jj_;					\
		for(i_ = 0; i_ < _G->n; i_++) {				\
			GMRFLib_ged_add(ged, i_, i_);			\
			for(jj_ = 0; jj_ < _G->nnbs[i_]; jj_++){	\
				j_ = _G->nbs[i_][jj_];			\
				GMRFLib_ged_add(ged, i_, j_);		\
				GMRFLib_ged_add(ged, j_, i_);		\
			}						\
		}							\
	}

	for (i = 0; i < 3; i++) {
		ADD_GRAPH(model->M[i]->graph);
	}

#undef ADD_GRAPH
	GMRFLib_graph_tp * tmp_graph = NULL;
	GMRFLib_ged_build(&tmp_graph,ged);

	GMRFLib_ged_tp *ged2 = NULL;
	GMRFLib_ged_init(&ged2, tmp_graph);
	GMRFLib_ged_append_graph(ged,tmp_graph);
	
	int i_,j_,jj_;
	for (i_ = 0; i_ < tmp_graph->n; i_++) {
	  for (jj_=0; jj_ < tmp_graph->nnbs[i_]; jj_++) {
	    j_ = tmp_graph->nbs[i_][jj_];
	    GMRFLib_ged_add(ged2,i_, j_+n);
	    GMRFLib_ged_add(ged2, i_+n, j_);  //I think this is correct.  Add in the top right and bottom left block.
	  }
	}

	GMRFLib_ged_build(&(model->graph), ged2);
	assert(model->n == model->graph->n);
	GMRFLib_ged_free(ged);
	GMRFLib_ged_frree(ged2);  

	model->Qfunc = inla_simplemvspde_Qfunction;
	model->Qfunc_arg = (void *) model;

	HYPER_NEW2(model->theta, 0.0, model->ntheta);
	*smodel = model;

	return INLA_OK;
}



}
