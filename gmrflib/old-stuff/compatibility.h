#ifndef __GMRFLib_COMPATBILITY_H__
#define __GMRFLib_COMPATBILITY_H__

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
   compability with versions < 2.0
*/
#define GMRFLib_KEEP_bchol      GMRFLib_KEEP_chol
#define GMRFLib_create_lattice  GMRFLib_graph_mk_lattice
#define GMRFLib_create_graph    GMRFLib_graph_mk_empty
#define GMRFLib_create_constr   GMRFLib_make_empty_constr

/*
  compability with version = 2.0
*/
#define GMRFLib_hidden_approx GMRFLib_init_problem_hidden
#define GMRFLib_blockupdate2  GMRFLib_blockupdate_hidden

/*
  compability with versions <= 3.0 and changes made in 3.0
*/
#define GMRFLib_ctime GMRFLib_cpu
#define GMRFLib_ai_marginal_hidden_int GMRFLib_ai_INLA
#define GMRFLib_AI_LINEAR_CORRECTION_EXACT GMRFLib_AI_LINEAR_CORRECTION_FAST

/* 
   changes made
 */
#define GMRFLib_optimise_reorder(g) GMRFLib_optimize_reorder(g, NULL)

/*
  a comment to prevent funny indentation
*/
    __END_DECLS
#endif
