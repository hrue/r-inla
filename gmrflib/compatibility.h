
/* compatibility.h
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
 *       H{\aa}vard Rue
 *       Department of Mathematical Sciences
 *       The Norwegian University of Science and Technology
 *       N-7491 Trondheim, Norway
 *       Voice: +47-7359-3533    URL  : http://www.math.ntnu.no/~hrue
 *       Fax  : +47-7359-3524    Email: havard.rue@math.ntnu.no
 *
 * RCSId: $Id: compatibility.h,v 1.21 2010/03/12 12:23:16 hrue Exp $
 *
 */

/*!
  \file compatibility.h
  \brief Include defines to ensure (partially) backwards-compatbility of GMRFLib.

  This file defines `old names' to the `new name' to ensure (partially) backwards-compatbility of GMRFLib.
*/

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
#define GMRFLib_create_lattice  GMRFLib_make_lattice_graph
#define GMRFLib_create_graph    GMRFLib_make_empty_graph
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
