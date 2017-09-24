
/* random.h
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
  \file random.h
  \brief Typedefs and defines for the RNG.
*/

#ifndef __GMRFLib_RANDOM_H__
#define __GMRFLib_RANDOM_H__

#if !defined(__FreeBSD__)
#include <malloc.h>
#endif

#include "GMRFLib/GMRFLib.h"

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
#include <stdlib.h>
#define GMRFLib_rng (GMRFLib_rng_ptr ? GMRFLib_rng_ptr :  GMRFLib_rng_init_default())
    gsl_rng * GMRFLib_rng_init_default(void);
double GMRFLib_rng_uniform(void);
int GMRFLib_rng_init(unsigned long int seed);
int GMRFLib_rng_set_default_seed(void);
int GMRFLib_rng_setstate(void *saved_state);
void *GMRFLib_rng_getstate(size_t *siz);


__END_DECLS
#endif
