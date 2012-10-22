
/* ar.h
 * 
 * Copyright (C) 2012 Havard Rue
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
#ifndef __INLA_AR_H__
#define __INLA_AR_H__
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
 *
 */

#include "inla.h"

typedef struct
{
	int n;
	int p;

	double **log_prec;
	double ***pacf_intern;
	
	/* 
	 * these are the stored values, all for prec = 1
	 */
	double **hold_pacf_intern;			       /* [i][id][0] */
	double **hold_Q;				       /* dim = 2*p + 1 */
	double **hold_Qmarg;				       /* dim = p */
}
	ar_def_tp;


int ar_pacf2phi(int p, double *pacf, double *phi);
int ar_phi2pacf(int p, double *phi, double *pacf);
int ar_test1();
int ar_marginal_distribution(int p, double *pacf, double *prec, double *Q);
double Qfunc_ar(int i, int j, void *arg);
double ar_map_pacf(double arg, map_arg_tp typ, void *param);


__END_DECLS
#endif
