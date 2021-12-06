
/* cgeneric.h
 * 
 * Copyright (C) 2021 Havard Rue
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
 */
#ifndef __INLA_CGENERIC_H__
#define __INLA_CGENERIC_H__

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

typedef enum {
	INLA_CGENERIC_VOID = 0, 
	INLA_CGENERIC_Q, 
	INLA_CGENERIC_GRAPH, 
	INLA_CGENERIC_MU, 
	INLA_CGENERIC_INITIAL, 
	INLA_CGENERIC_LOG_NORM_CONST, 
	INLA_CGENERIC_LOG_PRIOR, 
	INLA_CGENERIC_QUIT
}
	inla_cgeneric_cmd_tp;

#define INLA_CGENERIC_CMD_NAME(cmd_) ((cmd_) == INLA_CGENERIC_VOID ? "void" : \
				      ((cmd_) == INLA_CGENERIC_Q ? "Q" : \
				       ((cmd_) == INLA_CGENERIC_GRAPH ? "graph" : \
					((cmd_) == INLA_CGENERIC_MU ? "mu" : \
					 ((cmd_) == INLA_CGENERIC_INITIAL ? "initial" : \
					  ((cmd_) == INLA_CGENERIC_LOG_NORM_CONST ? "log_norm_const" : \
					   ((cmd_) == INLA_CGENERIC_LOG_PRIOR ? "log_prior" : \
					    ((cmd_) == INLA_CGENERIC_QUIT ? "quit" : "(***ERROR***)"))))))))

typedef struct 
{
	int n_ints;
	int **ints;

	int n_doubles;
	double **doubles;

	int n_chars;
	char **chars;
}
	inla_cgeneric_arg_tp;

typedef double * inla_cgeneric_func_tp(inla_cgeneric_cmd_tp cmd, double *theta, inla_cgeneric_arg_tp * args);
inla_cgeneric_func_tp inla_cgeneric_demo; 

__END_DECLS
#endif
