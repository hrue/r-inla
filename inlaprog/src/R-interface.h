
/* R-interface.h
 * 
 * Copyright (C) 2014 Havard Rue
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
#ifndef __INLA_R_INTERFACE_H__
#define __INLA_R_INTERFACE_H__
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
int inla_R_assign(const char *variable, int n, double *x);
int inla_R_funcall1(int *n_out, double **x_out, const char *function, int n, double *x);
int inla_R_funcall2(int *n_out, double **x_out, const char *function, const char *tag, int n, double *x);
int inla_R_get(int *n_out, double **x_out, const char *variable);
int inla_R_init(void);
int inla_R_inlaload(const char *filename);
int inla_R_library(const char *library);
int inla_R_load(const char *filename);
int inla_R_rgeneric(int *n_out, double **x_out, const char *cmd, const char *model, int n, double *theta);
int inla_R_source(const char *filename);
void inla_R_exit(void);

__END_DECLS
#endif
