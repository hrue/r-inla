
/* R-interface.h
 * 
 * Copyright (C) 2014-2024 Havard Rue
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
    typedef enum {
	INLA_R_INIT = 1,
	INLA_R_ASSIGN,
	INLA_R_FUNCALL1,
	INLA_R_FUNCALL2,
	INLA_R_FUNCALL_JP,
	INLA_R_GET,
	INLA_R_INLALOAD,
	INLA_R_LIBRARY,
	INLA_R_LOAD,
	INLA_R_RGENERIC,
	INLA_R_SOURCE,
	INLA_R_EXIT
} inla_R_cmd_tp;

void inla_set_R_home(char *home);

int inla_R_do_(inla_R_cmd_tp cmd, void *a1, void *a2, void *a3, void *a4, void *a5, void *a6);

#define inla_R_assign(a1, a2, a3) inla_R_do_(INLA_R_ASSIGN, (void *) (a1), (void *) (a2), (void *) (a3), (void *) NULL, (void *) NULL, (void *) NULL)
int inla_R_assign_(const char *variable, int *n, double *x);

#define inla_R_funcall1(a1, a2, a3, a4, a5) inla_R_do_(INLA_R_FUNCALL1, (void *) (a1), (void *) (a2), (void *) (a3), (void *) (a4), (void *) (a5), (void *) NULL)
int inla_R_funcall1_(int *n_out, double **x_out, const char *function, int *n, double *x);

#define inla_R_funcall2(a1, a2, a3, a4, a5, a6) inla_R_do_(INLA_R_FUNCALL2, (void *) (a1), (void *) (a2), (void *) (a3), (void *) (a4), (void *) (a5), (void *) (a6))
int inla_R_funcall2_(int *n_out, double **x_out, const char *function, const char *tag, int *n, double *x);

#define inla_R_funcall_jp(a1, a2, a3, a4, a5, a6) inla_R_do_(INLA_R_FUNCALL_JP, (void *) (a1), (void *) (a2), (void *) (a3), (void *) (a4), (void *) (a5), (void *) (a6))
int inla_R_funcall_jp_(int *n_out, double **x_out, const char *function, int *n, double *x, void *sexp);

#define inla_R_get(a1, a2, a3) inla_R_do_(INLA_R_GET, (void *) (a1), (void *) (a2), (void *) (a3), (void *) NULL, (void *) NULL, (void *) NULL)
int inla_R_get_(int *n_out, double **x_out, const char *variable);

#define inla_R_init() inla_R_do_(INLA_R_INIT, (void *) NULL, (void *) NULL, (void *) NULL, (void *) NULL, (void *) NULL, (void *) NULL)
int inla_R_init_(void);

#define inla_R_inlaload(a1) inla_R_do_(INLA_R_INLALOAD, (void *) (a1), (void *) NULL, (void *) NULL, (void *) NULL, (void *) NULL, (void *) NULL)
int inla_R_inlaload_(const char *filename);

#define inla_R_library(a1) inla_R_do_(INLA_R_LIBRARY, (void *) (a1), (void *) NULL, (void *) NULL, (void *) NULL, (void *) NULL, (void *) NULL)
int inla_R_library_(const char *filename);

#define inla_R_load(a1) inla_R_do_(INLA_R_LOAD, (void *) (a1), (void *) NULL, (void *) NULL, (void *) NULL, (void *) NULL, (void *) NULL)
int inla_R_load_(const char *filename);

#define inla_R_rgeneric(a1, a2, a3, a4, a5, a6) inla_R_do_(INLA_R_RGENERIC, (void *) (a1), (void *) (a2), (void *) (a3), (void *) (a4), (void *) (a5), (void *) (a6))
int inla_R_rgeneric_(int *n_out, double **x_out, const char *cmd, const char *model, int *n, double *theta);

#define inla_R_source(a1) inla_R_do_(INLA_R_SOURCE, (void *) (a1), (void *) NULL, (void *) NULL, (void *) NULL, (void *) NULL, (void *) NULL)
int inla_R_source_(const char *filename);
int inla_R_source_quiet_(const char *filename);

#define inla_R_exit() inla_R_do_(INLA_R_EXIT, (void *) NULL, (void *) NULL, (void *) NULL, (void *) NULL, (void *) NULL, (void *) NULL)
int inla_R_exit_(void);

void *inla_R_vector_of_strings(int n, char **s);

__END_DECLS
#endif
