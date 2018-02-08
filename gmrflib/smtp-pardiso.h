/* GMRFLib-smtp-pardiso.h
 * 
 * Copyright (C) 2018 Havard Rue & Alexander Litvinenko
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
  \file smtp-band.h
  \brief GMRFLib interface to the band-solver in the LAPACK library
*/

#ifndef __GMRFLib_SMTP_PARDISO_H__
#define __GMRFLib_SMTP_PARDISO_H__

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

typedef struct
{
	int n;
	int na;
	int *ia;
	int *ja;
	double *a;
}
	GMRFLib_csr_tp;

int GMRFLib_Q2csr(GMRFLib_csr_tp **csr, GMRFLib_graph_tp *graph, GMRFLib_Qfunc_tp *Qfunc, void *Qfunc_arg);
int GMRFLib_print_csr(FILE *fp, GMRFLib_csr_tp *csr);
int GMRFLib_csr2Q(GMRFLib_tabulate_Qfunc_tp **Qtab, GMRFLib_graph_tp **graph, GMRFLib_csr_tp *csr);

int pardiso_test();

__END_DECLS
#endif
