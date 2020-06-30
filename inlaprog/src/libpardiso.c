/* libpardiso.c
 * 
 * Copyright (C) 2018 Havard Rue
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
#include <stdio.h>
#if !defined(__FreeBSD__)
#include <malloc.h>
#endif
#include <stdlib.h>

/* 
   this creates an empty version of the pardiso functions with an appropriate wrapper to the relevant original metis function.
   Since the (pardiso modified) metis library is bundled into the libpardiso.so, I have to do the same.
 */



// do not change: also GMRFLib/smtp-pardiso.c uses this code
#define NOLIB_ECODE (270465)

#define NO_PARDISO_LIB						\
	{							\
		fprintf(stderr, "\n\n\t*** No PARDISO library is loaded. Exit.\n\n");	\
		exit(1);						\
	}

void pardisoinit(void *a, int *b, int *c, int *d, double *e, int *f)
{
	*f = NOLIB_ECODE;
	return;
}
void pardiso(void *a, int *b, int *c, int *d, int *e, int *f, double *g,
	     int *h, int *i, int *j, int *k, int *l, int *m, double *n, double *o, int *p, double *q) NO_PARDISO_LIB;
void pardiso_chkmatrix(int *a, int *s, double *d, int *f, int *g, int *h) NO_PARDISO_LIB;
void pardiso_chkvec(int *a, int *s, double *d, int *f) NO_PARDISO_LIB;
void pardiso_printstats(int *a, int *s, double *d, int *f, int *g, int *h, double *j, int *k) NO_PARDISO_LIB;
void pardiso_get_factor_csc(void **a, double *s, int *d, int *f, double *g, int *h, int *j, int *k, int l) NO_PARDISO_LIB;
void pardiso_get_inverse_factor_csc(void **a, double *s, int *d, int *f, int *g, int h) NO_PARDISO_LIB;

int METIS_NodeND(int *, int *, int *, int *, int *, int *, int *);
int METIS51_NodeND(int *nvtxs, int *xadj, int *adjncy, int *vwgt, int *options, int *perm, int *iperm)
{
	return METIS_NodeND(nvtxs, xadj, adjncy, vwgt, options, perm, iperm);
}

