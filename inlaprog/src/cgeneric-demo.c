
/* cgeneric-demo.c
 * 
 * Copyright (C) 2021-2024 Havard Rue
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


#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <strings.h>

#include "cgeneric.h"

#define Calloc(n_, type_)  (type_ *)calloc((n_), sizeof(type_))
#define SQR(x) ((x)*(x))

double *inla_cgeneric_iid_model(inla_cgeneric_cmd_tp cmd, double *theta, inla_cgeneric_data_tp *data)
{
	// this reimplement `inla.rgeneric.iid.model` using cgeneric

	double *ret = NULL, prec = (theta ? exp(theta[0]) : NAN), lprec = (theta ? theta[0] : NAN);

	assert(!strcasecmp(data->ints[0]->name, "n"));	       // this will always be the case
	int N = data->ints[0]->ints[0];			       // this will always be the case
	assert(N > 0);

	switch (cmd) {
	case INLA_CGENERIC_VOID:
	{
		assert(!(cmd == INLA_CGENERIC_VOID));
	}
		break;

	case INLA_CGENERIC_GRAPH:
	{
		// return a vector of indices with format
		// c(N, M, ii, jj)
		// where ii<=jj, ii is non-decreasing and jj is non-decreasing for the same ii
		// so like the loop
		// for i=0, ...
		// for j=i, ...
		// G_ij = 
		// and M is the total length while N is the dimension

		int M = N;
		ret = Calloc(2 + 2 * N, double);
		assert(ret);
		ret[0] = N;				       /* dimension */
		ret[1] = M;				       /* number of (i <= j) */
		for (int i = 0; i < M; i++) {
			ret[2 + i] = i;			       /* i */
			ret[2 + N + i] = i;		       /* j */
		}
	}
		break;

	case INLA_CGENERIC_Q:
	{
		// return c(-1, M, Qij) in the same order as defined in INLA_CGENERIC_GRAPH
		int M = N;
		ret = Calloc(2 + N, double);
		assert(ret);
		ret[0] = -1;				       /* REQUIRED! */
		ret[1] = M;				       /* number of (i <= j) */
		for (int i = 0; i < M; i++) {
			ret[2 + i] = prec;
		}
	}
		break;

	case INLA_CGENERIC_MU:
	{
		// return (N, mu)
		// if N==0 then mu is not needed as its taken to be mu[]==0
		ret = Calloc(1, double);
		assert(ret);
		ret[0] = 0;
	}
		break;

	case INLA_CGENERIC_INITIAL:
	{
		// return c(M, initials)
		// where M is the number of hyperparameters
		ret = Calloc(2, double);
		assert(ret);
		ret[0] = 1;
		ret[1] = 4.0;
	}
		break;

	case INLA_CGENERIC_LOG_NORM_CONST:
	{
		// return c(NORM_CONST) or a NULL-pointer if INLA should compute it by itself
		ret = Calloc(1, double);
		assert(ret);
		ret[0] = N * (-0.9189385332 + 0.5 * lprec);
	}
		break;

	case INLA_CGENERIC_LOG_PRIOR:
	{
		// return c(LOG_PRIOR)
		ret = Calloc(1, double);
		assert(ret);
		ret[0] = -prec + lprec;			       // prec ~ gamma(1,1)
	}
		break;

	case INLA_CGENERIC_QUIT:
	default:
		break;
	}

	return (ret);
}

double *inla_cgeneric_ar1_model(inla_cgeneric_cmd_tp cmd, double *theta, inla_cgeneric_data_tp *data)
{
	// this reimplement `inla.rgeneric.ar1.model` using cgeneric

	double *ret = NULL, prec, lprec, rho, rho_intern;

	if (theta) {
		lprec = theta[0];
		prec = exp(lprec);
		rho_intern = theta[1];
		rho = 2.0 * exp(rho_intern) / (1.0 + exp(rho_intern)) - 1.0;
	} else {
		prec = lprec = rho = rho_intern = NAN;
	}

	assert(!strcasecmp(data->ints[0]->name, "n"));	       // this will always be the case
	int N = data->ints[0]->ints[0];			       // this will always be the case
	assert(N > 0);

	switch (cmd) {
	case INLA_CGENERIC_VOID:
	{
		assert(!(cmd == INLA_CGENERIC_VOID));
	}
		break;

	case INLA_CGENERIC_GRAPH:
	{
		// return a vector of indices with format
		// c(N, M, ii, jj)
		// where ii<=jj, ii is non-decreasing and jj is non-decreasing for the same ii
		// so like the loop
		// for i=0, ...
		// for j=i, ...
		// G_ij = 
		// and M is the length of ii

		int M = N + N - 1, offset, i, k;
		ret = Calloc(2 + 2 * M, double);
		assert(ret);
		offset = 2;
		ret[0] = N;				       /* dimension */
		ret[1] = M;				       /* number of (i <= j) */
		for (k = i = 0; i < N; i++) {
			ret[offset + k] = i;		       /* i */
			ret[offset + M + k++] = i;	       /* j */
			if (i < N - 1) {
				ret[offset + k] = i;	       /* i */
				ret[offset + M + k++] = i + 1; /* j */
			}
		}
	}
		break;

	case INLA_CGENERIC_Q:
	{
		// optimized format
		// return c(-1, M, Qij) in the same order as defined in INLA_CGENERIC_GRAPH
		// where M is the length of Qij

		double param = prec / (1.0 - SQR(rho));
		int M = N + N - 1;
		int offset, i, k;
		ret = Calloc(2 + M, double);
		assert(ret);
		offset = 2;
		ret[0] = -1;				       /* REQUIRED */
		ret[1] = M;
		for (i = k = 0; i < N; i++) {
			ret[offset + k++] = param * (i == 0 || i == N - 1 ? 1.0 : (1.0 + SQR(rho)));
			if (i < N - 1) {
				ret[offset + k++] = -param * rho;
			}
		}
	}
		break;

	case INLA_CGENERIC_MU:
	{
		// return (N, mu)
		// if N==0 then mu is not needed as its taken to be mu[]==0

		ret = Calloc(1, double);
		assert(ret);
		ret[0] = 0;
	}
		break;

	case INLA_CGENERIC_INITIAL:
	{
		// return c(M, initials)
		// where M is the number of hyperparameters

		ret = Calloc(3, double);
		assert(ret);
		ret[0] = 2;
		ret[1] = 1.0;
		ret[2] = 1.0;
	}
		break;

	case INLA_CGENERIC_LOG_NORM_CONST:
	{
		// return c(NORM_CONST) or a NULL-pointer if INLA should compute it by itself

		double prec_innovation = prec / (1.0 - SQR(rho));
		ret = Calloc(1, double);
		assert(ret);
		ret[0] = N * (-0.5 * log(2.0 * M_PI) + 0.5 * log(prec_innovation)) + 0.5 * log(1.0 - SQR(rho));
	}
		break;

	case INLA_CGENERIC_LOG_PRIOR:
	{
		// return c(LOG_PRIOR)

		ret = Calloc(1, double);
		assert(ret);
		ret[0] = -prec + lprec - 0.5 * log(2.0 * M_PI) - 0.5 * SQR(rho_intern);
	}
		break;

	case INLA_CGENERIC_QUIT:
	default:
		break;
	}

	return (ret);
}

double *inla_cgeneric_generic0_model(inla_cgeneric_cmd_tp cmd, double *theta, inla_cgeneric_data_tp *data)
{
	double *ret = NULL, prec, lprec;

	if (theta) {
		lprec = theta[0];
		prec = exp(lprec);
	} else {
		prec = lprec = NAN;
	}

	assert(data->n_ints > 0);
	assert(!strcasecmp(data->ints[0]->name, "n"));	       // this will always be the case
	int N = data->ints[0]->ints[0];			       // this will always be the case
	assert(N > 0);

	/*
	 * we assume Cmatrix is stored column-wise, like
	 * for(i...)
	 *     for(j=i...)
	 *         Q_ij
	 * 
	 * so, upper half only
	 */
	assert(data->n_smats > 0);
	assert(!strcasecmp(data->smats[0]->name, "Cmatrix"));
	inla_cgeneric_smat_tp *Cmatrix = data->smats[0];
	assert(N == Cmatrix->nrow);
	assert(N == Cmatrix->ncol);

	switch (cmd) {
	case INLA_CGENERIC_VOID:
	{
		assert(!(cmd == INLA_CGENERIC_VOID));
	}
		break;

	case INLA_CGENERIC_GRAPH:
	{
		// return a vector of indices with format
		// c(N, M, ii, jj)
		// where ii<=jj, ii is non-decreasing and jj is non-decreasing for the same ii
		// so like the loop
		// for i=0, ...
		// for j=i, ...
		// G_ij = 
		// and M is the length of ii

		int M = Cmatrix->n, offset;
		ret = Calloc(2 + 2 * M, double);
		assert(ret);
		offset = 2;
		ret[0] = N;				       /* dimension */
		ret[1] = M;				       /* number of (i <= j) */
		for (int k = 0; k < M; k++) {
			ret[offset + k] = Cmatrix->i[k];       /* i */
			ret[offset + M + k] = Cmatrix->j[k];   /* j */
		}
	}
		break;

	case INLA_CGENERIC_Q:
	{
		// optimized format
		// return c(-1, M, Qij) in the same order as defined in INLA_CGENERIC_GRAPH
		// where M is the length of Qij

		int M = Cmatrix->n;
		int offset;
		ret = Calloc(2 + M, double);
		assert(ret);
		offset = 2;
		ret[0] = -1;				       /* REQUIRED */
		ret[1] = M;
		for (int k = 0; k < M; k++) {
			ret[offset + k] = prec * Cmatrix->x[k];
		}
	}
		break;

	case INLA_CGENERIC_MU:
	{
		// return (N, mu)
		// if N==0 then mu is not needed as its taken to be mu[]==0

		ret = Calloc(1, double);
		assert(ret);
		ret[0] = 0;
	}
		break;

	case INLA_CGENERIC_INITIAL:
	{
		// return c(M, initials)
		// where M is the number of hyperparameters

		ret = Calloc(2, double);
		assert(ret);
		ret[0] = 1;
		ret[1] = 4.0;
	}
		break;

	case INLA_CGENERIC_LOG_NORM_CONST:
	{
		// return c(NORM_CONST) or a NULL-pointer if INLA should compute it by itself. here we ignore the part that comes from the
		// 1/2*log(det(Q)), which could be added if needed later. this is how the 'generic0' model is implemented.
		ret = Calloc(1, double);
		assert(ret);
		ret[0] = N / 2.0 * (lprec - log(2.0 * M_PI));
	}
		break;

	case INLA_CGENERIC_LOG_PRIOR:
	{
		// return c(LOG_PRIOR). with a Gamma(1,1) for precision, this is the log prior for the log(precision).
		ret = Calloc(1, double);
		assert(ret);
		ret[0] = -prec + lprec;
	}
		break;

	case INLA_CGENERIC_QUIT:
	default:
		break;
	}

	return (ret);
}
