
/* GMRFLib-smtp-profile.c
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
 *        H{a}vard Rue
 *        CEMSE Division
 *        King Abdullah University of Science and Technology
 *        Thuwal 23955-6900, Saudi Arabia
 *        Email: haavard.rue@kaust.edu.sa
 *        Office: +966 (0)12 808 0640
 *
 */

/*!
  \file smtp-profile.c
  \brief Template for an interface to a new sparse-matrix library
*/

#include "GMRFLib/GMRFLib.h"
#include "GMRFLib/GMRFLibP.h"

#ifndef HGVERSION
#define HGVERSION
#endif
static const char RCSId[] = "file: " __FILE__ "  " HGVERSION;

/* Pre-hg-Id: $Id: smtp-profile.c,v 1.10 2010/02/26 17:54:28 hrue Exp $ */

int GMRFLib_compute_reordering_PROFILE(void)
{
	abort();
	return GMRFLib_SUCCESS;
}
int GMRFLib_build_sparse_matrix_PROFILE(void)
{
	abort();
	return GMRFLib_SUCCESS;
}
int GMRFLib_factorise_sparse_matrix_PROFILE(void)
{
	abort();
	return GMRFLib_SUCCESS;
}
int GMRFLib_free_fact_sparse_matrix_PROFILE(void)
{
	abort();
	return GMRFLib_SUCCESS;
}
int GMRFLib_solve_lt_sparse_matrix_PROFILE(void)
{
	abort();
	return GMRFLib_SUCCESS;
}
int GMRFLib_solve_llt_sparse_matrix_PROFILE(void)
{
	abort();
	return GMRFLib_SUCCESS;
}
int GMRFLib_solve_l_sparse_matrix_PROFILE(void)
{
	abort();
	return GMRFLib_SUCCESS;
}
int GMRFLib_solve_lt_sparse_matrix_special_PROFILE(void)
{
	abort();
	return GMRFLib_SUCCESS;
}
int GMRFLib_solve_l_sparse_matrix_special_PROFILE(void)
{
	abort();
	return GMRFLib_SUCCESS;
}
int GMRFLib_comp_cond_meansd_PROFILE(void)
{
	abort();
	return GMRFLib_SUCCESS;
}
int GMRFLib_log_determinant_PROFILE(void)
{
	abort();
	return GMRFLib_SUCCESS;
}
int GMRFLib_compute_Qinv_PROFILE(void)
{
	abort();
	return GMRFLib_SUCCESS;
}

int GMRFLib_bitmap_factorisation_PROFILE(void)
{
	abort();
	return GMRFLib_SUCCESS;
}
