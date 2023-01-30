
/* optimize.c
 * 
 * Copyright (C) 2001-2023 Havard Rue
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

/*!
  \file optimize.c
  \brief The optimising routines in GMRFLib.
*/

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>

#include "GMRFLib/GMRFLib.h"
#include "GMRFLib/GMRFLibP.h"

int GMRFLib_default_optimize_param(GMRFLib_optimize_param_tp ** optpar)
{
	/*
	 * define default values for the optimizer. methods are
	 */

	*optpar = Calloc(1, GMRFLib_optimize_param_tp);

	(*optpar)->fp = NULL;
	// (*optpar)->fp = stdout;FIXME("set fp=stdout");
	(*optpar)->opt_type = GMRFLib_OPTTYPE_NR;
	(*optpar)->nr_step_factor = 1.0;
	(*optpar)->nsearch_dir = 1;
	(*optpar)->restart_interval = 10;
	(*optpar)->max_iter = 25;
	(*optpar)->fixed_iter = 0;
	(*optpar)->max_linesearch_iter = 25;
	(*optpar)->step_len = GMRFLib_eps(0.25);
	(*optpar)->stencil = 5;				       /* 3,5,7 */
	(*optpar)->abserr_func = 0.005;
	(*optpar)->abserr_step = 0.005;

	return GMRFLib_SUCCESS;
}
