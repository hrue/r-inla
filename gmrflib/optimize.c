#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>

#include "GMRFLib/GMRFLib.h"
#include "GMRFLib/GMRFLibP.h"

int GMRFLib_default_optimize_param(GMRFLib_optimize_param_tp **optpar)
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
	(*optpar)->step_len = GSL_ROOT4_DBL_EPSILON;
	(*optpar)->stencil = 5;				       /* 3,5,7 */
	(*optpar)->abserr_func = 0.005;
	(*optpar)->abserr_step = 0.005;

	return GMRFLib_SUCCESS;
}
