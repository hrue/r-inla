#include <stddef.h>
#include <assert.h>
#include <strings.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#include "GMRFLib/GMRFLib.h"
#include "GMRFLib/GMRFLibP.h"

static GMRFLib_stiles_ctl_tp  stiles_ctl = {0, 0};

int GMRFLib_stiles_set_param(int verbose, int debug) 
{
	if (verbose >= 0) {
		stiles_ctl.verbose = verbose;
	}
	if (debug >= 0) {
		stiles_ctl.debug = debug;
	}
	//printf("\n\nstiles set verbose=%1d debug=%1d\n\n", verbose, debug);
	return GMRFLib_SUCCESS;
}

GMRFLib_ptr_tp *GMRFLib_stiles_get_graphs(void *mb) 
{
	// need to call it in the inla_ parts as need need access to inla_tp, which
	// is a mess to make available here
	GMRFLib_ptr_tp *inla_stiles_get_graphs(void *);
	return inla_stiles_get_graphs(mb);
}
