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
	returm GMRFLib_SUCCESS;
}
