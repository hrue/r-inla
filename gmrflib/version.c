#include <omp.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "GMRFLib/GMRFLib.h"

int GMRFLib_version(FILE *fp)
{
	fprintf((fp ? fp : stdout), "%s\n", GMRFLib_VERSION);
	return GMRFLib_SUCCESS;
}
