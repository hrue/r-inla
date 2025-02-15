#include <stddef.h>
#include <string.h>
#include <stdlib.h>

#include "GMRFLib/GMRFLib.h"
#include "GMRFLib/GMRFLibP.h"

int GMRFLib_version(FILE *fp)
{
	fprintf((fp ? fp : stdout), "%s\n", GMRFLib_VERSION);
	return GMRFLib_SUCCESS;
}
