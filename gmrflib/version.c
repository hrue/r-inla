#include <stddef.h>
#include <string.h>
#include <stdlib.h>

#include "GMRFLib/GMRFLib.h"
#include "GMRFLib/GMRFLibP.h"

int GMRFLib_version(FILE *fp)
{
	const char *GitIDs[] = {
#include "GitID.all"
		NULL
	};
	int i = 0;

	fp = (fp ? fp : stdout);

	fprintf(fp, "GMRFLib %s\n", GMRFLib_VERSION);
	while (GitIDs[i]) {
		fprintf(fp, "\t%s\n", GitIDs[i]);
		i++;
	}

	return GMRFLib_SUCCESS;
}
