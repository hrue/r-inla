#include <string.h>
#include <stdlib.h>
#include <stdio.h>

// need fortran_charlen_t
#include "GMRFLib/GMRFLibP.h"

void xerbla_(char *srname, int *info, fortran_charlen_t len)
{
	// fortran version
	char *name = (char *) calloc(len + 1, sizeof(char));
	assert(name);
	memcpy(name, srname, len * sizeof(char));
	name[len] = '\0';
	fprintf(stderr, "\nxerbla: %s %d\n", name, *info);
	exit(1);
}

void xerbla(char *srname, int *info)
{
	fprintf(stderr, "\nxerbla: %s %d\n", srname, *info);
	exit(1);
}
