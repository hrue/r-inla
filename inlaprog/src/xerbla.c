#include "GMRFLib/GMRFLib.h"

void xerbla_(char *srname, int *info, FORTRAN_CHARLEN_T len)
{
	// fortran version
	char *name = (char *) calloc(len + 1, sizeof(char));
	Assert(name);
	memcpy(name, srname, len * sizeof(char));
	name[len] = '\0';
	LOG_ERROR("name = %s info = %d", name, *info);
	exit(1);
}

void xerbla(char POSSIBLY_UNUSED(*srname), int POSSIBLY_UNUSED(*info))
{
	abort();
}
