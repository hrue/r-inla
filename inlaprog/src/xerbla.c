#include <stdio.h>


void xerbla_(char *srname, int *info)
{
	fprintf(stderr, "xerbla: %s %d\n", srname, *info);
}


void xerbla(char *srname, int *info)
{
	fprintf(stderr, "xerbla: %s %d\n", srname, *info);
}

