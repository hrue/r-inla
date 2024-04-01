#include <assert.h>
#include <stdlib.h>
#include <stdio.h>

#if 0
int main(int argc, char **argv)
{
	double x;
#if defined(WINDOWS)
	/*
	 * We cannot read from STDIN in binary format, so therefore we use the first argument 
	 */
	if (argc > 1) {
		FILE *fp = NULL;
		int i;

		for (i = 1; i < argc; i++) {
			fp = fopen(argv[i], "rb");
			assert(fp);
			while (fread(&x, sizeof(double), 1, fp) == 1) {
				printf("%g\n", x);
			}
			fclose(fp);
		}
	} else {
		fprintf(stderr, "Usage: %s BINARYFILE [BINARYFILE2..]\n", argv[0]);
		exit(0);
	}
#else
	while (fread(&x, sizeof(double), 1, stdin) == 1) {
		printf("%g\n", x);
	}
#endif
	return 0;
}
#endif
