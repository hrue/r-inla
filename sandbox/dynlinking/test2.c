#include <stdio.h>
#include <stdlib.h>
#include <ltdl.h>
int main(int argc, char **argv) {
	lt_dlhandle handle;
	lt_dlhandle handle2;
	double (*cosine)(double);
	double (*sine)(double);
	double (*sine2)(double);
	const char *error;
	
	lt_dlinit();
	
	handle = lt_dlopen ("libm.so.6");
	if (!handle) {
		fprintf (stderr, "%s\n", lt_dlerror());
		exit(1);
	}
	lt_dlerror();    /* Clear any existing error */

	handle2 = lt_dlopen("libm.so.6");
	if (!handle2) {
		fprintf (stderr, "%s\n", lt_dlerror());
		exit(1);
	}
	lt_dlerror();    /* Clear any existing error */

	cosine = lt_dlsym(handle, "cos");
	if ((error = lt_dlerror()) != NULL)  {
		fprintf (stderr, "%s\n", error);
		exit(1);
	}

	sine = lt_dlsym(handle2, "sin");
	if ((error = lt_dlerror()) != NULL)  {
		fprintf (stderr, "%s\n", error);
		exit(1);
	}

	sine2 = lt_dlsym(handle2, "sin");
	if ((error = lt_dlerror()) != NULL)  {
		fprintf (stderr, "%s\n", error);
		exit(1);
	}

	printf ("%f\n", (*cosine)(2.0));
	printf ("%f\n", (*sine)(2.0));
	printf ("%f\n", (*sine2)(2.0));

	lt_dlclose(handle);
	lt_dlclose(handle2);

	return 0;
}
