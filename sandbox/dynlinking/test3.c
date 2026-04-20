/// THIS DOES NOT WORK

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <dlfcn.h>

double real_function(double x) {
	return sin(x);
}


int main(int argc, char **argv) {
    void *handle;
    double (*fun)(double);
    char *error;

    handle = dlopen ("libelias.so", RTLD_NOW);
    if (!handle) {
        fprintf (stderr, "%s\n", dlerror());
        exit(1);
    }
    dlerror();    /* Clear any existing error */

    fun = dlsym(handle, "eliasfunc");
    if ((error = dlerror()) != NULL)  {
        fprintf (stderr, "%s\n", error);
        exit(1);
    }

    printf ("%f\n", (*fun)(1.0));
    printf ("%f\n", (*fun)(2.0));
    printf ("%f\n", (*fun)(3.0));

    dlclose(handle);

    return 0;
}
