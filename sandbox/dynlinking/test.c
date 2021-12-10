#include <stdio.h>
#include <stdlib.h>
#include <dlfcn.h>
int main(int argc, char **argv) {
    void *handle;
    void *handle2;
    double (*cosine)(double);
    double (*sine)(double);
    double (*sine2)(double);
    char *error;
    handle = dlopen ("libm.so.6", RTLD_NOW);
    if (!handle) {
        fprintf (stderr, "%s\n", dlerror());
        exit(1);
    }
    dlerror();    /* Clear any existing error */

    handle2 = dlopen ("libm.so.6", RTLD_NOW);
    if (!handle2) {
        fprintf (stderr, "%s\n", dlerror());
        exit(1);
    }
    dlerror();    /* Clear any existing error */

    cosine = dlsym(handle, "cos");
    if ((error = dlerror()) != NULL)  {
        fprintf (stderr, "%s\n", error);
        exit(1);
    }

    sine = dlsym(handle2, "sin");
    if ((error = dlerror()) != NULL)  {
        fprintf (stderr, "%s\n", error);
        exit(1);
    }

    sine2 = dlsym(handle2, "sin");
    if ((error = dlerror()) != NULL)  {
        fprintf (stderr, "%s\n", error);
        exit(1);
    }

    printf ("%f\n", (*cosine)(2.0));
    printf ("%f\n", (*sine)(2.0));
    printf ("%f\n", (*sine2)(2.0));

    dlclose(handle);
    dlclose(handle2);

    return 0;
}
