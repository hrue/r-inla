#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

#define MSG(fun_) fprintf(stderr, "\n\n*** Warning *** Fake function [%s] called. This may, or may not, go well.\n\n\n", fun_)
#define mkfun(fun_) void fun_(void) { static char first = 1;  if (first) MSG(# fun_); first = 0; return; }

#if defined(WINDOWS) && defined(INLA_WITH_MKL) && (defined(__MINGW32__) || defined(__MINGW64__))

mkfun(__GSHandlerCheck);
mkfun(__chkstk);

uint64_t __security_cookie;
void __security_init_cookie()
{
	__security_cookie = 0;
}
void __security_check_cookie(uint64_t retrieved)
{
	if (__security_cookie != retrieved) {
		// abort();
	}
}

#endif

#if defined(WINDOWS) && defined(INLA_WITH_OPENBLAS) && (defined(__MINGW32__) || defined(__MINGW64__))
mkfun(__imp__cprintf);
#endif
