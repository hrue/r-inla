#include <ctype.h>
#include <stdio.h>
#include <string.h>
#include "my-fix.h"
#include "iniparser.h"
#include "GMRFLib/GMRFLib.h"
#include "GMRFLib/GMRFLibP.h"

int my_is_double(char *str)
{
	/*
	 * return 1 if a double can be read from STR and 0 if not. 
	 */
	double x;

	return (sscanf(str, "%lf", &x) == 1 ? 1 : 0);
}

int my_is_int(char *str)
{
	/*
	 * return 1 if an int can be read from STR and 0 if not. 
	 */

	int x;

	return (sscanf(str, "%d", &x) == 1 ? 1 : 0);
}
char *my_strlwc(const char *s)
{
	/*
	 * str = A:B only lowercase part B. 
	 */
	long int i, start, debug = 0;
	char *f, *str;

	str = GMRFLib_strdup(s);
	if (debug)
		printf("str in %s\n", str);
	f = strchr(str, INIPARSER_SEP);
	start = (f ? (f - str) : 0);
	for (i = start; i < (long int) strlen(str); i++) {
		str[i] = (char) tolower((int) str[i]);
	}
	if (debug)
		printf("str out %s\n", str);

	return str;
}

#if defined(WINDOWS)
// provide just this functiononality. this is required for windows compilation, I think (without posix libs)
double drand48(void)
{
	return GMRFLib_uniform();
}
void srand48(long int seed)
{
	return;
}
#endif

#if defined(INLA_WINDOWS32_FIX)
void _mm_pause(void) {
  return;
}
#endif
