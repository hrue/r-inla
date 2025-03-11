
/*********************************************************/

/* TAUCS                                                 */

/* Author: Sivan Toledo                                  */

/*                                                       */

/* Simple complex arithmetic routines.                   */

/* They are called if the compiler does not support      */

/* complex. GCC supports complex, and so do all C99      */

/* compilers.                                            */

/*                                                       */

/*********************************************************/

#include <math.h>
#include "taucs.h"

#ifdef TAUCS_CORE_DOUBLE
double taucs_get_nan()
{
	double zero = 0.0;
	double inf = 1.0 / zero;
	double nan = inf - inf;

	return nan;
}
#endif
