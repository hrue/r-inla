
/* example-rng.c
 * 
 * Copyright (C) 2006 Havard Rue
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * The author's contact information:
 *
 *       H{\aa}vard Rue
 *       Department of Mathematical Sciences
 *       The Norwegian University of Science and Technology
 *       N-7491 Trondheim, Norway
 *       Voice: +47-7359-3533    URL  : http://www.math.ntnu.no/~hrue  
 *       Fax  : +47-7359-3524    Email: havard.rue@math.ntnu.no
 *
 */

#if !defined(__FreeBSD__)
#include <malloc.h>
#endif
#include <stdlib.h>
#include <stddef.h>
#include <stdio.h>

#include "GMRFLib/GMRFLib.h"

static const char RCSId[] = "$Id: example-rng.c,v 1.12 2007/04/02 19:53:20 hrue Exp $";

int main(int argc, char **argv)
{
	/*
	 * A simple illustration of the get/set-functiononality for the RNG, plus the use of GSL-routines requiring a RNG.
	 * 
	 * The seed is either given as the first argument or set by GMRFLib otherwise. 
	 */

	int i;
	unsigned long int seed;
	void *state;

	/*
	 * both options are fine... 
	 */
	if (1) {
		if (argc > 1) {
			seed = (unsigned long int) atol(argv[1]);
			printf("Use seed from the user: %lu\n", seed);
			GMRFLib_uniform_init(seed);
		} else {
			GMRFLib_uniform();		       /* force the seed to be set */
			seed = GMRFLib_rng_seed;
			printf("Use seed provided from GMRFLib: %lu\n", GMRFLib_rng_seed);
		}
	} else {
		printf("No initialisation....\n");
	}

	printf("first five U(0,1)'s\n");
	for (i = 0; i < 5; i++) {
		printf("\ti=%d Unif(0,1)=%f\n", i, GMRFLib_uniform());
	}

	printf("get state\n");
	state = GMRFLib_uniform_getstate();

	printf("...continue\n");
	for (i = 0; i < 5; i++) {
		printf("\ti=%d Unif(0,1)=%f\n", i, GMRFLib_uniform());
	}

	printf("change state: use the copy\n");
	GMRFLib_uniform_setstate(state);

	printf("...restart using the saved state\n");
	for (i = 0; i < 5; i++) {
		printf("\ti=%d Unif(0,1)=%f\n", i, GMRFLib_uniform());
	}

	printf("cleanup\n");
	free(state);

	/*
	 * Use the global pointer `GMRFLib_rng' to access the RNG-object used by GMRFLib, when calling GSL-routines requiring a
	 * `gsl_rng'-object. It is not needed to initialise the RNG-functions before GMRFLib_rng is used.
	 */

	printf("demonstrate the use of GSL routines using the GMRFLib's RNG-interface\n");
	printf("\tproduce some Poisson random variables using: gsl_ran_poisson(GMRFLib_rng, 1.0))\n");
	for (i = 0; i < 5; i++) {
		printf("\ti=%d Poisson(1)=%d\n", i, (int) gsl_ran_poisson(GMRFLib_rng, 1.0));
	}

	printf("get state\n");
	state = GMRFLib_uniform_getstate();

	printf("\tcontinue to produce some Poisson random variables\n");
	for (i = 0; i < 5; i++) {
		printf("\ti=%d Poisson(1)=%d\n", i, (int) gsl_ran_poisson(GMRFLib_rng, 1.0));
	}

	printf("change state: use the copy\n");
	GMRFLib_uniform_setstate(state);

	printf("generate some Poisson random variables restarting from the saved state\n");
	for (i = 0; i < 5; i++) {
		printf("\ti=%d Poisson(1)=%d\n", i, (int) gsl_ran_poisson(GMRFLib_rng, 1.0));
	}

	printf("cleanup\n");
	free(state);

	return 0;
}
