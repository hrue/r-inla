
/* random.c
 * 
 * Copyright (C) 2005-2006 Havard Rue
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
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 *
 * The author's contact information:
 *
 *        Haavard Rue
 *        CEMSE Division
 *        King Abdullah University of Science and Technology
 *        Thuwal 23955-6900, Saudi Arabia
 *        Email: haavard.rue@kaust.edu.sa
 *        Office: +966 (0)12 808 0640
 *
 */

#ifndef HGVERSION
#define HGVERSION
#endif
static const char RCSId[] = "file: " __FILE__ "  " HGVERSION;

/* Pre-hg-Id: $Id: random.c,v 1.28 2008/09/18 07:31:07 hrue Exp $ */

/*!
  \file random.c
  \brief Implement the interface towards the MT19937 random generator in GSL.

  The RNG-routines in GMRFLib are most easily available to the user as the following function pointers:

  - GMRFLib_uniform() used to return a Uniform(0,1) variable.
  - GMRFLib_uniform_init() to initialize the RNG
  - GMRFLib_uniform_getstate() which returns a pointer to a alloc'ed copy of the state of the RNG
  - GMRFLib_uniform_setstate() which set the state in the RNG

  These function points points to the functions \c GMRFLib_rng_uniform(), \c GMRFLib_rng_init(), \c GMRFLib_rng_getstate() and
  \c GMRFLib_rng_setstate() which implements the required functionality using the RNG-routines in GSL. GMRFLib use by default
  the MT19937 random generator in GSL. (This is the same RNG-generator always used by GMRFLib, but from version >= 2.3.0 it use
  the GSL-interface and their implementation.)
  
  The RNG-routines initialise itself by reading a (random) seed from \c /dev/urandom (if it exists), otherwise, it use a fixed
  seed. The behaviour can be overrided by the user, by setting a seed with the function \c GMRFLib_uniform_init().  The global
  variable (of type unsigned long int) \c GMRFLib_rng_seed contains the last seed used to initialise the RNG-routines.

  \sa GMRFLib_uniform, GMRFLib_uniform_init, GMRFLib_uniform_getstate, GMRFLib_uniform_setstate
  
  Example:
  \verbinclude example-doxygen-rng.txt
*/

#include <time.h>
#include <stdio.h>
#include <stddef.h>
#include <string.h>
#if !defined(__FreeBSD__)
#include <malloc.h>
#endif
#include <stdlib.h>
#include <unistd.h>
#include <sys/stat.h>
#include <fcntl.h>

#if defined(WINDOWS)
#include <windows.h>
#include <wincrypt.h>
#endif

#include "GMRFLib/GMRFLib.h"
#include "GMRFLib/GMRFLibP.h"

int GMRFLib_rng_set_default_seed(void)
{
	unsigned long int seed_default = (unsigned long int) time(NULL);
	unsigned long int seed;
	int fd, debug = 0;
	ssize_t nb;
	size_t len = sizeof(unsigned long int);
#pragma omp critical					       /* only one at the time */
	{
#if defined(WINDOWS)
		{
			// this is the eqv of /dev/random for Windows
			HCRYPTPROV prov;
			if (CryptAcquireContext(&prov, NULL, NULL, PROV_RSA_FULL, CRYPT_VERIFYCONTEXT)) {
				if (!CryptGenRandom(prov, (DWORD) len, (BYTE *) & seed)) {
					// error: fall back to default
					seed = seed_default;
				} else {
					CryptReleaseContext(prov, 0);
				}
			} else {
				// error: fall back to default
				seed = seed_default;
			}
		}
#else							       /* !defined(WINDOWS) */
		{
			fd = open("/dev/urandom", O_RDONLY);
			if (fd > 0) {
				nb = read(fd, (void *) &seed, len);
				if (nb != (ssize_t) len) {
					seed = seed_default;
				}
				close(fd);
			} else {
				// error: fall back to default
				seed = seed_default;
			}
		}
#endif							       /* defined(WINDOWS) */
	}

	if (debug)
		fprintf(stderr, "Init RNG with seed %lu\n", seed);

	GMRFLib_rng_init(seed);

	return GMRFLib_SUCCESS;
}
int GMRFLib_rng_init(unsigned long int seed)
{
	if (GMRFLib_rng_ptr) {
		gsl_rng_free(GMRFLib_rng_ptr);
	}
	GMRFLib_rng_ptr = gsl_rng_alloc(gsl_rng_mt19937);
	gsl_rng_set(GMRFLib_rng_ptr, seed);		       /* or GMRFLib_rng .... */
	GMRFLib_rng_seed = seed;			       /* keep a copy */

	return GMRFLib_SUCCESS;
}

gsl_rng *GMRFLib_rng_init_default(void)
{
	GMRFLib_rng_set_default_seed();
	return GMRFLib_rng_ptr;
}
double GMRFLib_rng_uniform(void)
{
	return gsl_rng_uniform_pos(GMRFLib_rng);
}
void *GMRFLib_rng_getstate(size_t * siz)
{
	size_t n;
	void *p, *pp;

	p = gsl_rng_state(GMRFLib_rng);
	n = gsl_rng_size(GMRFLib_rng);
	pp = Calloc(n, char);

	memcpy(pp, p, n);
	if (siz) {
		*siz = n;
	}

	return pp;
}
int GMRFLib_rng_setstate(void *saved_state)
{
	if (saved_state) {
		void *p;
		size_t n;

		p = gsl_rng_state(GMRFLib_rng);
		n = gsl_rng_size(GMRFLib_rng);
		memcpy(p, saved_state, n);
	}
	return GMRFLib_SUCCESS;
}
