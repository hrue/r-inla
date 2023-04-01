
/* random.c
 * 
 * Copyright (C) 2005-2023 Havard Rue
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

#include <time.h>
#include <stdio.h>
#include <stddef.h>
#include <string.h>
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

static unsigned long int GMRFLib_rng_seed;
#pragma omp threadprivate(GMRFLib_rng_seed)

int GMRFLib_rng_set_default_seed(void)
{
	unsigned long int seed_default = (unsigned long int) time(NULL);
	unsigned long int seed;
	const int debug = 0;
	size_t len = sizeof(unsigned long int);
#pragma omp critical (Name_96da5f632ecbd97ae1e5504794f8724fabfdee73)
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
			int fd = open("/dev/urandom", O_RDONLY);
			if (fd > 0) {
				ssize_t nb = read(fd, (void *) &seed, len);
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
		fprintf(stderr, "Init RNG with seed %zu\n", (size_t) seed);

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

void *GMRFLib_rng_getstate(size_t *siz)
{
	size_t n;
	void *p, *pp;

	p = gsl_rng_state(GMRFLib_rng);
	n = gsl_rng_size(GMRFLib_rng);
	pp = Calloc(n, char);

	Memcpy(pp, p, n);
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
		Memcpy(p, saved_state, n);
	}
	return GMRFLib_SUCCESS;
}
