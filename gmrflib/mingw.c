
/* mingw.c
 * 
 * Copyright (C) 2024-2024 Havard Rue
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

// workaround for mingw32 builds when linking with MKL


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
