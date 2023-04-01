
/* sha.h
 * 
 * Copyright (C) 2022-2023 Havard Rue
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
 *
 */

#ifndef __GMRFLib_SHA_H__
#define __GMRFLib_SHA_H__

/* Don't warn about deprecated functions, */
#ifndef OPENSSL_API_COMPAT
  // 0x10101000L == 1.1.1, 30000 == 3.0.0
#define OPENSSL_API_COMPAT 0x10101000L
#endif
#include <openssl/sha.h>

#undef __BEGIN_DECLS
#undef __END_DECLS
#ifdef __cplusplus
#define __BEGIN_DECLS extern "C" {
#define __END_DECLS }
#else
#define __BEGIN_DECLS					       /* empty */
#define __END_DECLS					       /* empty */
#endif

__BEGIN_DECLS
//
#define GMRFLib_SHA_TP         SHA256_CTX
#define GMRFLib_SHA_DIGEST_LEN SHA256_DIGEST_LENGTH
#define GMRFLib_SHA_Init       SHA256_Init
#define GMRFLib_SHA_Update     SHA256_Update
#define GMRFLib_SHA_Final      SHA256_Final
#define GMRFLib_SHA_UPDATE_LEN 64L
#define GMRFLib_SHA_UPDATE_CORE(_x, _len, _type) \
	if ((_len) > 0 && (_x)) {					\
		size_t len = (_len) * sizeof(_type);			\
		size_t n = (size_t) len / GMRFLib_SHA_UPDATE_LEN;	\
		size_t m = len - n * GMRFLib_SHA_UPDATE_LEN;		\
		unsigned char *xx = (unsigned char *) (_x);		\
		for(size_t i_ = 0; i_ < n; i_++) {			\
			GMRFLib_SHA_Update(&c, (const void *) (xx + i_ * GMRFLib_SHA_UPDATE_LEN), (size_t) GMRFLib_SHA_UPDATE_LEN); \
		}							\
		if (m) {						\
			GMRFLib_SHA_Update(&c, (const void *) (xx + n * GMRFLib_SHA_UPDATE_LEN), m); \
		}							\
	}
#define GMRFLib_SHA_IUPDATE(_x, _len) GMRFLib_SHA_UPDATE_CORE(_x, _len, int)
#define GMRFLib_SHA_DUPDATE(_x, _len) GMRFLib_SHA_UPDATE_CORE(_x, _len, double)
    __END_DECLS
#endif
