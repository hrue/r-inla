#ifndef __GMRFLib_SHA_H__
#       define __GMRFLib_SHA_H__

#       undef __BEGIN_DECLS
#       undef __END_DECLS
#       ifdef __cplusplus
#              define __BEGIN_DECLS extern "C" {
#              define __END_DECLS }
#       else
#              define __BEGIN_DECLS			       /* empty */
#              define __END_DECLS			       /* empty */
#       endif

#       include "GMRFLib/sha.h"

__BEGIN_DECLS
//
    typedef struct {
	uint32_t state[8];
	uint64_t bitlen;
	uint8_t buffer[64];
	uint32_t buflen;
} SHA256_CTX;

void sha256_transform(SHA256_CTX * ctx, const uint8_t data[]);
void sha256_init(SHA256_CTX * ctx);
void sha256_update(SHA256_CTX * ctx, const uint8_t data[], size_t len);
void sha256_final(SHA256_CTX * ctx, uint8_t hash[]);


#       define GMRFLib_SHA_TP         SHA256_CTX
#       define GMRFLib_SHA_DIGEST_LEN 32L
#       define GMRFLib_SHA_Init       sha256_init
#       define GMRFLib_SHA_Update     sha256_update
#       define GMRFLib_SHA_Final      sha256_final
#       define GMRFLib_SHA_UPDATE_LEN 64L
#       define GMRFLib_SHA_UPDATE_CORE(_x, _len, _type, _c) \
	if ((_len) > 0) {						\
		size_t len_ = (_len) * sizeof(_type);			\
		size_t n_ = (size_t) len_ / GMRFLib_SHA_UPDATE_LEN;	\
		size_t m_ = len_ - (n_) * GMRFLib_SHA_UPDATE_LEN;	\
		unsigned char *xx = (unsigned char *) (_x);		\
		for(size_t i_ = 0; i_ < (n_); i_++) {			\
			GMRFLib_SHA_Update(&(_c), (const void *) (xx + i_ * GMRFLib_SHA_UPDATE_LEN), (size_t) GMRFLib_SHA_UPDATE_LEN); \
		}							\
		if (m_) {						\
			GMRFLib_SHA_Update(&(_c), (const void *) (xx + (n_) * GMRFLib_SHA_UPDATE_LEN), m_); \
		}							\
	}
#       define GMRFLib_SHA_IUPDATE(_x, _len, _c) GMRFLib_SHA_UPDATE_CORE(_x, _len, int, _c)
#       define GMRFLib_SHA_DUPDATE(_x, _len, _c) GMRFLib_SHA_UPDATE_CORE(_x, _len, double, _c)
__END_DECLS
#endif
