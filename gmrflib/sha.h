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
