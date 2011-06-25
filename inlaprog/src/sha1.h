// I found this somewhere and I cannot find it again... It was under GPL, and I have modified the code to fit with inla.c
// Feb 11, 2009, Havard Rue

#ifndef __INLA_SHA1_H__
#define __INLA_SHA1_H__
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
#include <stdio.h>
    typedef struct {
	unsigned char data[64];
	unsigned int datalen;
	unsigned int bitlen[2];
	unsigned int state[5];
	unsigned int k[4];
} SHA_CTX;

#if !defined(SHA_DIGEST_LENGTH)
#define SHA_DIGEST_LENGTH 20
#endif

void SHA1_Transform(SHA_CTX * ctx, unsigned char data[]);
void SHA1_Init(SHA_CTX * ctx);
void SHA1_Update(SHA_CTX * ctx, unsigned char data[], unsigned long len);
void SHA1_Final(unsigned char hash[], SHA_CTX * ctx);


__END_DECLS
#endif
