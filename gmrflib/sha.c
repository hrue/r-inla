#include <stddef.h>
#include <string.h>
#include <stdint.h>

#include "GMRFLib/sha.h"

#define ROTR(x, n)    (((x) >> (n)) | ((x) << (32 - (n))))
#define Ch(x, y, z)   ((z) ^ ((x) & ((y) ^ (z))))
#define Maj(x, y, z)  (((x) & (y)) ^ ((z) & ((x) ^ (y))))
#define Sigma0(x)     (ROTR(x, 2)  ^ ROTR(x, 13) ^ ROTR(x, 22))
#define Sigma1(x)     (ROTR(x, 6)  ^ ROTR(x, 11) ^ ROTR(x, 25))
#define sigma0(x)     (ROTR(x, 7)  ^ ROTR(x, 18) ^ ((x) >> 3))
#define sigma1(x)     (ROTR(x, 17) ^ ROTR(x, 19) ^ ((x) >> 10))

static inline uint32_t read_be32(const uint8_t *p)
{
	return ((uint32_t) p[0] << 24) | ((uint32_t) p[1] << 16) | ((uint32_t) p[2] << 8) | p[3];
}

static const uint32_t K[64] = {
	0x428a2f98, 0x71374491, 0xb5c0fbcf, 0xe9b5dba5, 0x3956c25b, 0x59f111f1, 0x923f82a4, 0xab1c5ed5,
	0xd807aa98, 0x12835b01, 0x243185be, 0x550c7dc3, 0x72be5d74, 0x80deb1fe, 0x9bdc06a7, 0xc19bf174,
	0xe49b69c1, 0xefbe4786, 0x0fc19dc6, 0x240ca1cc, 0x2de92c6f, 0x4a7484aa, 0x5cb0a9dc, 0x76f988da,
	0x983e5152, 0xa831c66d, 0xb00327c8, 0xbf597fc7, 0xc6e00bf3, 0xd5a79147, 0x06ca6351, 0x14292967,
	0x27b70a85, 0x2e1b2138, 0x4d2c6dfc, 0x53380d13, 0x650a7354, 0x766a0abb, 0x81c2c92e, 0x92722c85,
	0xa2bfe8a1, 0xa81a664b, 0xc24b8b70, 0xc76c51a3, 0xd192e819, 0xd6990624, 0xf40e3585, 0x106aa070,
	0x19a4c116, 0x1e376c08, 0x2748774c, 0x34b0bcb5, 0x391c0cb3, 0x4ed8aa4a, 0x5b9cca4f, 0x682e6ff3,
	0x748f82ee, 0x78a5636f, 0x84c87814, 0x8cc70208, 0x90befffa, 0xa4506ceb, 0xbef9a3f7, 0xc67178f2
};

static void sha256_transform_opt(SHA256_CTX *ctx, const uint8_t *block)
{
	uint32_t a = ctx->state[0], b = ctx->state[1], c = ctx->state[2], d = ctx->state[3];
	uint32_t e = ctx->state[4], f = ctx->state[5], g = ctx->state[6], h = ctx->state[7];
	uint32_t W[64];

	for (int i = 0; i < 16; i++) {
		W[i] = read_be32(block + (i * 4));
	}
	for (int i = 16; i < 64; i++) {
		W[i] = sigma1(W[i - 2]) + W[i - 7] + sigma0(W[i - 15]) + W[i - 16];
	}

	for (int i = 0; i < 64; i++) {
		uint32_t t1 = h + Sigma1(e) + Ch(e, f, g) + K[i] + W[i];
		uint32_t t2 = Sigma0(a) + Maj(a, b, c);
		h = g;
		g = f;
		f = e;
		e = d + t1;
		d = c;
		c = b;
		b = a;
		a = t1 + t2;
	}

	ctx->state[0] += a;
	ctx->state[1] += b;
	ctx->state[2] += c;
	ctx->state[3] += d;
	ctx->state[4] += e;
	ctx->state[5] += f;
	ctx->state[6] += g;
	ctx->state[7] += h;
}

void sha256_init(SHA256_CTX *ctx)
{
	ctx->state[0] = 0x6a09e667;
	ctx->state[1] = 0xbb67ae85;
	ctx->state[2] = 0x3c6ef372;
	ctx->state[3] = 0xa54ff53a;
	ctx->state[4] = 0x510e527f;
	ctx->state[5] = 0x9b05688c;
	ctx->state[6] = 0x1f83d9ab;
	ctx->state[7] = 0x5be0cd19;
	ctx->bitlen = 0;
	ctx->buflen = 0;
}

void sha256_update(SHA256_CTX *ctx, const uint8_t *data, size_t len)
{
	size_t i = 0;

	if (ctx->buflen > 0) {
		size_t left = 64 - ctx->buflen;
		size_t fill = (len < left) ? len : left;
		memcpy(ctx->buffer + ctx->buflen, data, fill);
		ctx->buflen += fill;
		i += fill;
		if (ctx->buflen == 64) {
			sha256_transform_opt(ctx, ctx->buffer);
			ctx->bitlen += 512;
			ctx->buflen = 0;
		}
	}

	while (i + 64 <= len) {
		sha256_transform_opt(ctx, data + i);
		ctx->bitlen += 512;
		i += 64;
	}

	if (i < len) {
		memcpy(ctx->buffer + ctx->buflen, data + i, len - i);
		ctx->buflen += (len - i);
	}
}

void sha256_final(SHA256_CTX *ctx, uint8_t *hash)
{
	uint64_t total_bits = ctx->bitlen + (ctx->buflen * 8);

	ctx->buffer[ctx->buflen++] = 0x80;

	if (ctx->buflen > 56) {
		while (ctx->buflen < 64)
			ctx->buffer[ctx->buflen++] = 0x00;
		sha256_transform_opt(ctx, ctx->buffer);
		ctx->buflen = 0;
	}

	while (ctx->buflen < 56) {
		ctx->buffer[ctx->buflen++] = 0x00;
	}

	for (int i = 0; i < 8; i++) {
		ctx->buffer[56 + i] = (uint8_t) (total_bits >> (56 - i * 8));
	}
	sha256_transform_opt(ctx, ctx->buffer);

	for (int i = 0; i < 8; i++) {
		hash[i * 4] = (uint8_t) (ctx->state[i] >> 24);
		hash[i * 4 + 1] = (uint8_t) (ctx->state[i] >> 16);
		hash[i * 4 + 2] = (uint8_t) (ctx->state[i] >> 8);
		hash[i * 4 + 3] = (uint8_t) (ctx->state[i]);
	}
}


#if 0

/*
 * test program that gives the same result
*/
int main()
{
	const char *text = "hello world";
	uint8_t buf[32];
	SHA256_CTX ctx;

	// Use the functions exactly like OpenSSL legacy setup
	sha256_init(&ctx);
	sha256_update(&ctx, (const uint8_t *) text, strlen(text));
	sha256_final(&ctx, buf);

	// Output the resulting hex digest
	printf("SHA256: ");
	for (int i = 0; i < 32; i++) {
		printf("%02x", buf[i]);
	}
	printf("\n");

	return 0;
}

#       include <openssl/evp.h>
int main()
{
	const char *text = "hello world";
	unsigned char md_value[EVP_MAX_MD_SIZE];
	unsigned int md_len;

	// Allocate a message digest context
	EVP_MD_CTX *mdctx = EVP_MD_CTX_new();
	if (mdctx == NULL) {
		fprintf(stderr, "Failed to create EVP_MD_CTX\n");
		return 1;
	}
	// 1. MODERN INITIALIZATION (Replaces deprecated SHA256_Init)
	if (EVP_DigestInit_ex(mdctx, EVP_sha256(), NULL) != 1) {
		fprintf(stderr, "Digest initialization failed\n");
		EVP_MD_CTX_free(mdctx);
		return 1;
	}
	// 2. MODERN UPDATE (Replaces deprecated SHA256_Update)
	if (EVP_DigestUpdate(mdctx, text, strlen(text)) != 1) {
		fprintf(stderr, "Digest update failed\n");
		EVP_MD_CTX_free(mdctx);
		return 1;
	}
	// 3. MODERN FINAL (Replaces deprecated SHA256_Final)
	if (EVP_DigestFinal_ex(mdctx, md_value, &md_len) != 1) {
		fprintf(stderr, "Digest finalization failed\n");
		EVP_MD_CTX_free(mdctx);
		return 1;
	}
	// Clean up allocated context memory
	EVP_MD_CTX_free(mdctx);

	// Output the resulting hex digest (Should match the optimized standalone version)
	printf("SHA256 (OpenSSL EVP): ");
	for (unsigned int i = 0; i < md_len; i++) {
		printf("%02x", md_value[i]);
	}
	printf("\n");

	return 0;
}
#endif
