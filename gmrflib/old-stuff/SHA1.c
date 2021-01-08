


#if 0
// I found this somewhere and I cannot find it again... It was under GPL, and I have modified the code to fit with inla.c.
// can also use openssl

typedef struct {
	unsigned char data[64];
	unsigned int datalen;
	unsigned int bitlen[2];
	unsigned int state[5];
	unsigned int k[4];
} SHA_CTX;
#define SHA_DIGEST_LENGTH 20
void SHA1_Transform(SHA_CTX * ctx, unsigned char data[]);
void SHA1_Init(SHA_CTX * ctx);
void SHA1_Update(SHA_CTX * ctx, unsigned char data[], unsigned long len);
void SHA1_Final(unsigned char hash[], SHA_CTX * ctx);

#else




// DBL_INT_ADD treats two unsigned ints a and b as one 64-bit integer and adds c to it
#define ROTLEFT(a,b) ((a << b) | (a >> (32-b)))
#define DBL_INT_ADD(a,b,c) if (a > 0xffffffff - c) ++b; a += c;

void SHA1_Transform(SHA_CTX * ctx, unsigned char data[])
{
	unsigned int a, b, c, d, e, i, j, t, m[80];

	for (i = 0, j = 0; i < 16; ++i, j += 4)
		m[i] = (data[j] << 24) + (data[j + 1] << 16) + (data[j + 2] << 8) + (data[j + 3]);
	for (; i < 80; ++i) {
		m[i] = (m[i - 3] ^ m[i - 8] ^ m[i - 14] ^ m[i - 16]);
		m[i] = (m[i] << 1) | (m[i] >> 31);
	}

	a = ctx->state[0];
	b = ctx->state[1];
	c = ctx->state[2];
	d = ctx->state[3];
	e = ctx->state[4];

	for (i = 0; i < 20; ++i) {
		t = ROTLEFT(a, 5) + ((b & c) ^ (~b & d)) + e + ctx->k[0] + m[i];
		e = d;
		d = c;
		c = ROTLEFT(b, 30);
		b = a;
		a = t;
	}
	for (; i < 40; ++i) {
		t = ROTLEFT(a, 5) + (b ^ c ^ d) + e + ctx->k[1] + m[i];
		e = d;
		d = c;
		c = ROTLEFT(b, 30);
		b = a;
		a = t;
	}
	for (; i < 60; ++i) {
		t = ROTLEFT(a, 5) + ((b & c) ^ (b & d) ^ (c & d)) + e + ctx->k[2] + m[i];
		e = d;
		d = c;
		c = ROTLEFT(b, 30);
		b = a;
		a = t;
	}
	for (; i < 80; ++i) {
		t = ROTLEFT(a, 5) + (b ^ c ^ d) + e + ctx->k[3] + m[i];
		e = d;
		d = c;
		c = ROTLEFT(b, 30);
		b = a;
		a = t;
	}

	ctx->state[0] += a;
	ctx->state[1] += b;
	ctx->state[2] += c;
	ctx->state[3] += d;
	ctx->state[4] += e;
}
void SHA1_Init(SHA_CTX * ctx)
{
	ctx->datalen = 0;
	ctx->bitlen[0] = 0;
	ctx->bitlen[1] = 0;
	ctx->state[0] = 0x67452301;
	ctx->state[1] = 0xEFCDAB89;
	ctx->state[2] = 0x98BADCFE;
	ctx->state[3] = 0x10325476;
	ctx->state[4] = 0xc3d2e1f0;
	ctx->k[0] = 0x5a827999;
	ctx->k[1] = 0x6ed9eba1;
	ctx->k[2] = 0x8f1bbcdc;
	ctx->k[3] = 0xca62c1d6;
}
void SHA1_Update(SHA_CTX * ctx, unsigned char data[], unsigned long len)
{
	unsigned long i;

	for (i = 0; i < len; ++i) {
		ctx->data[ctx->datalen] = data[i];
		ctx->datalen++;
		if (ctx->datalen == 64) {
			SHA1_Transform(ctx, ctx->data);
			DBL_INT_ADD(ctx->bitlen[0], ctx->bitlen[1], 512);
			ctx->datalen = 0;
		}
	}
}
void SHA1_Final(unsigned char hash[], SHA_CTX * ctx)
{
	unsigned int i;

	i = ctx->datalen;

	// Pad whatever data is left in the buffer. 
	if (ctx->datalen < 56) {
		ctx->data[i++] = 0x80;
		while (i < 56)
			ctx->data[i++] = 0x00;
	} else {
		ctx->data[i++] = 0x80;
		while (i < 64)
			ctx->data[i++] = 0x00;
		SHA1_Transform(ctx, ctx->data);
		memset(ctx->data, 0, 56);
	}

	// Append to the padding the total message's length in bits and transform. 
	DBL_INT_ADD(ctx->bitlen[0], ctx->bitlen[1], 8 * ctx->datalen);
	ctx->data[63] = (unsigned char) (ctx->bitlen[0]);
	ctx->data[62] = (unsigned char) (ctx->bitlen[0] >> 8);
	ctx->data[61] = (unsigned char) (ctx->bitlen[0] >> 16);
	ctx->data[60] = (unsigned char) (ctx->bitlen[0] >> 24);
	ctx->data[59] = (unsigned char) (ctx->bitlen[1]);
	ctx->data[58] = (unsigned char) (ctx->bitlen[1] >> 8);
	ctx->data[57] = (unsigned char) (ctx->bitlen[1] >> 16);
	ctx->data[56] = (unsigned char) (ctx->bitlen[1] >> 24);
	SHA1_Transform(ctx, ctx->data);

	// Since this implementation uses little endian byte ordering and MD uses big endian, 
	// reverse all the bytes when copying the final state to the output hash. 
	for (i = 0; i < 4; ++i) {
		hash[i] = (ctx->state[0] >> (24 - i * 8)) & 0x000000ff;
		hash[i + 4] = (ctx->state[1] >> (24 - i * 8)) & 0x000000ff;
		hash[i + 8] = (ctx->state[2] >> (24 - i * 8)) & 0x000000ff;
		hash[i + 12] = (ctx->state[3] >> (24 - i * 8)) & 0x000000ff;
		hash[i + 16] = (ctx->state[4] >> (24 - i * 8)) & 0x000000ff;
	}
}
#undef ROTLEFT
#undef DBL_INT_ADD

#else

// dummy-functions for WINDOWS32 (cpp symbol INLA_WINDOWS32 defined)

int GMRFLib_graph_add_sha1(GMRFLib_graph_tp * g, int skip_sha1)
{
	g->sha1 = NULL;
	return GMRFLib_SUCCESS;
}
void SHA1_Transform(SHA_CTX * ctx, unsigned char data[])
{
}
void SHA1_Init(SHA_CTX * ctx)
{
}
void SHA1_Update(SHA_CTX * ctx, unsigned char data[], unsigned long len)
{
}
void SHA1_Final(unsigned char hash[], SHA_CTX * ctx)
{
}

#endif
