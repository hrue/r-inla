#include <stdio.h>

#include "GMRFLib/GMRFLib.h"
#include "GMRFLib/GMRFLibP.h"

static map_strvp remap_store;
static int remap_store_use = 0;
static int remap_store_must_init = 1;
static int remap_store_debug = 0;

int GMRFLib_remap_init_store(void)
{
	if (remap_store_use) {
		if (remap_store_must_init) {
			map_strvp_init_hint(&remap_store, 128);
			remap_store.alwaysdefault = 1;
			remap_store_must_init = 0;
			if (remap_store_debug) {
				printf("\tremap_store: init storage\n");
			}
		}
	}

	return GMRFLib_SUCCESS;
}

unsigned char *GMRFLib_remap_sha(int *remap, int n, int nrhs)
{
	GMRFLib_SHA_TP c;
	unsigned char *md = Malloc(GMRFLib_SHA_DIGEST_LEN + 1, unsigned char);
	GMRFLib_SHA_Init(&c);
	GMRFLib_SHA_IUPDATE(remap, n, c);
	GMRFLib_SHA_IUPDATE(&n, 1, c);
	GMRFLib_SHA_IUPDATE(&nrhs, 1, c);
	GMRFLib_SHA_Final(md, &c);
	md[GMRFLib_SHA_DIGEST_LEN] = '\0';

	return (md);
}

int *GMRFLib_remap_get(int *remap, int n, int nrhs)
{
	if (!remap_store_use) {
		return NULL;
	}

	unsigned char *sha = GMRFLib_remap_sha(remap, n, nrhs);
	void **p = map_strvp_ptr(&remap_store, (char *) sha);

	if (remap_store_debug) {
		if (p) {
			printf("\t[%1d] remap_store: remap in store\n", omp_get_thread_num());
		} else {
			printf("\t[%1d] remap_store: add remap into store\n", omp_get_thread_num());
		}
	}

	if (p) {
		// all good, we have it from before
		Free(sha);
		GMRFLib_remap_tp *r = *((GMRFLib_remap_tp **) p);
		return r->remap;
	}

	int *re = Malloc(n * nrhs, int);
	int *re1 = Malloc(n * nrhs, int);
	int *re2 = Malloc(n * nrhs, int);

	// two step mapping
	for (int j = 0; j < nrhs; j++) {
		int offset = j * n;
		for (int i = 0; i < n; i++) {
			re1[offset + i] = remap[i] + offset;
			re2[offset + i] = i * nrhs + j;
		}
	}
	for (int k = 0; k < n * nrhs; k++) {
		re[k] = re2[re1[k]];
	}

	GMRFLib_remap_tp *r = Calloc(1, GMRFLib_remap_tp);
	r->sha = sha;
	r->n = n;
	r->nrhs = nrhs;
	r->remap = re;

#pragma omp critical (Name_71dc250ae8a03e0bd798461c633f37625101e6b8)
	{
		map_strvp_set(&remap_store, (char *) r->sha, (void *) r);
	}

	Free(re1);
	Free(re2);

	return r->remap;
}
