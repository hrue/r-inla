#include <stdio.h>

#include "GMRFLib/GMRFLib.h"
#include "GMRFLib/GMRFLibP.h"

static map_strvp *remap_store = NULL;
static int remap_store_use = 1;
static int remap_store_must_init = 1;
static int remap_store_debug = 0;

int GMRFLib_remap_init_store(void)
{
	if (remap_store_use) {
#pragma omp critical (Name_925d0d4cc0c47ce7c7e2c91447db58147e515e4b)
		if (remap_store_must_init) {
			remap_store = Calloc(1, map_strvp);
			map_strvp_init_hint(remap_store, 64);
			remap_store->alwaysdefault = 1;
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

unsigned char *GMRFLib_remap_prettify_sha(unsigned char *sha)
{
	// THIS FUNCTION OVERWRITE SHA
	// we do a non-invertible compression to make it more pretty for output (do not need to be precise)
	int len = 'z' - 'a' + 1;
	for (int i = 0; i < GMRFLib_SHA_DIGEST_LEN; i++) {
		sha[i] = ((int) sha[i] % len) + 'a';
	}
	return sha;
}

int *GMRFLib_remap_get(int *remap, int n, int nrhs)
{
	if (!remap_store_use) {
		return NULL;
	}

	unsigned char *sha = GMRFLib_remap_sha(remap, n, nrhs);
	void **p = map_strvp_ptr(remap_store, (char *) sha);

	if (remap_store_debug) {
		unsigned char *sh = Malloc(GMRFLib_SHA_DIGEST_LEN + 1, unsigned char);
		Memcpy(sh, sha, GMRFLib_SHA_DIGEST_LEN + 1);
		sh = GMRFLib_remap_prettify_sha(sh);
		if (p) {
			printf("[%1d]{%s} remap_store: remap in store\n", omp_get_thread_num(), sh);
		} else {
			printf("[%1d]{%s} remap_store: add remap into store\n", omp_get_thread_num(), sh);
		}
		Free(sh);
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
		map_strvp_set(remap_store, (char *) r->sha, (void *) r);
	}

	Free(re1);
	Free(re2);

	return r->remap;
}

void GMRFLib_remap_print(FILE *fp)
{
	// write out the cache
	if (remap_store_use) {
		double tsiz = 0.0;
		fprintf(fp, "\nContents of remap_store: \n");
		int k = 0;
		map_strvp_storage *ptr = NULL;
		for (ptr = NULL; (ptr = map_strvp_nextptr(remap_store, ptr)) != NULL;) {
			GMRFLib_remap_tp *r = ((GMRFLib_remap_tp *) ptr->value);
			if (r && r->remap) {
				unsigned char *sh = Malloc(GMRFLib_SHA_DIGEST_LEN + 1, unsigned char);
				Memcpy(sh, r->sha, GMRFLib_SHA_DIGEST_LEN + 1);
				sh = GMRFLib_remap_prettify_sha(sh);
				sh = GMRFLib_remap_prettify_sha(sh);
				fprintf(fp, "\tSlot[%1d] sha[%s] n[%1d] rhs[%1d] remap[%1d %1d %1d...]\n",
					k, sh, r->n, r->nrhs, r->remap[0], r->remap[IMIN(r->n - 1, 1)], r->remap[IMIN(r->n - 1, 2)]);
				tsiz += (r->n * r->nrhs + 2) * sizeof(int);
				Free(sh);
			}
			k++;
		}
		fprintf(fp, "\tTotal size[%.2f]Mb\n", tsiz / SQR(1024.0));
	}
}

void GMRFLib_remap_reset(void)
{
	if (remap_store_use) {
		map_strvp_storage *ptr = NULL;
		for (ptr = NULL; (ptr = map_strvp_nextptr(remap_store, ptr)) != NULL;) {
			GMRFLib_remap_tp *r = ((GMRFLib_remap_tp *) ptr->value);
			if (r && r->remap) {
				Free(r->remap);
				Free(r->sha);
				Free(r);
			}
		}
		Free(remap_store);
		remap_store = NULL;
		remap_store_must_init = 1;
	}
}
