#include <stdio.h>

#include "GMRFLib/GMRFLib.h"
#include "GMRFLib/GMRFLibP.h"

static int remap_store_use = 1;
static map_strvp remap_store;
static int remap_store_must_init = 1;
static int remap_store_debug = 1;

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
	unsigned char *md = Calloc(GMRFLib_SHA_DIGEST_LEN + 1, unsigned char);

	Memset(md, 0, GMRFLib_SHA_DIGEST_LEN + 1);
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
	void **p = NULL;

	if (!sha) {
		return NULL;
	}
	
	p = map_strvp_ptr(&remap_store, (char *) sha);
	if (remap_store_debug) {
		if (p) {
			printf("\t[%1d] remap_store: remap is found in store: do not add\n", omp_get_thread_num());
		} else {
			printf("\t[%1d] remap_store: remap is not found in store: add\n", omp_get_thread_num());

			for(int i = 0; i < n; i++) printf("Enter with remap[%1d] = %d\n", i, remap[i]);
		}
	}

	if (p) {
		// all good, we have it from before
		Free(sha);
		GMRFLib_remap_tp *r = *((GMRFLib_remap_tp **) p);
		return r->remap;
	}
	
	int *r1 = Malloc(n * nrhs, int);
	int *r2 = Malloc(n * nrhs, int);

	Memcpy(r1, remap, n * sizeof(int));
	for(int j = 1; j < nrhs; j++) {
		int offset = j * n;
		for(int i = 0; i < n; i++) {
			r1[i + offset] = r1[i] + offset;
		}
	}

	if(1) {
		int *r3 = Malloc(n * nrhs, int);
		for (int j = 0; j < nrhs; j++) {
			for (int i = 0; i < n; i++) {
				r3[i * nrhs + j] = j * n + i;
				//r3[j * n + i] = i * nrhs + j;
			}
		}
		for (int k = 0; k < n * nrhs; k++) {
			r2[k] = r3[r1[k]];
		}
		Free(r3);
	} else {
		Memcpy(r2, r1, n * nrhs * sizeof(int));
	}

	GMRFLib_remap_tp * r = Calloc(1, GMRFLib_remap_tp);
	r->sha = sha;
	r->n = n;
	r->nrhs = nrhs;
	r->remap = r2;
	Free(r1);
	
#pragma omp critical (Name_71dc250ae8a03e0bd798461c633f37625101e6b8)
	{
		map_strvp_set(&remap_store, (char *) r->sha, (void *) r);
	}

	return r->remap;
}

void icopy(int *n, int *x, int *ix, int *y, int *iy) 
{
	// same as dcopy_(). seems like icopy_() is not std in blas
	for(int i = 0, iix = 0, iiy = 0;  i < *n; i++) {
		y[iiy] = x[iix];
		iix += *ix;
		iiy += *iy;
	}
}
