#include <stdio.h>
#include <omp.h>

#include "GMRFLib/GMRFLib.h"

static map_strvp *remap_store = NULL;
static int remap_store_use = 1;
static int remap_store_must_init = 1;
static int remap_store_debug = 0;

int GMRFLib_remap_init_store(void)
{
	if (remap_store_use) {
#pragma omp critical (Name_925d0d4cc0c47ce7c7e2c91447db58147e515e4b)
		if (remap_store_must_init) {
			if (remap_store_debug) {
				printf("\tremap_store: init storage\n");
			}
			remap_store = Calloc(1, map_strvp);
			map_strvp_init_hint(remap_store, 64);
			remap_store->alwaysdefault = 1;
			remap_store_must_init = 0;
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
	int numa_node = -1;
	GMRFLib_numa_get(NULL, &numa_node);
	GMRFLib_SHA_IUPDATE(&numa_node, 1, c);
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
	void **p = map_strvp_ptr(remap_store, (char *) sha);

	if (remap_store_debug) {

		unsigned char *sh = GMRFLib_prettify_sha(Strdup_sha(sha));
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
#pragma omp atomic
		r->count++;
		return r->remap;
	}

	int numa_node = -1;
	GMRFLib_numa_get(NULL, &numa_node);
	int *re = (int *) GMRFLib_numa_alloc_onnode(n * nrhs * sizeof(int), numa_node);
	int *re1 = Calloc(n * nrhs, int);
	int *re2 = Calloc(n * nrhs, int);

	// two step mapping
	for (int j = 0; j < nrhs; j++) {
		int offset = j * n;
#pragma omp simd
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
	r->numa_node = numa_node;
	r->count = 1;
	r->remap = re;

#pragma omp critical (Name_71dc250ae8a03e0bd798461c633f37625101e6b8)
	map_strvp_set(remap_store, (char *) r->sha, (void *) r);

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
				int nn = r->n * r->nrhs;
				fprintf(fp, "\tSlot[%2.2d] n[%1d] rhs[%1d] numa.node[%1d] count[%1d] remap[%1d %1d %1d...]\n",
					k, r->n, r->nrhs, r->numa_node, r->count, r->remap[0], r->remap[IMIN(nn - 1, 1)],
					r->remap[IMIN(nn - 1, 2)]);
				tsiz += (r->n * r->nrhs + 2) * sizeof(int);
			}
			k++;
		}
		fprintf(fp, "\tTotal size[%.2fMb]\n", tsiz / SQR(1024.0));
	}
}

void GMRFLib_remap_reset(void)
{
	if (remap_store_use) {
		map_strvp_storage *ptr = NULL;
		for (ptr = NULL; (ptr = map_strvp_nextptr(remap_store, ptr)) != NULL;) {
			GMRFLib_remap_tp *r = ((GMRFLib_remap_tp *) ptr->value);
			if (r && r->remap) {
				GMRFLib_numa_free((void *) r->remap, r->n * r->nrhs * sizeof(int));
				Free(r->sha);
				Free(r);
			}
		}
		Free(remap_store);
		remap_store = NULL;
		remap_store_must_init = 1;
	}
}
