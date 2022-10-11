#include <stdio.h>
#include <malloc.h>
#include <stdlib.h>
#include <omp.h>
#include "sparse-vec.h"

#define UNIFORM() (rand() / ((double) RAND_MAX + 1.0))
#define SEED(seed_) srand((unsigned int) (seed_))

int main(int argc, char **argv)
{
	int n, m;
	if (argc >= 3) {
		n = atoi(argv[1]);
		m = atoi(argv[2]);
	} else {
		printf("Usage:  length repititions\n");
		exit(0);
	}

	GMRFLib_sparse_vec_tp *sparse_vec = NULL;

	SEED(n);
	double *xx = (double *) calloc(n, sizeof(double));
	for (int i = 0; i < n; i++) {
		xx[i] = UNIFORM();
	}

	for (int i = 0, j = 0; i < n; i++) {
		j += 1 + (UNIFORM() < 0.9 ? 0 : 1 + (int) (UNIFORM() * 31));
		if (j >= n) {
			break;
		}
		GMRFLib_sparse_vec_add(&sparse_vec, j, UNIFORM());
	}
	GMRFLib_sparse_vec_prepare(sparse_vec);

	double sum1 = 0.0, sum2 = 0.0;
	double tref1 = 0.0, tref2 = 0.0;

	for (int k = 0; k < m; k++) {
		sum1 = sum2 = 0.0;
		tref1 -= omp_get_wtime();
		sum1 = GMRFLib_sparse_vec_dot_product(sparse_vec, xx);
		tref1 += omp_get_wtime();

		tref2 -= omp_get_wtime();
		for (int j = 0; j < sparse_vec->n; j++) {
			sum2 += sparse_vec->val[j] * xx[sparse_vec->idx[j]];
		}
		tref2 += omp_get_wtime();
		if (fabs(sum1 - sum2) > 1e-8) {
			printf("ERROR: sum1 %g sum2 %g\n", sum1, sum2);
			exit(1);
		}
	}
	printf("dot_product %.3f plain %.3f (%.3f, %.3f)\n", tref1, tref2, tref1 / (tref1 + tref2), tref2 / (tref1 + tref2));

	GMRFLib_sparse_vec_free(sparse_vec);

	exit(0);
}
