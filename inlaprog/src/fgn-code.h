//
// the code that interpolate in a static table
//

{
	int idx, i, len_par;
	double weight, tmp;

	assert(k == K);
	// make sure its in the range. Hack...
	H_intern = DMAX(H_intern, H_intern_start + H_intern_by);
	H_intern = DMIN(H_intern, H_intern_end - H_intern_by);

	assert(H_intern >= H_intern_start && H_intern <= H_intern_end);
	idx = (int) floor((H_intern - H_intern_start) / H_intern_by);	/* idx is the block-index */
	weight = (H_intern - (H_intern_start + idx * H_intern_by)) / H_intern_by;
	len_par = 2 * K - 1;
	idx *= len_par;					       /* and now the index in the table */

	double *fit_par = Calloc(len_par, double);
	for (i = 0; i < len_par; i++) {
		fit_par[i] = (1.0 - weight) * param[idx + i] + weight * param[idx + len_par + i];
	}

	// the first K are phi
	for (i = 0, tmp = 0.0; i < K; i++) {
		tmp += exp(-fit_par[i]);
		phi[i] = 1.0 / (1.0 + tmp);
	}

	// the remaining K-1 are the weights
	double psum, *par = Calloc(len_par, double);
	par[0] = psum = 1;
	for (i = 1; i < K; i++) {
		par[i] = exp(fit_par[K + (i - 1)]);
		psum += par[i];
	}
	for (i = 0; i < K; i++) {
		w[i] = par[i] / psum;
	}

	Free(fit_par);
	Free(par);

	return GMRFLib_SUCCESS;
}
