//  this is just scrap-code 

if (irow < n_row) {
	int low = irow;
	int high = n_row - 1;
	int mid;
	do {
		mid = low + (high - low) / 2;
		// printf("BEFORE low %d mid %d high %d n_row %d\n", low, mid, high, n_row);
		if (row_elm[mid].idx < kk) {
			low = mid + 1;
		} else if (row_elm[mid].idx > kk) {
			high = mid - 1;
		} else {
			high = low = mid;
		}
		// printf("AFTER low %d mid %d high %d n_row %d\n", low, mid, high, n_row);
	} while (high > low);
	irow = low;
	P(irow);
	P(n_row);
	P(row_elm[irow].idx);
	P(kk);
	assert(irow > n_row - 1 || row_elm[irow].idx == kk || row_elm[irow].idx > kk);
}
















if (0) {
	int ir;
	int ir2;
	int step;

	if (0) {
		step = 256;
		for (ir = irow + step; ir < nrow && row_elm[ir].idx < kk; ir += step);
		irow = ir - step;

		step = 32;
		for (ir = irow + step; ir < nrow && row_elm[ir].idx < kk; ir += step);
		irow = ir - step;

		step = 4;
		for (ir = irow + step; ir < nrow && row_elm[ir].idx < kk; ir += step);
		irow = ir - step;
	}

	step = 2;
	for (ir = irow + step; ir < nrow && row_elm[ir].idx < kk; ir += step);
	assert(ir - step >= irow);

	irow = ir - step;

	while (irow < nrow - 1 && row_elm[irow].idx < kk)
		irow++;
	// printf("ADD kk %d %d\n", kk, row_elm[irow].idx);
}

if (0) {
	int low = irow + 1;
	int high = n_row - 1;
	int mid;
	int idx_mid;
	int iter = 0;
	irow = low;
	while (high > low) {
		iter++;
		mid = low + (high - low) / 2;
		// printf("BEFORE iter low mid high %d %d %d %d\n", iter, low, mid, high);
		idx_mid = row_elm[mid].idx;
		if (idx_mid == kk) {
			irow = high = low = mid;
		} else if (idx_mid < kk) {
			low = mid + 1;
		} else {
			high = mid - 1;
		}
		// printf("AFTER iter low mid high %d %d %d %d\n", iter, low, mid, high);
	}
	irow = low;
	// printf("BREAK with irow %d\n", irow);
}















if (0) {
	// then add and accumate terms using '..._addto'
#pragma omp parallel for private (i, k, kk, j, jj) num_threads(nt)
	for (i = 0; i < nrow; i++) {
		GMRFLib_idxval_tp *row_idxval = NULL;
		GMRFLib_matrix_get_row_idxval(&row_idxval, i, pA);

		for (kk = 0; kk < pAA_pattern[i]->n; kk++) {   /* for(k = 0; k < N; k++) { */
			k = pAA_pattern[i]->idx[kk];
			for (jj = 0; jj < row_idxval->n; jj++) {
				j = row_idxval->store[jj].idx;
				double val = A_idxval[j]->store[kk].val;
				double val_row = row_idxval->store[jj].val;
				GMRFLib_idxval_addto(&(pAA_idxval[i]), k, val_row * val);
			}
		}
		GMRFLib_idxval_free(row_idxval);
	}
} else {
#pragma omp parallel for private (i, k, kk, j, jj) num_threads(nt)
	for (i = 0; i < nrow; i++) {
		GMRFLib_idxval_tp *row_idxval = NULL;
		GMRFLib_matrix_get_row_idxval(&row_idxval, i, pA);
		GMRFLib_idxval_elm_tp *row_elm = row_idxval->store;

		for (jj = 0; jj < pAA_pattern[i]->n; jj++) {
			j = pAA_pattern[i]->idx[jj];

			GMRFLib_idxval_elm_tp *At_elm = At_idxval[j]->store;

			printf("i j %d %d\n", i, j);
			int irow = 0, iAt = 0;
			while (irow < row_idxval->n && iAt < At_idxval[j]->n) {
				printf("\tirow iAt %d %d\n", irow, iAt);
				k = row_elm[irow].idx;
				kk = At_elm[iAt].idx;
				if (k == kk) {
					printf("\t\tEQ\n");
					GMRFLib_idxval_addto(&(pAA_idxval[i]), j, row_elm[irow].val * At_elm[iAt].val);
					irow++;
					iAt++;
				} else if (k < kk) {
					irow++;
				} else {
					iAt++;
				}
			}
		}
		GMRFLib_idxval_free(row_idxval);
	}
}
