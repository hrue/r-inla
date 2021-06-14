		if (0) {
			// then add and accumate terms using '..._addto'
#pragma omp parallel for private (i, k, kk, j, jj) num_threads(nt)
			for (i = 0; i < nrow; i++) {
				GMRFLib_idxval_tp *row_idxval = NULL;
				GMRFLib_matrix_get_row_idxval(&row_idxval, i, pA);

				for (kk = 0; kk < pAA_pattern[i]->n; kk++) {	/* for(k = 0; k < N; k++) { */
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
					while(irow < row_idxval->n && iAt < At_idxval[j]->n) {
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
