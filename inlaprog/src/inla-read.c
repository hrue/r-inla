
/* inla-read.c
 * 
 * Copyright (C) 2007-2023 Havard Rue
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or (at
 * your option) any later version.
 *  * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 *
 * The author's contact information:
 *
 *        Haavard Rue
 *        CEMSE Division
 *        King Abdullah University of Science and Technology
 *        Thuwal 23955-6900, Saudi Arabia
 *        Email: haavard.rue@kaust.edu.sa
 *        Office: +966 (0)12 808 0640
 *
 */

int inla_read_data_all(double **x, int *n, const char *filename, int *ncol_data_all)
{
	if (ncol_data_all)
		*ncol_data_all = -1;			       /* say that it is not defined */

	if (!filename) {
		/*
		 * useful for ini-files with no weight file. (backward compatability...)
		 */
		*n = 0;
		*x = NULL;
		return GMRFLib_SUCCESS;
	}

	int count = 0, err, len = 1024;
	GMRFLib_io_tp *io = NULL;

	if (GMRFLib_is_fmesher_file(filename, (long int) 0, -1) == GMRFLib_SUCCESS) {
		/*
		 * This is the binary-file interface 
		 */
		GMRFLib_matrix_tp *M = GMRFLib_read_fmesher_file(filename, (long int) 0, -1);
		assert(M->elems == M->nrow * M->ncol);	       /* no sparse matrix! */

		*n = M->nrow * M->ncol;
		*x = Calloc(*n, double);

		int i, j, k;
		for (i = k = 0; i < M->nrow; i++) {
			for (j = 0; j < M->ncol; j++) {
				(*x)[k++] = M->A[i + j * M->nrow];
			}
		}

		if (ncol_data_all)
			*ncol_data_all = M->ncol;

		GMRFLib_matrix_free(M);

		return INLA_OK;
	} else {
		double *c = Calloc(len, double);
		GMRFLib_EWRAP0(GMRFLib_io_open(&io, filename, "r"));
		{
			GMRFLib_error_handler_tp *old_handler = GMRFLib_set_error_handler_off();

			while (1) {
				err = GMRFLib_io_read_next(io, &c[count], "%lf");
				if (err == GMRFLib_SUCCESS) {
					count++;
					if (count >= len) {
						len += 1024;
						c = Realloc(c, len, double);
					}
				} else {
					break;
				}
			}
			GMRFLib_set_error_handler(old_handler);
		}
		GMRFLib_EWRAP0(GMRFLib_io_close(io));

		*n = count;
		if (count == 0) {
			Free(c);
			*x = NULL;
		} else {
			*x = c;
		}

		return INLA_OK;
	}
}

int inla_read_data_general(double **xx, int **ix, int *nndata, const char *filename, int n, int column, int n_columns, int verbose,
			   double default_value)
{
	/*
	 * read a column from file. the first (or first two) columns are indices and these are not counted.
	 *
	 * counting starts from 0. unread values are default set to DEFAULT_VALUE. n_columns is the total number of columns
	 * (except the indices).
	 *
	 * if xx exists, then read into xx. otherwise, read into ix. return number of reads, in *nx
	 */
	int nx = 0, ndata, i, j, ii, ncol_true;
	double *x = NULL;

	assert(xx || ix);
	/*
	 * first read all entries in the file 
	 */
	inla_read_data_all(&x, &nx, filename, NULL);
	if (nx) {
		assert(x);
	}
	if (verbose) {
		printf("\t\tread n=[%1d] entries from file=[%s]\n", nx, filename);
	}
	ncol_true = n_columns + 1;
	if (!inla_divisible(nx, ncol_true)) {
		inla_error_file_numelm(__GMRFLib_FuncName, filename, nx, ncol_true);
	}
	ndata = nx / ncol_true;
	if (xx) {
		*xx = Calloc(n, double);
		for (i = 0; i < n; i++) {
			(*xx)[i] = default_value;
		}
	} else {
		*ix = Calloc(n, int);
		for (i = 0; i < n; i++) {
			(*ix)[i] = (int) default_value;
		}
	}

	for (i = j = 0; i < nx; i += ncol_true, j++) {
		ii = (int) x[i];
		if (!LEGAL(ii, n)) {
			inla_error_file_error(__GMRFLib_FuncName, filename, nx, i, x[i]);
		}
		if (xx) {
			(*xx)[ii] = x[i + column + 1];
		} else {
			(*ix)[ii] = (int) x[i + column + 1];
		}
		if (verbose && j < PREVIEW) {
			printf("\t\tfile=[%s] %1d/%1d  (idx,y) = (%1d, %g)\n", filename, j, ndata, ii, x[i + column + 1]);
		}
	}
	if (nndata) {
		*nndata = ndata;
	}

	Free(x);
	return INLA_OK;
}

int inla_sread_str_int(char **tag, int *i, const char *str)
{
	char *strtok_ptr = NULL, *token;
	const char *delim = ":";
	int ok = 1;
	char *p = GMRFLib_strdup(str);

	token = GMRFLib_strtok_r(p, delim, &strtok_ptr);
	if (token) {
		p = NULL;
		*tag = GMRFLib_strdup(token);
		token = GMRFLib_strtok_r(p, delim, &strtok_ptr);
		if (sscanf(token, "%d", i) != 1)
			ok = 0;
	} else {
		ok = 0;
	}
	return (ok ? GMRFLib_SUCCESS : !GMRFLib_SUCCESS);
}

int inla_sread_colon_ints(int *i, int *j, const char *str)
{
	/*
	 * read integer I and J from STR using format I:J, I,J or I J
	 */
	if (sscanf(str, "%d:%d", i, j) == 2) {
		return INLA_OK;
	} else if (sscanf(str, "%d,%d", i, j) == 2) {
		return INLA_OK;
	} else if (sscanf(str, "%d %d", i, j) == 2) {
		return INLA_OK;
	} else {
		return INLA_FAIL;
	}
}

int inla_sread(void *x, int nx, const char *str, int code)
{
	/*
	 * code = 0: int. code = 1: double 
	 */

	char *strtok_ptr = NULL, *token, *p;
	const char *delim = " \t";
	double *dx = (double *) x;
	int *ix = (int *) x;
	int count = 0;
	const int debug = 0;
	int ok;

	if (debug)
		printf("read %d entries from %s\n", nx, str);

	assert(code == 0 || code == 1);
	p = GMRFLib_strdup(str);
	while ((token = GMRFLib_strtok_r(p, delim, &strtok_ptr))) {
		p = NULL;
		ok = 1;
		if (debug) {
			printf("strip [%s] into [%s]\n", str, token);
		}

		if (code == 0) {
			if (sscanf(token, "%d", &ix[count]) == 0) {
				ok = 0;
			}
		} else if (code == 1) {
			if (sscanf(token, "%lf", &dx[count]) == 0) {
				ok = 0;
			}
		}
		if (ok)
			count++;

		if (count == nx) {
			break;
		}
	}

	if (count != nx) {
		return INLA_FAIL;
	}
	Free(p);

	return INLA_OK;
}

int inla_sread_q(void **x, int *nx, const char *str, int code)
{
	/*
	 * code = 0: int. code = 1: double
	 * 
	 * this return the number of `ints' read in nx and in x 
	 */

	char *strtok_ptr = NULL, *token, *p;
	const char *delim = " \t";
	double *dx = NULL;
	double dx_try;
	int *ix = NULL;
	int count = 0;
	const int debug = 0;
	int ix_try;
	int ok;

	assert(code == 0 || code == 1);
	p = GMRFLib_strdup(str);

	while ((token = GMRFLib_strtok_r(p, delim, &strtok_ptr))) {
		p = NULL;
		ok = 0;
		if (debug) {
			printf("strip [%s] into [%s]\n", str, token);
		}
		if (code == 0) {
			if (sscanf(token, "%d", &ix_try) == 1)
				ok = 1;
		} else {
			if (sscanf(token, "%lf", &dx_try) == 1)
				ok = 1;
		}
		if (ok) {
			if (code == 0) {
				ix = Realloc(ix, count + 1, int);
				ix[count++] = ix_try;
			} else {
				dx = Realloc(dx, count + 1, double);
				dx[count++] = dx_try;
			}
		}
	}

	*nx = count;
	if (count == 0) {
		*x = NULL;
	} else {
		if (code == 0) {
			*x = (void *) ix;
		} else {
			*x = (void *) dx;
		}
	}

	if (debug) {
		int i;
		for (i = 0; i < *nx; i++) {
			if (code == 0) {
				printf("%s : %d %d\n", str, i, ix[i]);
			} else {
				printf("%s : %d %g\n", str, i, dx[i]);
			}
		}
	}

	Free(p);

	return INLA_OK;
}

int inla_sread_ints(int *x, int nx, const char *str)
{
	// read a fixed number of ints from str
	return inla_sread((void *) x, nx, str, 0);
}

int inla_sread_doubles(double *x, int nx, const char *str)
{
	// read a fixed number of doubles from str
	return inla_sread((void *) x, nx, str, 1);
}

int inla_sread_ints_q(int **x, int *nx, const char *str)
{
	// read an unknown number of ints from str
	return inla_sread_q((void **) x, nx, str, 0);
}

int inla_sread_doubles_q(double **x, int *nx, const char *str)
{
	// read an unknown number of doubles from str
	return inla_sread_q((void **) x, nx, str, 1);
}

int inla_is_NAs(int nx, const char *string)
{
	/*
	 * return GMRFLib_SUCCESS is string consists of nx NA's + whitespace separation 
	 */
	char *scopy, *p;
	const char *sep = " \t", *NA = "NA";
	int k = 0, nna = 0;
	const int debug = 0;

	if (debug) {
		printf("call inla_is_NAs: nx %d string %s\n", nx, string);
	}

	if (!string && nx)
		return !GMRFLib_SUCCESS;
	if (!string && !nx)
		return GMRFLib_SUCCESS;

	scopy = GMRFLib_strdup(string);
	p = strtok(scopy, sep);
	nna += (p && !strcmp(p, NA));

	if (debug)
		printf("get token %d : %s\n", k++, p);

	while (p) {
		p = strtok(NULL, sep);
		nna += (p && !strcmp(p, NA));

		if (debug)
			printf("get token %d : %s\n", k++, p);
	}

	Free(scopy);
	return (nna == nx ? GMRFLib_SUCCESS : !GMRFLib_SUCCESS);
}

const char *inla_string_join(const char *a, const char *b)
{
	/*
	 * join strings A and B into A:B. 
	 */
	char *ans = NULL;

	assert((a ? strlen(a) + 1 : 0) + (b ? strlen(b) + 1 : 0) < 1025);
	GMRFLib_sprintf(&ans, "%s%c%s", (a ? a : ""), INIPARSER_SEP, (b ? b : ""));
	return ans;
}

int inla_error_missing_required_field(const char *funcname, const char *secname, const char *field)
{
	fprintf(stderr, "\n\n*** ERROR *** \t%s: section [%s]: missing required field [%s]\n\n", funcname, secname, field);
	exit(EXIT_FAILURE);
	return INLA_OK;
}

int inla_error_open_file(const char *msg)
{
	fprintf(stderr, "\n\n*** ERROR *** \tfail to open file[%s] for writing. Exit...\n\n", msg);
	exit(EXIT_FAILURE);
	return INLA_OK;
}

int inla_error_general(const char *msg)
{
	fprintf(stderr, "\n\n*** ERROR *** \t%s\n\n", msg);
	exit(EXIT_FAILURE);
	return INLA_OK;
}

int inla_error_general2(const char *msg, const char *msg2)
{
	fprintf(stderr, "\n\n*** ERROR *** \t%s: %s\n\n", msg, msg2);
	exit(EXIT_FAILURE);
	return INLA_OK;
}

int inla_error_field_is_void(const char *funcname, const char *secname, const char *field, const char *value)
{
	fprintf(stderr, "\n\n*** ERROR *** \t%s: section [%s]: field [%s] is void: [%s]\n\n", funcname, secname, field, value);
	exit(EXIT_FAILURE);
	return INLA_OK;
}

int inla_error_file_numelm(const char *funcname, const char *filename, int n, int idiv)
{
	fprintf(stderr, "\n\n*** ERROR *** \t%s: file [%s] contains [%1d] elements, which is not a multiple of [%1d]\n\n", funcname,
		filename, n, idiv);
	exit(EXIT_FAILURE);
	return INLA_OK;
}

int inla_error_file_totnumelm(const char *funcname, const char *filename, int n, int total)
{
	fprintf(stderr, "\n\n*** ERROR *** \t%s: file [%s] contains [%1d] elements, which is different from [n] = [%1d]\n\n",
		funcname, filename, n, total);
	exit(EXIT_FAILURE);
	return INLA_OK;
}

int inla_error_file_error_sorted(const char *funcname, const char *filename, int n, int element_number, double val)
{
	fprintf(stderr,
		"\n\n*** ERROR *** \t%s: file [%s] contains [%1d] elements, but element [%1d] = [%g] is void. The 'values' needs to be increasing!\n\n",
		funcname, filename, n, element_number, val);
	abort();
	exit(EXIT_FAILURE);
	return INLA_OK;
}

int inla_error_file_error(const char *funcname, const char *filename, int n, int element_number, double val)
{
	fprintf(stderr, "\n\n*** ERROR *** \t%s: file [%s] contains [%1d] elements, but element [%1d] = [%g] is void.\n\n",
		funcname, filename, n, element_number, val);
	abort();
	exit(EXIT_FAILURE);
	return INLA_OK;
}

int inla_error_file_error2(const char *funcname, const char *filename, int n, int element_number, double val, int element_number2, double val2)
{
	fprintf(stderr,
		"\n\n*** ERROR *** \t%s: file [%s] contains [%1d] elements, but element [%1d] = [%g] or [%1d] = [%g] is void.\n\n",
		funcname, filename, n, element_number, val, element_number2, val2);
	exit(EXIT_FAILURE);
	return INLA_OK;
}

int inla_read_fileinfo(inla_tp *mb, dictionary *ini, int sec, File_tp *file, const char *FILENAME)
{
	char *tag = (FILENAME ? GMRFLib_strdup(FILENAME) : GMRFLib_strdup("FILENAME"));
	char *secname = GMRFLib_strdup(iniparser_getsecname(ini, sec));

	file->name = GMRFLib_strdup(iniparser_getstring(ini, inla_string_join(secname, tag), NULL));

	/*
	 * Only signal error if FILENAME is NULL. 
	 */
	if (!file->name && !FILENAME) {
		inla_error_missing_required_field(__GMRFLib_FuncName, secname, tag);
	}
	if (mb->verbose) {
		printf("\t\tfile->name=[%s]\n", file->name);
	}
	Free(tag);

	return INLA_OK;
}

int inla_read_weightsinfo(inla_tp *mb, dictionary *ini, int sec, File_tp *file)
{
	char *secname = GMRFLib_strdup(iniparser_getsecname(ini, sec));

	file->name = GMRFLib_strdup(iniparser_getstring(ini, inla_string_join(secname, "FILENAME"), NULL));
	if (!file->name) {
		inla_error_missing_required_field(__GMRFLib_FuncName, secname, "filename");
	}
	if (mb->verbose) {
		printf("\t\tfile->name=[%s]\n", file->name);
	}
	return INLA_OK;
}


int inla_trim_family(char *family)
{
	size_t i, j = 0;

	assert(family);
	for (i = 0, j = 0; i < strlen(family); i++) {
		if (family[i] != '.' && family[i] != '_' && family[i] != ' ' && family[i] != '\t') {
			family[j] = family[i];
			j++;
		}
	}
	family[j] = '\0';
	return GMRFLib_SUCCESS;
}

char *inla_make_tag(const char *string, int ds)
{
	return (inla_make_tag2(string, ds, NULL));
}

char *inla_make_tag2(const char *string, int ds, const char *estr)
{
	char *res;

	if (ds > 0) {
		GMRFLib_sprintf(&res, "%s%s[%1d]", string, (estr ? estr : ""), ds + 1);	/* yes, number these from 1...nds */
	} else {
		GMRFLib_sprintf(&res, "%s%s", string, (estr ? estr : ""));
	}

	return res;
}

GMRFLib_constr_tp *inla_read_constraint(const char *filename, int n)
{
	/*
	 * read constraints from file 
	 */
	double *x = NULL;
	int i, j, m, nc;

	inla_read_data_all(&x, &m, filename, NULL);
	nc = m / (n + 1);				       /* yes, integer division */
	if (nc * n + nc != m) {
		char *msg = NULL;
		GMRFLib_sprintf(&msg, "Number of elements[%1d] in file[%s] does is not a multiplum of n=[%1d]", m, filename, n + 1);
		inla_error_general(msg);
	}
	GMRFLib_constr_tp *c = NULL;

	GMRFLib_make_empty_constr(&c);
	c->nc = nc;
	c->a_matrix = Calloc(nc * n, double);
	c->e_vector = Calloc(nc, double);

	for (j = 0; j < nc; j++) {
		for (i = 0; i < n; i++) {
			c->a_matrix[i * nc + j] = x[i + j * n];
		}
		c->e_vector[j] = x[nc * n + j];
	}
	/*
	 * this is a stupid construction (i blame myself...): i have to pass the graph in order to pass `n'... 
	 */
	GMRFLib_graph_tp *g = NULL;

	GMRFLib_graph_mk_linear(&g, n, 0, 0);
	GMRFLib_prepare_constr(c, g, 0);
	GMRFLib_graph_free(g);
	Free(x);
	return c;
}



int inla_read_graph(const char *filename)
{
	/*
	 * Read a graph and print it on stdio. Compute also the connected components.
	 */
	GMRFLib_graph_tp *graph = NULL;

	GMRFLib_graph_read(&graph, filename);
	GMRFLib_graph_write2(stdout, graph);

	int *cc, i;
	cc = GMRFLib_graph_cc(graph);
	for (i = 0; i < graph->n; i++)
		printf("%d\n", cc[i]);
	Free(cc);

	return 0;
}

inla_file_contents_tp *inla_read_file_contents(const char *filename)
{
	/*
	 * just read the hole file into on long character vector 
	 */

	FILE *fp;
	long len;

	if (!filename) {
		return NULL;
	}
	fp = fopen(filename, "rb");
	if (!fp) {
		return NULL;
	}
	fseek(fp, 0L, SEEK_END);
	len = ftell(fp);
	assert(len > 0L);

	inla_file_contents_tp *fc = Calloc(1, inla_file_contents_tp);
	fc->contents = Calloc((size_t) len, char);

	rewind(fp);
	fc->len = fread(fc->contents, (size_t) 1, (size_t) len, fp);
	assert(fc->len == (size_t) len);
	fclose(fp);

	return fc;
}

int inla_write_file_contents(const char *filename, inla_file_contents_tp *fc)
{
	/*
	 * just dump the file contents to the new file 
	 */

	if (!fc) {
		return INLA_OK;
	}

	FILE *fp;
	size_t len;

	fp = fopen(filename, "wb");
	assert(fp);
	len = fwrite(fc->contents, (size_t) 1, fc->len, fp);
	assert(len == fc->len);

	fclose(fp);
	return INLA_OK;
}
