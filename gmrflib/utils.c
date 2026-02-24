#include <assert.h>
#include <float.h>
#include <limits.h>
#include <math.h>
#include <omp.h>
#include <signal.h>
#include <stdarg.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include <time.h>

#if defined(_WIN32)
#       include <windows.h>
#endif

#if defined(__APPLE__)
#       include <sys/types.h>
#       include <sys/sysctl.h>
#endif

#include "GMRFLib/GMRFLib.h"
#include "GMRFLib/hashP.h"

#define IDX_ALLOC_INITIAL 64
#define IDX_ALLOC_ADD     512

static int malloc_debug = -1;
static int malloc_debug_show_free = 0;			       /* print 'free' or not */

// better with function than macro...
char *Strdup(const char *s)
{
	return (s ? strdup(s) : (char *) NULL);
}

unsigned char *Strdup_sha(unsigned char *sha)
{
	if (!sha) {
		return NULL;
	}
	unsigned char *nnew = Calloc(GMRFLib_SHA_DIGEST_LEN + 1, unsigned char);
	Memcpy(nnew, sha, GMRFLib_SHA_DIGEST_LEN);
	nnew[GMRFLib_SHA_DIGEST_LEN] = '\0';

	return nnew;
}

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wattributes"
__attribute__((target_clones(INLA_CLONE_TARGETS "default")))
unsigned char *GMRFLib_prettify_sha(unsigned char *sha)
{
	if (!sha) {
		return NULL;
	}
	// THIS FUNCTION OVERWRITE SHA
	// we do a non-invertible compression to make it more pretty for output (do not need to be precise)
	int len = 'z' - 'a' + 1;
	for (int i = 0; i < GMRFLib_SHA_DIGEST_LEN; i++) {
		sha[i] = ((int) sha[i] % len) + 'a';
		if (sha[i] == '\0') {
			// otherwise strlen() etc cannot be used
			sha[i] = 'z';
		}
	}
	sha[GMRFLib_SHA_DIGEST_LEN] = '\0';

	return sha;
}
#pragma GCC diagnostic pop


/*
 * Measures the current (and peak) resident and virtual memories
 * usage of your linux C process, in kB
 *
 * taken from
 * https://stackoverflow.com/questions/1558402/memory-usage-of-current-process-in-c
 */
void GMRFLib_getMemory(int POSSIBLY_UNUSED(*currRealMem), int POSSIBLY_UNUSED(*peakRealMem),
		       int POSSIBLY_UNUSED(*currVirtMem), int POSSIBLY_UNUSED(*peakVirtMem))
{
#if defined(__linux__)
	// stores each word in status file
	char buffer[1024] = "";

	// linux file contains this-process info
	FILE *file = fopen("/proc/self/status", "r");

	if (file) {
		// read the entire file
		while (fscanf(file, " %1023s", buffer) == 1) {

			if (strcmp(buffer, "VmRSS:") == 0) {
				if (fscanf(file, " %d", currRealMem) == 0)
					currRealMem = 0;
			}
			if (strcmp(buffer, "VmHWM:") == 0) {
				if (fscanf(file, " %d", peakRealMem) == 0)
					peakRealMem = 0;
			}
			if (strcmp(buffer, "VmSize:") == 0) {
				if (fscanf(file, " %d", currVirtMem) == 0)
					currVirtMem = 0;
			}
			if (strcmp(buffer, "VmPeak:") == 0) {
				if (fscanf(file, " %d", peakVirtMem) == 0)
					peakVirtMem = 0;
			}
		}
		fclose(file);
	} else {
		currRealMem = peakRealMem = currVirtMem = peakVirtMem = 0;
	}
#endif
}

void GMRFLib_printMem_core(FILE POSSIBLY_UNUSED(*fp), const char POSSIBLY_UNUSED(*fnm), int POSSIBLY_UNUSED(lineno))
{
#if defined(__linux__)
	int crm = 0, prm = 0, cvm = 0, pvm = 0;
	FILE *ffp = (fp ? fp : stdout);
	GMRFLib_getMemory(&crm, &prm, &cvm, &pvm);
	fprintf(ffp, "%s:%d: {cur,peak}-Mem used: Real[%.1f, %.1f]Mb, Virt[%.1f, %.1f]Mb\n",
		fnm, lineno, crm / 1024.0, prm / 1024.0, cvm / 1024.0, pvm / 1024.0);
#endif
}

void GMRFLib_delay(int msec)
{
	long pause;
	clock_t now, then;

	pause = msec * (CLOCKS_PER_SEC / 1000);
	now = then = clock();
	while ((now - then) < pause) {
		now = clock();
	}
}

void GMRFLib_delay_random(int msec_low, int msec_high)
{
	GMRFLib_delay(msec_low + (int) ((msec_high - msec_low) * GMRFLib_uniform()));
}

char *GMRFLib_vec2char(double *x, int len)
{
	size_t estimated_size = IMAX(1, len * 32);	       // More conservative estimate
	char *str = Calloc(estimated_size, char);

	size_t offset = 0;
	for (int i = 0; i < len; i++) {
		int written = snprintf(str + offset, estimated_size - offset,
				       (i < len - 1 ? "%.8g," : "%.8g"), x[i]);
		if (written >= (int) (estimated_size - offset)) {
			estimated_size *= 2;
			str = Realloc(str, estimated_size, char);
		}
		offset += written;
	}

	size_t j = 0;
	for (size_t i = 0; i < strlen(str); i++) {
		if (str[i] != ' ')
			str[j++] = str[i];
	}
	str[j] = '\0';

	return (str);
}

int GMRFLib_sprintf(char **ptr, const char *fmt, ...)
{
	/*
	 * parts of this code is copied from the manual page of snprintf. 
	 */

	int n, size = 128 + 1;
	char *p = NULL;
	va_list ap;

	GMRFLib_ASSERT(ptr, GMRFLib_EINVARG);
	GMRFLib_ASSERT(fmt, GMRFLib_EINVARG);

	p = Calloc(size, char);

	while (1) {
		/*
		 * Try to print in the allocated space. 
		 */
		va_start(ap, fmt);
		n = vsnprintf(p, (unsigned int) size, fmt, ap);
		va_end(ap);

		/*
		 * if that worked, return the string, 
		 */
		if (n > -1 && n < size) {
			*ptr = p;
			return GMRFLib_SUCCESS;
		}

		/*
		 * ...else try again with more space 
		 */
		if (n > -1) {
			size = n + 1;
		} else {
			size *= 2;
		}
		p = Realloc(p, size, char);
	}

	return GMRFLib_SUCCESS;
}

void *GMRFLib_memcpy(void *dest, const void *src, size_t n)
{
	assert(n < PTRDIFF_MAX);
	memcpy(dest, src, n);
	return NULL;
}

void GMRFLib_malloc_debug_check(void)
{
	if (malloc_debug < 0) {
		char *def = getenv("INLA_MALLOC_DEBUG");
		if (def) {
			int val = atoi(def);
			malloc_debug = (val > 0 ? val : 0);
		} else {
			malloc_debug = 0;
		}
	}
}

void *GMRFLib_calloc(size_t nmemb, size_t size, const char *file, const char *funcname, int lineno)
{
	assert(nmemb * size < PTRDIFF_MAX);
	void *ptr = calloc(nmemb, size);

	if (malloc_debug > 0 && nmemb * size >= (size_t) malloc_debug) {
		printf(" *** MALLOC_DEBUG *** %s: %s: %1d: calloc nmemb = %zu size = %zu, got address %p\n", file, funcname, lineno, nmemb, size,
		       ptr);
	}

	if (ptr) {
		return ptr;
	}
	char *msg = NULL;
	GMRFLib_sprintf(&msg, "Failed to calloc nmemb=%1lu elements of size=%1lu bytes", nmemb, size);
	GMRFLib_handle_error(file, funcname, lineno, GMRFLib_EMEMORY, msg);
	abort();

	return NULL;
}

void *GMRFLib_malloc(size_t size, const char *file, const char *funcname, int lineno)
{
	assert(size < PTRDIFF_MAX);
	void *ptr = malloc(size);

	if (malloc_debug > 0 && size >= (size_t) malloc_debug) {
		printf(" *** MALLOC_DEBUG *** %s: %s: %1d: malloc size = %zu, got address %p\n", file, funcname, lineno, size, ptr);
	}

	if (ptr) {
		return ptr;
	}

	char *msg = NULL;
	GMRFLib_sprintf(&msg, "Failed to malloc size=%1lu bytes", size);
	GMRFLib_handle_error(file, funcname, lineno, GMRFLib_EMEMORY, msg);
	abort();

	return NULL;
}

void *GMRFLib_realloc(void *old_ptr, size_t size, const char *file, const char *funcname, int lineno)
{
	assert(size < PTRDIFF_MAX);
	void *ptr = realloc(old_ptr, size);

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wuse-after-free"
	if (malloc_debug > 0 && size >= (size_t) malloc_debug) {
		printf(" *** MALLOC_DEBUG *** %s: %s: %d: realloc size = %zu,  from address %p to address %p\n", file, funcname, lineno, size,
		       old_ptr, ptr);
	}
#pragma GCC diagnostic pop

	if (ptr) {
		return ptr;
	}
	char *msg = NULL;
	GMRFLib_sprintf(&msg, "Failed to realloc size=%1lu bytes", size);
	GMRFLib_handle_error(file, funcname, lineno, GMRFLib_EMEMORY, msg);
	abort();

	return NULL;
}

void GMRFLib_free(void *ptr, const char *file, const char *funcname, int lineno)
{
	if (ptr) {
		if (malloc_debug_show_free && malloc_debug > 0) {
			printf(" *** MALLOC_DEBUG *** %s: %s: %d: free address %p\n", file, funcname, lineno, ptr);
		}
		free(ptr);
	}
}

char *GMRFLib_rindex(const char *p, int ch)
{
	/*
	 * as Windows does not have it... 
	 */
	char *save = NULL, *pp = (char *) p;
	for (save = NULL;; ++pp) {
		if (*pp == ch) {
			save = pp;
		}
		if (!*pp) {
			return (save);
		}
	}
	abort();
	return NULL;
}

int GMRFLib_which(double val, double *array, int len)
{
	/*
	 * return the first index in array such that array[idx] == val, and -1 if not there 
	 */
	int i;

	for (i = 0; i < len; i++) {
		if (ISEQUAL(val, array[i])) {
			return i;
		}
	}
	return -1;
}

int GMRFLib_iwhich_sorted(int key, int *__restrict ix, unsigned int len)
{
	int *p = GMRFLib_bsearch(key, (int) len, ix);
	return (p ? (p - ix) : -1);
}

int GMRFLib_find_nonzero(double *array, int len, int direction)
{
	/*
	 * return the first/last index in array such that array[idx] != 0, and -1 if not there. direction > 0 : look for first. direction < 0 : look for last
	 */
	if (direction >= 0) {
		for (int i = 0; i < len; i++) {
			if (array[i] != 0.0)
				return i;
		}
		return -1;
	} else {
		for (int i = len - 1; i >= 0; i--) {
			if (array[i] != 0.0)
				return i;
		}
		return -1;
	}

	return -1;
}

int GMRFLib_find_value(double *array, int len, int direction, double value)
{
	/*
	 * return the first/last index in array such that array[idx] == value , and -1 if not there. direction > 0 : look for first. direction < 0 : look for last
	 */
	if (direction >= 0) {
		for (int i = 0; i < len; i++) {
			if (array[i] == value)
				return i;
		}
		return -1;
	} else {
		for (int i = len - 1; i >= 0; i--) {
			if (array[i] == value)
				return i;
		}
		return -1;
	}

	return -1;
}

int GMRFLib_find_ivalue(int *iarray, int len, int direction, int ivalue)
{
	/*
	 * return the first/last index in iarray such that iarray[idx] == ivalue , and -1 if not there. direction > 0 : look for first. direction < 0 : look for last
	 */
	if (direction >= 0) {
		for (int i = 0; i < len; i++) {
			if (iarray[i] == ivalue)
				return i;
		}
		return -1;
	} else {
		for (int i = len - 1; i >= 0; i--) {
			if (iarray[i] == ivalue)
				return i;
		}
		return -1;
	}

	return -1;
}

double GMRFLib_eps(double power)
{
	return (exp(GSL_LOG_DBL_EPSILON * power));
	// return (pow(DBL_EPSILON, power));
}

int GMRFLib_print_darray(FILE *fp, double *x, int n, const char *desc)
{
	int i;

	fp = (fp ? fp : stdout);
	fprintf(fp, "Double array with length %1d [%1s]\n", n, (desc ? desc : "no description given"));
	for (i = 0; i < n; i++) {
		fprintf(fp, "\telement %1d = %.10f\n", i, x[i]);
	}
	return GMRFLib_SUCCESS;
}

int GMRFLib_print_iarray(FILE *fp, int *x, int n, const char *desc)
{
	int i;

	fp = (fp ? fp : stdout);
	fprintf(fp, "Integer array with length %1d [%1s]\n", n, (desc ? desc : "no description given"));
	for (i = 0; i < n; i++) {
		fprintf(fp, "\telement %1d = %1d\n", i, x[i]);
	}
	return GMRFLib_SUCCESS;
}

int GMRFLib_icmp(const void *a, const void *b)
{
	const int *ia = (const int *) a;
	const int *ib = (const int *) b;

	return (*ia - *ib);
}

int GMRFLib_icmp_r(const void *a, const void *b)
{
	return (-GMRFLib_icmp(a, b));
}

int GMRFLib_dcmp(const void *a, const void *b)
{
	const double *da = (const double *) a;
	const double *db = (const double *) b;
	return (*da > *db ? 1 : (*da < *db ? -1 : 0));
}

int GMRFLib_dcmp_r(const void *a, const void *b)
{
	return (-GMRFLib_dcmp(a, b));
}

int GMRFLib_dcmp_abs(const void *a, const void *b)
{
	const double *da = NULL, *db = NULL;

	da = (const double *) a;
	db = (const double *) b;

	/*
	 * sort on ABS() 
	 */
	if (ABS(*da) > ABS(*db)) {
		return 1;
	}
	if (ABS(*da) < ABS(*db)) {
		return -1;
	}

	/*
	 * if they're equal, sort on sign 
	 */
	if ((*da) > (*db)) {
		return 1;
	}
	if ((*da) < (*db)) {
		return -1;
	}

	/*
	 * identical 
	 */
	return 0;
}

int GMRFLib_dcmp_abs_r(const void *a, const void *b)
{
	return (-GMRFLib_dcmp_abs(a, b));
}


double GMRFLib_log_apbex(double a, double b)
{
	/*
	 * evaluate log(a + exp(b))
	 */
	return (b + log1p(a / exp(b)));
}

int GMRFLib_normalize(int n, double *x)
{
	// scale x so the sum is 1

	double sum;
	sum = GMRFLib_dsum(n, x);
	GMRFLib_dscale(n, 1.0 / sum, x);

	return GMRFLib_SUCCESS;
}
int GMRFLib_unique_relative(int *n, double *x, double eps)
{
	/*
	 * assume x is sorted, remove ties and change *n accordingly.
	 * 
	 * ties are defined if relative error between x_i and x_j <= eps, roughly, by using the routine gsl_fcmp()
	 * 
	 */
	int i = 0, j = 0, jj = 0;

	while (jj < *n) {
		while (jj + 1 < *n && gsl_fcmp(x[jj + 1], x[j], eps) == 0) {
			jj++;
		}
		x[i] = x[jj];
		i++;
		jj++;
		j = jj;
	}
	*n = i;

	return GMRFLib_SUCCESS;
}

int GMRFLib_unique_relative2(int *n, double *x, double *y, double eps)
{
	/*
	 * assume x is sorted, remove ties and change *n accordingly. make the same changes to y.
	 * 
	 * ties are defined if relative error between x_i and x_j <= eps, roughly, by using the routine gsl_fcmp()
	 * 
	 */

	if (!y) {
		return GMRFLib_unique_relative(n, x, eps);
	}

	int i = 0, j = 0, jj = 0;

	while (jj < *n) {
		while (jj + 1 < *n && gsl_fcmp(x[jj + 1], x[j], eps) == 0) {
			jj++;
		}
		x[i] = x[jj];
		y[i] = y[jj];
		i++;
		jj++;
		j = jj;
	}
	*n = i;

	return GMRFLib_SUCCESS;
}

int GMRFLib_unique_additive(int *n, double *x, double eps)
{
	// assume x is sorted, remove ties and change *n accordingly. use the median in each bin.
	// ties are defined if |x_i and x_j| <= eps

	int ties = 0;
	for (int k = 0; k < *n - 1; k++) {
		if (ABS(x[k] - x[k + 1]) <= eps) {
			ties = 1;
			break;
		}
	}
	if (!ties) {
		return GMRFLib_SUCCESS;
	}

	int i = 0, j = 0, jj = 0;
	while (jj < *n) {
		while (jj + 1 < *n && ABS(x[jj + 1] - x[j]) <= eps) {
			jj++;
		}
		x[i] = x[jj];
		i++;
		jj++;
		j = jj;
	}
	*n = i;

	return GMRFLib_SUCCESS;
}

int GMRFLib_unique_additive2(int *n, double *x, double *y, double eps)
{
	// assume x is sorted, remove ties and change *n accordingly. use the median in each bin. make the same changes to y.
	// ties are defined if |x_i and x_j| <= eps

	if (!y) {
		return GMRFLib_unique_additive(n, x, eps);
	}

	int ties = 0;
	for (int k = 0; k < *n - 1; k++) {
		if (ABS(x[k] - x[k + 1]) <= eps) {
			ties = 1;
			break;
		}
	}
	if (!ties) {
		return GMRFLib_SUCCESS;
	}

	int i = 0, j = 0, jj = 0;
	while (jj < *n) {
		while (jj + 1 < *n && ABS(x[jj + 1] - x[j]) <= eps) {
			jj++;
		}
		x[i] = x[jj];
		y[i] = y[jj];
		i++;
		jj++;
		j = jj;
	}
	*n = i;

	return GMRFLib_SUCCESS;
}

int GMRFLib_printf_matrix(FILE *fp, double *A, int m, int n)
{
	// A is m x n matrix
#pragma omp critical (Name_bb051132870d1f0b90133946052e91194aa163a5)
	{
		fprintf(fp, "\n\n");
		for (int i = 0; i < m; i++) {
			fprintf(fp, "\t");
			for (int j = 0; j < n; j++)
				fprintf(fp, " %10.6f", A[i + j * m]);
			fprintf(fp, "\n");
		}
		fprintf(fp, "\n");
		fflush(fp);
	}
	return 0;
}

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wattributes"
__attribute__((target_clones(INLA_CLONE_TARGETS "default")))
int GMRFLib_gsl_matrix_count_eq(gsl_matrix *A, double value)
{
	int num = 0;
	for (size_t i = 0; i < A->size1; i++) {
		for (size_t j = 0; j < A->size2; j++) {
			num += (ISNAN(value) ? ISNAN(gsl_matrix_get(A, i, j)) : (gsl_matrix_get(A, i, j) == value));
		}
	}
	return num;
}
#pragma GCC diagnostic pop

int GMRFLib_printf_gsl_matrix(FILE *fp, gsl_matrix *matrix, const char *format)
{
	size_t i, j;

	for (i = 0; i < matrix->size1; i++) {
		for (j = 0; j < matrix->size2; j++) {
			fprintf(fp, (format ? format : " %g"), gsl_matrix_get(matrix, i, j));
		}
		fprintf(fp, "\n");
	}
	return GMRFLib_SUCCESS;
}

int GMRFLib_printf_gsl_matrix2(FILE *fp, gsl_matrix *matrix, const char *format, double cutoff)
{
	size_t i, j;

	for (i = 0; i < matrix->size1; i++) {
		for (j = 0; j < matrix->size2; j++) {
			double a = gsl_matrix_get(matrix, i, j);
			if (ABS(a) > cutoff) {
				fprintf(fp, (format ? format : " %g"), gsl_matrix_get(matrix, i, j));
			} else {
				fprintf(fp, "\t %s", "  . ");
			}
		}
		fprintf(fp, "\n");
	}
	return GMRFLib_SUCCESS;
}

int GMRFLib_printf_gsl_vector(FILE *fp, gsl_vector *vector, const char *format)
{
	size_t i;

	for (i = 0; i < vector->size; i++) {
		fprintf(fp, (format ? format : " %g"), gsl_vector_get(vector, i));
		fprintf(fp, "\n");
	}
	return GMRFLib_SUCCESS;
}

double GMRFLib_signed_pow(double x, double power)
{
	if (ISZERO(x)) {
		return 0.0;
	} else {
		return (x > 0.0 ? 1.0 : -1.0) * pow(ABS(x), power);
	}
}

int GMRFLib_2order_poleq(double *sol1, double *sol2, double a, double b, double c)
{
	/*
	 * solve the equation a*x^2 + b*x + c,
	 * 
	 * return the two solutions in sol1 and sol2. return GMRFLib_EMISC if there is complex solutions. 
	 */

	double tmp = SQR(b) - 4.0 * a * c;

	if (tmp < 0.0) {
		return GMRFLib_EMISC;
	}

	if (ISZERO(a)) {
		if (sol1) {
			*sol1 = -c / b;
		}
		if (sol2) {
			*sol2 = -c / b;
		}
	} else {
		if (sol1) {
			*sol1 = (-b + sqrt(tmp)) / (2.0 * a);
		}
		if (sol2) {
			*sol2 = (-b - sqrt(tmp)) / (2.0 * a);
		}
	}
	return GMRFLib_SUCCESS;
}

mapkit_size_t GMRFLib_nelm_map_ii(map_ii *hash)
{
	/*
	 * return the number of elements in HASH 
	 */

	mapkit_size_t i, nelm = 0;

	for (i = -1; (i = map_ii_next(hash, i)) != -1;) {
		nelm++;
	}
	return nelm;
}

mapkit_size_t GMRFLib_nelm_map_id(map_id *hash)
{
	/*
	 * return the number of elements in HASH 
	 */
	mapkit_size_t i, nelm = 0;

	for (i = -1; (i = map_id_next(hash, i)) != -1;) {
		nelm++;
	}
	return nelm;
}

map_ii *GMRFLib_duplicate_map_ii(map_ii *hash)
{
	/*
	 * return a copy of HASH 
	 */
	if (!hash) {
		return NULL;
	}

	mapkit_size_t i, nelm;
	map_ii *newhash = NULL;

	nelm = GMRFLib_nelm_map_ii(hash);
	newhash = Calloc(1, map_ii);
	map_ii_init_hint(newhash, nelm);

	for (i = -1; (i = map_ii_next(hash, i)) != -1;) {
		map_ii_set(newhash, hash->contents[i].key, hash->contents[i].value);
	}
	newhash->alwaysdefault = hash->alwaysdefault;
	newhash->defaultvalue = hash->defaultvalue;

	return newhash;
}

map_id *GMRFLib_duplicate_map_id(map_id *hash)
{
	/*
	 * return a copy of HASH 
	 */
	if (!hash) {
		return NULL;
	}

	mapkit_size_t i, nelm;
	map_id *newhash = NULL;

	nelm = GMRFLib_nelm_map_id(hash);
	newhash = Calloc(1, map_id);
	map_id_init_hint(newhash, nelm);

	for (i = -1; (i = map_id_next(hash, i)) != -1;) {
		map_id_set(newhash, hash->contents[i].key, hash->contents[i].value);
	}
	newhash->alwaysdefault = hash->alwaysdefault;
	newhash->defaultvalue = hash->defaultvalue;

	return newhash;
}

int GMRFLib_is_int(char *str, int *value)
{
	/*
	 * return 1 if an int can be read from STR and 0 if not. return value in *VALUE is non-NULL.
	 */
	if (!str) {
		return 0;
	}

	int x, err;
	err = (sscanf(str, "%d", &x) == 1 ? 1 : 0);
	if (value) {
		*value = x;
	}
	return err;
}

char *GMRFLib_strtok_r(char *s1, const char *s2, char **lasts)
{
	char *ret = NULL;

	if (*lasts == NULL && s1 == NULL) {		       /* added this: hrue */
		return NULL;
	}

	if (s1 == NULL) {
		s1 = *lasts;
	}
	while (*s1 && strchr(s2, *s1)) {
		++s1;
	}
	if (*s1 == '\0') {
		return NULL;
	}
	ret = s1;
	while (*s1 && !strchr(s2, *s1)) {
		++s1;
	}
	if (*s1) {
		*s1++ = '\0';
	}
	*lasts = s1;

	return ret;
}

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wattributes"
__attribute__((target_clones(INLA_CLONE_TARGETS "default")))
int GMRFLib_iuniques(int *nuniques, int **uniques, int *ix, int nx)
{
	/*
	 * return in number of unique entries in ix != 0 and list them in `uniques' 
	 */

	int nu, *un = NULL, i, j, *ixx = NULL;

	if (nx <= 0 || !ix) {
		*nuniques = 0;
		if (uniques) {
			*uniques = NULL;
		}
		return GMRFLib_SUCCESS;
	}

	ixx = Malloc(nx, int);
	Memcpy(ixx, ix, nx * sizeof(int));
	QSORT_FUN((void *) ixx, (size_t) nx, sizeof(int), GMRFLib_icmp);

	for (j = nu = i = 0; i < nx; i++) {
		if (ixx[i] && (!i || ixx[i] != ixx[j])) {
			nu++;
			j = i;
		}
	}
	un = Malloc(nu, int);

	for (j = nu = i = 0; i < nx; i++) {
		if (ixx[i] && (!i || ixx[i] != ixx[j])) {
			un[nu++] = ixx[i];
			j = i;
		}
	}

	*nuniques = nu;
	if (uniques) {
		*uniques = un;
	} else {
		Free(un);
	}
	Free(ixx);

	return GMRFLib_SUCCESS;
}
#pragma GCC diagnostic pop

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wattributes"
__attribute__((target_clones(INLA_CLONE_TARGETS "default")))
int GMRFLib_gsl_vec2plain(double **out, gsl_vector *vec)
{
	if (!vec || vec->size == 0) {
		*out = NULL;
	} else {
		*out = Malloc(vec->size, double);
		for (size_t i = 0; i < vec->size; i++) {
			(*out)[i] = gsl_vector_get(vec, i);
		}
	}
	return GMRFLib_SUCCESS;
}
#pragma GCC diagnostic pop

int GMRFLib_gsl_mat2plain(double **out, gsl_matrix *mat)
{
	if (!mat || mat->size1 == 0 || mat->size2 == 0) {
		*out = NULL;
	} else {
		*out = Malloc(mat->size1 * mat->size2, double);
		for (size_t j = 0; j < mat->size2; j++) {
			size_t off = j * mat->size1;
			for (size_t i = 0; i < mat->size1; i++) {
				(*out)[i + off] = gsl_matrix_get(mat, i, j);
			}
		}
	}
	return GMRFLib_SUCCESS;
}

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wattributes"
__attribute__((target_clones(INLA_CLONE_TARGETS "default")))
int GMRFLib_adjust_vector(double *x, int n)
{
	/*
	 * x := x - max(x[]) 
	 */
	double max_value;

	if (n <= 0 || !x) {
		return GMRFLib_SUCCESS;
	}

	max_value = GMRFLib_max_value(x, n, NULL);
#pragma omp simd
	for (int i = 0; i < n; i++) {
		x[i] -= max_value;
	}

	return GMRFLib_SUCCESS;
}
#pragma GCC diagnostic pop

int GMRFLib_scale_vector(double *x, int n)
{
	/*
	 * x := x/max(x)
	 */
	if (n <= 0 || !x) {
		return GMRFLib_SUCCESS;
	}

	double scale = GMRFLib_max_value(x, n, NULL);
	if (ISNONZERO(scale)) {
		scale = 1.0 / scale;
		GMRFLib_dscale(n, scale, x);
	}

	return GMRFLib_SUCCESS;
}

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wattributes"
__attribute__((target_clones(INLA_CLONE_TARGETS "default")))
int GMRFLib_is_zero(double *x, int n)
{
	// return 1 if x is a zero vector or zero-ptr, 0 otherwise
	if (!x) {
		return 1;
	}

	const int nstart = 32L;
	int m = IMIN(n, nstart);
	int nonzero = 0;

#pragma omp simd reduction(|:nonzero)
	for (int i = 0; i < m; i++) {
		nonzero |= ISNONZERO(x[i]);
	}
	if (nonzero)
		return 0;

	if (n > m) {
		for (int offset = nstart; offset < n; offset += offset) {
			int len = IMIN(offset, n - offset);
			if (len <= 0 || memcmp(x, x + offset, len * sizeof(double))) {
				return 0;
			}
		}
	}

	return 1;
}
#pragma GCC diagnostic pop

double GMRFLib_max_value(double *x, int n, int *idx)
{
	/*
	 * return the MAX(x[]), optional idx
	 */
	if (n <= 1) {
		if (n <= 0) {
			if (idx) {
				idx = NULL;
			}
			return NAN;
		} else {
			assert(x);
			if (idx) {
				*idx = 0;
			}
			return x[0];
		}
	}
	// sometimes the max is at the boundary
	assert(x);
	if (idx) {
		int imax;
		double max_val;
		if (x[0] > x[n - 1]) {
			imax = 0;
			max_val = x[0];
		} else {
			imax = n - 1;
			max_val = x[n - 1];
		}
		for (int i = 1; i < n - 1; i++) {
			if (unlikely(x[i] > max_val)) {
				max_val = x[i];
				imax = i;
			}
		}
		*idx = imax;
		return max_val;
	} else {
		double max_val = DMAX(x[0], x[n - 1]);
		for (int i = 1; i < n - 1; i++) {
			if (unlikely(x[i] > max_val)) {
				max_val = x[i];
			}
		}
		return max_val;
	}
}

double GMRFLib_min_value(double *x, int n, int *idx)
{
	/*
	 * return the MIN(x[]), optional idx
	 */

	if (n <= 1) {
		if (n <= 0) {
			if (idx) {
				idx = NULL;
			}
			return NAN;
		} else {
			assert(x);
			if (idx) {
				*idx = 0;
			}
			return x[0];
		}
	}
	// sometimes the min is at the boundary
	assert(x);
	if (idx) {
		int imin;
		double min_val;
		if (x[0] < x[n - 1]) {
			imin = 0;
			min_val = x[0];
		} else {
			imin = n - 1;
			min_val = x[n - 1];
		}
		for (int i = 1; i < n - 1; i++) {
			if (unlikely(x[i] < min_val)) {
				min_val = x[i];
				imin = i;
			}
		}
		*idx = imin;
		return min_val;
	} else {
		double min_val = DMIN(x[0], x[n - 1]);
		for (int i = 1; i < n - 1; i++) {
			if (unlikely(x[i] < min_val)) {
				min_val = x[i];
			}
		}
		return min_val;
	}
}

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wattributes"
__attribute__((target_clones(INLA_CLONE_TARGETS "default")))
int GMRFLib_imax_value(int *x, int n, int *idx)
{
	/*
	 * return the IMAX(x[]), optional idx
	 */

	if (n <= 1) {
		if (n <= 0) {
			if (idx) {
				idx = NULL;
			}
			return INT_MAX;
		} else {
			assert(x);
			if (idx) {
				*idx = 0;
			}
			return x[0];
		}
	}
	// sometimes the max is at the boundary, but we're just 'almost' sure
	assert(x);
	if (idx) {
		int imax;
		int max_val;
		if (x[0] > x[n - 1]) {
			imax = 0;
			max_val = x[0];
		} else {
			imax = n - 1;
			max_val = x[n - 1];
		}
		for (int i = 1; i < n - 1; i++) {
			if (unlikely(x[i] > max_val)) {
				max_val = x[i];
				imax = i;
			}
		}
		*idx = imax;
		return max_val;
	} else {
		int max_val = IMAX(x[0], x[n - 1]);
		for (int i = 1; i < n - 1; i++) {
			if (unlikely(x[i] > max_val)) {
				max_val = x[i];
			}
		}
		return max_val;
	}
}
#pragma GCC diagnostic pop

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wattributes"
__attribute__((target_clones(INLA_CLONE_TARGETS "default")))
int GMRFLib_imin_value(int *x, int n, int *idx)
{
	/*
	 * return the IMIN(x[]), optional idx
	 */

	if (n <= 1) {
		if (n <= 0) {
			if (idx) {
				idx = NULL;
			}
			return INT_MIN;
		} else {
			assert(x);
			if (idx) {
				*idx = 0;
			}
			return x[0];
		}
	}
	// sometimes the min is at the boundary, but we're just 'almost' sure
	assert(x);
	if (idx) {
		int imin;
		int min_val;
		if (x[0] < x[n - 1]) {
			imin = 0;
			min_val = x[0];
		} else {
			imin = n - 1;
			min_val = x[n - 1];
		}
		for (int i = 1; i < n - 1; i++) {
			if (unlikely(x[i] < min_val)) {
				min_val = x[i];
				imin = i;
			}
		}
		*idx = imin;
		return min_val;
	} else {
		int min_val = IMIN(x[0], x[n - 1]);
		for (int i = 1; i < n - 1; i++) {
			if (unlikely(x[i] < min_val)) {
				min_val = x[i];
			}
		}
		return min_val;
	}
}
#pragma GCC diagnostic pop

int GMRFLib_iamax_value(int *x, int n, int *idx)
{
	/*
	 * return IMAX(abs(x[])), optional idx
	 */
	int imax, max_val;

	max_val = IABS(x[0]);
	imax = 0;
	for (int i = 1; i < n; i++) {
		if (IABS(x[i]) > max_val) {
			max_val = IABS(x[i]);
			imax = i;
		}
	}

	if (idx) {
		*idx = imax;
	}
	return max_val;
}

double GMRFLib_logit(double p)
{
	// evaluate log(p/(1-p)) more safe than just log(p/(1-p))
	const double lim = 0.001;

	if (p > lim && p < 1.0 - lim) {
		return log(p / (1.0 - p));
	} else if (p < 0.5) {
		return (log(p) - log1p(-p));
	} else {
		double pp = 1.0 - p;
		return (-log(pp) + log1p(-pp));
	}
}

double GMRFLib_inv_logit(double x)
{
	// evaluate 1/(1+exp(-x))

	return 1.0 / (2.0 + expm1(-x));
}

const char *GMRFLib_function_name_strip(const char *name)
{
	char *s = (char *) name;
	if (!strncmp("GMRFLib_", s, 8)) {
		s += 8;
	}
	if (!strncmp("inla_", s, 5)) {
		s += 5;
	}
	return s;
}

int GMRFLib_debug_functions(const char *name)
{
	if (!name) {
		return 0;
	}

	static int not_defined = 0;
	if (not_defined) {
		return 0;
	}

	static map_stri **ddefs = NULL;
	static int *first = NULL;
	static int clen = 0;

	if (!ddefs) {
#pragma omp critical (Name_30c48b516c7b1cce1be137af0e429a5e3b52a645)
		if (!ddefs) {
			clen = GMRFLib_CACHE_LEN();
			first = Calloc(clen, int);
			map_stri **tmp = Calloc(clen, map_stri *);
			ddefs = tmp;
		}
	}
	int idx = 0;
	GMRFLib_CACHE_SET_IDX(idx);
	assert(idx < clen);

	if (!ddefs[idx]) {
#pragma omp critical (Name_c3afbb5a350a04cd0a2ad81d85df8cc44ff04279)
		if (!ddefs[idx]) {
			// format FUN[:N],...
			// prefix's GMRFLib_ and inla_ are removed automatically
			char *def = getenv("INLA_DEBUG");
			int verbose = 0;
			map_stri *tmp = NULL;

			if (def) {
				def = Strdup(def);
			}
			if (verbose) {
				printf("\t\tREAD %s\n", def);
			}

			if (!def) {
				not_defined = 1;
			} else {
				char sep1[] = ",;";

				tmp = Calloc(1, map_stri);
				map_stri_init_hint(tmp, 128);
				char *str = def;
				char *s = NULL;

				first[idx] = -1;
				while ((s = strtok(str, sep1))) {
					str = NULL;

					int val = 0;
					char *s2 = strchr(s, ':');
					char *ss = NULL;
					if (!s2) {
						ss = s;
						val = 1;
					} else {
						int len = s2 - s + 1;
						int len1 = len + 1;	/* to avoid compiler warning */
						assert(len >= 0);
						ss = Calloc(len + 1, char);
						ss[len1 - 1] = '\0';
						strncpy(ss, s, len - 1);
						val = atoi(s2 + 1);
						val = IMAX(val, 1);
					}
					// strip leading whitespace
					while (!strncmp(ss, " ", 1))
						ss++;
					// special option that override all others
					if (!strcmp(ss, "*")) {
						first[idx] = 2;
					}

					char *nm = NULL;
					if (strlen(ss)) {
						GMRFLib_sprintf(&nm, "%s", ss);
						map_stri_set(tmp, nm, val);
						GMRFLib_sprintf(&nm, "GMRFLib_%s", ss);
						map_stri_set(tmp, nm, val);
						GMRFLib_sprintf(&nm, "inla_%s", ss);
						map_stri_set(tmp, nm, val);
					}
					if (first[idx] != 2) {
						first[idx] = 0;
					}

					if (verbose) {
						printf("\t\t[%1d] debug init: ADD [%s]=%1d\n", omp_get_thread_num(), ss, val);
					}
				}
			}
			ddefs[idx] = tmp;
		}
	}

	if (!name || not_defined) {
		return 0;
	} else {
		int *p = map_stri_ptr(ddefs[idx], (char *) (first[idx] == 2 ? "*" : name));
		return (p ? *p : 0);
	}
}

int GMRFLib_trace_functions(const char *name)
{
	if (!name) {
		return 0;
	}

	static int not_defined = 0;
	if (not_defined) {
		return 0;
	}

	static map_stri **ddefs = NULL;
	static int *first = NULL;
	static int clen = 0;

	if (!ddefs) {
#pragma omp critical (Name_3a266edf254a33111bcf4ab49b3acc5833850a29)
		if (!ddefs) {
			clen = GMRFLib_CACHE_LEN();
			first = Calloc(clen, int);
			map_stri **tmp = Calloc(clen, map_stri *);
			ddefs = tmp;
		}
	}
	int idx = 0;
	GMRFLib_CACHE_SET_IDX(idx);
	assert(idx < clen);

	if (!ddefs[idx]) {
#pragma omp critical (Name_e9b04207643dde9dc8734f9ae0e41a3e03910f80)
		if (!ddefs[idx]) {
			// format FUN[:N],...
			// prefix's GMRFLib_ and inla_ are removed automatically
			char *def = getenv("INLA_TRACE");
			int verbose = 0;
			map_stri *tmp = NULL;

			if (def) {
				def = Strdup(def);
			}
			if (verbose) {
				printf("\t\tREAD %s\n", def);
			}

			if (!def) {
				not_defined = 1;
			} else {
				char sep1[] = ",;";

				tmp = Calloc(1, map_stri);
				map_stri_init_hint(tmp, 128);
				char *str = def;
				char *s = NULL;

				first[idx] = -1;
				while ((s = strtok(str, sep1))) {
					str = NULL;

					int val = 0;
					char *s2 = strchr(s, ':');
					char *ss = NULL;
					if (!s2) {
						ss = s;
						val = 1;
					} else {
						int len = s2 - s + 1;
						assert(len >= 0);
						ss = Calloc(len + 1, char);
						ss[len] = '\0';
						strncpy(ss, s, len - 1);
						val = atoi(s2 + 1);
						val = IMAX(val, 1);
					}
					// strip leading whitespace
					while (!strncmp(ss, " ", 1))
						ss++;
					// special option that override all others
					if (!strcmp(ss, "*")) {
						first[idx] = 2;
					}

					char *nm = NULL;
					if (strlen(ss)) {
						GMRFLib_sprintf(&nm, "%s", ss);
						map_stri_set(tmp, nm, val);
						GMRFLib_sprintf(&nm, "inla_%s", ss);
						map_stri_set(tmp, nm, val);
						GMRFLib_sprintf(&nm, "GMRFLib_%s", ss);
						map_stri_set(tmp, nm, val);
					}
					if (first[idx] != 2) {
						first[idx] = 0;
					}

					if (verbose) {
						printf("\t\t[%1d] debug init: ADD [%s]=%1d\n", omp_get_thread_num(), ss, val);
					}
				}
			}
			ddefs[idx] = tmp;
		}
	}

	if (!name || not_defined) {
		return 0;
	} else {
		int *p = map_stri_ptr(ddefs[idx], (char *) (first[idx] == 2 ? "*" : name));
		return (p ? *p : 0);
	}
}

int GMRFLib_trace_cache_hitmiss(const char *name)
{
	if (!name) {
		return 0;
	}

	static int not_defined = 0;
	if (not_defined) {
		return 0;
	}

	static map_stri **ddefs = NULL;
	static int *first = NULL;
	static int clen = 0;

	if (!ddefs) {
#pragma omp critical (Name_72150bb8d161e16549ba70e0a250eb5d4f572df6)
		if (!ddefs) {
			clen = GMRFLib_CACHE_LEN();
			first = Calloc(clen, int);
			map_stri **tmp = Calloc(clen, map_stri *);
			ddefs = tmp;
		}
	}
	int idx = 0;
	GMRFLib_CACHE_SET_IDX(idx);
	assert(idx < clen);

	if (!ddefs[idx]) {
#pragma omp critical (Name_4b0bcb2d4e2c1a81a1672358ca7320e389c962bc)
		if (!ddefs[idx]) {
			// format FUN[:N],...
			// prefix's GMRFLib_ and inla_ are removed automatically
			char *def = getenv("INLA_CACHE");
			int verbose = 0;
			map_stri *tmp = NULL;

			if (def) {
				def = Strdup(def);
			}
			if (verbose) {
				printf("\t\tREAD %s\n", def);
			}

			if (!def) {
				not_defined = 1;
			} else {
				char sep1[] = ",;";

				tmp = Calloc(1, map_stri);
				map_stri_init_hint(tmp, 128);
				char *str = def;
				char *s = NULL;

				first[idx] = -1;
				while ((s = strtok(str, sep1))) {
					str = NULL;

					int val = 0;
					char *s2 = strchr(s, ':');
					char *ss = NULL;
					if (!s2) {
						ss = s;
						val = 1;
					} else {
						int len = s2 - s + 1;
						assert(len >= 0);
						ss = Calloc(len + 1, char);
						ss[len] = '\0';
						strncpy(ss, s, len - 1);
						val = atoi(s2 + 1);
						val = IMAX(val, 1);
					}
					// strip leading whitespace
					while (!strncmp(ss, " ", 1))
						ss++;
					// special option that override all others
					if (!strcmp(ss, "*")) {
						first[idx] = 2;
					}

					char *nm = NULL;
					if (strlen(ss)) {
						GMRFLib_sprintf(&nm, "%s", ss);
						map_stri_set(tmp, nm, val);
						GMRFLib_sprintf(&nm, "inla_%s", ss);
						map_stri_set(tmp, nm, val);
						GMRFLib_sprintf(&nm, "GMRFLib_%s", ss);
						map_stri_set(tmp, nm, val);
					}
					if (first[idx] != 2) {
						first[idx] = 0;
					}

					if (verbose) {
						printf("\t\t[%1d] debug init: ADD [%s]=%1d\n", omp_get_thread_num(), ss, val);
					}
				}
			}
			ddefs[idx] = tmp;
		}
	}

	if (!name || not_defined) {
		return 0;
	} else {
		int *p = map_stri_ptr(ddefs[idx], (char *) (first[idx] == 2 ? "*" : name));
		return (p ? *p : 0);
	}
}

// ******************************************************************************************

int GMRFLib_vmatrix_init(GMRFLib_vmatrix_tp **vmatrix, int nrow, GMRFLib_graph_tp *graph)
{
	// graph is optional. If given, the initialise with lnnbs+1

	*vmatrix = Calloc(1, GMRFLib_vmatrix_tp);
	(*vmatrix)->nrow = nrow;
	(*vmatrix)->vmat = Calloc(nrow, map_ivp);
	if (graph) {
		for (int i = 0; i < nrow; i++) {
			map_ivp_init_hint(&((*vmatrix)->vmat[i]), (mapkit_size_t) (graph->lnnbs[i] + 1));
		}
	} else {
		for (int i = 0; i < nrow; i++) {
			map_ivp_init(&((*vmatrix)->vmat[i]));
		}
	}

	return GMRFLib_SUCCESS;
}

int GMRFLib_vmatrix_set(GMRFLib_vmatrix_tp *vmatrix, int i, int j, double *vec)
{
	map_ivp_set(&(vmatrix->vmat[i]), j, (void *) vec);
	return GMRFLib_SUCCESS;
}

double *GMRFLib_vmatrix_get(GMRFLib_vmatrix_tp *vmatrix, int i, int j)
{
	void *p = NULL;
	map_ivp_get(&(vmatrix->vmat[i]), j, &p);
	return ((double *) p);
}

int GMRFLib_vmatrix_free(GMRFLib_vmatrix_tp *vmatrix, int free_content)
{
	if (free_content) {
		for (int i = 0; i < vmatrix->nrow; i++) {
			for (int j = -1; (j = map_ivp_next(&(vmatrix->vmat[i]), j)) != -1;) {
				Free(vmatrix->vmat[i].contents[j].value);
			}
		}
	}

	for (int i = 0; i < vmatrix->nrow; i++) {
		map_ivp_free((map_ivp *) & (vmatrix->vmat[i]));
	}

	Free(vmatrix->vmat);
	Free(vmatrix);

	return GMRFLib_SUCCESS;
}

// ****************************************************************************************

/*
 * Implement Heap sort -- direct and indirect sorting
 * Based on descriptions in Sedgewick "Algorithms in C"
 *
 * Copyright (C) 1999  Thomas Walter
 */

void my_downheap2_id(int *__restrict data1, double *__restrict data2, const int N, int k)
{
	int v1 = data1[k];
	double v2 = data2[k];

	while (k <= N / 2) {
		int j = 2 * k;
		if (j < N && data1[j] < data1[j + 1]) {
			j++;
		}

		if (!(v1 < data1[j])) {
			break;
		}

		data1[k] = data1[j];
		data2[k] = data2[j];
		k = j;
	}
	data1[k] = v1;
	data2[k] = v2;
}

void gsl_sort2_id(int *__restrict data1, double *__restrict data2, const int n)
{
	int N, k;

	if (n == 0) {
		return;					       /* No data to sort */
	}

	/*
	 * We have n_data elements, last element is at 'n_data-1', first at '0' Set N to the last element number. 
	 */

	N = n - 1;
	k = N / 2;
	k++;						       /* Compensate the first use of 'k--' */
	do {
		k--;
		my_downheap2_id(data1, data2, N, k);
	} while (k > 0);

	while (N > 0) {
		int tmp1 = data1[0];
		data1[0] = data1[N];
		data1[N] = tmp1;

		double tmp2 = data2[0];
		data2[0] = data2[N];
		data2[N] = tmp2;

		/*
		 * then process the heap 
		 */
		N--;
		my_downheap2_id(data1, data2, N, 0);
	}
}

void my_downheap2_ii(int *__restrict data1, int *__restrict data2, const int N, int k)
{
	int v1 = data1[k];
	int v2 = data2[k];

	while (k <= N / 2) {
		int j = 2 * k;
		if (j < N && data1[j] < data1[j + 1]) {
			j++;
		}

		if (!(v1 < data1[j])) {
			break;
		}

		data1[k] = data1[j];
		data2[k] = data2[j];
		k = j;
	}
	data1[k] = v1;
	data2[k] = v2;
}

void gsl_sort2_ii(int *__restrict data1, int *__restrict data2, const int n)
{
	int N, k;

	if (n == 0) {
		return;					       /* No data to sort */
	}

	/*
	 * We have n_data elements, last element is at 'n_data-1', first at '0' Set N to the last element number. 
	 */

	N = n - 1;
	k = N / 2;
	k++;						       /* Compensate the first use of 'k--' */
	do {
		k--;
		my_downheap2_ii(data1, data2, N, k);
	} while (k > 0);

	while (N > 0) {
		int tmp1 = data1[0];
		data1[0] = data1[N];
		data1[N] = tmp1;

		int tmp2 = data2[0];
		data2[0] = data2[N];
		data2[N] = tmp2;

		/*
		 * then process the heap 
		 */
		N--;
		my_downheap2_ii(data1, data2, N, 0);
	}
}

void my_insertionSort_id(int *__restrict iarr, double *__restrict darr, int n)
{
	if (darr) {
		for (int i = 1; i < n; i++) {
			int key = iarr[i];
			double dkey = darr[i];
			int j = i - 1;
			while (j >= 0 && iarr[j] > key) {
				iarr[j + 1] = iarr[j];
				darr[j + 1] = darr[j];
				j--;
			}
			iarr[j + 1] = key;
			darr[j + 1] = dkey;
		}
	} else {
		for (int i = 1; i < n; i++) {
			int key = iarr[i];
			int j = i - 1;
			while (j >= 0 && iarr[j] > key) {
				iarr[j + 1] = iarr[j];
				j--;
			}
			iarr[j + 1] = key;
		}
	}
}

void my_insertionSort_ii(int *__restrict iarr, int *__restrict darr, int n)
{
	if (darr) {
		for (int i = 1; i < n; i++) {
			int key = iarr[i];
			int dkey = darr[i];
			int j = i - 1;
			while (j >= 0 && iarr[j] > key) {
				iarr[j + 1] = iarr[j];
				darr[j + 1] = darr[j];
				j--;
			}
			iarr[j + 1] = key;
			darr[j + 1] = dkey;
		}
	} else {
		for (int i = 1; i < n; i++) {
			int key = iarr[i];
			int j = i - 1;
			while (j >= 0 && iarr[j] > key) {
				iarr[j + 1] = iarr[j];
				j--;
			}
			iarr[j + 1] = key;
		}
	}
}

void my_insertionSort_dd(double *__restrict iarr, double *__restrict darr, int n)
{
	if (darr) {
		for (int i = 1; i < n; i++) {
			double key = iarr[i];
			double dkey = darr[i];
			int j = i - 1;
			while (j >= 0 && iarr[j] > key) {
				iarr[j + 1] = iarr[j];
				darr[j + 1] = darr[j];
				j--;
			}
			iarr[j + 1] = key;
			darr[j + 1] = dkey;
		}
	} else {
		for (int i = 1; i < n; i++) {
			double key = iarr[i];
			int j = i - 1;
			while (j >= 0 && iarr[j] > key) {
				iarr[j + 1] = iarr[j];
				j--;
			}
			iarr[j + 1] = key;
		}
	}
}

void my_insertionSort_i(int *__restrict iarr, int n)
{
	for (int i = 1; i < n; i++) {
		int key = iarr[i];
		int j = i - 1;
		while (j >= 0 && iarr[j] > key) {
			iarr[j + 1] = iarr[j];
			j--;
		}
		iarr[j + 1] = key;
	}
}

void my_insertionSort_d(double *__restrict iarr, int n)
{
	for (int i = 1; i < n; i++) {
		double key = iarr[i];
		int j = i - 1;
		while (j >= 0 && iarr[j] > key) {
			iarr[j + 1] = iarr[j];
			j--;
		}
		iarr[j + 1] = key;
	}
}

void gsl_sort2_dd(double *__restrict data1, double *__restrict data2, const int n)
{
	gsl_sort2(data1, (size_t) 1, data2, (size_t) 1, (size_t) n);
}

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wattributes"
__attribute__((target_clones(INLA_CLONE_TARGETS "default")))
void my_sort2_ii(int *__restrict ix, int *__restrict x, int n)
{
	if (n <= 1 || GMRFLib_is_sorted_iinc(n, ix))
		return;

	if (1) {
		// this one is now a better option (feb'2024). no need to initialize with 0's
		int *ixy = Malloc(n * 2, int);
#pragma omp simd
		for (int i = 0; i < n; i++) {
			int j = 2 * i;
			ixy[j] = ix[i];
			ixy[j + 1] = x[i];
		}
		QSORT_FUN((void *) ixy, (size_t) n, 2 * sizeof(int), GMRFLib_icmp);
#pragma omp simd
		for (int i = 0; i < n; i++) {
			int j = 2 * i;
			ix[i] = ixy[j];
			x[i] = ixy[j + 1];
		}
		Free(ixy);
		return;
	} else {
		if (n < GMRFLib_sort2_id_cut_off) {
			my_insertionSort_ii(ix, x, n);
		} else {
			gsl_sort2_ii(ix, x, n);
		}
	}
}
#pragma GCC diagnostic pop

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wattributes"
__attribute__((target_clones(INLA_CLONE_TARGETS "default")))
void my_sort2_id_work(int *__restrict ix, double *__restrict x, int n, double *work)
{
	// this does not go that well: see test 160

	for (int i = 0; i < n; i++) {
		int j = 2 * i;
		int *ip = (int *) (work + j);
		double *dp = (work + j + 1);
		*ip = ix[i];
		*dp = x[i];
	}

	QSORT_FUN((void *) work, (size_t) n, (size_t) (2 * sizeof(double)), GMRFLib_icmp);

	for (int i = 0; i < n; i++) {
		int j = 2 * i;
		int *ip = (int *) (work + j);
		double *dp = (work + j + 1);
		ix[i] = *ip;
		x[i] = *dp;
	}
}
#pragma GCC diagnostic pop

void my_sort2_id(int *__restrict ix, double *__restrict x, int n)
{
	return my_sort2_id_x(ix, x, n, NULL);
}

void my_sort2_id_x(int *__restrict ix, double *__restrict x, int n, void *UNUSED(work))
{
	if (n <= 1)
		return;

	// do not precheck if using insertionSort
	if (__builtin_expect(n < GMRFLib_sort2_id_cut_off, 0)) {
		my_insertionSort_id(ix, x, n);
	} else if (__builtin_expect(GMRFLib_is_sorted_iinc(n, ix), 1)) {
		// already sorted
	} else {
		gsl_sort2_id(ix, x, n);
	}

	return;
}

void my_sort2_dd(double *__restrict ix, double *__restrict x, int n)
{
	if (n <= 1)
		return;

	// do not precheck if using insertionSort
	if (__builtin_expect(n < GMRFLib_sort2_dd_cut_off, 0)) {
		my_insertionSort_dd(ix, x, n);
	} else if (__builtin_expect(GMRFLib_is_sorted_dinc(n, ix), 1)) {
		// already sorted
	} else {
		gsl_sort2_dd(ix, x, n);
	}
}

int my_sort2_id_test_cutoff(int verbose)
{
	const int nmax = 512;
	const int nmin = 64;
	const int nstep = 64;
	const int ntimes = 200;

	double time_used = 0.0;
	int *ix = Calloc(2 * nmax, int);
	double *x = Calloc(2 * nmax, double);

	double slope_xy = 0.0;
	double slope_xx = 0.0;
	double slope_x = 0.0;
	double slope_y = 0.0;
	double slope_n = 0.0;
	double cutoff = 1;
	double b;

	time_used -= GMRFLib_timer();

	for (int n = nmin; n <= nmax; n += nstep) {

		int *ixx = ix + nmax;
		double *xx = x + nmax;
		double time[2] = { 0.0, 0.0 };

		for (int times = 0; times < ntimes; times++) {

			for (int i = 0; i < n; i++) {
				ix[i] = (int) ((100 * nmax) * GMRFLib_uniform());
				x[i] = GMRFLib_uniform();
			}

			Memcpy(ixx, ix, n * sizeof(int));
			Memcpy(xx, x, n * sizeof(double));
			time[0] -= GMRFLib_timer();
			my_insertionSort_id(ixx, xx, n);
			time[0] += GMRFLib_timer();

			Memcpy(ixx, ix, n * sizeof(int));
			Memcpy(xx, x, n * sizeof(double));
			time[1] -= GMRFLib_timer();
			gsl_sort2_id(ixx, xx, n);
			time[1] += GMRFLib_timer();
		}

		slope_xx += SQR(n);
		slope_xy += n * (time[0] / time[1]);
		slope_x += n;
		slope_y += (time[0] / time[1]);
		slope_n++;

		b = (slope_xy / slope_n - (slope_x / slope_n) * (slope_y / slope_n)) / (slope_xx / slope_n - SQR(slope_x / slope_n));
		if (ISZERO(b))
			b = 1.0;
		cutoff = (slope_x / slope_n) + (1.0 - (slope_y / slope_n)) / b;

		if (verbose) {
			printf("sort-test n = %1d  time(insertSort/gsl_sort2) =  %.2f cutoff.est = %1d\n", n, time[0] / time[1], (int) cutoff);
		}
	}

	// this is a global variable
	GMRFLib_sort2_id_cut_off = IMAX(nmin, IMIN(nmax, (int) cutoff));

	time_used += GMRFLib_timer();
	if (verbose) {
		printf("sort-test took %.4f seconds\n", time_used);
	}

	Free(ix);
	Free(x);

	return GMRFLib_sort2_id_cut_off;
}

int my_sort2_dd_test_cutoff(int verbose)
{
	const int nmax = 448;
	const int nmin = 64;
	const int nstep = 64;
	const int ntimes = 100;

	double time_used = 0.0;
	double *ix = Calloc(2 * nmax, double);
	double *x = Calloc(2 * nmax, double);

	double slope_xy = 0.0;
	double slope_xx = 0.0;
	double slope_x = 0.0;
	double slope_y = 0.0;
	double slope_n = 0.0;
	double cutoff = 1;
	double b;

	time_used -= GMRFLib_timer();

	for (int n = nmin; n <= nmax; n += nstep) {

		double *ixx = ix + nmax;
		double *xx = x + nmax;
		double time[2] = { 0.0, 0.0 };

		for (int times = 0; times < ntimes; times++) {

			for (int i = 0; i < n; i++) {
				ix[i] = GMRFLib_uniform();
				x[i] = GMRFLib_uniform();
			}

			Memcpy(ixx, ix, n * sizeof(double));
			Memcpy(xx, x, n * sizeof(double));
			time[0] -= GMRFLib_timer();
			my_insertionSort_dd(ixx, xx, n);
			time[0] += GMRFLib_timer();

			Memcpy(ixx, ix, n * sizeof(double));
			Memcpy(xx, x, n * sizeof(double));
			time[1] -= GMRFLib_timer();
			gsl_sort2_dd(ixx, xx, n);
			time[1] += GMRFLib_timer();
		}

		slope_xx += SQR(n);
		slope_xy += n * (time[0] / time[1]);
		slope_x += n;
		slope_y += (time[0] / time[1]);
		slope_n++;

		b = (slope_xy / slope_n - (slope_x / slope_n) * (slope_y / slope_n)) / (slope_xx / slope_n - SQR(slope_x / slope_n));
		if (ISZERO(b))
			b = 1.0;
		cutoff = (slope_x / slope_n) + (1.0 - (slope_y / slope_n)) / b;

		if (verbose) {
			printf("sort-test n = %1d  time(insertSort/gsl_sort2) =  %.2f cutoff.est = %1d\n", n, time[0] / time[1], (int) cutoff);
		}
	}

	// this is a global variable
	GMRFLib_sort2_dd_cut_off = IMAX(nmin, IMIN(nmax, (int) cutoff));

	time_used += GMRFLib_timer();
	if (verbose) {
		printf("sort-test took %.4f seconds\n", time_used);
	}

	Free(ix);
	Free(x);

	return GMRFLib_sort2_dd_cut_off;
}

#define SOURCE_INCLUDE(CMP)			\
	for (int i = 0; i < n - 1; i++)		\
		if (a[i + 1] CMP a[i])		\
			return 0;		\
	return 1

int GMRFLib_is_sorted_iinc(int n, int *a)
{
	// increasing int's
	SOURCE_INCLUDE(<);
}

int GMRFLib_is_sorted_dinc(int n, double *a)
{
	// increasing double's
	SOURCE_INCLUDE(<);
}

int GMRFLib_is_sorted_idec(int n, int *a)
{
	// decreasing int's
	SOURCE_INCLUDE(>);
}

int GMRFLib_is_sorted_ddec(int n, double *a)
{
	// decreasing double's
	SOURCE_INCLUDE(>);
}

int GMRFLib_is_sorted_iinc_plain(int n, int *a)
{
	SOURCE_INCLUDE(<);
}

int GMRFLib_is_sorted_dinc_plain(int n, double *a)
{
	SOURCE_INCLUDE(<);
}

int GMRFLib_is_sorted_idec_plain(int n, int *a)
{
	SOURCE_INCLUDE(>);
}

int GMRFLib_is_sorted_ddec_plain(int n, double *a)
{
	SOURCE_INCLUDE(>);
}

#undef SOURCE_INCLUDE

int GMRFLib_is_sorted(void *a, size_t n, size_t size, int (*cmp)(const void *, const void *))
{
	if ( (cmp == (void *) GMRFLib_icmp) && size == sizeof(int)) {
		// increasing ints
		return GMRFLib_is_sorted_iinc(n, (int *) a);
	} else if (cmp == (void *) GMRFLib_icmp_r && size == sizeof(int)) {
		// decreasing ints
		return GMRFLib_is_sorted_idec(n, (int *) a);
	} else if (cmp == (void *) GMRFLib_dcmp && size == sizeof(double)) {
		// increasing doubles
		return GMRFLib_is_sorted_dinc(n, (double *) a);
	} else if (cmp == (void *) GMRFLib_dcmp_r && size == sizeof(double)) {
		// decreasing doubles
		return GMRFLib_is_sorted_ddec(n, (double *) a);
	} else {
		// by default not sorted
		return 0;
	}
	return 0;
}

void GMRFLib_qsort(void *a, size_t n, size_t size, int (*cmp)(const void *, const void *))
{
	// sort if not sorted
	if (n > 0 && !GMRFLib_is_sorted(a, n, size, cmp)) {
		QSORT_FUN(a, n, size, cmp);
	}
}

void GMRFLib_qsort2(void *x, size_t nmemb, size_t size_x, void *y, size_t size_y, int (*compar)(const void *, const void *))
{
	if (!y) {
		return(GMRFLib_qsort(x, nmemb, size_x, compar));
	}

	if (nmemb == 0) {
		return;
	}
	// there could be a test for GMRFLib_icmp_r but since I do not use it, I do not include it here
	if (compar == GMRFLib_icmp && size_x == size_y && size_x == sizeof(int)) {
		my_sort2_ii((int *) x, (int *) y, (int) nmemb);
		return;
	}

	if (compar == GMRFLib_dcmp && size_x == size_y && size_x == sizeof(double)) {
		my_sort2_dd((double *) x, (double *) y, (int) nmemb);
		return;
	}

	size_t siz = size_x + size_y;
	char *xy = Calloc(nmemb * siz, char);
	char *xx = (char *) x;
	char *yy = (char *) y;

	for (size_t i = 0, offset = 0; i < nmemb; i++) {
		Memcpy((void *) &xy[offset], (void *) &xx[i * size_x], size_x);
		offset += size_x;
		Memcpy((void *) &xy[offset], (void *) &yy[i * size_y], size_y);
		offset += size_y;
	}
	qsort((void *) xy, nmemb, siz, compar);
	for (size_t i = 0, offset = 0; i < nmemb; i++) {
		Memcpy((void *) &xx[i * size_x], (void *) &xy[offset], size_x);
		offset += size_x;
		Memcpy((void *) &yy[i * size_y], (void *) &xy[offset], size_y);
		offset += size_y;
	}
	Free(xy);
}


// easier interface to sort ints and doubles, increasingly
void GMRFLib_sort_i(int *ix, int n)
{
	return QSORT_FUN((void *) ix, (size_t) n, sizeof(int), GMRFLib_icmp);
}
void GMRFLib_sort_d(double *x, int n)
{
	if (n <= 32L) {
		my_insertionSort_d(x, n);
	} else {
		QSORT_FUN((void *) x, (size_t) n, sizeof(double), GMRFLib_dcmp);
	}
}


// 
double GMRFLib_cdfnorm_inv(double p)
{
	return (gsl_cdf_ugaussian_Pinv(p));
#if 0
	// https://arxiv.org/abs/0901.0638
	int sign = (p < 0.5 ? -1 : 1);
	double u = DMAX(p, 1.0 - p);
	double v = -log(2.0 * (1.0 - u));
	double P = 1.2533141359896652729 +
	    v * (3.0333178251950406994 +
		 v * (2.3884158540184385711 +
		      v * (0.73176759583280610539 +
			   v * (0.085838533424158257377 +
				v * (0.0034424140686962222423 + (0.000036313870818023761224 + 4.3304513840364031401e-8 * v) * v)))));
	double Q = 1 + v * (2.9202373175993672857 +
			    v * (2.9373357991677046357 +
				 v * (1.2356513216582148689 +
				      v * (0.2168237095066675527 +
					   v * (0.014494272424798068406 + (0.00030617264753008793976 + 1.3141263119543315917e-6 * v) * v)))));
	return (sign * v * P / Q);
#endif
}

double GMRFLib_cdfnorm(double x)
{
	return (0.5 * (1.0 + GMRFLib_erf(M_SQRT1_2 * x)));
}

double GMRFLib_erf(double x)
{
	return erf(x);
}

double GMRFLib_erfc(double x)
{
	return erfc(x);
}

double GMRFLib_erf_inv(double x)
{
	return (M_SQRT1_2 * GMRFLib_cdfnorm_inv((x + 1.0) * 0.5));
}

double GMRFLib_erfc_inv(double x)
{
	return (M_SQRT1_2 * GMRFLib_cdfnorm_inv(1.0 - x * 0.5));
}


void GMRFLib_sys_cache(GMRFLib_sys_cache_tp *l123)
{
	if (l123) {
		Memset(l123, 0, sizeof(GMRFLib_sys_cache_tp));
	}

	static GMRFLib_sys_cache_tp L123;
	static int first = 1;

	if (first) {
#pragma omp critical (Name_baeb0f7ca4f3dbd67a64a855c0a33b451124c7aa)
		if (first) {
#if defined(__linux__)
			long tmp;
			tmp = sysconf(_SC_LEVEL1_DCACHE_SIZE);
			L123.l1_data = (size_t) (tmp >= 0 ? tmp : 0);
			tmp = sysconf(_SC_LEVEL1_ICACHE_SIZE);
			L123.l1_inst = (size_t) (tmp >= 0 ? tmp : 0);
			tmp = sysconf(_SC_LEVEL2_CACHE_SIZE);
			L123.l2 = (size_t) (tmp >= 0 ? tmp : 0);
			tmp = sysconf(_SC_LEVEL3_CACHE_SIZE);
			L123.l3 = (size_t) (tmp >= 0 ? tmp : 0);
#elif defined(__APPLE__)
			size_t len = sizeof(size_t);
			sysctlbyname("hw.l1dcachesize", &(L123.l1_data), &len, NULL, 0);
			sysctlbyname("hw.l1icachesize", &(L123.l1_inst), &len, NULL, 0);
			sysctlbyname("hw.l2cachesize", &(L123.l2), &len, NULL, 0);
			sysctlbyname("hw.l3cachesize", &(L123.l3), &len, NULL, 0);
#elif defined(_MSC_VER) || defined(__MINGW32__) || defined (__MSVCRT__)
			DWORD len = 0;
			GetLogicalProcessorInformation(NULL, &len);
			SYSTEM_LOGICAL_PROCESSOR_INFORMATION *info = (SYSTEM_LOGICAL_PROCESSOR_INFORMATION *) malloc(len);
			if (!info) {
				goto label_end;
			}
			if (!GetLogicalProcessorInformation(info, &len)) {
				free(info);
				goto label_end;
			}
			for (DWORD i = 0; i < len / sizeof(SYSTEM_LOGICAL_PROCESSOR_INFORMATION); i++) {
				if (info[i].Relationship == RelationCache) {
					CACHE_DESCRIPTOR cache = info[i].Cache;
					if (cache.Level == 1 && cache.Type == CacheData) {
						L123.l1_data = (size_t) cache.Size;
					} else if (cache.Level == 1 && cache.Type == CacheInstruction) {
						L123.l1_inst = (size_t) cache.Size;
					} else if (cache.Level == 2) {
						L123.l2 = (size_t) cache.Size;
					} else if (cache.Level == 3) {
						L123.l3 = (size_t) cache.Size;
					}
				}
			}
		      label_end:
#else
			Memset(&L123, 0, sizeof(GMRFLib_sys_cache_tp));
#endif
			first = 0;
		}
	}

	if (l123) {
		Memcpy(l123, &L123, sizeof(GMRFLib_sys_cache_tp));
	}
}

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wattributes"
__attribute__((optimize("O3")))
    __attribute__((target_clones(INLA_CLONE_TARGETS "default")))
void GMRFLib_zero_small(int n, double eps, double *x)
{
	// if (ABS(x[i]) < eps) x[i]=0.0
#if defined(INLA_WITH_SIMDE_AVX512F_) && defined(__AVX512F__)
#       include "intrinsics/simde/zero-small-avx512f.h"
#elif defined(INLA_WITH_SIMDE_AVX2_) && (!defined(__x86_64__) || (defined(__x86_64__) && defined(__AVX2__)))
#       include "intrinsics/simde/zero-small-avx2.h"
#elif defined(INLA_WITH_SIMDE)
#       include "intrinsics/simde/zero-small-sse2.h"
#else
#       pragma omp simd
	for (int i = 0; i < n; i++) {
		if (ABS(x[i]) < eps) {
			x[i] = 0.0;
		}
	}
#endif
}
#pragma GCC diagnostic pop
