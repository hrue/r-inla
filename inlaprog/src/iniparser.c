
/**
   @file    iniparser.c
   @author  N. Devillard
   @date    Mar 2000
   @version $Revision: 1.26 $
   @brief   Parser for ini files.
*/

/*
  Partly rewritten by hrue

    $Id: iniparser.c,v 1.26 2009/03/30 16:47:04 hrue Exp $
    $Author: hrue $
    $Date: 2009/03/30 16:47:04 $
    $Revision: 1.26 $
*/

#include <stddef.h>
#include <string.h>
#include <stdio.h>
#if !defined(__FreeBSD__)
#include <malloc.h>
#endif
#include <stdlib.h>
#include "iniparser.h"
#include "gsl/gsl_math.h"
#include "my-fix.h"
#include "strlib.h"

static const char RCSId[] = "file: " __FILE__ "  " HGVERSION;
#define INI_INVALID_KEY     ((char*)-1)
#define MY_STRING_LOWERCASE(a) my_strlwc(a)

static void iniparser_add_entry(dictionary * d, char *sec, char *key, char *val)
{
	char *longkey = NULL;

	/*
	 * Make a key as section:keyword 
	 */
	longkey = (char *) calloc((size_t) ((sec ? strlen(sec) : 0) + (key ? strlen(key) : 0) + LEN_INIPARSER_SEP + 1), (size_t) 1);
	if (key != NULL) {
		sprintf(longkey, "%s%c%s", sec, INIPARSER_SEP, key);
	} else {
		strcpy(longkey, sec);
	}

	/*
	 * Add (key,val) to dictionary 
	 */
	dictionary_set(d, longkey, val);

	free(longkey);
	return;
}

/**
  @brief    Get number of sections in a dictionary
  @param    d   Dictionary to examine
  @return   int Number of sections found in dictionary

  This function returns the number of sections found in a dictionary.
  The test to recognize sections is done on the string stored in the
  dictionary: a section name is given as "section" whereas a key is
  stored as "section:key", thus the test looks for entries that do not
  contain a colon.

  This clearly fails in the case a section name contains a colon, but
  this should simply be avoided.

  This function returns -1 in case of error.
 */
int iniparser_getnsec(dictionary * d)
{
	int i;
	int nsec;

	if (d == NULL)
		return -1;
	nsec = 0;
	for (i = 0; i < d->size; i++) {
		if (d->key[i] == NULL)
			continue;
		if (strchr(d->key[i], INIPARSER_SEP) == NULL) {
			nsec++;
		}
	}
	return nsec;
}

/**
  @brief    Get name for section n in a dictionary.
  @param    d   Dictionary to examine
  @param    n   Section number (from 0 to nsec-1).
  @return   Pointer to char string

  This function locates the n-th section in a dictionary and returns
  its name as a pointer to a string statically allocated inside the
  dictionary. Do not free or modify the returned string!

  This function returns NULL in case of error.
 */
char *iniparser_getsecname(dictionary * d, int n)
{
	int i;
	int foundsec;

	if (d == NULL || n < 0)
		return NULL;
	foundsec = 0;
	for (i = 0; i < d->size; i++) {
		if (d->key[i] == NULL)
			continue;
		if (strchr(d->key[i], INIPARSER_SEP) == NULL) {
			foundsec++;
			if (foundsec > n)
				break;
		}
	}
	if (foundsec <= n) {
		return NULL;
	}
	return d->key[i];
}

/**
  @brief    Dump a dictionary to an opened file pointer.
  @param    d   Dictionary to dump.
  @param    f   Opened file pointer to dump to.
  @return   void

  This function prints out the contents of a dictionary, one element by
  line, onto the provided file pointer. It is OK to specify @c stderr
  or @c stdout as output files. This function is meant for debugging
  purposes mostly.
 */
void iniparser_dump(dictionary * d, FILE * f)
{
	int i;

	if (d == NULL || f == NULL)
		return;
	for (i = 0; i < d->size; i++) {
		if (d->key[i] == NULL)
			continue;
		if (d->val[i] != NULL) {
			fprintf(f, "[%s]=[%s]\n", d->key[i], d->val[i]);
		} else {
			fprintf(f, "[%s]=NULL\n", d->key[i]);
		}
	}
	return;
}

/**
  @brief    Save a dictionary to a loadable ini file
  @param    d   Dictionary to dump
  @param    f   Opened file pointer to dump to
  @return   void

  This function dumps a given dictionary into a loadable ini file.
  It is Ok to specify @c stderr or @c stdout as output files.
 */
void iniparser_dump_ini(dictionary * d, FILE * f)
{
	int i, j;
	char *keym = NULL;
	int nsec;
	char *secname = NULL;
	int seclen;

	if (d == NULL || f == NULL)
		return;

	nsec = iniparser_getnsec(d);
	if (nsec < 1) {
		/*
		 * No section in file: dump all keys as they are 
		 */
		for (i = 0; i < d->size; i++) {
			if (d->key[i] == NULL)
				continue;
			fprintf(f, "%s = %s\n", d->key[i], d->val[i]);
		}
		return;
	}
	for (i = 0; i < nsec; i++) {
		secname = iniparser_getsecname(d, i);
		seclen = (int) strlen(secname);
		fprintf(f, "\n[%s]\n", secname);

		keym = (char *) calloc(strlen(secname) + LEN_INIPARSER_SEP + 1, (size_t) 1);
		sprintf(keym, "%s%c", secname, INIPARSER_SEP);
		for (j = 0; j < d->size; j++) {
			if (d->key[j] == NULL)
				continue;
			if (!strncmp(d->key[j], keym, (size_t) (seclen + 1))) {
				fprintf(f, "%-30s = %s\n", d->key[j] + seclen + 1, d->val[j] ? d->val[j] : "");
			}
		}
		free(keym);
	}
	fprintf(f, "\n");
	return;
}

/**
  @brief	Get the string associated to a key, return NULL if not found
  @param    d   Dictionary to search
  @param    key Key string to look for
  @return   pointer to statically allocated character string, or NULL.

  This function queries a dictionary for a key. A key as read from an
  ini file is given as "section:key". If the key cannot be found,
  NULL is returned.
  The returned char pointer is pointing to a string allocated in
  the dictionary, do not free or modify it.

  This function is only provided for backwards compatibility with 
  previous versions of iniparser. It is recommended to use
  iniparser_getstring() instead.
 */
char *iniparser_getstr(dictionary * d, const char *key)
{
	return iniparser_getstring(d, key, NULL);
}

/**
  @brief    Get the string associated to a key
  @param    d       Dictionary to search
  @param    key     Key string to look for
  @param    def     Default value to return if key not found.
  @return   pointer to statically allocated character string

  This function queries a dictionary for a key. A key as read from an
  ini file is given as "section:key". If the key cannot be found,
  the pointer passed as 'def' is returned.
  The returned char pointer is pointing to a string allocated in
  the dictionary, do not free or modify it.
 */
char *iniparser_getstring(dictionary * d, const char *key, char *def)
{
	char *lc_key = NULL;
	char *sval = NULL;

	if (d == NULL || key == NULL)
		return def;

	if (!(lc_key = strdup(MY_STRING_LOWERCASE(key)))) {
		return NULL;
	}
	sval = dictionary_get(d, lc_key, def);
	free(lc_key);
	return sval;
}

/**
  @brief    Get the string associated to a key, convert to an int
  @param    d Dictionary to search
  @param    key Key string to look for
  @param    notfound Value to return in case of error
  @return   integer

  This function queries a dictionary for a key. A key as read from an
  ini file is given as "section:key". If the key cannot be found,
  the notfound value is returned.

 */
int iniparser_getint(dictionary * d, const char *key, int notfound)
{
	char *str = NULL;

	str = iniparser_getstring(d, key, INI_INVALID_KEY);
	if (str == INI_INVALID_KEY)
		return notfound;

	while (str[0] == ' ')
		str++;
	if (!strncasecmp(str, "Inf", (size_t) 3))
		return INT_MAX;
	if (!strncasecmp(str, "+Inf", (size_t) 4))
		return INT_MAX;
	if (!strncasecmp(str, "-Inf", (size_t) 4))
		return INT_MIN;
	if (!my_is_int(str))
		return notfound;
	return (int) strtol(str, NULL, 0);
}

/**
  @brief    Get the string associated to a key, convert to a double
  @param    d Dictionary to search
  @param    key Key string to look for
  @param    notfound Value to return in case of error
  @return   double

  This function queries a dictionary for a key. A key as read from an
  ini file is given as "section:key". If the key cannot be found,
  the notfound value is returned.
 */
double iniparser_getdouble(dictionary * d, const char *key, double notfound)
{
	char *str = NULL;

	str = iniparser_getstring(d, key, INI_INVALID_KEY);

	if (str == INI_INVALID_KEY)
		return notfound;

	while (str[0] == ' ')
		str++;
	if (!strncasecmp(str, "Inf", (size_t) 3))
		return GSL_POSINF;
	if (!strncasecmp(str, "+Inf", (size_t) 4))
		return GSL_POSINF;
	if (!strncasecmp(str, "-Inf", (size_t) 4))
		return GSL_NEGINF;
	if (!my_is_double(str))
		return notfound;
	return atof(str);
}

/**
  @brief    Get the string associated to a key, convert to a boolean
  @param    d Dictionary to search
  @param    key Key string to look for
  @param    notfound Value to return in case of error
  @return   integer

  This function queries a dictionary for a key. A key as read from an
  ini file is given as "section:key". If the key cannot be found,
  the notfound value is returned.

  A true boolean is found if one of the following is matched:

  - A string starting with 'y'
  - A string starting with 'Y'
  - A string starting with 't'
  - A string starting with 'T'
  - A string starting with '1'

  A false boolean is found if one of the following is matched:

  - A string starting with 'n'
  - A string starting with 'N'
  - A string starting with 'f'
  - A string starting with 'F'
  - A string starting with '0'

  The notfound value returned if no boolean is identified, does not
  necessarily have to be 0 or 1.
 */
int iniparser_getboolean(dictionary * d, const char *key, int notfound)
{
	char *c = NULL;
	int ret = notfound;

	c = iniparser_getstring(d, key, INI_INVALID_KEY);
	if (c == INI_INVALID_KEY)
		return notfound;
	if (c[0] == 'y' || c[0] == 'Y' || c[0] == '1' || c[0] == 't' || c[0] == 'T') {
		ret = 1;
	} else if (c[0] == 'n' || c[0] == 'N' || c[0] == '0' || c[0] == 'f' || c[0] == 'F') {
		ret = 0;
	} else {
		ret = notfound;
	}
	return ret;
}

/**
  @brief    Finds out if a given entry exists in a dictionary
  @param    ini     Dictionary to search
  @param    entry   Name of the entry to look for
  @return   integer 1 if entry exists, 0 otherwise

  Finds out if a given entry exists in the dictionary. Since sections
  are stored as keys with NULL associated values, this is the only way
  of querying for the presence of sections in a dictionary.
 */
int iniparser_find_entry(dictionary * ini, char *entry)
{
	int found = 0;

	if (iniparser_getstring(ini, entry, INI_INVALID_KEY) != INI_INVALID_KEY) {
		found = 1;
	}
	return found;
}

/**
  @brief    Set an entry in a dictionary.
  @param    ini     Dictionary to modify.
  @param    entry   Entry to modify (entry name)
  @param    val     New value to associate to the entry.
  @return   int 0 if Ok, -1 otherwise.

  If the given entry can be found in the dictionary, it is modified to
  contain the provided value. If it cannot be found, -1 is returned.
  It is Ok to set val to NULL.
 */
int iniparser_setstr(dictionary * ini, char *entry, char *val)
{
	dictionary_set(ini, MY_STRING_LOWERCASE(entry), val);
	return 0;
}

/**
  @brief    Delete an entry in a dictionary
  @param    ini     Dictionary to modify
  @param    entry   Entry to delete (entry name)
  @return   void

  If the given entry can be found, it is deleted from the dictionary.
 */
void iniparser_unset(dictionary * ini, char *entry)
{
	dictionary_unset(ini, MY_STRING_LOWERCASE(entry));
}

/**
  @brief    Parse an ini file and return an allocated dictionary object
  @param    ininame Name of the ini file to read.
  @return   Pointer to newly allocated dictionary

  This is the parser for ini files. This function is called, providing
  the name of the file to be read. It returns a dictionary object that
  should not be accessed directly, but through accessor functions
  instead.

  The returned dictionary must be freed using iniparser_freedict().
 */
char *iniparser_getline(FILE * fp)
{
	if (feof(fp)) {
		return NULL;
	}

	int debug = 0;
	size_t len = 0, len_buf = 0;
	char *buf = NULL;
	int c;

	len_buf = BUFSIZ;
	buf = Calloc(len_buf, char);
	buf[0] = '\0';

	do {
		c = fgetc(fp);
		if (c == '\r') {			       /* could be \r\n */
			int cc;

			cc = fgetc(fp);
			if (cc != '\n' && cc != EOF) {
				ungetc(cc, fp);
			}
		}
		if (c == EOF || c == '\n' || c == '\r') {
			if (debug) {
				printf("GETLINE RETURNS [%s]\n", buf);
			}
			return (buf);
		}
		buf[len] = c;
		buf[len + 1] = '\0';
		len++;

		if (len >= len_buf) {
			len_buf += BUFSIZ;
			buf = Realloc(buf, len_buf, char);
		}
	} while (1);

	return buf;
}
dictionary *iniparser_load(const char *ininame)
{
	dictionary *d = NULL;
	char *lin = NULL;
	char *sec = NULL;
	char *key = NULL;
	char *val = NULL;
	char *where = NULL;
	FILE *ini = NULL;
	int lineno = 0;
	size_t len_str = 0;

	if ((ini = fopen(ininame, "r")) == NULL) {
		return NULL;
	}

	/*
	 * Initialize a new dictionary entry
	 */
	if (!(d = dictionary_new(0))) {
		fclose(ini);
		return NULL;
	}
	while ((lin = iniparser_getline(ini)) != NULL) {

		if (len_str == 0) {
			sec = Calloc(strlen(lin) + 1, char);
			key = Calloc(strlen(lin) + 1, char);
			val = Calloc(strlen(lin) + 1, char);
			len_str = strlen(lin) + 1;
		} else {
			if (strlen(lin) + 1 > len_str) {
				len_str = strlen(lin) + 1;
				sec = (char *) realloc((void *) sec, len_str);
				key = (char *) realloc((void *) key, len_str);
				val = (char *) realloc((void *) val, len_str);
			}
		}

		lineno++;

		// if (!(lineno % 1000)) printf("lineno %d\n", lineno);

		where = strskp(lin);			       /* Skip leading spaces */
		if (*where == ';' || *where == '#' || *where == 0)
			continue;			       /* Comment lines */
		else {
			if (sscanf(where, INIPARSER_SECSEP "%[^" INIPARSER_SECSEP "]" INIPARSER_SECSEP, sec) == 1) {
				/*
				 * Valid section name 
				 */
				// printf("sec %s\n", sec);
				iniparser_add_entry(d, sec, NULL, NULL);
			} else if (sscanf(where, "%[^=] = \"%[^\"]\"", key, val) == 2
				   || sscanf(where, "%[^=] = '%[^\']'", key, val) == 2 || sscanf(where, "%[^=] = %[^#]", key, val) == 2) {
				strcpy(key, MY_STRING_LOWERCASE(strcrop(key)));
				/*
				 * sscanf cannot handle "" or '' as empty value,
				 * this is done here
				 */
				if (!strcmp(val, "\"\"") || !strcmp(val, "''")) {
					val[0] = (char) 0;
				} else {
					strcpy(val, strcrop(val));
				}
				iniparser_add_entry(d, sec, key, val);
			}
		}
		Free(lin);
	}

	Free(sec);
	Free(key);
	Free(val);
	fclose(ini);

	return d;
}

/**
  @brief    Free all memory associated to an ini dictionary
  @param    d Dictionary to free
  @return   void

  Free all memory associated to an ini dictionary.
  It is mandatory to call this function before the dictionary object
  gets out of the current context.
 */
void iniparser_freedict(dictionary * d)
{
	dictionary_del(d);
}
