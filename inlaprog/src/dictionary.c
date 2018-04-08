
/**
   @file	dictionary.c
   @author	N. Devillard
   @date	Aug 2000
   @version	$Revision: 1.30 $
   @brief	Implements a dictionary for string variables.

   This module implements a simple dictionary object, i.e. a list of string/string associations.
   This object is useful to store e.g.  informations retrieved from a
   configuration file (ini files).  */


/*
	$Id: dictionary.c,v 1.30 2009/03/30 16:47:04 hrue Exp $
	$Author: hrue $
	$Date: 2009/03/30 16:47:04 $
	$Revision: 1.30 $
*/

#include "dictionary.h"
#include "iniparser.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include "strlib.h"
#include "my-fix.h"
#include "GMRFLib/GMRFLib.h"

/** Maximum value size for integers and doubles. */
#define MAXVALSZ	80

/** Minimal allocated number of entries in a dictionary */
#define DICTMINSZ	262144

/** Invalid key token */
#define DICT_INVALID_KEY    ((char*)-1)

static void *mem_double(void *ptr, int size)
{
	void *newptr = NULL;

	newptr = calloc((size_t) (2 * size), 1);
	memcpy(newptr, ptr, (size_t) size);
	free(ptr);
	return newptr;
}

/**
  @brief	Compute the hash key for a string.
  @param	key		Character string to use for key.
  @return	1 unsigned int on at least 32 bits.

  This hash function has been taken from an Article in Dr Dobbs Journal.
  This is normally a collision-free function, distributing keys evenly.
  The key is stored anyway in the struct so that collision can be avoided
  by comparing the key itself in last resort.
 */
unsigned dictionary_hash(char *key)
{
	unsigned hash;
	size_t i, len;

	len = strlen(key);
	for (hash = 0, i = 0; i < len; i++) {
		hash += (unsigned) key[i];
		hash += (hash << 10);
		hash ^= (hash >> 6);
	}
	hash += (hash << 3);
	hash ^= (hash >> 11);
	hash += (hash << 15);
	return hash;
}

/**
  @brief	Create a new dictionary object.
  @param	size	Optional initial size of the dictionary.
  @return	1 newly allocated dictionary objet.

  This function allocates a new dictionary object of given size and returns
  it. If you do not know in advance (roughly) the number of entries in the
  dictionary, give size=0.
 */
dictionary *dictionary_new(int size)
{
	int i;
	dictionary *d = NULL;

	/*
	 * If no size was specified, allocate space for DICTMINSZ 
	 */
	if (size < DICTMINSZ)
		size = DICTMINSZ;

	if (!(d = (dictionary *) calloc(1, sizeof(dictionary)))) {
		return NULL;
	}
	d->size = size;
	d->used = (char *) calloc((size_t) size, sizeof(char));
	d->val = (char **) calloc((size_t) size, sizeof(char *));
	d->key = (char **) calloc((size_t) size, sizeof(char *));
	map_stri_init_hint(&(d->strihash), (size_t) size);
	map_ii_init_hint(&(d->iihash), (size_t) size);

	for (i = 0; i < size; i++) {
		map_ii_set(&(d->iihash), i, 1);
	}

	return d;
}

/**
  @brief	Delete a dictionary object
  @param	d	dictionary object to deallocate.
  @return	void

  Deallocate a dictionary object and all memory associated to it.
 */
void dictionary_del(dictionary * d)
{
	int i;

	if (d == NULL)
		return;
	for (i = 0; i < d->size; i++) {
		if (d->key[i] != NULL)
			free(d->key[i]);
		if (d->val[i] != NULL)
			free(d->val[i]);
	}
	free(d->val);
	free(d->key);
	free(d->used);
	map_stri_free(&(d->strihash));
	map_ii_free(&(d->iihash));
	free(d);

	return;
}

/**
  @brief	Get a value from a dictionary.
  @param	d		dictionary object to search.
  @param	key		Key to look for in the dictionary.
  @param    def     Default value to return if key not found.
  @return	1 pointer to internally allocated character string.

  This function locates a key in a dictionary and returns a pointer to its
  value, or the passed 'def' pointer if no such key can be found in
  dictionary. The returned character pointer points to data internal to the
  dictionary object, you should not try to free it or modify it.
 */
char *dictionary_get(dictionary * d, char *key, char *def)
{
	int i, *ip;

	ip = map_stri_ptr(&(d->strihash), key);
	if (ip) {
		i = *ip;
		d->used[i] = 1;
		return dictionary_replace_variables(d, d->val[i]);
	} else {
		return dictionary_replace_variables(d, def);
	}
}

/**
  @brief	Get a value from a dictionary, as a char.
  @param	d		dictionary object to search.
  @param	key		Key to look for in the dictionary.
  @param	def		Default value for the key if not found.
  @return 	char	

  This function locates a key in a dictionary using dictionary_get,
  and returns the first char of the found string.
 */
char dictionary_getchar(dictionary * d, char *key, char def)
{
	char *v = NULL;

	if ((v = dictionary_get(d, key, DICT_INVALID_KEY)) == DICT_INVALID_KEY) {
		return def;
	} else {
		return v[0];
	}
}

/**
  @brief	Get a value from a dictionary, as an int.
  @param	d		dictionary object to search.
  @param	key		Key to look for in the dictionary.
  @param	def		Default value for the key if not found.
  @return	int

  This function locates a key in a dictionary using dictionary_get,
  and applies atoi on it to return an int. If the value cannot be found
  in the dictionary, the default is returned.
 */
int dictionary_getint(dictionary * d, char *key, int def)
{
	char *v = NULL;

	if ((v = dictionary_get(d, key, DICT_INVALID_KEY)) == DICT_INVALID_KEY) {
		return def;
	} else {
		return atoi(v);
	}
}

/**
  @brief		Get a value from a dictionary, as a double.
  @param	d		dictionary object to search.
  @param	key		Key to look for in the dictionary.
  @param	def		Default value for the key if not found.
  @return	double

  This function locates a key in a dictionary using dictionary_get,
  and applies atof on it to return a double. If the value cannot be found
  in the dictionary, the default is returned.
 */
double dictionary_getdouble(dictionary * d, char *key, double def)
{
	char *v = NULL;

	if ((v = dictionary_get(d, key, DICT_INVALID_KEY)) == DICT_INVALID_KEY) {
		return def;
	} else {
		return atof(v);
	}
}

/**
  @brief	Set a value in a dictionary.
  @param	d		dictionary object to modify.
  @param	key		Key to modify or add.
  @param	val 	Value to add.
  @return	void

  If the given key is found in the dictionary, the associated value is
  replaced by the provided one. If the key cannot be found in the
  dictionary, it is added to it.

  It is Ok to provide a NULL value for val, but NULL values for the dictionary
  or the key are considered as errors: the function will return immediately
  in such a case.

  Notice that if you dictionary_set a variable to NULL, a call to
  dictionary_get will return a NULL value: the variable will be found, and
  its value (NULL) is returned. In other words, setting the variable
  content to NULL is equivalent to deleting the variable from the
  dictionary. It is not possible (in this implementation) to have a key in
  the dictionary without value.
 */
void dictionary_set(dictionary * d, char *key, char *val)
{
	int i = 0, *ip = NULL;

	if (d == NULL || key == NULL)
		return;

	ip = map_stri_ptr(&(d->strihash), key);
	if (ip) {
		i = *ip;

		/*
		 * Found a value: modify and return 
		 */
		if (d->val[i] != NULL)
			free(d->val[i]);
		d->val[i] = val ? strdup(val) : NULL;
		map_stri_set(&(d->strihash), d->val[i], 0);
	} else {
		/*
		 * Add a new value 
		 */
		/*
		 * See if dictionary needs to grow 
		 */
		if (d->n == d->size) {

			/*
			 * Reached maximum size: reallocate blackboard 
			 */
			d->used = (char *) mem_double(d->used, (int) (d->size * sizeof(char)));
			d->val = (char **) mem_double(d->val, (int) (d->size * sizeof(char *)));
			d->key = (char **) mem_double(d->key, (int) (d->size * sizeof(char *)));
			for (i = d->size; i < 2 * d->size; i++)
				map_ii_set(&(d->iihash), i, 1);
			d->size *= 2;
		}

		int j;
		for (j = -1; (j = (int) map_ii_next(&(d->iihash), j)) != -1;) {
			i = d->iihash.contents[j].key;
			break;
		}
		assert(d->key[i] == NULL);

		if (0) {
			/*
			 * Insert key in the first empty slot 
			 */
			for (i = 0; i < d->size; i++) {
				if (d->key[i] == NULL) {
					/*
					 * Add key here 
					 */
					break;
				}
			}
		}
		/*
		 * Copy key 
		 */
		d->key[i] = strdup(key);
		d->val[i] = val ? strdup(val) : NULL;
		d->used[i] = 0;
		map_stri_set(&(d->strihash), d->key[i], i);
		map_ii_remove(&(d->iihash), i);
		// printf("ADD VALUE [%s] = [%s]\n", d->key[i], d->val[i]);
		d->n++;
	}

	if (0) {
		printf("\n");
		dictionary_dump(d, stdout);
		printf("\n");
	}

	return;
}

/**
  @brief	Delete a key in a dictionary
  @param	d		dictionary object to modify.
  @param	key		Key to remove.
  @return   void

  This function deletes a key in a dictionary. Nothing is done if the
  key cannot be found.
 */
void dictionary_unset(dictionary * d, char *key)
{
	int i, *ip;

	if (key == NULL) {
		return;
	}

	ip = map_stri_ptr(&(d->strihash), key);

	if (!ip)
		return;

	map_stri_remove(&(d->strihash), key);

	i = *ip;
	map_ii_set(&(d->iihash), i, 1);

	free(d->key[i]);
	d->key[i] = NULL;
	if (d->val[i] != NULL) {
		free(d->val[i]);
		d->val[i] = NULL;
	}
	d->used[i] = 0;
	d->n--;
	return;
}

/**
  @brief	Set a key in a dictionary, providing an int.
  @param	d		Dictionary to update.
  @param	key		Key to modify or add
  @param	val		Integer value to store (will be stored as a string).
  @return	void

  This helper function calls dictionary_set() with the provided integer
  converted to a string using %d.
 */
void dictionary_setint(dictionary * d, char *key, int val)
{
	char sval[MAXVALSZ];

	sprintf(sval, "%d", val);
	dictionary_set(d, key, sval);
}

/**
  @brief	Set a key in a dictionary, providing a double.
  @param	d		Dictionary to update.
  @param	key		Key to modify or add
  @param	val		Double value to store (will be stored as a string).
  @return	void

  This helper function calls dictionary_set() with the provided double
  converted to a string using %g.
 */
void dictionary_setdouble(dictionary * d, char *key, double val)
{
	char sval[MAXVALSZ];

	sprintf(sval, "%g", val);
	dictionary_set(d, key, sval);
}

/**
  @brief	Dump a dictionary to an opened file pointer.
  @param	d	Dictionary to dump
  @param	f	Opened file pointer.
  @return	void

  Dumps a dictionary onto an opened file pointer. Key pairs are printed out
  as @c [Key]=[Value], one per line. It is Ok to provide stdout or stderr as
  output file pointers.
 */
void dictionary_dump(dictionary * d, FILE * out)
{
	int i;

	if (d == NULL || out == NULL)
		return;
	if (d->n < 1) {
		fprintf(out, "empty dictionary\n");
		return;
	}
	for (i = 0; i < d->size; i++) {
		if (d->key[i]) {
			fprintf(out, "%-40s\t[%s]\n", d->key[i], d->val[i] ? d->val[i] : "NULL");
		}
	}
	return;
}
char *dictionary_replace_variables(dictionary * d, char *str)
{
	/*
	 * variable-replacement! added by hrue@math.ntnu.no 
	 */

	if (!str || str == DICT_INVALID_KEY) {
		return str;
	}

	char *first = NULL, *last = NULL, *newstr = NULL, *envvar = NULL, *var = NULL, *newvar = NULL;
	int c, debug = 0;
	size_t len;
	int i;

	if (debug) {
		printf("ENTER with string = [%s] strlen=%1lu\n", str, (unsigned long) strlen(str));
	}
	first = strchr(str, '$');
	if (!first) {
		return str;
	}
	last = first + 1;
	c = *last;

	while (*last && ((c >= 'a' && c <= 'z') || (c >= 'A' && c <= 'Z'))) {
		c = *++last;
	}
	len = last - first;

	if (debug) {
		printf("first %s\n", first);
		printf("last %s\n", last);
	}

	var = (char *) calloc(len + 1, sizeof(char));
	var[0] = INIPARSER_SEP;
	strncpy(var + 1, first + 1, len);
	var[len] = '\0';
	if (debug) {
		printf("var is [%s]\n", var);
	}

	int *ip;
	ip = map_stri_ptr(&(d->strihash), var);
	if (ip) {
		i = *ip;
		if (debug) {
			printf("replace [%s] with [%s]\n", d->key[i], d->val[i]);
		}

		*first = '\0';
		GMRFLib_sprintf(&newstr, "%s%s%s", str, d->val[i], last);
		if (debug) {
			printf("new string [%s]\n", newstr);
		}
		return dictionary_replace_variables(d, newstr);

	} else {
		if (debug)
			printf("var not in strhash...\n");
	}

	if (debug) {
		printf("ask for envvar [%s]\n", &var[1]);
	}

	GMRFLib_sprintf(&newvar, "inla_%s", &var[1]);
	envvar = getenv(newvar);

	if (envvar) {
		*first = '\0';
		GMRFLib_sprintf(&newstr, "%s%s%s", str, envvar, last);
		if (debug) {
			printf("new string [%s]\n", newstr);
		}

		return dictionary_replace_variables(d, newstr);
	} else {
		if (debug)
			printf("getenv(%s) returns NULL\n", newvar);
	}

	return str;
}
int dictionary_dump_unused(dictionary * d, FILE * out)
{
	/*
	 * dump unused entries, added by hrue@math.ntnu.no
	 */

	int i, count = 0, first = 1;

	if (d == NULL || out == NULL)
		return count;
	if (d->n < 1) {
		fprintf(out, "empty dictionary\n");
		return count;
	}

	for (i = 0; i < d->size; i++) {
		if (d->key[i] && !d->used[i] && d->key[i][0] != INIPARSER_SEP && strchr(d->key[i], INIPARSER_SEP) != NULL) {

			if (first) {
				first = 0;
				fprintf(out, "\n\n");
			}
			fprintf(out, "\tUnused entry [%s]=[%s]\n", d->key[i], d->val[i] ? d->val[i] : "NULL");
			count++;
		}
	}
	return count;
}

/* Example code */
#ifdef TESTDIC
#define NVALS 20000
int main(int argc, char *argv[])
{
	dictionary *d = NULL;
	char *val = NULL;
	int i;
	char cval[90];

	/*
	 * allocate blackboard 
	 */
	printf("allocating...\n");
	d = dictionary_new(0);

	/*
	 * Set values in blackboard 
	 */
	printf("setting %d values...\n", NVALS);
	for (i = 0; i < NVALS; i++) {
		sprintf(cval, "%04d", i);
		dictionary_set(d, cval, "salut");
	}
	printf("getting %d values...\n", NVALS);
	for (i = 0; i < NVALS; i++) {
		sprintf(cval, "%04d", i);
		val = dictionary_get(d, cval, DICT_INVALID_KEY);
		if (val == DICT_INVALID_KEY) {
			printf("cannot get value for key [%s]\n", cval);
		}
	}
	printf("unsetting %d values...\n", NVALS);
	for (i = 0; i < NVALS; i++) {
		sprintf(cval, "%04d", i);
		dictionary_unset(d, cval);
	}
	if (d->n != 0) {
		printf("error deleting values\n");
	}

	printf("deallocating...\n");
	dictionary_del(d);
	return 0;
}
#endif
