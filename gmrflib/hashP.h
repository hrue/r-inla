
/* hashP.h
 * 
 * Copyright (C) 2001-2006 Havard Rue
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
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
 *
 */

/*!
  \file hashP.h
  \brief Include this file if you want to use the hash-functions.

  If you want to use the hash-functions which is used in the GMRFLib-library include \ref hashP.h. The documentation of the
  hash-function are as follows.

  \verbinclude hash_doc.txt

  The copyright information is available in the file \ref hashP.h
*/

/*

Copyright (c) 2002-2004, Jean-Sebastien Roy (js@jeannot.org)

Permission is hereby granted, free of charge, to any person obtaining a
copy of this software and associated documentation files (the
"Software"), to deal in the Software without restriction, including
without limitation the rights to use, copy, modify, merge, publish,
distribute, sublicense, and/or sell copies of the Software, and to
permit persons to whom the Software is furnished to do so, subject to
the following conditions:

The above copyright notice and this permission notice shall be included
in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

*/

/*

# MAPKIT, Version 1.4
# Copyright J.S. Roy (js@jeannot.org), 2002-2004
# See the LICENSE file for copyright information.
# @(#) $Jeannot: README,v 1.32 2004/04/16 13:48:50 js Exp $

Mapkit is a simple set of C functions to create and access maps (aka
dictionnaries, hash tables) and sparse structures (vectors, matrices).
The last version (and other software) is available at the URL :
http://www.jeannot.org/~js/code/index.en.html

This version provide the following maps :
map_ii     : int -> int (default value = 0)
map_id     : int -> double (default value = 0.0)
map_ivp    : int -> void * (default value = NULL)

map_vpi     : void* -> int (default value = 0)
map_vpd     : void* -> double (default value = 0.0)
map_vpvp    : void* -> void * (default value = NULL)

map_h_* : ditto, with a different hash function, 
useful when lots of collisions occurs when using the previous maps.

map_stri   : string -> int (default value = 0)
map_strd   : string -> double (default value = 0.0)
map_strvp  : string -> void * (default value = NULL)
map_strstr : string -> string (default value = "")
(strings are '\0' terminated char arrays).

spvector   : int(>0) -> double (default value = 0.0)
spmatrix   : int*int -> double (default value = 0.0)
By default, these two maps do not store elements whose value is the default
value, and always return the default value when an non existent key is queried.
(ie, the alwaysdefault field of the map is set. See below.)

For each map, the following basic functions are provided:
(prefix the function with the name of the map)
 _init(&map) : initialize the map (returns an error code)
 _(&map, key) : macro that returns *(a pointer to the value stored at 'key'
                (insert key if key is missing)). Can be used for reading
                and writing. May fail if there is not enough memory.
 _remove(&map, key) : remove key (key must exists. returns an error code)
 _free(&map) : free the map

The following functions are provided to tune the performance :
 _init_hint(&map, used) : initialize the map for 'used' elements
 _ensurecapacity(&map, used) : ensure at least 'used' elements can be stored
                               in the map
 _adjustcapacity(&map) : shrink the map to its smallest size
 _printstats(&map) : print statistics about the map's usage, and collisions 
                     if MAPKIT_COLLISIONS is defined during compilation.

Functions for more specific operations:
 _ptr(&map, key) : returns a pointer to the value stored at 'key'
                   or NULL if key is missing.
 _insertptr(&map, key) : returns a pointer to the value stored at 'key'
                         (insert if key is missing,
                         May return NULL if there is not enough memory.)
 _removeptr(&map, ptr) : remove the key pointing to value *ptr
 _value(&map, key) : returns the value at key. Fails if key is missing.

Iterators:
 _next(&map, index) : returns the next index to a full slot
 Typical use to scan all full slots of a map:
 for (i = -1 ; (i = map_next(&map, i)) != -1 ; )
 map.contents[i].key contains the key,
    map.contents[i].value contains the value 

 _nextptr(&map, ptr) : returns the next pointer to a full slot
 Typical use to scan all full slots of a map:
 for (ptr = NULL ; (ptr = map_nextptr(&map, ptr)) != NULL ; )
 ptr->key contains the key, ptr->value contains the value 

Query functions that return an error code:
 _set(&map, key, value) : sets the value at key
                          (insert as needed, returns an error code)
 _get(&map, key, &value) : returns the value at key (returns an error code)

Functions the manipulate the whole map:
 _copy : copy a map into a new uninitialized map
 _getall : allocate an array with all the map (key, value) pairs
 _getall_sorted : ditto, sorted.
 _setall : insert in a map all elements of an array with (key, value) pairs
 _removeall : remove from a map all keys from a key array
 _clean : remove all values equal to the defaultvalue.
(see mapkit.h for prototypes)

More advanced functions are also provided : see mapkit.h for the detailled list
of available functions, and mapkit_generic.h for the definition of error codes.

The default value of a map can be changed by setting the defaultvalue field of
the map. The behavior of the map in the case of missing keys can be changed by
setting the alwaysdefault field of the map to 1. In this case, functions never
fail even if the key is missing, and values equal to the defaultvalue are not
stored in the map. Queries return the defaultvalue when the key is missing.

Performance is optimized for large structures and consecutive keys. 
Resizing the hash table is costly, so using hints about its size is encouraged
(e.g. map_init_hint, map_setall and map_removeall functions).
Moreover, like most hash tables, the worst-case behavior (all keys having the
same hash) results in a very low performance (statistics (_printstats) should
be used to detect this behavior). In this case, the hash function should be
changed to a more robust one.

The maps are implemented using open-addressing and double hashing.
Feeback on design and implementation is welcome !

Adjust the Makefile to reflect the correct options for your setup.

Macro definitions:
  MAPKIT_EXITONERROR : if defined, will exit upon failure.
    The return code is the MAPKIT error number.
  MAPKIT_DEBUG : if defined, will print information about errors
    and various events. Messages get printed on stderr.
  MAPKIT_COLLISIONS : if defined, statistics about hash collisions will be
    gathered.

The following targets are provided :
 all     : compiles everything
 test    : test the library performs correcly on your setup
 bench   : benchmark the library
           !! may take a long time and require large amounts of memory !!
           (you should unset MAPKIT_DEBUG before benchmarking.)
 example : an example which counts the occurences of the lines of a file.

If your compiler does not support the C99 'inline' keyword (rare) define an
empty inline macro (eg. -Dinline= ).
(Under Visual C inline is defined to __inline)
If your system does not conform to Single UNIX v3 (rare) and lacks the random() 
call, define the NO_RANDOM macro (-DNO_RANDOM). (This is the default on Windows)

Simple example:

  spvector v;
  double a, b;
  
  spvector_init(&v);
  spvector_set(&v, 4, 7.0); 
  spvector_(&v, 4) = 7.0; 
  a = spvector_value(&v, 4);
  a = spvector_(&v, 4); 
  b = spvector_value(&v, 12);
  spvector_free(&v);

(See mapkit_example.c for a more complete example.)

*/

/*
  THIS FILE CONTAINS the source mapkit.h and mapkit_generic.h

*/

/* MapKit, Generic functions */

/*
 * Copyright (c) 2002-2004, Jean-Sebastien Roy (js@jeannot.org)
 * 
 * Permission is hereby granted, free of charge, to any person obtaining a
 * copy of this software and associated documentation files (the
 * "Software"), to deal in the Software without restriction, including
 * without limitation the rights to use, copy, modify, merge, publish,
 * distribute, sublicense, and/or sell copies of the Software, and to
 * permit persons to whom the Software is furnished to do so, subject to
 * the following conditions:
 * 
 * The above copyright notice and this permission notice shall be included
 * in all copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
 * OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
 * IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
 * CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
 * TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
 * SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

/* @(#) $Jeannot: mapkit_generic.h,v 1.31 2004/04/10 14:17:03 js Exp $ */

/*
Important defines for MapKit:
  MAPKIT_EXITONERROR : if defined, will exit upon failure.
    The return code is the MAPKIT error number.
  MAPKIT_DEBUG : if defined, will print information about errors
    and various events. Messages get printed on stderr.
  MAPKIT_COLLISIONS : if defined, statistics about hash collisions will be
    gathered.
*/

#ifndef _MAPKIT_GENERIC_
#define _MAPKIT_GENERIC_

#include <time.h>
#include <strings.h>
#if !defined(__FreeBSD__)
#include <malloc.h>
#endif
#include <stdlib.h>
#include <stddef.h>
#include <stdio.h>
#include <limits.h>
#include <string.h>

#include <stdint.h>
#include <inttypes.h>

#define MAPKIT_VERSION "1.4"

/* MS Visual C does not recognize 'inline' */
#if defined(_MSC_VER) && !defined(__cplusplus) && !defined(inline)
#define inline __inline
#endif

/*
  added from: see http://www.greenend.org.uk/rjk/2003/03/inline.html

  hrue added __STRICT_ANSI__ test
*/
#ifdef __STRICT_ANSI__
#define INLINE static
#else
#ifndef INLINE
#if __GNUC__
#define INLINE static inline
#else
#define INLINE inline
#endif
#endif
#endif

#if defined(__cplusplus)
extern "C" {
#endif

	typedef enum {
		MAPKIT_OK = 0,
		MAPKIT_ETOOBIG = 1,
		MAPKIT_ENOMEM = 2,
		MAPKIT_EBADKEY = 3,
		MAPKIT_EKEYNOTFOUND = 4
	} mapkit_error;

/* Constants */
	extern const char *mapkit_error_msg[5];

/* Macros */

/* print an appropriate error message if _err != MAPKIT_OK */
#define CHKMAPKIT(_err) { mapkit_error __err = (_err); if (__err) fprintf(stderr, "MAPKIT: ERROR: %s (%s, line %d)\n", mapkit_error_msg[__err], __FILE__, __LINE__); }

/* if MAPKIT_DEBUG is set, most mapkit functions will output various
   informational & error messages */
//#ifdef MAPKIT_DEBUG
#if 1
#define MAPKIT_MSG(_err) fprintf(stderr, "MAPKIT: ERROR: %s (%s, line %d)\n", mapkit_error_msg[_err], __FILE__, __LINE__)
#else
#define MAPKIT_MSG(_err)
#endif

#define MAPKIT_FATAL_ERROR(_err) { MAPKIT_MSG(_err); exit(_err); }

/* if MAPKIT_EXITONERROR is set, the program will quit instead of returning
   an error message */
//#ifdef MAPKIT_EXITONERROR
#if 1
#define MAPKIT_ERROR(_err) MAPKIT_FATAL_ERROR(_err)
#define MAPKIT_ERROR_NORET(_err) MAPKIT_FATAL_ERROR(_err)
#else
#define MAPKIT_ERROR(_err) { MAPKIT_MSG(_err); return (_err); }
#define MAPKIT_ERROR_NORET(_err) MAPKIT_MSG(_err)
#endif

/* return index value if key not found */
#define MAPKIT_KEYNOTFOUND (-1)

/* special key values */
#define MAPKIT_FULLSLOT (0)
#define MAPKIT_FREESLOT (-1)
#define MAPKIT_DELETEDSLOT (-2)

/* default initial size */
#define MAPKIT_DEFAULT_EXPECTEDUSED (100)


/* 
   hrue: need to make sure that the SIZE_T is the same size as a pointer!
   
   MAPKIT defines this default as 'long' but long is 4 bit for Windows 64... so we use the new c99 invention instead.
*/
	
#define mapkit_size_t intptr_t /* sizes, signed to accomodate MAPKIT_KEYNOTFOUND */
#define mapkit_hash_t uintptr_t /* hash values, must be unsigned (for index, decrement) */

#define MAPKIT_SIZE_T_MAX __INTPTR_MAX__
	
/* find the lowest number in mapkit_primes greater than n */
	extern mapkit_size_t mapkit_nextprime(const mapkit_size_t n);

/* string hash function */
	extern mapkit_hash_t mapkit_strhash(const char *string);

/* memory hash function */
	extern mapkit_hash_t mapkit_memhash(const void *data, const size_t len);

/* integer hash function */
	INLINE mapkit_hash_t mapkit_hash(mapkit_hash_t key);

/* Thomas Wang's integer hash function */
#if MAPKIT_SIZE_T_MAX <= 0x7fffffff

/* 32 bits version */
	mapkit_hash_t mapkit_hash(mapkit_hash_t key) {
		key += ~(key << 15);
		key ^= (key >> 10);
		key += (key << 3);
		key ^= (key >> 6);
		key += ~(key << 11);
		key ^= (key >> 16);
		return key;
	}
#else

/* 64 bits version */
	mapkit_hash_t mapkit_hash(mapkit_hash_t key) {
		key += ~(key << 32);
		key ^= (key >> 22);
		key += ~(key << 13);
		key ^= (key >> 8);
		key += (key << 3);
		key ^= (key >> 15);
		key += ~(key << 27);
		key ^= (key >> 31);
		return key;
	}
#endif
#if defined(__cplusplus)
}
#endif
#endif							       /* _MAPKIT_GENERIC_ */
   /*
    *  
       *//*
       * MapKit 
		      *//*
		      * Copyright (c) 2002-2004, Jean-Sebastien Roy (js@jeannot.org)
		      * 
		      * Permission is hereby granted, free of charge, to any person obtaining a
		      * copy of this software and associated documentation files (the
		      * "Software"), to deal in the Software without restriction, including
		      * without limitation the rights to use, copy, modify, merge, publish,
		      * distribute, sublicense, and/or sell copies of the Software, and to
		      * permit persons to whom the Software is furnished to do so, subject to
		      * the following conditions:
		      * 
		      * The above copyright notice and this permission notice shall be included
		      * in all copies or substantial portions of the Software.
		      * 
		      * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
		      * OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
		      * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
		      * IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
		      * CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
		      * TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
		      * SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
		      
	  *//*
	  * @(#) $Jeannot: mapkit.m4h,v 1.30 2004/04/05 18:43:41 js Exp $ 
									     *//*
									     * MapKit Copyright (c) 2002-2004, Jean-Sebastien Roy (js@jeannot.org)
									     * 
									     * @(#) $Jeannot: mapkit.m4,v 1.83 2004/04/10 14:17:03 js Exp $ 
	 *//*
	 * @(#) $Jeannot: mapkit_defs.m4,v 1.20 2004/03/31 19:12:55 js Exp $ 
	 */
#ifndef _MAPKIT_
#define _MAPKIT_

/* MS Visual C does not recognize 'inline' */
#if defined(_MSC_VER) && !defined(__cplusplus) && !defined(inline)
#define inline __inline
#endif
#if defined(__cplusplus)
extern "C" {
#endif

/* Map structures */

/*
  Prototypes for map_ii (int -> int)
  Default value : 0
  Uses a state field.
*/

	typedef struct {
		signed char state;
		int key;
		int value;
	} map_ii_storage;

	typedef struct {
		int key;
		int value;
	} map_ii_element;

	typedef struct {
		mapkit_size_t size,			       /* size of the hash table */
		 fill,					       /* # of non-free slots in the hash table */
		 used,					       /* # of full slots */
		 maxfill,				       /* max # of non-free slots before rehash */
		 minused /* min # of full slots before rehash */ ;

		double maxfillfactor,			       /* maxfill = maxfillfactor * size */
		 minusedfactor;				       /* minused = minusedfactor * size */

#ifdef MAPKIT_COLLISIONS
		unsigned long insertionindexs, insertionindex_collisions;
		unsigned long keyindexs, keyindex_collisions;
#endif

		int defaultvalue;

		int alwaysdefault;

		map_ii_storage *contents;
	} map_ii;

/* Macros */

/* set to inform the type is available */
#define MAPKIT_map_ii

/* read/write access to the value at key (insert as needed) */
#define map_ii_(spm, key) (*(map_ii_insertptr(spm, key)))

/* true if key exists */
#define map_ii_haskey(spm, key) ((map_ii_ptr(spm, key)) != NULL)

/* Prototypes */

/* initialize the structure */
	extern mapkit_error map_ii_init(map_ii * spm);

/* free contents */
	extern void map_ii_free(map_ii * spm);

/* copy the "from" map into a new "to" map */
	extern mapkit_error map_ii_copy(map_ii * to, map_ii * from);

/* initialize the structure for "used" elements */
	extern mapkit_error map_ii_init_hint(map_ii * spm, mapkit_size_t used);

/* ensure at least "used" slots (including currently occupied slots) */
	extern mapkit_error map_ii_ensurecapacity(map_ii * spm, mapkit_size_t used);

/* restore maxfill & minused and ensure these constraints */
	extern mapkit_error map_ii_adjustcapacity(map_ii * spm);

/* reallocate contents for "size" slots, keeping previous values */
	extern mapkit_error map_ii_reallocate(map_ii * spm, mapkit_size_t size);

/* return an adequate size for growing from "used" elements */
	extern mapkit_size_t map_ii_growsize(map_ii * spm, mapkit_size_t used);

/* return an adequate size for shrinking from "used" elements */
	extern mapkit_size_t map_ii_shrinksize(map_ii * spm, mapkit_size_t used);

/* return an adequate size for "used" elements */
	extern mapkit_size_t map_ii_meansize(map_ii * spm, mapkit_size_t used);

/* pointer to the value at key (insert as needed) */
	INLINE int *map_ii_insertptr(map_ii * spm, int key);

/* pointer to the value at key or NULL */
	INLINE int *map_ii_ptr(map_ii * spm, int key);

/* remove the key pointing to the value *ptr. No checking is done */
	INLINE mapkit_error map_ii_removeptr(map_ii * spm, int *ptr);

/* set the value at key (insert as needed) */
	INLINE mapkit_error map_ii_set(map_ii * spm, int key, int value);

/* return the value at key (key must exists if alwaysdefaultvalue is not set) */
	INLINE int map_ii_value(map_ii * spm, int key);

/* get the value at key (return an error code) */
	INLINE mapkit_error map_ii_get(map_ii * spm, int key, int *value);

/* remove key (key must exists) */
	INLINE mapkit_error map_ii_remove(map_ii * spm, int key);

/*
 * return the next index to a full slot of the map, or -1 if none
 * if index == -1, returns the first index, or -1 if none
 * the contents of the slot can be accessed by spm->contents[index].key and
 * spm->contents[index].value
 */
	extern mapkit_size_t map_ii_next(map_ii * spm, mapkit_size_t iindex);

/*
 * return the next ptr to a full slot of the map, or NULL if none
 * if pos_contents == NULL, returns the first full slot, or NULL if none
 * the contents of the slot can be accessed by pos_contents->key and
 * pos_contents->value
 */
	extern map_ii_storage *map_ii_nextptr(map_ii * spm, map_ii_storage * pos_contents);

/* allocate an array of map_ii_element with all "count" (key,value) pairs */
	extern mapkit_error map_ii_getall(map_ii * spm, map_ii_element ** array, mapkit_size_t * count);

/* remove all values equal to the defaultvalue */
	extern mapkit_error map_ii_clean(map_ii * spm);

/* compare the key of two map_ii_element (for qsort) */
	extern int map_ii_compare(const void *e1, const void *e2);

/* allocate an array of map_ii_element with all "count" (key,value) pairs */
	extern mapkit_error map_ii_getall_sorted(map_ii * spm, map_ii_element ** array, mapkit_size_t * count);

/* insert all count map_ii_element into spm */
	extern mapkit_error map_ii_setall(map_ii * spm, map_ii_element * array, mapkit_size_t count);

/* remove all "count" elements whose keys are in "array" */
	extern mapkit_error map_ii_removeall(map_ii * spm, int *array, mapkit_size_t count);

/* print statistics */
	extern void map_ii_printstats(map_ii * spm);

/* Private functions */

	extern mapkit_error map_ii_remove_s(map_ii * spm, int key);
	extern int map_ii_value_s(map_ii * spm, int key);
	extern mapkit_error map_ii_get_s(map_ii * spm, int key, int *value);
	extern mapkit_error map_ii_set_s(map_ii * spm, int key, int value);
	extern int *map_ii_insertptr_s(map_ii * spm, int key);
	extern int *map_ii_ptr_s(map_ii * spm, int key);

/* Implementation */

/* INLINEd functions */
	mapkit_error map_ii_remove(map_ii * spm, int key) {
		map_ii_storage *contents;

		 contents = &(spm->contents[((mapkit_hash_t) key) % spm->size]);

		if ((contents->state == MAPKIT_FULLSLOT)
		    && ((contents->key) == (key))) {
			contents->state = MAPKIT_DELETEDSLOT;
#ifdef MAPKIT_COLLISIONS
			spm->keyindexs++;
#endif
			spm->used--;
			if (spm->used < spm->minused) {
#ifdef MAPKIT_DEBUG
				fprintf(stderr, "MAPKIT: used < minused\n");
#endif
				return map_ii_reallocate(spm, map_ii_shrinksize(spm, spm->used));
			}
		} else
			 return map_ii_remove_s(spm, key);
		return MAPKIT_OK;
	}

	mapkit_error map_ii_removeptr(map_ii * spm, int *ptr) {
		map_ii_storage *sptr = (map_ii_storage *) ((char *) ptr - offsetof(map_ii_storage, value));

		sptr->state = MAPKIT_DELETEDSLOT;
		spm->used--;
		if (spm->used < spm->minused) {
#ifdef MAPKIT_DEBUG
			fprintf(stderr, "MAPKIT: used < minused\n");
#endif
			return map_ii_reallocate(spm, map_ii_shrinksize(spm, spm->used));
		}

		return MAPKIT_OK;
	}

	mapkit_error map_ii_set(map_ii * spm, int key, int value) {
		map_ii_storage *contents;

		if (spm->alwaysdefault && ((value) == (spm->defaultvalue)))
			return map_ii_remove(spm, key);

		contents = &(spm->contents[((mapkit_hash_t) key) % spm->size]);

		if (contents->state == MAPKIT_FREESLOT) {
			contents->state = MAPKIT_FULLSLOT;
			contents->key = key;
			contents->value = value;
#ifdef MAPKIT_COLLISIONS
			spm->insertionindexs++;
#endif
			spm->used++;
			spm->fill++;
			if (spm->fill > spm->maxfill) {
#ifdef MAPKIT_DEBUG
				fprintf(stderr, "MAPKIT: fill > maxfill\n");
#endif
				return map_ii_reallocate(spm, map_ii_growsize(spm, spm->used));
			}
			return MAPKIT_OK;
		} else if ((contents->state == MAPKIT_FULLSLOT)
			   && ((contents->key) == (key))) {
			contents->value = value;
#ifdef MAPKIT_COLLISIONS
			spm->insertionindexs++;
#endif
			return MAPKIT_OK;
		} else
			return map_ii_set_s(spm, key, value);
	}

	int map_ii_value(map_ii * spm, int key) {
		map_ii_storage *contents;

		contents = &(spm->contents[((mapkit_hash_t) key) % spm->size]);

		if ((contents->state == MAPKIT_FULLSLOT)
		    && ((contents->key) == (key))) {
#ifdef MAPKIT_COLLISIONS
			spm->keyindexs++;
#endif
			return contents->value;
		} else if (spm->alwaysdefault && (contents->state == MAPKIT_FREESLOT)) {
#ifdef MAPKIT_COLLISIONS
			spm->keyindexs++;
#endif
			return spm->defaultvalue;
		} else
			return map_ii_value_s(spm, key);
	}

	mapkit_error map_ii_get(map_ii * spm, int key, int *value) {
		map_ii_storage *contents;

		contents = &(spm->contents[((mapkit_hash_t) key) % spm->size]);

		if ((contents->state == MAPKIT_FULLSLOT)
		    && ((contents->key) == (key))) {
#ifdef MAPKIT_COLLISIONS
			spm->keyindexs++;
#endif
			*value = contents->value;
		} else if (spm->alwaysdefault && (contents->state == MAPKIT_FREESLOT)) {
#ifdef MAPKIT_COLLISIONS
			spm->keyindexs++;
#endif
			*value = spm->defaultvalue;
		} else
			return map_ii_get_s(spm, key, value);

		return MAPKIT_OK;
	}

	int *map_ii_insertptr(map_ii * spm, int key) {
		map_ii_storage *contents;

		contents = &(spm->contents[((mapkit_hash_t) key) % spm->size]);

		if (contents->state == MAPKIT_FREESLOT) {
			if (spm->fill >= spm->maxfill)
				return map_ii_insertptr_s(spm, key);

			contents->state = MAPKIT_FULLSLOT;
			contents->key = key;
			contents->value = spm->defaultvalue;
#ifdef MAPKIT_COLLISIONS
			spm->insertionindexs++;
#endif
			spm->used++;
			spm->fill++;
			return &(contents->value);
		} else if ((contents->state == MAPKIT_FULLSLOT)
			   && ((contents->key) == (key))) {
#ifdef MAPKIT_COLLISIONS
			spm->insertionindexs++;
#endif
			return &(contents->value);
		} else
			return map_ii_insertptr_s(spm, key);
	}

	int *map_ii_ptr(map_ii * spm, int key) {
		map_ii_storage *contents;

		contents = &(spm->contents[((mapkit_hash_t) key) % spm->size]);

		if ((contents->state == MAPKIT_FULLSLOT)
		    && ((contents->key) == (key))) {
#ifdef MAPKIT_COLLISIONS
			spm->insertionindexs++;
#endif
			if ((!spm->alwaysdefault) || (!((contents->value) == (spm->defaultvalue))))
				return &(contents->value);
			else
				return NULL;
		} else if (contents->state == MAPKIT_FREESLOT) {
#ifdef MAPKIT_COLLISIONS
			spm->insertionindexs++;
#endif
			return NULL;
		} else
			return map_ii_ptr_s(spm, key);
	}

/*
  Prototypes for map_id (int -> double)
  Default value : 0.0
  Uses a state field.
*/

	typedef struct {
		signed char state;
		int key;
		double value;
	} map_id_storage;

	typedef struct {
		int key;
		double value;
	} map_id_element;

	typedef struct {
		mapkit_size_t size,			       /* size of the hash table */
		 fill,					       /* # of non-free slots in the hash table */
		 used,					       /* # of full slots */
		 maxfill,				       /* max # of non-free slots before rehash */
		 minused /* min # of full slots before rehash */ ;

		double maxfillfactor,			       /* maxfill = maxfillfactor * size */
		 minusedfactor;				       /* minused = minusedfactor * size */

#ifdef MAPKIT_COLLISIONS
		unsigned long insertionindexs, insertionindex_collisions;
		unsigned long keyindexs, keyindex_collisions;
#endif

		double defaultvalue;

		int alwaysdefault;

		map_id_storage *contents;
	} map_id;

/* Macros */

/* set to inform the type is available */
#define MAPKIT_map_id

/* read/write access to the value at key (insert as needed) */
#define map_id_(spm, key) (*(map_id_insertptr(spm, key)))

/* true if key exists */
#define map_id_haskey(spm, key) ((map_id_ptr(spm, key)) != NULL)

/* Prototypes */

/* initialize the structure */
	extern mapkit_error map_id_init(map_id * spm);

/* free contents */
	extern void map_id_free(map_id * spm);

/* copy the "from" map into a new "to" map */
	extern mapkit_error map_id_copy(map_id * to, map_id * from);

/* initialize the structure for "used" elements */
	extern mapkit_error map_id_init_hint(map_id * spm, mapkit_size_t used);

/* ensure at least "used" slots (including currently occupied slots) */
	extern mapkit_error map_id_ensurecapacity(map_id * spm, mapkit_size_t used);

/* restore maxfill & minused and ensure these constraints */
	extern mapkit_error map_id_adjustcapacity(map_id * spm);

/* reallocate contents for "size" slots, keeping previous values */
	extern mapkit_error map_id_reallocate(map_id * spm, mapkit_size_t size);

/* return an adequate size for growing from "used" elements */
	extern mapkit_size_t map_id_growsize(map_id * spm, mapkit_size_t used);

/* return an adequate size for shrinking from "used" elements */
	extern mapkit_size_t map_id_shrinksize(map_id * spm, mapkit_size_t used);

/* return an adequate size for "used" elements */
	extern mapkit_size_t map_id_meansize(map_id * spm, mapkit_size_t used);

/* pointer to the value at key (insert as needed) */
	INLINE double *map_id_insertptr(map_id * spm, int key);

/* pointer to the value at key or NULL */
	INLINE double *map_id_ptr(map_id * spm, int key);

/* remove the key pointing to the value *ptr. No checking is done */
	INLINE mapkit_error map_id_removeptr(map_id * spm, double *ptr);

/* set the value at key (insert as needed) */
	INLINE mapkit_error map_id_set(map_id * spm, int key, double value);

/* return the value at key (key must exists if alwaysdefaultvalue is not set) */
	INLINE double map_id_value(map_id * spm, int key);

/* get the value at key (return an error code) */
	INLINE mapkit_error map_id_get(map_id * spm, int key, double *value);

/* remove key (key must exists) */
	INLINE mapkit_error map_id_remove(map_id * spm, int key);

/*
 * return the next index to a full slot of the map, or -1 if none
 * if index == -1, returns the first index, or -1 if none
 * the contents of the slot can be accessed by spm->contents[index].key and
 * spm->contents[index].value
 */
	extern mapkit_size_t map_id_next(map_id * spm, mapkit_size_t iindex);

/*
 * return the next ptr to a full slot of the map, or NULL if none
 * if pos_contents == NULL, returns the first full slot, or NULL if none
 * the contents of the slot can be accessed by pos_contents->key and
 * pos_contents->value
 */
	extern map_id_storage *map_id_nextptr(map_id * spm, map_id_storage * pos_contents);

/* allocate an array of map_id_element with all "count" (key,value) pairs */
	extern mapkit_error map_id_getall(map_id * spm, map_id_element ** array, mapkit_size_t * count);

/* remove all values equal to the defaultvalue */
	extern mapkit_error map_id_clean(map_id * spm);

/* compare the key of two map_id_element (for qsort) */
	extern int map_id_compare(const void *e1, const void *e2);

/* allocate an array of map_id_element with all "count" (key,value) pairs */
	extern mapkit_error map_id_getall_sorted(map_id * spm, map_id_element ** array, mapkit_size_t * count);

/* insert all count map_id_element into spm */
	extern mapkit_error map_id_setall(map_id * spm, map_id_element * array, mapkit_size_t count);

/* remove all "count" elements whose keys are in "array" */
	extern mapkit_error map_id_removeall(map_id * spm, int *array, mapkit_size_t count);

/* print statistics */
	extern void map_id_printstats(map_id * spm);

/* Private functions */

	extern mapkit_error map_id_remove_s(map_id * spm, int key);
	extern double map_id_value_s(map_id * spm, int key);
	extern mapkit_error map_id_get_s(map_id * spm, int key, double *value);
	extern mapkit_error map_id_set_s(map_id * spm, int key, double value);
	extern double *map_id_insertptr_s(map_id * spm, int key);
	extern double *map_id_ptr_s(map_id * spm, int key);

/* Implementation */

/* INLINEd functions */
	mapkit_error map_id_remove(map_id * spm, int key) {
		map_id_storage *contents;

		contents = &(spm->contents[((mapkit_hash_t) key) % spm->size]);

		if ((contents->state == MAPKIT_FULLSLOT)
		    && ((contents->key) == (key))) {
			contents->state = MAPKIT_DELETEDSLOT;
#ifdef MAPKIT_COLLISIONS
			spm->keyindexs++;
#endif
			spm->used--;
			if (spm->used < spm->minused) {
#ifdef MAPKIT_DEBUG
				fprintf(stderr, "MAPKIT: used < minused\n");
#endif
				return map_id_reallocate(spm, map_id_shrinksize(spm, spm->used));
			}
		} else
			return map_id_remove_s(spm, key);
		return MAPKIT_OK;
	}

	mapkit_error map_id_removeptr(map_id * spm, double *ptr) {
		map_id_storage *sptr = (map_id_storage *) ((char *) ptr - offsetof(map_id_storage, value));

		sptr->state = MAPKIT_DELETEDSLOT;
		spm->used--;
		if (spm->used < spm->minused) {
#ifdef MAPKIT_DEBUG
			fprintf(stderr, "MAPKIT: used < minused\n");
#endif
			return map_id_reallocate(spm, map_id_shrinksize(spm, spm->used));
		}

		return MAPKIT_OK;
	}

	mapkit_error map_id_set(map_id * spm, int key, double value) {
		map_id_storage *contents;

		if (spm->alwaysdefault && ((value) == (spm->defaultvalue)))
			return map_id_remove(spm, key);

		contents = &(spm->contents[((mapkit_hash_t) key) % spm->size]);

		if (contents->state == MAPKIT_FREESLOT) {
			contents->state = MAPKIT_FULLSLOT;
			contents->key = key;
			contents->value = value;
#ifdef MAPKIT_COLLISIONS
			spm->insertionindexs++;
#endif
			spm->used++;
			spm->fill++;
			if (spm->fill > spm->maxfill) {
#ifdef MAPKIT_DEBUG
				fprintf(stderr, "MAPKIT: fill > maxfill\n");
#endif
				return map_id_reallocate(spm, map_id_growsize(spm, spm->used));
			}
			return MAPKIT_OK;
		} else if ((contents->state == MAPKIT_FULLSLOT)
			   && ((contents->key) == (key))) {
			contents->value = value;
#ifdef MAPKIT_COLLISIONS
			spm->insertionindexs++;
#endif
			return MAPKIT_OK;
		} else
			return map_id_set_s(spm, key, value);
	}

	double map_id_value(map_id * spm, int key) {
		map_id_storage *contents;

		contents = &(spm->contents[((mapkit_hash_t) key) % spm->size]);

		if ((contents->state == MAPKIT_FULLSLOT)
		    && ((contents->key) == (key))) {
#ifdef MAPKIT_COLLISIONS
			spm->keyindexs++;
#endif
			return contents->value;
		} else if (spm->alwaysdefault && (contents->state == MAPKIT_FREESLOT)) {
#ifdef MAPKIT_COLLISIONS
			spm->keyindexs++;
#endif
			return spm->defaultvalue;
		} else
			return map_id_value_s(spm, key);
	}

	mapkit_error map_id_get(map_id * spm, int key, double *value) {
		map_id_storage *contents;

		contents = &(spm->contents[((mapkit_hash_t) key) % spm->size]);

		if ((contents->state == MAPKIT_FULLSLOT)
		    && ((contents->key) == (key))) {
#ifdef MAPKIT_COLLISIONS
			spm->keyindexs++;
#endif
			*value = contents->value;
		} else if (spm->alwaysdefault && (contents->state == MAPKIT_FREESLOT)) {
#ifdef MAPKIT_COLLISIONS
			spm->keyindexs++;
#endif
			*value = spm->defaultvalue;
		} else
			return map_id_get_s(spm, key, value);

		return MAPKIT_OK;
	}

	double *map_id_insertptr(map_id * spm, int key) {
		map_id_storage *contents;

		contents = &(spm->contents[((mapkit_hash_t) key) % spm->size]);

		if (contents->state == MAPKIT_FREESLOT) {
			if (spm->fill >= spm->maxfill)
				return map_id_insertptr_s(spm, key);

			contents->state = MAPKIT_FULLSLOT;
			contents->key = key;
			contents->value = spm->defaultvalue;
#ifdef MAPKIT_COLLISIONS
			spm->insertionindexs++;
#endif
			spm->used++;
			spm->fill++;
			return &(contents->value);
		} else if ((contents->state == MAPKIT_FULLSLOT)
			   && ((contents->key) == (key))) {
#ifdef MAPKIT_COLLISIONS
			spm->insertionindexs++;
#endif
			return &(contents->value);
		} else
			return map_id_insertptr_s(spm, key);
	}

	double *map_id_ptr(map_id * spm, int key) {
		map_id_storage *contents;

		contents = &(spm->contents[((mapkit_hash_t) key) % spm->size]);

		if ((contents->state == MAPKIT_FULLSLOT)
		    && ((contents->key) == (key))) {
#ifdef MAPKIT_COLLISIONS
			spm->insertionindexs++;
#endif
			if ((!spm->alwaysdefault) || (!((contents->value) == (spm->defaultvalue))))
				return &(contents->value);
			else
				return NULL;
		} else if (contents->state == MAPKIT_FREESLOT) {
#ifdef MAPKIT_COLLISIONS
			spm->insertionindexs++;
#endif
			return NULL;
		} else
			return map_id_ptr_s(spm, key);
	}

/*
  Prototypes for map_ivp (int -> void*)
  Default value : NULL
  Uses a state field.
*/

	typedef struct {
		signed char state;
		int key;
		void *value;
	} map_ivp_storage;

	typedef struct {
		int key;
		void *value;
	} map_ivp_element;

	typedef struct {
		mapkit_size_t size,			       /* size of the hash table */
		 fill,					       /* # of non-free slots in the hash table */
		 used,					       /* # of full slots */
		 maxfill,				       /* max # of non-free slots before rehash */
		 minused /* min # of full slots before rehash */ ;

		double maxfillfactor,			       /* maxfill = maxfillfactor * size */
		 minusedfactor;				       /* minused = minusedfactor * size */

#ifdef MAPKIT_COLLISIONS
		unsigned long insertionindexs, insertionindex_collisions;
		unsigned long keyindexs, keyindex_collisions;
#endif

		void *defaultvalue;

		int alwaysdefault;

		map_ivp_storage *contents;
	} map_ivp;

/* Macros */

/* set to inform the type is available */
#define MAPKIT_map_ivp

/* read/write access to the value at key (insert as needed) */
#define map_ivp_(spm, key) (*(map_ivp_insertptr(spm, key)))

/* true if key exists */
#define map_ivp_haskey(spm, key) ((map_ivp_ptr(spm, key)) != NULL)

/* Prototypes */

/* initialize the structure */
	extern mapkit_error map_ivp_init(map_ivp * spm);

/* free contents */
	extern void map_ivp_free(map_ivp * spm);

/* copy the "from" map into a new "to" map */
	extern mapkit_error map_ivp_copy(map_ivp * to, map_ivp * from);

/* initialize the structure for "used" elements */
	extern mapkit_error map_ivp_init_hint(map_ivp * spm, mapkit_size_t used);

/* ensure at least "used" slots (including currently occupied slots) */
	extern mapkit_error map_ivp_ensurecapacity(map_ivp * spm, mapkit_size_t used);

/* restore maxfill & minused and ensure these constraints */
	extern mapkit_error map_ivp_adjustcapacity(map_ivp * spm);

/* reallocate contents for "size" slots, keeping previous values */
	extern mapkit_error map_ivp_reallocate(map_ivp * spm, mapkit_size_t size);

/* return an adequate size for growing from "used" elements */
	extern mapkit_size_t map_ivp_growsize(map_ivp * spm, mapkit_size_t used);

/* return an adequate size for shrinking from "used" elements */
	extern mapkit_size_t map_ivp_shrinksize(map_ivp * spm, mapkit_size_t used);

/* return an adequate size for "used" elements */
	extern mapkit_size_t map_ivp_meansize(map_ivp * spm, mapkit_size_t used);

/* pointer to the value at key (insert as needed) */
	INLINE void **map_ivp_insertptr(map_ivp * spm, int key);

/* pointer to the value at key or NULL */
	INLINE void **map_ivp_ptr(map_ivp * spm, int key);

/* remove the key pointing to the value *ptr. No checking is done */
	INLINE mapkit_error map_ivp_removeptr(map_ivp * spm, void **ptr);

/* set the value at key (insert as needed) */
	INLINE mapkit_error map_ivp_set(map_ivp * spm, int key, void *value);

/* return the value at key (key must exists if alwaysdefaultvalue is not set) */
	INLINE void *map_ivp_value(map_ivp * spm, int key);

/* get the value at key (return an error code) */
	INLINE mapkit_error map_ivp_get(map_ivp * spm, int key, void **value);

/* remove key (key must exists) */
	INLINE mapkit_error map_ivp_remove(map_ivp * spm, int key);

/*
 * return the next index to a full slot of the map, or -1 if none
 * if index == -1, returns the first index, or -1 if none
 * the contents of the slot can be accessed by spm->contents[index].key and
 * spm->contents[index].value
 */
	extern mapkit_size_t map_ivp_next(map_ivp * spm, mapkit_size_t iindex);

/*
 * return the next ptr to a full slot of the map, or NULL if none
 * if pos_contents == NULL, returns the first full slot, or NULL if none
 * the contents of the slot can be accessed by pos_contents->key and
 * pos_contents->value
 */
	extern map_ivp_storage *map_ivp_nextptr(map_ivp * spm, map_ivp_storage * pos_contents);

/* allocate an array of map_ivp_element with all "count" (key,value) pairs */
	extern mapkit_error map_ivp_getall(map_ivp * spm, map_ivp_element ** array, mapkit_size_t * count);

/* remove all values equal to the defaultvalue */
	extern mapkit_error map_ivp_clean(map_ivp * spm);

/* compare the key of two map_ivp_element (for qsort) */
	extern int map_ivp_compare(const void *e1, const void *e2);

/* allocate an array of map_ivp_element with all "count" (key,value) pairs */
	extern mapkit_error map_ivp_getall_sorted(map_ivp * spm, map_ivp_element ** array, mapkit_size_t * count);

/* insert all count map_ivp_element into spm */
	extern mapkit_error map_ivp_setall(map_ivp * spm, map_ivp_element * array, mapkit_size_t count);

/* remove all "count" elements whose keys are in "array" */
	extern mapkit_error map_ivp_removeall(map_ivp * spm, int *array, mapkit_size_t count);

/* print statistics */
	extern void map_ivp_printstats(map_ivp * spm);

/* Private functions */

	extern mapkit_error map_ivp_remove_s(map_ivp * spm, int key);
	extern void *map_ivp_value_s(map_ivp * spm, int key);
	extern mapkit_error map_ivp_get_s(map_ivp * spm, int key, void **value);
	extern mapkit_error map_ivp_set_s(map_ivp * spm, int key, void *value);
	extern void **map_ivp_insertptr_s(map_ivp * spm, int key);
	extern void **map_ivp_ptr_s(map_ivp * spm, int key);

/* Implementation */

/* INLINEd functions */
	mapkit_error map_ivp_remove(map_ivp * spm, int key) {
		map_ivp_storage *contents;

		contents = &(spm->contents[((mapkit_hash_t) key) % spm->size]);

		if ((contents->state == MAPKIT_FULLSLOT)
		    && ((contents->key) == (key))) {
			contents->state = MAPKIT_DELETEDSLOT;
#ifdef MAPKIT_COLLISIONS
			spm->keyindexs++;
#endif
			spm->used--;
			if (spm->used < spm->minused) {
#ifdef MAPKIT_DEBUG
				fprintf(stderr, "MAPKIT: used < minused\n");
#endif
				return map_ivp_reallocate(spm, map_ivp_shrinksize(spm, spm->used));
			}
		} else
			return map_ivp_remove_s(spm, key);
		return MAPKIT_OK;
	}

	mapkit_error map_ivp_removeptr(map_ivp * spm, void **ptr) {
		map_ivp_storage *sptr = (map_ivp_storage *) ((char *) ptr - offsetof(map_ivp_storage, value));

		sptr->state = MAPKIT_DELETEDSLOT;
		spm->used--;
		if (spm->used < spm->minused) {
#ifdef MAPKIT_DEBUG
			fprintf(stderr, "MAPKIT: used < minused\n");
#endif
			return map_ivp_reallocate(spm, map_ivp_shrinksize(spm, spm->used));
		}

		return MAPKIT_OK;
	}

	mapkit_error map_ivp_set(map_ivp * spm, int key, void *value) {
		map_ivp_storage *contents;

		if (spm->alwaysdefault && ((value) == (spm->defaultvalue)))
			return map_ivp_remove(spm, key);

		contents = &(spm->contents[((mapkit_hash_t) key) % spm->size]);

		if (contents->state == MAPKIT_FREESLOT) {
			contents->state = MAPKIT_FULLSLOT;
			contents->key = key;
			contents->value = value;
#ifdef MAPKIT_COLLISIONS
			spm->insertionindexs++;
#endif
			spm->used++;
			spm->fill++;
			if (spm->fill > spm->maxfill) {
#ifdef MAPKIT_DEBUG
				fprintf(stderr, "MAPKIT: fill > maxfill\n");
#endif
				return map_ivp_reallocate(spm, map_ivp_growsize(spm, spm->used));
			}
			return MAPKIT_OK;
		} else if ((contents->state == MAPKIT_FULLSLOT)
			   && ((contents->key) == (key))) {
			contents->value = value;
#ifdef MAPKIT_COLLISIONS
			spm->insertionindexs++;
#endif
			return MAPKIT_OK;
		} else
			return map_ivp_set_s(spm, key, value);
	}

	void *map_ivp_value(map_ivp * spm, int key) {
		map_ivp_storage *contents;

		contents = &(spm->contents[((mapkit_hash_t) key) % spm->size]);

		if ((contents->state == MAPKIT_FULLSLOT)
		    && ((contents->key) == (key))) {
#ifdef MAPKIT_COLLISIONS
			spm->keyindexs++;
#endif
			return contents->value;
		} else if (spm->alwaysdefault && (contents->state == MAPKIT_FREESLOT)) {
#ifdef MAPKIT_COLLISIONS
			spm->keyindexs++;
#endif
			return spm->defaultvalue;
		} else
			return map_ivp_value_s(spm, key);
	}

	mapkit_error map_ivp_get(map_ivp * spm, int key, void **value) {
		map_ivp_storage *contents;

		contents = &(spm->contents[((mapkit_hash_t) key) % spm->size]);

		if ((contents->state == MAPKIT_FULLSLOT)
		    && ((contents->key) == (key))) {
#ifdef MAPKIT_COLLISIONS
			spm->keyindexs++;
#endif
			*value = contents->value;
		} else if (spm->alwaysdefault && (contents->state == MAPKIT_FREESLOT)) {
#ifdef MAPKIT_COLLISIONS
			spm->keyindexs++;
#endif
			*value = spm->defaultvalue;
		} else
			return map_ivp_get_s(spm, key, value);

		return MAPKIT_OK;
	}

	void **map_ivp_insertptr(map_ivp * spm, int key) {
		map_ivp_storage *contents;

		contents = &(spm->contents[((mapkit_hash_t) key) % spm->size]);

		if (contents->state == MAPKIT_FREESLOT) {
			if (spm->fill >= spm->maxfill)
				return map_ivp_insertptr_s(spm, key);

			contents->state = MAPKIT_FULLSLOT;
			contents->key = key;
			contents->value = spm->defaultvalue;
#ifdef MAPKIT_COLLISIONS
			spm->insertionindexs++;
#endif
			spm->used++;
			spm->fill++;
			return &(contents->value);
		} else if ((contents->state == MAPKIT_FULLSLOT)
			   && ((contents->key) == (key))) {
#ifdef MAPKIT_COLLISIONS
			spm->insertionindexs++;
#endif
			return &(contents->value);
		} else
			return map_ivp_insertptr_s(spm, key);
	}

	void **map_ivp_ptr(map_ivp * spm, int key) {
		map_ivp_storage *contents;

		contents = &(spm->contents[((mapkit_hash_t) key) % spm->size]);

		if ((contents->state == MAPKIT_FULLSLOT)
		    && ((contents->key) == (key))) {
#ifdef MAPKIT_COLLISIONS
			spm->insertionindexs++;
#endif
			if ((!spm->alwaysdefault) || (!((contents->value) == (spm->defaultvalue))))
				return &(contents->value);
			else
				return NULL;
		} else if (contents->state == MAPKIT_FREESLOT) {
#ifdef MAPKIT_COLLISIONS
			spm->insertionindexs++;
#endif
			return NULL;
		} else
			return map_ivp_ptr_s(spm, key);
	}

/*
  Prototypes for map_h_ii (int -> int)
  Default value : 0
  Uses a state field.
*/

	typedef struct {
		signed char state;
		int key;
		int value;
	} map_h_ii_storage;

	typedef struct {
		int key;
		int value;
	} map_h_ii_element;

	typedef struct {
		mapkit_size_t size,			       /* size of the hash table */
		 fill,					       /* # of non-free slots in the hash table */
		 used,					       /* # of full slots */
		 maxfill,				       /* max # of non-free slots before rehash */
		 minused /* min # of full slots before rehash */ ;

		double maxfillfactor,			       /* maxfill = maxfillfactor * size */
		 minusedfactor;				       /* minused = minusedfactor * size */

#ifdef MAPKIT_COLLISIONS
		unsigned long insertionindexs, insertionindex_collisions;
		unsigned long keyindexs, keyindex_collisions;
#endif

		int defaultvalue;

		int alwaysdefault;

		map_h_ii_storage *contents;
	} map_h_ii;

/* Macros */

/* set to inform the type is available */
#define MAPKIT_map_h_ii

/* read/write access to the value at key (insert as needed) */
#define map_h_ii_(spm, key) (*(map_h_ii_insertptr(spm, key)))

/* true if key exists */
#define map_h_ii_haskey(spm, key) ((map_h_ii_ptr(spm, key)) != NULL)

/* Prototypes */

/* initialize the structure */
	extern mapkit_error map_h_ii_init(map_h_ii * spm);

/* free contents */
	extern void map_h_ii_free(map_h_ii * spm);

/* copy the "from" map into a new "to" map */
	extern mapkit_error map_h_ii_copy(map_h_ii * to, map_h_ii * from);

/* initialize the structure for "used" elements */
	extern mapkit_error map_h_ii_init_hint(map_h_ii * spm, mapkit_size_t used);

/* ensure at least "used" slots (including currently occupied slots) */
	extern mapkit_error map_h_ii_ensurecapacity(map_h_ii * spm, mapkit_size_t used);

/* restore maxfill & minused and ensure these constraints */
	extern mapkit_error map_h_ii_adjustcapacity(map_h_ii * spm);

/* reallocate contents for "size" slots, keeping previous values */
	extern mapkit_error map_h_ii_reallocate(map_h_ii * spm, mapkit_size_t size);

/* return an adequate size for growing from "used" elements */
	extern mapkit_size_t map_h_ii_growsize(map_h_ii * spm, mapkit_size_t used);

/* return an adequate size for shrinking from "used" elements */
	extern mapkit_size_t map_h_ii_shrinksize(map_h_ii * spm, mapkit_size_t used);

/* return an adequate size for "used" elements */
	extern mapkit_size_t map_h_ii_meansize(map_h_ii * spm, mapkit_size_t used);

/* pointer to the value at key (insert as needed) */
	INLINE int *map_h_ii_insertptr(map_h_ii * spm, int key);

/* pointer to the value at key or NULL */
	INLINE int *map_h_ii_ptr(map_h_ii * spm, int key);

/* remove the key pointing to the value *ptr. No checking is done */
	INLINE mapkit_error map_h_ii_removeptr(map_h_ii * spm, int *ptr);

/* set the value at key (insert as needed) */
	INLINE mapkit_error map_h_ii_set(map_h_ii * spm, int key, int value);

/* return the value at key (key must exists if alwaysdefaultvalue is not set) */
	INLINE int map_h_ii_value(map_h_ii * spm, int key);

/* get the value at key (return an error code) */
	INLINE mapkit_error map_h_ii_get(map_h_ii * spm, int key, int *value);

/* remove key (key must exists) */
	INLINE mapkit_error map_h_ii_remove(map_h_ii * spm, int key);

/*
 * return the next index to a full slot of the map, or -1 if none
 * if index == -1, returns the first index, or -1 if none
 * the contents of the slot can be accessed by spm->contents[index].key and
 * spm->contents[index].value
 */
	extern mapkit_size_t map_h_ii_next(map_h_ii * spm, mapkit_size_t iindex);

/*
 * return the next ptr to a full slot of the map, or NULL if none
 * if pos_contents == NULL, returns the first full slot, or NULL if none
 * the contents of the slot can be accessed by pos_contents->key and
 * pos_contents->value
 */
	extern map_h_ii_storage *map_h_ii_nextptr(map_h_ii * spm, map_h_ii_storage * pos_contents);

/* allocate an array of map_h_ii_element with all "count" (key,value) pairs */
	extern mapkit_error map_h_ii_getall(map_h_ii * spm, map_h_ii_element ** array, mapkit_size_t * count);

/* remove all values equal to the defaultvalue */
	extern mapkit_error map_h_ii_clean(map_h_ii * spm);

/* compare the key of two map_h_ii_element (for qsort) */
	extern int map_h_ii_compare(const void *e1, const void *e2);

/* allocate an array of map_h_ii_element with all "count" (key,value) pairs */
	extern mapkit_error map_h_ii_getall_sorted(map_h_ii * spm, map_h_ii_element ** array, mapkit_size_t * count);

/* insert all count map_h_ii_element into spm */
	extern mapkit_error map_h_ii_setall(map_h_ii * spm, map_h_ii_element * array, mapkit_size_t count);

/* remove all "count" elements whose keys are in "array" */
	extern mapkit_error map_h_ii_removeall(map_h_ii * spm, int *array, mapkit_size_t count);

/* print statistics */
	extern void map_h_ii_printstats(map_h_ii * spm);

/* Private functions */

	extern mapkit_error map_h_ii_remove_s(map_h_ii * spm, int key, mapkit_hash_t hash);
	extern int map_h_ii_value_s(map_h_ii * spm, int key, mapkit_hash_t hash);
	extern mapkit_error map_h_ii_get_s(map_h_ii * spm, int key, int *value, mapkit_hash_t hash);
	extern mapkit_error map_h_ii_set_s(map_h_ii * spm, int key, int value, mapkit_hash_t hash);
	extern int *map_h_ii_insertptr_s(map_h_ii * spm, int key, mapkit_hash_t hash);
	extern int *map_h_ii_ptr_s(map_h_ii * spm, int key, mapkit_hash_t hash);

/* Implementation */

/* INLINEd functions */
	mapkit_error map_h_ii_remove(map_h_ii * spm, int key) {
		mapkit_hash_t hash;
		map_h_ii_storage *contents;

		contents = &(spm->contents[(hash = mapkit_hash((mapkit_hash_t) key)) % spm->size]);

		if ((contents->state == MAPKIT_FULLSLOT)
		    && ((contents->key) == (key))) {
			contents->state = MAPKIT_DELETEDSLOT;
#ifdef MAPKIT_COLLISIONS
			spm->keyindexs++;
#endif
			spm->used--;
			if (spm->used < spm->minused) {
#ifdef MAPKIT_DEBUG
				fprintf(stderr, "MAPKIT: used < minused\n");
#endif
				return map_h_ii_reallocate(spm, map_h_ii_shrinksize(spm, spm->used));
			}
		} else
			return map_h_ii_remove_s(spm, key, hash);
		return MAPKIT_OK;
	}

	mapkit_error map_h_ii_removeptr(map_h_ii * spm, int *ptr) {
		map_h_ii_storage *sptr = (map_h_ii_storage *) ((char *) ptr - offsetof(map_h_ii_storage, value));

		sptr->state = MAPKIT_DELETEDSLOT;
		spm->used--;
		if (spm->used < spm->minused) {
#ifdef MAPKIT_DEBUG
			fprintf(stderr, "MAPKIT: used < minused\n");
#endif
			return map_h_ii_reallocate(spm, map_h_ii_shrinksize(spm, spm->used));
		}

		return MAPKIT_OK;
	}

	mapkit_error map_h_ii_set(map_h_ii * spm, int key, int value) {
		mapkit_hash_t hash;
		map_h_ii_storage *contents;

		if (spm->alwaysdefault && ((value) == (spm->defaultvalue)))
			return map_h_ii_remove(spm, key);

		contents = &(spm->contents[(hash = mapkit_hash((mapkit_hash_t) key)) % spm->size]);

		if (contents->state == MAPKIT_FREESLOT) {
			contents->state = MAPKIT_FULLSLOT;
			contents->key = key;
			contents->value = value;
#ifdef MAPKIT_COLLISIONS
			spm->insertionindexs++;
#endif
			spm->used++;
			spm->fill++;
			if (spm->fill > spm->maxfill) {
#ifdef MAPKIT_DEBUG
				fprintf(stderr, "MAPKIT: fill > maxfill\n");
#endif
				return map_h_ii_reallocate(spm, map_h_ii_growsize(spm, spm->used));
			}
			return MAPKIT_OK;
		} else if ((contents->state == MAPKIT_FULLSLOT)
			   && ((contents->key) == (key))) {
			contents->value = value;
#ifdef MAPKIT_COLLISIONS
			spm->insertionindexs++;
#endif
			return MAPKIT_OK;
		} else
			return map_h_ii_set_s(spm, key, value, hash);
	}

	int map_h_ii_value(map_h_ii * spm, int key) {
		mapkit_hash_t hash;
		map_h_ii_storage *contents;

		contents = &(spm->contents[(hash = mapkit_hash((mapkit_hash_t) key)) % spm->size]);

		if ((contents->state == MAPKIT_FULLSLOT)
		    && ((contents->key) == (key))) {
#ifdef MAPKIT_COLLISIONS
			spm->keyindexs++;
#endif
			return contents->value;
		} else if (spm->alwaysdefault && (contents->state == MAPKIT_FREESLOT)) {
#ifdef MAPKIT_COLLISIONS
			spm->keyindexs++;
#endif
			return spm->defaultvalue;
		} else
			return map_h_ii_value_s(spm, key, hash);
	}

	mapkit_error map_h_ii_get(map_h_ii * spm, int key, int *value) {
		mapkit_hash_t hash;
		map_h_ii_storage *contents;

		contents = &(spm->contents[(hash = mapkit_hash((mapkit_hash_t) key)) % spm->size]);

		if ((contents->state == MAPKIT_FULLSLOT)
		    && ((contents->key) == (key))) {
#ifdef MAPKIT_COLLISIONS
			spm->keyindexs++;
#endif
			*value = contents->value;
		} else if (spm->alwaysdefault && (contents->state == MAPKIT_FREESLOT)) {
#ifdef MAPKIT_COLLISIONS
			spm->keyindexs++;
#endif
			*value = spm->defaultvalue;
		} else
			return map_h_ii_get_s(spm, key, value, hash);

		return MAPKIT_OK;
	}

	int *map_h_ii_insertptr(map_h_ii * spm, int key) {
		mapkit_hash_t hash;
		map_h_ii_storage *contents;

		contents = &(spm->contents[(hash = mapkit_hash((mapkit_hash_t) key)) % spm->size]);

		if (contents->state == MAPKIT_FREESLOT) {
			if (spm->fill >= spm->maxfill)
				return map_h_ii_insertptr_s(spm, key, hash);

			contents->state = MAPKIT_FULLSLOT;
			contents->key = key;
			contents->value = spm->defaultvalue;
#ifdef MAPKIT_COLLISIONS
			spm->insertionindexs++;
#endif
			spm->used++;
			spm->fill++;
			return &(contents->value);
		} else if ((contents->state == MAPKIT_FULLSLOT)
			   && ((contents->key) == (key))) {
#ifdef MAPKIT_COLLISIONS
			spm->insertionindexs++;
#endif
			return &(contents->value);
		} else
			return map_h_ii_insertptr_s(spm, key, hash);
	}

	int *map_h_ii_ptr(map_h_ii * spm, int key) {
		mapkit_hash_t hash;
		map_h_ii_storage *contents;

		contents = &(spm->contents[(hash = mapkit_hash((mapkit_hash_t) key)) % spm->size]);

		if ((contents->state == MAPKIT_FULLSLOT)
		    && ((contents->key) == (key))) {
#ifdef MAPKIT_COLLISIONS
			spm->insertionindexs++;
#endif
			if ((!spm->alwaysdefault) || (!((contents->value) == (spm->defaultvalue))))
				return &(contents->value);
			else
				return NULL;
		} else if (contents->state == MAPKIT_FREESLOT) {
#ifdef MAPKIT_COLLISIONS
			spm->insertionindexs++;
#endif
			return NULL;
		} else
			return map_h_ii_ptr_s(spm, key, hash);
	}

/*
  Prototypes for map_h_id (int -> double)
  Default value : 0.0
  Uses a state field.
*/

	typedef struct {
		signed char state;
		int key;
		double value;
	} map_h_id_storage;

	typedef struct {
		int key;
		double value;
	} map_h_id_element;

	typedef struct {
		mapkit_size_t size,			       /* size of the hash table */
		 fill,					       /* # of non-free slots in the hash table */
		 used,					       /* # of full slots */
		 maxfill,				       /* max # of non-free slots before rehash */
		 minused /* min # of full slots before rehash */ ;

		double maxfillfactor,			       /* maxfill = maxfillfactor * size */
		 minusedfactor;				       /* minused = minusedfactor * size */

#ifdef MAPKIT_COLLISIONS
		unsigned long insertionindexs, insertionindex_collisions;
		unsigned long keyindexs, keyindex_collisions;
#endif

		double defaultvalue;

		int alwaysdefault;

		map_h_id_storage *contents;
	} map_h_id;

/* Macros */

/* set to inform the type is available */
#define MAPKIT_map_h_id

/* read/write access to the value at key (insert as needed) */
#define map_h_id_(spm, key) (*(map_h_id_insertptr(spm, key)))

/* true if key exists */
#define map_h_id_haskey(spm, key) ((map_h_id_ptr(spm, key)) != NULL)

/* Prototypes */

/* initialize the structure */
	extern mapkit_error map_h_id_init(map_h_id * spm);

/* free contents */
	extern void map_h_id_free(map_h_id * spm);

/* copy the "from" map into a new "to" map */
	extern mapkit_error map_h_id_copy(map_h_id * to, map_h_id * from);

/* initialize the structure for "used" elements */
	extern mapkit_error map_h_id_init_hint(map_h_id * spm, mapkit_size_t used);

/* ensure at least "used" slots (including currently occupied slots) */
	extern mapkit_error map_h_id_ensurecapacity(map_h_id * spm, mapkit_size_t used);

/* restore maxfill & minused and ensure these constraints */
	extern mapkit_error map_h_id_adjustcapacity(map_h_id * spm);

/* reallocate contents for "size" slots, keeping previous values */
	extern mapkit_error map_h_id_reallocate(map_h_id * spm, mapkit_size_t size);

/* return an adequate size for growing from "used" elements */
	extern mapkit_size_t map_h_id_growsize(map_h_id * spm, mapkit_size_t used);

/* return an adequate size for shrinking from "used" elements */
	extern mapkit_size_t map_h_id_shrinksize(map_h_id * spm, mapkit_size_t used);

/* return an adequate size for "used" elements */
	extern mapkit_size_t map_h_id_meansize(map_h_id * spm, mapkit_size_t used);

/* pointer to the value at key (insert as needed) */
	INLINE double *map_h_id_insertptr(map_h_id * spm, int key);

/* pointer to the value at key or NULL */
	INLINE double *map_h_id_ptr(map_h_id * spm, int key);

/* remove the key pointing to the value *ptr. No checking is done */
	INLINE mapkit_error map_h_id_removeptr(map_h_id * spm, double *ptr);

/* set the value at key (insert as needed) */
	INLINE mapkit_error map_h_id_set(map_h_id * spm, int key, double value);

/* return the value at key (key must exists if alwaysdefaultvalue is not set) */
	INLINE double map_h_id_value(map_h_id * spm, int key);

/* get the value at key (return an error code) */
	INLINE mapkit_error map_h_id_get(map_h_id * spm, int key, double *value);

/* remove key (key must exists) */
	INLINE mapkit_error map_h_id_remove(map_h_id * spm, int key);

/*
 * return the next index to a full slot of the map, or -1 if none
 * if index == -1, returns the first index, or -1 if none
 * the contents of the slot can be accessed by spm->contents[index].key and
 * spm->contents[index].value
 */
	extern mapkit_size_t map_h_id_next(map_h_id * spm, mapkit_size_t iindex);

/*
 * return the next ptr to a full slot of the map, or NULL if none
 * if pos_contents == NULL, returns the first full slot, or NULL if none
 * the contents of the slot can be accessed by pos_contents->key and
 * pos_contents->value
 */
	extern map_h_id_storage *map_h_id_nextptr(map_h_id * spm, map_h_id_storage * pos_contents);

/* allocate an array of map_h_id_element with all "count" (key,value) pairs */
	extern mapkit_error map_h_id_getall(map_h_id * spm, map_h_id_element ** array, mapkit_size_t * count);

/* remove all values equal to the defaultvalue */
	extern mapkit_error map_h_id_clean(map_h_id * spm);

/* compare the key of two map_h_id_element (for qsort) */
	extern int map_h_id_compare(const void *e1, const void *e2);

/* allocate an array of map_h_id_element with all "count" (key,value) pairs */
	extern mapkit_error map_h_id_getall_sorted(map_h_id * spm, map_h_id_element ** array, mapkit_size_t * count);

/* insert all count map_h_id_element into spm */
	extern mapkit_error map_h_id_setall(map_h_id * spm, map_h_id_element * array, mapkit_size_t count);

/* remove all "count" elements whose keys are in "array" */
	extern mapkit_error map_h_id_removeall(map_h_id * spm, int *array, mapkit_size_t count);

/* print statistics */
	extern void map_h_id_printstats(map_h_id * spm);

/* Private functions */

	extern mapkit_error map_h_id_remove_s(map_h_id * spm, int key, mapkit_hash_t hash);
	extern double map_h_id_value_s(map_h_id * spm, int key, mapkit_hash_t hash);
	extern mapkit_error map_h_id_get_s(map_h_id * spm, int key, double *value, mapkit_hash_t hash);
	extern mapkit_error map_h_id_set_s(map_h_id * spm, int key, double value, mapkit_hash_t hash);
	extern double *map_h_id_insertptr_s(map_h_id * spm, int key, mapkit_hash_t hash);
	extern double *map_h_id_ptr_s(map_h_id * spm, int key, mapkit_hash_t hash);

/* Implementation */

/* INLINEd functions */
	mapkit_error map_h_id_remove(map_h_id * spm, int key) {
		mapkit_hash_t hash;
		map_h_id_storage *contents;

		contents = &(spm->contents[(hash = mapkit_hash((mapkit_hash_t) key)) % spm->size]);

		if ((contents->state == MAPKIT_FULLSLOT)
		    && ((contents->key) == (key))) {
			contents->state = MAPKIT_DELETEDSLOT;
#ifdef MAPKIT_COLLISIONS
			spm->keyindexs++;
#endif
			spm->used--;
			if (spm->used < spm->minused) {
#ifdef MAPKIT_DEBUG
				fprintf(stderr, "MAPKIT: used < minused\n");
#endif
				return map_h_id_reallocate(spm, map_h_id_shrinksize(spm, spm->used));
			}
		} else
			return map_h_id_remove_s(spm, key, hash);
		return MAPKIT_OK;
	}

	mapkit_error map_h_id_removeptr(map_h_id * spm, double *ptr) {
		map_h_id_storage *sptr = (map_h_id_storage *) ((char *) ptr - offsetof(map_h_id_storage, value));

		sptr->state = MAPKIT_DELETEDSLOT;
		spm->used--;
		if (spm->used < spm->minused) {
#ifdef MAPKIT_DEBUG
			fprintf(stderr, "MAPKIT: used < minused\n");
#endif
			return map_h_id_reallocate(spm, map_h_id_shrinksize(spm, spm->used));
		}

		return MAPKIT_OK;
	}

	mapkit_error map_h_id_set(map_h_id * spm, int key, double value) {
		mapkit_hash_t hash;
		map_h_id_storage *contents;

		if (spm->alwaysdefault && ((value) == (spm->defaultvalue)))
			return map_h_id_remove(spm, key);

		contents = &(spm->contents[(hash = mapkit_hash((mapkit_hash_t) key)) % spm->size]);

		if (contents->state == MAPKIT_FREESLOT) {
			contents->state = MAPKIT_FULLSLOT;
			contents->key = key;
			contents->value = value;
#ifdef MAPKIT_COLLISIONS
			spm->insertionindexs++;
#endif
			spm->used++;
			spm->fill++;
			if (spm->fill > spm->maxfill) {
#ifdef MAPKIT_DEBUG
				fprintf(stderr, "MAPKIT: fill > maxfill\n");
#endif
				return map_h_id_reallocate(spm, map_h_id_growsize(spm, spm->used));
			}
			return MAPKIT_OK;
		} else if ((contents->state == MAPKIT_FULLSLOT)
			   && ((contents->key) == (key))) {
			contents->value = value;
#ifdef MAPKIT_COLLISIONS
			spm->insertionindexs++;
#endif
			return MAPKIT_OK;
		} else
			return map_h_id_set_s(spm, key, value, hash);
	}

	double map_h_id_value(map_h_id * spm, int key) {
		mapkit_hash_t hash;
		map_h_id_storage *contents;

		contents = &(spm->contents[(hash = mapkit_hash((mapkit_hash_t) key)) % spm->size]);

		if ((contents->state == MAPKIT_FULLSLOT)
		    && ((contents->key) == (key))) {
#ifdef MAPKIT_COLLISIONS
			spm->keyindexs++;
#endif
			return contents->value;
		} else if (spm->alwaysdefault && (contents->state == MAPKIT_FREESLOT)) {
#ifdef MAPKIT_COLLISIONS
			spm->keyindexs++;
#endif
			return spm->defaultvalue;
		} else
			return map_h_id_value_s(spm, key, hash);
	}

	mapkit_error map_h_id_get(map_h_id * spm, int key, double *value) {
		mapkit_hash_t hash;
		map_h_id_storage *contents;

		contents = &(spm->contents[(hash = mapkit_hash((mapkit_hash_t) key)) % spm->size]);

		if ((contents->state == MAPKIT_FULLSLOT)
		    && ((contents->key) == (key))) {
#ifdef MAPKIT_COLLISIONS
			spm->keyindexs++;
#endif
			*value = contents->value;
		} else if (spm->alwaysdefault && (contents->state == MAPKIT_FREESLOT)) {
#ifdef MAPKIT_COLLISIONS
			spm->keyindexs++;
#endif
			*value = spm->defaultvalue;
		} else
			return map_h_id_get_s(spm, key, value, hash);

		return MAPKIT_OK;
	}

	double *map_h_id_insertptr(map_h_id * spm, int key) {
		mapkit_hash_t hash;
		map_h_id_storage *contents;

		contents = &(spm->contents[(hash = mapkit_hash((mapkit_hash_t) key)) % spm->size]);

		if (contents->state == MAPKIT_FREESLOT) {
			if (spm->fill >= spm->maxfill)
				return map_h_id_insertptr_s(spm, key, hash);

			contents->state = MAPKIT_FULLSLOT;
			contents->key = key;
			contents->value = spm->defaultvalue;
#ifdef MAPKIT_COLLISIONS
			spm->insertionindexs++;
#endif
			spm->used++;
			spm->fill++;
			return &(contents->value);
		} else if ((contents->state == MAPKIT_FULLSLOT)
			   && ((contents->key) == (key))) {
#ifdef MAPKIT_COLLISIONS
			spm->insertionindexs++;
#endif
			return &(contents->value);
		} else
			return map_h_id_insertptr_s(spm, key, hash);
	}

	double *map_h_id_ptr(map_h_id * spm, int key) {
		mapkit_hash_t hash;
		map_h_id_storage *contents;

		contents = &(spm->contents[(hash = mapkit_hash((mapkit_hash_t) key)) % spm->size]);

		if ((contents->state == MAPKIT_FULLSLOT)
		    && ((contents->key) == (key))) {
#ifdef MAPKIT_COLLISIONS
			spm->insertionindexs++;
#endif
			if ((!spm->alwaysdefault) || (!((contents->value) == (spm->defaultvalue))))
				return &(contents->value);
			else
				return NULL;
		} else if (contents->state == MAPKIT_FREESLOT) {
#ifdef MAPKIT_COLLISIONS
			spm->insertionindexs++;
#endif
			return NULL;
		} else
			return map_h_id_ptr_s(spm, key, hash);
	}

/*
  Prototypes for map_h_ivp (int -> void*)
  Default value : NULL
  Uses a state field.
*/

	typedef struct {
		signed char state;
		int key;
		void *value;
	} map_h_ivp_storage;

	typedef struct {
		int key;
		void *value;
	} map_h_ivp_element;

	typedef struct {
		mapkit_size_t size,			       /* size of the hash table */
		 fill,					       /* # of non-free slots in the hash table */
		 used,					       /* # of full slots */
		 maxfill,				       /* max # of non-free slots before rehash */
		 minused /* min # of full slots before rehash */ ;

		double maxfillfactor,			       /* maxfill = maxfillfactor * size */
		 minusedfactor;				       /* minused = minusedfactor * size */

#ifdef MAPKIT_COLLISIONS
		unsigned long insertionindexs, insertionindex_collisions;
		unsigned long keyindexs, keyindex_collisions;
#endif

		void *defaultvalue;

		int alwaysdefault;

		map_h_ivp_storage *contents;
	} map_h_ivp;

/* Macros */

/* set to inform the type is available */
#define MAPKIT_map_h_ivp

/* read/write access to the value at key (insert as needed) */
#define map_h_ivp_(spm, key) (*(map_h_ivp_insertptr(spm, key)))

/* true if key exists */
#define map_h_ivp_haskey(spm, key) ((map_h_ivp_ptr(spm, key)) != NULL)

/* Prototypes */

/* initialize the structure */
	extern mapkit_error map_h_ivp_init(map_h_ivp * spm);

/* free contents */
	extern void map_h_ivp_free(map_h_ivp * spm);

/* copy the "from" map into a new "to" map */
	extern mapkit_error map_h_ivp_copy(map_h_ivp * to, map_h_ivp * from);

/* initialize the structure for "used" elements */
	extern mapkit_error map_h_ivp_init_hint(map_h_ivp * spm, mapkit_size_t used);

/* ensure at least "used" slots (including currently occupied slots) */
	extern mapkit_error map_h_ivp_ensurecapacity(map_h_ivp * spm, mapkit_size_t used);

/* restore maxfill & minused and ensure these constraints */
	extern mapkit_error map_h_ivp_adjustcapacity(map_h_ivp * spm);

/* reallocate contents for "size" slots, keeping previous values */
	extern mapkit_error map_h_ivp_reallocate(map_h_ivp * spm, mapkit_size_t size);

/* return an adequate size for growing from "used" elements */
	extern mapkit_size_t map_h_ivp_growsize(map_h_ivp * spm, mapkit_size_t used);

/* return an adequate size for shrinking from "used" elements */
	extern mapkit_size_t map_h_ivp_shrinksize(map_h_ivp * spm, mapkit_size_t used);

/* return an adequate size for "used" elements */
	extern mapkit_size_t map_h_ivp_meansize(map_h_ivp * spm, mapkit_size_t used);

/* pointer to the value at key (insert as needed) */
	INLINE void **map_h_ivp_insertptr(map_h_ivp * spm, int key);

/* pointer to the value at key or NULL */
	INLINE void **map_h_ivp_ptr(map_h_ivp * spm, int key);

/* remove the key pointing to the value *ptr. No checking is done */
	INLINE mapkit_error map_h_ivp_removeptr(map_h_ivp * spm, void **ptr);

/* set the value at key (insert as needed) */
	INLINE mapkit_error map_h_ivp_set(map_h_ivp * spm, int key, void *value);

/* return the value at key (key must exists if alwaysdefaultvalue is not set) */
	INLINE void *map_h_ivp_value(map_h_ivp * spm, int key);

/* get the value at key (return an error code) */
	INLINE mapkit_error map_h_ivp_get(map_h_ivp * spm, int key, void **value);

/* remove key (key must exists) */
	INLINE mapkit_error map_h_ivp_remove(map_h_ivp * spm, int key);

/*
 * return the next index to a full slot of the map, or -1 if none
 * if index == -1, returns the first index, or -1 if none
 * the contents of the slot can be accessed by spm->contents[index].key and
 * spm->contents[index].value
 */
	extern mapkit_size_t map_h_ivp_next(map_h_ivp * spm, mapkit_size_t iindex);

/*
 * return the next ptr to a full slot of the map, or NULL if none
 * if pos_contents == NULL, returns the first full slot, or NULL if none
 * the contents of the slot can be accessed by pos_contents->key and
 * pos_contents->value
 */
	extern map_h_ivp_storage *map_h_ivp_nextptr(map_h_ivp * spm, map_h_ivp_storage * pos_contents);

/* allocate an array of map_h_ivp_element with all "count" (key,value) pairs */
	extern mapkit_error map_h_ivp_getall(map_h_ivp * spm, map_h_ivp_element ** array, mapkit_size_t * count);

/* remove all values equal to the defaultvalue */
	extern mapkit_error map_h_ivp_clean(map_h_ivp * spm);

/* compare the key of two map_h_ivp_element (for qsort) */
	extern int map_h_ivp_compare(const void *e1, const void *e2);

/* allocate an array of map_h_ivp_element with all "count" (key,value) pairs */
	extern mapkit_error map_h_ivp_getall_sorted(map_h_ivp * spm, map_h_ivp_element ** array, mapkit_size_t * count);

/* insert all count map_h_ivp_element into spm */
	extern mapkit_error map_h_ivp_setall(map_h_ivp * spm, map_h_ivp_element * array, mapkit_size_t count);

/* remove all "count" elements whose keys are in "array" */
	extern mapkit_error map_h_ivp_removeall(map_h_ivp * spm, int *array, mapkit_size_t count);

/* print statistics */
	extern void map_h_ivp_printstats(map_h_ivp * spm);

/* Private functions */

	extern mapkit_error map_h_ivp_remove_s(map_h_ivp * spm, int key, mapkit_hash_t hash);
	extern void *map_h_ivp_value_s(map_h_ivp * spm, int key, mapkit_hash_t hash);
	extern mapkit_error map_h_ivp_get_s(map_h_ivp * spm, int key, void **value, mapkit_hash_t hash);
	extern mapkit_error map_h_ivp_set_s(map_h_ivp * spm, int key, void *value, mapkit_hash_t hash);
	extern void **map_h_ivp_insertptr_s(map_h_ivp * spm, int key, mapkit_hash_t hash);
	extern void **map_h_ivp_ptr_s(map_h_ivp * spm, int key, mapkit_hash_t hash);

/* Implementation */

/* INLINEd functions */
	mapkit_error map_h_ivp_remove(map_h_ivp * spm, int key) {
		mapkit_hash_t hash;
		map_h_ivp_storage *contents;

		contents = &(spm->contents[(hash = mapkit_hash((mapkit_hash_t) key)) % spm->size]);

		if ((contents->state == MAPKIT_FULLSLOT)
		    && ((contents->key) == (key))) {
			contents->state = MAPKIT_DELETEDSLOT;
#ifdef MAPKIT_COLLISIONS
			spm->keyindexs++;
#endif
			spm->used--;
			if (spm->used < spm->minused) {
#ifdef MAPKIT_DEBUG
				fprintf(stderr, "MAPKIT: used < minused\n");
#endif
				return map_h_ivp_reallocate(spm, map_h_ivp_shrinksize(spm, spm->used));
			}
		} else
			return map_h_ivp_remove_s(spm, key, hash);
		return MAPKIT_OK;
	}

	mapkit_error map_h_ivp_removeptr(map_h_ivp * spm, void **ptr) {
		map_h_ivp_storage *sptr = (map_h_ivp_storage *) ((char *) ptr - offsetof(map_h_ivp_storage, value));

		sptr->state = MAPKIT_DELETEDSLOT;
		spm->used--;
		if (spm->used < spm->minused) {
#ifdef MAPKIT_DEBUG
			fprintf(stderr, "MAPKIT: used < minused\n");
#endif
			return map_h_ivp_reallocate(spm, map_h_ivp_shrinksize(spm, spm->used));
		}

		return MAPKIT_OK;
	}

	mapkit_error map_h_ivp_set(map_h_ivp * spm, int key, void *value) {
		mapkit_hash_t hash;
		map_h_ivp_storage *contents;

		if (spm->alwaysdefault && ((value) == (spm->defaultvalue)))
			return map_h_ivp_remove(spm, key);

		contents = &(spm->contents[(hash = mapkit_hash((mapkit_hash_t) key)) % spm->size]);

		if (contents->state == MAPKIT_FREESLOT) {
			contents->state = MAPKIT_FULLSLOT;
			contents->key = key;
			contents->value = value;
#ifdef MAPKIT_COLLISIONS
			spm->insertionindexs++;
#endif
			spm->used++;
			spm->fill++;
			if (spm->fill > spm->maxfill) {
#ifdef MAPKIT_DEBUG
				fprintf(stderr, "MAPKIT: fill > maxfill\n");
#endif
				return map_h_ivp_reallocate(spm, map_h_ivp_growsize(spm, spm->used));
			}
			return MAPKIT_OK;
		} else if ((contents->state == MAPKIT_FULLSLOT)
			   && ((contents->key) == (key))) {
			contents->value = value;
#ifdef MAPKIT_COLLISIONS
			spm->insertionindexs++;
#endif
			return MAPKIT_OK;
		} else
			return map_h_ivp_set_s(spm, key, value, hash);
	}

	void *map_h_ivp_value(map_h_ivp * spm, int key) {
		mapkit_hash_t hash;
		map_h_ivp_storage *contents;

		contents = &(spm->contents[(hash = mapkit_hash((mapkit_hash_t) key)) % spm->size]);

		if ((contents->state == MAPKIT_FULLSLOT)
		    && ((contents->key) == (key))) {
#ifdef MAPKIT_COLLISIONS
			spm->keyindexs++;
#endif
			return contents->value;
		} else if (spm->alwaysdefault && (contents->state == MAPKIT_FREESLOT)) {
#ifdef MAPKIT_COLLISIONS
			spm->keyindexs++;
#endif
			return spm->defaultvalue;
		} else
			return map_h_ivp_value_s(spm, key, hash);
	}

	mapkit_error map_h_ivp_get(map_h_ivp * spm, int key, void **value) {
		mapkit_hash_t hash;
		map_h_ivp_storage *contents;

		contents = &(spm->contents[(hash = mapkit_hash((mapkit_hash_t) key)) % spm->size]);

		if ((contents->state == MAPKIT_FULLSLOT)
		    && ((contents->key) == (key))) {
#ifdef MAPKIT_COLLISIONS
			spm->keyindexs++;
#endif
			*value = contents->value;
		} else if (spm->alwaysdefault && (contents->state == MAPKIT_FREESLOT)) {
#ifdef MAPKIT_COLLISIONS
			spm->keyindexs++;
#endif
			*value = spm->defaultvalue;
		} else
			return map_h_ivp_get_s(spm, key, value, hash);

		return MAPKIT_OK;
	}

	void **map_h_ivp_insertptr(map_h_ivp * spm, int key) {
		mapkit_hash_t hash;
		map_h_ivp_storage *contents;

		contents = &(spm->contents[(hash = mapkit_hash((mapkit_hash_t) key)) % spm->size]);

		if (contents->state == MAPKIT_FREESLOT) {
			if (spm->fill >= spm->maxfill)
				return map_h_ivp_insertptr_s(spm, key, hash);

			contents->state = MAPKIT_FULLSLOT;
			contents->key = key;
			contents->value = spm->defaultvalue;
#ifdef MAPKIT_COLLISIONS
			spm->insertionindexs++;
#endif
			spm->used++;
			spm->fill++;
			return &(contents->value);
		} else if ((contents->state == MAPKIT_FULLSLOT)
			   && ((contents->key) == (key))) {
#ifdef MAPKIT_COLLISIONS
			spm->insertionindexs++;
#endif
			return &(contents->value);
		} else
			return map_h_ivp_insertptr_s(spm, key, hash);
	}

	void **map_h_ivp_ptr(map_h_ivp * spm, int key) {
		mapkit_hash_t hash;
		map_h_ivp_storage *contents;

		contents = &(spm->contents[(hash = mapkit_hash((mapkit_hash_t) key)) % spm->size]);

		if ((contents->state == MAPKIT_FULLSLOT)
		    && ((contents->key) == (key))) {
#ifdef MAPKIT_COLLISIONS
			spm->insertionindexs++;
#endif
			if ((!spm->alwaysdefault) || (!((contents->value) == (spm->defaultvalue))))
				return &(contents->value);
			else
				return NULL;
		} else if (contents->state == MAPKIT_FREESLOT) {
#ifdef MAPKIT_COLLISIONS
			spm->insertionindexs++;
#endif
			return NULL;
		} else
			return map_h_ivp_ptr_s(spm, key, hash);
	}

/*
  Prototypes for map_vpi (void* -> int)
  Default value : 0
  Uses a state field.
*/

	typedef struct {
		signed char state;
		void *key;
		int value;
	} map_vpi_storage;

	typedef struct {
		void *key;
		int value;
	} map_vpi_element;

	typedef struct {
		mapkit_size_t size,			       /* size of the hash table */
		 fill,					       /* # of non-free slots in the hash table */
		 used,					       /* # of full slots */
		 maxfill,				       /* max # of non-free slots before rehash */
		 minused /* min # of full slots before rehash */ ;

		double maxfillfactor,			       /* maxfill = maxfillfactor * size */
		 minusedfactor;				       /* minused = minusedfactor * size */

#ifdef MAPKIT_COLLISIONS
		unsigned long insertionindexs, insertionindex_collisions;
		unsigned long keyindexs, keyindex_collisions;
#endif

		int defaultvalue;

		int alwaysdefault;

		map_vpi_storage *contents;
	} map_vpi;

/* Macros */

/* set to inform the type is available */
#define MAPKIT_map_vpi

/* read/write access to the value at key (insert as needed) */
#define map_vpi_(spm, key) (*(map_vpi_insertptr(spm, key)))

/* true if key exists */
#define map_vpi_haskey(spm, key) ((map_vpi_ptr(spm, key)) != NULL)

/* Prototypes */

/* initialize the structure */
	extern mapkit_error map_vpi_init(map_vpi * spm);

/* free contents */
	extern void map_vpi_free(map_vpi * spm);

/* copy the "from" map into a new "to" map */
	extern mapkit_error map_vpi_copy(map_vpi * to, map_vpi * from);

/* initialize the structure for "used" elements */
	extern mapkit_error map_vpi_init_hint(map_vpi * spm, mapkit_size_t used);

/* ensure at least "used" slots (including currently occupied slots) */
	extern mapkit_error map_vpi_ensurecapacity(map_vpi * spm, mapkit_size_t used);

/* restore maxfill & minused and ensure these constraints */
	extern mapkit_error map_vpi_adjustcapacity(map_vpi * spm);

/* reallocate contents for "size" slots, keeping previous values */
	extern mapkit_error map_vpi_reallocate(map_vpi * spm, mapkit_size_t size);

/* return an adequate size for growing from "used" elements */
	extern mapkit_size_t map_vpi_growsize(map_vpi * spm, mapkit_size_t used);

/* return an adequate size for shrinking from "used" elements */
	extern mapkit_size_t map_vpi_shrinksize(map_vpi * spm, mapkit_size_t used);

/* return an adequate size for "used" elements */
	extern mapkit_size_t map_vpi_meansize(map_vpi * spm, mapkit_size_t used);

/* pointer to the value at key (insert as needed) */
	INLINE int *map_vpi_insertptr(map_vpi * spm, void *key);

/* pointer to the value at key or NULL */
	INLINE int *map_vpi_ptr(map_vpi * spm, void *key);

/* remove the key pointing to the value *ptr. No checking is done */
	INLINE mapkit_error map_vpi_removeptr(map_vpi * spm, int *ptr);

/* set the value at key (insert as needed) */
	INLINE mapkit_error map_vpi_set(map_vpi * spm, void *key, int value);

/* return the value at key (key must exists if alwaysdefaultvalue is not set) */
	INLINE int map_vpi_value(map_vpi * spm, void *key);

/* get the value at key (return an error code) */
	INLINE mapkit_error map_vpi_get(map_vpi * spm, void *key, int *value);

/* remove key (key must exists) */
	INLINE mapkit_error map_vpi_remove(map_vpi * spm, void *key);

/*
 * return the next index to a full slot of the map, or -1 if none
 * if index == -1, returns the first index, or -1 if none
 * the contents of the slot can be accessed by spm->contents[index].key and
 * spm->contents[index].value
 */
	extern mapkit_size_t map_vpi_next(map_vpi * spm, mapkit_size_t iindex);

/*
 * return the next ptr to a full slot of the map, or NULL if none
 * if pos_contents == NULL, returns the first full slot, or NULL if none
 * the contents of the slot can be accessed by pos_contents->key and
 * pos_contents->value
 */
	extern map_vpi_storage *map_vpi_nextptr(map_vpi * spm, map_vpi_storage * pos_contents);

/* allocate an array of map_vpi_element with all "count" (key,value) pairs */
	extern mapkit_error map_vpi_getall(map_vpi * spm, map_vpi_element ** array, mapkit_size_t * count);

/* remove all values equal to the defaultvalue */
	extern mapkit_error map_vpi_clean(map_vpi * spm);

/* compare the key of two map_vpi_element (for qsort) */
	extern int map_vpi_compare(const void *e1, const void *e2);

/* allocate an array of map_vpi_element with all "count" (key,value) pairs */
	extern mapkit_error map_vpi_getall_sorted(map_vpi * spm, map_vpi_element ** array, mapkit_size_t * count);

/* insert all count map_vpi_element into spm */
	extern mapkit_error map_vpi_setall(map_vpi * spm, map_vpi_element * array, mapkit_size_t count);

/* remove all "count" elements whose keys are in "array" */
	extern mapkit_error map_vpi_removeall(map_vpi * spm, void **array, mapkit_size_t count);

/* print statistics */
	extern void map_vpi_printstats(map_vpi * spm);

/* Private functions */

	extern mapkit_error map_vpi_remove_s(map_vpi * spm, void *key);
	extern int map_vpi_value_s(map_vpi * spm, void *key);
	extern mapkit_error map_vpi_get_s(map_vpi * spm, void *key, int *value);
	extern mapkit_error map_vpi_set_s(map_vpi * spm, void *key, int value);
	extern int *map_vpi_insertptr_s(map_vpi * spm, void *key);
	extern int *map_vpi_ptr_s(map_vpi * spm, void *key);

/* Implementation */

/* INLINEd functions */
	mapkit_error map_vpi_remove(map_vpi * spm, void *key) {
		map_vpi_storage *contents;

		contents = &(spm->contents[((mapkit_hash_t) key) % spm->size]);

		if ((contents->state == MAPKIT_FULLSLOT)
		    && ((contents->key) == (key))) {
			contents->state = MAPKIT_DELETEDSLOT;
#ifdef MAPKIT_COLLISIONS
			spm->keyindexs++;
#endif
			spm->used--;
			if (spm->used < spm->minused) {
#ifdef MAPKIT_DEBUG
				fprintf(stderr, "MAPKIT: used < minused\n");
#endif
				return map_vpi_reallocate(spm, map_vpi_shrinksize(spm, spm->used));
			}
		} else
			return map_vpi_remove_s(spm, key);
		return MAPKIT_OK;
	}

	mapkit_error map_vpi_removeptr(map_vpi * spm, int *ptr) {
		map_vpi_storage *sptr = (map_vpi_storage *) ((char *) ptr - offsetof(map_vpi_storage, value));

		sptr->state = MAPKIT_DELETEDSLOT;
		spm->used--;
		if (spm->used < spm->minused) {
#ifdef MAPKIT_DEBUG
			fprintf(stderr, "MAPKIT: used < minused\n");
#endif
			return map_vpi_reallocate(spm, map_vpi_shrinksize(spm, spm->used));
		}

		return MAPKIT_OK;
	}

	mapkit_error map_vpi_set(map_vpi * spm, void *key, int value) {
		map_vpi_storage *contents;

		if (spm->alwaysdefault && ((value) == (spm->defaultvalue)))
			return map_vpi_remove(spm, key);

		contents = &(spm->contents[((mapkit_hash_t) key) % spm->size]);

		if (contents->state == MAPKIT_FREESLOT) {
			contents->state = MAPKIT_FULLSLOT;
			contents->key = key;
			contents->value = value;
#ifdef MAPKIT_COLLISIONS
			spm->insertionindexs++;
#endif
			spm->used++;
			spm->fill++;
			if (spm->fill > spm->maxfill) {
#ifdef MAPKIT_DEBUG
				fprintf(stderr, "MAPKIT: fill > maxfill\n");
#endif
				return map_vpi_reallocate(spm, map_vpi_growsize(spm, spm->used));
			}
			return MAPKIT_OK;
		} else if ((contents->state == MAPKIT_FULLSLOT)
			   && ((contents->key) == (key))) {
			contents->value = value;
#ifdef MAPKIT_COLLISIONS
			spm->insertionindexs++;
#endif
			return MAPKIT_OK;
		} else
			return map_vpi_set_s(spm, key, value);
	}

	int map_vpi_value(map_vpi * spm, void *key) {
		map_vpi_storage *contents;

		contents = &(spm->contents[((mapkit_hash_t) key) % spm->size]);

		if ((contents->state == MAPKIT_FULLSLOT)
		    && ((contents->key) == (key))) {
#ifdef MAPKIT_COLLISIONS
			spm->keyindexs++;
#endif
			return contents->value;
		} else if (spm->alwaysdefault && (contents->state == MAPKIT_FREESLOT)) {
#ifdef MAPKIT_COLLISIONS
			spm->keyindexs++;
#endif
			return spm->defaultvalue;
		} else
			return map_vpi_value_s(spm, key);
	}

	mapkit_error map_vpi_get(map_vpi * spm, void *key, int *value) {
		map_vpi_storage *contents;

		contents = &(spm->contents[((mapkit_hash_t) key) % spm->size]);

		if ((contents->state == MAPKIT_FULLSLOT)
		    && ((contents->key) == (key))) {
#ifdef MAPKIT_COLLISIONS
			spm->keyindexs++;
#endif
			*value = contents->value;
		} else if (spm->alwaysdefault && (contents->state == MAPKIT_FREESLOT)) {
#ifdef MAPKIT_COLLISIONS
			spm->keyindexs++;
#endif
			*value = spm->defaultvalue;
		} else
			return map_vpi_get_s(spm, key, value);

		return MAPKIT_OK;
	}

	int *map_vpi_insertptr(map_vpi * spm, void *key) {
		map_vpi_storage *contents;

		contents = &(spm->contents[((mapkit_hash_t) key) % spm->size]);

		if (contents->state == MAPKIT_FREESLOT) {
			if (spm->fill >= spm->maxfill)
				return map_vpi_insertptr_s(spm, key);

			contents->state = MAPKIT_FULLSLOT;
			contents->key = key;
			contents->value = spm->defaultvalue;
#ifdef MAPKIT_COLLISIONS
			spm->insertionindexs++;
#endif
			spm->used++;
			spm->fill++;
			return &(contents->value);
		} else if ((contents->state == MAPKIT_FULLSLOT)
			   && ((contents->key) == (key))) {
#ifdef MAPKIT_COLLISIONS
			spm->insertionindexs++;
#endif
			return &(contents->value);
		} else
			return map_vpi_insertptr_s(spm, key);
	}

	int *map_vpi_ptr(map_vpi * spm, void *key) {
		map_vpi_storage *contents;

		contents = &(spm->contents[((mapkit_hash_t) key) % spm->size]);

		if ((contents->state == MAPKIT_FULLSLOT)
		    && ((contents->key) == (key))) {
#ifdef MAPKIT_COLLISIONS
			spm->insertionindexs++;
#endif
			if ((!spm->alwaysdefault) || (!((contents->value) == (spm->defaultvalue))))
				return &(contents->value);
			else
				return NULL;
		} else if (contents->state == MAPKIT_FREESLOT) {
#ifdef MAPKIT_COLLISIONS
			spm->insertionindexs++;
#endif
			return NULL;
		} else
			return map_vpi_ptr_s(spm, key);
	}

/*
  Prototypes for map_vpd (void* -> double)
  Default value : 0.0
  Uses a state field.
*/

	typedef struct {
		signed char state;
		void *key;
		double value;
	} map_vpd_storage;

	typedef struct {
		void *key;
		double value;
	} map_vpd_element;

	typedef struct {
		mapkit_size_t size,			       /* size of the hash table */
		 fill,					       /* # of non-free slots in the hash table */
		 used,					       /* # of full slots */
		 maxfill,				       /* max # of non-free slots before rehash */
		 minused /* min # of full slots before rehash */ ;

		double maxfillfactor,			       /* maxfill = maxfillfactor * size */
		 minusedfactor;				       /* minused = minusedfactor * size */

#ifdef MAPKIT_COLLISIONS
		unsigned long insertionindexs, insertionindex_collisions;
		unsigned long keyindexs, keyindex_collisions;
#endif

		double defaultvalue;

		int alwaysdefault;

		map_vpd_storage *contents;
	} map_vpd;

/* Macros */

/* set to inform the type is available */
#define MAPKIT_map_vpd

/* read/write access to the value at key (insert as needed) */
#define map_vpd_(spm, key) (*(map_vpd_insertptr(spm, key)))

/* true if key exists */
#define map_vpd_haskey(spm, key) ((map_vpd_ptr(spm, key)) != NULL)

/* Prototypes */

/* initialize the structure */
	extern mapkit_error map_vpd_init(map_vpd * spm);

/* free contents */
	extern void map_vpd_free(map_vpd * spm);

/* copy the "from" map into a new "to" map */
	extern mapkit_error map_vpd_copy(map_vpd * to, map_vpd * from);

/* initialize the structure for "used" elements */
	extern mapkit_error map_vpd_init_hint(map_vpd * spm, mapkit_size_t used);

/* ensure at least "used" slots (including currently occupied slots) */
	extern mapkit_error map_vpd_ensurecapacity(map_vpd * spm, mapkit_size_t used);

/* restore maxfill & minused and ensure these constraints */
	extern mapkit_error map_vpd_adjustcapacity(map_vpd * spm);

/* reallocate contents for "size" slots, keeping previous values */
	extern mapkit_error map_vpd_reallocate(map_vpd * spm, mapkit_size_t size);

/* return an adequate size for growing from "used" elements */
	extern mapkit_size_t map_vpd_growsize(map_vpd * spm, mapkit_size_t used);

/* return an adequate size for shrinking from "used" elements */
	extern mapkit_size_t map_vpd_shrinksize(map_vpd * spm, mapkit_size_t used);

/* return an adequate size for "used" elements */
	extern mapkit_size_t map_vpd_meansize(map_vpd * spm, mapkit_size_t used);

/* pointer to the value at key (insert as needed) */
	INLINE double *map_vpd_insertptr(map_vpd * spm, void *key);

/* pointer to the value at key or NULL */
	INLINE double *map_vpd_ptr(map_vpd * spm, void *key);

/* remove the key pointing to the value *ptr. No checking is done */
	INLINE mapkit_error map_vpd_removeptr(map_vpd * spm, double *ptr);

/* set the value at key (insert as needed) */
	INLINE mapkit_error map_vpd_set(map_vpd * spm, void *key, double value);

/* return the value at key (key must exists if alwaysdefaultvalue is not set) */
	INLINE double map_vpd_value(map_vpd * spm, void *key);

/* get the value at key (return an error code) */
	INLINE mapkit_error map_vpd_get(map_vpd * spm, void *key, double *value);

/* remove key (key must exists) */
	INLINE mapkit_error map_vpd_remove(map_vpd * spm, void *key);

/*
 * return the next index to a full slot of the map, or -1 if none
 * if index == -1, returns the first index, or -1 if none
 * the contents of the slot can be accessed by spm->contents[index].key and
 * spm->contents[index].value
 */
	extern mapkit_size_t map_vpd_next(map_vpd * spm, mapkit_size_t iindex);

/*
 * return the next ptr to a full slot of the map, or NULL if none
 * if pos_contents == NULL, returns the first full slot, or NULL if none
 * the contents of the slot can be accessed by pos_contents->key and
 * pos_contents->value
 */
	extern map_vpd_storage *map_vpd_nextptr(map_vpd * spm, map_vpd_storage * pos_contents);

/* allocate an array of map_vpd_element with all "count" (key,value) pairs */
	extern mapkit_error map_vpd_getall(map_vpd * spm, map_vpd_element ** array, mapkit_size_t * count);

/* remove all values equal to the defaultvalue */
	extern mapkit_error map_vpd_clean(map_vpd * spm);

/* compare the key of two map_vpd_element (for qsort) */
	extern int map_vpd_compare(const void *e1, const void *e2);

/* allocate an array of map_vpd_element with all "count" (key,value) pairs */
	extern mapkit_error map_vpd_getall_sorted(map_vpd * spm, map_vpd_element ** array, mapkit_size_t * count);

/* insert all count map_vpd_element into spm */
	extern mapkit_error map_vpd_setall(map_vpd * spm, map_vpd_element * array, mapkit_size_t count);

/* remove all "count" elements whose keys are in "array" */
	extern mapkit_error map_vpd_removeall(map_vpd * spm, void **array, mapkit_size_t count);

/* print statistics */
	extern void map_vpd_printstats(map_vpd * spm);

/* Private functions */

	extern mapkit_error map_vpd_remove_s(map_vpd * spm, void *key);
	extern double map_vpd_value_s(map_vpd * spm, void *key);
	extern mapkit_error map_vpd_get_s(map_vpd * spm, void *key, double *value);
	extern mapkit_error map_vpd_set_s(map_vpd * spm, void *key, double value);
	extern double *map_vpd_insertptr_s(map_vpd * spm, void *key);
	extern double *map_vpd_ptr_s(map_vpd * spm, void *key);

/* Implementation */

/* INLINEd functions */
	mapkit_error map_vpd_remove(map_vpd * spm, void *key) {
		map_vpd_storage *contents;

		contents = &(spm->contents[((mapkit_hash_t) key) % spm->size]);

		if ((contents->state == MAPKIT_FULLSLOT)
		    && ((contents->key) == (key))) {
			contents->state = MAPKIT_DELETEDSLOT;
#ifdef MAPKIT_COLLISIONS
			spm->keyindexs++;
#endif
			spm->used--;
			if (spm->used < spm->minused) {
#ifdef MAPKIT_DEBUG
				fprintf(stderr, "MAPKIT: used < minused\n");
#endif
				return map_vpd_reallocate(spm, map_vpd_shrinksize(spm, spm->used));
			}
		} else
			return map_vpd_remove_s(spm, key);
		return MAPKIT_OK;
	}

	mapkit_error map_vpd_removeptr(map_vpd * spm, double *ptr) {
		map_vpd_storage *sptr = (map_vpd_storage *) ((char *) ptr - offsetof(map_vpd_storage, value));

		sptr->state = MAPKIT_DELETEDSLOT;
		spm->used--;
		if (spm->used < spm->minused) {
#ifdef MAPKIT_DEBUG
			fprintf(stderr, "MAPKIT: used < minused\n");
#endif
			return map_vpd_reallocate(spm, map_vpd_shrinksize(spm, spm->used));
		}

		return MAPKIT_OK;
	}

	mapkit_error map_vpd_set(map_vpd * spm, void *key, double value) {
		map_vpd_storage *contents;

		if (spm->alwaysdefault && ((value) == (spm->defaultvalue)))
			return map_vpd_remove(spm, key);

		contents = &(spm->contents[((mapkit_hash_t) key) % spm->size]);

		if (contents->state == MAPKIT_FREESLOT) {
			contents->state = MAPKIT_FULLSLOT;
			contents->key = key;
			contents->value = value;
#ifdef MAPKIT_COLLISIONS
			spm->insertionindexs++;
#endif
			spm->used++;
			spm->fill++;
			if (spm->fill > spm->maxfill) {
#ifdef MAPKIT_DEBUG
				fprintf(stderr, "MAPKIT: fill > maxfill\n");
#endif
				return map_vpd_reallocate(spm, map_vpd_growsize(spm, spm->used));
			}
			return MAPKIT_OK;
		} else if ((contents->state == MAPKIT_FULLSLOT)
			   && ((contents->key) == (key))) {
			contents->value = value;
#ifdef MAPKIT_COLLISIONS
			spm->insertionindexs++;
#endif
			return MAPKIT_OK;
		} else
			return map_vpd_set_s(spm, key, value);
	}

	double map_vpd_value(map_vpd * spm, void *key) {
		map_vpd_storage *contents;

		contents = &(spm->contents[((mapkit_hash_t) key) % spm->size]);

		if ((contents->state == MAPKIT_FULLSLOT)
		    && ((contents->key) == (key))) {
#ifdef MAPKIT_COLLISIONS
			spm->keyindexs++;
#endif
			return contents->value;
		} else if (spm->alwaysdefault && (contents->state == MAPKIT_FREESLOT)) {
#ifdef MAPKIT_COLLISIONS
			spm->keyindexs++;
#endif
			return spm->defaultvalue;
		} else
			return map_vpd_value_s(spm, key);
	}

	mapkit_error map_vpd_get(map_vpd * spm, void *key, double *value) {
		map_vpd_storage *contents;

		contents = &(spm->contents[((mapkit_hash_t) key) % spm->size]);

		if ((contents->state == MAPKIT_FULLSLOT)
		    && ((contents->key) == (key))) {
#ifdef MAPKIT_COLLISIONS
			spm->keyindexs++;
#endif
			*value = contents->value;
		} else if (spm->alwaysdefault && (contents->state == MAPKIT_FREESLOT)) {
#ifdef MAPKIT_COLLISIONS
			spm->keyindexs++;
#endif
			*value = spm->defaultvalue;
		} else
			return map_vpd_get_s(spm, key, value);

		return MAPKIT_OK;
	}

	double *map_vpd_insertptr(map_vpd * spm, void *key) {
		map_vpd_storage *contents;

		contents = &(spm->contents[((mapkit_hash_t) key) % spm->size]);

		if (contents->state == MAPKIT_FREESLOT) {
			if (spm->fill >= spm->maxfill)
				return map_vpd_insertptr_s(spm, key);

			contents->state = MAPKIT_FULLSLOT;
			contents->key = key;
			contents->value = spm->defaultvalue;
#ifdef MAPKIT_COLLISIONS
			spm->insertionindexs++;
#endif
			spm->used++;
			spm->fill++;
			return &(contents->value);
		} else if ((contents->state == MAPKIT_FULLSLOT)
			   && ((contents->key) == (key))) {
#ifdef MAPKIT_COLLISIONS
			spm->insertionindexs++;
#endif
			return &(contents->value);
		} else
			return map_vpd_insertptr_s(spm, key);
	}

	double *map_vpd_ptr(map_vpd * spm, void *key) {
		map_vpd_storage *contents;

		contents = &(spm->contents[((mapkit_hash_t) key) % spm->size]);

		if ((contents->state == MAPKIT_FULLSLOT)
		    && ((contents->key) == (key))) {
#ifdef MAPKIT_COLLISIONS
			spm->insertionindexs++;
#endif
			if ((!spm->alwaysdefault) || (!((contents->value) == (spm->defaultvalue))))
				return &(contents->value);
			else
				return NULL;
		} else if (contents->state == MAPKIT_FREESLOT) {
#ifdef MAPKIT_COLLISIONS
			spm->insertionindexs++;
#endif
			return NULL;
		} else
			return map_vpd_ptr_s(spm, key);
	}

/*
  Prototypes for map_vpvp (void* -> void*)
  Default value : NULL
  Uses a state field.
*/

	typedef struct {
		signed char state;
		void *key;
		void *value;
	} map_vpvp_storage;

	typedef struct {
		void *key;
		void *value;
	} map_vpvp_element;

	typedef struct {
		mapkit_size_t size,			       /* size of the hash table */
		 fill,					       /* # of non-free slots in the hash table */
		 used,					       /* # of full slots */
		 maxfill,				       /* max # of non-free slots before rehash */
		 minused /* min # of full slots before rehash */ ;

		double maxfillfactor,			       /* maxfill = maxfillfactor * size */
		 minusedfactor;				       /* minused = minusedfactor * size */

#ifdef MAPKIT_COLLISIONS
		unsigned long insertionindexs, insertionindex_collisions;
		unsigned long keyindexs, keyindex_collisions;
#endif

		void *defaultvalue;

		int alwaysdefault;

		map_vpvp_storage *contents;
	} map_vpvp;

/* Macros */

/* set to inform the type is available */
#define MAPKIT_map_vpvp

/* read/write access to the value at key (insert as needed) */
#define map_vpvp_(spm, key) (*(map_vpvp_insertptr(spm, key)))

/* true if key exists */
#define map_vpvp_haskey(spm, key) ((map_vpvp_ptr(spm, key)) != NULL)

/* Prototypes */

/* initialize the structure */
	extern mapkit_error map_vpvp_init(map_vpvp * spm);

/* free contents */
	extern void map_vpvp_free(map_vpvp * spm);

/* copy the "from" map into a new "to" map */
	extern mapkit_error map_vpvp_copy(map_vpvp * to, map_vpvp * from);

/* initialize the structure for "used" elements */
	extern mapkit_error map_vpvp_init_hint(map_vpvp * spm, mapkit_size_t used);

/* ensure at least "used" slots (including currently occupied slots) */
	extern mapkit_error map_vpvp_ensurecapacity(map_vpvp * spm, mapkit_size_t used);

/* restore maxfill & minused and ensure these constraints */
	extern mapkit_error map_vpvp_adjustcapacity(map_vpvp * spm);

/* reallocate contents for "size" slots, keeping previous values */
	extern mapkit_error map_vpvp_reallocate(map_vpvp * spm, mapkit_size_t size);

/* return an adequate size for growing from "used" elements */
	extern mapkit_size_t map_vpvp_growsize(map_vpvp * spm, mapkit_size_t used);

/* return an adequate size for shrinking from "used" elements */
	extern mapkit_size_t map_vpvp_shrinksize(map_vpvp * spm, mapkit_size_t used);

/* return an adequate size for "used" elements */
	extern mapkit_size_t map_vpvp_meansize(map_vpvp * spm, mapkit_size_t used);

/* pointer to the value at key (insert as needed) */
	INLINE void **map_vpvp_insertptr(map_vpvp * spm, void *key);

/* pointer to the value at key or NULL */
	INLINE void **map_vpvp_ptr(map_vpvp * spm, void *key);

/* remove the key pointing to the value *ptr. No checking is done */
	INLINE mapkit_error map_vpvp_removeptr(map_vpvp * spm, void **ptr);

/* set the value at key (insert as needed) */
	INLINE mapkit_error map_vpvp_set(map_vpvp * spm, void *key, void *value);

/* return the value at key (key must exists if alwaysdefaultvalue is not set) */
	INLINE void *map_vpvp_value(map_vpvp * spm, void *key);

/* get the value at key (return an error code) */
	INLINE mapkit_error map_vpvp_get(map_vpvp * spm, void *key, void **value);

/* remove key (key must exists) */
	INLINE mapkit_error map_vpvp_remove(map_vpvp * spm, void *key);

/*
 * return the next index to a full slot of the map, or -1 if none
 * if index == -1, returns the first index, or -1 if none
 * the contents of the slot can be accessed by spm->contents[index].key and
 * spm->contents[index].value
 */
	extern mapkit_size_t map_vpvp_next(map_vpvp * spm, mapkit_size_t iindex);

/*
 * return the next ptr to a full slot of the map, or NULL if none
 * if pos_contents == NULL, returns the first full slot, or NULL if none
 * the contents of the slot can be accessed by pos_contents->key and
 * pos_contents->value
 */
	extern map_vpvp_storage *map_vpvp_nextptr(map_vpvp * spm, map_vpvp_storage * pos_contents);

/* allocate an array of map_vpvp_element with all "count" (key,value) pairs */
	extern mapkit_error map_vpvp_getall(map_vpvp * spm, map_vpvp_element ** array, mapkit_size_t * count);

/* remove all values equal to the defaultvalue */
	extern mapkit_error map_vpvp_clean(map_vpvp * spm);

/* compare the key of two map_vpvp_element (for qsort) */
	extern int map_vpvp_compare(const void *e1, const void *e2);

/* allocate an array of map_vpvp_element with all "count" (key,value) pairs */
	extern mapkit_error map_vpvp_getall_sorted(map_vpvp * spm, map_vpvp_element ** array, mapkit_size_t * count);

/* insert all count map_vpvp_element into spm */
	extern mapkit_error map_vpvp_setall(map_vpvp * spm, map_vpvp_element * array, mapkit_size_t count);

/* remove all "count" elements whose keys are in "array" */
	extern mapkit_error map_vpvp_removeall(map_vpvp * spm, void **array, mapkit_size_t count);

/* print statistics */
	extern void map_vpvp_printstats(map_vpvp * spm);

/* Private functions */

	extern mapkit_error map_vpvp_remove_s(map_vpvp * spm, void *key);
	extern void *map_vpvp_value_s(map_vpvp * spm, void *key);
	extern mapkit_error map_vpvp_get_s(map_vpvp * spm, void *key, void **value);
	extern mapkit_error map_vpvp_set_s(map_vpvp * spm, void *key, void *value);
	extern void **map_vpvp_insertptr_s(map_vpvp * spm, void *key);
	extern void **map_vpvp_ptr_s(map_vpvp * spm, void *key);

/* Implementation */

/* INLINEd functions */
	mapkit_error map_vpvp_remove(map_vpvp * spm, void *key) {
		map_vpvp_storage *contents;

		contents = &(spm->contents[((mapkit_hash_t) key) % spm->size]);

		if ((contents->state == MAPKIT_FULLSLOT)
		    && ((contents->key) == (key))) {
			contents->state = MAPKIT_DELETEDSLOT;
#ifdef MAPKIT_COLLISIONS
			spm->keyindexs++;
#endif
			spm->used--;
			if (spm->used < spm->minused) {
#ifdef MAPKIT_DEBUG
				fprintf(stderr, "MAPKIT: used < minused\n");
#endif
				return map_vpvp_reallocate(spm, map_vpvp_shrinksize(spm, spm->used));
			}
		} else
			return map_vpvp_remove_s(spm, key);
		return MAPKIT_OK;
	}

	mapkit_error map_vpvp_removeptr(map_vpvp * spm, void **ptr) {
		map_vpvp_storage *sptr = (map_vpvp_storage *) ((char *) ptr - offsetof(map_vpvp_storage, value));

		sptr->state = MAPKIT_DELETEDSLOT;
		spm->used--;
		if (spm->used < spm->minused) {
#ifdef MAPKIT_DEBUG
			fprintf(stderr, "MAPKIT: used < minused\n");
#endif
			return map_vpvp_reallocate(spm, map_vpvp_shrinksize(spm, spm->used));
		}

		return MAPKIT_OK;
	}

	mapkit_error map_vpvp_set(map_vpvp * spm, void *key, void *value) {
		map_vpvp_storage *contents;

		if (spm->alwaysdefault && ((value) == (spm->defaultvalue)))
			return map_vpvp_remove(spm, key);

		contents = &(spm->contents[((mapkit_hash_t) key) % spm->size]);

		if (contents->state == MAPKIT_FREESLOT) {
			contents->state = MAPKIT_FULLSLOT;
			contents->key = key;
			contents->value = value;
#ifdef MAPKIT_COLLISIONS
			spm->insertionindexs++;
#endif
			spm->used++;
			spm->fill++;
			if (spm->fill > spm->maxfill) {
#ifdef MAPKIT_DEBUG
				fprintf(stderr, "MAPKIT: fill > maxfill\n");
#endif
				return map_vpvp_reallocate(spm, map_vpvp_growsize(spm, spm->used));
			}
			return MAPKIT_OK;
		} else if ((contents->state == MAPKIT_FULLSLOT)
			   && ((contents->key) == (key))) {
			contents->value = value;
#ifdef MAPKIT_COLLISIONS
			spm->insertionindexs++;
#endif
			return MAPKIT_OK;
		} else
			return map_vpvp_set_s(spm, key, value);
	}

	void *map_vpvp_value(map_vpvp * spm, void *key) {
		map_vpvp_storage *contents;

		contents = &(spm->contents[((mapkit_hash_t) key) % spm->size]);

		if ((contents->state == MAPKIT_FULLSLOT)
		    && ((contents->key) == (key))) {
#ifdef MAPKIT_COLLISIONS
			spm->keyindexs++;
#endif
			return contents->value;
		} else if (spm->alwaysdefault && (contents->state == MAPKIT_FREESLOT)) {
#ifdef MAPKIT_COLLISIONS
			spm->keyindexs++;
#endif
			return spm->defaultvalue;
		} else
			return map_vpvp_value_s(spm, key);
	}

	mapkit_error map_vpvp_get(map_vpvp * spm, void *key, void **value) {
		map_vpvp_storage *contents;

		contents = &(spm->contents[((mapkit_hash_t) key) % spm->size]);

		if ((contents->state == MAPKIT_FULLSLOT)
		    && ((contents->key) == (key))) {
#ifdef MAPKIT_COLLISIONS
			spm->keyindexs++;
#endif
			*value = contents->value;
		} else if (spm->alwaysdefault && (contents->state == MAPKIT_FREESLOT)) {
#ifdef MAPKIT_COLLISIONS
			spm->keyindexs++;
#endif
			*value = spm->defaultvalue;
		} else
			return map_vpvp_get_s(spm, key, value);

		return MAPKIT_OK;
	}

	void **map_vpvp_insertptr(map_vpvp * spm, void *key) {
		map_vpvp_storage *contents;

		contents = &(spm->contents[((mapkit_hash_t) key) % spm->size]);

		if (contents->state == MAPKIT_FREESLOT) {
			if (spm->fill >= spm->maxfill)
				return map_vpvp_insertptr_s(spm, key);

			contents->state = MAPKIT_FULLSLOT;
			contents->key = key;
			contents->value = spm->defaultvalue;
#ifdef MAPKIT_COLLISIONS
			spm->insertionindexs++;
#endif
			spm->used++;
			spm->fill++;
			return &(contents->value);
		} else if ((contents->state == MAPKIT_FULLSLOT)
			   && ((contents->key) == (key))) {
#ifdef MAPKIT_COLLISIONS
			spm->insertionindexs++;
#endif
			return &(contents->value);
		} else
			return map_vpvp_insertptr_s(spm, key);
	}

	void **map_vpvp_ptr(map_vpvp * spm, void *key) {
		map_vpvp_storage *contents;

		contents = &(spm->contents[((mapkit_hash_t) key) % spm->size]);

		if ((contents->state == MAPKIT_FULLSLOT)
		    && ((contents->key) == (key))) {
#ifdef MAPKIT_COLLISIONS
			spm->insertionindexs++;
#endif
			if ((!spm->alwaysdefault) || (!((contents->value) == (spm->defaultvalue))))
				return &(contents->value);
			else
				return NULL;
		} else if (contents->state == MAPKIT_FREESLOT) {
#ifdef MAPKIT_COLLISIONS
			spm->insertionindexs++;
#endif
			return NULL;
		} else
			return map_vpvp_ptr_s(spm, key);
	}

/*
  Prototypes for map_h_vpi (void* -> int)
  Default value : 0
  Uses a state field.
*/

	typedef struct {
		signed char state;
		void *key;
		int value;
	} map_h_vpi_storage;

	typedef struct {
		void *key;
		int value;
	} map_h_vpi_element;

	typedef struct {
		mapkit_size_t size,			       /* size of the hash table */
		 fill,					       /* # of non-free slots in the hash table */
		 used,					       /* # of full slots */
		 maxfill,				       /* max # of non-free slots before rehash */
		 minused /* min # of full slots before rehash */ ;

		double maxfillfactor,			       /* maxfill = maxfillfactor * size */
		 minusedfactor;				       /* minused = minusedfactor * size */

#ifdef MAPKIT_COLLISIONS
		unsigned long insertionindexs, insertionindex_collisions;
		unsigned long keyindexs, keyindex_collisions;
#endif

		int defaultvalue;

		int alwaysdefault;

		map_h_vpi_storage *contents;
	} map_h_vpi;

/* Macros */

/* set to inform the type is available */
#define MAPKIT_map_h_vpi

/* read/write access to the value at key (insert as needed) */
#define map_h_vpi_(spm, key) (*(map_h_vpi_insertptr(spm, key)))

/* true if key exists */
#define map_h_vpi_haskey(spm, key) ((map_h_vpi_ptr(spm, key)) != NULL)

/* Prototypes */

/* initialize the structure */
	extern mapkit_error map_h_vpi_init(map_h_vpi * spm);

/* free contents */
	extern void map_h_vpi_free(map_h_vpi * spm);

/* copy the "from" map into a new "to" map */
	extern mapkit_error map_h_vpi_copy(map_h_vpi * to, map_h_vpi * from);

/* initialize the structure for "used" elements */
	extern mapkit_error map_h_vpi_init_hint(map_h_vpi * spm, mapkit_size_t used);

/* ensure at least "used" slots (including currently occupied slots) */
	extern mapkit_error map_h_vpi_ensurecapacity(map_h_vpi * spm, mapkit_size_t used);

/* restore maxfill & minused and ensure these constraints */
	extern mapkit_error map_h_vpi_adjustcapacity(map_h_vpi * spm);

/* reallocate contents for "size" slots, keeping previous values */
	extern mapkit_error map_h_vpi_reallocate(map_h_vpi * spm, mapkit_size_t size);

/* return an adequate size for growing from "used" elements */
	extern mapkit_size_t map_h_vpi_growsize(map_h_vpi * spm, mapkit_size_t used);

/* return an adequate size for shrinking from "used" elements */
	extern mapkit_size_t map_h_vpi_shrinksize(map_h_vpi * spm, mapkit_size_t used);

/* return an adequate size for "used" elements */
	extern mapkit_size_t map_h_vpi_meansize(map_h_vpi * spm, mapkit_size_t used);

/* pointer to the value at key (insert as needed) */
	INLINE int *map_h_vpi_insertptr(map_h_vpi * spm, void *key);

/* pointer to the value at key or NULL */
	INLINE int *map_h_vpi_ptr(map_h_vpi * spm, void *key);

/* remove the key pointing to the value *ptr. No checking is done */
	INLINE mapkit_error map_h_vpi_removeptr(map_h_vpi * spm, int *ptr);

/* set the value at key (insert as needed) */
	INLINE mapkit_error map_h_vpi_set(map_h_vpi * spm, void *key, int value);

/* return the value at key (key must exists if alwaysdefaultvalue is not set) */
	INLINE int map_h_vpi_value(map_h_vpi * spm, void *key);

/* get the value at key (return an error code) */
	INLINE mapkit_error map_h_vpi_get(map_h_vpi * spm, void *key, int *value);

/* remove key (key must exists) */
	INLINE mapkit_error map_h_vpi_remove(map_h_vpi * spm, void *key);

/*
 * return the next index to a full slot of the map, or -1 if none
 * if index == -1, returns the first index, or -1 if none
 * the contents of the slot can be accessed by spm->contents[index].key and
 * spm->contents[index].value
 */
	extern mapkit_size_t map_h_vpi_next(map_h_vpi * spm, mapkit_size_t iindex);

/*
 * return the next ptr to a full slot of the map, or NULL if none
 * if pos_contents == NULL, returns the first full slot, or NULL if none
 * the contents of the slot can be accessed by pos_contents->key and
 * pos_contents->value
 */
	extern map_h_vpi_storage *map_h_vpi_nextptr(map_h_vpi * spm, map_h_vpi_storage * pos_contents);

/* allocate an array of map_h_vpi_element with all "count" (key,value) pairs */
	extern mapkit_error map_h_vpi_getall(map_h_vpi * spm, map_h_vpi_element ** array, mapkit_size_t * count);

/* remove all values equal to the defaultvalue */
	extern mapkit_error map_h_vpi_clean(map_h_vpi * spm);

/* compare the key of two map_h_vpi_element (for qsort) */
	extern int map_h_vpi_compare(const void *e1, const void *e2);

/* allocate an array of map_h_vpi_element with all "count" (key,value) pairs */
	extern mapkit_error map_h_vpi_getall_sorted(map_h_vpi * spm, map_h_vpi_element ** array, mapkit_size_t * count);

/* insert all count map_h_vpi_element into spm */
	extern mapkit_error map_h_vpi_setall(map_h_vpi * spm, map_h_vpi_element * array, mapkit_size_t count);

/* remove all "count" elements whose keys are in "array" */
	extern mapkit_error map_h_vpi_removeall(map_h_vpi * spm, void **array, mapkit_size_t count);

/* print statistics */
	extern void map_h_vpi_printstats(map_h_vpi * spm);

/* Private functions */

	extern mapkit_error map_h_vpi_remove_s(map_h_vpi * spm, void *key, mapkit_hash_t hash);
	extern int map_h_vpi_value_s(map_h_vpi * spm, void *key, mapkit_hash_t hash);
	extern mapkit_error map_h_vpi_get_s(map_h_vpi * spm, void *key, int *value, mapkit_hash_t hash);
	extern mapkit_error map_h_vpi_set_s(map_h_vpi * spm, void *key, int value, mapkit_hash_t hash);
	extern int *map_h_vpi_insertptr_s(map_h_vpi * spm, void *key, mapkit_hash_t hash);
	extern int *map_h_vpi_ptr_s(map_h_vpi * spm, void *key, mapkit_hash_t hash);

/* Implementation */

/* INLINEd functions */
	mapkit_error map_h_vpi_remove(map_h_vpi * spm, void *key) {
		mapkit_hash_t hash;
		map_h_vpi_storage *contents;

		contents = &(spm->contents[(hash = mapkit_hash((mapkit_hash_t) key)) % spm->size]);

		if ((contents->state == MAPKIT_FULLSLOT)
		    && ((contents->key) == (key))) {
			contents->state = MAPKIT_DELETEDSLOT;
#ifdef MAPKIT_COLLISIONS
			spm->keyindexs++;
#endif
			spm->used--;
			if (spm->used < spm->minused) {
#ifdef MAPKIT_DEBUG
				fprintf(stderr, "MAPKIT: used < minused\n");
#endif
				return map_h_vpi_reallocate(spm, map_h_vpi_shrinksize(spm, spm->used));
			}
		} else
			return map_h_vpi_remove_s(spm, key, hash);
		return MAPKIT_OK;
	}

	mapkit_error map_h_vpi_removeptr(map_h_vpi * spm, int *ptr) {
		map_h_vpi_storage *sptr = (map_h_vpi_storage *) ((char *) ptr - offsetof(map_h_vpi_storage, value));

		sptr->state = MAPKIT_DELETEDSLOT;
		spm->used--;
		if (spm->used < spm->minused) {
#ifdef MAPKIT_DEBUG
			fprintf(stderr, "MAPKIT: used < minused\n");
#endif
			return map_h_vpi_reallocate(spm, map_h_vpi_shrinksize(spm, spm->used));
		}

		return MAPKIT_OK;
	}

	mapkit_error map_h_vpi_set(map_h_vpi * spm, void *key, int value) {
		mapkit_hash_t hash;
		map_h_vpi_storage *contents;

		if (spm->alwaysdefault && ((value) == (spm->defaultvalue)))
			return map_h_vpi_remove(spm, key);

		contents = &(spm->contents[(hash = mapkit_hash((mapkit_hash_t) key)) % spm->size]);

		if (contents->state == MAPKIT_FREESLOT) {
			contents->state = MAPKIT_FULLSLOT;
			contents->key = key;
			contents->value = value;
#ifdef MAPKIT_COLLISIONS
			spm->insertionindexs++;
#endif
			spm->used++;
			spm->fill++;
			if (spm->fill > spm->maxfill) {
#ifdef MAPKIT_DEBUG
				fprintf(stderr, "MAPKIT: fill > maxfill\n");
#endif
				return map_h_vpi_reallocate(spm, map_h_vpi_growsize(spm, spm->used));
			}
			return MAPKIT_OK;
		} else if ((contents->state == MAPKIT_FULLSLOT)
			   && ((contents->key) == (key))) {
			contents->value = value;
#ifdef MAPKIT_COLLISIONS
			spm->insertionindexs++;
#endif
			return MAPKIT_OK;
		} else
			return map_h_vpi_set_s(spm, key, value, hash);
	}

	int map_h_vpi_value(map_h_vpi * spm, void *key) {
		mapkit_hash_t hash;
		map_h_vpi_storage *contents;

		contents = &(spm->contents[(hash = mapkit_hash((mapkit_hash_t) key)) % spm->size]);

		if ((contents->state == MAPKIT_FULLSLOT)
		    && ((contents->key) == (key))) {
#ifdef MAPKIT_COLLISIONS
			spm->keyindexs++;
#endif
			return contents->value;
		} else if (spm->alwaysdefault && (contents->state == MAPKIT_FREESLOT)) {
#ifdef MAPKIT_COLLISIONS
			spm->keyindexs++;
#endif
			return spm->defaultvalue;
		} else
			return map_h_vpi_value_s(spm, key, hash);
	}

	mapkit_error map_h_vpi_get(map_h_vpi * spm, void *key, int *value) {
		mapkit_hash_t hash;
		map_h_vpi_storage *contents;

		contents = &(spm->contents[(hash = mapkit_hash((mapkit_hash_t) key)) % spm->size]);

		if ((contents->state == MAPKIT_FULLSLOT)
		    && ((contents->key) == (key))) {
#ifdef MAPKIT_COLLISIONS
			spm->keyindexs++;
#endif
			*value = contents->value;
		} else if (spm->alwaysdefault && (contents->state == MAPKIT_FREESLOT)) {
#ifdef MAPKIT_COLLISIONS
			spm->keyindexs++;
#endif
			*value = spm->defaultvalue;
		} else
			return map_h_vpi_get_s(spm, key, value, hash);

		return MAPKIT_OK;
	}

	int *map_h_vpi_insertptr(map_h_vpi * spm, void *key) {
		mapkit_hash_t hash;
		map_h_vpi_storage *contents;

		contents = &(spm->contents[(hash = mapkit_hash((mapkit_hash_t) key)) % spm->size]);

		if (contents->state == MAPKIT_FREESLOT) {
			if (spm->fill >= spm->maxfill)
				return map_h_vpi_insertptr_s(spm, key, hash);

			contents->state = MAPKIT_FULLSLOT;
			contents->key = key;
			contents->value = spm->defaultvalue;
#ifdef MAPKIT_COLLISIONS
			spm->insertionindexs++;
#endif
			spm->used++;
			spm->fill++;
			return &(contents->value);
		} else if ((contents->state == MAPKIT_FULLSLOT)
			   && ((contents->key) == (key))) {
#ifdef MAPKIT_COLLISIONS
			spm->insertionindexs++;
#endif
			return &(contents->value);
		} else
			return map_h_vpi_insertptr_s(spm, key, hash);
	}

	int *map_h_vpi_ptr(map_h_vpi * spm, void *key) {
		mapkit_hash_t hash;
		map_h_vpi_storage *contents;

		contents = &(spm->contents[(hash = mapkit_hash((mapkit_hash_t) key)) % spm->size]);

		if ((contents->state == MAPKIT_FULLSLOT)
		    && ((contents->key) == (key))) {
#ifdef MAPKIT_COLLISIONS
			spm->insertionindexs++;
#endif
			if ((!spm->alwaysdefault) || (!((contents->value) == (spm->defaultvalue))))
				return &(contents->value);
			else
				return NULL;
		} else if (contents->state == MAPKIT_FREESLOT) {
#ifdef MAPKIT_COLLISIONS
			spm->insertionindexs++;
#endif
			return NULL;
		} else
			return map_h_vpi_ptr_s(spm, key, hash);
	}

/*
  Prototypes for map_h_vpd (void* -> double)
  Default value : 0.0
  Uses a state field.
*/

	typedef struct {
		signed char state;
		void *key;
		double value;
	} map_h_vpd_storage;

	typedef struct {
		void *key;
		double value;
	} map_h_vpd_element;

	typedef struct {
		mapkit_size_t size,			       /* size of the hash table */
		 fill,					       /* # of non-free slots in the hash table */
		 used,					       /* # of full slots */
		 maxfill,				       /* max # of non-free slots before rehash */
		 minused /* min # of full slots before rehash */ ;

		double maxfillfactor,			       /* maxfill = maxfillfactor * size */
		 minusedfactor;				       /* minused = minusedfactor * size */

#ifdef MAPKIT_COLLISIONS
		unsigned long insertionindexs, insertionindex_collisions;
		unsigned long keyindexs, keyindex_collisions;
#endif

		double defaultvalue;

		int alwaysdefault;

		map_h_vpd_storage *contents;
	} map_h_vpd;

/* Macros */

/* set to inform the type is available */
#define MAPKIT_map_h_vpd

/* read/write access to the value at key (insert as needed) */
#define map_h_vpd_(spm, key) (*(map_h_vpd_insertptr(spm, key)))

/* true if key exists */
#define map_h_vpd_haskey(spm, key) ((map_h_vpd_ptr(spm, key)) != NULL)

/* Prototypes */

/* initialize the structure */
	extern mapkit_error map_h_vpd_init(map_h_vpd * spm);

/* free contents */
	extern void map_h_vpd_free(map_h_vpd * spm);

/* copy the "from" map into a new "to" map */
	extern mapkit_error map_h_vpd_copy(map_h_vpd * to, map_h_vpd * from);

/* initialize the structure for "used" elements */
	extern mapkit_error map_h_vpd_init_hint(map_h_vpd * spm, mapkit_size_t used);

/* ensure at least "used" slots (including currently occupied slots) */
	extern mapkit_error map_h_vpd_ensurecapacity(map_h_vpd * spm, mapkit_size_t used);

/* restore maxfill & minused and ensure these constraints */
	extern mapkit_error map_h_vpd_adjustcapacity(map_h_vpd * spm);

/* reallocate contents for "size" slots, keeping previous values */
	extern mapkit_error map_h_vpd_reallocate(map_h_vpd * spm, mapkit_size_t size);

/* return an adequate size for growing from "used" elements */
	extern mapkit_size_t map_h_vpd_growsize(map_h_vpd * spm, mapkit_size_t used);

/* return an adequate size for shrinking from "used" elements */
	extern mapkit_size_t map_h_vpd_shrinksize(map_h_vpd * spm, mapkit_size_t used);

/* return an adequate size for "used" elements */
	extern mapkit_size_t map_h_vpd_meansize(map_h_vpd * spm, mapkit_size_t used);

/* pointer to the value at key (insert as needed) */
	INLINE double *map_h_vpd_insertptr(map_h_vpd * spm, void *key);

/* pointer to the value at key or NULL */
	INLINE double *map_h_vpd_ptr(map_h_vpd * spm, void *key);

/* remove the key pointing to the value *ptr. No checking is done */
	INLINE mapkit_error map_h_vpd_removeptr(map_h_vpd * spm, double *ptr);

/* set the value at key (insert as needed) */
	INLINE mapkit_error map_h_vpd_set(map_h_vpd * spm, void *key, double value);

/* return the value at key (key must exists if alwaysdefaultvalue is not set) */
	INLINE double map_h_vpd_value(map_h_vpd * spm, void *key);

/* get the value at key (return an error code) */
	INLINE mapkit_error map_h_vpd_get(map_h_vpd * spm, void *key, double *value);

/* remove key (key must exists) */
	INLINE mapkit_error map_h_vpd_remove(map_h_vpd * spm, void *key);

/*
 * return the next index to a full slot of the map, or -1 if none
 * if index == -1, returns the first index, or -1 if none
 * the contents of the slot can be accessed by spm->contents[index].key and
 * spm->contents[index].value
 */
	extern mapkit_size_t map_h_vpd_next(map_h_vpd * spm, mapkit_size_t iindex);

/*
 * return the next ptr to a full slot of the map, or NULL if none
 * if pos_contents == NULL, returns the first full slot, or NULL if none
 * the contents of the slot can be accessed by pos_contents->key and
 * pos_contents->value
 */
	extern map_h_vpd_storage *map_h_vpd_nextptr(map_h_vpd * spm, map_h_vpd_storage * pos_contents);

/* allocate an array of map_h_vpd_element with all "count" (key,value) pairs */
	extern mapkit_error map_h_vpd_getall(map_h_vpd * spm, map_h_vpd_element ** array, mapkit_size_t * count);

/* remove all values equal to the defaultvalue */
	extern mapkit_error map_h_vpd_clean(map_h_vpd * spm);

/* compare the key of two map_h_vpd_element (for qsort) */
	extern int map_h_vpd_compare(const void *e1, const void *e2);

/* allocate an array of map_h_vpd_element with all "count" (key,value) pairs */
	extern mapkit_error map_h_vpd_getall_sorted(map_h_vpd * spm, map_h_vpd_element ** array, mapkit_size_t * count);

/* insert all count map_h_vpd_element into spm */
	extern mapkit_error map_h_vpd_setall(map_h_vpd * spm, map_h_vpd_element * array, mapkit_size_t count);

/* remove all "count" elements whose keys are in "array" */
	extern mapkit_error map_h_vpd_removeall(map_h_vpd * spm, void **array, mapkit_size_t count);

/* print statistics */
	extern void map_h_vpd_printstats(map_h_vpd * spm);

/* Private functions */

	extern mapkit_error map_h_vpd_remove_s(map_h_vpd * spm, void *key, mapkit_hash_t hash);
	extern double map_h_vpd_value_s(map_h_vpd * spm, void *key, mapkit_hash_t hash);
	extern mapkit_error map_h_vpd_get_s(map_h_vpd * spm, void *key, double *value, mapkit_hash_t hash);
	extern mapkit_error map_h_vpd_set_s(map_h_vpd * spm, void *key, double value, mapkit_hash_t hash);
	extern double *map_h_vpd_insertptr_s(map_h_vpd * spm, void *key, mapkit_hash_t hash);
	extern double *map_h_vpd_ptr_s(map_h_vpd * spm, void *key, mapkit_hash_t hash);

/* Implementation */

/* INLINEd functions */
	mapkit_error map_h_vpd_remove(map_h_vpd * spm, void *key) {
		mapkit_hash_t hash;
		map_h_vpd_storage *contents;

		contents = &(spm->contents[(hash = mapkit_hash((mapkit_hash_t) key)) % spm->size]);

		if ((contents->state == MAPKIT_FULLSLOT)
		    && ((contents->key) == (key))) {
			contents->state = MAPKIT_DELETEDSLOT;
#ifdef MAPKIT_COLLISIONS
			spm->keyindexs++;
#endif
			spm->used--;
			if (spm->used < spm->minused) {
#ifdef MAPKIT_DEBUG
				fprintf(stderr, "MAPKIT: used < minused\n");
#endif
				return map_h_vpd_reallocate(spm, map_h_vpd_shrinksize(spm, spm->used));
			}
		} else
			return map_h_vpd_remove_s(spm, key, hash);
		return MAPKIT_OK;
	}

	mapkit_error map_h_vpd_removeptr(map_h_vpd * spm, double *ptr) {
		map_h_vpd_storage *sptr = (map_h_vpd_storage *) ((char *) ptr - offsetof(map_h_vpd_storage, value));

		sptr->state = MAPKIT_DELETEDSLOT;
		spm->used--;
		if (spm->used < spm->minused) {
#ifdef MAPKIT_DEBUG
			fprintf(stderr, "MAPKIT: used < minused\n");
#endif
			return map_h_vpd_reallocate(spm, map_h_vpd_shrinksize(spm, spm->used));
		}

		return MAPKIT_OK;
	}

	mapkit_error map_h_vpd_set(map_h_vpd * spm, void *key, double value) {
		mapkit_hash_t hash;
		map_h_vpd_storage *contents;

		if (spm->alwaysdefault && ((value) == (spm->defaultvalue)))
			return map_h_vpd_remove(spm, key);

		contents = &(spm->contents[(hash = mapkit_hash((mapkit_hash_t) key)) % spm->size]);

		if (contents->state == MAPKIT_FREESLOT) {
			contents->state = MAPKIT_FULLSLOT;
			contents->key = key;
			contents->value = value;
#ifdef MAPKIT_COLLISIONS
			spm->insertionindexs++;
#endif
			spm->used++;
			spm->fill++;
			if (spm->fill > spm->maxfill) {
#ifdef MAPKIT_DEBUG
				fprintf(stderr, "MAPKIT: fill > maxfill\n");
#endif
				return map_h_vpd_reallocate(spm, map_h_vpd_growsize(spm, spm->used));
			}
			return MAPKIT_OK;
		} else if ((contents->state == MAPKIT_FULLSLOT)
			   && ((contents->key) == (key))) {
			contents->value = value;
#ifdef MAPKIT_COLLISIONS
			spm->insertionindexs++;
#endif
			return MAPKIT_OK;
		} else
			return map_h_vpd_set_s(spm, key, value, hash);
	}

	double map_h_vpd_value(map_h_vpd * spm, void *key) {
		mapkit_hash_t hash;
		map_h_vpd_storage *contents;

		contents = &(spm->contents[(hash = mapkit_hash((mapkit_hash_t) key)) % spm->size]);

		if ((contents->state == MAPKIT_FULLSLOT)
		    && ((contents->key) == (key))) {
#ifdef MAPKIT_COLLISIONS
			spm->keyindexs++;
#endif
			return contents->value;
		} else if (spm->alwaysdefault && (contents->state == MAPKIT_FREESLOT)) {
#ifdef MAPKIT_COLLISIONS
			spm->keyindexs++;
#endif
			return spm->defaultvalue;
		} else
			return map_h_vpd_value_s(spm, key, hash);
	}

	mapkit_error map_h_vpd_get(map_h_vpd * spm, void *key, double *value) {
		mapkit_hash_t hash;
		map_h_vpd_storage *contents;

		contents = &(spm->contents[(hash = mapkit_hash((mapkit_hash_t) key)) % spm->size]);

		if ((contents->state == MAPKIT_FULLSLOT)
		    && ((contents->key) == (key))) {
#ifdef MAPKIT_COLLISIONS
			spm->keyindexs++;
#endif
			*value = contents->value;
		} else if (spm->alwaysdefault && (contents->state == MAPKIT_FREESLOT)) {
#ifdef MAPKIT_COLLISIONS
			spm->keyindexs++;
#endif
			*value = spm->defaultvalue;
		} else
			return map_h_vpd_get_s(spm, key, value, hash);

		return MAPKIT_OK;
	}

	double *map_h_vpd_insertptr(map_h_vpd * spm, void *key) {
		mapkit_hash_t hash;
		map_h_vpd_storage *contents;

		contents = &(spm->contents[(hash = mapkit_hash((mapkit_hash_t) key)) % spm->size]);

		if (contents->state == MAPKIT_FREESLOT) {
			if (spm->fill >= spm->maxfill)
				return map_h_vpd_insertptr_s(spm, key, hash);

			contents->state = MAPKIT_FULLSLOT;
			contents->key = key;
			contents->value = spm->defaultvalue;
#ifdef MAPKIT_COLLISIONS
			spm->insertionindexs++;
#endif
			spm->used++;
			spm->fill++;
			return &(contents->value);
		} else if ((contents->state == MAPKIT_FULLSLOT)
			   && ((contents->key) == (key))) {
#ifdef MAPKIT_COLLISIONS
			spm->insertionindexs++;
#endif
			return &(contents->value);
		} else
			return map_h_vpd_insertptr_s(spm, key, hash);
	}

	double *map_h_vpd_ptr(map_h_vpd * spm, void *key) {
		mapkit_hash_t hash;
		map_h_vpd_storage *contents;

		contents = &(spm->contents[(hash = mapkit_hash((mapkit_hash_t) key)) % spm->size]);

		if ((contents->state == MAPKIT_FULLSLOT)
		    && ((contents->key) == (key))) {
#ifdef MAPKIT_COLLISIONS
			spm->insertionindexs++;
#endif
			if ((!spm->alwaysdefault) || (!((contents->value) == (spm->defaultvalue))))
				return &(contents->value);
			else
				return NULL;
		} else if (contents->state == MAPKIT_FREESLOT) {
#ifdef MAPKIT_COLLISIONS
			spm->insertionindexs++;
#endif
			return NULL;
		} else
			return map_h_vpd_ptr_s(spm, key, hash);
	}

/*
  Prototypes for map_h_vpvp (void* -> void*)
  Default value : NULL
  Uses a state field.
*/

	typedef struct {
		signed char state;
		void *key;
		void *value;
	} map_h_vpvp_storage;

	typedef struct {
		void *key;
		void *value;
	} map_h_vpvp_element;

	typedef struct {
		mapkit_size_t size,			       /* size of the hash table */
		 fill,					       /* # of non-free slots in the hash table */
		 used,					       /* # of full slots */
		 maxfill,				       /* max # of non-free slots before rehash */
		 minused /* min # of full slots before rehash */ ;

		double maxfillfactor,			       /* maxfill = maxfillfactor * size */
		 minusedfactor;				       /* minused = minusedfactor * size */

#ifdef MAPKIT_COLLISIONS
		unsigned long insertionindexs, insertionindex_collisions;
		unsigned long keyindexs, keyindex_collisions;
#endif

		void *defaultvalue;

		int alwaysdefault;

		map_h_vpvp_storage *contents;
	} map_h_vpvp;

/* Macros */

/* set to inform the type is available */
#define MAPKIT_map_h_vpvp

/* read/write access to the value at key (insert as needed) */
#define map_h_vpvp_(spm, key) (*(map_h_vpvp_insertptr(spm, key)))

/* true if key exists */
#define map_h_vpvp_haskey(spm, key) ((map_h_vpvp_ptr(spm, key)) != NULL)

/* Prototypes */

/* initialize the structure */
	extern mapkit_error map_h_vpvp_init(map_h_vpvp * spm);

/* free contents */
	extern void map_h_vpvp_free(map_h_vpvp * spm);

/* copy the "from" map into a new "to" map */
	extern mapkit_error map_h_vpvp_copy(map_h_vpvp * to, map_h_vpvp * from);

/* initialize the structure for "used" elements */
	extern mapkit_error map_h_vpvp_init_hint(map_h_vpvp * spm, mapkit_size_t used);

/* ensure at least "used" slots (including currently occupied slots) */
	extern mapkit_error map_h_vpvp_ensurecapacity(map_h_vpvp * spm, mapkit_size_t used);

/* restore maxfill & minused and ensure these constraints */
	extern mapkit_error map_h_vpvp_adjustcapacity(map_h_vpvp * spm);

/* reallocate contents for "size" slots, keeping previous values */
	extern mapkit_error map_h_vpvp_reallocate(map_h_vpvp * spm, mapkit_size_t size);

/* return an adequate size for growing from "used" elements */
	extern mapkit_size_t map_h_vpvp_growsize(map_h_vpvp * spm, mapkit_size_t used);

/* return an adequate size for shrinking from "used" elements */
	extern mapkit_size_t map_h_vpvp_shrinksize(map_h_vpvp * spm, mapkit_size_t used);

/* return an adequate size for "used" elements */
	extern mapkit_size_t map_h_vpvp_meansize(map_h_vpvp * spm, mapkit_size_t used);

/* pointer to the value at key (insert as needed) */
	INLINE void **map_h_vpvp_insertptr(map_h_vpvp * spm, void *key);

/* pointer to the value at key or NULL */
	INLINE void **map_h_vpvp_ptr(map_h_vpvp * spm, void *key);

/* remove the key pointing to the value *ptr. No checking is done */
	INLINE mapkit_error map_h_vpvp_removeptr(map_h_vpvp * spm, void **ptr);

/* set the value at key (insert as needed) */
	INLINE mapkit_error map_h_vpvp_set(map_h_vpvp * spm, void *key, void *value);

/* return the value at key (key must exists if alwaysdefaultvalue is not set) */
	INLINE void *map_h_vpvp_value(map_h_vpvp * spm, void *key);

/* get the value at key (return an error code) */
	INLINE mapkit_error map_h_vpvp_get(map_h_vpvp * spm, void *key, void **value);

/* remove key (key must exists) */
	INLINE mapkit_error map_h_vpvp_remove(map_h_vpvp * spm, void *key);

/*
 * return the next index to a full slot of the map, or -1 if none
 * if index == -1, returns the first index, or -1 if none
 * the contents of the slot can be accessed by spm->contents[index].key and
 * spm->contents[index].value
 */
	extern mapkit_size_t map_h_vpvp_next(map_h_vpvp * spm, mapkit_size_t iindex);

/*
 * return the next ptr to a full slot of the map, or NULL if none
 * if pos_contents == NULL, returns the first full slot, or NULL if none
 * the contents of the slot can be accessed by pos_contents->key and
 * pos_contents->value
 */
	extern map_h_vpvp_storage *map_h_vpvp_nextptr(map_h_vpvp * spm, map_h_vpvp_storage * pos_contents);

/* allocate an array of map_h_vpvp_element with all "count" (key,value) pairs */
	extern mapkit_error map_h_vpvp_getall(map_h_vpvp * spm, map_h_vpvp_element ** array, mapkit_size_t * count);

/* remove all values equal to the defaultvalue */
	extern mapkit_error map_h_vpvp_clean(map_h_vpvp * spm);

/* compare the key of two map_h_vpvp_element (for qsort) */
	extern int map_h_vpvp_compare(const void *e1, const void *e2);

/* allocate an array of map_h_vpvp_element with all "count" (key,value) pairs */
	extern mapkit_error map_h_vpvp_getall_sorted(map_h_vpvp * spm, map_h_vpvp_element ** array, mapkit_size_t * count);

/* insert all count map_h_vpvp_element into spm */
	extern mapkit_error map_h_vpvp_setall(map_h_vpvp * spm, map_h_vpvp_element * array, mapkit_size_t count);

/* remove all "count" elements whose keys are in "array" */
	extern mapkit_error map_h_vpvp_removeall(map_h_vpvp * spm, void **array, mapkit_size_t count);

/* print statistics */
	extern void map_h_vpvp_printstats(map_h_vpvp * spm);

/* Private functions */

	extern mapkit_error map_h_vpvp_remove_s(map_h_vpvp * spm, void *key, mapkit_hash_t hash);
	extern void *map_h_vpvp_value_s(map_h_vpvp * spm, void *key, mapkit_hash_t hash);
	extern mapkit_error map_h_vpvp_get_s(map_h_vpvp * spm, void *key, void **value, mapkit_hash_t hash);
	extern mapkit_error map_h_vpvp_set_s(map_h_vpvp * spm, void *key, void *value, mapkit_hash_t hash);
	extern void **map_h_vpvp_insertptr_s(map_h_vpvp * spm, void *key, mapkit_hash_t hash);
	extern void **map_h_vpvp_ptr_s(map_h_vpvp * spm, void *key, mapkit_hash_t hash);

/* Implementation */

/* INLINEd functions */
	mapkit_error map_h_vpvp_remove(map_h_vpvp * spm, void *key) {
		mapkit_hash_t hash;
		map_h_vpvp_storage *contents;

		contents = &(spm->contents[(hash = mapkit_hash((mapkit_hash_t) key)) % spm->size]);

		if ((contents->state == MAPKIT_FULLSLOT)
		    && ((contents->key) == (key))) {
			contents->state = MAPKIT_DELETEDSLOT;
#ifdef MAPKIT_COLLISIONS
			spm->keyindexs++;
#endif
			spm->used--;
			if (spm->used < spm->minused) {
#ifdef MAPKIT_DEBUG
				fprintf(stderr, "MAPKIT: used < minused\n");
#endif
				return map_h_vpvp_reallocate(spm, map_h_vpvp_shrinksize(spm, spm->used));
			}
		} else
			return map_h_vpvp_remove_s(spm, key, hash);
		return MAPKIT_OK;
	}

	mapkit_error map_h_vpvp_removeptr(map_h_vpvp * spm, void **ptr) {
		map_h_vpvp_storage *sptr = (map_h_vpvp_storage *) ((char *) ptr - offsetof(map_h_vpvp_storage, value));

		sptr->state = MAPKIT_DELETEDSLOT;
		spm->used--;
		if (spm->used < spm->minused) {
#ifdef MAPKIT_DEBUG
			fprintf(stderr, "MAPKIT: used < minused\n");
#endif
			return map_h_vpvp_reallocate(spm, map_h_vpvp_shrinksize(spm, spm->used));
		}

		return MAPKIT_OK;
	}

	mapkit_error map_h_vpvp_set(map_h_vpvp * spm, void *key, void *value) {
		mapkit_hash_t hash;
		map_h_vpvp_storage *contents;

		if (spm->alwaysdefault && ((value) == (spm->defaultvalue)))
			return map_h_vpvp_remove(spm, key);

		contents = &(spm->contents[(hash = mapkit_hash((mapkit_hash_t) key)) % spm->size]);

		if (contents->state == MAPKIT_FREESLOT) {
			contents->state = MAPKIT_FULLSLOT;
			contents->key = key;
			contents->value = value;
#ifdef MAPKIT_COLLISIONS
			spm->insertionindexs++;
#endif
			spm->used++;
			spm->fill++;
			if (spm->fill > spm->maxfill) {
#ifdef MAPKIT_DEBUG
				fprintf(stderr, "MAPKIT: fill > maxfill\n");
#endif
				return map_h_vpvp_reallocate(spm, map_h_vpvp_growsize(spm, spm->used));
			}
			return MAPKIT_OK;
		} else if ((contents->state == MAPKIT_FULLSLOT)
			   && ((contents->key) == (key))) {
			contents->value = value;
#ifdef MAPKIT_COLLISIONS
			spm->insertionindexs++;
#endif
			return MAPKIT_OK;
		} else
			return map_h_vpvp_set_s(spm, key, value, hash);
	}

	void *map_h_vpvp_value(map_h_vpvp * spm, void *key) {
		mapkit_hash_t hash;
		map_h_vpvp_storage *contents;

		contents = &(spm->contents[(hash = mapkit_hash((mapkit_hash_t) key)) % spm->size]);

		if ((contents->state == MAPKIT_FULLSLOT)
		    && ((contents->key) == (key))) {
#ifdef MAPKIT_COLLISIONS
			spm->keyindexs++;
#endif
			return contents->value;
		} else if (spm->alwaysdefault && (contents->state == MAPKIT_FREESLOT)) {
#ifdef MAPKIT_COLLISIONS
			spm->keyindexs++;
#endif
			return spm->defaultvalue;
		} else
			return map_h_vpvp_value_s(spm, key, hash);
	}

	mapkit_error map_h_vpvp_get(map_h_vpvp * spm, void *key, void **value) {
		mapkit_hash_t hash;
		map_h_vpvp_storage *contents;

		contents = &(spm->contents[(hash = mapkit_hash((mapkit_hash_t) key)) % spm->size]);

		if ((contents->state == MAPKIT_FULLSLOT)
		    && ((contents->key) == (key))) {
#ifdef MAPKIT_COLLISIONS
			spm->keyindexs++;
#endif
			*value = contents->value;
		} else if (spm->alwaysdefault && (contents->state == MAPKIT_FREESLOT)) {
#ifdef MAPKIT_COLLISIONS
			spm->keyindexs++;
#endif
			*value = spm->defaultvalue;
		} else
			return map_h_vpvp_get_s(spm, key, value, hash);

		return MAPKIT_OK;
	}

	void **map_h_vpvp_insertptr(map_h_vpvp * spm, void *key) {
		mapkit_hash_t hash;
		map_h_vpvp_storage *contents;

		contents = &(spm->contents[(hash = mapkit_hash((mapkit_hash_t) key)) % spm->size]);

		if (contents->state == MAPKIT_FREESLOT) {
			if (spm->fill >= spm->maxfill)
				return map_h_vpvp_insertptr_s(spm, key, hash);

			contents->state = MAPKIT_FULLSLOT;
			contents->key = key;
			contents->value = spm->defaultvalue;
#ifdef MAPKIT_COLLISIONS
			spm->insertionindexs++;
#endif
			spm->used++;
			spm->fill++;
			return &(contents->value);
		} else if ((contents->state == MAPKIT_FULLSLOT)
			   && ((contents->key) == (key))) {
#ifdef MAPKIT_COLLISIONS
			spm->insertionindexs++;
#endif
			return &(contents->value);
		} else
			return map_h_vpvp_insertptr_s(spm, key, hash);
	}

	void **map_h_vpvp_ptr(map_h_vpvp * spm, void *key) {
		mapkit_hash_t hash;
		map_h_vpvp_storage *contents;

		contents = &(spm->contents[(hash = mapkit_hash((mapkit_hash_t) key)) % spm->size]);

		if ((contents->state == MAPKIT_FULLSLOT)
		    && ((contents->key) == (key))) {
#ifdef MAPKIT_COLLISIONS
			spm->insertionindexs++;
#endif
			if ((!spm->alwaysdefault) || (!((contents->value) == (spm->defaultvalue))))
				return &(contents->value);
			else
				return NULL;
		} else if (contents->state == MAPKIT_FREESLOT) {
#ifdef MAPKIT_COLLISIONS
			spm->insertionindexs++;
#endif
			return NULL;
		} else
			return map_h_vpvp_ptr_s(spm, key, hash);
	}

/*
  Prototypes for map_stri (char* -> int)
  Default value : 0
  Uses a state field.
*/

	typedef struct {
		signed char state;
		char *key;
		int value;
	} map_stri_storage;

	typedef struct {
		char *key;
		int value;
	} map_stri_element;

	typedef struct {
		mapkit_size_t size,			       /* size of the hash table */
		 fill,					       /* # of non-free slots in the hash table */
		 used,					       /* # of full slots */
		 maxfill,				       /* max # of non-free slots before rehash */
		 minused /* min # of full slots before rehash */ ;

		double maxfillfactor,			       /* maxfill = maxfillfactor * size */
		 minusedfactor;				       /* minused = minusedfactor * size */

#ifdef MAPKIT_COLLISIONS
		unsigned long insertionindexs, insertionindex_collisions;
		unsigned long keyindexs, keyindex_collisions;
#endif

		int defaultvalue;

		int alwaysdefault;

		map_stri_storage *contents;
	} map_stri;

/* Macros */

/* set to inform the type is available */
#define MAPKIT_map_stri

/* read/write access to the value at key (insert as needed) */
#define map_stri_(spm, key) (*(map_stri_insertptr(spm, key)))

/* true if key exists */
#define map_stri_haskey(spm, key) ((map_stri_ptr(spm, key)) != NULL)

/* Prototypes */

/* initialize the structure */
	extern mapkit_error map_stri_init(map_stri * spm);

/* free contents */
	extern void map_stri_free(map_stri * spm);

/* copy the "from" map into a new "to" map */
	extern mapkit_error map_stri_copy(map_stri * to, map_stri * from);

/* initialize the structure for "used" elements */
	extern mapkit_error map_stri_init_hint(map_stri * spm, mapkit_size_t used);

/* ensure at least "used" slots (including currently occupied slots) */
	extern mapkit_error map_stri_ensurecapacity(map_stri * spm, mapkit_size_t used);

/* restore maxfill & minused and ensure these constraints */
	extern mapkit_error map_stri_adjustcapacity(map_stri * spm);

/* reallocate contents for "size" slots, keeping previous values */
	extern mapkit_error map_stri_reallocate(map_stri * spm, mapkit_size_t size);

/* return an adequate size for growing from "used" elements */
	extern mapkit_size_t map_stri_growsize(map_stri * spm, mapkit_size_t used);

/* return an adequate size for shrinking from "used" elements */
	extern mapkit_size_t map_stri_shrinksize(map_stri * spm, mapkit_size_t used);

/* return an adequate size for "used" elements */
	extern mapkit_size_t map_stri_meansize(map_stri * spm, mapkit_size_t used);

/* pointer to the value at key (insert as needed) */
	INLINE int *map_stri_insertptr(map_stri * spm, char *key);

/* pointer to the value at key or NULL */
	INLINE int *map_stri_ptr(map_stri * spm, char *key);

/* remove the key pointing to the value *ptr. No checking is done */
	INLINE mapkit_error map_stri_removeptr(map_stri * spm, int *ptr);

/* set the value at key (insert as needed) */
	INLINE mapkit_error map_stri_set(map_stri * spm, char *key, int value);

/* return the value at key (key must exists if alwaysdefaultvalue is not set) */
	INLINE int map_stri_value(map_stri * spm, char *key);

/* get the value at key (return an error code) */
	INLINE mapkit_error map_stri_get(map_stri * spm, char *key, int *value);

/* remove key (key must exists) */
	INLINE mapkit_error map_stri_remove(map_stri * spm, char *key);

/*
 * return the next index to a full slot of the map, or -1 if none
 * if index == -1, returns the first index, or -1 if none
 * the contents of the slot can be accessed by spm->contents[index].key and
 * spm->contents[index].value
 */
	extern mapkit_size_t map_stri_next(map_stri * spm, mapkit_size_t iindex);

/*
 * return the next ptr to a full slot of the map, or NULL if none
 * if pos_contents == NULL, returns the first full slot, or NULL if none
 * the contents of the slot can be accessed by pos_contents->key and
 * pos_contents->value
 */
	extern map_stri_storage *map_stri_nextptr(map_stri * spm, map_stri_storage * pos_contents);

/* allocate an array of map_stri_element with all "count" (key,value) pairs */
	extern mapkit_error map_stri_getall(map_stri * spm, map_stri_element ** array, mapkit_size_t * count);

/* remove all values equal to the defaultvalue */
	extern mapkit_error map_stri_clean(map_stri * spm);

/* compare the key of two map_stri_element (for qsort) */
	extern int map_stri_compare(const void *e1, const void *e2);

/* allocate an array of map_stri_element with all "count" (key,value) pairs */
	extern mapkit_error map_stri_getall_sorted(map_stri * spm, map_stri_element ** array, mapkit_size_t * count);

/* insert all count map_stri_element into spm */
	extern mapkit_error map_stri_setall(map_stri * spm, map_stri_element * array, mapkit_size_t count);

/* remove all "count" elements whose keys are in "array" */
	extern mapkit_error map_stri_removeall(map_stri * spm, char **array, mapkit_size_t count);

/* print statistics */
	extern void map_stri_printstats(map_stri * spm);

/* Private functions */

	extern mapkit_error map_stri_remove_s(map_stri * spm, char *key, mapkit_hash_t hash);
	extern int map_stri_value_s(map_stri * spm, char *key, mapkit_hash_t hash);
	extern mapkit_error map_stri_get_s(map_stri * spm, char *key, int *value, mapkit_hash_t hash);
	extern mapkit_error map_stri_set_s(map_stri * spm, char *key, int value, mapkit_hash_t hash);
	extern int *map_stri_insertptr_s(map_stri * spm, char *key, mapkit_hash_t hash);
	extern int *map_stri_ptr_s(map_stri * spm, char *key, mapkit_hash_t hash);

/* Implementation */

/* INLINEd functions */
	mapkit_error map_stri_remove(map_stri * spm, char *key) {
		mapkit_hash_t hash;
		map_stri_storage *contents;

		contents = &(spm->contents[(hash = mapkit_strhash(key)) % spm->size]);

		if ((contents->state == MAPKIT_FULLSLOT)
		    && (strcmp(contents->key, key) == 0)) {
			contents->state = MAPKIT_DELETEDSLOT;
#ifdef MAPKIT_COLLISIONS
			spm->keyindexs++;
#endif
			spm->used--;
			if (spm->used < spm->minused) {
#ifdef MAPKIT_DEBUG
				fprintf(stderr, "MAPKIT: used < minused\n");
#endif
				return map_stri_reallocate(spm, map_stri_shrinksize(spm, spm->used));
			}
		} else
			return map_stri_remove_s(spm, key, hash);
		return MAPKIT_OK;
	}

	mapkit_error map_stri_removeptr(map_stri * spm, int *ptr) {
		map_stri_storage *sptr = (map_stri_storage *) ((char *) ptr - offsetof(map_stri_storage, value));

		sptr->state = MAPKIT_DELETEDSLOT;
		spm->used--;
		if (spm->used < spm->minused) {
#ifdef MAPKIT_DEBUG
			fprintf(stderr, "MAPKIT: used < minused\n");
#endif
			return map_stri_reallocate(spm, map_stri_shrinksize(spm, spm->used));
		}

		return MAPKIT_OK;
	}

	mapkit_error map_stri_set(map_stri * spm, char *key, int value) {
		mapkit_hash_t hash;
		map_stri_storage *contents;

		if (spm->alwaysdefault && ((value) == (spm->defaultvalue)))
			return map_stri_remove(spm, key);

		contents = &(spm->contents[(hash = mapkit_strhash(key)) % spm->size]);

		if (contents->state == MAPKIT_FREESLOT) {
			contents->state = MAPKIT_FULLSLOT;
			contents->key = key;
			contents->value = value;
#ifdef MAPKIT_COLLISIONS
			spm->insertionindexs++;
#endif
			spm->used++;
			spm->fill++;
			if (spm->fill > spm->maxfill) {
#ifdef MAPKIT_DEBUG
				fprintf(stderr, "MAPKIT: fill > maxfill\n");
#endif
				return map_stri_reallocate(spm, map_stri_growsize(spm, spm->used));
			}
			return MAPKIT_OK;
		} else if ((contents->state == MAPKIT_FULLSLOT)
			   && (strcmp(contents->key, key) == 0)) {
			contents->value = value;
#ifdef MAPKIT_COLLISIONS
			spm->insertionindexs++;
#endif
			return MAPKIT_OK;
		} else
			return map_stri_set_s(spm, key, value, hash);
	}

	int map_stri_value(map_stri * spm, char *key) {
		mapkit_hash_t hash;
		map_stri_storage *contents;

		contents = &(spm->contents[(hash = mapkit_strhash(key)) % spm->size]);

		if ((contents->state == MAPKIT_FULLSLOT)
		    && (strcmp(contents->key, key) == 0)) {
#ifdef MAPKIT_COLLISIONS
			spm->keyindexs++;
#endif
			return contents->value;
		} else if (spm->alwaysdefault && (contents->state == MAPKIT_FREESLOT)) {
#ifdef MAPKIT_COLLISIONS
			spm->keyindexs++;
#endif
			return spm->defaultvalue;
		} else
			return map_stri_value_s(spm, key, hash);
	}

	mapkit_error map_stri_get(map_stri * spm, char *key, int *value) {
		mapkit_hash_t hash;
		map_stri_storage *contents;

		contents = &(spm->contents[(hash = mapkit_strhash(key)) % spm->size]);

		if ((contents->state == MAPKIT_FULLSLOT)
		    && (strcmp(contents->key, key) == 0)) {
#ifdef MAPKIT_COLLISIONS
			spm->keyindexs++;
#endif
			*value = contents->value;
		} else if (spm->alwaysdefault && (contents->state == MAPKIT_FREESLOT)) {
#ifdef MAPKIT_COLLISIONS
			spm->keyindexs++;
#endif
			*value = spm->defaultvalue;
		} else
			return map_stri_get_s(spm, key, value, hash);

		return MAPKIT_OK;
	}

	int *map_stri_insertptr(map_stri * spm, char *key) {
		mapkit_hash_t hash;
		map_stri_storage *contents;

		contents = &(spm->contents[(hash = mapkit_strhash(key)) % spm->size]);

		if (contents->state == MAPKIT_FREESLOT) {
			if (spm->fill >= spm->maxfill)
				return map_stri_insertptr_s(spm, key, hash);

			contents->state = MAPKIT_FULLSLOT;
			contents->key = key;
			contents->value = spm->defaultvalue;
#ifdef MAPKIT_COLLISIONS
			spm->insertionindexs++;
#endif
			spm->used++;
			spm->fill++;
			return &(contents->value);
		} else if ((contents->state == MAPKIT_FULLSLOT)
			   && (strcmp(contents->key, key) == 0)) {
#ifdef MAPKIT_COLLISIONS
			spm->insertionindexs++;
#endif
			return &(contents->value);
		} else
			return map_stri_insertptr_s(spm, key, hash);
	}

	int *map_stri_ptr(map_stri * spm, char *key) {
		mapkit_hash_t hash;
		map_stri_storage *contents;

		contents = &(spm->contents[(hash = mapkit_strhash(key)) % spm->size]);

		if ((contents->state == MAPKIT_FULLSLOT)
		    && (strcmp(contents->key, key) == 0)) {
#ifdef MAPKIT_COLLISIONS
			spm->insertionindexs++;
#endif
			if ((!spm->alwaysdefault) || (!((contents->value) == (spm->defaultvalue))))
				return &(contents->value);
			else
				return NULL;
		} else if (contents->state == MAPKIT_FREESLOT) {
#ifdef MAPKIT_COLLISIONS
			spm->insertionindexs++;
#endif
			return NULL;
		} else
			return map_stri_ptr_s(spm, key, hash);
	}

/*
  Prototypes for map_strd (char* -> double)
  Default value : 0.0
  Uses a state field.
*/

	typedef struct {
		signed char state;
		char *key;
		double value;
	} map_strd_storage;

	typedef struct {
		char *key;
		double value;
	} map_strd_element;

	typedef struct {
		mapkit_size_t size,			       /* size of the hash table */
		 fill,					       /* # of non-free slots in the hash table */
		 used,					       /* # of full slots */
		 maxfill,				       /* max # of non-free slots before rehash */
		 minused /* min # of full slots before rehash */ ;

		double maxfillfactor,			       /* maxfill = maxfillfactor * size */
		 minusedfactor;				       /* minused = minusedfactor * size */

#ifdef MAPKIT_COLLISIONS
		unsigned long insertionindexs, insertionindex_collisions;
		unsigned long keyindexs, keyindex_collisions;
#endif

		double defaultvalue;

		int alwaysdefault;

		map_strd_storage *contents;
	} map_strd;

/* Macros */

/* set to inform the type is available */
#define MAPKIT_map_strd

/* read/write access to the value at key (insert as needed) */
#define map_strd_(spm, key) (*(map_strd_insertptr(spm, key)))

/* true if key exists */
#define map_strd_haskey(spm, key) ((map_strd_ptr(spm, key)) != NULL)

/* Prototypes */

/* initialize the structure */
	extern mapkit_error map_strd_init(map_strd * spm);

/* free contents */
	extern void map_strd_free(map_strd * spm);

/* copy the "from" map into a new "to" map */
	extern mapkit_error map_strd_copy(map_strd * to, map_strd * from);

/* initialize the structure for "used" elements */
	extern mapkit_error map_strd_init_hint(map_strd * spm, mapkit_size_t used);

/* ensure at least "used" slots (including currently occupied slots) */
	extern mapkit_error map_strd_ensurecapacity(map_strd * spm, mapkit_size_t used);

/* restore maxfill & minused and ensure these constraints */
	extern mapkit_error map_strd_adjustcapacity(map_strd * spm);

/* reallocate contents for "size" slots, keeping previous values */
	extern mapkit_error map_strd_reallocate(map_strd * spm, mapkit_size_t size);

/* return an adequate size for growing from "used" elements */
	extern mapkit_size_t map_strd_growsize(map_strd * spm, mapkit_size_t used);

/* return an adequate size for shrinking from "used" elements */
	extern mapkit_size_t map_strd_shrinksize(map_strd * spm, mapkit_size_t used);

/* return an adequate size for "used" elements */
	extern mapkit_size_t map_strd_meansize(map_strd * spm, mapkit_size_t used);

/* pointer to the value at key (insert as needed) */
	INLINE double *map_strd_insertptr(map_strd * spm, char *key);

/* pointer to the value at key or NULL */
	INLINE double *map_strd_ptr(map_strd * spm, char *key);

/* remove the key pointing to the value *ptr. No checking is done */
	INLINE mapkit_error map_strd_removeptr(map_strd * spm, double *ptr);

/* set the value at key (insert as needed) */
	INLINE mapkit_error map_strd_set(map_strd * spm, char *key, double value);

/* return the value at key (key must exists if alwaysdefaultvalue is not set) */
	INLINE double map_strd_value(map_strd * spm, char *key);

/* get the value at key (return an error code) */
	INLINE mapkit_error map_strd_get(map_strd * spm, char *key, double *value);

/* remove key (key must exists) */
	INLINE mapkit_error map_strd_remove(map_strd * spm, char *key);

/*
 * return the next index to a full slot of the map, or -1 if none
 * if index == -1, returns the first index, or -1 if none
 * the contents of the slot can be accessed by spm->contents[index].key and
 * spm->contents[index].value
 */
	extern mapkit_size_t map_strd_next(map_strd * spm, mapkit_size_t iindex);

/*
 * return the next ptr to a full slot of the map, or NULL if none
 * if pos_contents == NULL, returns the first full slot, or NULL if none
 * the contents of the slot can be accessed by pos_contents->key and
 * pos_contents->value
 */
	extern map_strd_storage *map_strd_nextptr(map_strd * spm, map_strd_storage * pos_contents);

/* allocate an array of map_strd_element with all "count" (key,value) pairs */
	extern mapkit_error map_strd_getall(map_strd * spm, map_strd_element ** array, mapkit_size_t * count);

/* remove all values equal to the defaultvalue */
	extern mapkit_error map_strd_clean(map_strd * spm);

/* compare the key of two map_strd_element (for qsort) */
	extern int map_strd_compare(const void *e1, const void *e2);

/* allocate an array of map_strd_element with all "count" (key,value) pairs */
	extern mapkit_error map_strd_getall_sorted(map_strd * spm, map_strd_element ** array, mapkit_size_t * count);

/* insert all count map_strd_element into spm */
	extern mapkit_error map_strd_setall(map_strd * spm, map_strd_element * array, mapkit_size_t count);

/* remove all "count" elements whose keys are in "array" */
	extern mapkit_error map_strd_removeall(map_strd * spm, char **array, mapkit_size_t count);

/* print statistics */
	extern void map_strd_printstats(map_strd * spm);

/* Private functions */

	extern mapkit_error map_strd_remove_s(map_strd * spm, char *key, mapkit_hash_t hash);
	extern double map_strd_value_s(map_strd * spm, char *key, mapkit_hash_t hash);
	extern mapkit_error map_strd_get_s(map_strd * spm, char *key, double *value, mapkit_hash_t hash);
	extern mapkit_error map_strd_set_s(map_strd * spm, char *key, double value, mapkit_hash_t hash);
	extern double *map_strd_insertptr_s(map_strd * spm, char *key, mapkit_hash_t hash);
	extern double *map_strd_ptr_s(map_strd * spm, char *key, mapkit_hash_t hash);

/* Implementation */

/* INLINEd functions */
	mapkit_error map_strd_remove(map_strd * spm, char *key) {
		mapkit_hash_t hash;
		map_strd_storage *contents;

		contents = &(spm->contents[(hash = mapkit_strhash(key)) % spm->size]);

		if ((contents->state == MAPKIT_FULLSLOT)
		    && (strcmp(contents->key, key) == 0)) {
			contents->state = MAPKIT_DELETEDSLOT;
#ifdef MAPKIT_COLLISIONS
			spm->keyindexs++;
#endif
			spm->used--;
			if (spm->used < spm->minused) {
#ifdef MAPKIT_DEBUG
				fprintf(stderr, "MAPKIT: used < minused\n");
#endif
				return map_strd_reallocate(spm, map_strd_shrinksize(spm, spm->used));
			}
		} else
			return map_strd_remove_s(spm, key, hash);
		return MAPKIT_OK;
	}

	mapkit_error map_strd_removeptr(map_strd * spm, double *ptr) {
		map_strd_storage *sptr = (map_strd_storage *) ((char *) ptr - offsetof(map_strd_storage, value));

		sptr->state = MAPKIT_DELETEDSLOT;
		spm->used--;
		if (spm->used < spm->minused) {
#ifdef MAPKIT_DEBUG
			fprintf(stderr, "MAPKIT: used < minused\n");
#endif
			return map_strd_reallocate(spm, map_strd_shrinksize(spm, spm->used));
		}

		return MAPKIT_OK;
	}

	mapkit_error map_strd_set(map_strd * spm, char *key, double value) {
		mapkit_hash_t hash;
		map_strd_storage *contents;

		if (spm->alwaysdefault && ((value) == (spm->defaultvalue)))
			return map_strd_remove(spm, key);

		contents = &(spm->contents[(hash = mapkit_strhash(key)) % spm->size]);

		if (contents->state == MAPKIT_FREESLOT) {
			contents->state = MAPKIT_FULLSLOT;
			contents->key = key;
			contents->value = value;
#ifdef MAPKIT_COLLISIONS
			spm->insertionindexs++;
#endif
			spm->used++;
			spm->fill++;
			if (spm->fill > spm->maxfill) {
#ifdef MAPKIT_DEBUG
				fprintf(stderr, "MAPKIT: fill > maxfill\n");
#endif
				return map_strd_reallocate(spm, map_strd_growsize(spm, spm->used));
			}
			return MAPKIT_OK;
		} else if ((contents->state == MAPKIT_FULLSLOT)
			   && (strcmp(contents->key, key) == 0)) {
			contents->value = value;
#ifdef MAPKIT_COLLISIONS
			spm->insertionindexs++;
#endif
			return MAPKIT_OK;
		} else
			return map_strd_set_s(spm, key, value, hash);
	}

	double map_strd_value(map_strd * spm, char *key) {
		mapkit_hash_t hash;
		map_strd_storage *contents;

		contents = &(spm->contents[(hash = mapkit_strhash(key)) % spm->size]);

		if ((contents->state == MAPKIT_FULLSLOT)
		    && (strcmp(contents->key, key) == 0)) {
#ifdef MAPKIT_COLLISIONS
			spm->keyindexs++;
#endif
			return contents->value;
		} else if (spm->alwaysdefault && (contents->state == MAPKIT_FREESLOT)) {
#ifdef MAPKIT_COLLISIONS
			spm->keyindexs++;
#endif
			return spm->defaultvalue;
		} else
			return map_strd_value_s(spm, key, hash);
	}

	mapkit_error map_strd_get(map_strd * spm, char *key, double *value) {
		mapkit_hash_t hash;
		map_strd_storage *contents;

		contents = &(spm->contents[(hash = mapkit_strhash(key)) % spm->size]);

		if ((contents->state == MAPKIT_FULLSLOT)
		    && (strcmp(contents->key, key) == 0)) {
#ifdef MAPKIT_COLLISIONS
			spm->keyindexs++;
#endif
			*value = contents->value;
		} else if (spm->alwaysdefault && (contents->state == MAPKIT_FREESLOT)) {
#ifdef MAPKIT_COLLISIONS
			spm->keyindexs++;
#endif
			*value = spm->defaultvalue;
		} else
			return map_strd_get_s(spm, key, value, hash);

		return MAPKIT_OK;
	}

	double *map_strd_insertptr(map_strd * spm, char *key) {
		mapkit_hash_t hash;
		map_strd_storage *contents;

		contents = &(spm->contents[(hash = mapkit_strhash(key)) % spm->size]);

		if (contents->state == MAPKIT_FREESLOT) {
			if (spm->fill >= spm->maxfill)
				return map_strd_insertptr_s(spm, key, hash);

			contents->state = MAPKIT_FULLSLOT;
			contents->key = key;
			contents->value = spm->defaultvalue;
#ifdef MAPKIT_COLLISIONS
			spm->insertionindexs++;
#endif
			spm->used++;
			spm->fill++;
			return &(contents->value);
		} else if ((contents->state == MAPKIT_FULLSLOT)
			   && (strcmp(contents->key, key) == 0)) {
#ifdef MAPKIT_COLLISIONS
			spm->insertionindexs++;
#endif
			return &(contents->value);
		} else
			return map_strd_insertptr_s(spm, key, hash);
	}

	double *map_strd_ptr(map_strd * spm, char *key) {
		mapkit_hash_t hash;
		map_strd_storage *contents;

		contents = &(spm->contents[(hash = mapkit_strhash(key)) % spm->size]);

		if ((contents->state == MAPKIT_FULLSLOT)
		    && (strcmp(contents->key, key) == 0)) {
#ifdef MAPKIT_COLLISIONS
			spm->insertionindexs++;
#endif
			if ((!spm->alwaysdefault) || (!((contents->value) == (spm->defaultvalue))))
				return &(contents->value);
			else
				return NULL;
		} else if (contents->state == MAPKIT_FREESLOT) {
#ifdef MAPKIT_COLLISIONS
			spm->insertionindexs++;
#endif
			return NULL;
		} else
			return map_strd_ptr_s(spm, key, hash);
	}

/*
  Prototypes for map_strvp (char* -> void*)
  Default value : NULL
  Uses a state field.
*/

	typedef struct {
		signed char state;
		char *key;
		void *value;
	} map_strvp_storage;

	typedef struct {
		char *key;
		void *value;
	} map_strvp_element;

	typedef struct {
		mapkit_size_t size,			       /* size of the hash table */
		 fill,					       /* # of non-free slots in the hash table */
		 used,					       /* # of full slots */
		 maxfill,				       /* max # of non-free slots before rehash */
		 minused /* min # of full slots before rehash */ ;

		double maxfillfactor,			       /* maxfill = maxfillfactor * size */
		 minusedfactor;				       /* minused = minusedfactor * size */

#ifdef MAPKIT_COLLISIONS
		unsigned long insertionindexs, insertionindex_collisions;
		unsigned long keyindexs, keyindex_collisions;
#endif

		void *defaultvalue;

		int alwaysdefault;

		map_strvp_storage *contents;
	} map_strvp;

/* Macros */

/* set to inform the type is available */
#define MAPKIT_map_strvp

/* read/write access to the value at key (insert as needed) */
#define map_strvp_(spm, key) (*(map_strvp_insertptr(spm, key)))

/* true if key exists */
#define map_strvp_haskey(spm, key) ((map_strvp_ptr(spm, key)) != NULL)

/* Prototypes */

/* initialize the structure */
	extern mapkit_error map_strvp_init(map_strvp * spm);

/* free contents */
	extern void map_strvp_free(map_strvp * spm);

/* copy the "from" map into a new "to" map */
	extern mapkit_error map_strvp_copy(map_strvp * to, map_strvp * from);

/* initialize the structure for "used" elements */
	extern mapkit_error map_strvp_init_hint(map_strvp * spm, mapkit_size_t used);

/* ensure at least "used" slots (including currently occupied slots) */
	extern mapkit_error map_strvp_ensurecapacity(map_strvp * spm, mapkit_size_t used);

/* restore maxfill & minused and ensure these constraints */
	extern mapkit_error map_strvp_adjustcapacity(map_strvp * spm);

/* reallocate contents for "size" slots, keeping previous values */
	extern mapkit_error map_strvp_reallocate(map_strvp * spm, mapkit_size_t size);

/* return an adequate size for growing from "used" elements */
	extern mapkit_size_t map_strvp_growsize(map_strvp * spm, mapkit_size_t used);

/* return an adequate size for shrinking from "used" elements */
	extern mapkit_size_t map_strvp_shrinksize(map_strvp * spm, mapkit_size_t used);

/* return an adequate size for "used" elements */
	extern mapkit_size_t map_strvp_meansize(map_strvp * spm, mapkit_size_t used);

/* pointer to the value at key (insert as needed) */
	INLINE void **map_strvp_insertptr(map_strvp * spm, char *key);

/* pointer to the value at key or NULL */
	INLINE void **map_strvp_ptr(map_strvp * spm, char *key);

/* remove the key pointing to the value *ptr. No checking is done */
	INLINE mapkit_error map_strvp_removeptr(map_strvp * spm, void **ptr);

/* set the value at key (insert as needed) */
	INLINE mapkit_error map_strvp_set(map_strvp * spm, char *key, void *value);

/* return the value at key (key must exists if alwaysdefaultvalue is not set) */
	INLINE void *map_strvp_value(map_strvp * spm, char *key);

/* get the value at key (return an error code) */
	INLINE mapkit_error map_strvp_get(map_strvp * spm, char *key, void **value);

/* remove key (key must exists) */
	INLINE mapkit_error map_strvp_remove(map_strvp * spm, char *key);

/*
 * return the next index to a full slot of the map, or -1 if none
 * if index == -1, returns the first index, or -1 if none
 * the contents of the slot can be accessed by spm->contents[index].key and
 * spm->contents[index].value
 */
	extern mapkit_size_t map_strvp_next(map_strvp * spm, mapkit_size_t iindex);

/*
 * return the next ptr to a full slot of the map, or NULL if none
 * if pos_contents == NULL, returns the first full slot, or NULL if none
 * the contents of the slot can be accessed by pos_contents->key and
 * pos_contents->value
 */
	extern map_strvp_storage *map_strvp_nextptr(map_strvp * spm, map_strvp_storage * pos_contents);

/* allocate an array of map_strvp_element with all "count" (key,value) pairs */
	extern mapkit_error map_strvp_getall(map_strvp * spm, map_strvp_element ** array, mapkit_size_t * count);

/* remove all values equal to the defaultvalue */
	extern mapkit_error map_strvp_clean(map_strvp * spm);

/* compare the key of two map_strvp_element (for qsort) */
	extern int map_strvp_compare(const void *e1, const void *e2);

/* allocate an array of map_strvp_element with all "count" (key,value) pairs */
	extern mapkit_error map_strvp_getall_sorted(map_strvp * spm, map_strvp_element ** array, mapkit_size_t * count);

/* insert all count map_strvp_element into spm */
	extern mapkit_error map_strvp_setall(map_strvp * spm, map_strvp_element * array, mapkit_size_t count);

/* remove all "count" elements whose keys are in "array" */
	extern mapkit_error map_strvp_removeall(map_strvp * spm, char **array, mapkit_size_t count);

/* print statistics */
	extern void map_strvp_printstats(map_strvp * spm);

/* Private functions */

	extern mapkit_error map_strvp_remove_s(map_strvp * spm, char *key, mapkit_hash_t hash);
	extern void *map_strvp_value_s(map_strvp * spm, char *key, mapkit_hash_t hash);
	extern mapkit_error map_strvp_get_s(map_strvp * spm, char *key, void **value, mapkit_hash_t hash);
	extern mapkit_error map_strvp_set_s(map_strvp * spm, char *key, void *value, mapkit_hash_t hash);
	extern void **map_strvp_insertptr_s(map_strvp * spm, char *key, mapkit_hash_t hash);
	extern void **map_strvp_ptr_s(map_strvp * spm, char *key, mapkit_hash_t hash);

/* Implementation */

/* INLINEd functions */
	mapkit_error map_strvp_remove(map_strvp * spm, char *key) {
		mapkit_hash_t hash;
		map_strvp_storage *contents;

		contents = &(spm->contents[(hash = mapkit_strhash(key)) % spm->size]);

		if ((contents->state == MAPKIT_FULLSLOT)
		    && (strcmp(contents->key, key) == 0)) {
			contents->state = MAPKIT_DELETEDSLOT;
#ifdef MAPKIT_COLLISIONS
			spm->keyindexs++;
#endif
			spm->used--;
			if (spm->used < spm->minused) {
#ifdef MAPKIT_DEBUG
				fprintf(stderr, "MAPKIT: used < minused\n");
#endif
				return map_strvp_reallocate(spm, map_strvp_shrinksize(spm, spm->used));
			}
		} else
			return map_strvp_remove_s(spm, key, hash);
		return MAPKIT_OK;
	}

	mapkit_error map_strvp_removeptr(map_strvp * spm, void **ptr) {
		map_strvp_storage *sptr = (map_strvp_storage *) ((char *) ptr - offsetof(map_strvp_storage, value));

		sptr->state = MAPKIT_DELETEDSLOT;
		spm->used--;
		if (spm->used < spm->minused) {
#ifdef MAPKIT_DEBUG
			fprintf(stderr, "MAPKIT: used < minused\n");
#endif
			return map_strvp_reallocate(spm, map_strvp_shrinksize(spm, spm->used));
		}

		return MAPKIT_OK;
	}

	mapkit_error map_strvp_set(map_strvp * spm, char *key, void *value) {
		mapkit_hash_t hash;
		map_strvp_storage *contents;

		if (spm->alwaysdefault && ((value) == (spm->defaultvalue)))
			return map_strvp_remove(spm, key);

		contents = &(spm->contents[(hash = mapkit_strhash(key)) % spm->size]);

		if (contents->state == MAPKIT_FREESLOT) {
			contents->state = MAPKIT_FULLSLOT;
			contents->key = key;
			contents->value = value;
#ifdef MAPKIT_COLLISIONS
			spm->insertionindexs++;
#endif
			spm->used++;
			spm->fill++;
			if (spm->fill > spm->maxfill) {
#ifdef MAPKIT_DEBUG
				fprintf(stderr, "MAPKIT: fill > maxfill\n");
#endif
				return map_strvp_reallocate(spm, map_strvp_growsize(spm, spm->used));
			}
			return MAPKIT_OK;
		} else if ((contents->state == MAPKIT_FULLSLOT)
			   && (strcmp(contents->key, key) == 0)) {
			contents->value = value;
#ifdef MAPKIT_COLLISIONS
			spm->insertionindexs++;
#endif
			return MAPKIT_OK;
		} else
			return map_strvp_set_s(spm, key, value, hash);
	}

	void *map_strvp_value(map_strvp * spm, char *key) {
		mapkit_hash_t hash;
		map_strvp_storage *contents;

		contents = &(spm->contents[(hash = mapkit_strhash(key)) % spm->size]);

		if ((contents->state == MAPKIT_FULLSLOT)
		    && (strcmp(contents->key, key) == 0)) {
#ifdef MAPKIT_COLLISIONS
			spm->keyindexs++;
#endif
			return contents->value;
		} else if (spm->alwaysdefault && (contents->state == MAPKIT_FREESLOT)) {
#ifdef MAPKIT_COLLISIONS
			spm->keyindexs++;
#endif
			return spm->defaultvalue;
		} else
			return map_strvp_value_s(spm, key, hash);
	}

	mapkit_error map_strvp_get(map_strvp * spm, char *key, void **value) {
		mapkit_hash_t hash;
		map_strvp_storage *contents;

		contents = &(spm->contents[(hash = mapkit_strhash(key)) % spm->size]);

		if ((contents->state == MAPKIT_FULLSLOT)
		    && (strcmp(contents->key, key) == 0)) {
#ifdef MAPKIT_COLLISIONS
			spm->keyindexs++;
#endif
			*value = contents->value;
		} else if (spm->alwaysdefault && (contents->state == MAPKIT_FREESLOT)) {
#ifdef MAPKIT_COLLISIONS
			spm->keyindexs++;
#endif
			*value = spm->defaultvalue;
		} else
			return map_strvp_get_s(spm, key, value, hash);

		return MAPKIT_OK;
	}

	void **map_strvp_insertptr(map_strvp * spm, char *key) {
		mapkit_hash_t hash;
		map_strvp_storage *contents;

		contents = &(spm->contents[(hash = mapkit_strhash(key)) % spm->size]);

		if (contents->state == MAPKIT_FREESLOT) {
			if (spm->fill >= spm->maxfill)
				return map_strvp_insertptr_s(spm, key, hash);

			contents->state = MAPKIT_FULLSLOT;
			contents->key = key;
			contents->value = spm->defaultvalue;
#ifdef MAPKIT_COLLISIONS
			spm->insertionindexs++;
#endif
			spm->used++;
			spm->fill++;
			return &(contents->value);
		} else if ((contents->state == MAPKIT_FULLSLOT)
			   && (strcmp(contents->key, key) == 0)) {
#ifdef MAPKIT_COLLISIONS
			spm->insertionindexs++;
#endif
			return &(contents->value);
		} else
			return map_strvp_insertptr_s(spm, key, hash);
	}
	void **map_strvp_ptr(map_strvp * spm, char *key) {
		mapkit_hash_t hash;
		map_strvp_storage *contents;

		contents = &(spm->contents[(hash = mapkit_strhash(key)) % spm->size]);

		if ((contents->state == MAPKIT_FULLSLOT)
		    && (strcmp(contents->key, key) == 0)) {
#ifdef MAPKIT_COLLISIONS
			spm->insertionindexs++;
#endif
			if ((!spm->alwaysdefault) || (!((contents->value) == (spm->defaultvalue))))
				return &(contents->value);
			else
				return NULL;
		} else if (contents->state == MAPKIT_FREESLOT) {
#ifdef MAPKIT_COLLISIONS
			spm->insertionindexs++;
#endif
			return NULL;
		} else
			return map_strvp_ptr_s(spm, key, hash);
	}

/*
  Prototypes for map_strstr (char* -> char*)
  Default value : ""
  Uses a state field.
*/

	typedef struct {
		signed char state;
		char *key;
		char *value;
	} map_strstr_storage;

	typedef struct {
		char *key;
		char *value;
	} map_strstr_element;

	typedef struct {
		mapkit_size_t size,			       /* size of the hash table */
		 fill,					       /* # of non-free slots in the hash table */
		 used,					       /* # of full slots */
		 maxfill,				       /* max # of non-free slots before rehash */
		 minused /* min # of full slots before rehash */ ;

		double maxfillfactor,			       /* maxfill = maxfillfactor * size */
		 minusedfactor;				       /* minused = minusedfactor * size */

#ifdef MAPKIT_COLLISIONS
		unsigned long insertionindexs, insertionindex_collisions;
		unsigned long keyindexs, keyindex_collisions;
#endif

		char *defaultvalue;

		int alwaysdefault;

		map_strstr_storage *contents;
	} map_strstr;

/* Macros */

/* set to inform the type is available */
#define MAPKIT_map_strstr

/* read/write access to the value at key (insert as needed) */
#define map_strstr_(spm, key) (*(map_strstr_insertptr(spm, key)))

/* true if key exists */
#define map_strstr_haskey(spm, key) ((map_strstr_ptr(spm, key)) != NULL)

/* Prototypes */

/* initialize the structure */
	extern mapkit_error map_strstr_init(map_strstr * spm);

/* free contents */
	extern void map_strstr_free(map_strstr * spm);

/* copy the "from" map into a new "to" map */
	extern mapkit_error map_strstr_copy(map_strstr * to, map_strstr * from);

/* initialize the structure for "used" elements */
	extern mapkit_error map_strstr_init_hint(map_strstr * spm, mapkit_size_t used);

/* ensure at least "used" slots (including currently occupied slots) */
	extern mapkit_error map_strstr_ensurecapacity(map_strstr * spm, mapkit_size_t used);

/* restore maxfill & minused and ensure these constraints */
	extern mapkit_error map_strstr_adjustcapacity(map_strstr * spm);

/* reallocate contents for "size" slots, keeping previous values */
	extern mapkit_error map_strstr_reallocate(map_strstr * spm, mapkit_size_t size);

/* return an adequate size for growing from "used" elements */
	extern mapkit_size_t map_strstr_growsize(map_strstr * spm, mapkit_size_t used);

/* return an adequate size for shrinking from "used" elements */
	extern mapkit_size_t map_strstr_shrinksize(map_strstr * spm, mapkit_size_t used);

/* return an adequate size for "used" elements */
	extern mapkit_size_t map_strstr_meansize(map_strstr * spm, mapkit_size_t used);

/* pointer to the value at key (insert as needed) */
	INLINE char **map_strstr_insertptr(map_strstr * spm, char *key);

/* pointer to the value at key or NULL */
	INLINE char **map_strstr_ptr(map_strstr * spm, char *key);

/* remove the key pointing to the value *ptr. No checking is done */
	INLINE mapkit_error map_strstr_removeptr(map_strstr * spm, char **ptr);

/* set the value at key (insert as needed) */
	INLINE mapkit_error map_strstr_set(map_strstr * spm, char *key, char *value);

/* return the value at key (key must exists if alwaysdefaultvalue is not set) */
	INLINE char *map_strstr_value(map_strstr * spm, char *key);

/* get the value at key (return an error code) */
	INLINE mapkit_error map_strstr_get(map_strstr * spm, char *key, char **value);

/* remove key (key must exists) */
	INLINE mapkit_error map_strstr_remove(map_strstr * spm, char *key);

/*
 * return the next index to a full slot of the map, or -1 if none
 * if index == -1, returns the first index, or -1 if none
 * the contents of the slot can be accessed by spm->contents[index].key and
 * spm->contents[index].value
 */
	extern mapkit_size_t map_strstr_next(map_strstr * spm, mapkit_size_t iindex);

/*
 * return the next ptr to a full slot of the map, or NULL if none
 * if pos_contents == NULL, returns the first full slot, or NULL if none
 * the contents of the slot can be accessed by pos_contents->key and
 * pos_contents->value
 */
	extern map_strstr_storage *map_strstr_nextptr(map_strstr * spm, map_strstr_storage * pos_contents);

/* allocate an array of map_strstr_element with all "count" (key,value) pairs */
	extern mapkit_error map_strstr_getall(map_strstr * spm, map_strstr_element ** array, mapkit_size_t * count);

/* remove all values equal to the defaultvalue */
	extern mapkit_error map_strstr_clean(map_strstr * spm);

/* compare the key of two map_strstr_element (for qsort) */
	extern int map_strstr_compare(const void *e1, const void *e2);

/* allocate an array of map_strstr_element with all "count" (key,value) pairs */
	extern mapkit_error map_strstr_getall_sorted(map_strstr * spm, map_strstr_element ** array, mapkit_size_t * count);

/* insert all count map_strstr_element into spm */
	extern mapkit_error map_strstr_setall(map_strstr * spm, map_strstr_element * array, mapkit_size_t count);

/* remove all "count" elements whose keys are in "array" */
	extern mapkit_error map_strstr_removeall(map_strstr * spm, char **array, mapkit_size_t count);

/* print statistics */
	extern void map_strstr_printstats(map_strstr * spm);

/* Private functions */

	extern mapkit_error map_strstr_remove_s(map_strstr * spm, char *key, mapkit_hash_t hash);
	extern char *map_strstr_value_s(map_strstr * spm, char *key, mapkit_hash_t hash);
	extern mapkit_error map_strstr_get_s(map_strstr * spm, char *key, char **value, mapkit_hash_t hash);
	extern mapkit_error map_strstr_set_s(map_strstr * spm, char *key, char *value, mapkit_hash_t hash);
	extern char **map_strstr_insertptr_s(map_strstr * spm, char *key, mapkit_hash_t hash);
	extern char **map_strstr_ptr_s(map_strstr * spm, char *key, mapkit_hash_t hash);

/* Implementation */

/* INLINEd functions */
	mapkit_error map_strstr_remove(map_strstr * spm, char *key) {
		mapkit_hash_t hash;
		map_strstr_storage *contents;

		contents = &(spm->contents[(hash = mapkit_strhash(key)) % spm->size]);

		if ((contents->state == MAPKIT_FULLSLOT)
		    && (strcmp(contents->key, key) == 0)) {
			contents->state = MAPKIT_DELETEDSLOT;
#ifdef MAPKIT_COLLISIONS
			spm->keyindexs++;
#endif
			spm->used--;
			if (spm->used < spm->minused) {
#ifdef MAPKIT_DEBUG
				fprintf(stderr, "MAPKIT: used < minused\n");
#endif
				return map_strstr_reallocate(spm, map_strstr_shrinksize(spm, spm->used));
			}
		} else
			return map_strstr_remove_s(spm, key, hash);
		return MAPKIT_OK;
	}

	mapkit_error map_strstr_removeptr(map_strstr * spm, char **ptr) {
		map_strstr_storage *sptr = (map_strstr_storage *) ((char *) ptr - offsetof(map_strstr_storage, value));

		sptr->state = MAPKIT_DELETEDSLOT;
		spm->used--;
		if (spm->used < spm->minused) {
#ifdef MAPKIT_DEBUG
			fprintf(stderr, "MAPKIT: used < minused\n");
#endif
			return map_strstr_reallocate(spm, map_strstr_shrinksize(spm, spm->used));
		}

		return MAPKIT_OK;
	}

	mapkit_error map_strstr_set(map_strstr * spm, char *key, char *value) {
		mapkit_hash_t hash;
		map_strstr_storage *contents;

		if (spm->alwaysdefault && (strcmp(value, spm->defaultvalue) == 0))
			return map_strstr_remove(spm, key);

		contents = &(spm->contents[(hash = mapkit_strhash(key)) % spm->size]);

		if (contents->state == MAPKIT_FREESLOT) {
			contents->state = MAPKIT_FULLSLOT;
			contents->key = key;
			contents->value = value;
#ifdef MAPKIT_COLLISIONS
			spm->insertionindexs++;
#endif
			spm->used++;
			spm->fill++;
			if (spm->fill > spm->maxfill) {
#ifdef MAPKIT_DEBUG
				fprintf(stderr, "MAPKIT: fill > maxfill\n");
#endif
				return map_strstr_reallocate(spm, map_strstr_growsize(spm, spm->used));
			}
			return MAPKIT_OK;
		} else if ((contents->state == MAPKIT_FULLSLOT)
			   && (strcmp(contents->key, key) == 0)) {
			contents->value = value;
#ifdef MAPKIT_COLLISIONS
			spm->insertionindexs++;
#endif
			return MAPKIT_OK;
		} else
			return map_strstr_set_s(spm, key, value, hash);
	}

	char *map_strstr_value(map_strstr * spm, char *key) {
		mapkit_hash_t hash;
		map_strstr_storage *contents;

		contents = &(spm->contents[(hash = mapkit_strhash(key)) % spm->size]);

		if ((contents->state == MAPKIT_FULLSLOT)
		    && (strcmp(contents->key, key) == 0)) {
#ifdef MAPKIT_COLLISIONS
			spm->keyindexs++;
#endif
			return contents->value;
		} else if (spm->alwaysdefault && (contents->state == MAPKIT_FREESLOT)) {
#ifdef MAPKIT_COLLISIONS
			spm->keyindexs++;
#endif
			return spm->defaultvalue;
		} else
			return map_strstr_value_s(spm, key, hash);
	}

	mapkit_error map_strstr_get(map_strstr * spm, char *key, char **value) {
		mapkit_hash_t hash;
		map_strstr_storage *contents;

		contents = &(spm->contents[(hash = mapkit_strhash(key)) % spm->size]);

		if ((contents->state == MAPKIT_FULLSLOT)
		    && (strcmp(contents->key, key) == 0)) {
#ifdef MAPKIT_COLLISIONS
			spm->keyindexs++;
#endif
			*value = contents->value;
		} else if (spm->alwaysdefault && (contents->state == MAPKIT_FREESLOT)) {
#ifdef MAPKIT_COLLISIONS
			spm->keyindexs++;
#endif
			*value = spm->defaultvalue;
		} else
			return map_strstr_get_s(spm, key, value, hash);

		return MAPKIT_OK;
	}

	char **map_strstr_insertptr(map_strstr * spm, char *key) {
		mapkit_hash_t hash;
		map_strstr_storage *contents;

		contents = &(spm->contents[(hash = mapkit_strhash(key)) % spm->size]);

		if (contents->state == MAPKIT_FREESLOT) {
			if (spm->fill >= spm->maxfill)
				return map_strstr_insertptr_s(spm, key, hash);

			contents->state = MAPKIT_FULLSLOT;
			contents->key = key;
			contents->value = spm->defaultvalue;
#ifdef MAPKIT_COLLISIONS
			spm->insertionindexs++;
#endif
			spm->used++;
			spm->fill++;
			return &(contents->value);
		} else if ((contents->state == MAPKIT_FULLSLOT)
			   && (strcmp(contents->key, key) == 0)) {
#ifdef MAPKIT_COLLISIONS
			spm->insertionindexs++;
#endif
			return &(contents->value);
		} else
			return map_strstr_insertptr_s(spm, key, hash);
	}

	char **map_strstr_ptr(map_strstr * spm, char *key) {
		mapkit_hash_t hash;
		map_strstr_storage *contents;

		contents = &(spm->contents[(hash = mapkit_strhash(key)) % spm->size]);

		if ((contents->state == MAPKIT_FULLSLOT)
		    && (strcmp(contents->key, key) == 0)) {
#ifdef MAPKIT_COLLISIONS
			spm->insertionindexs++;
#endif
			if ((!spm->alwaysdefault) || (!(strcmp(contents->value, spm->defaultvalue) == 0)))
				return &(contents->value);
			else
				return NULL;
		} else if (contents->state == MAPKIT_FREESLOT) {
#ifdef MAPKIT_COLLISIONS
			spm->insertionindexs++;
#endif
			return NULL;
		} else
			return map_strstr_ptr_s(spm, key, hash);
	}

/* Sparse structures */

/*
  Prototypes for spvector (int -> double)
  Default value : 0.0
  By default, unassigned values always default to the default value.
*/

	typedef struct {
		int key;
		double value;
	} spvector_storage;

	typedef spvector_storage spvector_element;

	typedef struct {
		mapkit_size_t size,			       /* size of the hash table */
		 fill,					       /* # of non-free slots in the hash table */
		 used,					       /* # of full slots */
		 maxfill,				       /* max # of non-free slots before rehash */
		 minused /* min # of full slots before rehash */ ;

		double maxfillfactor,			       /* maxfill = maxfillfactor * size */
		 minusedfactor;				       /* minused = minusedfactor * size */

#ifdef MAPKIT_COLLISIONS
		unsigned long insertionindexs, insertionindex_collisions;
		unsigned long keyindexs, keyindex_collisions;
#endif

		double defaultvalue;

		int alwaysdefault;

		spvector_storage *contents;
	} spvector;

/* Macros */

/* set to inform the type is available */
#define MAPKIT_spvector

/* read/write access to the value at key (insert as needed) */
#define spvector_(spm, key) (*(spvector_insertptr(spm, key)))

/* true if key exists */
#define spvector_haskey(spm, key) ((spvector_ptr(spm, key)) != NULL)

/* Prototypes */

/* initialize the structure */
	extern mapkit_error spvector_init(spvector * spm);

/* free contents */
	extern void spvector_free(spvector * spm);

/* copy the "from" map into a new "to" map */
	extern mapkit_error spvector_copy(spvector * to, spvector * from);

/* initialize the structure for "used" elements */
	extern mapkit_error spvector_init_hint(spvector * spm, mapkit_size_t used);

/* ensure at least "used" slots (including currently occupied slots) */
	extern mapkit_error spvector_ensurecapacity(spvector * spm, mapkit_size_t used);

/* restore maxfill & minused and ensure these constraints */
	extern mapkit_error spvector_adjustcapacity(spvector * spm);

/* reallocate contents for "size" slots, keeping previous values */
	extern mapkit_error spvector_reallocate(spvector * spm, mapkit_size_t size);

/* return an adequate size for growing from "used" elements */
	extern mapkit_size_t spvector_growsize(spvector * spm, mapkit_size_t used);

/* return an adequate size for shrinking from "used" elements */
	extern mapkit_size_t spvector_shrinksize(spvector * spm, mapkit_size_t used);

/* return an adequate size for "used" elements */
	extern mapkit_size_t spvector_meansize(spvector * spm, mapkit_size_t used);

/* pointer to the value at key (insert as needed) */
	INLINE double *spvector_insertptr(spvector * spm, int key);

/* pointer to the value at key or NULL */
	INLINE double *spvector_ptr(spvector * spm, int key);

/* remove the key pointing to the value *ptr. No checking is done */
	INLINE mapkit_error spvector_removeptr(spvector * spm, double *ptr);

/* set the value at key (insert as needed) */
	INLINE mapkit_error spvector_set(spvector * spm, int key, double value);

/* return the value at key (key must exists if alwaysdefaultvalue is not set) */
	INLINE double spvector_value(spvector * spm, int key);

/* get the value at key (return an error code) */
	INLINE mapkit_error spvector_get(spvector * spm, int key, double *value);

/* remove key (key must exists) */
	INLINE mapkit_error spvector_remove(spvector * spm, int key);

/*
 * return the next index to a full slot of the map, or -1 if none
 * if index == -1, returns the first index, or -1 if none
 * the contents of the slot can be accessed by spm->contents[index].key and
 * spm->contents[index].value
 */
	extern mapkit_size_t spvector_next(spvector * spm, mapkit_size_t iindex);

/*
 * return the next ptr to a full slot of the map, or NULL if none
 * if pos_contents == NULL, returns the first full slot, or NULL if none
 * the contents of the slot can be accessed by pos_contents->key and
 * pos_contents->value
 */
	extern spvector_storage *spvector_nextptr(spvector * spm, spvector_storage * pos_contents);

/* allocate an array of spvector_element with all "count" (key,value) pairs */
	extern mapkit_error spvector_getall(spvector * spm, spvector_element ** array, mapkit_size_t * count);

/* remove all values equal to the defaultvalue */
	extern mapkit_error spvector_clean(spvector * spm);

/* compare the key of two spvector_element (for qsort) */
	extern int spvector_compare(const void *e1, const void *e2);

/* allocate an array of spvector_element with all "count" (key,value) pairs */
	extern mapkit_error spvector_getall_sorted(spvector * spm, spvector_element ** array, mapkit_size_t * count);

/* insert all count spvector_element into spm */
	extern mapkit_error spvector_setall(spvector * spm, spvector_element * array, mapkit_size_t count);

/* remove all "count" elements whose keys are in "array" */
	extern mapkit_error spvector_removeall(spvector * spm, int *array, mapkit_size_t count);

/* print statistics */
	extern void spvector_printstats(spvector * spm);

/* Private functions */

	extern mapkit_error spvector_remove_s(spvector * spm, int key);
	extern double spvector_value_s(spvector * spm, int key);
	extern mapkit_error spvector_get_s(spvector * spm, int key, double *value);
	extern mapkit_error spvector_set_s(spvector * spm, int key, double value);
	extern double *spvector_insertptr_s(spvector * spm, int key);
	extern double *spvector_ptr_s(spvector * spm, int key);

/* Implementation */

/* INLINEd functions */
	mapkit_error spvector_remove(spvector * spm, int key) {
		spvector_storage *contents;

		/*
		 * key must be >= MAPKIT_FULLSLOT 
		 */
		if (key < MAPKIT_FULLSLOT)
			MAPKIT_ERROR(MAPKIT_EBADKEY);

		contents = &(spm->contents[((mapkit_hash_t) key) % spm->size]);

		if (contents->key == key) {
			contents->key = MAPKIT_DELETEDSLOT;
#ifdef MAPKIT_COLLISIONS
			spm->keyindexs++;
#endif
			spm->used--;
			if (spm->used < spm->minused) {
#ifdef MAPKIT_DEBUG
				fprintf(stderr, "MAPKIT: used < minused\n");
#endif
				return spvector_reallocate(spm, spvector_shrinksize(spm, spm->used));
			}
		} else
			return spvector_remove_s(spm, key);
		return MAPKIT_OK;
	}

	mapkit_error spvector_removeptr(spvector * spm, double *ptr) {
		spvector_storage *sptr = (spvector_storage *) ((char *) ptr - offsetof(spvector_storage, value));

		sptr->key = MAPKIT_DELETEDSLOT;
		spm->used--;
		if (spm->used < spm->minused) {
#ifdef MAPKIT_DEBUG
			fprintf(stderr, "MAPKIT: used < minused\n");
#endif
			return spvector_reallocate(spm, spvector_shrinksize(spm, spm->used));
		}

		return MAPKIT_OK;
	}

	mapkit_error spvector_set(spvector * spm, int key, double value) {
		int ckey;
		spvector_storage *contents;

		if (spm->alwaysdefault && ((value) == (spm->defaultvalue)))
			return spvector_remove(spm, key);

		/*
		 * key must be >= MAPKIT_FULLSLOT 
		 */
		if (key < MAPKIT_FULLSLOT)
			MAPKIT_ERROR(MAPKIT_EBADKEY);

		contents = &(spm->contents[((mapkit_hash_t) key) % spm->size]);

		if ((ckey = contents->key) == MAPKIT_FREESLOT) {
			contents->key = key;
			contents->value = value;
#ifdef MAPKIT_COLLISIONS
			spm->insertionindexs++;
#endif
			spm->used++;
			spm->fill++;
			if (spm->fill > spm->maxfill) {
#ifdef MAPKIT_DEBUG
				fprintf(stderr, "MAPKIT: fill > maxfill\n");
#endif
				return spvector_reallocate(spm, spvector_growsize(spm, spm->used));
			}
			return MAPKIT_OK;
		} else if (ckey == key) {
			contents->value = value;
#ifdef MAPKIT_COLLISIONS
			spm->insertionindexs++;
#endif
			return MAPKIT_OK;
		} else
			return spvector_set_s(spm, key, value);
	}

	double spvector_value(spvector * spm, int key) {
		spvector_storage *contents;
		int ckey;

		/*
		 * key must be >= MAPKIT_FULLSLOT 
		 */
		if (key < MAPKIT_FULLSLOT)
			MAPKIT_FATAL_ERROR(MAPKIT_EBADKEY);

		contents = &(spm->contents[((mapkit_hash_t) key) % spm->size]);

		if ((ckey = contents->key) == key) {
#ifdef MAPKIT_COLLISIONS
			spm->keyindexs++;
#endif
			return contents->value;
		} else if (spm->alwaysdefault && (ckey == MAPKIT_FREESLOT)) {
#ifdef MAPKIT_COLLISIONS
			spm->keyindexs++;
#endif
			return spm->defaultvalue;
		} else
			return spvector_value_s(spm, key);
	}

	mapkit_error spvector_get(spvector * spm, int key, double *value) {
		spvector_storage *contents;
		int ckey;

		/*
		 * key must be >= MAPKIT_FULLSLOT 
		 */
		if (key < MAPKIT_FULLSLOT)
			MAPKIT_ERROR(MAPKIT_EBADKEY);

		contents = &(spm->contents[((mapkit_hash_t) key) % spm->size]);

		if ((ckey = contents->key) == key) {
#ifdef MAPKIT_COLLISIONS
			spm->keyindexs++;
#endif
			*value = contents->value;
		} else if (spm->alwaysdefault && (ckey == MAPKIT_FREESLOT)) {
#ifdef MAPKIT_COLLISIONS
			spm->keyindexs++;
#endif
			*value = spm->defaultvalue;
		} else
			return spvector_get_s(spm, key, value);

		return MAPKIT_OK;
	}

	double *spvector_insertptr(spvector * spm, int key) {
		int ckey;
		spvector_storage *contents;

		/*
		 * key must be >= MAPKIT_FULLSLOT 
		 */
		if (key < MAPKIT_FULLSLOT) {
			MAPKIT_ERROR_NORET(MAPKIT_EBADKEY);
			return NULL;
		}

		contents = &(spm->contents[((mapkit_hash_t) key) % spm->size]);

		if ((ckey = contents->key) == MAPKIT_FREESLOT) {
			if (spm->fill >= spm->maxfill)
				return spvector_insertptr_s(spm, key);

			contents->key = key;
			contents->value = spm->defaultvalue;
#ifdef MAPKIT_COLLISIONS
			spm->insertionindexs++;
#endif
			spm->used++;
			spm->fill++;
			return &(contents->value);
		} else if (ckey == key) {
#ifdef MAPKIT_COLLISIONS
			spm->insertionindexs++;
#endif
			return &(contents->value);
		} else
			return spvector_insertptr_s(spm, key);
	}

	double *spvector_ptr(spvector * spm, int key) {
		int ckey;
		spvector_storage *contents;

		/*
		 * key must be >= MAPKIT_FULLSLOT 
		 */
		if (key < MAPKIT_FULLSLOT) {
			MAPKIT_ERROR_NORET(MAPKIT_EBADKEY);
			return NULL;
		}

		contents = &(spm->contents[((mapkit_hash_t) key) % spm->size]);

		if ((ckey = contents->key) == key) {
#ifdef MAPKIT_COLLISIONS
			spm->insertionindexs++;
#endif
			if ((!spm->alwaysdefault) || (!((contents->value) == (spm->defaultvalue))))
				return &(contents->value);
			else
				return NULL;
		} else if (ckey == MAPKIT_FREESLOT) {
#ifdef MAPKIT_COLLISIONS
			spm->insertionindexs++;
#endif
			return NULL;
		} else
			return spvector_ptr_s(spm, key);
	}

	typedef struct spmatrix_key_pair_struct {
		int key1, key2;
	} spmatrix_key_pair;

/*
  Prototypes for _spmatrix (spmatrix_key_pair -> double)
  Default value : 0.0
  Uses a state field.
  By default, unassigned values always default to the default value.
*/

	typedef struct {
		signed char state;
		spmatrix_key_pair key;
		double value;
	} _spmatrix_storage;

	typedef struct {
		spmatrix_key_pair key;
		double value;
	} _spmatrix_element;

	typedef struct {
		mapkit_size_t size,			       /* size of the hash table */
		 fill,					       /* # of non-free slots in the hash table */
		 used,					       /* # of full slots */
		 maxfill,				       /* max # of non-free slots before rehash */
		 minused /* min # of full slots before rehash */ ;

		double maxfillfactor,			       /* maxfill = maxfillfactor * size */
		 minusedfactor;				       /* minused = minusedfactor * size */

#ifdef MAPKIT_COLLISIONS
		unsigned long insertionindexs, insertionindex_collisions;
		unsigned long keyindexs, keyindex_collisions;
#endif

		double defaultvalue;

		int alwaysdefault;

		_spmatrix_storage *contents;
	} _spmatrix;

/* Macros */

/* set to inform the type is available */
#define MAPKIT__spmatrix

/* read/write access to the value at key (insert as needed) */
#define _spmatrix_(spm, key) (*(_spmatrix_insertptr(spm, key)))

/* true if key exists */
#define _spmatrix_haskey(spm, key) ((_spmatrix_ptr(spm, key)) != NULL)

/* Prototypes */

/* initialize the structure */
	extern mapkit_error _spmatrix_init(_spmatrix * spm);

/* free contents */
	extern void _spmatrix_free(_spmatrix * spm);

/* copy the "from" map into a new "to" map */
	extern mapkit_error _spmatrix_copy(_spmatrix * to, _spmatrix * from);

/* initialize the structure for "used" elements */
	extern mapkit_error _spmatrix_init_hint(_spmatrix * spm, mapkit_size_t used);

/* ensure at least "used" slots (including currently occupied slots) */
	extern mapkit_error _spmatrix_ensurecapacity(_spmatrix * spm, mapkit_size_t used);

/* restore maxfill & minused and ensure these constraints */
	extern mapkit_error _spmatrix_adjustcapacity(_spmatrix * spm);

/* reallocate contents for "size" slots, keeping previous values */
	extern mapkit_error _spmatrix_reallocate(_spmatrix * spm, mapkit_size_t size);

/* return an adequate size for growing from "used" elements */
	extern mapkit_size_t _spmatrix_growsize(_spmatrix * spm, mapkit_size_t used);

/* return an adequate size for shrinking from "used" elements */
	extern mapkit_size_t _spmatrix_shrinksize(_spmatrix * spm, mapkit_size_t used);

/* return an adequate size for "used" elements */
	extern mapkit_size_t _spmatrix_meansize(_spmatrix * spm, mapkit_size_t used);

/* pointer to the value at key (insert as needed) */
	INLINE double *_spmatrix_insertptr(_spmatrix * spm, spmatrix_key_pair key);

/* pointer to the value at key or NULL */
	INLINE double *_spmatrix_ptr(_spmatrix * spm, spmatrix_key_pair key);

/* remove the key pointing to the value *ptr. No checking is done */
	INLINE mapkit_error _spmatrix_removeptr(_spmatrix * spm, double *ptr);

/* set the value at key (insert as needed) */
	INLINE mapkit_error _spmatrix_set(_spmatrix * spm, spmatrix_key_pair key, double value);

/* return the value at key (key must exists if alwaysdefaultvalue is not set) */
	INLINE double _spmatrix_value(_spmatrix * spm, spmatrix_key_pair key);

/* get the value at key (return an error code) */
	INLINE mapkit_error _spmatrix_get(_spmatrix * spm, spmatrix_key_pair key, double *value);

/* remove key (key must exists) */
	INLINE mapkit_error _spmatrix_remove(_spmatrix * spm, spmatrix_key_pair key);

/*
 * return the next index to a full slot of the map, or -1 if none
 * if index == -1, returns the first index, or -1 if none
 * the contents of the slot can be accessed by spm->contents[index].key and
 * spm->contents[index].value
 */
	extern mapkit_size_t _spmatrix_next(_spmatrix * spm, mapkit_size_t iindex);

/*
 * return the next ptr to a full slot of the map, or NULL if none
 * if pos_contents == NULL, returns the first full slot, or NULL if none
 * the contents of the slot can be accessed by pos_contents->key and
 * pos_contents->value
 */
	extern _spmatrix_storage *_spmatrix_nextptr(_spmatrix * spm, _spmatrix_storage * pos_contents);

/* allocate an array of _spmatrix_element with all "count" (key,value) pairs */
	extern mapkit_error _spmatrix_getall(_spmatrix * spm, _spmatrix_element ** array, mapkit_size_t * count);

/* remove all values equal to the defaultvalue */
	extern mapkit_error _spmatrix_clean(_spmatrix * spm);

/* compare the key of two _spmatrix_element (for qsort) */
	extern int _spmatrix_compare(const void *e1, const void *e2);

/* allocate an array of _spmatrix_element with all "count" (key,value) pairs */
	extern mapkit_error _spmatrix_getall_sorted(_spmatrix * spm, _spmatrix_element ** array, mapkit_size_t * count);

/* insert all count _spmatrix_element into spm */
	extern mapkit_error _spmatrix_setall(_spmatrix * spm, _spmatrix_element * array, mapkit_size_t count);

/* remove all "count" elements whose keys are in "array" */
	extern mapkit_error _spmatrix_removeall(_spmatrix * spm, spmatrix_key_pair * array, mapkit_size_t count);

/* print statistics */
	extern void _spmatrix_printstats(_spmatrix * spm);

/* Private functions */

	extern mapkit_error _spmatrix_remove_s(_spmatrix * spm, spmatrix_key_pair key, mapkit_hash_t hash);
	extern double _spmatrix_value_s(_spmatrix * spm, spmatrix_key_pair key, mapkit_hash_t hash);
	extern mapkit_error _spmatrix_get_s(_spmatrix * spm, spmatrix_key_pair key, double *value, mapkit_hash_t hash);
	extern mapkit_error _spmatrix_set_s(_spmatrix * spm, spmatrix_key_pair key, double value, mapkit_hash_t hash);
	extern double *_spmatrix_insertptr_s(_spmatrix * spm, spmatrix_key_pair key, mapkit_hash_t hash);
	extern double *_spmatrix_ptr_s(_spmatrix * spm, spmatrix_key_pair key, mapkit_hash_t hash);

/* Implementation */

/* INLINEd functions */
	mapkit_error _spmatrix_remove(_spmatrix * spm, spmatrix_key_pair key) {
		mapkit_hash_t hash;
		_spmatrix_storage *contents;

		contents = &(spm->contents[(hash = (mapkit_hash_t) ((key).key1 + (key).key2 * 65519)) % spm->size]);

		if ((contents->state == MAPKIT_FULLSLOT)
		    && (((contents->key).key1 == (key).key1) && ((contents->key).key2 == (key).key2))) {
			contents->state = MAPKIT_DELETEDSLOT;
#ifdef MAPKIT_COLLISIONS
			spm->keyindexs++;
#endif
			spm->used--;
			if (spm->used < spm->minused) {
#ifdef MAPKIT_DEBUG
				fprintf(stderr, "MAPKIT: used < minused\n");
#endif
				return _spmatrix_reallocate(spm, _spmatrix_shrinksize(spm, spm->used));
			}
		} else
			return _spmatrix_remove_s(spm, key, hash);
		return MAPKIT_OK;
	}

	mapkit_error _spmatrix_removeptr(_spmatrix * spm, double *ptr) {
		_spmatrix_storage *sptr = (_spmatrix_storage *) ((char *) ptr - offsetof(_spmatrix_storage, value));

		sptr->state = MAPKIT_DELETEDSLOT;
		spm->used--;
		if (spm->used < spm->minused) {
#ifdef MAPKIT_DEBUG
			fprintf(stderr, "MAPKIT: used < minused\n");
#endif
			return _spmatrix_reallocate(spm, _spmatrix_shrinksize(spm, spm->used));
		}

		return MAPKIT_OK;
	}

	mapkit_error _spmatrix_set(_spmatrix * spm, spmatrix_key_pair key, double value) {
		mapkit_hash_t hash;
		_spmatrix_storage *contents;

		if (spm->alwaysdefault && ((value) == (spm->defaultvalue)))
			return _spmatrix_remove(spm, key);

		contents = &(spm->contents[(hash = (mapkit_hash_t) ((key).key1 + (key).key2 * 65519)) % spm->size]);

		if (contents->state == MAPKIT_FREESLOT) {
			contents->state = MAPKIT_FULLSLOT;
			contents->key = key;
			contents->value = value;
#ifdef MAPKIT_COLLISIONS
			spm->insertionindexs++;
#endif
			spm->used++;
			spm->fill++;
			if (spm->fill > spm->maxfill) {
#ifdef MAPKIT_DEBUG
				fprintf(stderr, "MAPKIT: fill > maxfill\n");
#endif
				return _spmatrix_reallocate(spm, _spmatrix_growsize(spm, spm->used));
			}
			return MAPKIT_OK;
		} else if ((contents->state == MAPKIT_FULLSLOT)
			   && (((contents->key).key1 == (key).key1) && ((contents->key).key2 == (key).key2))) {
			contents->value = value;
#ifdef MAPKIT_COLLISIONS
			spm->insertionindexs++;
#endif
			return MAPKIT_OK;
		} else
			return _spmatrix_set_s(spm, key, value, hash);
	}

	double _spmatrix_value(_spmatrix * spm, spmatrix_key_pair key) {
		mapkit_hash_t hash;
		_spmatrix_storage *contents;

		contents = &(spm->contents[(hash = (mapkit_hash_t) ((key).key1 + (key).key2 * 65519)) % spm->size]);

		if ((contents->state == MAPKIT_FULLSLOT)
		    && (((contents->key).key1 == (key).key1) && ((contents->key).key2 == (key).key2))) {
#ifdef MAPKIT_COLLISIONS
			spm->keyindexs++;
#endif
			return contents->value;
		} else if (spm->alwaysdefault && (contents->state == MAPKIT_FREESLOT)) {
#ifdef MAPKIT_COLLISIONS
			spm->keyindexs++;
#endif
			return spm->defaultvalue;
		} else
			return _spmatrix_value_s(spm, key, hash);
	}

	mapkit_error _spmatrix_get(_spmatrix * spm, spmatrix_key_pair key, double *value) {
		mapkit_hash_t hash;
		_spmatrix_storage *contents;

		contents = &(spm->contents[(hash = (mapkit_hash_t) ((key).key1 + (key).key2 * 65519)) % spm->size]);

		if ((contents->state == MAPKIT_FULLSLOT)
		    && (((contents->key).key1 == (key).key1) && ((contents->key).key2 == (key).key2))) {
#ifdef MAPKIT_COLLISIONS
			spm->keyindexs++;
#endif
			*value = contents->value;
		} else if (spm->alwaysdefault && (contents->state == MAPKIT_FREESLOT)) {
#ifdef MAPKIT_COLLISIONS
			spm->keyindexs++;
#endif
			*value = spm->defaultvalue;
		} else
			return _spmatrix_get_s(spm, key, value, hash);

		return MAPKIT_OK;
	}

	double *_spmatrix_insertptr(_spmatrix * spm, spmatrix_key_pair key) {
		mapkit_hash_t hash;
		_spmatrix_storage *contents;

		contents = &(spm->contents[(hash = (mapkit_hash_t) ((key).key1 + (key).key2 * 65519)) % spm->size]);

		if (contents->state == MAPKIT_FREESLOT) {
			if (spm->fill >= spm->maxfill)
				return _spmatrix_insertptr_s(spm, key, hash);

			contents->state = MAPKIT_FULLSLOT;
			contents->key = key;
			contents->value = spm->defaultvalue;
#ifdef MAPKIT_COLLISIONS
			spm->insertionindexs++;
#endif
			spm->used++;
			spm->fill++;
			return &(contents->value);
		} else if ((contents->state == MAPKIT_FULLSLOT)
			   && (((contents->key).key1 == (key).key1) && ((contents->key).key2 == (key).key2))) {
#ifdef MAPKIT_COLLISIONS
			spm->insertionindexs++;
#endif
			return &(contents->value);
		} else
			return _spmatrix_insertptr_s(spm, key, hash);
	}

	double *_spmatrix_ptr(_spmatrix * spm, spmatrix_key_pair key) {
		mapkit_hash_t hash;
		_spmatrix_storage *contents;

		contents = &(spm->contents[(hash = (mapkit_hash_t) ((key).key1 + (key).key2 * 65519)) % spm->size]);

		if ((contents->state == MAPKIT_FULLSLOT)
		    && (((contents->key).key1 == (key).key1) && ((contents->key).key2 == (key).key2))) {
#ifdef MAPKIT_COLLISIONS
			spm->insertionindexs++;
#endif
			if ((!spm->alwaysdefault) || (!((contents->value) == (spm->defaultvalue))))
				return &(contents->value);
			else
				return NULL;
		} else if (contents->state == MAPKIT_FREESLOT) {
#ifdef MAPKIT_COLLISIONS
			spm->insertionindexs++;
#endif
			return NULL;
		} else
			return _spmatrix_ptr_s(spm, key, hash);
	}

/* Aliases */
	typedef _spmatrix spmatrix;
	typedef _spmatrix_storage spmatrix_storage;

	typedef struct {
		int key1;
		int key2;
		double value;
	} spmatrix_element;

	typedef struct {
		int key1;
		int key2;
	} spmatrix_key;

/* Macros */

/* set to inform the type is available */
#define MAPKIT_spmatrix

/* read/write access to the value at (key1, key2) (insert as needed) */
#define spmatrix_(spm, key1, key2) (*(spmatrix_insertptr(spm, key1, key2)))

/* true if (key1, key2) exists */
#define spmatrix_haskey(spm, key1, key2) ((spmatrix_ptr(spm, key1, key2)) != NULL)

/* Basic functions */

/* set default values */
	extern mapkit_error spmatrix_init(spmatrix * spm);

/* free contents */
	extern void spmatrix_free(spmatrix * spm);

/* copy the "from" map into a new "to" map */
	extern mapkit_error spmatrix_copy(spmatrix * to, spmatrix * from);

/* Advanced functions */

/* initialize the structure for "used" elements */
	extern mapkit_error spmatrix_init_hint(spmatrix * spm, mapkit_size_t used);

/* ensure at least "used" slots (including currently occupied slots) */
	extern mapkit_error spmatrix_ensurecapacity(spmatrix * spm, mapkit_size_t used);

/* restore maxfill & minused and ensure these constraints */
	extern mapkit_error spmatrix_adjustcapacity(spmatrix * spm);

/* Matrix functions */

/* return a pointer to the value at (key1,key2) (insert as needed) */
	INLINE double *spmatrix_insertptr(spmatrix * spm, int key1, int key2);

/* return a pointer to thevalue at (key1,key2) or NULL */
	INLINE double *spmatrix_ptr(spmatrix * spm, int key1, int key2);

/* remove the key pointing to the value *ptr. No checking is done */
	INLINE mapkit_error spmatrix_removeptr(spmatrix * spm, double *ptr);

/* return the value at (key1,key2) */
	INLINE double spmatrix_value(spmatrix * spm, int key1, int key2);

/* get the value at (key1,key2), return an error code */
	INLINE mapkit_error spmatrix_get(spmatrix * spm, int key1, int key2, double *value);

/* set the value at (key1,key2) */
	INLINE mapkit_error spmatrix_set(spmatrix * spm, int key1, int key2, double value);

/* remove the value at (key1,key2) */
	INLINE mapkit_error spmatrix_remove(spmatrix * spm, int key1, int key2);

/*
 * return the next index to a full slot of the map, or -1 if none
 * if index == -1, returns the first index, or -1 if none
 * the contents of the slot can be accessed by spm->contents[index].key and
 * spm->contents[index].value
 * the key, of type spmatrix_key_pair contains the key pair.
 */
	INLINE mapkit_size_t spmatrix_next(spmatrix * spm, mapkit_size_t iindex);

/*
 * return the next ptr to a full slot of the map, or NULL if none
 * If pos_contents == NULL, returns the first full slot, or NULL if none
 * The contents of the slot can be accessed by pos_contents->key and
 * pos_contents->value
 * the key, of type spmatrix_key_pair contains the key pair.
 */
	INLINE spmatrix_storage *spmatrix_nextptr(spmatrix * spm, spmatrix_storage * pos_contents);

/* allocate an array of spmatrix_element with all
  (key1, key2, value) triplets */
	extern mapkit_error spmatrix_getall(spmatrix * spm, spmatrix_element ** array, mapkit_size_t * count);

/* remove all values equal to the defaultvalue */
	extern mapkit_error spmatrix_clean(spmatrix * spm);

/* allocate an array of spmatrix_element with all
  (key1, key2, value) triplets, sorted by key1, key2 */
	extern mapkit_error spmatrix_getall_sorted1(spmatrix * spm, spmatrix_element ** array, mapkit_size_t * count);

/* allocate an array au spmatrix_element with all
  (key1, key2, value) triplets, sorted by key2, key1 */
	extern mapkit_error spmatrix_getall_sorted2(spmatrix * spm, spmatrix_element ** array, mapkit_size_t * count);

/* insert all count spmatrix_element into spm */
	extern mapkit_error spmatrix_setall(spmatrix * spm, spmatrix_element * array, mapkit_size_t count);

/* remove all "count" elements whose keys are in "array" */
	extern mapkit_error spmatrix_removeall(spmatrix * spm, spmatrix_key * array, mapkit_size_t count);

/* print statistics */
	extern void spmatrix_printstats(spmatrix * spm);

/* Implementation */

/* Inlined functions */

	double spmatrix_value(spmatrix * spm, int key1, int key2) {
		spmatrix_key_pair key;

		(((key).key1 = key1), ((key).key2 = key2));
		return _spmatrix_value(spm, key);
	}

	mapkit_error spmatrix_get(spmatrix * spm, int key1, int key2, double *value) {
		spmatrix_key_pair key;

		(((key).key1 = key1), ((key).key2 = key2));
		return _spmatrix_get(spm, key, value);
	}

	mapkit_error spmatrix_set(spmatrix * spm, int key1, int key2, double value) {
		spmatrix_key_pair key;

		(((key).key1 = key1), ((key).key2 = key2));
		return _spmatrix_set(spm, key, value);
	}

	mapkit_error spmatrix_remove(spmatrix * spm, int key1, int key2) {
		spmatrix_key_pair key;

		(((key).key1 = key1), ((key).key2 = key2));
		return _spmatrix_remove(spm, key);
	}

	double *spmatrix_insertptr(spmatrix * spm, int key1, int key2) {
		spmatrix_key_pair key;

		(((key).key1 = key1), ((key).key2 = key2));
		return _spmatrix_insertptr(spm, key);
	}

	double *spmatrix_ptr(spmatrix * spm, int key1, int key2) {
		spmatrix_key_pair key;

		(((key).key1 = key1), ((key).key2 = key2));
		return _spmatrix_ptr(spm, key);
	}

	mapkit_error spmatrix_removeptr(spmatrix * spm, double *ptr) {
		return _spmatrix_removeptr(spm, ptr);
	}

	mapkit_size_t spmatrix_next(spmatrix * spm, mapkit_size_t iindex) {
		return _spmatrix_next(spm, iindex);
	}

	spmatrix_storage *spmatrix_nextptr(spmatrix * spm, spmatrix_storage * pos_contents) {
		return _spmatrix_nextptr(spm, pos_contents);
	}

#if defined(__cplusplus)
}
#endif

#endif							       /* _MAPKIT_ */
