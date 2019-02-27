
/* hash.c
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
 */

/*!
  \file hash.c
  \brief The hash-library used in GMRFLib, see \ref hashP.h if you want to use it.
*/

#ifndef HGVERSION
#define HGVERSION
#endif
//static const char RCSId[] = "file: " __FILE__ "  " HGVERSION;

/* Pre-hg-Id: $Id: hash.c,v 1.13 2009/06/03 09:04:54 hrue Exp $ */

#include <string.h>
#if !defined(__FreeBSD__)
#include <malloc.h>
#endif
#include <stdlib.h>
#include <stdio.h>
#include <limits.h>

#include "GMRFLib/GMRFLib.h"
#include "GMRFLib/GMRFLibP.h"
#include "GMRFLib/hashP.h"

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
  THIS FILE CONTAINS the source mapkit.c and mapkit_generic.c

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

*/

/*  */

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

/* static char const rcsid[] =
   "@(#) $Jeannot: mapkit_generic.c,v 1.16 2004/03/05 16:06:06 js Exp $"; */

const char *mapkit_error_msg[5] = {
	"ok",
	"fill >= max. admissible size",
	"malloc failed",
	"bad key",
	"key not found"
};

/* Static constants */

#if LONG_MAX <= 0x7fffffff

/* 32 bits longs */
#define MAPKIT_PRIMES 343

static mapkit_size_t mapkit_primes[MAPKIT_PRIMES] = { 5L, 7L, 13L, 19L, 31L, 43L, 61L, 73L, 103L, 109L, 139L, 151L,
	181L, 193L, 199L, 229L, 241L, 271L, 283L, 313L, 349L, 421L, 433L,
	463L, 523L, 571L, 601L, 619L, 643L, 661L, 811L, 829L, 859L, 883L,
	1021L, 1063L, 1093L, 1153L, 1231L, 1291L, 1321L, 1429L, 1489L, 1609L,
	1669L, 1723L, 1789L, 1873L, 1951L, 2029L, 2113L, 2143L, 2239L, 2341L,
	2383L, 2551L, 2659L, 2791L, 2803L, 2971L, 3001L, 3121L, 3259L, 3391L,
	3559L, 3673L, 3853L, 4021L, 4219L, 4423L, 4639L, 4801L, 5023L, 5233L,
	5479L, 5743L, 5881L, 6133L, 6361L, 6661L, 6961L, 7309L, 7591L, 7951L,
	8293L, 8629L, 9043L, 9463L, 9931L, 10333L, 10711L, 11173L, 11719L,
	12253L, 12823L, 13399L, 14011L, 14629L, 15331L, 16069L, 16831L, 17659L,
	18541L, 19429L, 20359L, 21319L, 22369L, 23371L, 24421L, 25603L, 26881L,
	28183L, 29569L, 30871L, 32413L, 34033L, 35731L, 37363L, 39229L, 41179L,
	43051L, 45181L, 47419L, 49789L, 52183L, 54631L, 57349L, 60169L, 63031L,
	66109L, 69403L, 72871L, 76423L, 80233L, 84223L, 88339L, 92683L, 97303L,
	102103L, 107101L, 112363L, 117979L, 123733L, 129919L, 136399L, 143113L,
	150223L, 157669L, 165313L, 173431L, 182101L, 191143L, 200383L, 210361L,
	220879L, 231841L, 243433L, 255589L, 267961L, 281251L, 295201L, 309931L,
	325309L, 341461L, 358429L, 376099L, 394819L, 414433L, 435109L, 456811L,
	479431L, 503383L, 528511L, 554893L, 582511L, 611551L, 642079L, 674161L,
	707671L, 742993L, 780049L, 819031L, 859801L, 902599L, 947413L, 994711L,
	1044289L, 1096351L, 1150873L, 1208299L, 1268623L, 1331989L, 1398559L,
	1468459L, 1541821L, 1618831L, 1699741L, 1784581L, 1873771L, 1967419L,
	2065729L, 2168989L, 2277001L, 2390473L, 2509963L, 2635099L, 2766793L,
	2905033L, 3050263L, 3202483L, 3362329L, 3530071L, 3706489L, 3891781L,
	4086289L, 4290469L, 4504963L, 4730179L, 4966531L, 5214823L, 5475523L,
	5748961L, 6036271L, 6338029L, 6654649L, 6987109L, 7336423L, 7703239L,
	8088331L, 8492593L, 8916913L, 9362671L, 9830803L, 10322161L, 10838119L,
	11379961L, 11948791L, 12546013L, 13173031L, 13831669L, 14523193L,
	15248671L, 16011031L, 16811131L, 17651659L, 18534169L, 19460659L,
	20433691L, 21455281L, 22527469L, 23653813L, 24836311L, 26077591L,
	27381469L, 28750261L, 30187609L, 31696939L, 33281659L, 34945681L,
	36692881L, 38527303L, 40453549L, 42476143L, 44599741L, 46829149L,
	49170421L, 51628693L, 54210073L, 56920441L, 59765779L, 62753881L,
	65890753L, 69184879L, 72643801L, 76275961L, 80089663L, 84094039L,
	88298521L, 92713099L, 97348243L, 102215539L, 107326033L, 112691731L,
	118326181L, 124242331L, 130454101L, 136976053L, 143824789L, 151015939L,
	158566651L, 166494901L, 174819583L, 183560551L, 192738463L, 202375339L,
	212493961L, 223118629L, 234274321L, 245987983L, 258287221L, 271201321L,
	284761159L, 298999201L, 313948993L, 329646313L, 346128283L, 363434053L,
	381605641L, 400685893L, 420720019L, 441755893L, 463843339L, 487035169L,
	511386781L, 536955961L, 563803063L, 591992419L, 621591913L, 652670593L,
	685303783L, 719567743L, 755545711L, 793322989L, 832989109L, 874638211L,
	918369919L, 964287913L, 1012502221L, 1063127293L, 1116283603L,
	1172097763L, 1230702463L, 1292237251L, 1356848263L, 1424690611L,
	1495924909L, 1570720981L, 1649256319L, 1731718699L, 1818304471L,
	1909219603L, 2004680263L, 2147483647L
};
#else

/* 64 bits longs */
#define MAPKIT_PRIMES 798

static mapkit_size_t mapkit_primes[MAPKIT_PRIMES] = { 5L, 7L, 13L, 19L, 31L, 43L, 61L, 73L, 103L, 109L, 139L, 151L,
	181L, 193L, 199L, 229L, 241L, 271L, 283L, 313L, 349L, 421L, 433L,
	463L, 523L, 571L, 601L, 619L, 643L, 661L, 811L, 829L, 859L, 883L,
	1021L, 1063L, 1093L, 1153L, 1231L, 1291L, 1321L, 1429L, 1489L, 1609L,
	1669L, 1723L, 1789L, 1873L, 1951L, 2029L, 2113L, 2143L, 2239L, 2341L,
	2383L, 2551L, 2659L, 2791L, 2803L, 2971L, 3001L, 3121L, 3259L, 3391L,
	3559L, 3673L, 3853L, 4021L, 4219L, 4423L, 4639L, 4801L, 5023L, 5233L,
	5479L, 5743L, 5881L, 6133L, 6361L, 6661L, 6961L, 7309L, 7591L, 7951L,
	8293L, 8629L, 9043L, 9463L, 9931L, 10333L, 10711L, 11173L, 11719L,
	12253L, 12823L, 13399L, 14011L, 14629L, 15331L, 16069L, 16831L, 17659L,
	18541L, 19429L, 20359L, 21319L, 22369L, 23371L, 24421L, 25603L, 26881L,
	28183L, 29569L, 30871L, 32413L, 34033L, 35731L, 37363L, 39229L, 41179L,
	43051L, 45181L, 47419L, 49789L, 52183L, 54631L, 57349L, 60169L, 63031L,
	66109L, 69403L, 72871L, 76423L, 80233L, 84223L, 88339L, 92683L, 97303L,
	102103L, 107101L, 112363L, 117979L, 123733L, 129919L, 136399L, 143113L,
	150223L, 157669L, 165313L, 173431L, 182101L, 191143L, 200383L, 210361L,
	220879L, 231841L, 243433L, 255589L, 267961L, 281251L, 295201L, 309931L,
	325309L, 341461L, 358429L, 376099L, 394819L, 414433L, 435109L, 456811L,
	479431L, 503383L, 528511L, 554893L, 582511L, 611551L, 642079L, 674161L,
	707671L, 742993L, 780049L, 819031L, 859801L, 902599L, 947413L, 994711L,
	1044289L, 1096351L, 1150873L, 1208299L, 1268623L, 1331989L, 1398559L,
	1468459L, 1541821L, 1618831L, 1699741L, 1784581L, 1873771L, 1967419L,
	2065729L, 2168989L, 2277001L, 2390473L, 2509963L, 2635099L, 2766793L,
	2905033L, 3050263L, 3202483L, 3362329L, 3530071L, 3706489L, 3891781L,
	4086289L, 4290469L, 4504963L, 4730179L, 4966531L, 5214823L, 5475523L,
	5748961L, 6036271L, 6338029L, 6654649L, 6987109L, 7336423L, 7703239L,
	8088331L, 8492593L, 8916913L, 9362671L, 9830803L, 10322161L, 10838119L,
	11379961L, 11948791L, 12546013L, 13173031L, 13831669L, 14523193L,
	15248671L, 16011031L, 16811131L, 17651659L, 18534169L, 19460659L,
	20433691L, 21455281L, 22527469L, 23653813L, 24836311L, 26077591L,
	27381469L, 28750261L, 30187609L, 31696939L, 33281659L, 34945681L,
	36692881L, 38527303L, 40453549L, 42476143L, 44599741L, 46829149L,
	49170421L, 51628693L, 54210073L, 56920441L, 59765779L, 62753881L,
	65890753L, 69184879L, 72643801L, 76275961L, 80089663L, 84094039L,
	88298521L, 92713099L, 97348243L, 102215539L, 107326033L, 112691731L,
	118326181L, 124242331L, 130454101L, 136976053L, 143824789L, 151015939L,
	158566651L, 166494901L, 174819583L, 183560551L, 192738463L, 202375339L,
	212493961L, 223118629L, 234274321L, 245987983L, 258287221L, 271201321L,
	284761159L, 298999201L, 313948993L, 329646313L, 346128283L, 363434053L,
	381605641L, 400685893L, 420720019L, 441755893L, 463843339L, 487035169L,
	511386781L, 536955961L, 563803063L, 591992419L, 621591913L, 652670593L,
	685303783L, 719567743L, 755545711L, 793322989L, 832989109L, 874638211L,
	918369919L, 964287913L, 1012502221L, 1063127293L, 1116283603L,
	1172097763L, 1230702463L, 1292237251L, 1356848263L, 1424690611L,
	1495924909L, 1570720981L, 1649256319L, 1731718699L, 1818304471L,
	1909219603L, 2004680263L, 2104914271L, 2210159113L, 2320667053L,
	2436700369L, 2558534863L, 2686461301L, 2820784321L, 2961823519L,
	3109914463L, 3265409671L, 3428679931L, 3600113089L, 3780118621L,
	3969124501L, 4167580723L, 4375959571L, 4594757509L, 4824495043L,
	5065719649L, 5319004963L, 5584954621L, 5864201491L, 6157411219L,
	6465281743L, 6788545669L, 7127972773L, 7484371333L, 7858589383L,
	8251518469L, 8664094291L, 9097298143L, 9552162799L, 10029770899L,
	10531258921L, 11057821561L, 11610712021L, 12191247619L, 12800809999L,
	13440850483L, 14112891943L, 14818535863L, 15559462363L, 16337435383L,
	17154304303L, 18012019033L, 18912619603L, 19858249549L, 20851162021L,
	21893719921L, 22988405833L, 24137825713L, 25344716401L, 26611952023L,
	27942548563L, 29339675203L, 30806658811L, 32346991741L, 33964340653L,
	35662557061L, 37445682481L, 39317966071L, 41283863683L, 43348056571L,
	45515459221L, 47791229623L, 50180790763L, 52689830191L, 55324321471L,
	58090535023L, 60995061511L, 64044814111L, 67247053933L, 70609405513L,
	74139875629L, 77846867719L, 81739211053L, 85826171569L, 90117479719L,
	94623353611L, 99354520999L, 104322246649L, 109538358781L, 115015276639L,
	120766040323L, 126804341893L, 133144558123L, 139801785943L, 146791874503L,
	154131467131L, 161838040339L, 169929942253L, 178426439299L, 187347761059L,
	196715147611L, 206550902629L, 216878447329L, 227722369339L, 239108487523L,
	251063911393L, 263617106239L, 276797961259L, 290637858673L, 305169751459L,
	320428238989L, 336449650219L, 353272132573L, 370935738499L, 389482524943L,
	408956649469L, 429404481559L, 450874704883L, 473418439609L, 497089361581L,
	521943829573L, 548041020721L, 575443071409L, 604215224551L, 634425985759L,
	666147282601L, 699454645993L, 734427378103L, 771148746229L, 809706182809L,
	850191491911L, 892701066403L, 937336119691L, 984202925611L, 1033413071611L,
	1085083725139L, 1139337910231L, 1196304805513L, 1256120045179L,
	1318926046771L, 1384872349033L, 1454115966343L, 1526821763179L,
	1603162851163L, 1683320993701L, 1767487043281L, 1855861394839L,
	1948654464571L, 2046087187561L, 2148391546201L, 2255811123229L,
	2368601679049L, 2487031762351L, 2611383349123L, 2741952516463L,
	2879050142233L, 3023002649209L, 3174152780683L, 3332860418863L,
	3499503438721L, 3674478609619L, 3858202539553L, 4051112665363L,
	4253668298569L, 4466351712343L, 4689669297121L, 4924152760951L,
	5170360398871L, 5428878418111L, 5700322338679L, 5985338454673L,
	6284605376983L, 6598835645251L, 6928777426579L, 7275216297553L,
	7638977112073L, 8020925966743L, 8421972264949L, 8843070877009L,
	9285224420533L, 9749485640683L, 10236959922619L, 10748807917999L,
	11286248313631L, 11850560728873L, 12443088765229L, 13065243203029L,
	13718505362413L, 14404430629621L, 15124652159923L, 15880884767119L,
	16674929004691L, 17508675454561L, 18384109226773L, 19303314687661L,
	20268480422041L, 21281904441859L, 22345999663861L, 23463299646301L,
	24636464626651L, 25868287857613L, 27161702249623L, 28519787361073L,
	29945776728871L, 31443065564923L, 33015218842561L, 34665979784323L,
	36399278773369L, 38219242711423L, 40130204846401L, 42136715087683L,
	44243550841573L, 46455728382823L, 48778514801713L, 51217440541681L,
	53778312568273L, 56467228196209L, 59290589605441L, 62255119083403L,
	65367875036911L, 68636268788083L, 72068082226771L, 75671486337751L,
	79455060651391L, 83427813682459L, 87599204365963L, 91979164583389L,
	96578122811629L, 101407028951749L, 106477380398581L, 111801249418381L,
	117391311889171L, 123260877482989L, 129423921356869L, 135895117423759L,
	142689873294391L, 149824366957723L, 157315585305553L, 165181364570713L,
	173440432797391L, 182112454437103L, 191218077158161L, 200778981015229L,
	210817930064869L, 221358826564963L, 232426767892591L, 244048106287189L,
	256250511601291L, 269063037181111L, 282516189040069L, 296641998489379L,
	311474098412341L, 327047803332493L, 343400193498883L, 360570203173129L,
	378598713329569L, 397528648995091L, 417405081444811L, 438275335516903L,
	460189102292461L, 483198557406979L, 507358485275689L, 532726409539063L,
	559362730015711L, 587330866516411L, 616697409841309L, 647532280332391L,
	679908894347161L, 713904339063823L, 749599556016823L, 787079533813429L,
	826433510503831L, 867755186028253L, 911142945327949L, 956700092593771L,
	1004535097222273L, 1054761852083113L, 1107499944685039L, 1162874941918513L,
	1221018689014429L, 1282069623462961L, 1346173104633829L, 1413481759865509L,
	1484155847858683L, 1558363640250409L, 1636281822261421L, 1718095913373613L,
	1804000709041543L, 1894200744492733L, 1988910781716979L, 2088356320802581L,
	2192774136837901L, 2302412843679721L, 2417533485863611L, 2538410160156331L,
	2665330668163621L, 2798597201570449L, 2938527061648813L, 3085453414731223L,
	3239726085467581L, 3401712389739493L, 3571798009224793L, 3750387909685933L,
	3937907305169503L, 4134802670427643L, 4341542803947763L, 4558619944144753L,
	4786550941350793L, 5025878488418239L, 5277172412838241L, 5541031033479973L,
	5818082585151829L, 6108986714409409L, 6414436050127759L, 6735157852632241L,
	7071915745260913L, 7425511532522983L, 7796787109147669L, 8186626464604573L,
	8595957787833721L, 9025755677222701L, 9477043461083581L, 9950895634134841L,
	10448440415839693L, 10970862436630939L, 11519405558461063L,
	12095375836382221L, 12700144628200483L, 13335151859609851L,
	14001909452590141L, 14702004925217923L, 15437105171477401L,
	16208960430048979L, 17019408451549663L, 17870378874126133L,
	18763897817830393L, 19702092708721561L, 20687197344157579L,
	21721557211362139L, 22807635071930101L, 23948016825524251L,
	25145417666799979L, 26402688550138993L, 27722822977645669L,
	29108964126526519L, 30564412332851533L, 32092632949493869L,
	33697264596968521L, 35382127826815741L, 37151234218156453L,
	39008795929062781L, 40959235725511513L, 43007197511786173L,
	45157557387373681L, 47415435256740379L, 49786207019575321L,
	52275517370552743L, 54889293239080033L, 57633757901033443L,
	60515445796082611L, 63541218085884829L, 66718278990176071L,
	70054192939684843L, 73556902586668771L, 77234747716002061L,
	81096485101801909L, 85151309356891693L, 89408874824728471L,
	93879318565964671L, 98573284494261781L, 103501948718973733L,
	108677046154922281L, 114110898462666361L, 119816443385799541L,
	125807265555085729L, 132097628832838831L, 138702510274479409L,
	145637635788203203L, 152919517577610211L, 160565493456489889L,
	168593768129314393L, 177023456535779353L, 185874629362567699L,
	195168360830694109L, 204926778872228263L, 215173117815838381L,
	225931773706628389L, 237228362391959599L, 249089780511553453L,
	261544269537130579L, 274621483013986711L, 288352557164685031L,
	302770185022918213L, 317908694274062653L, 333804128987765479L,
	350494335437151013L, 368019052209007741L, 386420004819456331L,
	405741005060428693L, 426028055313449341L, 447329458079118319L,
	469695930983072839L, 493180727532224659L, 517839763908833833L,
	543731752104274861L, 570918339709487131L, 599464256694960871L,
	629437469529708223L, 660909343006190899L, 693954810156499303L,
	728652550664320351L, 765085178197535851L, 803339437107410461L,
	843506408962777591L, 885681729410916439L, 929965815881462293L,
	976464106675535479L, 1025287312009307773L, 1076551677609771973L,
	1130379261490259239L, 1186898224564770163L, 1246243135793007493L,
	1308555292582657939L, 1373983057211790409L, 1442682210072376321L,
	1514816320575993133L, 1590557136604791001L, 1670084993435030329L,
	1753589243106781699L, 1841268705262119643L, 1933332140525224759L,
	2029998747551483719L, 2131498684929056983L, 2238073619175509449L,
	2349977300134282729L, 2467476165140995279L, 2590849973398042621L,
	2720392472067944623L, 2856412095671341201L, 2999232700454907349L,
	3149194335477651103L, 3306654052251531673L, 3471986754864108193L,
	3645586092607310851L, 3827865397237676521L, 4019258667099560083L,
	4220221600454537311L, 4431232680477263053L, 4652794314501125203L,
	4885434030226180201L, 5129705731737487051L, 5386191018324357391L,
	5655500569240573411L, 5938275597702599353L, 6235189377587729959L,
	6546948846467114779L, 6874296288790467079L, 7218011103229988773L,
	7578911658391486813L, 7957857241311059629L, 8355750103376609851L,
	8773537608545440321L, 9212214488972709751L
};
#endif

/* Prototypes */

/* Implementation */

mapkit_size_t mapkit_nextprime(const mapkit_size_t n)
{
	long low, high, mid;

	if (n >= mapkit_primes[MAPKIT_PRIMES - 1])
		return mapkit_primes[MAPKIT_PRIMES - 1];

	for (low = -1, high = MAPKIT_PRIMES - 1; high - low > 1;) {
		mid = (high + low) >> 1;
		if (n <= mapkit_primes[mid])
			high = mid;
		else
			low = mid;
	}

	return mapkit_primes[high];
}

/* similar to the hash function used in PERL */
mapkit_hash_t mapkit_strhash(const char *string)
{
	unsigned char c;
	const unsigned char *ustring = (const unsigned char *) string;
	mapkit_hash_t hash = 5381;			       /* bernstein constant */

	for (; (c = *ustring++);)
		hash = hash * 33 + c;

	return (hash + (hash >> 5));
}

/* similar to the hash function used in PERL */
mapkit_hash_t mapkit_memhash(const void *data, const size_t len)
{
	size_t i;
	const unsigned char *ustring = (const unsigned char *) data;
	mapkit_hash_t hash = 5381;			       /* bernstein constant */

	for (i = 0; i < len; i++)
		hash = hash * 33 + *ustring++;

	return (hash + (hash >> 5));
}

/*  */

/* MapKit */

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

/* static char const rcsid[] =
   "@(#) $Jeannot: mapkit.m4c,v 1.23 2004/03/31 19:12:55 js Exp $"; */

#include <stdlib.h>
#include <stdio.h>
#include <limits.h>

/*
  MapKit
  Copyright (c) 2002-2004, Jean-Sebastien Roy (js@jeannot.org)
  @(#) $Jeannot: mapkit.m4,v 1.83 2004/04/10 14:17:03 js Exp $
*/

/* @(#) $Jeannot: mapkit_defs.m4,v 1.20 2004/03/31 19:12:55 js Exp $ */

/* Map structures */

#ifdef MAPKIT_map_ii

/*
  Implementation for map_ii (int -> int)
  Default value : 0
  Uses a state field.
*/

/* Static prototypes */

/* Return the index of key, or MAPKIT_KEYNOTFOUND */
static mapkit_size_t map_ii_keyindex(map_ii * spm, int key);

/* Return the index of key or -(insertion index)-1 if key not found */
static mapkit_size_t map_ii_insertionindex(map_ii * spm, int key);

/* Implementation */

mapkit_error map_ii_init(map_ii * spm)
{
	return map_ii_init_hint(spm, MAPKIT_DEFAULT_EXPECTEDUSED);
}

mapkit_error map_ii_init_hint(map_ii * spm, mapkit_size_t used)
{
#ifdef MAPKIT_DEBUG
	fprintf(stderr, "MAPKIT: init\n");
#endif

	spm->size = 0;
	spm->fill = 0;
	spm->used = 0;
	spm->maxfillfactor = 0.5;
	spm->minusedfactor = 0.2;
	spm->contents = NULL;
#ifdef MAPKIT_COLLISIONS
	spm->insertionindexs = spm->insertionindex_collisions = 0;
	spm->keyindexs = spm->keyindex_collisions = 0;
#endif
	spm->defaultvalue = 0;
	spm->alwaysdefault = 0;

	return map_ii_reallocate(spm, map_ii_meansize(spm, used));
}

mapkit_error map_ii_ensurecapacity(map_ii * spm, mapkit_size_t used)
{
	if (used > (spm->used + spm->maxfill - spm->fill)) {
#ifdef MAPKIT_DEBUG
		fprintf(stderr, "MAPKIT: need more capacity\n");
#endif
		return map_ii_reallocate(spm, map_ii_meansize(spm, used));
	} else
		return MAPKIT_OK;
}

mapkit_error map_ii_adjustcapacity(map_ii * spm)
{
	spm->minused = (mapkit_size_t) (spm->size * spm->minusedfactor);
	spm->maxfill = (mapkit_size_t) (spm->size * spm->maxfillfactor);

	if (spm->used < spm->minused) {
#ifdef MAPKIT_DEBUG
		fprintf(stderr, "MAPKIT: used < minused\n");
#endif
		return map_ii_reallocate(spm, map_ii_meansize(spm, spm->used));
	} else if (spm->fill > spm->maxfill) {
#ifdef MAPKIT_DEBUG
		fprintf(stderr, "MAPKIT: fill > maxfill\n");
#endif
		return map_ii_reallocate(spm, map_ii_meansize(spm, spm->used));
	} else
		return MAPKIT_OK;
}

void map_ii_free(map_ii * spm)
{
	free(spm->contents);
	spm->contents = NULL;
	spm->size = 0;
	spm->fill = 0;
	spm->used = 0;
	spm->maxfill = 0;
	spm->minused = 0;
#ifdef MAPKIT_COLLISIONS
	spm->insertionindexs = spm->insertionindex_collisions = 0;
	spm->keyindexs = spm->keyindex_collisions = 0;
#endif
}

mapkit_error map_ii_copy(map_ii * to, map_ii * from)
{
	map_ii_storage *contentscopy;
	size_t size = from->size * sizeof(*from->contents);

	contentscopy = (map_ii_storage *) malloc(size);
	if (contentscopy == NULL)
		MAPKIT_ERROR(MAPKIT_ENOMEM);

	memcpy(to, from, sizeof(map_ii));
	to->contents = contentscopy;
	memcpy(to->contents, from->contents, size);

	return MAPKIT_OK;
}

mapkit_size_t map_ii_growsize(map_ii * spm, mapkit_size_t used)
{
	return (mapkit_size_t) (4.0 * used / (3 * spm->minusedfactor + spm->maxfillfactor));
}

mapkit_size_t map_ii_shrinksize(map_ii * spm, mapkit_size_t used)
{
	return (mapkit_size_t) (4.0 * used / (spm->minusedfactor + 3 * spm->maxfillfactor));
}

mapkit_size_t map_ii_meansize(map_ii * spm, mapkit_size_t used)
{
	return (mapkit_size_t) (2.0 * used / (spm->minusedfactor + spm->maxfillfactor));
}

mapkit_error map_ii_reallocate(map_ii * spm, mapkit_size_t newsize)
{
	mapkit_size_t iindex;
	mapkit_size_t oldsize;
	map_ii_storage *newcontents, *oldcontents;

	/*
	 * At least one free entry 
	 */
	if (newsize <= spm->used)
		newsize = spm->used + 1;
	newsize = mapkit_nextprime(newsize);
	if (newsize <= spm->used)
		MAPKIT_ERROR(MAPKIT_ETOOBIG);

	newcontents = (map_ii_storage *) malloc(newsize * sizeof(*spm->contents));
	if (newcontents == NULL)
		MAPKIT_ERROR(MAPKIT_ENOMEM);

	/*
	 * Initialize all entries to "free" 
	 */
	for (iindex = 0; iindex < newsize; iindex++)
		newcontents[iindex].state = MAPKIT_FREESLOT;

	oldcontents = spm->contents;
	oldsize = spm->size;
	spm->contents = newcontents;
	spm->size = newsize;

#ifdef MAPKIT_DEBUG
	fprintf(stderr, "MAPKIT: reallocate %ld -> %ld\n", (long) oldsize, (long) newsize);
#endif

	spm->maxfill = (mapkit_size_t) (newsize * spm->maxfillfactor);
	/*
	 * At least one free entry 
	 */
	if (spm->maxfill >= newsize)
		spm->maxfill = newsize - 1;
	spm->minused = (mapkit_size_t) (newsize * spm->minusedfactor);
	spm->used = 0;

	if (oldcontents != NULL) {
		int used = 0;
		int key;
		int defaultvalue = spm->defaultvalue;
		int notalwaysdefault = !spm->alwaysdefault;

		/*
		 * Copy all entries from old to new 
		 */
		for (iindex = 0; iindex < oldsize; iindex++)
			if (oldcontents[iindex].state == MAPKIT_FULLSLOT) {
				mapkit_size_t ins_iindex;
				map_ii_storage *contents;

				key = oldcontents[iindex].key;

				/*
				 * Fast path 
				 */
				ins_iindex = ((mapkit_hash_t) key) % spm->size;
				contents = &(newcontents[ins_iindex]);

				if (contents->state != MAPKIT_FREESLOT) {
					ins_iindex = map_ii_insertionindex(spm, key);
					contents = &(newcontents[ins_iindex]);
				}
#ifdef MAPKIT_COLLISIONS
				else
					spm->insertionindexs++;
#endif
				if (notalwaysdefault || (!((oldcontents[iindex].value) == (defaultvalue)))) {
					contents->value = oldcontents[iindex].value;
					contents->state = MAPKIT_FULLSLOT;
					contents->key = key;
					used++;
				}
			}
		free(oldcontents);
		spm->used = used;
	}
	spm->fill = spm->used;

	return MAPKIT_OK;
}

int map_ii_value_s(map_ii * spm, int key)
{
	mapkit_size_t iindex;

	iindex = map_ii_keyindex(spm, key);

	if (iindex < 0) {
		if (spm->alwaysdefault)
			return spm->defaultvalue;
		else
			MAPKIT_FATAL_ERROR(MAPKIT_EKEYNOTFOUND);
	}

	return spm->contents[iindex].value;
}

mapkit_error map_ii_get_s(map_ii * spm, int key, int *value)
{
	mapkit_size_t iindex;

	iindex = map_ii_keyindex(spm, key);

	if (iindex < 0) {
		if (spm->alwaysdefault) {
			*value = spm->defaultvalue;
			return MAPKIT_OK;
		} else {
			MAPKIT_ERROR(MAPKIT_EKEYNOTFOUND);
		}
	}

	*value = spm->contents[iindex].value;
	return MAPKIT_OK;
}

mapkit_error map_ii_set_s(map_ii * spm, int key, int value)
{
	mapkit_size_t iindex;

	iindex = map_ii_insertionindex(spm, key);

	if (iindex < 0)
		/*
		 * FULLSLOT 
		 */
		spm->contents[-iindex - 1].value = value;
	else {
		map_ii_storage *element = &(spm->contents[iindex]);
		int ffree = element->state == MAPKIT_FREESLOT;

		element->state = MAPKIT_FULLSLOT;
		element->value = value;
		element->key = key;
		spm->used++;

		if (ffree && ((++spm->fill) > spm->maxfill)) {
#ifdef MAPKIT_DEBUG
			fprintf(stderr, "MAPKIT: fill > maxfill\n");
#endif
			return map_ii_reallocate(spm, map_ii_growsize(spm, spm->used));
		}
	}
	return MAPKIT_OK;
}

int *map_ii_insertptr_s(map_ii * spm, int key)
{
	mapkit_size_t iindex;

	iindex = map_ii_insertionindex(spm, key);

	if (iindex < 0)
		return &(spm->contents[-iindex - 1].value);
	else {
		map_ii_storage *element = &(spm->contents[iindex]);

		if (element->state == MAPKIT_FREESLOT) {
			/*
			 * FREESLOT 
			 */
			if (spm->fill >= spm->maxfill) {
				mapkit_error err;

#ifdef MAPKIT_DEBUG
				fprintf(stderr, "MAPKIT: fill => maxfill before insert\n");
#endif
				/*
				 * Must reallocate -before- inserting defaultvalue 
				 */
				err = map_ii_reallocate(spm, map_ii_growsize(spm, spm->used + 1));
				if (err) {
					MAPKIT_ERROR_NORET(err);
					return NULL;
				}

				iindex = map_ii_insertionindex(spm, key);
				/*
				 * FREESLOT 
				 */
				spm->contents[iindex].state = MAPKIT_FULLSLOT;
				spm->contents[iindex].key = key;
				spm->contents[iindex].value = spm->defaultvalue;
				spm->used++;
				spm->fill++;
				return &(spm->contents[iindex].value);
			} else
				spm->fill++;
		}

		element->state = MAPKIT_FULLSLOT;
		element->key = key;
		element->value = spm->defaultvalue;
		spm->used++;
		return &(element->value);
	}
}

int *map_ii_ptr_s(map_ii * spm, int key)
{
	mapkit_size_t iindex;

	iindex = map_ii_keyindex(spm, key);

	if ((iindex >= 0) && ((!spm->alwaysdefault) || (!((spm->contents[iindex].value) == (spm->defaultvalue)))))
		return &(spm->contents[iindex].value);
	else
		return NULL;
}

mapkit_error map_ii_remove_s(map_ii * spm, int key)
{
	mapkit_size_t iindex;

	iindex = map_ii_keyindex(spm, key);

	if (iindex < 0) {
		if (spm->alwaysdefault)
			return MAPKIT_OK;
		else
			return MAPKIT_EKEYNOTFOUND;
	}

	spm->contents[iindex].state = MAPKIT_DELETEDSLOT;
	spm->used--;
	if (spm->used < spm->minused) {
#ifdef MAPKIT_DEBUG
		fprintf(stderr, "MAPKIT: used < minused\n");
#endif
		return map_ii_reallocate(spm, map_ii_shrinksize(spm, spm->used));
	}

	return MAPKIT_OK;
}

mapkit_size_t map_ii_keyindex(map_ii * spm, int key)
{
	mapkit_size_t iindex, decrement;

	signed char state;

	iindex = ((mapkit_hash_t) key) % spm->size;
	decrement = (((mapkit_hash_t) key) % (spm->size - 2));
	decrement += (decrement == 0);

#ifdef MAPKIT_COLLISIONS
	spm->keyindexs++;
#endif

	while ((state = spm->contents[iindex].state) != MAPKIT_FREESLOT
	       && (state == MAPKIT_DELETEDSLOT || (!((spm->contents[iindex].key) == (key))))) {
#ifdef MAPKIT_COLLISIONS
		spm->keyindex_collisions++;
#endif
		iindex -= decrement;
		if (iindex < 0)
			iindex += spm->size;
	}

	if (state == MAPKIT_FREESLOT)
		return MAPKIT_KEYNOTFOUND;
	return iindex;
}

mapkit_size_t map_ii_insertionindex(map_ii * spm, int key)
{
	mapkit_size_t iindex, decrement;
	signed char state;

#ifdef MAPKIT_COLLISIONS
	spm->insertionindexs++;
#endif

	iindex = ((mapkit_hash_t) key) % spm->size;

	/*
	 * Fast path (largely superfluous) 
	 */
	if ((state = spm->contents[iindex].state) == MAPKIT_FREESLOT)
		return iindex;
	if ((state == MAPKIT_FULLSLOT)
	    && ((spm->contents[iindex].key) == (key)))
		return -iindex - 1;

	decrement = (((mapkit_hash_t) key) % (spm->size - 2));
	decrement += (decrement == 0);

	while ((state == MAPKIT_FULLSLOT)
	       && (!((spm->contents[iindex].key) == (key)))) {
#ifdef MAPKIT_COLLISIONS
		spm->insertionindex_collisions++;
#endif
		iindex -= decrement;
		if (iindex < 0)
			iindex += spm->size;
		state = spm->contents[iindex].state;
	}

	if (state == MAPKIT_DELETEDSLOT) {
		mapkit_size_t index2 = iindex;

		while ((state = spm->contents[iindex].state) != MAPKIT_FREESLOT && ((state == MAPKIT_DELETEDSLOT)
										    || (!((spm->contents[iindex].key) == (key))))) {
			iindex -= decrement;
			if (iindex < 0)
				iindex += spm->size;
			state = spm->contents[iindex].state;
		}
		if (state == MAPKIT_FREESLOT)
			return index2;
	}

	if (state == MAPKIT_FULLSLOT)
		return -iindex - 1;
	return iindex;
}

mapkit_size_t map_ii_next(map_ii * spm, mapkit_size_t iindex)
{
	mapkit_size_t size = spm->size;
	map_ii_storage *pos_contents;
	int defaultvalue = spm->defaultvalue;
	int notalwaysdefault = !spm->alwaysdefault;

	pos_contents = &(spm->contents[++iindex]);

	for (; iindex < size; iindex++, pos_contents++)
		if ((pos_contents->state == MAPKIT_FULLSLOT)
		    && (notalwaysdefault || (!((pos_contents->value) == (defaultvalue)))))
			return iindex;

	return -1;
}

map_ii_storage *map_ii_nextptr(map_ii * spm, map_ii_storage * pos_contents)
{
	map_ii_storage *end = &(spm->contents[spm->size]);
	int defaultvalue = spm->defaultvalue;
	int notalwaysdefault = !spm->alwaysdefault;

	if (pos_contents == NULL)
		pos_contents = spm->contents;
	else {
		pos_contents++;
		if (pos_contents <= spm->contents)
			return NULL;
	}

	for (; pos_contents < end; pos_contents++)
		if ((pos_contents->state == MAPKIT_FULLSLOT)
		    && (notalwaysdefault || (!((pos_contents->value) == (defaultvalue)))))
			return pos_contents;

	return NULL;
}

mapkit_error map_ii_getall(map_ii * spm, map_ii_element ** array, mapkit_size_t * count)
{
	mapkit_size_t iindex;
	mapkit_size_t size = spm->size, vcount = 0;
	map_ii_element *pos_array;
	map_ii_storage *pos_contents;
	int defaultvalue = spm->defaultvalue;
	int notalwaysdefault = !spm->alwaysdefault;

	pos_array = *array = (map_ii_element *) malloc(sizeof(**array) * spm->used);
	if (*array == NULL)
		MAPKIT_ERROR(MAPKIT_ENOMEM);

	pos_contents = spm->contents;

	for (iindex = 0; iindex < size; iindex++, pos_contents++)
		if ((pos_contents->state == MAPKIT_FULLSLOT)
		    && (notalwaysdefault || (!((pos_contents->value) == (defaultvalue))))) {
			pos_array->key = pos_contents->key;
			pos_array->value = pos_contents->value;
			pos_array++;
			vcount++;
		}
	*count = vcount;

	return MAPKIT_OK;
}

mapkit_error map_ii_clean(map_ii * spm)
{
	mapkit_size_t iindex, count = 0;
	mapkit_size_t size = spm->size;
	map_ii_storage *pos_contents;
	int defaultvalue = spm->defaultvalue;

	pos_contents = spm->contents;
	for (iindex = 0; iindex < size; iindex++, pos_contents++)
		if ((pos_contents->state == MAPKIT_FULLSLOT)
		    && ((pos_contents->value) == (defaultvalue))) {
			pos_contents->state = MAPKIT_DELETEDSLOT;
			count++;
		}

	spm->used -= count;
	if (spm->used < spm->minused) {
#ifdef MAPKIT_DEBUG
		fprintf(stderr, "MAPKIT: used < minused\n");
#endif
		return map_ii_reallocate(spm, map_ii_meansize(spm, spm->used));
	}

	return MAPKIT_OK;
}

int map_ii_compare(const void *e1, const void *e2)
{
	int key1 = ((const map_ii_element *) e1)->key;
	int key2 = ((const map_ii_element *) e2)->key;

	return ((key1) < (key2) ? -1 : ((key1) == (key2) ? 0 : 1));
}

mapkit_error map_ii_getall_sorted(map_ii * spm, map_ii_element ** array, mapkit_size_t * count)
{
	mapkit_error err;

	err = map_ii_getall(spm, array, count);
	if (err)
		MAPKIT_ERROR(err);

	qsort(*array, *count, sizeof(**array), map_ii_compare);

	return MAPKIT_OK;
}

mapkit_error map_ii_setall(map_ii * spm, map_ii_element * array, mapkit_size_t count)
{
	mapkit_size_t array_iindex;
	mapkit_error err;

	err = map_ii_ensurecapacity(spm, spm->used + count);
	if (err)
		MAPKIT_ERROR(err);

	if (spm->alwaysdefault)
		/*
		 * Prevent shrinking 
		 */
		spm->minused = 0;

	for (array_iindex = 0; array_iindex < count; array_iindex++) {
		err = map_ii_set(spm, array[array_iindex].key, array[array_iindex].value);
		if (err)
			MAPKIT_ERROR(err);
	}

	if (spm->alwaysdefault)
		map_ii_adjustcapacity(spm);

	return MAPKIT_OK;
}

mapkit_error map_ii_removeall(map_ii * spm, int *array, mapkit_size_t count)
{
	mapkit_size_t array_iindex;
	mapkit_error err;

	/*
	 * Prevent shrinking 
	 */
	spm->minused = 0;

	for (array_iindex = 0; array_iindex < count; array_iindex++) {
		err = map_ii_remove(spm, array[array_iindex]);
		if (err)
			MAPKIT_ERROR(err);
	}

	map_ii_adjustcapacity(spm);

	return MAPKIT_OK;
}

void map_ii_printstats(map_ii * spm)
{
	fprintf(stderr, "MAPKIT: map_ii statistics\n");
	fprintf(stderr, "MAPKIT: alwaysdefault = %d\n", spm->alwaysdefault);
	fprintf(stderr, "MAPKIT: minused = %ld, maxfill = %ld\n", (long) spm->minused, (long) spm->maxfill);
	fprintf(stderr, "MAPKIT: minusedfactor = %g, maxfillfactor = %g\n", spm->minusedfactor, spm->maxfillfactor);
#ifdef MAPKIT_COLLISIONS
	fprintf(stderr, "MAPKIT: insertionindexs = %lu, collisions = %lu\n", (unsigned long) spm->insertionindexs,
		(unsigned long) spm->insertionindex_collisions);
	fprintf(stderr, "MAPKIT: keyindexs = %lu, collisions = %lu\n", (unsigned long) spm->keyindexs, (unsigned long) spm->keyindex_collisions);
#endif
}

#endif							       /* MAPKIT_map_ii */

#ifdef MAPKIT_map_id

/*
  Implementation for map_id (int -> double)
  Default value : 0.0
  Uses a state field.
*/

/* Static prototypes */

/* Return the iindex of key, or MAPKIT_KEYNOTFOUND */
static mapkit_size_t map_id_keyindex(map_id * spm, int key);

/* Return the iindex of key or -(insertion iindex)-1 if key not found */
static mapkit_size_t map_id_insertionindex(map_id * spm, int key);

/* Implementation */

mapkit_error map_id_init(map_id * spm)
{
	return map_id_init_hint(spm, MAPKIT_DEFAULT_EXPECTEDUSED);
}

mapkit_error map_id_init_hint(map_id * spm, mapkit_size_t used)
{
#ifdef MAPKIT_DEBUG
	fprintf(stderr, "MAPKIT: init\n");
#endif

	spm->size = 0;
	spm->fill = 0;
	spm->used = 0;
	spm->maxfillfactor = 0.5;
	spm->minusedfactor = 0.2;
	spm->contents = NULL;
#ifdef MAPKIT_COLLISIONS
	spm->insertionindexs = spm->insertionindex_collisions = 0;
	spm->keyindexs = spm->keyindex_collisions = 0;
#endif
	spm->defaultvalue = 0.0;
	spm->alwaysdefault = 0;

	return map_id_reallocate(spm, map_id_meansize(spm, used));
}

mapkit_error map_id_ensurecapacity(map_id * spm, mapkit_size_t used)
{
	if (used > (spm->used + spm->maxfill - spm->fill)) {
#ifdef MAPKIT_DEBUG
		fprintf(stderr, "MAPKIT: need more capacity\n");
#endif
		return map_id_reallocate(spm, map_id_meansize(spm, used));
	} else
		return MAPKIT_OK;
}

mapkit_error map_id_adjustcapacity(map_id * spm)
{
	spm->minused = (mapkit_size_t) (spm->size * spm->minusedfactor);
	spm->maxfill = (mapkit_size_t) (spm->size * spm->maxfillfactor);

	if (spm->used < spm->minused) {
#ifdef MAPKIT_DEBUG
		fprintf(stderr, "MAPKIT: used < minused\n");
#endif
		return map_id_reallocate(spm, map_id_meansize(spm, spm->used));
	} else if (spm->fill > spm->maxfill) {
#ifdef MAPKIT_DEBUG
		fprintf(stderr, "MAPKIT: fill > maxfill\n");
#endif
		return map_id_reallocate(spm, map_id_meansize(spm, spm->used));
	} else
		return MAPKIT_OK;
}

void map_id_free(map_id * spm)
{
	free(spm->contents);
	spm->contents = NULL;
	spm->size = 0;
	spm->fill = 0;
	spm->used = 0;
	spm->maxfill = 0;
	spm->minused = 0;
#ifdef MAPKIT_COLLISIONS
	spm->insertionindexs = spm->insertionindex_collisions = 0;
	spm->keyindexs = spm->keyindex_collisions = 0;
#endif
}

mapkit_error map_id_copy(map_id * to, map_id * from)
{
	map_id_storage *contentscopy;
	size_t size = from->size * sizeof(*from->contents);

	contentscopy = (map_id_storage *) malloc(size);
	if (contentscopy == NULL)
		MAPKIT_ERROR(MAPKIT_ENOMEM);

	memcpy(to, from, sizeof(map_id));
	to->contents = contentscopy;
	memcpy(to->contents, from->contents, size);

	return MAPKIT_OK;
}

mapkit_size_t map_id_growsize(map_id * spm, mapkit_size_t used)
{
	return (mapkit_size_t) (4.0 * used / (3 * spm->minusedfactor + spm->maxfillfactor));
}

mapkit_size_t map_id_shrinksize(map_id * spm, mapkit_size_t used)
{
	return (mapkit_size_t) (4.0 * used / (spm->minusedfactor + 3 * spm->maxfillfactor));
}

mapkit_size_t map_id_meansize(map_id * spm, mapkit_size_t used)
{
	return (mapkit_size_t) (2.0 * used / (spm->minusedfactor + spm->maxfillfactor));
}

mapkit_error map_id_reallocate(map_id * spm, mapkit_size_t newsize)
{
	mapkit_size_t iindex;
	mapkit_size_t oldsize;
	map_id_storage *newcontents, *oldcontents;

	/*
	 * At least one free entry 
	 */
	if (newsize <= spm->used)
		newsize = spm->used + 1;
	newsize = mapkit_nextprime(newsize);
	if (newsize <= spm->used)
		MAPKIT_ERROR(MAPKIT_ETOOBIG);

	newcontents = (map_id_storage *) malloc(newsize * sizeof(*spm->contents));
	if (newcontents == NULL)
		MAPKIT_ERROR(MAPKIT_ENOMEM);

	/*
	 * Initialize all entries to "free" 
	 */
	for (iindex = 0; iindex < newsize; iindex++)
		newcontents[iindex].state = MAPKIT_FREESLOT;

	oldcontents = spm->contents;
	oldsize = spm->size;
	spm->contents = newcontents;
	spm->size = newsize;

#ifdef MAPKIT_DEBUG
	fprintf(stderr, "MAPKIT: reallocate %ld -> %ld\n", (long) oldsize, (long) newsize);
#endif

	spm->maxfill = (mapkit_size_t) (newsize * spm->maxfillfactor);
	/*
	 * At least one free entry 
	 */
	if (spm->maxfill >= newsize)
		spm->maxfill = newsize - 1;
	spm->minused = (mapkit_size_t) (newsize * spm->minusedfactor);
	spm->used = 0;

	if (oldcontents != NULL) {
		int used = 0;
		int key;
		double defaultvalue = spm->defaultvalue;
		int notalwaysdefault = !spm->alwaysdefault;

		/*
		 * Copy all entries from old to new 
		 */
		for (iindex = 0; iindex < oldsize; iindex++)
			if (oldcontents[iindex].state == MAPKIT_FULLSLOT) {
				mapkit_size_t ins_iindex;
				map_id_storage *contents;

				key = oldcontents[iindex].key;

				/*
				 * Fast path 
				 */
				ins_iindex = ((mapkit_hash_t) key) % spm->size;
				contents = &(newcontents[ins_iindex]);

				if (contents->state != MAPKIT_FREESLOT) {
					ins_iindex = map_id_insertionindex(spm, key);
					contents = &(newcontents[ins_iindex]);
				}
#ifdef MAPKIT_COLLISIONS
				else
					spm->insertionindexs++;
#endif
				if (notalwaysdefault || (!((oldcontents[iindex].value) == (defaultvalue)))) {
					contents->value = oldcontents[iindex].value;
					contents->state = MAPKIT_FULLSLOT;
					contents->key = key;
					used++;
				}
			}
		free(oldcontents);
		spm->used = used;
	}
	spm->fill = spm->used;

	return MAPKIT_OK;
}

double map_id_value_s(map_id * spm, int key)
{
	mapkit_size_t iindex;

	iindex = map_id_keyindex(spm, key);

	if (iindex < 0) {
		if (spm->alwaysdefault)
			return spm->defaultvalue;
		else
			MAPKIT_FATAL_ERROR(MAPKIT_EKEYNOTFOUND);
	}

	return spm->contents[iindex].value;
}

mapkit_error map_id_get_s(map_id * spm, int key, double *value)
{
	mapkit_size_t iindex;

	iindex = map_id_keyindex(spm, key);

	if (iindex < 0) {
		if (spm->alwaysdefault) {
			*value = spm->defaultvalue;
			return MAPKIT_OK;
		} else {
			MAPKIT_ERROR(MAPKIT_EKEYNOTFOUND);
		}
	}

	*value = spm->contents[iindex].value;
	return MAPKIT_OK;
}

mapkit_error map_id_set_s(map_id * spm, int key, double value)
{
	mapkit_size_t iindex;

	iindex = map_id_insertionindex(spm, key);

	if (iindex < 0)
		/*
		 * FULLSLOT 
		 */
		spm->contents[-iindex - 1].value = value;
	else {
		map_id_storage *element = &(spm->contents[iindex]);
		int ffree = element->state == MAPKIT_FREESLOT;

		element->state = MAPKIT_FULLSLOT;
		element->value = value;
		element->key = key;
		spm->used++;

		if (ffree && ((++spm->fill) > spm->maxfill)) {
#ifdef MAPKIT_DEBUG
			fprintf(stderr, "MAPKIT: fill > maxfill\n");
#endif
			return map_id_reallocate(spm, map_id_growsize(spm, spm->used));
		}
	}
	return MAPKIT_OK;
}

double *map_id_insertptr_s(map_id * spm, int key)
{
	mapkit_size_t iindex;

	iindex = map_id_insertionindex(spm, key);

	if (iindex < 0)
		return &(spm->contents[-iindex - 1].value);
	else {
		map_id_storage *element = &(spm->contents[iindex]);

		if (element->state == MAPKIT_FREESLOT) {
			/*
			 * FREESLOT 
			 */
			if (spm->fill >= spm->maxfill) {
				mapkit_error err;

#ifdef MAPKIT_DEBUG
				fprintf(stderr, "MAPKIT: fill => maxfill before insert\n");
#endif
				/*
				 * Must reallocate -before- inserting defaultvalue 
				 */
				err = map_id_reallocate(spm, map_id_growsize(spm, spm->used + 1));
				if (err) {
					MAPKIT_ERROR_NORET(err);
					return NULL;
				}

				iindex = map_id_insertionindex(spm, key);
				/*
				 * FREESLOT 
				 */
				spm->contents[iindex].state = MAPKIT_FULLSLOT;
				spm->contents[iindex].key = key;
				spm->contents[iindex].value = spm->defaultvalue;
				spm->used++;
				spm->fill++;
				return &(spm->contents[iindex].value);
			} else
				spm->fill++;
		}

		element->state = MAPKIT_FULLSLOT;
		element->key = key;
		element->value = spm->defaultvalue;
		spm->used++;
		return &(element->value);
	}
}

double *map_id_ptr_s(map_id * spm, int key)
{
	mapkit_size_t iindex;

	iindex = map_id_keyindex(spm, key);

	if ((iindex >= 0) && ((!spm->alwaysdefault) || (!((spm->contents[iindex].value) == (spm->defaultvalue)))))
		return &(spm->contents[iindex].value);
	else
		return NULL;
}

mapkit_error map_id_remove_s(map_id * spm, int key)
{
	mapkit_size_t iindex;

	iindex = map_id_keyindex(spm, key);

	if (iindex < 0) {
		if (spm->alwaysdefault)
			return MAPKIT_OK;
		else
			return MAPKIT_EKEYNOTFOUND;
	}

	spm->contents[iindex].state = MAPKIT_DELETEDSLOT;
	spm->used--;
	if (spm->used < spm->minused) {
#ifdef MAPKIT_DEBUG
		fprintf(stderr, "MAPKIT: used < minused\n");
#endif
		return map_id_reallocate(spm, map_id_shrinksize(spm, spm->used));
	}

	return MAPKIT_OK;
}

mapkit_size_t map_id_keyindex(map_id * spm, int key)
{
	mapkit_size_t iindex, decrement;

	signed char state;

	iindex = ((mapkit_hash_t) key) % spm->size;
	decrement = (((mapkit_hash_t) key) % (spm->size - 2));
	decrement += (decrement == 0);

#ifdef MAPKIT_COLLISIONS
	spm->keyindexs++;
#endif

	while ((state = spm->contents[iindex].state) != MAPKIT_FREESLOT
	       && (state == MAPKIT_DELETEDSLOT || (!((spm->contents[iindex].key) == (key))))) {
#ifdef MAPKIT_COLLISIONS
		spm->keyindex_collisions++;
#endif
		iindex -= decrement;
		if (iindex < 0)
			iindex += spm->size;
	}

	if (state == MAPKIT_FREESLOT)
		return MAPKIT_KEYNOTFOUND;
	return iindex;
}

mapkit_size_t map_id_insertionindex(map_id * spm, int key)
{
	mapkit_size_t iindex, decrement;
	signed char state;

#ifdef MAPKIT_COLLISIONS
	spm->insertionindexs++;
#endif

	iindex = ((mapkit_hash_t) key) % spm->size;

	/*
	 * Fast path (largely superfluous) 
	 */
	if ((state = spm->contents[iindex].state) == MAPKIT_FREESLOT)
		return iindex;
	if ((state == MAPKIT_FULLSLOT)
	    && ((spm->contents[iindex].key) == (key)))
		return -iindex - 1;

	decrement = (((mapkit_hash_t) key) % (spm->size - 2));
	decrement += (decrement == 0);

	while ((state == MAPKIT_FULLSLOT)
	       && (!((spm->contents[iindex].key) == (key)))) {
#ifdef MAPKIT_COLLISIONS
		spm->insertionindex_collisions++;
#endif
		iindex -= decrement;
		if (iindex < 0)
			iindex += spm->size;
		state = spm->contents[iindex].state;
	}

	if (state == MAPKIT_DELETEDSLOT) {
		mapkit_size_t index2 = iindex;

		while ((state = spm->contents[iindex].state) != MAPKIT_FREESLOT && ((state == MAPKIT_DELETEDSLOT)
										    || (!((spm->contents[iindex].key) == (key))))) {
			iindex -= decrement;
			if (iindex < 0)
				iindex += spm->size;
			state = spm->contents[iindex].state;
		}
		if (state == MAPKIT_FREESLOT)
			return index2;
	}

	if (state == MAPKIT_FULLSLOT)
		return -iindex - 1;
	return iindex;
}

mapkit_size_t map_id_next(map_id * spm, mapkit_size_t iindex)
{
	mapkit_size_t size = spm->size;
	map_id_storage *pos_contents;
	double defaultvalue = spm->defaultvalue;
	int notalwaysdefault = !spm->alwaysdefault;

	pos_contents = &(spm->contents[++iindex]);

	for (; iindex < size; iindex++, pos_contents++)
		if ((pos_contents->state == MAPKIT_FULLSLOT)
		    && (notalwaysdefault || (!((pos_contents->value) == (defaultvalue)))))
			return iindex;

	return -1;
}

map_id_storage *map_id_nextptr(map_id * spm, map_id_storage * pos_contents)
{
	map_id_storage *end = &(spm->contents[spm->size]);
	double defaultvalue = spm->defaultvalue;
	int notalwaysdefault = !spm->alwaysdefault;

	if (pos_contents == NULL)
		pos_contents = spm->contents;
	else {
		pos_contents++;
		if (pos_contents <= spm->contents)
			return NULL;
	}

	for (; pos_contents < end; pos_contents++)
		if ((pos_contents->state == MAPKIT_FULLSLOT)
		    && (notalwaysdefault || (!((pos_contents->value) == (defaultvalue)))))
			return pos_contents;

	return NULL;
}

mapkit_error map_id_getall(map_id * spm, map_id_element ** array, mapkit_size_t * count)
{
	mapkit_size_t iindex;
	mapkit_size_t size = spm->size, vcount = 0;
	map_id_element *pos_array;
	map_id_storage *pos_contents;
	double defaultvalue = spm->defaultvalue;
	int notalwaysdefault = !spm->alwaysdefault;

	pos_array = *array = (map_id_element *) malloc(sizeof(**array) * spm->used);
	if (*array == NULL)
		MAPKIT_ERROR(MAPKIT_ENOMEM);

	pos_contents = spm->contents;

	for (iindex = 0; iindex < size; iindex++, pos_contents++)
		if ((pos_contents->state == MAPKIT_FULLSLOT)
		    && (notalwaysdefault || (!((pos_contents->value) == (defaultvalue))))) {
			pos_array->key = pos_contents->key;
			pos_array->value = pos_contents->value;
			pos_array++;
			vcount++;
		}
	*count = vcount;

	return MAPKIT_OK;
}

mapkit_error map_id_clean(map_id * spm)
{
	mapkit_size_t iindex, count = 0;
	mapkit_size_t size = spm->size;
	map_id_storage *pos_contents;
	double defaultvalue = spm->defaultvalue;

	pos_contents = spm->contents;
	for (iindex = 0; iindex < size; iindex++, pos_contents++)
		if ((pos_contents->state == MAPKIT_FULLSLOT)
		    && ((pos_contents->value) == (defaultvalue))) {
			pos_contents->state = MAPKIT_DELETEDSLOT;
			count++;
		}

	spm->used -= count;
	if (spm->used < spm->minused) {
#ifdef MAPKIT_DEBUG
		fprintf(stderr, "MAPKIT: used < minused\n");
#endif
		return map_id_reallocate(spm, map_id_meansize(spm, spm->used));
	}

	return MAPKIT_OK;
}

int map_id_compare(const void *e1, const void *e2)
{
	int key1 = ((const map_id_element *) e1)->key;
	int key2 = ((const map_id_element *) e2)->key;

	return ((key1) < (key2) ? -1 : ((key1) == (key2) ? 0 : 1));
}

mapkit_error map_id_getall_sorted(map_id * spm, map_id_element ** array, mapkit_size_t * count)
{
	mapkit_error err;

	err = map_id_getall(spm, array, count);
	if (err)
		MAPKIT_ERROR(err);

	qsort(*array, *count, sizeof(**array), map_id_compare);

	return MAPKIT_OK;
}

mapkit_error map_id_setall(map_id * spm, map_id_element * array, mapkit_size_t count)
{
	mapkit_size_t array_iindex;
	mapkit_error err;

	err = map_id_ensurecapacity(spm, spm->used + count);
	if (err)
		MAPKIT_ERROR(err);

	if (spm->alwaysdefault)
		/*
		 * Prevent shrinking 
		 */
		spm->minused = 0;

	for (array_iindex = 0; array_iindex < count; array_iindex++) {
		err = map_id_set(spm, array[array_iindex].key, array[array_iindex].value);
		if (err)
			MAPKIT_ERROR(err);
	}

	if (spm->alwaysdefault)
		map_id_adjustcapacity(spm);

	return MAPKIT_OK;
}

mapkit_error map_id_removeall(map_id * spm, int *array, mapkit_size_t count)
{
	mapkit_size_t array_iindex;
	mapkit_error err;

	/*
	 * Prevent shrinking 
	 */
	spm->minused = 0;

	for (array_iindex = 0; array_iindex < count; array_iindex++) {
		err = map_id_remove(spm, array[array_iindex]);
		if (err)
			MAPKIT_ERROR(err);
	}

	map_id_adjustcapacity(spm);

	return MAPKIT_OK;
}

void map_id_printstats(map_id * spm)
{
	fprintf(stderr, "MAPKIT: map_id statistics\n");
	fprintf(stderr, "MAPKIT: alwaysdefault = %d\n", spm->alwaysdefault);
	fprintf(stderr, "MAPKIT: minused = %ld, maxfill = %ld\n", (long) spm->minused, (long) spm->maxfill);
	fprintf(stderr, "MAPKIT: minusedfactor = %g, maxfillfactor = %g\n", spm->minusedfactor, spm->maxfillfactor);
#ifdef MAPKIT_COLLISIONS
	fprintf(stderr, "MAPKIT: insertionindexs = %lu, collisions = %lu\n", (unsigned long) spm->insertionindexs,
		(unsigned long) spm->insertionindex_collisions);
	fprintf(stderr, "MAPKIT: keyindexs = %lu, collisions = %lu\n", (unsigned long) spm->keyindexs, (unsigned long) spm->keyindex_collisions);
#endif
}

#endif							       /* MAPKIT_map_id */

#ifdef MAPKIT_map_ivp

/*
  Implementation for map_ivp (int -> void*)
  Default value : NULL
  Uses a state field.
*/

/* Static prototypes */

/* Return the iindex of key, or MAPKIT_KEYNOTFOUND */
static mapkit_size_t map_ivp_keyindex(map_ivp * spm, int key);

/* Return the iindex of key or -(insertion iindex)-1 if key not found */
static mapkit_size_t map_ivp_insertionindex(map_ivp * spm, int key);

/* Implementation */

mapkit_error map_ivp_init(map_ivp * spm)
{
	return map_ivp_init_hint(spm, MAPKIT_DEFAULT_EXPECTEDUSED);
}

mapkit_error map_ivp_init_hint(map_ivp * spm, mapkit_size_t used)
{
#ifdef MAPKIT_DEBUG
	fprintf(stderr, "MAPKIT: init\n");
#endif

	spm->size = 0;
	spm->fill = 0;
	spm->used = 0;
	spm->maxfillfactor = 0.5;
	spm->minusedfactor = 0.2;
	spm->contents = NULL;
#ifdef MAPKIT_COLLISIONS
	spm->insertionindexs = spm->insertionindex_collisions = 0;
	spm->keyindexs = spm->keyindex_collisions = 0;
#endif
	spm->defaultvalue = NULL;
	spm->alwaysdefault = 0;

	return map_ivp_reallocate(spm, map_ivp_meansize(spm, used));
}

mapkit_error map_ivp_ensurecapacity(map_ivp * spm, mapkit_size_t used)
{
	if (used > (spm->used + spm->maxfill - spm->fill)) {
#ifdef MAPKIT_DEBUG
		fprintf(stderr, "MAPKIT: need more capacity\n");
#endif
		return map_ivp_reallocate(spm, map_ivp_meansize(spm, used));
	} else
		return MAPKIT_OK;
}

mapkit_error map_ivp_adjustcapacity(map_ivp * spm)
{
	spm->minused = (mapkit_size_t) (spm->size * spm->minusedfactor);
	spm->maxfill = (mapkit_size_t) (spm->size * spm->maxfillfactor);

	if (spm->used < spm->minused) {
#ifdef MAPKIT_DEBUG
		fprintf(stderr, "MAPKIT: used < minused\n");
#endif
		return map_ivp_reallocate(spm, map_ivp_meansize(spm, spm->used));
	} else if (spm->fill > spm->maxfill) {
#ifdef MAPKIT_DEBUG
		fprintf(stderr, "MAPKIT: fill > maxfill\n");
#endif
		return map_ivp_reallocate(spm, map_ivp_meansize(spm, spm->used));
	} else
		return MAPKIT_OK;
}

void map_ivp_free(map_ivp * spm)
{
	free(spm->contents);
	spm->contents = NULL;
	spm->size = 0;
	spm->fill = 0;
	spm->used = 0;
	spm->maxfill = 0;
	spm->minused = 0;
#ifdef MAPKIT_COLLISIONS
	spm->insertionindexs = spm->insertionindex_collisions = 0;
	spm->keyindexs = spm->keyindex_collisions = 0;
#endif
}

mapkit_error map_ivp_copy(map_ivp * to, map_ivp * from)
{
	map_ivp_storage *contentscopy;
	size_t size = from->size * sizeof(*from->contents);

	contentscopy = (map_ivp_storage *) malloc(size);
	if (contentscopy == NULL)
		MAPKIT_ERROR(MAPKIT_ENOMEM);

	memcpy(to, from, sizeof(map_ivp));
	to->contents = contentscopy;
	memcpy(to->contents, from->contents, size);

	return MAPKIT_OK;
}

mapkit_size_t map_ivp_growsize(map_ivp * spm, mapkit_size_t used)
{
	return (mapkit_size_t) (4.0 * used / (3 * spm->minusedfactor + spm->maxfillfactor));
}

mapkit_size_t map_ivp_shrinksize(map_ivp * spm, mapkit_size_t used)
{
	return (mapkit_size_t) (4.0 * used / (spm->minusedfactor + 3 * spm->maxfillfactor));
}

mapkit_size_t map_ivp_meansize(map_ivp * spm, mapkit_size_t used)
{
	return (mapkit_size_t) (2.0 * used / (spm->minusedfactor + spm->maxfillfactor));
}

mapkit_error map_ivp_reallocate(map_ivp * spm, mapkit_size_t newsize)
{
	mapkit_size_t iindex;
	mapkit_size_t oldsize;
	map_ivp_storage *newcontents, *oldcontents;

	/*
	 * At least one free entry 
	 */
	if (newsize <= spm->used)
		newsize = spm->used + 1;
	newsize = mapkit_nextprime(newsize);
	if (newsize <= spm->used)
		MAPKIT_ERROR(MAPKIT_ETOOBIG);

	newcontents = (map_ivp_storage *) malloc(newsize * sizeof(*spm->contents));
	if (newcontents == NULL)
		MAPKIT_ERROR(MAPKIT_ENOMEM);

	/*
	 * Initialize all entries to "free" 
	 */
	for (iindex = 0; iindex < newsize; iindex++)
		newcontents[iindex].state = MAPKIT_FREESLOT;

	oldcontents = spm->contents;
	oldsize = spm->size;
	spm->contents = newcontents;
	spm->size = newsize;

#ifdef MAPKIT_DEBUG
	fprintf(stderr, "MAPKIT: reallocate %ld -> %ld\n", (long) oldsize, (long) newsize);
#endif

	spm->maxfill = (mapkit_size_t) (newsize * spm->maxfillfactor);
	/*
	 * At least one free entry 
	 */
	if (spm->maxfill >= newsize)
		spm->maxfill = newsize - 1;
	spm->minused = (mapkit_size_t) (newsize * spm->minusedfactor);
	spm->used = 0;

	if (oldcontents != NULL) {
		int used = 0;
		int key;
		void *defaultvalue = spm->defaultvalue;
		int notalwaysdefault = !spm->alwaysdefault;

		/*
		 * Copy all entries from old to new 
		 */
		for (iindex = 0; iindex < oldsize; iindex++)
			if (oldcontents[iindex].state == MAPKIT_FULLSLOT) {
				mapkit_size_t ins_iindex;
				map_ivp_storage *contents;

				key = oldcontents[iindex].key;

				/*
				 * Fast path 
				 */
				ins_iindex = ((mapkit_hash_t) key) % spm->size;
				contents = &(newcontents[ins_iindex]);

				if (contents->state != MAPKIT_FREESLOT) {
					ins_iindex = map_ivp_insertionindex(spm, key);
					contents = &(newcontents[ins_iindex]);
				}
#ifdef MAPKIT_COLLISIONS
				else
					spm->insertionindexs++;
#endif
				if (notalwaysdefault || (!((oldcontents[iindex].value) == (defaultvalue)))) {
					contents->value = oldcontents[iindex].value;
					contents->state = MAPKIT_FULLSLOT;
					contents->key = key;
					used++;
				}
			}
		free(oldcontents);
		spm->used = used;
	}
	spm->fill = spm->used;

	return MAPKIT_OK;
}

void *map_ivp_value_s(map_ivp * spm, int key)
{
	mapkit_size_t iindex;

	iindex = map_ivp_keyindex(spm, key);

	if (iindex < 0) {
		if (spm->alwaysdefault)
			return spm->defaultvalue;
		else
			MAPKIT_FATAL_ERROR(MAPKIT_EKEYNOTFOUND);
	}

	return spm->contents[iindex].value;
}

mapkit_error map_ivp_get_s(map_ivp * spm, int key, void **value)
{
	mapkit_size_t iindex;

	iindex = map_ivp_keyindex(spm, key);

	if (iindex < 0) {
		if (spm->alwaysdefault) {
			*value = spm->defaultvalue;
			return MAPKIT_OK;
		} else
			MAPKIT_ERROR(MAPKIT_EKEYNOTFOUND);
	}

	*value = spm->contents[iindex].value;
	return MAPKIT_OK;
}

mapkit_error map_ivp_set_s(map_ivp * spm, int key, void *value)
{
	mapkit_size_t iindex;

	iindex = map_ivp_insertionindex(spm, key);

	if (iindex < 0)
		/*
		 * FULLSLOT 
		 */
		spm->contents[-iindex - 1].value = value;
	else {
		map_ivp_storage *element = &(spm->contents[iindex]);
		int ffree = element->state == MAPKIT_FREESLOT;

		element->state = MAPKIT_FULLSLOT;
		element->value = value;
		element->key = key;
		spm->used++;

		if (ffree && ((++spm->fill) > spm->maxfill)) {
#ifdef MAPKIT_DEBUG
			fprintf(stderr, "MAPKIT: fill > maxfill\n");
#endif
			return map_ivp_reallocate(spm, map_ivp_growsize(spm, spm->used));
		}
	}
	return MAPKIT_OK;
}

void **map_ivp_insertptr_s(map_ivp * spm, int key)
{
	mapkit_size_t iindex;

	iindex = map_ivp_insertionindex(spm, key);

	if (iindex < 0)
		return &(spm->contents[-iindex - 1].value);
	else {
		map_ivp_storage *element = &(spm->contents[iindex]);

		if (element->state == MAPKIT_FREESLOT) {
			/*
			 * FREESLOT 
			 */
			if (spm->fill >= spm->maxfill) {
				mapkit_error err;

#ifdef MAPKIT_DEBUG
				fprintf(stderr, "MAPKIT: fill => maxfill before insert\n");
#endif
				/*
				 * Must reallocate -before- inserting defaultvalue 
				 */
				err = map_ivp_reallocate(spm, map_ivp_growsize(spm, spm->used + 1));
				if (err) {
					MAPKIT_ERROR_NORET(err);
					return NULL;
				}

				iindex = map_ivp_insertionindex(spm, key);
				/*
				 * FREESLOT 
				 */
				spm->contents[iindex].state = MAPKIT_FULLSLOT;
				spm->contents[iindex].key = key;
				spm->contents[iindex].value = spm->defaultvalue;
				spm->used++;
				spm->fill++;
				return &(spm->contents[iindex].value);
			} else
				spm->fill++;
		}

		element->state = MAPKIT_FULLSLOT;
		element->key = key;
		element->value = spm->defaultvalue;
		spm->used++;
		return &(element->value);
	}
}

void **map_ivp_ptr_s(map_ivp * spm, int key)
{
	mapkit_size_t iindex;

	iindex = map_ivp_keyindex(spm, key);

	if ((iindex >= 0) && ((!spm->alwaysdefault) || (!((spm->contents[iindex].value) == (spm->defaultvalue)))))
		return &(spm->contents[iindex].value);
	else
		return NULL;
}

mapkit_error map_ivp_remove_s(map_ivp * spm, int key)
{
	mapkit_size_t iindex;

	iindex = map_ivp_keyindex(spm, key);

	if (iindex < 0) {
		if (spm->alwaysdefault)
			return MAPKIT_OK;
		else
			return MAPKIT_EKEYNOTFOUND;
	}

	spm->contents[iindex].state = MAPKIT_DELETEDSLOT;
	spm->used--;
	if (spm->used < spm->minused) {
#ifdef MAPKIT_DEBUG
		fprintf(stderr, "MAPKIT: used < minused\n");
#endif
		return map_ivp_reallocate(spm, map_ivp_shrinksize(spm, spm->used));
	}

	return MAPKIT_OK;
}

mapkit_size_t map_ivp_keyindex(map_ivp * spm, int key)
{
	mapkit_size_t iindex, decrement;

	signed char state;

	iindex = ((mapkit_hash_t) key) % spm->size;
	decrement = (((mapkit_hash_t) key) % (spm->size - 2));
	decrement += (decrement == 0);

#ifdef MAPKIT_COLLISIONS
	spm->keyindexs++;
#endif

	while ((state = spm->contents[iindex].state) != MAPKIT_FREESLOT
	       && (state == MAPKIT_DELETEDSLOT || (!((spm->contents[iindex].key) == (key))))) {
#ifdef MAPKIT_COLLISIONS
		spm->keyindex_collisions++;
#endif
		iindex -= decrement;
		if (iindex < 0)
			iindex += spm->size;
	}

	if (state == MAPKIT_FREESLOT)
		return MAPKIT_KEYNOTFOUND;
	return iindex;
}

mapkit_size_t map_ivp_insertionindex(map_ivp * spm, int key)
{
	mapkit_size_t iindex, decrement;
	signed char state;

#ifdef MAPKIT_COLLISIONS
	spm->insertionindexs++;
#endif

	iindex = ((mapkit_hash_t) key) % spm->size;

	/*
	 * Fast path (largely superfluous) 
	 */
	if ((state = spm->contents[iindex].state) == MAPKIT_FREESLOT)
		return iindex;
	if ((state == MAPKIT_FULLSLOT)
	    && ((spm->contents[iindex].key) == (key)))
		return -iindex - 1;

	decrement = (((mapkit_hash_t) key) % (spm->size - 2));
	decrement += (decrement == 0);

	while ((state == MAPKIT_FULLSLOT)
	       && (!((spm->contents[iindex].key) == (key)))) {
#ifdef MAPKIT_COLLISIONS
		spm->insertionindex_collisions++;
#endif
		iindex -= decrement;
		if (iindex < 0)
			iindex += spm->size;
		state = spm->contents[iindex].state;
	}

	if (state == MAPKIT_DELETEDSLOT) {
		mapkit_size_t index2 = iindex;

		while ((state = spm->contents[iindex].state) != MAPKIT_FREESLOT && ((state == MAPKIT_DELETEDSLOT)
										    || (!((spm->contents[iindex].key) == (key))))) {
			iindex -= decrement;
			if (iindex < 0)
				iindex += spm->size;
			state = spm->contents[iindex].state;
		}
		if (state == MAPKIT_FREESLOT)
			return index2;
	}

	if (state == MAPKIT_FULLSLOT)
		return -iindex - 1;
	return iindex;
}

mapkit_size_t map_ivp_next(map_ivp * spm, mapkit_size_t iindex)
{
	mapkit_size_t size = spm->size;
	map_ivp_storage *pos_contents;
	void *defaultvalue = spm->defaultvalue;
	int notalwaysdefault = !spm->alwaysdefault;

	pos_contents = &(spm->contents[++iindex]);

	for (; iindex < size; iindex++, pos_contents++)
		if ((pos_contents->state == MAPKIT_FULLSLOT)
		    && (notalwaysdefault || (!((pos_contents->value) == (defaultvalue)))))
			return iindex;

	return -1;
}

map_ivp_storage *map_ivp_nextptr(map_ivp * spm, map_ivp_storage * pos_contents)
{
	map_ivp_storage *end = &(spm->contents[spm->size]);
	void *defaultvalue = spm->defaultvalue;
	int notalwaysdefault = !spm->alwaysdefault;

	if (pos_contents == NULL)
		pos_contents = spm->contents;
	else {
		pos_contents++;
		if (pos_contents <= spm->contents)
			return NULL;
	}

	for (; pos_contents < end; pos_contents++)
		if ((pos_contents->state == MAPKIT_FULLSLOT)
		    && (notalwaysdefault || (!((pos_contents->value) == (defaultvalue)))))
			return pos_contents;

	return NULL;
}

mapkit_error map_ivp_getall(map_ivp * spm, map_ivp_element ** array, mapkit_size_t * count)
{
	mapkit_size_t iindex;
	mapkit_size_t size = spm->size, vcount = 0;
	map_ivp_element *pos_array;
	map_ivp_storage *pos_contents;
	void *defaultvalue = spm->defaultvalue;
	int notalwaysdefault = !spm->alwaysdefault;

	pos_array = *array = (map_ivp_element *) malloc(sizeof(**array) * spm->used);
	if (*array == NULL)
		MAPKIT_ERROR(MAPKIT_ENOMEM);

	pos_contents = spm->contents;

	for (iindex = 0; iindex < size; iindex++, pos_contents++)
		if ((pos_contents->state == MAPKIT_FULLSLOT)
		    && (notalwaysdefault || (!((pos_contents->value) == (defaultvalue))))) {
			pos_array->key = pos_contents->key;
			pos_array->value = pos_contents->value;
			pos_array++;
			vcount++;
		}
	*count = vcount;

	return MAPKIT_OK;
}

mapkit_error map_ivp_clean(map_ivp * spm)
{
	mapkit_size_t iindex, count = 0;
	mapkit_size_t size = spm->size;
	map_ivp_storage *pos_contents;
	void *defaultvalue = spm->defaultvalue;

	pos_contents = spm->contents;
	for (iindex = 0; iindex < size; iindex++, pos_contents++)
		if ((pos_contents->state == MAPKIT_FULLSLOT)
		    && ((pos_contents->value) == (defaultvalue))) {
			pos_contents->state = MAPKIT_DELETEDSLOT;
			count++;
		}

	spm->used -= count;
	if (spm->used < spm->minused) {
#ifdef MAPKIT_DEBUG
		fprintf(stderr, "MAPKIT: used < minused\n");
#endif
		return map_ivp_reallocate(spm, map_ivp_meansize(spm, spm->used));
	}

	return MAPKIT_OK;
}

int map_ivp_compare(const void *e1, const void *e2)
{
	int key1 = ((const map_ivp_element *) e1)->key;
	int key2 = ((const map_ivp_element *) e2)->key;

	return ((key1) < (key2) ? -1 : ((key1) == (key2) ? 0 : 1));
}

mapkit_error map_ivp_getall_sorted(map_ivp * spm, map_ivp_element ** array, mapkit_size_t * count)
{
	mapkit_error err;

	err = map_ivp_getall(spm, array, count);
	if (err)
		MAPKIT_ERROR(err);

	qsort(*array, *count, sizeof(**array), map_ivp_compare);

	return MAPKIT_OK;
}

mapkit_error map_ivp_setall(map_ivp * spm, map_ivp_element * array, mapkit_size_t count)
{
	mapkit_size_t array_iindex;
	mapkit_error err;

	err = map_ivp_ensurecapacity(spm, spm->used + count);
	if (err)
		MAPKIT_ERROR(err);

	if (spm->alwaysdefault)
		/*
		 * Prevent shrinking 
		 */
		spm->minused = 0;

	for (array_iindex = 0; array_iindex < count; array_iindex++) {
		err = map_ivp_set(spm, array[array_iindex].key, array[array_iindex].value);
		if (err)
			MAPKIT_ERROR(err);
	}

	if (spm->alwaysdefault)
		map_ivp_adjustcapacity(spm);

	return MAPKIT_OK;
}

mapkit_error map_ivp_removeall(map_ivp * spm, int *array, mapkit_size_t count)
{
	mapkit_size_t array_iindex;
	mapkit_error err;

	/*
	 * Prevent shrinking 
	 */
	spm->minused = 0;

	for (array_iindex = 0; array_iindex < count; array_iindex++) {
		err = map_ivp_remove(spm, array[array_iindex]);
		if (err)
			MAPKIT_ERROR(err);
	}

	map_ivp_adjustcapacity(spm);

	return MAPKIT_OK;
}

void map_ivp_printstats(map_ivp * spm)
{
	fprintf(stderr, "MAPKIT: map_ivp statistics\n");
	fprintf(stderr, "MAPKIT: alwaysdefault = %d\n", spm->alwaysdefault);
	fprintf(stderr, "MAPKIT: minused = %ld, maxfill = %ld\n", (long) spm->minused, (long) spm->maxfill);
	fprintf(stderr, "MAPKIT: minusedfactor = %g, maxfillfactor = %g\n", spm->minusedfactor, spm->maxfillfactor);
#ifdef MAPKIT_COLLISIONS
	fprintf(stderr, "MAPKIT: insertionindexs = %lu, collisions = %lu\n", (unsigned long) spm->insertionindexs,
		(unsigned long) spm->insertionindex_collisions);
	fprintf(stderr, "MAPKIT: keyindexs = %lu, collisions = %lu\n", (unsigned long) spm->keyindexs, (unsigned long) spm->keyindex_collisions);
#endif
}

#endif							       /* MAPKIT_map_ivp */

#ifdef MAPKIT_map_h_ii

/*
  Implementation for map_h_ii (int -> int)
  Default value : 0
  Uses a state field.
*/

/* Static prototypes */

/* Return the iindex of key, or MAPKIT_KEYNOTFOUND */
static mapkit_size_t map_h_ii_keyindex(map_h_ii * spm, int key, mapkit_hash_t hash);

/* Return the iindex of key or -(insertion iindex)-1 if key not found */
static mapkit_size_t map_h_ii_insertionindex(map_h_ii * spm, int key, mapkit_hash_t hash);

/* Implementation */

mapkit_error map_h_ii_init(map_h_ii * spm)
{
	return map_h_ii_init_hint(spm, MAPKIT_DEFAULT_EXPECTEDUSED);
}

mapkit_error map_h_ii_init_hint(map_h_ii * spm, mapkit_size_t used)
{
#ifdef MAPKIT_DEBUG
	fprintf(stderr, "MAPKIT: init\n");
#endif

	spm->size = 0;
	spm->fill = 0;
	spm->used = 0;
	spm->maxfillfactor = 0.5;
	spm->minusedfactor = 0.2;
	spm->contents = NULL;
#ifdef MAPKIT_COLLISIONS
	spm->insertionindexs = spm->insertionindex_collisions = 0;
	spm->keyindexs = spm->keyindex_collisions = 0;
#endif
	spm->defaultvalue = 0;
	spm->alwaysdefault = 0;

	return map_h_ii_reallocate(spm, map_h_ii_meansize(spm, used));
}

mapkit_error map_h_ii_ensurecapacity(map_h_ii * spm, mapkit_size_t used)
{
	if (used > (spm->used + spm->maxfill - spm->fill)) {
#ifdef MAPKIT_DEBUG
		fprintf(stderr, "MAPKIT: need more capacity\n");
#endif
		return map_h_ii_reallocate(spm, map_h_ii_meansize(spm, used));
	} else
		return MAPKIT_OK;
}

mapkit_error map_h_ii_adjustcapacity(map_h_ii * spm)
{
	spm->minused = (mapkit_size_t) (spm->size * spm->minusedfactor);
	spm->maxfill = (mapkit_size_t) (spm->size * spm->maxfillfactor);

	if (spm->used < spm->minused) {
#ifdef MAPKIT_DEBUG
		fprintf(stderr, "MAPKIT: used < minused\n");
#endif
		return map_h_ii_reallocate(spm, map_h_ii_meansize(spm, spm->used));
	} else if (spm->fill > spm->maxfill) {
#ifdef MAPKIT_DEBUG
		fprintf(stderr, "MAPKIT: fill > maxfill\n");
#endif
		return map_h_ii_reallocate(spm, map_h_ii_meansize(spm, spm->used));
	} else
		return MAPKIT_OK;
}

void map_h_ii_free(map_h_ii * spm)
{
	free(spm->contents);
	spm->contents = NULL;
	spm->size = 0;
	spm->fill = 0;
	spm->used = 0;
	spm->maxfill = 0;
	spm->minused = 0;
#ifdef MAPKIT_COLLISIONS
	spm->insertionindexs = spm->insertionindex_collisions = 0;
	spm->keyindexs = spm->keyindex_collisions = 0;
#endif
}

mapkit_error map_h_ii_copy(map_h_ii * to, map_h_ii * from)
{
	map_h_ii_storage *contentscopy;
	size_t size = from->size * sizeof(*from->contents);

	contentscopy = (map_h_ii_storage *) malloc(size);
	if (contentscopy == NULL)
		MAPKIT_ERROR(MAPKIT_ENOMEM);

	memcpy(to, from, sizeof(map_h_ii));
	to->contents = contentscopy;
	memcpy(to->contents, from->contents, size);

	return MAPKIT_OK;
}

mapkit_size_t map_h_ii_growsize(map_h_ii * spm, mapkit_size_t used)
{
	return (mapkit_size_t) (4.0 * used / (3 * spm->minusedfactor + spm->maxfillfactor));
}

mapkit_size_t map_h_ii_shrinksize(map_h_ii * spm, mapkit_size_t used)
{
	return (mapkit_size_t) (4.0 * used / (spm->minusedfactor + 3 * spm->maxfillfactor));
}

mapkit_size_t map_h_ii_meansize(map_h_ii * spm, mapkit_size_t used)
{
	return (mapkit_size_t) (2.0 * used / (spm->minusedfactor + spm->maxfillfactor));
}

mapkit_error map_h_ii_reallocate(map_h_ii * spm, mapkit_size_t newsize)
{
	mapkit_size_t iindex;
	mapkit_size_t oldsize;
	map_h_ii_storage *newcontents, *oldcontents;

	/*
	 * At least one free entry 
	 */
	if (newsize <= spm->used)
		newsize = spm->used + 1;
	newsize = mapkit_nextprime(newsize);
	if (newsize <= spm->used)
		MAPKIT_ERROR(MAPKIT_ETOOBIG);

	newcontents = (map_h_ii_storage *) malloc(newsize * sizeof(*spm->contents));
	if (newcontents == NULL)
		MAPKIT_ERROR(MAPKIT_ENOMEM);

	/*
	 * Initialize all entries to "free" 
	 */
	for (iindex = 0; iindex < newsize; iindex++)
		newcontents[iindex].state = MAPKIT_FREESLOT;

	oldcontents = spm->contents;
	oldsize = spm->size;
	spm->contents = newcontents;
	spm->size = newsize;

#ifdef MAPKIT_DEBUG
	fprintf(stderr, "MAPKIT: reallocate %ld -> %ld\n", (long) oldsize, (long) newsize);
#endif

	spm->maxfill = (mapkit_size_t) (newsize * spm->maxfillfactor);
	/*
	 * At least one free entry 
	 */
	if (spm->maxfill >= newsize)
		spm->maxfill = newsize - 1;
	spm->minused = (mapkit_size_t) (newsize * spm->minusedfactor);
	spm->used = 0;

	if (oldcontents != NULL) {
		int used = 0;
		int key;
		int defaultvalue = spm->defaultvalue;
		int notalwaysdefault = !spm->alwaysdefault;

		/*
		 * Copy all entries from old to new 
		 */
		for (iindex = 0; iindex < oldsize; iindex++)
			if (oldcontents[iindex].state == MAPKIT_FULLSLOT) {
				mapkit_size_t ins_iindex;
				mapkit_hash_t hash;
				map_h_ii_storage *contents;

				key = oldcontents[iindex].key;

				/*
				 * Fast path 
				 */
				ins_iindex = (hash = mapkit_hash((mapkit_hash_t) key)) % spm->size;
				contents = &(newcontents[ins_iindex]);

				if (contents->state != MAPKIT_FREESLOT) {
					ins_iindex = map_h_ii_insertionindex(spm, key, hash);
					contents = &(newcontents[ins_iindex]);
				}
#ifdef MAPKIT_COLLISIONS
				else
					spm->insertionindexs++;
#endif
				if (notalwaysdefault || (!((oldcontents[iindex].value) == (defaultvalue)))) {
					contents->value = oldcontents[iindex].value;
					contents->state = MAPKIT_FULLSLOT;
					contents->key = key;
					used++;
				}
			}
		free(oldcontents);
		spm->used = used;
	}
	spm->fill = spm->used;

	return MAPKIT_OK;
}

int map_h_ii_value_s(map_h_ii * spm, int key, mapkit_hash_t hash)
{
	mapkit_size_t iindex;

	iindex = map_h_ii_keyindex(spm, key, hash);

	if (iindex < 0) {
		if (spm->alwaysdefault)
			return spm->defaultvalue;
		else
			MAPKIT_FATAL_ERROR(MAPKIT_EKEYNOTFOUND);
	}

	return spm->contents[iindex].value;
}

mapkit_error map_h_ii_get_s(map_h_ii * spm, int key, int *value, mapkit_hash_t hash)
{
	mapkit_size_t iindex;

	iindex = map_h_ii_keyindex(spm, key, hash);

	if (iindex < 0) {
		if (spm->alwaysdefault) {
			*value = spm->defaultvalue;
			return MAPKIT_OK;
		} else
			MAPKIT_ERROR(MAPKIT_EKEYNOTFOUND);
	}

	*value = spm->contents[iindex].value;
	return MAPKIT_OK;
}

mapkit_error map_h_ii_set_s(map_h_ii * spm, int key, int value, mapkit_hash_t hash)
{
	mapkit_size_t iindex;

	iindex = map_h_ii_insertionindex(spm, key, hash);

	if (iindex < 0)
		/*
		 * FULLSLOT 
		 */
		spm->contents[-iindex - 1].value = value;
	else {
		map_h_ii_storage *element = &(spm->contents[iindex]);
		int ffree = element->state == MAPKIT_FREESLOT;

		element->state = MAPKIT_FULLSLOT;
		element->value = value;
		element->key = key;
		spm->used++;

		if (ffree && ((++spm->fill) > spm->maxfill)) {
#ifdef MAPKIT_DEBUG
			fprintf(stderr, "MAPKIT: fill > maxfill\n");
#endif
			return map_h_ii_reallocate(spm, map_h_ii_growsize(spm, spm->used));
		}
	}
	return MAPKIT_OK;
}

int *map_h_ii_insertptr_s(map_h_ii * spm, int key, mapkit_hash_t hash)
{
	mapkit_size_t iindex;

	iindex = map_h_ii_insertionindex(spm, key, hash);

	if (iindex < 0)
		return &(spm->contents[-iindex - 1].value);
	else {
		map_h_ii_storage *element = &(spm->contents[iindex]);

		if (element->state == MAPKIT_FREESLOT) {
			/*
			 * FREESLOT 
			 */
			if (spm->fill >= spm->maxfill) {
				mapkit_error err;

#ifdef MAPKIT_DEBUG
				fprintf(stderr, "MAPKIT: fill => maxfill before insert\n");
#endif
				/*
				 * Must reallocate -before- inserting defaultvalue 
				 */
				err = map_h_ii_reallocate(spm, map_h_ii_growsize(spm, spm->used + 1));
				if (err) {
					MAPKIT_ERROR_NORET(err);
					return NULL;
				}

				iindex = map_h_ii_insertionindex(spm, key, hash);
				/*
				 * FREESLOT 
				 */
				spm->contents[iindex].state = MAPKIT_FULLSLOT;
				spm->contents[iindex].key = key;
				spm->contents[iindex].value = spm->defaultvalue;
				spm->used++;
				spm->fill++;
				return &(spm->contents[iindex].value);
			} else
				spm->fill++;
		}

		element->state = MAPKIT_FULLSLOT;
		element->key = key;
		element->value = spm->defaultvalue;
		spm->used++;
		return &(element->value);
	}
}

int *map_h_ii_ptr_s(map_h_ii * spm, int key, mapkit_hash_t hash)
{
	mapkit_size_t iindex;

	iindex = map_h_ii_keyindex(spm, key, hash);

	if ((iindex >= 0) && ((!spm->alwaysdefault) || (!((spm->contents[iindex].value) == (spm->defaultvalue)))))
		return &(spm->contents[iindex].value);
	else
		return NULL;
}

mapkit_error map_h_ii_remove_s(map_h_ii * spm, int key, mapkit_hash_t hash)
{
	mapkit_size_t iindex;

	iindex = map_h_ii_keyindex(spm, key, hash);

	if (iindex < 0) {
		if (spm->alwaysdefault)
			return MAPKIT_OK;
		else
			return MAPKIT_EKEYNOTFOUND;
	}

	spm->contents[iindex].state = MAPKIT_DELETEDSLOT;
	spm->used--;
	if (spm->used < spm->minused) {
#ifdef MAPKIT_DEBUG
		fprintf(stderr, "MAPKIT: used < minused\n");
#endif
		return map_h_ii_reallocate(spm, map_h_ii_shrinksize(spm, spm->used));
	}

	return MAPKIT_OK;
}

mapkit_size_t map_h_ii_keyindex(map_h_ii * spm, int key, mapkit_hash_t hash)
{
	mapkit_size_t iindex, decrement;

	signed char state;

	iindex = hash % spm->size;
	decrement = (hash % (spm->size - 2));
	decrement += (decrement == 0);

#ifdef MAPKIT_COLLISIONS
	spm->keyindexs++;
#endif

	while ((state = spm->contents[iindex].state) != MAPKIT_FREESLOT
	       && (state == MAPKIT_DELETEDSLOT || (!((spm->contents[iindex].key) == (key))))) {
#ifdef MAPKIT_COLLISIONS
		spm->keyindex_collisions++;
#endif
		iindex -= decrement;
		if (iindex < 0)
			iindex += spm->size;
	}

	if (state == MAPKIT_FREESLOT)
		return MAPKIT_KEYNOTFOUND;
	return iindex;
}

mapkit_size_t map_h_ii_insertionindex(map_h_ii * spm, int key, mapkit_hash_t hash)
{
	mapkit_size_t iindex, decrement;
	signed char state;

#ifdef MAPKIT_COLLISIONS
	spm->insertionindexs++;
#endif

	iindex = hash % spm->size;

	/*
	 * Fast path (largely superfluous) 
	 */
	if ((state = spm->contents[iindex].state) == MAPKIT_FREESLOT)
		return iindex;
	if ((state == MAPKIT_FULLSLOT)
	    && ((spm->contents[iindex].key) == (key)))
		return -iindex - 1;

	decrement = (hash % (spm->size - 2));
	decrement += (decrement == 0);

	while ((state == MAPKIT_FULLSLOT)
	       && (!((spm->contents[iindex].key) == (key)))) {
#ifdef MAPKIT_COLLISIONS
		spm->insertionindex_collisions++;
#endif
		iindex -= decrement;
		if (iindex < 0)
			iindex += spm->size;
		state = spm->contents[iindex].state;
	}

	if (state == MAPKIT_DELETEDSLOT) {
		mapkit_size_t index2 = iindex;

		while ((state = spm->contents[iindex].state) != MAPKIT_FREESLOT && ((state == MAPKIT_DELETEDSLOT)
										    || (!((spm->contents[iindex].key) == (key))))) {
			iindex -= decrement;
			if (iindex < 0)
				iindex += spm->size;
			state = spm->contents[iindex].state;
		}
		if (state == MAPKIT_FREESLOT)
			return index2;
	}

	if (state == MAPKIT_FULLSLOT)
		return -iindex - 1;
	return iindex;
}

mapkit_size_t map_h_ii_next(map_h_ii * spm, mapkit_size_t iindex)
{
	mapkit_size_t size = spm->size;
	map_h_ii_storage *pos_contents;
	int defaultvalue = spm->defaultvalue;
	int notalwaysdefault = !spm->alwaysdefault;

	pos_contents = &(spm->contents[++iindex]);

	for (; iindex < size; iindex++, pos_contents++)
		if ((pos_contents->state == MAPKIT_FULLSLOT)
		    && (notalwaysdefault || (!((pos_contents->value) == (defaultvalue)))))
			return iindex;

	return -1;
}

map_h_ii_storage *map_h_ii_nextptr(map_h_ii * spm, map_h_ii_storage * pos_contents)
{
	map_h_ii_storage *end = &(spm->contents[spm->size]);
	int defaultvalue = spm->defaultvalue;
	int notalwaysdefault = !spm->alwaysdefault;

	if (pos_contents == NULL)
		pos_contents = spm->contents;
	else {
		pos_contents++;
		if (pos_contents <= spm->contents)
			return NULL;
	}

	for (; pos_contents < end; pos_contents++)
		if ((pos_contents->state == MAPKIT_FULLSLOT)
		    && (notalwaysdefault || (!((pos_contents->value) == (defaultvalue)))))
			return pos_contents;

	return NULL;
}

mapkit_error map_h_ii_getall(map_h_ii * spm, map_h_ii_element ** array, mapkit_size_t * count)
{
	mapkit_size_t iindex;
	mapkit_size_t size = spm->size, vcount = 0;
	map_h_ii_element *pos_array;
	map_h_ii_storage *pos_contents;
	int defaultvalue = spm->defaultvalue;
	int notalwaysdefault = !spm->alwaysdefault;

	pos_array = *array = (map_h_ii_element *) malloc(sizeof(**array) * spm->used);
	if (*array == NULL)
		MAPKIT_ERROR(MAPKIT_ENOMEM);

	pos_contents = spm->contents;

	for (iindex = 0; iindex < size; iindex++, pos_contents++)
		if ((pos_contents->state == MAPKIT_FULLSLOT)
		    && (notalwaysdefault || (!((pos_contents->value) == (defaultvalue))))) {
			pos_array->key = pos_contents->key;
			pos_array->value = pos_contents->value;
			pos_array++;
			vcount++;
		}
	*count = vcount;

	return MAPKIT_OK;
}

mapkit_error map_h_ii_clean(map_h_ii * spm)
{
	mapkit_size_t iindex, count = 0;
	mapkit_size_t size = spm->size;
	map_h_ii_storage *pos_contents;
	int defaultvalue = spm->defaultvalue;

	pos_contents = spm->contents;
	for (iindex = 0; iindex < size; iindex++, pos_contents++)
		if ((pos_contents->state == MAPKIT_FULLSLOT)
		    && ((pos_contents->value) == (defaultvalue))) {
			pos_contents->state = MAPKIT_DELETEDSLOT;
			count++;
		}

	spm->used -= count;
	if (spm->used < spm->minused) {
#ifdef MAPKIT_DEBUG
		fprintf(stderr, "MAPKIT: used < minused\n");
#endif
		return map_h_ii_reallocate(spm, map_h_ii_meansize(spm, spm->used));
	}

	return MAPKIT_OK;
}

int map_h_ii_compare(const void *e1, const void *e2)
{
	int key1 = ((const map_h_ii_element *) e1)->key;
	int key2 = ((const map_h_ii_element *) e2)->key;

	return ((key1) < (key2) ? -1 : ((key1) == (key2) ? 0 : 1));
}

mapkit_error map_h_ii_getall_sorted(map_h_ii * spm, map_h_ii_element ** array, mapkit_size_t * count)
{
	mapkit_error err;

	err = map_h_ii_getall(spm, array, count);
	if (err)
		MAPKIT_ERROR(err);

	qsort(*array, *count, sizeof(**array), map_h_ii_compare);

	return MAPKIT_OK;
}

mapkit_error map_h_ii_setall(map_h_ii * spm, map_h_ii_element * array, mapkit_size_t count)
{
	mapkit_size_t array_iindex;
	mapkit_error err;

	err = map_h_ii_ensurecapacity(spm, spm->used + count);
	if (err)
		MAPKIT_ERROR(err);

	if (spm->alwaysdefault)
		/*
		 * Prevent shrinking 
		 */
		spm->minused = 0;

	for (array_iindex = 0; array_iindex < count; array_iindex++) {
		err = map_h_ii_set(spm, array[array_iindex].key, array[array_iindex].value);
		if (err)
			MAPKIT_ERROR(err);
	}

	if (spm->alwaysdefault)
		map_h_ii_adjustcapacity(spm);

	return MAPKIT_OK;
}

mapkit_error map_h_ii_removeall(map_h_ii * spm, int *array, mapkit_size_t count)
{
	mapkit_size_t array_iindex;
	mapkit_error err;

	/*
	 * Prevent shrinking 
	 */
	spm->minused = 0;

	for (array_iindex = 0; array_iindex < count; array_iindex++) {
		err = map_h_ii_remove(spm, array[array_iindex]);
		if (err)
			MAPKIT_ERROR(err);
	}

	map_h_ii_adjustcapacity(spm);

	return MAPKIT_OK;
}

void map_h_ii_printstats(map_h_ii * spm)
{
	fprintf(stderr, "MAPKIT: map_h_ii statistics\n");
	fprintf(stderr, "MAPKIT: alwaysdefault = %d\n", spm->alwaysdefault);
	fprintf(stderr, "MAPKIT: minused = %ld, maxfill = %ld\n", (long) spm->minused, (long) spm->maxfill);
	fprintf(stderr, "MAPKIT: minusedfactor = %g, maxfillfactor = %g\n", spm->minusedfactor, spm->maxfillfactor);
#ifdef MAPKIT_COLLISIONS
	fprintf(stderr, "MAPKIT: insertionindexs = %lu, collisions = %lu\n", (unsigned long) spm->insertionindexs,
		(unsigned long) spm->insertionindex_collisions);
	fprintf(stderr, "MAPKIT: keyindexs = %lu, collisions = %lu\n", (unsigned long) spm->keyindexs, (unsigned long) spm->keyindex_collisions);
#endif
}

#endif							       /* MAPKIT_map_h_ii */

#ifdef MAPKIT_map_h_id

/*
  Implementation for map_h_id (int -> double)
  Default value : 0.0
  Uses a state field.
*/

/* Static prototypes */

/* Return the iindex of key, or MAPKIT_KEYNOTFOUND */
static mapkit_size_t map_h_id_keyindex(map_h_id * spm, int key, mapkit_hash_t hash);

/* Return the iindex of key or -(insertion iindex)-1 if key not found */
static mapkit_size_t map_h_id_insertionindex(map_h_id * spm, int key, mapkit_hash_t hash);

/* Implementation */

mapkit_error map_h_id_init(map_h_id * spm)
{
	return map_h_id_init_hint(spm, MAPKIT_DEFAULT_EXPECTEDUSED);
}

mapkit_error map_h_id_init_hint(map_h_id * spm, mapkit_size_t used)
{
#ifdef MAPKIT_DEBUG
	fprintf(stderr, "MAPKIT: init\n");
#endif

	spm->size = 0;
	spm->fill = 0;
	spm->used = 0;
	spm->maxfillfactor = 0.5;
	spm->minusedfactor = 0.2;
	spm->contents = NULL;
#ifdef MAPKIT_COLLISIONS
	spm->insertionindexs = spm->insertionindex_collisions = 0;
	spm->keyindexs = spm->keyindex_collisions = 0;
#endif
	spm->defaultvalue = 0.0;
	spm->alwaysdefault = 0;

	return map_h_id_reallocate(spm, map_h_id_meansize(spm, used));
}

mapkit_error map_h_id_ensurecapacity(map_h_id * spm, mapkit_size_t used)
{
	if (used > (spm->used + spm->maxfill - spm->fill)) {
#ifdef MAPKIT_DEBUG
		fprintf(stderr, "MAPKIT: need more capacity\n");
#endif
		return map_h_id_reallocate(spm, map_h_id_meansize(spm, used));
	} else
		return MAPKIT_OK;
}

mapkit_error map_h_id_adjustcapacity(map_h_id * spm)
{
	spm->minused = (mapkit_size_t) (spm->size * spm->minusedfactor);
	spm->maxfill = (mapkit_size_t) (spm->size * spm->maxfillfactor);

	if (spm->used < spm->minused) {
#ifdef MAPKIT_DEBUG
		fprintf(stderr, "MAPKIT: used < minused\n");
#endif
		return map_h_id_reallocate(spm, map_h_id_meansize(spm, spm->used));
	} else if (spm->fill > spm->maxfill) {
#ifdef MAPKIT_DEBUG
		fprintf(stderr, "MAPKIT: fill > maxfill\n");
#endif
		return map_h_id_reallocate(spm, map_h_id_meansize(spm, spm->used));
	} else
		return MAPKIT_OK;
}

void map_h_id_free(map_h_id * spm)
{
	free(spm->contents);
	spm->contents = NULL;
	spm->size = 0;
	spm->fill = 0;
	spm->used = 0;
	spm->maxfill = 0;
	spm->minused = 0;
#ifdef MAPKIT_COLLISIONS
	spm->insertionindexs = spm->insertionindex_collisions = 0;
	spm->keyindexs = spm->keyindex_collisions = 0;
#endif
}

mapkit_error map_h_id_copy(map_h_id * to, map_h_id * from)
{
	map_h_id_storage *contentscopy;
	size_t size = from->size * sizeof(*from->contents);

	contentscopy = (map_h_id_storage *) malloc(size);
	if (contentscopy == NULL)
		MAPKIT_ERROR(MAPKIT_ENOMEM);

	memcpy(to, from, sizeof(map_h_id));
	to->contents = contentscopy;
	memcpy(to->contents, from->contents, size);

	return MAPKIT_OK;
}

mapkit_size_t map_h_id_growsize(map_h_id * spm, mapkit_size_t used)
{
	return (mapkit_size_t) (4.0 * used / (3 * spm->minusedfactor + spm->maxfillfactor));
}

mapkit_size_t map_h_id_shrinksize(map_h_id * spm, mapkit_size_t used)
{
	return (mapkit_size_t) (4.0 * used / (spm->minusedfactor + 3 * spm->maxfillfactor));
}

mapkit_size_t map_h_id_meansize(map_h_id * spm, mapkit_size_t used)
{
	return (mapkit_size_t) (2.0 * used / (spm->minusedfactor + spm->maxfillfactor));
}

mapkit_error map_h_id_reallocate(map_h_id * spm, mapkit_size_t newsize)
{
	mapkit_size_t iindex;
	mapkit_size_t oldsize;
	map_h_id_storage *newcontents, *oldcontents;

	/*
	 * At least one free entry 
	 */
	if (newsize <= spm->used)
		newsize = spm->used + 1;
	newsize = mapkit_nextprime(newsize);
	if (newsize <= spm->used)
		MAPKIT_ERROR(MAPKIT_ETOOBIG);

	newcontents = (map_h_id_storage *) malloc(newsize * sizeof(*spm->contents));
	if (newcontents == NULL)
		MAPKIT_ERROR(MAPKIT_ENOMEM);

	/*
	 * Initialize all entries to "free" 
	 */
	for (iindex = 0; iindex < newsize; iindex++)
		newcontents[iindex].state = MAPKIT_FREESLOT;

	oldcontents = spm->contents;
	oldsize = spm->size;
	spm->contents = newcontents;
	spm->size = newsize;

#ifdef MAPKIT_DEBUG
	fprintf(stderr, "MAPKIT: reallocate %ld -> %ld\n", (long) oldsize, (long) newsize);
#endif

	spm->maxfill = (mapkit_size_t) (newsize * spm->maxfillfactor);
	/*
	 * At least one free entry 
	 */
	if (spm->maxfill >= newsize)
		spm->maxfill = newsize - 1;
	spm->minused = (mapkit_size_t) (newsize * spm->minusedfactor);
	spm->used = 0;

	if (oldcontents != NULL) {
		int used = 0;
		int key;
		double defaultvalue = spm->defaultvalue;
		int notalwaysdefault = !spm->alwaysdefault;

		/*
		 * Copy all entries from old to new 
		 */
		for (iindex = 0; iindex < oldsize; iindex++)
			if (oldcontents[iindex].state == MAPKIT_FULLSLOT) {
				mapkit_size_t ins_iindex;
				mapkit_hash_t hash;
				map_h_id_storage *contents;

				key = oldcontents[iindex].key;

				/*
				 * Fast path 
				 */
				ins_iindex = (hash = mapkit_hash((mapkit_hash_t) key)) % spm->size;
				contents = &(newcontents[ins_iindex]);

				if (contents->state != MAPKIT_FREESLOT) {
					ins_iindex = map_h_id_insertionindex(spm, key, hash);
					contents = &(newcontents[ins_iindex]);
				}
#ifdef MAPKIT_COLLISIONS
				else
					spm->insertionindexs++;
#endif
				if (notalwaysdefault || (!((oldcontents[iindex].value) == (defaultvalue)))) {
					contents->value = oldcontents[iindex].value;
					contents->state = MAPKIT_FULLSLOT;
					contents->key = key;
					used++;
				}
			}
		free(oldcontents);
		spm->used = used;
	}
	spm->fill = spm->used;

	return MAPKIT_OK;
}

double map_h_id_value_s(map_h_id * spm, int key, mapkit_hash_t hash)
{
	mapkit_size_t iindex;

	iindex = map_h_id_keyindex(spm, key, hash);

	if (iindex < 0) {
		if (spm->alwaysdefault)
			return spm->defaultvalue;
		else
			MAPKIT_FATAL_ERROR(MAPKIT_EKEYNOTFOUND);
	}

	return spm->contents[iindex].value;
}

mapkit_error map_h_id_get_s(map_h_id * spm, int key, double *value, mapkit_hash_t hash)
{
	mapkit_size_t iindex;

	iindex = map_h_id_keyindex(spm, key, hash);

	if (iindex < 0) {
		if (spm->alwaysdefault) {
			*value = spm->defaultvalue;
			return MAPKIT_OK;
		} else
			MAPKIT_ERROR(MAPKIT_EKEYNOTFOUND);
	}

	*value = spm->contents[iindex].value;
	return MAPKIT_OK;
}

mapkit_error map_h_id_set_s(map_h_id * spm, int key, double value, mapkit_hash_t hash)
{
	mapkit_size_t iindex;

	iindex = map_h_id_insertionindex(spm, key, hash);

	if (iindex < 0)
		/*
		 * FULLSLOT 
		 */
		spm->contents[-iindex - 1].value = value;
	else {
		map_h_id_storage *element = &(spm->contents[iindex]);
		int ffree = element->state == MAPKIT_FREESLOT;

		element->state = MAPKIT_FULLSLOT;
		element->value = value;
		element->key = key;
		spm->used++;

		if (ffree && ((++spm->fill) > spm->maxfill)) {
#ifdef MAPKIT_DEBUG
			fprintf(stderr, "MAPKIT: fill > maxfill\n");
#endif
			return map_h_id_reallocate(spm, map_h_id_growsize(spm, spm->used));
		}
	}
	return MAPKIT_OK;
}

double *map_h_id_insertptr_s(map_h_id * spm, int key, mapkit_hash_t hash)
{
	mapkit_size_t iindex;

	iindex = map_h_id_insertionindex(spm, key, hash);

	if (iindex < 0)
		return &(spm->contents[-iindex - 1].value);
	else {
		map_h_id_storage *element = &(spm->contents[iindex]);

		if (element->state == MAPKIT_FREESLOT) {
			/*
			 * FREESLOT 
			 */
			if (spm->fill >= spm->maxfill) {
				mapkit_error err;

#ifdef MAPKIT_DEBUG
				fprintf(stderr, "MAPKIT: fill => maxfill before insert\n");
#endif
				/*
				 * Must reallocate -before- inserting defaultvalue 
				 */
				err = map_h_id_reallocate(spm, map_h_id_growsize(spm, spm->used + 1));
				if (err) {
					MAPKIT_ERROR_NORET(err);
					return NULL;
				}

				iindex = map_h_id_insertionindex(spm, key, hash);
				/*
				 * FREESLOT 
				 */
				spm->contents[iindex].state = MAPKIT_FULLSLOT;
				spm->contents[iindex].key = key;
				spm->contents[iindex].value = spm->defaultvalue;
				spm->used++;
				spm->fill++;
				return &(spm->contents[iindex].value);
			} else
				spm->fill++;
		}

		element->state = MAPKIT_FULLSLOT;
		element->key = key;
		element->value = spm->defaultvalue;
		spm->used++;
		return &(element->value);
	}
}

double *map_h_id_ptr_s(map_h_id * spm, int key, mapkit_hash_t hash)
{
	mapkit_size_t iindex;

	iindex = map_h_id_keyindex(spm, key, hash);

	if ((iindex >= 0) && ((!spm->alwaysdefault) || (!((spm->contents[iindex].value) == (spm->defaultvalue)))))
		return &(spm->contents[iindex].value);
	else
		return NULL;
}

mapkit_error map_h_id_remove_s(map_h_id * spm, int key, mapkit_hash_t hash)
{
	mapkit_size_t iindex;

	iindex = map_h_id_keyindex(spm, key, hash);

	if (iindex < 0) {
		if (spm->alwaysdefault)
			return MAPKIT_OK;
		else
			return MAPKIT_EKEYNOTFOUND;
	}

	spm->contents[iindex].state = MAPKIT_DELETEDSLOT;
	spm->used--;
	if (spm->used < spm->minused) {
#ifdef MAPKIT_DEBUG
		fprintf(stderr, "MAPKIT: used < minused\n");
#endif
		return map_h_id_reallocate(spm, map_h_id_shrinksize(spm, spm->used));
	}

	return MAPKIT_OK;
}

mapkit_size_t map_h_id_keyindex(map_h_id * spm, int key, mapkit_hash_t hash)
{
	mapkit_size_t iindex, decrement;

	signed char state;

	iindex = hash % spm->size;
	decrement = (hash % (spm->size - 2));
	decrement += (decrement == 0);

#ifdef MAPKIT_COLLISIONS
	spm->keyindexs++;
#endif

	while ((state = spm->contents[iindex].state) != MAPKIT_FREESLOT
	       && (state == MAPKIT_DELETEDSLOT || (!((spm->contents[iindex].key) == (key))))) {
#ifdef MAPKIT_COLLISIONS
		spm->keyindex_collisions++;
#endif
		iindex -= decrement;
		if (iindex < 0)
			iindex += spm->size;
	}

	if (state == MAPKIT_FREESLOT)
		return MAPKIT_KEYNOTFOUND;
	return iindex;
}

mapkit_size_t map_h_id_insertionindex(map_h_id * spm, int key, mapkit_hash_t hash)
{
	mapkit_size_t iindex, decrement;
	signed char state;

#ifdef MAPKIT_COLLISIONS
	spm->insertionindexs++;
#endif

	iindex = hash % spm->size;

	/*
	 * Fast path (largely superfluous) 
	 */
	if ((state = spm->contents[iindex].state) == MAPKIT_FREESLOT)
		return iindex;
	if ((state == MAPKIT_FULLSLOT)
	    && ((spm->contents[iindex].key) == (key)))
		return -iindex - 1;

	decrement = (hash % (spm->size - 2));
	decrement += (decrement == 0);

	while ((state == MAPKIT_FULLSLOT)
	       && (!((spm->contents[iindex].key) == (key)))) {
#ifdef MAPKIT_COLLISIONS
		spm->insertionindex_collisions++;
#endif
		iindex -= decrement;
		if (iindex < 0)
			iindex += spm->size;
		state = spm->contents[iindex].state;
	}

	if (state == MAPKIT_DELETEDSLOT) {
		mapkit_size_t index2 = iindex;

		while ((state = spm->contents[iindex].state) != MAPKIT_FREESLOT && ((state == MAPKIT_DELETEDSLOT)
										    || (!((spm->contents[iindex].key) == (key))))) {
			iindex -= decrement;
			if (iindex < 0)
				iindex += spm->size;
			state = spm->contents[iindex].state;
		}
		if (state == MAPKIT_FREESLOT)
			return index2;
	}

	if (state == MAPKIT_FULLSLOT)
		return -iindex - 1;
	return iindex;
}

mapkit_size_t map_h_id_next(map_h_id * spm, mapkit_size_t iindex)
{
	mapkit_size_t size = spm->size;
	map_h_id_storage *pos_contents;
	double defaultvalue = spm->defaultvalue;
	int notalwaysdefault = !spm->alwaysdefault;

	pos_contents = &(spm->contents[++iindex]);

	for (; iindex < size; iindex++, pos_contents++)
		if ((pos_contents->state == MAPKIT_FULLSLOT)
		    && (notalwaysdefault || (!((pos_contents->value) == (defaultvalue)))))
			return iindex;

	return -1;
}

map_h_id_storage *map_h_id_nextptr(map_h_id * spm, map_h_id_storage * pos_contents)
{
	map_h_id_storage *end = &(spm->contents[spm->size]);
	double defaultvalue = spm->defaultvalue;
	int notalwaysdefault = !spm->alwaysdefault;

	if (pos_contents == NULL)
		pos_contents = spm->contents;
	else {
		pos_contents++;
		if (pos_contents <= spm->contents)
			return NULL;
	}

	for (; pos_contents < end; pos_contents++)
		if ((pos_contents->state == MAPKIT_FULLSLOT)
		    && (notalwaysdefault || (!((pos_contents->value) == (defaultvalue)))))
			return pos_contents;

	return NULL;
}

mapkit_error map_h_id_getall(map_h_id * spm, map_h_id_element ** array, mapkit_size_t * count)
{
	mapkit_size_t iindex;
	mapkit_size_t size = spm->size, vcount = 0;
	map_h_id_element *pos_array;
	map_h_id_storage *pos_contents;
	double defaultvalue = spm->defaultvalue;
	int notalwaysdefault = !spm->alwaysdefault;

	pos_array = *array = (map_h_id_element *) malloc(sizeof(**array) * spm->used);
	if (*array == NULL)
		MAPKIT_ERROR(MAPKIT_ENOMEM);

	pos_contents = spm->contents;

	for (iindex = 0; iindex < size; iindex++, pos_contents++)
		if ((pos_contents->state == MAPKIT_FULLSLOT)
		    && (notalwaysdefault || (!((pos_contents->value) == (defaultvalue))))) {
			pos_array->key = pos_contents->key;
			pos_array->value = pos_contents->value;
			pos_array++;
			vcount++;
		}
	*count = vcount;

	return MAPKIT_OK;
}

mapkit_error map_h_id_clean(map_h_id * spm)
{
	mapkit_size_t iindex, count = 0;
	mapkit_size_t size = spm->size;
	map_h_id_storage *pos_contents;
	double defaultvalue = spm->defaultvalue;

	pos_contents = spm->contents;
	for (iindex = 0; iindex < size; iindex++, pos_contents++)
		if ((pos_contents->state == MAPKIT_FULLSLOT)
		    && ((pos_contents->value) == (defaultvalue))) {
			pos_contents->state = MAPKIT_DELETEDSLOT;
			count++;
		}

	spm->used -= count;
	if (spm->used < spm->minused) {
#ifdef MAPKIT_DEBUG
		fprintf(stderr, "MAPKIT: used < minused\n");
#endif
		return map_h_id_reallocate(spm, map_h_id_meansize(spm, spm->used));
	}

	return MAPKIT_OK;
}

int map_h_id_compare(const void *e1, const void *e2)
{
	int key1 = ((const map_h_id_element *) e1)->key;
	int key2 = ((const map_h_id_element *) e2)->key;

	return ((key1) < (key2) ? -1 : ((key1) == (key2) ? 0 : 1));
}

mapkit_error map_h_id_getall_sorted(map_h_id * spm, map_h_id_element ** array, mapkit_size_t * count)
{
	mapkit_error err;

	err = map_h_id_getall(spm, array, count);
	if (err)
		MAPKIT_ERROR(err);

	qsort(*array, *count, sizeof(**array), map_h_id_compare);

	return MAPKIT_OK;
}

mapkit_error map_h_id_setall(map_h_id * spm, map_h_id_element * array, mapkit_size_t count)
{
	mapkit_size_t array_iindex;
	mapkit_error err;

	err = map_h_id_ensurecapacity(spm, spm->used + count);
	if (err)
		MAPKIT_ERROR(err);

	if (spm->alwaysdefault)
		/*
		 * Prevent shrinking 
		 */
		spm->minused = 0;

	for (array_iindex = 0; array_iindex < count; array_iindex++) {
		err = map_h_id_set(spm, array[array_iindex].key, array[array_iindex].value);
		if (err)
			MAPKIT_ERROR(err);
	}

	if (spm->alwaysdefault)
		map_h_id_adjustcapacity(spm);

	return MAPKIT_OK;
}

mapkit_error map_h_id_removeall(map_h_id * spm, int *array, mapkit_size_t count)
{
	mapkit_size_t array_iindex;
	mapkit_error err;

	/*
	 * Prevent shrinking 
	 */
	spm->minused = 0;

	for (array_iindex = 0; array_iindex < count; array_iindex++) {
		err = map_h_id_remove(spm, array[array_iindex]);
		if (err)
			MAPKIT_ERROR(err);
	}

	map_h_id_adjustcapacity(spm);

	return MAPKIT_OK;
}

void map_h_id_printstats(map_h_id * spm)
{
	fprintf(stderr, "MAPKIT: map_h_id statistics\n");
	fprintf(stderr, "MAPKIT: alwaysdefault = %d\n", spm->alwaysdefault);
	fprintf(stderr, "MAPKIT: minused = %ld, maxfill = %ld\n", (long) spm->minused, (long) spm->maxfill);
	fprintf(stderr, "MAPKIT: minusedfactor = %g, maxfillfactor = %g\n", spm->minusedfactor, spm->maxfillfactor);
#ifdef MAPKIT_COLLISIONS
	fprintf(stderr, "MAPKIT: insertionindexs = %lu, collisions = %lu\n", (unsigned long) spm->insertionindexs,
		(unsigned long) spm->insertionindex_collisions);
	fprintf(stderr, "MAPKIT: keyindexs = %lu, collisions = %lu\n", (unsigned long) spm->keyindexs, (unsigned long) spm->keyindex_collisions);
#endif
}

#endif							       /* MAPKIT_map_h_id */

#ifdef MAPKIT_map_h_ivp

/*
  Implementation for map_h_ivp (int -> void*)
  Default value : NULL
  Uses a state field.
*/

/* Static prototypes */

/* Return the iindex of key, or MAPKIT_KEYNOTFOUND */
static mapkit_size_t map_h_ivp_keyindex(map_h_ivp * spm, int key, mapkit_hash_t hash);

/* Return the iindex of key or -(insertion iindex)-1 if key not found */
static mapkit_size_t map_h_ivp_insertionindex(map_h_ivp * spm, int key, mapkit_hash_t hash);

/* Implementation */

mapkit_error map_h_ivp_init(map_h_ivp * spm)
{
	return map_h_ivp_init_hint(spm, MAPKIT_DEFAULT_EXPECTEDUSED);
}

mapkit_error map_h_ivp_init_hint(map_h_ivp * spm, mapkit_size_t used)
{
#ifdef MAPKIT_DEBUG
	fprintf(stderr, "MAPKIT: init\n");
#endif

	spm->size = 0;
	spm->fill = 0;
	spm->used = 0;
	spm->maxfillfactor = 0.5;
	spm->minusedfactor = 0.2;
	spm->contents = NULL;
#ifdef MAPKIT_COLLISIONS
	spm->insertionindexs = spm->insertionindex_collisions = 0;
	spm->keyindexs = spm->keyindex_collisions = 0;
#endif
	spm->defaultvalue = NULL;
	spm->alwaysdefault = 0;

	return map_h_ivp_reallocate(spm, map_h_ivp_meansize(spm, used));
}

mapkit_error map_h_ivp_ensurecapacity(map_h_ivp * spm, mapkit_size_t used)
{
	if (used > (spm->used + spm->maxfill - spm->fill)) {
#ifdef MAPKIT_DEBUG
		fprintf(stderr, "MAPKIT: need more capacity\n");
#endif
		return map_h_ivp_reallocate(spm, map_h_ivp_meansize(spm, used));
	} else
		return MAPKIT_OK;
}

mapkit_error map_h_ivp_adjustcapacity(map_h_ivp * spm)
{
	spm->minused = (mapkit_size_t) (spm->size * spm->minusedfactor);
	spm->maxfill = (mapkit_size_t) (spm->size * spm->maxfillfactor);

	if (spm->used < spm->minused) {
#ifdef MAPKIT_DEBUG
		fprintf(stderr, "MAPKIT: used < minused\n");
#endif
		return map_h_ivp_reallocate(spm, map_h_ivp_meansize(spm, spm->used));
	} else if (spm->fill > spm->maxfill) {
#ifdef MAPKIT_DEBUG
		fprintf(stderr, "MAPKIT: fill > maxfill\n");
#endif
		return map_h_ivp_reallocate(spm, map_h_ivp_meansize(spm, spm->used));
	} else
		return MAPKIT_OK;
}

void map_h_ivp_free(map_h_ivp * spm)
{
	free(spm->contents);
	spm->contents = NULL;
	spm->size = 0;
	spm->fill = 0;
	spm->used = 0;
	spm->maxfill = 0;
	spm->minused = 0;
#ifdef MAPKIT_COLLISIONS
	spm->insertionindexs = spm->insertionindex_collisions = 0;
	spm->keyindexs = spm->keyindex_collisions = 0;
#endif
}

mapkit_error map_h_ivp_copy(map_h_ivp * to, map_h_ivp * from)
{
	map_h_ivp_storage *contentscopy;
	size_t size = from->size * sizeof(*from->contents);

	contentscopy = (map_h_ivp_storage *) malloc(size);
	if (contentscopy == NULL)
		MAPKIT_ERROR(MAPKIT_ENOMEM);

	memcpy(to, from, sizeof(map_h_ivp));
	to->contents = contentscopy;
	memcpy(to->contents, from->contents, size);

	return MAPKIT_OK;
}

mapkit_size_t map_h_ivp_growsize(map_h_ivp * spm, mapkit_size_t used)
{
	return (mapkit_size_t) (4.0 * used / (3 * spm->minusedfactor + spm->maxfillfactor));
}

mapkit_size_t map_h_ivp_shrinksize(map_h_ivp * spm, mapkit_size_t used)
{
	return (mapkit_size_t) (4.0 * used / (spm->minusedfactor + 3 * spm->maxfillfactor));
}

mapkit_size_t map_h_ivp_meansize(map_h_ivp * spm, mapkit_size_t used)
{
	return (mapkit_size_t) (2.0 * used / (spm->minusedfactor + spm->maxfillfactor));
}

mapkit_error map_h_ivp_reallocate(map_h_ivp * spm, mapkit_size_t newsize)
{
	mapkit_size_t iindex;
	mapkit_size_t oldsize;
	map_h_ivp_storage *newcontents, *oldcontents;

	/*
	 * At least one free entry 
	 */
	if (newsize <= spm->used)
		newsize = spm->used + 1;
	newsize = mapkit_nextprime(newsize);
	if (newsize <= spm->used)
		MAPKIT_ERROR(MAPKIT_ETOOBIG);

	newcontents = (map_h_ivp_storage *) malloc(newsize * sizeof(*spm->contents));
	if (newcontents == NULL)
		MAPKIT_ERROR(MAPKIT_ENOMEM);

	/*
	 * Initialize all entries to "free" 
	 */
	for (iindex = 0; iindex < newsize; iindex++)
		newcontents[iindex].state = MAPKIT_FREESLOT;

	oldcontents = spm->contents;
	oldsize = spm->size;
	spm->contents = newcontents;
	spm->size = newsize;

#ifdef MAPKIT_DEBUG
	fprintf(stderr, "MAPKIT: reallocate %ld -> %ld\n", (long) oldsize, (long) newsize);
#endif

	spm->maxfill = (mapkit_size_t) (newsize * spm->maxfillfactor);
	/*
	 * At least one free entry 
	 */
	if (spm->maxfill >= newsize)
		spm->maxfill = newsize - 1;
	spm->minused = (mapkit_size_t) (newsize * spm->minusedfactor);
	spm->used = 0;

	if (oldcontents != NULL) {
		int used = 0;
		int key;
		void *defaultvalue = spm->defaultvalue;
		int notalwaysdefault = !spm->alwaysdefault;

		/*
		 * Copy all entries from old to new 
		 */
		for (iindex = 0; iindex < oldsize; iindex++)
			if (oldcontents[iindex].state == MAPKIT_FULLSLOT) {
				mapkit_size_t ins_iindex;
				mapkit_hash_t hash;
				map_h_ivp_storage *contents;

				key = oldcontents[iindex].key;

				/*
				 * Fast path 
				 */
				ins_iindex = (hash = mapkit_hash((mapkit_hash_t) key)) % spm->size;
				contents = &(newcontents[ins_iindex]);

				if (contents->state != MAPKIT_FREESLOT) {
					ins_iindex = map_h_ivp_insertionindex(spm, key, hash);
					contents = &(newcontents[ins_iindex]);
				}
#ifdef MAPKIT_COLLISIONS
				else
					spm->insertionindexs++;
#endif
				if (notalwaysdefault || (!((oldcontents[iindex].value) == (defaultvalue)))) {
					contents->value = oldcontents[iindex].value;
					contents->state = MAPKIT_FULLSLOT;
					contents->key = key;
					used++;
				}
			}
		free(oldcontents);
		spm->used = used;
	}
	spm->fill = spm->used;

	return MAPKIT_OK;
}

void *map_h_ivp_value_s(map_h_ivp * spm, int key, mapkit_hash_t hash)
{
	mapkit_size_t iindex;

	iindex = map_h_ivp_keyindex(spm, key, hash);

	if (iindex < 0) {
		if (spm->alwaysdefault)
			return spm->defaultvalue;
		else
			MAPKIT_FATAL_ERROR(MAPKIT_EKEYNOTFOUND);
	}

	return spm->contents[iindex].value;
}

mapkit_error map_h_ivp_get_s(map_h_ivp * spm, int key, void **value, mapkit_hash_t hash)
{
	mapkit_size_t iindex;

	iindex = map_h_ivp_keyindex(spm, key, hash);

	if (iindex < 0) {
		if (spm->alwaysdefault) {
			*value = spm->defaultvalue;
			return MAPKIT_OK;
		} else
			MAPKIT_ERROR(MAPKIT_EKEYNOTFOUND);
	}

	*value = spm->contents[iindex].value;
	return MAPKIT_OK;
}

mapkit_error map_h_ivp_set_s(map_h_ivp * spm, int key, void *value, mapkit_hash_t hash)
{
	mapkit_size_t iindex;

	iindex = map_h_ivp_insertionindex(spm, key, hash);

	if (iindex < 0)
		/*
		 * FULLSLOT 
		 */
		spm->contents[-iindex - 1].value = value;
	else {
		map_h_ivp_storage *element = &(spm->contents[iindex]);
		int ffree = element->state == MAPKIT_FREESLOT;

		element->state = MAPKIT_FULLSLOT;
		element->value = value;
		element->key = key;
		spm->used++;

		if (ffree && ((++spm->fill) > spm->maxfill)) {
#ifdef MAPKIT_DEBUG
			fprintf(stderr, "MAPKIT: fill > maxfill\n");
#endif
			return map_h_ivp_reallocate(spm, map_h_ivp_growsize(spm, spm->used));
		}
	}
	return MAPKIT_OK;
}

void **map_h_ivp_insertptr_s(map_h_ivp * spm, int key, mapkit_hash_t hash)
{
	mapkit_size_t iindex;

	iindex = map_h_ivp_insertionindex(spm, key, hash);

	if (iindex < 0)
		return &(spm->contents[-iindex - 1].value);
	else {
		map_h_ivp_storage *element = &(spm->contents[iindex]);

		if (element->state == MAPKIT_FREESLOT) {
			/*
			 * FREESLOT 
			 */
			if (spm->fill >= spm->maxfill) {
				mapkit_error err;

#ifdef MAPKIT_DEBUG
				fprintf(stderr, "MAPKIT: fill => maxfill before insert\n");
#endif
				/*
				 * Must reallocate -before- inserting defaultvalue 
				 */
				err = map_h_ivp_reallocate(spm, map_h_ivp_growsize(spm, spm->used + 1));
				if (err) {
					MAPKIT_ERROR_NORET(err);
					return NULL;
				}

				iindex = map_h_ivp_insertionindex(spm, key, hash);
				/*
				 * FREESLOT 
				 */
				spm->contents[iindex].state = MAPKIT_FULLSLOT;
				spm->contents[iindex].key = key;
				spm->contents[iindex].value = spm->defaultvalue;
				spm->used++;
				spm->fill++;
				return &(spm->contents[iindex].value);
			} else
				spm->fill++;
		}

		element->state = MAPKIT_FULLSLOT;
		element->key = key;
		element->value = spm->defaultvalue;
		spm->used++;
		return &(element->value);
	}
}

void **map_h_ivp_ptr_s(map_h_ivp * spm, int key, mapkit_hash_t hash)
{
	mapkit_size_t iindex;

	iindex = map_h_ivp_keyindex(spm, key, hash);

	if ((iindex >= 0) && ((!spm->alwaysdefault) || (!((spm->contents[iindex].value) == (spm->defaultvalue)))))
		return &(spm->contents[iindex].value);
	else
		return NULL;
}

mapkit_error map_h_ivp_remove_s(map_h_ivp * spm, int key, mapkit_hash_t hash)
{
	mapkit_size_t iindex;

	iindex = map_h_ivp_keyindex(spm, key, hash);

	if (iindex < 0) {
		if (spm->alwaysdefault)
			return MAPKIT_OK;
		else
			return MAPKIT_EKEYNOTFOUND;
	}

	spm->contents[iindex].state = MAPKIT_DELETEDSLOT;
	spm->used--;
	if (spm->used < spm->minused) {
#ifdef MAPKIT_DEBUG
		fprintf(stderr, "MAPKIT: used < minused\n");
#endif
		return map_h_ivp_reallocate(spm, map_h_ivp_shrinksize(spm, spm->used));
	}

	return MAPKIT_OK;
}

mapkit_size_t map_h_ivp_keyindex(map_h_ivp * spm, int key, mapkit_hash_t hash)
{
	mapkit_size_t iindex, decrement;

	signed char state;

	iindex = hash % spm->size;
	decrement = (hash % (spm->size - 2));
	decrement += (decrement == 0);

#ifdef MAPKIT_COLLISIONS
	spm->keyindexs++;
#endif

	while ((state = spm->contents[iindex].state) != MAPKIT_FREESLOT
	       && (state == MAPKIT_DELETEDSLOT || (!((spm->contents[iindex].key) == (key))))) {
#ifdef MAPKIT_COLLISIONS
		spm->keyindex_collisions++;
#endif
		iindex -= decrement;
		if (iindex < 0)
			iindex += spm->size;
	}

	if (state == MAPKIT_FREESLOT)
		return MAPKIT_KEYNOTFOUND;
	return iindex;
}

mapkit_size_t map_h_ivp_insertionindex(map_h_ivp * spm, int key, mapkit_hash_t hash)
{
	mapkit_size_t iindex, decrement;
	signed char state;

#ifdef MAPKIT_COLLISIONS
	spm->insertionindexs++;
#endif

	iindex = hash % spm->size;

	/*
	 * Fast path (largely superfluous) 
	 */
	if ((state = spm->contents[iindex].state) == MAPKIT_FREESLOT)
		return iindex;
	if ((state == MAPKIT_FULLSLOT)
	    && ((spm->contents[iindex].key) == (key)))
		return -iindex - 1;

	decrement = (hash % (spm->size - 2));
	decrement += (decrement == 0);

	while ((state == MAPKIT_FULLSLOT)
	       && (!((spm->contents[iindex].key) == (key)))) {
#ifdef MAPKIT_COLLISIONS
		spm->insertionindex_collisions++;
#endif
		iindex -= decrement;
		if (iindex < 0)
			iindex += spm->size;
		state = spm->contents[iindex].state;
	}

	if (state == MAPKIT_DELETEDSLOT) {
		mapkit_size_t index2 = iindex;

		while ((state = spm->contents[iindex].state) != MAPKIT_FREESLOT && ((state == MAPKIT_DELETEDSLOT)
										    || (!((spm->contents[iindex].key) == (key))))) {
			iindex -= decrement;
			if (iindex < 0)
				iindex += spm->size;
			state = spm->contents[iindex].state;
		}
		if (state == MAPKIT_FREESLOT)
			return index2;
	}

	if (state == MAPKIT_FULLSLOT)
		return -iindex - 1;
	return iindex;
}

mapkit_size_t map_h_ivp_next(map_h_ivp * spm, mapkit_size_t iindex)
{
	mapkit_size_t size = spm->size;
	map_h_ivp_storage *pos_contents;
	void *defaultvalue = spm->defaultvalue;
	int notalwaysdefault = !spm->alwaysdefault;

	pos_contents = &(spm->contents[++iindex]);

	for (; iindex < size; iindex++, pos_contents++)
		if ((pos_contents->state == MAPKIT_FULLSLOT)
		    && (notalwaysdefault || (!((pos_contents->value) == (defaultvalue)))))
			return iindex;

	return -1;
}

map_h_ivp_storage *map_h_ivp_nextptr(map_h_ivp * spm, map_h_ivp_storage * pos_contents)
{
	map_h_ivp_storage *end = &(spm->contents[spm->size]);
	void *defaultvalue = spm->defaultvalue;
	int notalwaysdefault = !spm->alwaysdefault;

	if (pos_contents == NULL)
		pos_contents = spm->contents;
	else {
		pos_contents++;
		if (pos_contents <= spm->contents)
			return NULL;
	}

	for (; pos_contents < end; pos_contents++)
		if ((pos_contents->state == MAPKIT_FULLSLOT)
		    && (notalwaysdefault || (!((pos_contents->value) == (defaultvalue)))))
			return pos_contents;

	return NULL;
}

mapkit_error map_h_ivp_getall(map_h_ivp * spm, map_h_ivp_element ** array, mapkit_size_t * count)
{
	mapkit_size_t iindex;
	mapkit_size_t size = spm->size, vcount = 0;
	map_h_ivp_element *pos_array;
	map_h_ivp_storage *pos_contents;
	void *defaultvalue = spm->defaultvalue;
	int notalwaysdefault = !spm->alwaysdefault;

	pos_array = *array = (map_h_ivp_element *) malloc(sizeof(**array) * spm->used);
	if (*array == NULL)
		MAPKIT_ERROR(MAPKIT_ENOMEM);

	pos_contents = spm->contents;

	for (iindex = 0; iindex < size; iindex++, pos_contents++)
		if ((pos_contents->state == MAPKIT_FULLSLOT)
		    && (notalwaysdefault || (!((pos_contents->value) == (defaultvalue))))) {
			pos_array->key = pos_contents->key;
			pos_array->value = pos_contents->value;
			pos_array++;
			vcount++;
		}
	*count = vcount;

	return MAPKIT_OK;
}

mapkit_error map_h_ivp_clean(map_h_ivp * spm)
{
	mapkit_size_t iindex, count = 0;
	mapkit_size_t size = spm->size;
	map_h_ivp_storage *pos_contents;
	void *defaultvalue = spm->defaultvalue;

	pos_contents = spm->contents;
	for (iindex = 0; iindex < size; iindex++, pos_contents++)
		if ((pos_contents->state == MAPKIT_FULLSLOT)
		    && ((pos_contents->value) == (defaultvalue))) {
			pos_contents->state = MAPKIT_DELETEDSLOT;
			count++;
		}

	spm->used -= count;
	if (spm->used < spm->minused) {
#ifdef MAPKIT_DEBUG
		fprintf(stderr, "MAPKIT: used < minused\n");
#endif
		return map_h_ivp_reallocate(spm, map_h_ivp_meansize(spm, spm->used));
	}

	return MAPKIT_OK;
}

int map_h_ivp_compare(const void *e1, const void *e2)
{
	int key1 = ((const map_h_ivp_element *) e1)->key;
	int key2 = ((const map_h_ivp_element *) e2)->key;

	return ((key1) < (key2) ? -1 : ((key1) == (key2) ? 0 : 1));
}

mapkit_error map_h_ivp_getall_sorted(map_h_ivp * spm, map_h_ivp_element ** array, mapkit_size_t * count)
{
	mapkit_error err;

	err = map_h_ivp_getall(spm, array, count);
	if (err)
		MAPKIT_ERROR(err);

	qsort(*array, *count, sizeof(**array), map_h_ivp_compare);

	return MAPKIT_OK;
}

mapkit_error map_h_ivp_setall(map_h_ivp * spm, map_h_ivp_element * array, mapkit_size_t count)
{
	mapkit_size_t array_iindex;
	mapkit_error err;

	err = map_h_ivp_ensurecapacity(spm, spm->used + count);
	if (err)
		MAPKIT_ERROR(err);

	if (spm->alwaysdefault)
		/*
		 * Prevent shrinking 
		 */
		spm->minused = 0;

	for (array_iindex = 0; array_iindex < count; array_iindex++) {
		err = map_h_ivp_set(spm, array[array_iindex].key, array[array_iindex].value);
		if (err)
			MAPKIT_ERROR(err);
	}

	if (spm->alwaysdefault)
		map_h_ivp_adjustcapacity(spm);

	return MAPKIT_OK;
}

mapkit_error map_h_ivp_removeall(map_h_ivp * spm, int *array, mapkit_size_t count)
{
	mapkit_size_t array_iindex;
	mapkit_error err;

	/*
	 * Prevent shrinking 
	 */
	spm->minused = 0;

	for (array_iindex = 0; array_iindex < count; array_iindex++) {
		err = map_h_ivp_remove(spm, array[array_iindex]);
		if (err)
			MAPKIT_ERROR(err);
	}

	map_h_ivp_adjustcapacity(spm);

	return MAPKIT_OK;
}

void map_h_ivp_printstats(map_h_ivp * spm)
{
	fprintf(stderr, "MAPKIT: map_h_ivp statistics\n");
	fprintf(stderr, "MAPKIT: alwaysdefault = %d\n", spm->alwaysdefault);
	fprintf(stderr, "MAPKIT: minused = %ld, maxfill = %ld\n", (long) spm->minused, (long) spm->maxfill);
	fprintf(stderr, "MAPKIT: minusedfactor = %g, maxfillfactor = %g\n", spm->minusedfactor, spm->maxfillfactor);
#ifdef MAPKIT_COLLISIONS
	fprintf(stderr, "MAPKIT: insertionindexs = %lu, collisions = %lu\n", (unsigned long) spm->insertionindexs,
		(unsigned long) spm->insertionindex_collisions);
	fprintf(stderr, "MAPKIT: keyindexs = %lu, collisions = %lu\n", (unsigned long) spm->keyindexs, (unsigned long) spm->keyindex_collisions);
#endif
}

#endif							       /* MAPKIT_map_h_ivp */

#ifdef MAPKIT_map_vpi

/*
  Implementation for map_vpi (void* -> int)
  Default value : 0
  Uses a state field.
*/

/* Static prototypes */

/* Return the iindex of key, or MAPKIT_KEYNOTFOUND */
static mapkit_size_t map_vpi_keyindex(map_vpi * spm, void *key);

/* Return the iindex of key or -(insertion iindex)-1 if key not found */
static mapkit_size_t map_vpi_insertionindex(map_vpi * spm, void *key);

/* Implementation */

mapkit_error map_vpi_init(map_vpi * spm)
{
	return map_vpi_init_hint(spm, MAPKIT_DEFAULT_EXPECTEDUSED);
}

mapkit_error map_vpi_init_hint(map_vpi * spm, mapkit_size_t used)
{
#ifdef MAPKIT_DEBUG
	fprintf(stderr, "MAPKIT: init\n");
#endif

	spm->size = 0;
	spm->fill = 0;
	spm->used = 0;
	spm->maxfillfactor = 0.5;
	spm->minusedfactor = 0.2;
	spm->contents = NULL;
#ifdef MAPKIT_COLLISIONS
	spm->insertionindexs = spm->insertionindex_collisions = 0;
	spm->keyindexs = spm->keyindex_collisions = 0;
#endif
	spm->defaultvalue = 0;
	spm->alwaysdefault = 0;

	return map_vpi_reallocate(spm, map_vpi_meansize(spm, used));
}

mapkit_error map_vpi_ensurecapacity(map_vpi * spm, mapkit_size_t used)
{
	if (used > (spm->used + spm->maxfill - spm->fill)) {
#ifdef MAPKIT_DEBUG
		fprintf(stderr, "MAPKIT: need more capacity\n");
#endif
		return map_vpi_reallocate(spm, map_vpi_meansize(spm, used));
	} else
		return MAPKIT_OK;
}

mapkit_error map_vpi_adjustcapacity(map_vpi * spm)
{
	spm->minused = (mapkit_size_t) (spm->size * spm->minusedfactor);
	spm->maxfill = (mapkit_size_t) (spm->size * spm->maxfillfactor);

	if (spm->used < spm->minused) {
#ifdef MAPKIT_DEBUG
		fprintf(stderr, "MAPKIT: used < minused\n");
#endif
		return map_vpi_reallocate(spm, map_vpi_meansize(spm, spm->used));
	} else if (spm->fill > spm->maxfill) {
#ifdef MAPKIT_DEBUG
		fprintf(stderr, "MAPKIT: fill > maxfill\n");
#endif
		return map_vpi_reallocate(spm, map_vpi_meansize(spm, spm->used));
	} else
		return MAPKIT_OK;
}

void map_vpi_free(map_vpi * spm)
{
	free(spm->contents);
	spm->contents = NULL;
	spm->size = 0;
	spm->fill = 0;
	spm->used = 0;
	spm->maxfill = 0;
	spm->minused = 0;
#ifdef MAPKIT_COLLISIONS
	spm->insertionindexs = spm->insertionindex_collisions = 0;
	spm->keyindexs = spm->keyindex_collisions = 0;
#endif
}

mapkit_error map_vpi_copy(map_vpi * to, map_vpi * from)
{
	map_vpi_storage *contentscopy;
	size_t size = from->size * sizeof(*from->contents);

	contentscopy = (map_vpi_storage *) malloc(size);
	if (contentscopy == NULL)
		MAPKIT_ERROR(MAPKIT_ENOMEM);

	memcpy(to, from, sizeof(map_vpi));
	to->contents = contentscopy;
	memcpy(to->contents, from->contents, size);

	return MAPKIT_OK;
}

mapkit_size_t map_vpi_growsize(map_vpi * spm, mapkit_size_t used)
{
	return (mapkit_size_t) (4.0 * used / (3 * spm->minusedfactor + spm->maxfillfactor));
}

mapkit_size_t map_vpi_shrinksize(map_vpi * spm, mapkit_size_t used)
{
	return (mapkit_size_t) (4.0 * used / (spm->minusedfactor + 3 * spm->maxfillfactor));
}

mapkit_size_t map_vpi_meansize(map_vpi * spm, mapkit_size_t used)
{
	return (mapkit_size_t) (2.0 * used / (spm->minusedfactor + spm->maxfillfactor));
}

mapkit_error map_vpi_reallocate(map_vpi * spm, mapkit_size_t newsize)
{
	mapkit_size_t iindex;
	mapkit_size_t oldsize;
	map_vpi_storage *newcontents, *oldcontents;

	/*
	 * At least one free entry 
	 */
	if (newsize <= spm->used)
		newsize = spm->used + 1;
	newsize = mapkit_nextprime(newsize);
	if (newsize <= spm->used)
		MAPKIT_ERROR(MAPKIT_ETOOBIG);

	newcontents = (map_vpi_storage *) malloc(newsize * sizeof(*spm->contents));
	if (newcontents == NULL)
		MAPKIT_ERROR(MAPKIT_ENOMEM);

	/*
	 * Initialize all entries to "free" 
	 */
	for (iindex = 0; iindex < newsize; iindex++)
		newcontents[iindex].state = MAPKIT_FREESLOT;

	oldcontents = spm->contents;
	oldsize = spm->size;
	spm->contents = newcontents;
	spm->size = newsize;

#ifdef MAPKIT_DEBUG
	fprintf(stderr, "MAPKIT: reallocate %ld -> %ld\n", (long) oldsize, (long) newsize);
#endif

	spm->maxfill = (mapkit_size_t) (newsize * spm->maxfillfactor);
	/*
	 * At least one free entry 
	 */
	if (spm->maxfill >= newsize)
		spm->maxfill = newsize - 1;
	spm->minused = (mapkit_size_t) (newsize * spm->minusedfactor);
	spm->used = 0;

	if (oldcontents != NULL) {
		int used = 0;
		void *key;
		int defaultvalue = spm->defaultvalue;
		int notalwaysdefault = !spm->alwaysdefault;

		/*
		 * Copy all entries from old to new 
		 */
		for (iindex = 0; iindex < oldsize; iindex++)
			if (oldcontents[iindex].state == MAPKIT_FULLSLOT) {
				mapkit_size_t ins_iindex;
				map_vpi_storage *contents;

				key = oldcontents[iindex].key;

				/*
				 * Fast path 
				 */
				ins_iindex = ((mapkit_hash_t) key) % spm->size;
				contents = &(newcontents[ins_iindex]);

				if (contents->state != MAPKIT_FREESLOT) {
					ins_iindex = map_vpi_insertionindex(spm, key);
					contents = &(newcontents[ins_iindex]);
				}
#ifdef MAPKIT_COLLISIONS
				else
					spm->insertionindexs++;
#endif
				if (notalwaysdefault || (!((oldcontents[iindex].value) == (defaultvalue)))) {
					contents->value = oldcontents[iindex].value;
					contents->state = MAPKIT_FULLSLOT;
					contents->key = key;
					used++;
				}
			}
		free(oldcontents);
		spm->used = used;
	}
	spm->fill = spm->used;

	return MAPKIT_OK;
}

int map_vpi_value_s(map_vpi * spm, void *key)
{
	mapkit_size_t iindex;

	iindex = map_vpi_keyindex(spm, key);

	if (iindex < 0) {
		if (spm->alwaysdefault)
			return spm->defaultvalue;
		else
			MAPKIT_FATAL_ERROR(MAPKIT_EKEYNOTFOUND);
	}

	return spm->contents[iindex].value;
}

mapkit_error map_vpi_get_s(map_vpi * spm, void *key, int *value)
{
	mapkit_size_t iindex;

	iindex = map_vpi_keyindex(spm, key);

	if (iindex < 0) {
		if (spm->alwaysdefault) {
			*value = spm->defaultvalue;
			return MAPKIT_OK;
		} else
			MAPKIT_ERROR(MAPKIT_EKEYNOTFOUND);
	}

	*value = spm->contents[iindex].value;
	return MAPKIT_OK;
}

mapkit_error map_vpi_set_s(map_vpi * spm, void *key, int value)
{
	mapkit_size_t iindex;

	iindex = map_vpi_insertionindex(spm, key);

	if (iindex < 0)
		/*
		 * FULLSLOT 
		 */
		spm->contents[-iindex - 1].value = value;
	else {
		map_vpi_storage *element = &(spm->contents[iindex]);
		int ffree = element->state == MAPKIT_FREESLOT;

		element->state = MAPKIT_FULLSLOT;
		element->value = value;
		element->key = key;
		spm->used++;

		if (ffree && ((++spm->fill) > spm->maxfill)) {
#ifdef MAPKIT_DEBUG
			fprintf(stderr, "MAPKIT: fill > maxfill\n");
#endif
			return map_vpi_reallocate(spm, map_vpi_growsize(spm, spm->used));
		}
	}
	return MAPKIT_OK;
}

int *map_vpi_insertptr_s(map_vpi * spm, void *key)
{
	mapkit_size_t iindex;

	iindex = map_vpi_insertionindex(spm, key);

	if (iindex < 0)
		return &(spm->contents[-iindex - 1].value);
	else {
		map_vpi_storage *element = &(spm->contents[iindex]);

		if (element->state == MAPKIT_FREESLOT) {
			/*
			 * FREESLOT 
			 */
			if (spm->fill >= spm->maxfill) {
				mapkit_error err;

#ifdef MAPKIT_DEBUG
				fprintf(stderr, "MAPKIT: fill => maxfill before insert\n");
#endif
				/*
				 * Must reallocate -before- inserting defaultvalue 
				 */
				err = map_vpi_reallocate(spm, map_vpi_growsize(spm, spm->used + 1));
				if (err) {
					MAPKIT_ERROR_NORET(err);
					return NULL;
				}

				iindex = map_vpi_insertionindex(spm, key);
				/*
				 * FREESLOT 
				 */
				spm->contents[iindex].state = MAPKIT_FULLSLOT;
				spm->contents[iindex].key = key;
				spm->contents[iindex].value = spm->defaultvalue;
				spm->used++;
				spm->fill++;
				return &(spm->contents[iindex].value);
			} else
				spm->fill++;
		}

		element->state = MAPKIT_FULLSLOT;
		element->key = key;
		element->value = spm->defaultvalue;
		spm->used++;
		return &(element->value);
	}
}

int *map_vpi_ptr_s(map_vpi * spm, void *key)
{
	mapkit_size_t iindex;

	iindex = map_vpi_keyindex(spm, key);

	if ((iindex >= 0) && ((!spm->alwaysdefault) || (!((spm->contents[iindex].value) == (spm->defaultvalue)))))
		return &(spm->contents[iindex].value);
	else
		return NULL;
}

mapkit_error map_vpi_remove_s(map_vpi * spm, void *key)
{
	mapkit_size_t iindex;

	iindex = map_vpi_keyindex(spm, key);

	if (iindex < 0) {
		if (spm->alwaysdefault)
			return MAPKIT_OK;
		else
			return MAPKIT_EKEYNOTFOUND;
	}

	spm->contents[iindex].state = MAPKIT_DELETEDSLOT;
	spm->used--;
	if (spm->used < spm->minused) {
#ifdef MAPKIT_DEBUG
		fprintf(stderr, "MAPKIT: used < minused\n");
#endif
		return map_vpi_reallocate(spm, map_vpi_shrinksize(spm, spm->used));
	}

	return MAPKIT_OK;
}

mapkit_size_t map_vpi_keyindex(map_vpi * spm, void *key)
{
	mapkit_size_t iindex, decrement;

	signed char state;

	iindex = ((mapkit_hash_t) key) % spm->size;
	decrement = (((mapkit_hash_t) key) % (spm->size - 2));
	decrement += (decrement == 0);

#ifdef MAPKIT_COLLISIONS
	spm->keyindexs++;
#endif

	while ((state = spm->contents[iindex].state) != MAPKIT_FREESLOT
	       && (state == MAPKIT_DELETEDSLOT || (!((spm->contents[iindex].key) == (key))))) {
#ifdef MAPKIT_COLLISIONS
		spm->keyindex_collisions++;
#endif
		iindex -= decrement;
		if (iindex < 0)
			iindex += spm->size;
	}

	if (state == MAPKIT_FREESLOT)
		return MAPKIT_KEYNOTFOUND;
	return iindex;
}

mapkit_size_t map_vpi_insertionindex(map_vpi * spm, void *key)
{
	mapkit_size_t iindex, decrement;
	signed char state;

#ifdef MAPKIT_COLLISIONS
	spm->insertionindexs++;
#endif

	iindex = ((mapkit_hash_t) key) % spm->size;

	/*
	 * Fast path (largely superfluous) 
	 */
	if ((state = spm->contents[iindex].state) == MAPKIT_FREESLOT)
		return iindex;
	if ((state == MAPKIT_FULLSLOT)
	    && ((spm->contents[iindex].key) == (key)))
		return -iindex - 1;

	decrement = (((mapkit_hash_t) key) % (spm->size - 2));
	decrement += (decrement == 0);

	while ((state == MAPKIT_FULLSLOT)
	       && (!((spm->contents[iindex].key) == (key)))) {
#ifdef MAPKIT_COLLISIONS
		spm->insertionindex_collisions++;
#endif
		iindex -= decrement;
		if (iindex < 0)
			iindex += spm->size;
		state = spm->contents[iindex].state;
	}

	if (state == MAPKIT_DELETEDSLOT) {
		mapkit_size_t index2 = iindex;

		while ((state = spm->contents[iindex].state) != MAPKIT_FREESLOT && ((state == MAPKIT_DELETEDSLOT)
										    || (!((spm->contents[iindex].key) == (key))))) {
			iindex -= decrement;
			if (iindex < 0)
				iindex += spm->size;
			state = spm->contents[iindex].state;
		}
		if (state == MAPKIT_FREESLOT)
			return index2;
	}

	if (state == MAPKIT_FULLSLOT)
		return -iindex - 1;
	return iindex;
}

mapkit_size_t map_vpi_next(map_vpi * spm, mapkit_size_t iindex)
{
	mapkit_size_t size = spm->size;
	map_vpi_storage *pos_contents;
	int defaultvalue = spm->defaultvalue;
	int notalwaysdefault = !spm->alwaysdefault;

	pos_contents = &(spm->contents[++iindex]);

	for (; iindex < size; iindex++, pos_contents++)
		if ((pos_contents->state == MAPKIT_FULLSLOT)
		    && (notalwaysdefault || (!((pos_contents->value) == (defaultvalue)))))
			return iindex;

	return -1;
}

map_vpi_storage *map_vpi_nextptr(map_vpi * spm, map_vpi_storage * pos_contents)
{
	map_vpi_storage *end = &(spm->contents[spm->size]);
	int defaultvalue = spm->defaultvalue;
	int notalwaysdefault = !spm->alwaysdefault;

	if (pos_contents == NULL)
		pos_contents = spm->contents;
	else {
		pos_contents++;
		if (pos_contents <= spm->contents)
			return NULL;
	}

	for (; pos_contents < end; pos_contents++)
		if ((pos_contents->state == MAPKIT_FULLSLOT)
		    && (notalwaysdefault || (!((pos_contents->value) == (defaultvalue)))))
			return pos_contents;

	return NULL;
}

mapkit_error map_vpi_getall(map_vpi * spm, map_vpi_element ** array, mapkit_size_t * count)
{
	mapkit_size_t iindex;
	mapkit_size_t size = spm->size, vcount = 0;
	map_vpi_element *pos_array;
	map_vpi_storage *pos_contents;
	int defaultvalue = spm->defaultvalue;
	int notalwaysdefault = !spm->alwaysdefault;

	pos_array = *array = (map_vpi_element *) malloc(sizeof(**array) * spm->used);
	if (*array == NULL)
		MAPKIT_ERROR(MAPKIT_ENOMEM);

	pos_contents = spm->contents;

	for (iindex = 0; iindex < size; iindex++, pos_contents++)
		if ((pos_contents->state == MAPKIT_FULLSLOT)
		    && (notalwaysdefault || (!((pos_contents->value) == (defaultvalue))))) {
			pos_array->key = pos_contents->key;
			pos_array->value = pos_contents->value;
			pos_array++;
			vcount++;
		}
	*count = vcount;

	return MAPKIT_OK;
}

mapkit_error map_vpi_clean(map_vpi * spm)
{
	mapkit_size_t iindex, count = 0;
	mapkit_size_t size = spm->size;
	map_vpi_storage *pos_contents;
	int defaultvalue = spm->defaultvalue;

	pos_contents = spm->contents;
	for (iindex = 0; iindex < size; iindex++, pos_contents++)
		if ((pos_contents->state == MAPKIT_FULLSLOT)
		    && ((pos_contents->value) == (defaultvalue))) {
			pos_contents->state = MAPKIT_DELETEDSLOT;
			count++;
		}

	spm->used -= count;
	if (spm->used < spm->minused) {
#ifdef MAPKIT_DEBUG
		fprintf(stderr, "MAPKIT: used < minused\n");
#endif
		return map_vpi_reallocate(spm, map_vpi_meansize(spm, spm->used));
	}

	return MAPKIT_OK;
}

int map_vpi_compare(const void *e1, const void *e2)
{
	void *key1 = ((const map_vpi_element *) e1)->key;
	void *key2 = ((const map_vpi_element *) e2)->key;

	return ((key1) < (key2) ? -1 : ((key1) == (key2) ? 0 : 1));
}

mapkit_error map_vpi_getall_sorted(map_vpi * spm, map_vpi_element ** array, mapkit_size_t * count)
{
	mapkit_error err;

	err = map_vpi_getall(spm, array, count);
	if (err)
		MAPKIT_ERROR(err);

	qsort(*array, *count, sizeof(**array), map_vpi_compare);

	return MAPKIT_OK;
}

mapkit_error map_vpi_setall(map_vpi * spm, map_vpi_element * array, mapkit_size_t count)
{
	mapkit_size_t array_iindex;
	mapkit_error err;

	err = map_vpi_ensurecapacity(spm, spm->used + count);
	if (err)
		MAPKIT_ERROR(err);

	if (spm->alwaysdefault)
		/*
		 * Prevent shrinking 
		 */
		spm->minused = 0;

	for (array_iindex = 0; array_iindex < count; array_iindex++) {
		err = map_vpi_set(spm, array[array_iindex].key, array[array_iindex].value);
		if (err)
			MAPKIT_ERROR(err);
	}

	if (spm->alwaysdefault)
		map_vpi_adjustcapacity(spm);

	return MAPKIT_OK;
}

mapkit_error map_vpi_removeall(map_vpi * spm, void **array, mapkit_size_t count)
{
	mapkit_size_t array_iindex;
	mapkit_error err;

	/*
	 * Prevent shrinking 
	 */
	spm->minused = 0;

	for (array_iindex = 0; array_iindex < count; array_iindex++) {
		err = map_vpi_remove(spm, array[array_iindex]);
		if (err)
			MAPKIT_ERROR(err);
	}

	map_vpi_adjustcapacity(spm);

	return MAPKIT_OK;
}

void map_vpi_printstats(map_vpi * spm)
{
	fprintf(stderr, "MAPKIT: map_vpi statistics\n");
	fprintf(stderr, "MAPKIT: alwaysdefault = %d\n", spm->alwaysdefault);
	fprintf(stderr, "MAPKIT: minused = %ld, maxfill = %ld\n", (long) spm->minused, (long) spm->maxfill);
	fprintf(stderr, "MAPKIT: minusedfactor = %g, maxfillfactor = %g\n", spm->minusedfactor, spm->maxfillfactor);
#ifdef MAPKIT_COLLISIONS
	fprintf(stderr, "MAPKIT: insertionindexs = %lu, collisions = %lu\n", (unsigned long) spm->insertionindexs,
		(unsigned long) spm->insertionindex_collisions);
	fprintf(stderr, "MAPKIT: keyindexs = %lu, collisions = %lu\n", (unsigned long) spm->keyindexs, (unsigned long) spm->keyindex_collisions);
#endif
}

#endif							       /* MAPKIT_map_vpi */

#ifdef MAPKIT_map_vpd

/*
  Implementation for map_vpd (void* -> double)
  Default value : 0.0
  Uses a state field.
*/

/* Static prototypes */

/* Return the iindex of key, or MAPKIT_KEYNOTFOUND */
static mapkit_size_t map_vpd_keyindex(map_vpd * spm, void *key);

/* Return the iindex of key or -(insertion iindex)-1 if key not found */
static mapkit_size_t map_vpd_insertionindex(map_vpd * spm, void *key);

/* Implementation */

mapkit_error map_vpd_init(map_vpd * spm)
{
	return map_vpd_init_hint(spm, MAPKIT_DEFAULT_EXPECTEDUSED);
}

mapkit_error map_vpd_init_hint(map_vpd * spm, mapkit_size_t used)
{
#ifdef MAPKIT_DEBUG
	fprintf(stderr, "MAPKIT: init\n");
#endif

	spm->size = 0;
	spm->fill = 0;
	spm->used = 0;
	spm->maxfillfactor = 0.5;
	spm->minusedfactor = 0.2;
	spm->contents = NULL;
#ifdef MAPKIT_COLLISIONS
	spm->insertionindexs = spm->insertionindex_collisions = 0;
	spm->keyindexs = spm->keyindex_collisions = 0;
#endif
	spm->defaultvalue = 0.0;
	spm->alwaysdefault = 0;

	return map_vpd_reallocate(spm, map_vpd_meansize(spm, used));
}

mapkit_error map_vpd_ensurecapacity(map_vpd * spm, mapkit_size_t used)
{
	if (used > (spm->used + spm->maxfill - spm->fill)) {
#ifdef MAPKIT_DEBUG
		fprintf(stderr, "MAPKIT: need more capacity\n");
#endif
		return map_vpd_reallocate(spm, map_vpd_meansize(spm, used));
	} else
		return MAPKIT_OK;
}

mapkit_error map_vpd_adjustcapacity(map_vpd * spm)
{
	spm->minused = (mapkit_size_t) (spm->size * spm->minusedfactor);
	spm->maxfill = (mapkit_size_t) (spm->size * spm->maxfillfactor);

	if (spm->used < spm->minused) {
#ifdef MAPKIT_DEBUG
		fprintf(stderr, "MAPKIT: used < minused\n");
#endif
		return map_vpd_reallocate(spm, map_vpd_meansize(spm, spm->used));
	} else if (spm->fill > spm->maxfill) {
#ifdef MAPKIT_DEBUG
		fprintf(stderr, "MAPKIT: fill > maxfill\n");
#endif
		return map_vpd_reallocate(spm, map_vpd_meansize(spm, spm->used));
	} else
		return MAPKIT_OK;
}

void map_vpd_free(map_vpd * spm)
{
	free(spm->contents);
	spm->contents = NULL;
	spm->size = 0;
	spm->fill = 0;
	spm->used = 0;
	spm->maxfill = 0;
	spm->minused = 0;
#ifdef MAPKIT_COLLISIONS
	spm->insertionindexs = spm->insertionindex_collisions = 0;
	spm->keyindexs = spm->keyindex_collisions = 0;
#endif
}

mapkit_error map_vpd_copy(map_vpd * to, map_vpd * from)
{
	map_vpd_storage *contentscopy;
	size_t size = from->size * sizeof(*from->contents);

	contentscopy = (map_vpd_storage *) malloc(size);
	if (contentscopy == NULL)
		MAPKIT_ERROR(MAPKIT_ENOMEM);

	memcpy(to, from, sizeof(map_vpd));
	to->contents = contentscopy;
	memcpy(to->contents, from->contents, size);

	return MAPKIT_OK;
}

mapkit_size_t map_vpd_growsize(map_vpd * spm, mapkit_size_t used)
{
	return (mapkit_size_t) (4.0 * used / (3 * spm->minusedfactor + spm->maxfillfactor));
}

mapkit_size_t map_vpd_shrinksize(map_vpd * spm, mapkit_size_t used)
{
	return (mapkit_size_t) (4.0 * used / (spm->minusedfactor + 3 * spm->maxfillfactor));
}

mapkit_size_t map_vpd_meansize(map_vpd * spm, mapkit_size_t used)
{
	return (mapkit_size_t) (2.0 * used / (spm->minusedfactor + spm->maxfillfactor));
}

mapkit_error map_vpd_reallocate(map_vpd * spm, mapkit_size_t newsize)
{
	mapkit_size_t iindex;
	mapkit_size_t oldsize;
	map_vpd_storage *newcontents, *oldcontents;

	/*
	 * At least one free entry 
	 */
	if (newsize <= spm->used)
		newsize = spm->used + 1;
	newsize = mapkit_nextprime(newsize);
	if (newsize <= spm->used)
		MAPKIT_ERROR(MAPKIT_ETOOBIG);

	newcontents = (map_vpd_storage *) malloc(newsize * sizeof(*spm->contents));
	if (newcontents == NULL)
		MAPKIT_ERROR(MAPKIT_ENOMEM);

	/*
	 * Initialize all entries to "free" 
	 */
	for (iindex = 0; iindex < newsize; iindex++)
		newcontents[iindex].state = MAPKIT_FREESLOT;

	oldcontents = spm->contents;
	oldsize = spm->size;
	spm->contents = newcontents;
	spm->size = newsize;

#ifdef MAPKIT_DEBUG
	fprintf(stderr, "MAPKIT: reallocate %ld -> %ld\n", (long) oldsize, (long) newsize);
#endif

	spm->maxfill = (mapkit_size_t) (newsize * spm->maxfillfactor);
	/*
	 * At least one free entry 
	 */
	if (spm->maxfill >= newsize)
		spm->maxfill = newsize - 1;
	spm->minused = (mapkit_size_t) (newsize * spm->minusedfactor);
	spm->used = 0;

	if (oldcontents != NULL) {
		int used = 0;
		void *key;
		double defaultvalue = spm->defaultvalue;
		int notalwaysdefault = !spm->alwaysdefault;

		/*
		 * Copy all entries from old to new 
		 */
		for (iindex = 0; iindex < oldsize; iindex++)
			if (oldcontents[iindex].state == MAPKIT_FULLSLOT) {
				mapkit_size_t ins_iindex;
				map_vpd_storage *contents;

				key = oldcontents[iindex].key;

				/*
				 * Fast path 
				 */
				ins_iindex = ((mapkit_hash_t) key) % spm->size;
				contents = &(newcontents[ins_iindex]);

				if (contents->state != MAPKIT_FREESLOT) {
					ins_iindex = map_vpd_insertionindex(spm, key);
					contents = &(newcontents[ins_iindex]);
				}
#ifdef MAPKIT_COLLISIONS
				else
					spm->insertionindexs++;
#endif
				if (notalwaysdefault || (!((oldcontents[iindex].value) == (defaultvalue)))) {
					contents->value = oldcontents[iindex].value;
					contents->state = MAPKIT_FULLSLOT;
					contents->key = key;
					used++;
				}
			}
		free(oldcontents);
		spm->used = used;
	}
	spm->fill = spm->used;

	return MAPKIT_OK;
}

double map_vpd_value_s(map_vpd * spm, void *key)
{
	mapkit_size_t iindex;

	iindex = map_vpd_keyindex(spm, key);

	if (iindex < 0) {
		if (spm->alwaysdefault)
			return spm->defaultvalue;
		else
			MAPKIT_FATAL_ERROR(MAPKIT_EKEYNOTFOUND);
	}

	return spm->contents[iindex].value;
}

mapkit_error map_vpd_get_s(map_vpd * spm, void *key, double *value)
{
	mapkit_size_t iindex;

	iindex = map_vpd_keyindex(spm, key);

	if (iindex < 0) {
		if (spm->alwaysdefault) {
			*value = spm->defaultvalue;
			return MAPKIT_OK;
		} else
			MAPKIT_ERROR(MAPKIT_EKEYNOTFOUND);
	}

	*value = spm->contents[iindex].value;
	return MAPKIT_OK;
}

mapkit_error map_vpd_set_s(map_vpd * spm, void *key, double value)
{
	mapkit_size_t iindex;

	iindex = map_vpd_insertionindex(spm, key);

	if (iindex < 0)
		/*
		 * FULLSLOT 
		 */
		spm->contents[-iindex - 1].value = value;
	else {
		map_vpd_storage *element = &(spm->contents[iindex]);
		int ffree = element->state == MAPKIT_FREESLOT;

		element->state = MAPKIT_FULLSLOT;
		element->value = value;
		element->key = key;
		spm->used++;

		if (ffree && ((++spm->fill) > spm->maxfill)) {
#ifdef MAPKIT_DEBUG
			fprintf(stderr, "MAPKIT: fill > maxfill\n");
#endif
			return map_vpd_reallocate(spm, map_vpd_growsize(spm, spm->used));
		}
	}
	return MAPKIT_OK;
}

double *map_vpd_insertptr_s(map_vpd * spm, void *key)
{
	mapkit_size_t iindex;

	iindex = map_vpd_insertionindex(spm, key);

	if (iindex < 0)
		return &(spm->contents[-iindex - 1].value);
	else {
		map_vpd_storage *element = &(spm->contents[iindex]);

		if (element->state == MAPKIT_FREESLOT) {
			/*
			 * FREESLOT 
			 */
			if (spm->fill >= spm->maxfill) {
				mapkit_error err;

#ifdef MAPKIT_DEBUG
				fprintf(stderr, "MAPKIT: fill => maxfill before insert\n");
#endif
				/*
				 * Must reallocate -before- inserting defaultvalue 
				 */
				err = map_vpd_reallocate(spm, map_vpd_growsize(spm, spm->used + 1));
				if (err) {
					MAPKIT_ERROR_NORET(err);
					return NULL;
				}

				iindex = map_vpd_insertionindex(spm, key);
				/*
				 * FREESLOT 
				 */
				spm->contents[iindex].state = MAPKIT_FULLSLOT;
				spm->contents[iindex].key = key;
				spm->contents[iindex].value = spm->defaultvalue;
				spm->used++;
				spm->fill++;
				return &(spm->contents[iindex].value);
			} else
				spm->fill++;
		}

		element->state = MAPKIT_FULLSLOT;
		element->key = key;
		element->value = spm->defaultvalue;
		spm->used++;
		return &(element->value);
	}
}

double *map_vpd_ptr_s(map_vpd * spm, void *key)
{
	mapkit_size_t iindex;

	iindex = map_vpd_keyindex(spm, key);

	if ((iindex >= 0) && ((!spm->alwaysdefault) || (!((spm->contents[iindex].value) == (spm->defaultvalue)))))
		return &(spm->contents[iindex].value);
	else
		return NULL;
}

mapkit_error map_vpd_remove_s(map_vpd * spm, void *key)
{
	mapkit_size_t iindex;

	iindex = map_vpd_keyindex(spm, key);

	if (iindex < 0) {
		if (spm->alwaysdefault)
			return MAPKIT_OK;
		else
			return MAPKIT_EKEYNOTFOUND;
	}

	spm->contents[iindex].state = MAPKIT_DELETEDSLOT;
	spm->used--;
	if (spm->used < spm->minused) {
#ifdef MAPKIT_DEBUG
		fprintf(stderr, "MAPKIT: used < minused\n");
#endif
		return map_vpd_reallocate(spm, map_vpd_shrinksize(spm, spm->used));
	}

	return MAPKIT_OK;
}

mapkit_size_t map_vpd_keyindex(map_vpd * spm, void *key)
{
	mapkit_size_t iindex, decrement;

	signed char state;

	iindex = ((mapkit_hash_t) key) % spm->size;
	decrement = (((mapkit_hash_t) key) % (spm->size - 2));
	decrement += (decrement == 0);

#ifdef MAPKIT_COLLISIONS
	spm->keyindexs++;
#endif

	while ((state = spm->contents[iindex].state) != MAPKIT_FREESLOT
	       && (state == MAPKIT_DELETEDSLOT || (!((spm->contents[iindex].key) == (key))))) {
#ifdef MAPKIT_COLLISIONS
		spm->keyindex_collisions++;
#endif
		iindex -= decrement;
		if (iindex < 0)
			iindex += spm->size;
	}

	if (state == MAPKIT_FREESLOT)
		return MAPKIT_KEYNOTFOUND;
	return iindex;
}

mapkit_size_t map_vpd_insertionindex(map_vpd * spm, void *key)
{
	mapkit_size_t iindex, decrement;
	signed char state;

#ifdef MAPKIT_COLLISIONS
	spm->insertionindexs++;
#endif

	iindex = ((mapkit_hash_t) key) % spm->size;

	/*
	 * Fast path (largely superfluous) 
	 */
	if ((state = spm->contents[iindex].state) == MAPKIT_FREESLOT)
		return iindex;
	if ((state == MAPKIT_FULLSLOT)
	    && ((spm->contents[iindex].key) == (key)))
		return -iindex - 1;

	decrement = (((mapkit_hash_t) key) % (spm->size - 2));
	decrement += (decrement == 0);

	while ((state == MAPKIT_FULLSLOT)
	       && (!((spm->contents[iindex].key) == (key)))) {
#ifdef MAPKIT_COLLISIONS
		spm->insertionindex_collisions++;
#endif
		iindex -= decrement;
		if (iindex < 0)
			iindex += spm->size;
		state = spm->contents[iindex].state;
	}

	if (state == MAPKIT_DELETEDSLOT) {
		mapkit_size_t index2 = iindex;

		while ((state = spm->contents[iindex].state) != MAPKIT_FREESLOT && ((state == MAPKIT_DELETEDSLOT)
										    || (!((spm->contents[iindex].key) == (key))))) {
			iindex -= decrement;
			if (iindex < 0)
				iindex += spm->size;
			state = spm->contents[iindex].state;
		}
		if (state == MAPKIT_FREESLOT)
			return index2;
	}

	if (state == MAPKIT_FULLSLOT)
		return -iindex - 1;
	return iindex;
}

mapkit_size_t map_vpd_next(map_vpd * spm, mapkit_size_t iindex)
{
	mapkit_size_t size = spm->size;
	map_vpd_storage *pos_contents;
	double defaultvalue = spm->defaultvalue;
	int notalwaysdefault = !spm->alwaysdefault;

	pos_contents = &(spm->contents[++iindex]);

	for (; iindex < size; iindex++, pos_contents++)
		if ((pos_contents->state == MAPKIT_FULLSLOT)
		    && (notalwaysdefault || (!((pos_contents->value) == (defaultvalue)))))
			return iindex;

	return -1;
}

map_vpd_storage *map_vpd_nextptr(map_vpd * spm, map_vpd_storage * pos_contents)
{
	map_vpd_storage *end = &(spm->contents[spm->size]);
	double defaultvalue = spm->defaultvalue;
	int notalwaysdefault = !spm->alwaysdefault;

	if (pos_contents == NULL)
		pos_contents = spm->contents;
	else {
		pos_contents++;
		if (pos_contents <= spm->contents)
			return NULL;
	}

	for (; pos_contents < end; pos_contents++)
		if ((pos_contents->state == MAPKIT_FULLSLOT)
		    && (notalwaysdefault || (!((pos_contents->value) == (defaultvalue)))))
			return pos_contents;

	return NULL;
}

mapkit_error map_vpd_getall(map_vpd * spm, map_vpd_element ** array, mapkit_size_t * count)
{
	mapkit_size_t iindex;
	mapkit_size_t size = spm->size, vcount = 0;
	map_vpd_element *pos_array;
	map_vpd_storage *pos_contents;
	double defaultvalue = spm->defaultvalue;
	int notalwaysdefault = !spm->alwaysdefault;

	pos_array = *array = (map_vpd_element *) malloc(sizeof(**array) * spm->used);
	if (*array == NULL)
		MAPKIT_ERROR(MAPKIT_ENOMEM);

	pos_contents = spm->contents;

	for (iindex = 0; iindex < size; iindex++, pos_contents++)
		if ((pos_contents->state == MAPKIT_FULLSLOT)
		    && (notalwaysdefault || (!((pos_contents->value) == (defaultvalue))))) {
			pos_array->key = pos_contents->key;
			pos_array->value = pos_contents->value;
			pos_array++;
			vcount++;
		}
	*count = vcount;

	return MAPKIT_OK;
}

mapkit_error map_vpd_clean(map_vpd * spm)
{
	mapkit_size_t iindex, count = 0;
	mapkit_size_t size = spm->size;
	map_vpd_storage *pos_contents;
	double defaultvalue = spm->defaultvalue;

	pos_contents = spm->contents;
	for (iindex = 0; iindex < size; iindex++, pos_contents++)
		if ((pos_contents->state == MAPKIT_FULLSLOT)
		    && ((pos_contents->value) == (defaultvalue))) {
			pos_contents->state = MAPKIT_DELETEDSLOT;
			count++;
		}

	spm->used -= count;
	if (spm->used < spm->minused) {
#ifdef MAPKIT_DEBUG
		fprintf(stderr, "MAPKIT: used < minused\n");
#endif
		return map_vpd_reallocate(spm, map_vpd_meansize(spm, spm->used));
	}

	return MAPKIT_OK;
}

int map_vpd_compare(const void *e1, const void *e2)
{
	void *key1 = ((const map_vpd_element *) e1)->key;
	void *key2 = ((const map_vpd_element *) e2)->key;

	return ((key1) < (key2) ? -1 : ((key1) == (key2) ? 0 : 1));
}

mapkit_error map_vpd_getall_sorted(map_vpd * spm, map_vpd_element ** array, mapkit_size_t * count)
{
	mapkit_error err;

	err = map_vpd_getall(spm, array, count);
	if (err)
		MAPKIT_ERROR(err);

	qsort(*array, *count, sizeof(**array), map_vpd_compare);

	return MAPKIT_OK;
}

mapkit_error map_vpd_setall(map_vpd * spm, map_vpd_element * array, mapkit_size_t count)
{
	mapkit_size_t array_iindex;
	mapkit_error err;

	err = map_vpd_ensurecapacity(spm, spm->used + count);
	if (err)
		MAPKIT_ERROR(err);

	if (spm->alwaysdefault)
		/*
		 * Prevent shrinking 
		 */
		spm->minused = 0;

	for (array_iindex = 0; array_iindex < count; array_iindex++) {
		err = map_vpd_set(spm, array[array_iindex].key, array[array_iindex].value);
		if (err)
			MAPKIT_ERROR(err);
	}

	if (spm->alwaysdefault)
		map_vpd_adjustcapacity(spm);

	return MAPKIT_OK;
}

mapkit_error map_vpd_removeall(map_vpd * spm, void **array, mapkit_size_t count)
{
	mapkit_size_t array_iindex;
	mapkit_error err;

	/*
	 * Prevent shrinking 
	 */
	spm->minused = 0;

	for (array_iindex = 0; array_iindex < count; array_iindex++) {
		err = map_vpd_remove(spm, array[array_iindex]);
		if (err)
			MAPKIT_ERROR(err);
	}

	map_vpd_adjustcapacity(spm);

	return MAPKIT_OK;
}

void map_vpd_printstats(map_vpd * spm)
{
	fprintf(stderr, "MAPKIT: map_vpd statistics\n");
	fprintf(stderr, "MAPKIT: alwaysdefault = %d\n", spm->alwaysdefault);
	fprintf(stderr, "MAPKIT: minused = %ld, maxfill = %ld\n", (long) spm->minused, (long) spm->maxfill);
	fprintf(stderr, "MAPKIT: minusedfactor = %g, maxfillfactor = %g\n", spm->minusedfactor, spm->maxfillfactor);
#ifdef MAPKIT_COLLISIONS
	fprintf(stderr, "MAPKIT: insertionindexs = %lu, collisions = %lu\n", (unsigned long) spm->insertionindexs,
		(unsigned long) spm->insertionindex_collisions);
	fprintf(stderr, "MAPKIT: keyindexs = %lu, collisions = %lu\n", (unsigned long) spm->keyindexs, (unsigned long) spm->keyindex_collisions);
#endif
}

#endif							       /* MAPKIT_map_vpd */

#ifdef MAPKIT_map_vpvp

/*
  Implementation for map_vpvp (void* -> void*)
  Default value : NULL
  Uses a state field.
*/

/* Static prototypes */

/* Return the iindex of key, or MAPKIT_KEYNOTFOUND */
static mapkit_size_t map_vpvp_keyindex(map_vpvp * spm, void *key);

/* Return the iindex of key or -(insertion iindex)-1 if key not found */
static mapkit_size_t map_vpvp_insertionindex(map_vpvp * spm, void *key);

/* Implementation */

mapkit_error map_vpvp_init(map_vpvp * spm)
{
	return map_vpvp_init_hint(spm, MAPKIT_DEFAULT_EXPECTEDUSED);
}

mapkit_error map_vpvp_init_hint(map_vpvp * spm, mapkit_size_t used)
{
#ifdef MAPKIT_DEBUG
	fprintf(stderr, "MAPKIT: init\n");
#endif

	spm->size = 0;
	spm->fill = 0;
	spm->used = 0;
	spm->maxfillfactor = 0.5;
	spm->minusedfactor = 0.2;
	spm->contents = NULL;
#ifdef MAPKIT_COLLISIONS
	spm->insertionindexs = spm->insertionindex_collisions = 0;
	spm->keyindexs = spm->keyindex_collisions = 0;
#endif
	spm->defaultvalue = NULL;
	spm->alwaysdefault = 0;

	return map_vpvp_reallocate(spm, map_vpvp_meansize(spm, used));
}

mapkit_error map_vpvp_ensurecapacity(map_vpvp * spm, mapkit_size_t used)
{
	if (used > (spm->used + spm->maxfill - spm->fill)) {
#ifdef MAPKIT_DEBUG
		fprintf(stderr, "MAPKIT: need more capacity\n");
#endif
		return map_vpvp_reallocate(spm, map_vpvp_meansize(spm, used));
	} else
		return MAPKIT_OK;
}

mapkit_error map_vpvp_adjustcapacity(map_vpvp * spm)
{
	spm->minused = (mapkit_size_t) (spm->size * spm->minusedfactor);
	spm->maxfill = (mapkit_size_t) (spm->size * spm->maxfillfactor);

	if (spm->used < spm->minused) {
#ifdef MAPKIT_DEBUG
		fprintf(stderr, "MAPKIT: used < minused\n");
#endif
		return map_vpvp_reallocate(spm, map_vpvp_meansize(spm, spm->used));
	} else if (spm->fill > spm->maxfill) {
#ifdef MAPKIT_DEBUG
		fprintf(stderr, "MAPKIT: fill > maxfill\n");
#endif
		return map_vpvp_reallocate(spm, map_vpvp_meansize(spm, spm->used));
	} else
		return MAPKIT_OK;
}

void map_vpvp_free(map_vpvp * spm)
{
	free(spm->contents);
	spm->contents = NULL;
	spm->size = 0;
	spm->fill = 0;
	spm->used = 0;
	spm->maxfill = 0;
	spm->minused = 0;
#ifdef MAPKIT_COLLISIONS
	spm->insertionindexs = spm->insertionindex_collisions = 0;
	spm->keyindexs = spm->keyindex_collisions = 0;
#endif
}

mapkit_error map_vpvp_copy(map_vpvp * to, map_vpvp * from)
{
	map_vpvp_storage *contentscopy;
	size_t size = from->size * sizeof(*from->contents);

	contentscopy = (map_vpvp_storage *) malloc(size);
	if (contentscopy == NULL)
		MAPKIT_ERROR(MAPKIT_ENOMEM);

	memcpy(to, from, sizeof(map_vpvp));
	to->contents = contentscopy;
	memcpy(to->contents, from->contents, size);

	return MAPKIT_OK;
}

mapkit_size_t map_vpvp_growsize(map_vpvp * spm, mapkit_size_t used)
{
	return (mapkit_size_t) (4.0 * used / (3 * spm->minusedfactor + spm->maxfillfactor));
}

mapkit_size_t map_vpvp_shrinksize(map_vpvp * spm, mapkit_size_t used)
{
	return (mapkit_size_t) (4.0 * used / (spm->minusedfactor + 3 * spm->maxfillfactor));
}

mapkit_size_t map_vpvp_meansize(map_vpvp * spm, mapkit_size_t used)
{
	return (mapkit_size_t) (2.0 * used / (spm->minusedfactor + spm->maxfillfactor));
}

mapkit_error map_vpvp_reallocate(map_vpvp * spm, mapkit_size_t newsize)
{
	mapkit_size_t iindex;
	mapkit_size_t oldsize;
	map_vpvp_storage *newcontents, *oldcontents;

	/*
	 * At least one free entry 
	 */
	if (newsize <= spm->used)
		newsize = spm->used + 1;
	newsize = mapkit_nextprime(newsize);
	if (newsize <= spm->used)
		MAPKIT_ERROR(MAPKIT_ETOOBIG);

	newcontents = (map_vpvp_storage *) malloc(newsize * sizeof(*spm->contents));
	if (newcontents == NULL)
		MAPKIT_ERROR(MAPKIT_ENOMEM);

	/*
	 * Initialize all entries to "free" 
	 */
	for (iindex = 0; iindex < newsize; iindex++)
		newcontents[iindex].state = MAPKIT_FREESLOT;

	oldcontents = spm->contents;
	oldsize = spm->size;
	spm->contents = newcontents;
	spm->size = newsize;

#ifdef MAPKIT_DEBUG
	fprintf(stderr, "MAPKIT: reallocate %ld -> %ld\n", (long) oldsize, (long) newsize);
#endif

	spm->maxfill = (mapkit_size_t) (newsize * spm->maxfillfactor);
	/*
	 * At least one free entry 
	 */
	if (spm->maxfill >= newsize)
		spm->maxfill = newsize - 1;
	spm->minused = (mapkit_size_t) (newsize * spm->minusedfactor);
	spm->used = 0;

	if (oldcontents != NULL) {
		int used = 0;
		void *key;
		void *defaultvalue = spm->defaultvalue;
		int notalwaysdefault = !spm->alwaysdefault;

		/*
		 * Copy all entries from old to new 
		 */
		for (iindex = 0; iindex < oldsize; iindex++)
			if (oldcontents[iindex].state == MAPKIT_FULLSLOT) {
				mapkit_size_t ins_iindex;
				map_vpvp_storage *contents;

				key = oldcontents[iindex].key;

				/*
				 * Fast path 
				 */
				ins_iindex = ((mapkit_hash_t) key) % spm->size;
				contents = &(newcontents[ins_iindex]);

				if (contents->state != MAPKIT_FREESLOT) {
					ins_iindex = map_vpvp_insertionindex(spm, key);
					contents = &(newcontents[ins_iindex]);
				}
#ifdef MAPKIT_COLLISIONS
				else
					spm->insertionindexs++;
#endif
				if (notalwaysdefault || (!((oldcontents[iindex].value) == (defaultvalue)))) {
					contents->value = oldcontents[iindex].value;
					contents->state = MAPKIT_FULLSLOT;
					contents->key = key;
					used++;
				}
			}
		free(oldcontents);
		spm->used = used;
	}
	spm->fill = spm->used;

	return MAPKIT_OK;
}

void *map_vpvp_value_s(map_vpvp * spm, void *key)
{
	mapkit_size_t iindex;

	iindex = map_vpvp_keyindex(spm, key);

	if (iindex < 0) {
		if (spm->alwaysdefault)
			return spm->defaultvalue;
		else
			MAPKIT_FATAL_ERROR(MAPKIT_EKEYNOTFOUND);
	}

	return spm->contents[iindex].value;
}

mapkit_error map_vpvp_get_s(map_vpvp * spm, void *key, void **value)
{
	mapkit_size_t iindex;

	iindex = map_vpvp_keyindex(spm, key);

	if (iindex < 0) {
		if (spm->alwaysdefault) {
			*value = spm->defaultvalue;
			return MAPKIT_OK;
		} else
			MAPKIT_ERROR(MAPKIT_EKEYNOTFOUND);
	}

	*value = spm->contents[iindex].value;
	return MAPKIT_OK;
}

mapkit_error map_vpvp_set_s(map_vpvp * spm, void *key, void *value)
{
	mapkit_size_t iindex;

	iindex = map_vpvp_insertionindex(spm, key);

	if (iindex < 0)
		/*
		 * FULLSLOT 
		 */
		spm->contents[-iindex - 1].value = value;
	else {
		map_vpvp_storage *element = &(spm->contents[iindex]);
		int ffree = element->state == MAPKIT_FREESLOT;

		element->state = MAPKIT_FULLSLOT;
		element->value = value;
		element->key = key;
		spm->used++;

		if (ffree && ((++spm->fill) > spm->maxfill)) {
#ifdef MAPKIT_DEBUG
			fprintf(stderr, "MAPKIT: fill > maxfill\n");
#endif
			return map_vpvp_reallocate(spm, map_vpvp_growsize(spm, spm->used));
		}
	}
	return MAPKIT_OK;
}

void **map_vpvp_insertptr_s(map_vpvp * spm, void *key)
{
	mapkit_size_t iindex;

	iindex = map_vpvp_insertionindex(spm, key);

	if (iindex < 0)
		return &(spm->contents[-iindex - 1].value);
	else {
		map_vpvp_storage *element = &(spm->contents[iindex]);

		if (element->state == MAPKIT_FREESLOT) {
			/*
			 * FREESLOT 
			 */
			if (spm->fill >= spm->maxfill) {
				mapkit_error err;

#ifdef MAPKIT_DEBUG
				fprintf(stderr, "MAPKIT: fill => maxfill before insert\n");
#endif
				/*
				 * Must reallocate -before- inserting defaultvalue 
				 */
				err = map_vpvp_reallocate(spm, map_vpvp_growsize(spm, spm->used + 1));
				if (err) {
					MAPKIT_ERROR_NORET(err);
					return NULL;
				}

				iindex = map_vpvp_insertionindex(spm, key);
				/*
				 * FREESLOT 
				 */
				spm->contents[iindex].state = MAPKIT_FULLSLOT;
				spm->contents[iindex].key = key;
				spm->contents[iindex].value = spm->defaultvalue;
				spm->used++;
				spm->fill++;
				return &(spm->contents[iindex].value);
			} else
				spm->fill++;
		}

		element->state = MAPKIT_FULLSLOT;
		element->key = key;
		element->value = spm->defaultvalue;
		spm->used++;
		return &(element->value);
	}
}

void **map_vpvp_ptr_s(map_vpvp * spm, void *key)
{
	mapkit_size_t iindex;

	iindex = map_vpvp_keyindex(spm, key);

	if ((iindex >= 0) && ((!spm->alwaysdefault) || (!((spm->contents[iindex].value) == (spm->defaultvalue)))))
		return &(spm->contents[iindex].value);
	else
		return NULL;
}

mapkit_error map_vpvp_remove_s(map_vpvp * spm, void *key)
{
	mapkit_size_t iindex;

	iindex = map_vpvp_keyindex(spm, key);

	if (iindex < 0) {
		if (spm->alwaysdefault)
			return MAPKIT_OK;
		else
			return MAPKIT_EKEYNOTFOUND;
	}

	spm->contents[iindex].state = MAPKIT_DELETEDSLOT;
	spm->used--;
	if (spm->used < spm->minused) {
#ifdef MAPKIT_DEBUG
		fprintf(stderr, "MAPKIT: used < minused\n");
#endif
		return map_vpvp_reallocate(spm, map_vpvp_shrinksize(spm, spm->used));
	}

	return MAPKIT_OK;
}

mapkit_size_t map_vpvp_keyindex(map_vpvp * spm, void *key)
{
	mapkit_size_t iindex, decrement;

	signed char state;

	iindex = ((mapkit_hash_t) key) % spm->size;
	decrement = (((mapkit_hash_t) key) % (spm->size - 2));
	decrement += (decrement == 0);

#ifdef MAPKIT_COLLISIONS
	spm->keyindexs++;
#endif

	while ((state = spm->contents[iindex].state) != MAPKIT_FREESLOT
	       && (state == MAPKIT_DELETEDSLOT || (!((spm->contents[iindex].key) == (key))))) {
#ifdef MAPKIT_COLLISIONS
		spm->keyindex_collisions++;
#endif
		iindex -= decrement;
		if (iindex < 0)
			iindex += spm->size;
	}

	if (state == MAPKIT_FREESLOT)
		return MAPKIT_KEYNOTFOUND;
	return iindex;
}

mapkit_size_t map_vpvp_insertionindex(map_vpvp * spm, void *key)
{
	mapkit_size_t iindex, decrement;
	signed char state;

#ifdef MAPKIT_COLLISIONS
	spm->insertionindexs++;
#endif

	iindex = ((mapkit_hash_t) key) % spm->size;

	/*
	 * Fast path (largely superfluous) 
	 */
	if ((state = spm->contents[iindex].state) == MAPKIT_FREESLOT)
		return iindex;
	if ((state == MAPKIT_FULLSLOT)
	    && ((spm->contents[iindex].key) == (key)))
		return -iindex - 1;

	decrement = (((mapkit_hash_t) key) % (spm->size - 2));
	decrement += (decrement == 0);

	while ((state == MAPKIT_FULLSLOT)
	       && (!((spm->contents[iindex].key) == (key)))) {
#ifdef MAPKIT_COLLISIONS
		spm->insertionindex_collisions++;
#endif
		iindex -= decrement;
		if (iindex < 0)
			iindex += spm->size;
		state = spm->contents[iindex].state;
	}

	if (state == MAPKIT_DELETEDSLOT) {
		mapkit_size_t index2 = iindex;

		while ((state = spm->contents[iindex].state) != MAPKIT_FREESLOT && ((state == MAPKIT_DELETEDSLOT)
										    || (!((spm->contents[iindex].key) == (key))))) {
			iindex -= decrement;
			if (iindex < 0)
				iindex += spm->size;
			state = spm->contents[iindex].state;
		}
		if (state == MAPKIT_FREESLOT)
			return index2;
	}

	if (state == MAPKIT_FULLSLOT)
		return -iindex - 1;
	return iindex;
}

mapkit_size_t map_vpvp_next(map_vpvp * spm, mapkit_size_t iindex)
{
	mapkit_size_t size = spm->size;
	map_vpvp_storage *pos_contents;
	void *defaultvalue = spm->defaultvalue;
	int notalwaysdefault = !spm->alwaysdefault;

	pos_contents = &(spm->contents[++iindex]);

	for (; iindex < size; iindex++, pos_contents++)
		if ((pos_contents->state == MAPKIT_FULLSLOT)
		    && (notalwaysdefault || (!((pos_contents->value) == (defaultvalue)))))
			return iindex;

	return -1;
}

map_vpvp_storage *map_vpvp_nextptr(map_vpvp * spm, map_vpvp_storage * pos_contents)
{
	map_vpvp_storage *end = &(spm->contents[spm->size]);
	void *defaultvalue = spm->defaultvalue;
	int notalwaysdefault = !spm->alwaysdefault;

	if (pos_contents == NULL)
		pos_contents = spm->contents;
	else {
		pos_contents++;
		if (pos_contents <= spm->contents)
			return NULL;
	}

	for (; pos_contents < end; pos_contents++)
		if ((pos_contents->state == MAPKIT_FULLSLOT)
		    && (notalwaysdefault || (!((pos_contents->value) == (defaultvalue)))))
			return pos_contents;

	return NULL;
}

mapkit_error map_vpvp_getall(map_vpvp * spm, map_vpvp_element ** array, mapkit_size_t * count)
{
	mapkit_size_t iindex;
	mapkit_size_t size = spm->size, vcount = 0;
	map_vpvp_element *pos_array;
	map_vpvp_storage *pos_contents;
	void *defaultvalue = spm->defaultvalue;
	int notalwaysdefault = !spm->alwaysdefault;

	pos_array = *array = (map_vpvp_element *) malloc(sizeof(**array) * spm->used);
	if (*array == NULL)
		MAPKIT_ERROR(MAPKIT_ENOMEM);

	pos_contents = spm->contents;

	for (iindex = 0; iindex < size; iindex++, pos_contents++)
		if ((pos_contents->state == MAPKIT_FULLSLOT)
		    && (notalwaysdefault || (!((pos_contents->value) == (defaultvalue))))) {
			pos_array->key = pos_contents->key;
			pos_array->value = pos_contents->value;
			pos_array++;
			vcount++;
		}
	*count = vcount;

	return MAPKIT_OK;
}

mapkit_error map_vpvp_clean(map_vpvp * spm)
{
	mapkit_size_t iindex, count = 0;
	mapkit_size_t size = spm->size;
	map_vpvp_storage *pos_contents;
	void *defaultvalue = spm->defaultvalue;

	pos_contents = spm->contents;
	for (iindex = 0; iindex < size; iindex++, pos_contents++)
		if ((pos_contents->state == MAPKIT_FULLSLOT)
		    && ((pos_contents->value) == (defaultvalue))) {
			pos_contents->state = MAPKIT_DELETEDSLOT;
			count++;
		}

	spm->used -= count;
	if (spm->used < spm->minused) {
#ifdef MAPKIT_DEBUG
		fprintf(stderr, "MAPKIT: used < minused\n");
#endif
		return map_vpvp_reallocate(spm, map_vpvp_meansize(spm, spm->used));
	}

	return MAPKIT_OK;
}

int map_vpvp_compare(const void *e1, const void *e2)
{
	void *key1 = ((const map_vpvp_element *) e1)->key;
	void *key2 = ((const map_vpvp_element *) e2)->key;

	return ((key1) < (key2) ? -1 : ((key1) == (key2) ? 0 : 1));
}

mapkit_error map_vpvp_getall_sorted(map_vpvp * spm, map_vpvp_element ** array, mapkit_size_t * count)
{
	mapkit_error err;

	err = map_vpvp_getall(spm, array, count);
	if (err)
		MAPKIT_ERROR(err);

	qsort(*array, *count, sizeof(**array), map_vpvp_compare);

	return MAPKIT_OK;
}

mapkit_error map_vpvp_setall(map_vpvp * spm, map_vpvp_element * array, mapkit_size_t count)
{
	mapkit_size_t array_iindex;
	mapkit_error err;

	err = map_vpvp_ensurecapacity(spm, spm->used + count);
	if (err)
		MAPKIT_ERROR(err);

	if (spm->alwaysdefault)
		/*
		 * Prevent shrinking 
		 */
		spm->minused = 0;

	for (array_iindex = 0; array_iindex < count; array_iindex++) {
		err = map_vpvp_set(spm, array[array_iindex].key, array[array_iindex].value);
		if (err)
			MAPKIT_ERROR(err);
	}

	if (spm->alwaysdefault)
		map_vpvp_adjustcapacity(spm);

	return MAPKIT_OK;
}

mapkit_error map_vpvp_removeall(map_vpvp * spm, void **array, mapkit_size_t count)
{
	mapkit_size_t array_iindex;
	mapkit_error err;

	/*
	 * Prevent shrinking 
	 */
	spm->minused = 0;

	for (array_iindex = 0; array_iindex < count; array_iindex++) {
		err = map_vpvp_remove(spm, array[array_iindex]);
		if (err)
			MAPKIT_ERROR(err);
	}

	map_vpvp_adjustcapacity(spm);

	return MAPKIT_OK;
}

void map_vpvp_printstats(map_vpvp * spm)
{
	fprintf(stderr, "MAPKIT: map_vpvp statistics\n");
	fprintf(stderr, "MAPKIT: alwaysdefault = %d\n", spm->alwaysdefault);
	fprintf(stderr, "MAPKIT: minused = %ld, maxfill = %ld\n", (long) spm->minused, (long) spm->maxfill);
	fprintf(stderr, "MAPKIT: minusedfactor = %g, maxfillfactor = %g\n", spm->minusedfactor, spm->maxfillfactor);
#ifdef MAPKIT_COLLISIONS
	fprintf(stderr, "MAPKIT: insertionindexs = %lu, collisions = %lu\n", (unsigned long) spm->insertionindexs,
		(unsigned long) spm->insertionindex_collisions);
	fprintf(stderr, "MAPKIT: keyindexs = %lu, collisions = %lu\n", (unsigned long) spm->keyindexs, (unsigned long) spm->keyindex_collisions);
#endif
}

#endif							       /* MAPKIT_map_vpvp */

#ifdef MAPKIT_map_h_vpi

/*
  Implementation for map_h_vpi (void* -> int)
  Default value : 0
  Uses a state field.
*/

/* Static prototypes */

/* Return the iindex of key, or MAPKIT_KEYNOTFOUND */
static mapkit_size_t map_h_vpi_keyindex(map_h_vpi * spm, void *key, mapkit_hash_t hash);

/* Return the iindex of key or -(insertion iindex)-1 if key not found */
static mapkit_size_t map_h_vpi_insertionindex(map_h_vpi * spm, void *key, mapkit_hash_t hash);

/* Implementation */

mapkit_error map_h_vpi_init(map_h_vpi * spm)
{
	return map_h_vpi_init_hint(spm, MAPKIT_DEFAULT_EXPECTEDUSED);
}

mapkit_error map_h_vpi_init_hint(map_h_vpi * spm, mapkit_size_t used)
{
#ifdef MAPKIT_DEBUG
	fprintf(stderr, "MAPKIT: init\n");
#endif

	spm->size = 0;
	spm->fill = 0;
	spm->used = 0;
	spm->maxfillfactor = 0.5;
	spm->minusedfactor = 0.2;
	spm->contents = NULL;
#ifdef MAPKIT_COLLISIONS
	spm->insertionindexs = spm->insertionindex_collisions = 0;
	spm->keyindexs = spm->keyindex_collisions = 0;
#endif
	spm->defaultvalue = 0;
	spm->alwaysdefault = 0;

	return map_h_vpi_reallocate(spm, map_h_vpi_meansize(spm, used));
}

mapkit_error map_h_vpi_ensurecapacity(map_h_vpi * spm, mapkit_size_t used)
{
	if (used > (spm->used + spm->maxfill - spm->fill)) {
#ifdef MAPKIT_DEBUG
		fprintf(stderr, "MAPKIT: need more capacity\n");
#endif
		return map_h_vpi_reallocate(spm, map_h_vpi_meansize(spm, used));
	} else
		return MAPKIT_OK;
}

mapkit_error map_h_vpi_adjustcapacity(map_h_vpi * spm)
{
	spm->minused = (mapkit_size_t) (spm->size * spm->minusedfactor);
	spm->maxfill = (mapkit_size_t) (spm->size * spm->maxfillfactor);

	if (spm->used < spm->minused) {
#ifdef MAPKIT_DEBUG
		fprintf(stderr, "MAPKIT: used < minused\n");
#endif
		return map_h_vpi_reallocate(spm, map_h_vpi_meansize(spm, spm->used));
	} else if (spm->fill > spm->maxfill) {
#ifdef MAPKIT_DEBUG
		fprintf(stderr, "MAPKIT: fill > maxfill\n");
#endif
		return map_h_vpi_reallocate(spm, map_h_vpi_meansize(spm, spm->used));
	} else
		return MAPKIT_OK;
}

void map_h_vpi_free(map_h_vpi * spm)
{
	free(spm->contents);
	spm->contents = NULL;
	spm->size = 0;
	spm->fill = 0;
	spm->used = 0;
	spm->maxfill = 0;
	spm->minused = 0;
#ifdef MAPKIT_COLLISIONS
	spm->insertionindexs = spm->insertionindex_collisions = 0;
	spm->keyindexs = spm->keyindex_collisions = 0;
#endif
}

mapkit_error map_h_vpi_copy(map_h_vpi * to, map_h_vpi * from)
{
	map_h_vpi_storage *contentscopy;
	size_t size = from->size * sizeof(*from->contents);

	contentscopy = (map_h_vpi_storage *) malloc(size);
	if (contentscopy == NULL)
		MAPKIT_ERROR(MAPKIT_ENOMEM);

	memcpy(to, from, sizeof(map_h_vpi));
	to->contents = contentscopy;
	memcpy(to->contents, from->contents, size);

	return MAPKIT_OK;
}

mapkit_size_t map_h_vpi_growsize(map_h_vpi * spm, mapkit_size_t used)
{
	return (mapkit_size_t) (4.0 * used / (3 * spm->minusedfactor + spm->maxfillfactor));
}

mapkit_size_t map_h_vpi_shrinksize(map_h_vpi * spm, mapkit_size_t used)
{
	return (mapkit_size_t) (4.0 * used / (spm->minusedfactor + 3 * spm->maxfillfactor));
}

mapkit_size_t map_h_vpi_meansize(map_h_vpi * spm, mapkit_size_t used)
{
	return (mapkit_size_t) (2.0 * used / (spm->minusedfactor + spm->maxfillfactor));
}

mapkit_error map_h_vpi_reallocate(map_h_vpi * spm, mapkit_size_t newsize)
{
	mapkit_size_t iindex;
	mapkit_size_t oldsize;
	map_h_vpi_storage *newcontents, *oldcontents;

	/*
	 * At least one free entry 
	 */
	if (newsize <= spm->used)
		newsize = spm->used + 1;
	newsize = mapkit_nextprime(newsize);
	if (newsize <= spm->used)
		MAPKIT_ERROR(MAPKIT_ETOOBIG);

	newcontents = (map_h_vpi_storage *) malloc(newsize * sizeof(*spm->contents));
	if (newcontents == NULL)
		MAPKIT_ERROR(MAPKIT_ENOMEM);

	/*
	 * Initialize all entries to "free" 
	 */
	for (iindex = 0; iindex < newsize; iindex++)
		newcontents[iindex].state = MAPKIT_FREESLOT;

	oldcontents = spm->contents;
	oldsize = spm->size;
	spm->contents = newcontents;
	spm->size = newsize;

#ifdef MAPKIT_DEBUG
	fprintf(stderr, "MAPKIT: reallocate %ld -> %ld\n", (long) oldsize, (long) newsize);
#endif

	spm->maxfill = (mapkit_size_t) (newsize * spm->maxfillfactor);
	/*
	 * At least one free entry 
	 */
	if (spm->maxfill >= newsize)
		spm->maxfill = newsize - 1;
	spm->minused = (mapkit_size_t) (newsize * spm->minusedfactor);
	spm->used = 0;

	if (oldcontents != NULL) {
		int used = 0;
		void *key;
		int defaultvalue = spm->defaultvalue;
		int notalwaysdefault = !spm->alwaysdefault;

		/*
		 * Copy all entries from old to new 
		 */
		for (iindex = 0; iindex < oldsize; iindex++)
			if (oldcontents[iindex].state == MAPKIT_FULLSLOT) {
				mapkit_size_t ins_iindex;
				mapkit_hash_t hash;
				map_h_vpi_storage *contents;

				key = oldcontents[iindex].key;

				/*
				 * Fast path 
				 */
				ins_iindex = (hash = mapkit_hash((mapkit_hash_t) key)) % spm->size;
				contents = &(newcontents[ins_iindex]);

				if (contents->state != MAPKIT_FREESLOT) {
					ins_iindex = map_h_vpi_insertionindex(spm, key, hash);
					contents = &(newcontents[ins_iindex]);
				}
#ifdef MAPKIT_COLLISIONS
				else
					spm->insertionindexs++;
#endif
				if (notalwaysdefault || (!((oldcontents[iindex].value) == (defaultvalue)))) {
					contents->value = oldcontents[iindex].value;
					contents->state = MAPKIT_FULLSLOT;
					contents->key = key;
					used++;
				}
			}
		free(oldcontents);
		spm->used = used;
	}
	spm->fill = spm->used;

	return MAPKIT_OK;
}

int map_h_vpi_value_s(map_h_vpi * spm, void *key, mapkit_hash_t hash)
{
	mapkit_size_t iindex;

	iindex = map_h_vpi_keyindex(spm, key, hash);

	if (iindex < 0) {
		if (spm->alwaysdefault)
			return spm->defaultvalue;
		else
			MAPKIT_FATAL_ERROR(MAPKIT_EKEYNOTFOUND);
	}

	return spm->contents[iindex].value;
}

mapkit_error map_h_vpi_get_s(map_h_vpi * spm, void *key, int *value, mapkit_hash_t hash)
{
	mapkit_size_t iindex;

	iindex = map_h_vpi_keyindex(spm, key, hash);

	if (iindex < 0) {
		if (spm->alwaysdefault) {
			*value = spm->defaultvalue;
			return MAPKIT_OK;
		} else
			MAPKIT_ERROR(MAPKIT_EKEYNOTFOUND);
	}

	*value = spm->contents[iindex].value;
	return MAPKIT_OK;
}

mapkit_error map_h_vpi_set_s(map_h_vpi * spm, void *key, int value, mapkit_hash_t hash)
{
	mapkit_size_t iindex;

	iindex = map_h_vpi_insertionindex(spm, key, hash);

	if (iindex < 0)
		/*
		 * FULLSLOT 
		 */
		spm->contents[-iindex - 1].value = value;
	else {
		map_h_vpi_storage *element = &(spm->contents[iindex]);
		int ffree = element->state == MAPKIT_FREESLOT;

		element->state = MAPKIT_FULLSLOT;
		element->value = value;
		element->key = key;
		spm->used++;

		if (ffree && ((++spm->fill) > spm->maxfill)) {
#ifdef MAPKIT_DEBUG
			fprintf(stderr, "MAPKIT: fill > maxfill\n");
#endif
			return map_h_vpi_reallocate(spm, map_h_vpi_growsize(spm, spm->used));
		}
	}
	return MAPKIT_OK;
}

int *map_h_vpi_insertptr_s(map_h_vpi * spm, void *key, mapkit_hash_t hash)
{
	mapkit_size_t iindex;

	iindex = map_h_vpi_insertionindex(spm, key, hash);

	if (iindex < 0)
		return &(spm->contents[-iindex - 1].value);
	else {
		map_h_vpi_storage *element = &(spm->contents[iindex]);

		if (element->state == MAPKIT_FREESLOT) {
			/*
			 * FREESLOT 
			 */
			if (spm->fill >= spm->maxfill) {
				mapkit_error err;

#ifdef MAPKIT_DEBUG
				fprintf(stderr, "MAPKIT: fill => maxfill before insert\n");
#endif
				/*
				 * Must reallocate -before- inserting defaultvalue 
				 */
				err = map_h_vpi_reallocate(spm, map_h_vpi_growsize(spm, spm->used + 1));
				if (err) {
					MAPKIT_ERROR_NORET(err);
					return NULL;
				}

				iindex = map_h_vpi_insertionindex(spm, key, hash);
				/*
				 * FREESLOT 
				 */
				spm->contents[iindex].state = MAPKIT_FULLSLOT;
				spm->contents[iindex].key = key;
				spm->contents[iindex].value = spm->defaultvalue;
				spm->used++;
				spm->fill++;
				return &(spm->contents[iindex].value);
			} else
				spm->fill++;
		}

		element->state = MAPKIT_FULLSLOT;
		element->key = key;
		element->value = spm->defaultvalue;
		spm->used++;
		return &(element->value);
	}
}

int *map_h_vpi_ptr_s(map_h_vpi * spm, void *key, mapkit_hash_t hash)
{
	mapkit_size_t iindex;

	iindex = map_h_vpi_keyindex(spm, key, hash);

	if ((iindex >= 0) && ((!spm->alwaysdefault) || (!((spm->contents[iindex].value) == (spm->defaultvalue)))))
		return &(spm->contents[iindex].value);
	else
		return NULL;
}

mapkit_error map_h_vpi_remove_s(map_h_vpi * spm, void *key, mapkit_hash_t hash)
{
	mapkit_size_t iindex;

	iindex = map_h_vpi_keyindex(spm, key, hash);

	if (iindex < 0) {
		if (spm->alwaysdefault)
			return MAPKIT_OK;
		else
			return MAPKIT_EKEYNOTFOUND;
	}

	spm->contents[iindex].state = MAPKIT_DELETEDSLOT;
	spm->used--;
	if (spm->used < spm->minused) {
#ifdef MAPKIT_DEBUG
		fprintf(stderr, "MAPKIT: used < minused\n");
#endif
		return map_h_vpi_reallocate(spm, map_h_vpi_shrinksize(spm, spm->used));
	}

	return MAPKIT_OK;
}

mapkit_size_t map_h_vpi_keyindex(map_h_vpi * spm, void *key, mapkit_hash_t hash)
{
	mapkit_size_t iindex, decrement;

	signed char state;

	iindex = hash % spm->size;
	decrement = (hash % (spm->size - 2));
	decrement += (decrement == 0);

#ifdef MAPKIT_COLLISIONS
	spm->keyindexs++;
#endif

	while ((state = spm->contents[iindex].state) != MAPKIT_FREESLOT
	       && (state == MAPKIT_DELETEDSLOT || (!((spm->contents[iindex].key) == (key))))) {
#ifdef MAPKIT_COLLISIONS
		spm->keyindex_collisions++;
#endif
		iindex -= decrement;
		if (iindex < 0)
			iindex += spm->size;
	}

	if (state == MAPKIT_FREESLOT)
		return MAPKIT_KEYNOTFOUND;
	return iindex;
}

mapkit_size_t map_h_vpi_insertionindex(map_h_vpi * spm, void *key, mapkit_hash_t hash)
{
	mapkit_size_t iindex, decrement;
	signed char state;

#ifdef MAPKIT_COLLISIONS
	spm->insertionindexs++;
#endif

	iindex = hash % spm->size;

	/*
	 * Fast path (largely superfluous) 
	 */
	if ((state = spm->contents[iindex].state) == MAPKIT_FREESLOT)
		return iindex;
	if ((state == MAPKIT_FULLSLOT)
	    && ((spm->contents[iindex].key) == (key)))
		return -iindex - 1;

	decrement = (hash % (spm->size - 2));
	decrement += (decrement == 0);

	while ((state == MAPKIT_FULLSLOT)
	       && (!((spm->contents[iindex].key) == (key)))) {
#ifdef MAPKIT_COLLISIONS
		spm->insertionindex_collisions++;
#endif
		iindex -= decrement;
		if (iindex < 0)
			iindex += spm->size;
		state = spm->contents[iindex].state;
	}

	if (state == MAPKIT_DELETEDSLOT) {
		mapkit_size_t index2 = iindex;

		while ((state = spm->contents[iindex].state) != MAPKIT_FREESLOT && ((state == MAPKIT_DELETEDSLOT)
										    || (!((spm->contents[iindex].key) == (key))))) {
			iindex -= decrement;
			if (iindex < 0)
				iindex += spm->size;
			state = spm->contents[iindex].state;
		}
		if (state == MAPKIT_FREESLOT)
			return index2;
	}

	if (state == MAPKIT_FULLSLOT)
		return -iindex - 1;
	return iindex;
}

mapkit_size_t map_h_vpi_next(map_h_vpi * spm, mapkit_size_t iindex)
{
	mapkit_size_t size = spm->size;
	map_h_vpi_storage *pos_contents;
	int defaultvalue = spm->defaultvalue;
	int notalwaysdefault = !spm->alwaysdefault;

	pos_contents = &(spm->contents[++iindex]);

	for (; iindex < size; iindex++, pos_contents++)
		if ((pos_contents->state == MAPKIT_FULLSLOT)
		    && (notalwaysdefault || (!((pos_contents->value) == (defaultvalue)))))
			return iindex;

	return -1;
}

map_h_vpi_storage *map_h_vpi_nextptr(map_h_vpi * spm, map_h_vpi_storage * pos_contents)
{
	map_h_vpi_storage *end = &(spm->contents[spm->size]);
	int defaultvalue = spm->defaultvalue;
	int notalwaysdefault = !spm->alwaysdefault;

	if (pos_contents == NULL)
		pos_contents = spm->contents;
	else {
		pos_contents++;
		if (pos_contents <= spm->contents)
			return NULL;
	}

	for (; pos_contents < end; pos_contents++)
		if ((pos_contents->state == MAPKIT_FULLSLOT)
		    && (notalwaysdefault || (!((pos_contents->value) == (defaultvalue)))))
			return pos_contents;

	return NULL;
}

mapkit_error map_h_vpi_getall(map_h_vpi * spm, map_h_vpi_element ** array, mapkit_size_t * count)
{
	mapkit_size_t iindex;
	mapkit_size_t size = spm->size, vcount = 0;
	map_h_vpi_element *pos_array;
	map_h_vpi_storage *pos_contents;
	int defaultvalue = spm->defaultvalue;
	int notalwaysdefault = !spm->alwaysdefault;

	pos_array = *array = (map_h_vpi_element *) malloc(sizeof(**array) * spm->used);
	if (*array == NULL)
		MAPKIT_ERROR(MAPKIT_ENOMEM);

	pos_contents = spm->contents;

	for (iindex = 0; iindex < size; iindex++, pos_contents++)
		if ((pos_contents->state == MAPKIT_FULLSLOT)
		    && (notalwaysdefault || (!((pos_contents->value) == (defaultvalue))))) {
			pos_array->key = pos_contents->key;
			pos_array->value = pos_contents->value;
			pos_array++;
			vcount++;
		}
	*count = vcount;

	return MAPKIT_OK;
}

mapkit_error map_h_vpi_clean(map_h_vpi * spm)
{
	mapkit_size_t iindex, count = 0;
	mapkit_size_t size = spm->size;
	map_h_vpi_storage *pos_contents;
	int defaultvalue = spm->defaultvalue;

	pos_contents = spm->contents;
	for (iindex = 0; iindex < size; iindex++, pos_contents++)
		if ((pos_contents->state == MAPKIT_FULLSLOT)
		    && ((pos_contents->value) == (defaultvalue))) {
			pos_contents->state = MAPKIT_DELETEDSLOT;
			count++;
		}

	spm->used -= count;
	if (spm->used < spm->minused) {
#ifdef MAPKIT_DEBUG
		fprintf(stderr, "MAPKIT: used < minused\n");
#endif
		return map_h_vpi_reallocate(spm, map_h_vpi_meansize(spm, spm->used));
	}

	return MAPKIT_OK;
}

int map_h_vpi_compare(const void *e1, const void *e2)
{
	void *key1 = ((const map_h_vpi_element *) e1)->key;
	void *key2 = ((const map_h_vpi_element *) e2)->key;

	return ((key1) < (key2) ? -1 : ((key1) == (key2) ? 0 : 1));
}

mapkit_error map_h_vpi_getall_sorted(map_h_vpi * spm, map_h_vpi_element ** array, mapkit_size_t * count)
{
	mapkit_error err;

	err = map_h_vpi_getall(spm, array, count);
	if (err)
		MAPKIT_ERROR(err);

	qsort(*array, *count, sizeof(**array), map_h_vpi_compare);

	return MAPKIT_OK;
}

mapkit_error map_h_vpi_setall(map_h_vpi * spm, map_h_vpi_element * array, mapkit_size_t count)
{
	mapkit_size_t array_iindex;
	mapkit_error err;

	err = map_h_vpi_ensurecapacity(spm, spm->used + count);
	if (err)
		MAPKIT_ERROR(err);

	if (spm->alwaysdefault)
		/*
		 * Prevent shrinking 
		 */
		spm->minused = 0;

	for (array_iindex = 0; array_iindex < count; array_iindex++) {
		err = map_h_vpi_set(spm, array[array_iindex].key, array[array_iindex].value);
		if (err)
			MAPKIT_ERROR(err);
	}

	if (spm->alwaysdefault)
		map_h_vpi_adjustcapacity(spm);

	return MAPKIT_OK;
}

mapkit_error map_h_vpi_removeall(map_h_vpi * spm, void **array, mapkit_size_t count)
{
	mapkit_size_t array_iindex;
	mapkit_error err;

	/*
	 * Prevent shrinking 
	 */
	spm->minused = 0;

	for (array_iindex = 0; array_iindex < count; array_iindex++) {
		err = map_h_vpi_remove(spm, array[array_iindex]);
		if (err)
			MAPKIT_ERROR(err);
	}

	map_h_vpi_adjustcapacity(spm);

	return MAPKIT_OK;
}

void map_h_vpi_printstats(map_h_vpi * spm)
{
	fprintf(stderr, "MAPKIT: map_h_vpi statistics\n");
	fprintf(stderr, "MAPKIT: alwaysdefault = %d\n", spm->alwaysdefault);
	fprintf(stderr, "MAPKIT: minused = %ld, maxfill = %ld\n", (long) spm->minused, (long) spm->maxfill);
	fprintf(stderr, "MAPKIT: minusedfactor = %g, maxfillfactor = %g\n", spm->minusedfactor, spm->maxfillfactor);
#ifdef MAPKIT_COLLISIONS
	fprintf(stderr, "MAPKIT: insertionindexs = %lu, collisions = %lu\n", (unsigned long) spm->insertionindexs,
		(unsigned long) spm->insertionindex_collisions);
	fprintf(stderr, "MAPKIT: keyindexs = %lu, collisions = %lu\n", (unsigned long) spm->keyindexs, (unsigned long) spm->keyindex_collisions);
#endif
}

#endif							       /* MAPKIT_map_h_vpi */

#ifdef MAPKIT_map_h_vpd

/*
  Implementation for map_h_vpd (void* -> double)
  Default value : 0.0
  Uses a state field.
*/

/* Static prototypes */

/* Return the iindex of key, or MAPKIT_KEYNOTFOUND */
static mapkit_size_t map_h_vpd_keyindex(map_h_vpd * spm, void *key, mapkit_hash_t hash);

/* Return the iindex of key or -(insertion iindex)-1 if key not found */
static mapkit_size_t map_h_vpd_insertionindex(map_h_vpd * spm, void *key, mapkit_hash_t hash);

/* Implementation */

mapkit_error map_h_vpd_init(map_h_vpd * spm)
{
	return map_h_vpd_init_hint(spm, MAPKIT_DEFAULT_EXPECTEDUSED);
}

mapkit_error map_h_vpd_init_hint(map_h_vpd * spm, mapkit_size_t used)
{
#ifdef MAPKIT_DEBUG
	fprintf(stderr, "MAPKIT: init\n");
#endif

	spm->size = 0;
	spm->fill = 0;
	spm->used = 0;
	spm->maxfillfactor = 0.5;
	spm->minusedfactor = 0.2;
	spm->contents = NULL;
#ifdef MAPKIT_COLLISIONS
	spm->insertionindexs = spm->insertionindex_collisions = 0;
	spm->keyindexs = spm->keyindex_collisions = 0;
#endif
	spm->defaultvalue = 0.0;
	spm->alwaysdefault = 0;

	return map_h_vpd_reallocate(spm, map_h_vpd_meansize(spm, used));
}

mapkit_error map_h_vpd_ensurecapacity(map_h_vpd * spm, mapkit_size_t used)
{
	if (used > (spm->used + spm->maxfill - spm->fill)) {
#ifdef MAPKIT_DEBUG
		fprintf(stderr, "MAPKIT: need more capacity\n");
#endif
		return map_h_vpd_reallocate(spm, map_h_vpd_meansize(spm, used));
	} else
		return MAPKIT_OK;
}

mapkit_error map_h_vpd_adjustcapacity(map_h_vpd * spm)
{
	spm->minused = (mapkit_size_t) (spm->size * spm->minusedfactor);
	spm->maxfill = (mapkit_size_t) (spm->size * spm->maxfillfactor);

	if (spm->used < spm->minused) {
#ifdef MAPKIT_DEBUG
		fprintf(stderr, "MAPKIT: used < minused\n");
#endif
		return map_h_vpd_reallocate(spm, map_h_vpd_meansize(spm, spm->used));
	} else if (spm->fill > spm->maxfill) {
#ifdef MAPKIT_DEBUG
		fprintf(stderr, "MAPKIT: fill > maxfill\n");
#endif
		return map_h_vpd_reallocate(spm, map_h_vpd_meansize(spm, spm->used));
	} else
		return MAPKIT_OK;
}

void map_h_vpd_free(map_h_vpd * spm)
{
	free(spm->contents);
	spm->contents = NULL;
	spm->size = 0;
	spm->fill = 0;
	spm->used = 0;
	spm->maxfill = 0;
	spm->minused = 0;
#ifdef MAPKIT_COLLISIONS
	spm->insertionindexs = spm->insertionindex_collisions = 0;
	spm->keyindexs = spm->keyindex_collisions = 0;
#endif
}

mapkit_error map_h_vpd_copy(map_h_vpd * to, map_h_vpd * from)
{
	map_h_vpd_storage *contentscopy;
	size_t size = from->size * sizeof(*from->contents);

	contentscopy = (map_h_vpd_storage *) malloc(size);
	if (contentscopy == NULL)
		MAPKIT_ERROR(MAPKIT_ENOMEM);

	memcpy(to, from, sizeof(map_h_vpd));
	to->contents = contentscopy;
	memcpy(to->contents, from->contents, size);

	return MAPKIT_OK;
}

mapkit_size_t map_h_vpd_growsize(map_h_vpd * spm, mapkit_size_t used)
{
	return (mapkit_size_t) (4.0 * used / (3 * spm->minusedfactor + spm->maxfillfactor));
}

mapkit_size_t map_h_vpd_shrinksize(map_h_vpd * spm, mapkit_size_t used)
{
	return (mapkit_size_t) (4.0 * used / (spm->minusedfactor + 3 * spm->maxfillfactor));
}

mapkit_size_t map_h_vpd_meansize(map_h_vpd * spm, mapkit_size_t used)
{
	return (mapkit_size_t) (2.0 * used / (spm->minusedfactor + spm->maxfillfactor));
}

mapkit_error map_h_vpd_reallocate(map_h_vpd * spm, mapkit_size_t newsize)
{
	mapkit_size_t iindex;
	mapkit_size_t oldsize;
	map_h_vpd_storage *newcontents, *oldcontents;

	/*
	 * At least one free entry 
	 */
	if (newsize <= spm->used)
		newsize = spm->used + 1;
	newsize = mapkit_nextprime(newsize);
	if (newsize <= spm->used)
		MAPKIT_ERROR(MAPKIT_ETOOBIG);

	newcontents = (map_h_vpd_storage *) malloc(newsize * sizeof(*spm->contents));
	if (newcontents == NULL)
		MAPKIT_ERROR(MAPKIT_ENOMEM);

	/*
	 * Initialize all entries to "free" 
	 */
	for (iindex = 0; iindex < newsize; iindex++)
		newcontents[iindex].state = MAPKIT_FREESLOT;

	oldcontents = spm->contents;
	oldsize = spm->size;
	spm->contents = newcontents;
	spm->size = newsize;

#ifdef MAPKIT_DEBUG
	fprintf(stderr, "MAPKIT: reallocate %ld -> %ld\n", (long) oldsize, (long) newsize);
#endif

	spm->maxfill = (mapkit_size_t) (newsize * spm->maxfillfactor);
	/*
	 * At least one free entry 
	 */
	if (spm->maxfill >= newsize)
		spm->maxfill = newsize - 1;
	spm->minused = (mapkit_size_t) (newsize * spm->minusedfactor);
	spm->used = 0;

	if (oldcontents != NULL) {
		int used = 0;
		void *key;
		double defaultvalue = spm->defaultvalue;
		int notalwaysdefault = !spm->alwaysdefault;

		/*
		 * Copy all entries from old to new 
		 */
		for (iindex = 0; iindex < oldsize; iindex++)
			if (oldcontents[iindex].state == MAPKIT_FULLSLOT) {
				mapkit_size_t ins_iindex;
				mapkit_hash_t hash;
				map_h_vpd_storage *contents;

				key = oldcontents[iindex].key;

				/*
				 * Fast path 
				 */
				ins_iindex = (hash = mapkit_hash((mapkit_hash_t) key)) % spm->size;
				contents = &(newcontents[ins_iindex]);

				if (contents->state != MAPKIT_FREESLOT) {
					ins_iindex = map_h_vpd_insertionindex(spm, key, hash);
					contents = &(newcontents[ins_iindex]);
				}
#ifdef MAPKIT_COLLISIONS
				else
					spm->insertionindexs++;
#endif
				if (notalwaysdefault || (!((oldcontents[iindex].value) == (defaultvalue)))) {
					contents->value = oldcontents[iindex].value;
					contents->state = MAPKIT_FULLSLOT;
					contents->key = key;
					used++;
				}
			}
		free(oldcontents);
		spm->used = used;
	}
	spm->fill = spm->used;

	return MAPKIT_OK;
}

double map_h_vpd_value_s(map_h_vpd * spm, void *key, mapkit_hash_t hash)
{
	mapkit_size_t iindex;

	iindex = map_h_vpd_keyindex(spm, key, hash);

	if (iindex < 0) {
		if (spm->alwaysdefault)
			return spm->defaultvalue;
		else
			MAPKIT_FATAL_ERROR(MAPKIT_EKEYNOTFOUND);
	}

	return spm->contents[iindex].value;
}

mapkit_error map_h_vpd_get_s(map_h_vpd * spm, void *key, double *value, mapkit_hash_t hash)
{
	mapkit_size_t iindex;

	iindex = map_h_vpd_keyindex(spm, key, hash);

	if (iindex < 0) {
		if (spm->alwaysdefault) {
			*value = spm->defaultvalue;
			return MAPKIT_OK;
		} else
			MAPKIT_ERROR(MAPKIT_EKEYNOTFOUND);
	}

	*value = spm->contents[iindex].value;
	return MAPKIT_OK;
}

mapkit_error map_h_vpd_set_s(map_h_vpd * spm, void *key, double value, mapkit_hash_t hash)
{
	mapkit_size_t iindex;

	iindex = map_h_vpd_insertionindex(spm, key, hash);

	if (iindex < 0)
		/*
		 * FULLSLOT 
		 */
		spm->contents[-iindex - 1].value = value;
	else {
		map_h_vpd_storage *element = &(spm->contents[iindex]);
		int ffree = element->state == MAPKIT_FREESLOT;

		element->state = MAPKIT_FULLSLOT;
		element->value = value;
		element->key = key;
		spm->used++;

		if (ffree && ((++spm->fill) > spm->maxfill)) {
#ifdef MAPKIT_DEBUG
			fprintf(stderr, "MAPKIT: fill > maxfill\n");
#endif
			return map_h_vpd_reallocate(spm, map_h_vpd_growsize(spm, spm->used));
		}
	}
	return MAPKIT_OK;
}

double *map_h_vpd_insertptr_s(map_h_vpd * spm, void *key, mapkit_hash_t hash)
{
	mapkit_size_t iindex;

	iindex = map_h_vpd_insertionindex(spm, key, hash);

	if (iindex < 0)
		return &(spm->contents[-iindex - 1].value);
	else {
		map_h_vpd_storage *element = &(spm->contents[iindex]);

		if (element->state == MAPKIT_FREESLOT) {
			/*
			 * FREESLOT 
			 */
			if (spm->fill >= spm->maxfill) {
				mapkit_error err;

#ifdef MAPKIT_DEBUG
				fprintf(stderr, "MAPKIT: fill => maxfill before insert\n");
#endif
				/*
				 * Must reallocate -before- inserting defaultvalue 
				 */
				err = map_h_vpd_reallocate(spm, map_h_vpd_growsize(spm, spm->used + 1));
				if (err) {
					MAPKIT_ERROR_NORET(err);
					return NULL;
				}

				iindex = map_h_vpd_insertionindex(spm, key, hash);
				/*
				 * FREESLOT 
				 */
				spm->contents[iindex].state = MAPKIT_FULLSLOT;
				spm->contents[iindex].key = key;
				spm->contents[iindex].value = spm->defaultvalue;
				spm->used++;
				spm->fill++;
				return &(spm->contents[iindex].value);
			} else
				spm->fill++;
		}

		element->state = MAPKIT_FULLSLOT;
		element->key = key;
		element->value = spm->defaultvalue;
		spm->used++;
		return &(element->value);
	}
}

double *map_h_vpd_ptr_s(map_h_vpd * spm, void *key, mapkit_hash_t hash)
{
	mapkit_size_t iindex;

	iindex = map_h_vpd_keyindex(spm, key, hash);

	if ((iindex >= 0) && ((!spm->alwaysdefault) || (!((spm->contents[iindex].value) == (spm->defaultvalue)))))
		return &(spm->contents[iindex].value);
	else
		return NULL;
}

mapkit_error map_h_vpd_remove_s(map_h_vpd * spm, void *key, mapkit_hash_t hash)
{
	mapkit_size_t iindex;

	iindex = map_h_vpd_keyindex(spm, key, hash);

	if (iindex < 0) {
		if (spm->alwaysdefault)
			return MAPKIT_OK;
		else
			return MAPKIT_EKEYNOTFOUND;
	}

	spm->contents[iindex].state = MAPKIT_DELETEDSLOT;
	spm->used--;
	if (spm->used < spm->minused) {
#ifdef MAPKIT_DEBUG
		fprintf(stderr, "MAPKIT: used < minused\n");
#endif
		return map_h_vpd_reallocate(spm, map_h_vpd_shrinksize(spm, spm->used));
	}

	return MAPKIT_OK;
}

mapkit_size_t map_h_vpd_keyindex(map_h_vpd * spm, void *key, mapkit_hash_t hash)
{
	mapkit_size_t iindex, decrement;

	signed char state;

	iindex = hash % spm->size;
	decrement = (hash % (spm->size - 2));
	decrement += (decrement == 0);

#ifdef MAPKIT_COLLISIONS
	spm->keyindexs++;
#endif

	while ((state = spm->contents[iindex].state) != MAPKIT_FREESLOT
	       && (state == MAPKIT_DELETEDSLOT || (!((spm->contents[iindex].key) == (key))))) {
#ifdef MAPKIT_COLLISIONS
		spm->keyindex_collisions++;
#endif
		iindex -= decrement;
		if (iindex < 0)
			iindex += spm->size;
	}

	if (state == MAPKIT_FREESLOT)
		return MAPKIT_KEYNOTFOUND;
	return iindex;
}

mapkit_size_t map_h_vpd_insertionindex(map_h_vpd * spm, void *key, mapkit_hash_t hash)
{
	mapkit_size_t iindex, decrement;
	signed char state;

#ifdef MAPKIT_COLLISIONS
	spm->insertionindexs++;
#endif

	iindex = hash % spm->size;

	/*
	 * Fast path (largely superfluous) 
	 */
	if ((state = spm->contents[iindex].state) == MAPKIT_FREESLOT)
		return iindex;
	if ((state == MAPKIT_FULLSLOT)
	    && ((spm->contents[iindex].key) == (key)))
		return -iindex - 1;

	decrement = (hash % (spm->size - 2));
	decrement += (decrement == 0);

	while ((state == MAPKIT_FULLSLOT)
	       && (!((spm->contents[iindex].key) == (key)))) {
#ifdef MAPKIT_COLLISIONS
		spm->insertionindex_collisions++;
#endif
		iindex -= decrement;
		if (iindex < 0)
			iindex += spm->size;
		state = spm->contents[iindex].state;
	}

	if (state == MAPKIT_DELETEDSLOT) {
		mapkit_size_t index2 = iindex;

		while ((state = spm->contents[iindex].state) != MAPKIT_FREESLOT && ((state == MAPKIT_DELETEDSLOT)
										    || (!((spm->contents[iindex].key) == (key))))) {
			iindex -= decrement;
			if (iindex < 0)
				iindex += spm->size;
			state = spm->contents[iindex].state;
		}
		if (state == MAPKIT_FREESLOT)
			return index2;
	}

	if (state == MAPKIT_FULLSLOT)
		return -iindex - 1;
	return iindex;
}

mapkit_size_t map_h_vpd_next(map_h_vpd * spm, mapkit_size_t iindex)
{
	mapkit_size_t size = spm->size;
	map_h_vpd_storage *pos_contents;
	double defaultvalue = spm->defaultvalue;
	int notalwaysdefault = !spm->alwaysdefault;

	pos_contents = &(spm->contents[++iindex]);

	for (; iindex < size; iindex++, pos_contents++)
		if ((pos_contents->state == MAPKIT_FULLSLOT)
		    && (notalwaysdefault || (!((pos_contents->value) == (defaultvalue)))))
			return iindex;

	return -1;
}

map_h_vpd_storage *map_h_vpd_nextptr(map_h_vpd * spm, map_h_vpd_storage * pos_contents)
{
	map_h_vpd_storage *end = &(spm->contents[spm->size]);
	double defaultvalue = spm->defaultvalue;
	int notalwaysdefault = !spm->alwaysdefault;

	if (pos_contents == NULL)
		pos_contents = spm->contents;
	else {
		pos_contents++;
		if (pos_contents <= spm->contents)
			return NULL;
	}

	for (; pos_contents < end; pos_contents++)
		if ((pos_contents->state == MAPKIT_FULLSLOT)
		    && (notalwaysdefault || (!((pos_contents->value) == (defaultvalue)))))
			return pos_contents;

	return NULL;
}

mapkit_error map_h_vpd_getall(map_h_vpd * spm, map_h_vpd_element ** array, mapkit_size_t * count)
{
	mapkit_size_t iindex;
	mapkit_size_t size = spm->size, vcount = 0;
	map_h_vpd_element *pos_array;
	map_h_vpd_storage *pos_contents;
	double defaultvalue = spm->defaultvalue;
	int notalwaysdefault = !spm->alwaysdefault;

	pos_array = *array = (map_h_vpd_element *) malloc(sizeof(**array) * spm->used);
	if (*array == NULL)
		MAPKIT_ERROR(MAPKIT_ENOMEM);

	pos_contents = spm->contents;

	for (iindex = 0; iindex < size; iindex++, pos_contents++)
		if ((pos_contents->state == MAPKIT_FULLSLOT)
		    && (notalwaysdefault || (!((pos_contents->value) == (defaultvalue))))) {
			pos_array->key = pos_contents->key;
			pos_array->value = pos_contents->value;
			pos_array++;
			vcount++;
		}
	*count = vcount;

	return MAPKIT_OK;
}

mapkit_error map_h_vpd_clean(map_h_vpd * spm)
{
	mapkit_size_t iindex, count = 0;
	mapkit_size_t size = spm->size;
	map_h_vpd_storage *pos_contents;
	double defaultvalue = spm->defaultvalue;

	pos_contents = spm->contents;
	for (iindex = 0; iindex < size; iindex++, pos_contents++)
		if ((pos_contents->state == MAPKIT_FULLSLOT)
		    && ((pos_contents->value) == (defaultvalue))) {
			pos_contents->state = MAPKIT_DELETEDSLOT;
			count++;
		}

	spm->used -= count;
	if (spm->used < spm->minused) {
#ifdef MAPKIT_DEBUG
		fprintf(stderr, "MAPKIT: used < minused\n");
#endif
		return map_h_vpd_reallocate(spm, map_h_vpd_meansize(spm, spm->used));
	}

	return MAPKIT_OK;
}

int map_h_vpd_compare(const void *e1, const void *e2)
{
	void *key1 = ((const map_h_vpd_element *) e1)->key;
	void *key2 = ((const map_h_vpd_element *) e2)->key;

	return ((key1) < (key2) ? -1 : ((key1) == (key2) ? 0 : 1));
}

mapkit_error map_h_vpd_getall_sorted(map_h_vpd * spm, map_h_vpd_element ** array, mapkit_size_t * count)
{
	mapkit_error err;

	err = map_h_vpd_getall(spm, array, count);
	if (err)
		MAPKIT_ERROR(err);

	qsort(*array, *count, sizeof(**array), map_h_vpd_compare);

	return MAPKIT_OK;
}

mapkit_error map_h_vpd_setall(map_h_vpd * spm, map_h_vpd_element * array, mapkit_size_t count)
{
	mapkit_size_t array_iindex;
	mapkit_error err;

	err = map_h_vpd_ensurecapacity(spm, spm->used + count);
	if (err)
		MAPKIT_ERROR(err);

	if (spm->alwaysdefault)
		/*
		 * Prevent shrinking 
		 */
		spm->minused = 0;

	for (array_iindex = 0; array_iindex < count; array_iindex++) {
		err = map_h_vpd_set(spm, array[array_iindex].key, array[array_iindex].value);
		if (err)
			MAPKIT_ERROR(err);
	}

	if (spm->alwaysdefault)
		map_h_vpd_adjustcapacity(spm);

	return MAPKIT_OK;
}

mapkit_error map_h_vpd_removeall(map_h_vpd * spm, void **array, mapkit_size_t count)
{
	mapkit_size_t array_iindex;
	mapkit_error err;

	/*
	 * Prevent shrinking 
	 */
	spm->minused = 0;

	for (array_iindex = 0; array_iindex < count; array_iindex++) {
		err = map_h_vpd_remove(spm, array[array_iindex]);
		if (err)
			MAPKIT_ERROR(err);
	}

	map_h_vpd_adjustcapacity(spm);

	return MAPKIT_OK;
}

void map_h_vpd_printstats(map_h_vpd * spm)
{
	fprintf(stderr, "MAPKIT: map_h_vpd statistics\n");
	fprintf(stderr, "MAPKIT: alwaysdefault = %d\n", spm->alwaysdefault);
	fprintf(stderr, "MAPKIT: minused = %ld, maxfill = %ld\n", (long) spm->minused, (long) spm->maxfill);
	fprintf(stderr, "MAPKIT: minusedfactor = %g, maxfillfactor = %g\n", spm->minusedfactor, spm->maxfillfactor);
#ifdef MAPKIT_COLLISIONS
	fprintf(stderr, "MAPKIT: insertionindexs = %lu, collisions = %lu\n", (unsigned long) spm->insertionindexs,
		(unsigned long) spm->insertionindex_collisions);
	fprintf(stderr, "MAPKIT: keyindexs = %lu, collisions = %lu\n", (unsigned long) spm->keyindexs, (unsigned long) spm->keyindex_collisions);
#endif
}

#endif							       /* MAPKIT_map_h_vpd */

#ifdef MAPKIT_map_h_vpvp

/*
  Implementation for map_h_vpvp (void* -> void*)
  Default value : NULL
  Uses a state field.
*/

/* Static prototypes */

/* Return the iindex of key, or MAPKIT_KEYNOTFOUND */
static mapkit_size_t map_h_vpvp_keyindex(map_h_vpvp * spm, void *key, mapkit_hash_t hash);

/* Return the iindex of key or -(insertion iindex)-1 if key not found */
static mapkit_size_t map_h_vpvp_insertionindex(map_h_vpvp * spm, void *key, mapkit_hash_t hash);

/* Implementation */

mapkit_error map_h_vpvp_init(map_h_vpvp * spm)
{
	return map_h_vpvp_init_hint(spm, MAPKIT_DEFAULT_EXPECTEDUSED);
}

mapkit_error map_h_vpvp_init_hint(map_h_vpvp * spm, mapkit_size_t used)
{
#ifdef MAPKIT_DEBUG
	fprintf(stderr, "MAPKIT: init\n");
#endif

	spm->size = 0;
	spm->fill = 0;
	spm->used = 0;
	spm->maxfillfactor = 0.5;
	spm->minusedfactor = 0.2;
	spm->contents = NULL;
#ifdef MAPKIT_COLLISIONS
	spm->insertionindexs = spm->insertionindex_collisions = 0;
	spm->keyindexs = spm->keyindex_collisions = 0;
#endif
	spm->defaultvalue = NULL;
	spm->alwaysdefault = 0;

	return map_h_vpvp_reallocate(spm, map_h_vpvp_meansize(spm, used));
}

mapkit_error map_h_vpvp_ensurecapacity(map_h_vpvp * spm, mapkit_size_t used)
{
	if (used > (spm->used + spm->maxfill - spm->fill)) {
#ifdef MAPKIT_DEBUG
		fprintf(stderr, "MAPKIT: need more capacity\n");
#endif
		return map_h_vpvp_reallocate(spm, map_h_vpvp_meansize(spm, used));
	} else
		return MAPKIT_OK;
}

mapkit_error map_h_vpvp_adjustcapacity(map_h_vpvp * spm)
{
	spm->minused = (mapkit_size_t) (spm->size * spm->minusedfactor);
	spm->maxfill = (mapkit_size_t) (spm->size * spm->maxfillfactor);

	if (spm->used < spm->minused) {
#ifdef MAPKIT_DEBUG
		fprintf(stderr, "MAPKIT: used < minused\n");
#endif
		return map_h_vpvp_reallocate(spm, map_h_vpvp_meansize(spm, spm->used));
	} else if (spm->fill > spm->maxfill) {
#ifdef MAPKIT_DEBUG
		fprintf(stderr, "MAPKIT: fill > maxfill\n");
#endif
		return map_h_vpvp_reallocate(spm, map_h_vpvp_meansize(spm, spm->used));
	} else
		return MAPKIT_OK;
}

void map_h_vpvp_free(map_h_vpvp * spm)
{
	free(spm->contents);
	spm->contents = NULL;
	spm->size = 0;
	spm->fill = 0;
	spm->used = 0;
	spm->maxfill = 0;
	spm->minused = 0;
#ifdef MAPKIT_COLLISIONS
	spm->insertionindexs = spm->insertionindex_collisions = 0;
	spm->keyindexs = spm->keyindex_collisions = 0;
#endif
}

mapkit_error map_h_vpvp_copy(map_h_vpvp * to, map_h_vpvp * from)
{
	map_h_vpvp_storage *contentscopy;
	size_t size = from->size * sizeof(*from->contents);

	contentscopy = (map_h_vpvp_storage *) malloc(size);
	if (contentscopy == NULL)
		MAPKIT_ERROR(MAPKIT_ENOMEM);

	memcpy(to, from, sizeof(map_h_vpvp));
	to->contents = contentscopy;
	memcpy(to->contents, from->contents, size);

	return MAPKIT_OK;
}

mapkit_size_t map_h_vpvp_growsize(map_h_vpvp * spm, mapkit_size_t used)
{
	return (mapkit_size_t) (4.0 * used / (3 * spm->minusedfactor + spm->maxfillfactor));
}

mapkit_size_t map_h_vpvp_shrinksize(map_h_vpvp * spm, mapkit_size_t used)
{
	return (mapkit_size_t) (4.0 * used / (spm->minusedfactor + 3 * spm->maxfillfactor));
}

mapkit_size_t map_h_vpvp_meansize(map_h_vpvp * spm, mapkit_size_t used)
{
	return (mapkit_size_t) (2.0 * used / (spm->minusedfactor + spm->maxfillfactor));
}

mapkit_error map_h_vpvp_reallocate(map_h_vpvp * spm, mapkit_size_t newsize)
{
	mapkit_size_t iindex;
	mapkit_size_t oldsize;
	map_h_vpvp_storage *newcontents, *oldcontents;

	/*
	 * At least one free entry 
	 */
	if (newsize <= spm->used)
		newsize = spm->used + 1;
	newsize = mapkit_nextprime(newsize);
	if (newsize <= spm->used)
		MAPKIT_ERROR(MAPKIT_ETOOBIG);

	newcontents = (map_h_vpvp_storage *) malloc(newsize * sizeof(*spm->contents));
	if (newcontents == NULL)
		MAPKIT_ERROR(MAPKIT_ENOMEM);

	/*
	 * Initialize all entries to "free" 
	 */
	for (iindex = 0; iindex < newsize; iindex++)
		newcontents[iindex].state = MAPKIT_FREESLOT;

	oldcontents = spm->contents;
	oldsize = spm->size;
	spm->contents = newcontents;
	spm->size = newsize;

#ifdef MAPKIT_DEBUG
	fprintf(stderr, "MAPKIT: reallocate %ld -> %ld\n", (long) oldsize, (long) newsize);
#endif

	spm->maxfill = (mapkit_size_t) (newsize * spm->maxfillfactor);
	/*
	 * At least one free entry 
	 */
	if (spm->maxfill >= newsize)
		spm->maxfill = newsize - 1;
	spm->minused = (mapkit_size_t) (newsize * spm->minusedfactor);
	spm->used = 0;

	if (oldcontents != NULL) {
		int used = 0;
		void *key;
		void *defaultvalue = spm->defaultvalue;
		int notalwaysdefault = !spm->alwaysdefault;

		/*
		 * Copy all entries from old to new 
		 */
		for (iindex = 0; iindex < oldsize; iindex++)
			if (oldcontents[iindex].state == MAPKIT_FULLSLOT) {
				mapkit_size_t ins_iindex;
				mapkit_hash_t hash;
				map_h_vpvp_storage *contents;

				key = oldcontents[iindex].key;

				/*
				 * Fast path 
				 */
				ins_iindex = (hash = mapkit_hash((mapkit_hash_t) key)) % spm->size;
				contents = &(newcontents[ins_iindex]);

				if (contents->state != MAPKIT_FREESLOT) {
					ins_iindex = map_h_vpvp_insertionindex(spm, key, hash);
					contents = &(newcontents[ins_iindex]);
				}
#ifdef MAPKIT_COLLISIONS
				else
					spm->insertionindexs++;
#endif
				if (notalwaysdefault || (!((oldcontents[iindex].value) == (defaultvalue)))) {
					contents->value = oldcontents[iindex].value;
					contents->state = MAPKIT_FULLSLOT;
					contents->key = key;
					used++;
				}
			}
		free(oldcontents);
		spm->used = used;
	}
	spm->fill = spm->used;

	return MAPKIT_OK;
}

void *map_h_vpvp_value_s(map_h_vpvp * spm, void *key, mapkit_hash_t hash)
{
	mapkit_size_t iindex;

	iindex = map_h_vpvp_keyindex(spm, key, hash);

	if (iindex < 0) {
		if (spm->alwaysdefault)
			return spm->defaultvalue;
		else
			MAPKIT_FATAL_ERROR(MAPKIT_EKEYNOTFOUND);
	}

	return spm->contents[iindex].value;
}

mapkit_error map_h_vpvp_get_s(map_h_vpvp * spm, void *key, void **value, mapkit_hash_t hash)
{
	mapkit_size_t iindex;

	iindex = map_h_vpvp_keyindex(spm, key, hash);

	if (iindex < 0) {
		if (spm->alwaysdefault) {
			*value = spm->defaultvalue;
			return MAPKIT_OK;
		} else
			MAPKIT_ERROR(MAPKIT_EKEYNOTFOUND);
	}

	*value = spm->contents[iindex].value;
	return MAPKIT_OK;
}

mapkit_error map_h_vpvp_set_s(map_h_vpvp * spm, void *key, void *value, mapkit_hash_t hash)
{
	mapkit_size_t iindex;

	iindex = map_h_vpvp_insertionindex(spm, key, hash);

	if (iindex < 0)
		/*
		 * FULLSLOT 
		 */
		spm->contents[-iindex - 1].value = value;
	else {
		map_h_vpvp_storage *element = &(spm->contents[iindex]);
		int ffree = element->state == MAPKIT_FREESLOT;

		element->state = MAPKIT_FULLSLOT;
		element->value = value;
		element->key = key;
		spm->used++;

		if (ffree && ((++spm->fill) > spm->maxfill)) {
#ifdef MAPKIT_DEBUG
			fprintf(stderr, "MAPKIT: fill > maxfill\n");
#endif
			return map_h_vpvp_reallocate(spm, map_h_vpvp_growsize(spm, spm->used));
		}
	}
	return MAPKIT_OK;
}

void **map_h_vpvp_insertptr_s(map_h_vpvp * spm, void *key, mapkit_hash_t hash)
{
	mapkit_size_t iindex;

	iindex = map_h_vpvp_insertionindex(spm, key, hash);

	if (iindex < 0)
		return &(spm->contents[-iindex - 1].value);
	else {
		map_h_vpvp_storage *element = &(spm->contents[iindex]);

		if (element->state == MAPKIT_FREESLOT) {
			/*
			 * FREESLOT 
			 */
			if (spm->fill >= spm->maxfill) {
				mapkit_error err;

#ifdef MAPKIT_DEBUG
				fprintf(stderr, "MAPKIT: fill => maxfill before insert\n");
#endif
				/*
				 * Must reallocate -before- inserting defaultvalue 
				 */
				err = map_h_vpvp_reallocate(spm, map_h_vpvp_growsize(spm, spm->used + 1));
				if (err) {
					MAPKIT_ERROR_NORET(err);
					return NULL;
				}

				iindex = map_h_vpvp_insertionindex(spm, key, hash);
				/*
				 * FREESLOT 
				 */
				spm->contents[iindex].state = MAPKIT_FULLSLOT;
				spm->contents[iindex].key = key;
				spm->contents[iindex].value = spm->defaultvalue;
				spm->used++;
				spm->fill++;
				return &(spm->contents[iindex].value);
			} else
				spm->fill++;
		}

		element->state = MAPKIT_FULLSLOT;
		element->key = key;
		element->value = spm->defaultvalue;
		spm->used++;
		return &(element->value);
	}
}

void **map_h_vpvp_ptr_s(map_h_vpvp * spm, void *key, mapkit_hash_t hash)
{
	mapkit_size_t iindex;

	iindex = map_h_vpvp_keyindex(spm, key, hash);

	if ((iindex >= 0) && ((!spm->alwaysdefault) || (!((spm->contents[iindex].value) == (spm->defaultvalue)))))
		return &(spm->contents[iindex].value);
	else
		return NULL;
}

mapkit_error map_h_vpvp_remove_s(map_h_vpvp * spm, void *key, mapkit_hash_t hash)
{
	mapkit_size_t iindex;

	iindex = map_h_vpvp_keyindex(spm, key, hash);

	if (iindex < 0) {
		if (spm->alwaysdefault)
			return MAPKIT_OK;
		else
			return MAPKIT_EKEYNOTFOUND;
	}

	spm->contents[iindex].state = MAPKIT_DELETEDSLOT;
	spm->used--;
	if (spm->used < spm->minused) {
#ifdef MAPKIT_DEBUG
		fprintf(stderr, "MAPKIT: used < minused\n");
#endif
		return map_h_vpvp_reallocate(spm, map_h_vpvp_shrinksize(spm, spm->used));
	}

	return MAPKIT_OK;
}

mapkit_size_t map_h_vpvp_keyindex(map_h_vpvp * spm, void *key, mapkit_hash_t hash)
{
	mapkit_size_t iindex, decrement;

	signed char state;

	iindex = hash % spm->size;
	decrement = (hash % (spm->size - 2));
	decrement += (decrement == 0);

#ifdef MAPKIT_COLLISIONS
	spm->keyindexs++;
#endif

	while ((state = spm->contents[iindex].state) != MAPKIT_FREESLOT
	       && (state == MAPKIT_DELETEDSLOT || (!((spm->contents[iindex].key) == (key))))) {
#ifdef MAPKIT_COLLISIONS
		spm->keyindex_collisions++;
#endif
		iindex -= decrement;
		if (iindex < 0)
			iindex += spm->size;
	}

	if (state == MAPKIT_FREESLOT)
		return MAPKIT_KEYNOTFOUND;
	return iindex;
}

mapkit_size_t map_h_vpvp_insertionindex(map_h_vpvp * spm, void *key, mapkit_hash_t hash)
{
	mapkit_size_t iindex, decrement;
	signed char state;

#ifdef MAPKIT_COLLISIONS
	spm->insertionindexs++;
#endif

	iindex = hash % spm->size;

	/*
	 * Fast path (largely superfluous) 
	 */
	if ((state = spm->contents[iindex].state) == MAPKIT_FREESLOT)
		return iindex;
	if ((state == MAPKIT_FULLSLOT)
	    && ((spm->contents[iindex].key) == (key)))
		return -iindex - 1;

	decrement = (hash % (spm->size - 2));
	decrement += (decrement == 0);

	while ((state == MAPKIT_FULLSLOT)
	       && (!((spm->contents[iindex].key) == (key)))) {
#ifdef MAPKIT_COLLISIONS
		spm->insertionindex_collisions++;
#endif
		iindex -= decrement;
		if (iindex < 0)
			iindex += spm->size;
		state = spm->contents[iindex].state;
	}

	if (state == MAPKIT_DELETEDSLOT) {
		mapkit_size_t index2 = iindex;

		while ((state = spm->contents[iindex].state) != MAPKIT_FREESLOT && ((state == MAPKIT_DELETEDSLOT)
										    || (!((spm->contents[iindex].key) == (key))))) {
			iindex -= decrement;
			if (iindex < 0)
				iindex += spm->size;
			state = spm->contents[iindex].state;
		}
		if (state == MAPKIT_FREESLOT)
			return index2;
	}

	if (state == MAPKIT_FULLSLOT)
		return -iindex - 1;
	return iindex;
}

mapkit_size_t map_h_vpvp_next(map_h_vpvp * spm, mapkit_size_t iindex)
{
	mapkit_size_t size = spm->size;
	map_h_vpvp_storage *pos_contents;
	void *defaultvalue = spm->defaultvalue;
	int notalwaysdefault = !spm->alwaysdefault;

	pos_contents = &(spm->contents[++iindex]);

	for (; iindex < size; iindex++, pos_contents++)
		if ((pos_contents->state == MAPKIT_FULLSLOT)
		    && (notalwaysdefault || (!((pos_contents->value) == (defaultvalue)))))
			return iindex;

	return -1;
}

map_h_vpvp_storage *map_h_vpvp_nextptr(map_h_vpvp * spm, map_h_vpvp_storage * pos_contents)
{
	map_h_vpvp_storage *end = &(spm->contents[spm->size]);
	void *defaultvalue = spm->defaultvalue;
	int notalwaysdefault = !spm->alwaysdefault;

	if (pos_contents == NULL)
		pos_contents = spm->contents;
	else {
		pos_contents++;
		if (pos_contents <= spm->contents)
			return NULL;
	}

	for (; pos_contents < end; pos_contents++)
		if ((pos_contents->state == MAPKIT_FULLSLOT)
		    && (notalwaysdefault || (!((pos_contents->value) == (defaultvalue)))))
			return pos_contents;

	return NULL;
}

mapkit_error map_h_vpvp_getall(map_h_vpvp * spm, map_h_vpvp_element ** array, mapkit_size_t * count)
{
	mapkit_size_t iindex;
	mapkit_size_t size = spm->size, vcount = 0;
	map_h_vpvp_element *pos_array;
	map_h_vpvp_storage *pos_contents;
	void *defaultvalue = spm->defaultvalue;
	int notalwaysdefault = !spm->alwaysdefault;

	pos_array = *array = (map_h_vpvp_element *) malloc(sizeof(**array) * spm->used);
	if (*array == NULL)
		MAPKIT_ERROR(MAPKIT_ENOMEM);

	pos_contents = spm->contents;

	for (iindex = 0; iindex < size; iindex++, pos_contents++)
		if ((pos_contents->state == MAPKIT_FULLSLOT)
		    && (notalwaysdefault || (!((pos_contents->value) == (defaultvalue))))) {
			pos_array->key = pos_contents->key;
			pos_array->value = pos_contents->value;
			pos_array++;
			vcount++;
		}
	*count = vcount;

	return MAPKIT_OK;
}

mapkit_error map_h_vpvp_clean(map_h_vpvp * spm)
{
	mapkit_size_t iindex, count = 0;
	mapkit_size_t size = spm->size;
	map_h_vpvp_storage *pos_contents;
	void *defaultvalue = spm->defaultvalue;

	pos_contents = spm->contents;
	for (iindex = 0; iindex < size; iindex++, pos_contents++)
		if ((pos_contents->state == MAPKIT_FULLSLOT)
		    && ((pos_contents->value) == (defaultvalue))) {
			pos_contents->state = MAPKIT_DELETEDSLOT;
			count++;
		}

	spm->used -= count;
	if (spm->used < spm->minused) {
#ifdef MAPKIT_DEBUG
		fprintf(stderr, "MAPKIT: used < minused\n");
#endif
		return map_h_vpvp_reallocate(spm, map_h_vpvp_meansize(spm, spm->used));
	}

	return MAPKIT_OK;
}

int map_h_vpvp_compare(const void *e1, const void *e2)
{
	void *key1 = ((const map_h_vpvp_element *) e1)->key;
	void *key2 = ((const map_h_vpvp_element *) e2)->key;

	return ((key1) < (key2) ? -1 : ((key1) == (key2) ? 0 : 1));
}

mapkit_error map_h_vpvp_getall_sorted(map_h_vpvp * spm, map_h_vpvp_element ** array, mapkit_size_t * count)
{
	mapkit_error err;

	err = map_h_vpvp_getall(spm, array, count);
	if (err)
		MAPKIT_ERROR(err);

	qsort(*array, *count, sizeof(**array), map_h_vpvp_compare);

	return MAPKIT_OK;
}

mapkit_error map_h_vpvp_setall(map_h_vpvp * spm, map_h_vpvp_element * array, mapkit_size_t count)
{
	mapkit_size_t array_iindex;
	mapkit_error err;

	err = map_h_vpvp_ensurecapacity(spm, spm->used + count);
	if (err)
		MAPKIT_ERROR(err);

	if (spm->alwaysdefault)
		/*
		 * Prevent shrinking 
		 */
		spm->minused = 0;

	for (array_iindex = 0; array_iindex < count; array_iindex++) {
		err = map_h_vpvp_set(spm, array[array_iindex].key, array[array_iindex].value);
		if (err)
			MAPKIT_ERROR(err);
	}

	if (spm->alwaysdefault)
		map_h_vpvp_adjustcapacity(spm);

	return MAPKIT_OK;
}

mapkit_error map_h_vpvp_removeall(map_h_vpvp * spm, void **array, mapkit_size_t count)
{
	mapkit_size_t array_iindex;
	mapkit_error err;

	/*
	 * Prevent shrinking 
	 */
	spm->minused = 0;

	for (array_iindex = 0; array_iindex < count; array_iindex++) {
		err = map_h_vpvp_remove(spm, array[array_iindex]);
		if (err)
			MAPKIT_ERROR(err);
	}

	map_h_vpvp_adjustcapacity(spm);

	return MAPKIT_OK;
}

void map_h_vpvp_printstats(map_h_vpvp * spm)
{
	fprintf(stderr, "MAPKIT: map_h_vpvp statistics\n");
	fprintf(stderr, "MAPKIT: alwaysdefault = %d\n", spm->alwaysdefault);
	fprintf(stderr, "MAPKIT: minused = %ld, maxfill = %ld\n", (long) spm->minused, (long) spm->maxfill);
	fprintf(stderr, "MAPKIT: minusedfactor = %g, maxfillfactor = %g\n", spm->minusedfactor, spm->maxfillfactor);
#ifdef MAPKIT_COLLISIONS
	fprintf(stderr, "MAPKIT: insertionindexs = %lu, collisions = %lu\n", (unsigned long) spm->insertionindexs,
		(unsigned long) spm->insertionindex_collisions);
	fprintf(stderr, "MAPKIT: keyindexs = %lu, collisions = %lu\n", (unsigned long) spm->keyindexs, (unsigned long) spm->keyindex_collisions);
#endif
}

#endif							       /* MAPKIT_map_h_vpvp */

#ifdef MAPKIT_map_stri

/*
  Implementation for map_stri (char* -> int)
  Default value : 0
  Uses a state field.
*/

/* Static prototypes */

/* Return the iindex of key, or MAPKIT_KEYNOTFOUND */
static mapkit_size_t map_stri_keyindex(map_stri * spm, char *key, mapkit_hash_t hash);

/* Return the iindex of key or -(insertion iindex)-1 if key not found */
static mapkit_size_t map_stri_insertionindex(map_stri * spm, char *key, mapkit_hash_t hash);

/* Implementation */

mapkit_error map_stri_init(map_stri * spm)
{
	return map_stri_init_hint(spm, MAPKIT_DEFAULT_EXPECTEDUSED);
}

mapkit_error map_stri_init_hint(map_stri * spm, mapkit_size_t used)
{
#ifdef MAPKIT_DEBUG
	fprintf(stderr, "MAPKIT: init\n");
#endif

	spm->size = 0;
	spm->fill = 0;
	spm->used = 0;
	spm->maxfillfactor = 0.5;
	spm->minusedfactor = 0.2;
	spm->contents = NULL;
#ifdef MAPKIT_COLLISIONS
	spm->insertionindexs = spm->insertionindex_collisions = 0;
	spm->keyindexs = spm->keyindex_collisions = 0;
#endif
	spm->defaultvalue = 0;
	spm->alwaysdefault = 0;

	return map_stri_reallocate(spm, map_stri_meansize(spm, used));
}

mapkit_error map_stri_ensurecapacity(map_stri * spm, mapkit_size_t used)
{
	if (used > (spm->used + spm->maxfill - spm->fill)) {
#ifdef MAPKIT_DEBUG
		fprintf(stderr, "MAPKIT: need more capacity\n");
#endif
		return map_stri_reallocate(spm, map_stri_meansize(spm, used));
	} else
		return MAPKIT_OK;
}

mapkit_error map_stri_adjustcapacity(map_stri * spm)
{
	spm->minused = (mapkit_size_t) (spm->size * spm->minusedfactor);
	spm->maxfill = (mapkit_size_t) (spm->size * spm->maxfillfactor);

	if (spm->used < spm->minused) {
#ifdef MAPKIT_DEBUG
		fprintf(stderr, "MAPKIT: used < minused\n");
#endif
		return map_stri_reallocate(spm, map_stri_meansize(spm, spm->used));
	} else if (spm->fill > spm->maxfill) {
#ifdef MAPKIT_DEBUG
		fprintf(stderr, "MAPKIT: fill > maxfill\n");
#endif
		return map_stri_reallocate(spm, map_stri_meansize(spm, spm->used));
	} else
		return MAPKIT_OK;
}

void map_stri_free(map_stri * spm)
{
	free(spm->contents);
	spm->contents = NULL;
	spm->size = 0;
	spm->fill = 0;
	spm->used = 0;
	spm->maxfill = 0;
	spm->minused = 0;
#ifdef MAPKIT_COLLISIONS
	spm->insertionindexs = spm->insertionindex_collisions = 0;
	spm->keyindexs = spm->keyindex_collisions = 0;
#endif
}

mapkit_error map_stri_copy(map_stri * to, map_stri * from)
{
	map_stri_storage *contentscopy;
	size_t size = from->size * sizeof(*from->contents);

	contentscopy = (map_stri_storage *) malloc(size);
	if (contentscopy == NULL)
		MAPKIT_ERROR(MAPKIT_ENOMEM);

	memcpy(to, from, sizeof(map_stri));
	to->contents = contentscopy;
	memcpy(to->contents, from->contents, size);

	return MAPKIT_OK;
}

mapkit_size_t map_stri_growsize(map_stri * spm, mapkit_size_t used)
{
	return (mapkit_size_t) (4.0 * used / (3 * spm->minusedfactor + spm->maxfillfactor));
}

mapkit_size_t map_stri_shrinksize(map_stri * spm, mapkit_size_t used)
{
	return (mapkit_size_t) (4.0 * used / (spm->minusedfactor + 3 * spm->maxfillfactor));
}

mapkit_size_t map_stri_meansize(map_stri * spm, mapkit_size_t used)
{
	return (mapkit_size_t) (2.0 * used / (spm->minusedfactor + spm->maxfillfactor));
}

mapkit_error map_stri_reallocate(map_stri * spm, mapkit_size_t newsize)
{
	mapkit_size_t iindex;
	mapkit_size_t oldsize;
	map_stri_storage *newcontents, *oldcontents;

	/*
	 * At least one free entry 
	 */
	if (newsize <= spm->used)
		newsize = spm->used + 1;
	newsize = mapkit_nextprime(newsize);
	if (newsize <= spm->used)
		MAPKIT_ERROR(MAPKIT_ETOOBIG);

	newcontents = (map_stri_storage *) malloc(newsize * sizeof(*spm->contents));
	if (newcontents == NULL)
		MAPKIT_ERROR(MAPKIT_ENOMEM);

	/*
	 * Initialize all entries to "free" 
	 */
	for (iindex = 0; iindex < newsize; iindex++)
		newcontents[iindex].state = MAPKIT_FREESLOT;

	oldcontents = spm->contents;
	oldsize = spm->size;
	spm->contents = newcontents;
	spm->size = newsize;

#ifdef MAPKIT_DEBUG
	fprintf(stderr, "MAPKIT: reallocate %ld -> %ld\n", (long) oldsize, (long) newsize);
#endif

	spm->maxfill = (mapkit_size_t) (newsize * spm->maxfillfactor);
	/*
	 * At least one free entry 
	 */
	if (spm->maxfill >= newsize)
		spm->maxfill = newsize - 1;
	spm->minused = (mapkit_size_t) (newsize * spm->minusedfactor);
	spm->used = 0;

	if (oldcontents != NULL) {
		int used = 0;
		char *key;
		int defaultvalue = spm->defaultvalue;
		int notalwaysdefault = !spm->alwaysdefault;

		/*
		 * Copy all entries from old to new 
		 */
		for (iindex = 0; iindex < oldsize; iindex++)
			if (oldcontents[iindex].state == MAPKIT_FULLSLOT) {
				mapkit_size_t ins_iindex;
				mapkit_hash_t hash;
				map_stri_storage *contents;

				key = oldcontents[iindex].key;

				/*
				 * Fast path 
				 */
				ins_iindex = (hash = mapkit_strhash(key)) % spm->size;
				contents = &(newcontents[ins_iindex]);

				if (contents->state != MAPKIT_FREESLOT) {
					ins_iindex = map_stri_insertionindex(spm, key, hash);
					contents = &(newcontents[ins_iindex]);
				}
#ifdef MAPKIT_COLLISIONS
				else
					spm->insertionindexs++;
#endif
				if (notalwaysdefault || (!((oldcontents[iindex].value) == (defaultvalue)))) {
					contents->value = oldcontents[iindex].value;
					contents->state = MAPKIT_FULLSLOT;
					contents->key = key;
					used++;
				}
			}
		free(oldcontents);
		spm->used = used;
	}
	spm->fill = spm->used;

	return MAPKIT_OK;
}

int map_stri_value_s(map_stri * spm, char *key, mapkit_hash_t hash)
{
	mapkit_size_t iindex;

	iindex = map_stri_keyindex(spm, key, hash);

	if (iindex < 0) {
		if (spm->alwaysdefault)
			return spm->defaultvalue;
		else
			MAPKIT_FATAL_ERROR(MAPKIT_EKEYNOTFOUND);
	}

	return spm->contents[iindex].value;
}

mapkit_error map_stri_get_s(map_stri * spm, char *key, int *value, mapkit_hash_t hash)
{
	mapkit_size_t iindex;

	iindex = map_stri_keyindex(spm, key, hash);

	if (iindex < 0) {
		if (spm->alwaysdefault) {
			*value = spm->defaultvalue;
			return MAPKIT_OK;
		} else
			MAPKIT_ERROR(MAPKIT_EKEYNOTFOUND);
	}

	*value = spm->contents[iindex].value;
	return MAPKIT_OK;
}

mapkit_error map_stri_set_s(map_stri * spm, char *key, int value, mapkit_hash_t hash)
{
	mapkit_size_t iindex;

	iindex = map_stri_insertionindex(spm, key, hash);

	if (iindex < 0)
		/*
		 * FULLSLOT 
		 */
		spm->contents[-iindex - 1].value = value;
	else {
		map_stri_storage *element = &(spm->contents[iindex]);
		int ffree = element->state == MAPKIT_FREESLOT;

		element->state = MAPKIT_FULLSLOT;
		element->value = value;
		element->key = key;
		spm->used++;

		if (ffree && ((++spm->fill) > spm->maxfill)) {
#ifdef MAPKIT_DEBUG
			fprintf(stderr, "MAPKIT: fill > maxfill\n");
#endif
			return map_stri_reallocate(spm, map_stri_growsize(spm, spm->used));
		}
	}
	return MAPKIT_OK;
}

int *map_stri_insertptr_s(map_stri * spm, char *key, mapkit_hash_t hash)
{
	mapkit_size_t iindex;

	iindex = map_stri_insertionindex(spm, key, hash);

	if (iindex < 0)
		return &(spm->contents[-iindex - 1].value);
	else {
		map_stri_storage *element = &(spm->contents[iindex]);

		if (element->state == MAPKIT_FREESLOT) {
			/*
			 * FREESLOT 
			 */
			if (spm->fill >= spm->maxfill) {
				mapkit_error err;

#ifdef MAPKIT_DEBUG
				fprintf(stderr, "MAPKIT: fill => maxfill before insert\n");
#endif
				/*
				 * Must reallocate -before- inserting defaultvalue 
				 */
				err = map_stri_reallocate(spm, map_stri_growsize(spm, spm->used + 1));
				if (err) {
					MAPKIT_ERROR_NORET(err);
					return NULL;
				}

				iindex = map_stri_insertionindex(spm, key, hash);
				/*
				 * FREESLOT 
				 */
				spm->contents[iindex].state = MAPKIT_FULLSLOT;
				spm->contents[iindex].key = key;
				spm->contents[iindex].value = spm->defaultvalue;
				spm->used++;
				spm->fill++;
				return &(spm->contents[iindex].value);
			} else
				spm->fill++;
		}

		element->state = MAPKIT_FULLSLOT;
		element->key = key;
		element->value = spm->defaultvalue;
		spm->used++;
		return &(element->value);
	}
}

int *map_stri_ptr_s(map_stri * spm, char *key, mapkit_hash_t hash)
{
	mapkit_size_t iindex;

	iindex = map_stri_keyindex(spm, key, hash);

	if ((iindex >= 0) && ((!spm->alwaysdefault) || (!((spm->contents[iindex].value) == (spm->defaultvalue)))))
		return &(spm->contents[iindex].value);
	else
		return NULL;
}

mapkit_error map_stri_remove_s(map_stri * spm, char *key, mapkit_hash_t hash)
{
	mapkit_size_t iindex;

	iindex = map_stri_keyindex(spm, key, hash);

	if (iindex < 0) {
		if (spm->alwaysdefault)
			return MAPKIT_OK;
		else
			return MAPKIT_EKEYNOTFOUND;
	}

	spm->contents[iindex].state = MAPKIT_DELETEDSLOT;
	spm->used--;
	if (spm->used < spm->minused) {
#ifdef MAPKIT_DEBUG
		fprintf(stderr, "MAPKIT: used < minused\n");
#endif
		return map_stri_reallocate(spm, map_stri_shrinksize(spm, spm->used));
	}

	return MAPKIT_OK;
}

mapkit_size_t map_stri_keyindex(map_stri * spm, char *key, mapkit_hash_t hash)
{
	mapkit_size_t iindex, decrement;

	signed char state;

	iindex = hash % spm->size;
	decrement = (hash % (spm->size - 2));
	decrement += (decrement == 0);

#ifdef MAPKIT_COLLISIONS
	spm->keyindexs++;
#endif

	while ((state = spm->contents[iindex].state) != MAPKIT_FREESLOT
	       && (state == MAPKIT_DELETEDSLOT || (!(strcmp(spm->contents[iindex].key, key) == 0)))) {
#ifdef MAPKIT_COLLISIONS
		spm->keyindex_collisions++;
#endif
		iindex -= decrement;
		if (iindex < 0)
			iindex += spm->size;
	}

	if (state == MAPKIT_FREESLOT)
		return MAPKIT_KEYNOTFOUND;
	return iindex;
}

mapkit_size_t map_stri_insertionindex(map_stri * spm, char *key, mapkit_hash_t hash)
{
	mapkit_size_t iindex, decrement;
	signed char state;

#ifdef MAPKIT_COLLISIONS
	spm->insertionindexs++;
#endif

	iindex = hash % spm->size;

	/*
	 * Fast path (largely superfluous) 
	 */
	if ((state = spm->contents[iindex].state) == MAPKIT_FREESLOT)
		return iindex;
	if ((state == MAPKIT_FULLSLOT)
	    && (strcmp(spm->contents[iindex].key, key) == 0))
		return -iindex - 1;

	decrement = (hash % (spm->size - 2));
	decrement += (decrement == 0);

	while ((state == MAPKIT_FULLSLOT)
	       && (!(strcmp(spm->contents[iindex].key, key) == 0))) {
#ifdef MAPKIT_COLLISIONS
		spm->insertionindex_collisions++;
#endif
		iindex -= decrement;
		if (iindex < 0)
			iindex += spm->size;
		state = spm->contents[iindex].state;
	}

	if (state == MAPKIT_DELETEDSLOT) {
		mapkit_size_t index2 = iindex;

		while ((state = spm->contents[iindex].state) != MAPKIT_FREESLOT && ((state == MAPKIT_DELETEDSLOT)
										    || (!(strcmp(spm->contents[iindex].key, key) == 0)))) {
			iindex -= decrement;
			if (iindex < 0)
				iindex += spm->size;
			state = spm->contents[iindex].state;
		}
		if (state == MAPKIT_FREESLOT)
			return index2;
	}

	if (state == MAPKIT_FULLSLOT)
		return -iindex - 1;
	return iindex;
}

mapkit_size_t map_stri_next(map_stri * spm, mapkit_size_t iindex)
{
	mapkit_size_t size = spm->size;
	map_stri_storage *pos_contents;
	int defaultvalue = spm->defaultvalue;
	int notalwaysdefault = !spm->alwaysdefault;

	pos_contents = &(spm->contents[++iindex]);

	for (; iindex < size; iindex++, pos_contents++)
		if ((pos_contents->state == MAPKIT_FULLSLOT)
		    && (notalwaysdefault || (!((pos_contents->value) == (defaultvalue)))))
			return iindex;

	return -1;
}

map_stri_storage *map_stri_nextptr(map_stri * spm, map_stri_storage * pos_contents)
{
	map_stri_storage *end = &(spm->contents[spm->size]);
	int defaultvalue = spm->defaultvalue;
	int notalwaysdefault = !spm->alwaysdefault;

	if (pos_contents == NULL)
		pos_contents = spm->contents;
	else {
		pos_contents++;
		if (pos_contents <= spm->contents)
			return NULL;
	}

	for (; pos_contents < end; pos_contents++)
		if ((pos_contents->state == MAPKIT_FULLSLOT)
		    && (notalwaysdefault || (!((pos_contents->value) == (defaultvalue)))))
			return pos_contents;

	return NULL;
}

mapkit_error map_stri_getall(map_stri * spm, map_stri_element ** array, mapkit_size_t * count)
{
	mapkit_size_t iindex;
	mapkit_size_t size = spm->size, vcount = 0;
	map_stri_element *pos_array;
	map_stri_storage *pos_contents;
	int defaultvalue = spm->defaultvalue;
	int notalwaysdefault = !spm->alwaysdefault;

	pos_array = *array = (map_stri_element *) malloc(sizeof(**array) * spm->used);
	if (*array == NULL)
		MAPKIT_ERROR(MAPKIT_ENOMEM);

	pos_contents = spm->contents;

	for (iindex = 0; iindex < size; iindex++, pos_contents++)
		if ((pos_contents->state == MAPKIT_FULLSLOT)
		    && (notalwaysdefault || (!((pos_contents->value) == (defaultvalue))))) {
			pos_array->key = pos_contents->key;
			pos_array->value = pos_contents->value;
			pos_array++;
			vcount++;
		}
	*count = vcount;

	return MAPKIT_OK;
}

mapkit_error map_stri_clean(map_stri * spm)
{
	mapkit_size_t iindex, count = 0;
	mapkit_size_t size = spm->size;
	map_stri_storage *pos_contents;
	int defaultvalue = spm->defaultvalue;

	pos_contents = spm->contents;
	for (iindex = 0; iindex < size; iindex++, pos_contents++)
		if ((pos_contents->state == MAPKIT_FULLSLOT)
		    && ((pos_contents->value) == (defaultvalue))) {
			pos_contents->state = MAPKIT_DELETEDSLOT;
			count++;
		}

	spm->used -= count;
	if (spm->used < spm->minused) {
#ifdef MAPKIT_DEBUG
		fprintf(stderr, "MAPKIT: used < minused\n");
#endif
		return map_stri_reallocate(spm, map_stri_meansize(spm, spm->used));
	}

	return MAPKIT_OK;
}

int map_stri_compare(const void *e1, const void *e2)
{
	char *key1 = ((const map_stri_element *) e1)->key;
	char *key2 = ((const map_stri_element *) e2)->key;

	return strcmp(key1, key2);
}

mapkit_error map_stri_getall_sorted(map_stri * spm, map_stri_element ** array, mapkit_size_t * count)
{
	mapkit_error err;

	err = map_stri_getall(spm, array, count);
	if (err)
		MAPKIT_ERROR(err);

	qsort(*array, *count, sizeof(**array), map_stri_compare);

	return MAPKIT_OK;
}

mapkit_error map_stri_setall(map_stri * spm, map_stri_element * array, mapkit_size_t count)
{
	mapkit_size_t array_iindex;
	mapkit_error err;

	err = map_stri_ensurecapacity(spm, spm->used + count);
	if (err)
		MAPKIT_ERROR(err);

	if (spm->alwaysdefault)
		/*
		 * Prevent shrinking 
		 */
		spm->minused = 0;

	for (array_iindex = 0; array_iindex < count; array_iindex++) {
		err = map_stri_set(spm, array[array_iindex].key, array[array_iindex].value);
		if (err)
			MAPKIT_ERROR(err);
	}

	if (spm->alwaysdefault)
		map_stri_adjustcapacity(spm);

	return MAPKIT_OK;
}

mapkit_error map_stri_removeall(map_stri * spm, char **array, mapkit_size_t count)
{
	mapkit_size_t array_iindex;
	mapkit_error err;

	/*
	 * Prevent shrinking 
	 */
	spm->minused = 0;

	for (array_iindex = 0; array_iindex < count; array_iindex++) {
		err = map_stri_remove(spm, array[array_iindex]);
		if (err)
			MAPKIT_ERROR(err);
	}

	map_stri_adjustcapacity(spm);

	return MAPKIT_OK;
}

void map_stri_printstats(map_stri * spm)
{
	fprintf(stderr, "MAPKIT: map_stri statistics\n");
	fprintf(stderr, "MAPKIT: alwaysdefault = %d\n", spm->alwaysdefault);
	fprintf(stderr, "MAPKIT: minused = %ld, maxfill = %ld\n", (long) spm->minused, (long) spm->maxfill);
	fprintf(stderr, "MAPKIT: minusedfactor = %g, maxfillfactor = %g\n", spm->minusedfactor, spm->maxfillfactor);
#ifdef MAPKIT_COLLISIONS
	fprintf(stderr, "MAPKIT: insertionindexs = %lu, collisions = %lu\n", (unsigned long) spm->insertionindexs,
		(unsigned long) spm->insertionindex_collisions);
	fprintf(stderr, "MAPKIT: keyindexs = %lu, collisions = %lu\n", (unsigned long) spm->keyindexs, (unsigned long) spm->keyindex_collisions);
#endif
}

#endif							       /* MAPKIT_map_stri */

#ifdef MAPKIT_map_strd

/*
  Implementation for map_strd (char* -> double)
  Default value : 0.0
  Uses a state field.
*/

/* Static prototypes */

/* Return the iindex of key, or MAPKIT_KEYNOTFOUND */
static mapkit_size_t map_strd_keyindex(map_strd * spm, char *key, mapkit_hash_t hash);

/* Return the iindex of key or -(insertion iindex)-1 if key not found */
static mapkit_size_t map_strd_insertionindex(map_strd * spm, char *key, mapkit_hash_t hash);

/* Implementation */

mapkit_error map_strd_init(map_strd * spm)
{
	return map_strd_init_hint(spm, MAPKIT_DEFAULT_EXPECTEDUSED);
}

mapkit_error map_strd_init_hint(map_strd * spm, mapkit_size_t used)
{
#ifdef MAPKIT_DEBUG
	fprintf(stderr, "MAPKIT: init\n");
#endif

	spm->size = 0;
	spm->fill = 0;
	spm->used = 0;
	spm->maxfillfactor = 0.5;
	spm->minusedfactor = 0.2;
	spm->contents = NULL;
#ifdef MAPKIT_COLLISIONS
	spm->insertionindexs = spm->insertionindex_collisions = 0;
	spm->keyindexs = spm->keyindex_collisions = 0;
#endif
	spm->defaultvalue = 0.0;
	spm->alwaysdefault = 0;

	return map_strd_reallocate(spm, map_strd_meansize(spm, used));
}

mapkit_error map_strd_ensurecapacity(map_strd * spm, mapkit_size_t used)
{
	if (used > (spm->used + spm->maxfill - spm->fill)) {
#ifdef MAPKIT_DEBUG
		fprintf(stderr, "MAPKIT: need more capacity\n");
#endif
		return map_strd_reallocate(spm, map_strd_meansize(spm, used));
	} else
		return MAPKIT_OK;
}

mapkit_error map_strd_adjustcapacity(map_strd * spm)
{
	spm->minused = (mapkit_size_t) (spm->size * spm->minusedfactor);
	spm->maxfill = (mapkit_size_t) (spm->size * spm->maxfillfactor);

	if (spm->used < spm->minused) {
#ifdef MAPKIT_DEBUG
		fprintf(stderr, "MAPKIT: used < minused\n");
#endif
		return map_strd_reallocate(spm, map_strd_meansize(spm, spm->used));
	} else if (spm->fill > spm->maxfill) {
#ifdef MAPKIT_DEBUG
		fprintf(stderr, "MAPKIT: fill > maxfill\n");
#endif
		return map_strd_reallocate(spm, map_strd_meansize(spm, spm->used));
	} else
		return MAPKIT_OK;
}

void map_strd_free(map_strd * spm)
{
	free(spm->contents);
	spm->contents = NULL;
	spm->size = 0;
	spm->fill = 0;
	spm->used = 0;
	spm->maxfill = 0;
	spm->minused = 0;
#ifdef MAPKIT_COLLISIONS
	spm->insertionindexs = spm->insertionindex_collisions = 0;
	spm->keyindexs = spm->keyindex_collisions = 0;
#endif
}

mapkit_error map_strd_copy(map_strd * to, map_strd * from)
{
	map_strd_storage *contentscopy;
	size_t size = from->size * sizeof(*from->contents);

	contentscopy = (map_strd_storage *) malloc(size);
	if (contentscopy == NULL)
		MAPKIT_ERROR(MAPKIT_ENOMEM);

	memcpy(to, from, sizeof(map_strd));
	to->contents = contentscopy;
	memcpy(to->contents, from->contents, size);

	return MAPKIT_OK;
}

mapkit_size_t map_strd_growsize(map_strd * spm, mapkit_size_t used)
{
	return (mapkit_size_t) (4.0 * used / (3 * spm->minusedfactor + spm->maxfillfactor));
}

mapkit_size_t map_strd_shrinksize(map_strd * spm, mapkit_size_t used)
{
	return (mapkit_size_t) (4.0 * used / (spm->minusedfactor + 3 * spm->maxfillfactor));
}

mapkit_size_t map_strd_meansize(map_strd * spm, mapkit_size_t used)
{
	return (mapkit_size_t) (2.0 * used / (spm->minusedfactor + spm->maxfillfactor));
}

mapkit_error map_strd_reallocate(map_strd * spm, mapkit_size_t newsize)
{
	mapkit_size_t iindex;
	mapkit_size_t oldsize;
	map_strd_storage *newcontents, *oldcontents;

	/*
	 * At least one free entry 
	 */
	if (newsize <= spm->used)
		newsize = spm->used + 1;
	newsize = mapkit_nextprime(newsize);
	if (newsize <= spm->used)
		MAPKIT_ERROR(MAPKIT_ETOOBIG);

	newcontents = (map_strd_storage *) malloc(newsize * sizeof(*spm->contents));
	if (newcontents == NULL)
		MAPKIT_ERROR(MAPKIT_ENOMEM);

	/*
	 * Initialize all entries to "free" 
	 */
	for (iindex = 0; iindex < newsize; iindex++)
		newcontents[iindex].state = MAPKIT_FREESLOT;

	oldcontents = spm->contents;
	oldsize = spm->size;
	spm->contents = newcontents;
	spm->size = newsize;

#ifdef MAPKIT_DEBUG
	fprintf(stderr, "MAPKIT: reallocate %ld -> %ld\n", (long) oldsize, (long) newsize);
#endif

	spm->maxfill = (mapkit_size_t) (newsize * spm->maxfillfactor);
	/*
	 * At least one free entry 
	 */
	if (spm->maxfill >= newsize)
		spm->maxfill = newsize - 1;
	spm->minused = (mapkit_size_t) (newsize * spm->minusedfactor);
	spm->used = 0;

	if (oldcontents != NULL) {
		int used = 0;
		char *key;
		double defaultvalue = spm->defaultvalue;
		int notalwaysdefault = !spm->alwaysdefault;

		/*
		 * Copy all entries from old to new 
		 */
		for (iindex = 0; iindex < oldsize; iindex++)
			if (oldcontents[iindex].state == MAPKIT_FULLSLOT) {
				mapkit_size_t ins_iindex;
				mapkit_hash_t hash;
				map_strd_storage *contents;

				key = oldcontents[iindex].key;

				/*
				 * Fast path 
				 */
				ins_iindex = (hash = mapkit_strhash(key)) % spm->size;
				contents = &(newcontents[ins_iindex]);

				if (contents->state != MAPKIT_FREESLOT) {
					ins_iindex = map_strd_insertionindex(spm, key, hash);
					contents = &(newcontents[ins_iindex]);
				}
#ifdef MAPKIT_COLLISIONS
				else
					spm->insertionindexs++;
#endif
				if (notalwaysdefault || (!((oldcontents[iindex].value) == (defaultvalue)))) {
					contents->value = oldcontents[iindex].value;
					contents->state = MAPKIT_FULLSLOT;
					contents->key = key;
					used++;
				}
			}
		free(oldcontents);
		spm->used = used;
	}
	spm->fill = spm->used;

	return MAPKIT_OK;
}

double map_strd_value_s(map_strd * spm, char *key, mapkit_hash_t hash)
{
	mapkit_size_t iindex;

	iindex = map_strd_keyindex(spm, key, hash);

	if (iindex < 0) {
		if (spm->alwaysdefault)
			return spm->defaultvalue;
		else
			MAPKIT_FATAL_ERROR(MAPKIT_EKEYNOTFOUND);
	}

	return spm->contents[iindex].value;
}

mapkit_error map_strd_get_s(map_strd * spm, char *key, double *value, mapkit_hash_t hash)
{
	mapkit_size_t iindex;

	iindex = map_strd_keyindex(spm, key, hash);

	if (iindex < 0) {
		if (spm->alwaysdefault) {
			*value = spm->defaultvalue;
			return MAPKIT_OK;
		} else
			MAPKIT_ERROR(MAPKIT_EKEYNOTFOUND);
	}

	*value = spm->contents[iindex].value;
	return MAPKIT_OK;
}

mapkit_error map_strd_set_s(map_strd * spm, char *key, double value, mapkit_hash_t hash)
{
	mapkit_size_t iindex;

	iindex = map_strd_insertionindex(spm, key, hash);

	if (iindex < 0)
		/*
		 * FULLSLOT 
		 */
		spm->contents[-iindex - 1].value = value;
	else {
		map_strd_storage *element = &(spm->contents[iindex]);
		int ffree = element->state == MAPKIT_FREESLOT;

		element->state = MAPKIT_FULLSLOT;
		element->value = value;
		element->key = key;
		spm->used++;

		if (ffree && ((++spm->fill) > spm->maxfill)) {
#ifdef MAPKIT_DEBUG
			fprintf(stderr, "MAPKIT: fill > maxfill\n");
#endif
			return map_strd_reallocate(spm, map_strd_growsize(spm, spm->used));
		}
	}
	return MAPKIT_OK;
}

double *map_strd_insertptr_s(map_strd * spm, char *key, mapkit_hash_t hash)
{
	mapkit_size_t iindex;

	iindex = map_strd_insertionindex(spm, key, hash);

	if (iindex < 0)
		return &(spm->contents[-iindex - 1].value);
	else {
		map_strd_storage *element = &(spm->contents[iindex]);

		if (element->state == MAPKIT_FREESLOT) {
			/*
			 * FREESLOT 
			 */
			if (spm->fill >= spm->maxfill) {
				mapkit_error err;

#ifdef MAPKIT_DEBUG
				fprintf(stderr, "MAPKIT: fill => maxfill before insert\n");
#endif
				/*
				 * Must reallocate -before- inserting defaultvalue 
				 */
				err = map_strd_reallocate(spm, map_strd_growsize(spm, spm->used + 1));
				if (err) {
					MAPKIT_ERROR_NORET(err);
					return NULL;
				}

				iindex = map_strd_insertionindex(spm, key, hash);
				/*
				 * FREESLOT 
				 */
				spm->contents[iindex].state = MAPKIT_FULLSLOT;
				spm->contents[iindex].key = key;
				spm->contents[iindex].value = spm->defaultvalue;
				spm->used++;
				spm->fill++;
				return &(spm->contents[iindex].value);
			} else
				spm->fill++;
		}

		element->state = MAPKIT_FULLSLOT;
		element->key = key;
		element->value = spm->defaultvalue;
		spm->used++;
		return &(element->value);
	}
}

double *map_strd_ptr_s(map_strd * spm, char *key, mapkit_hash_t hash)
{
	mapkit_size_t iindex;

	iindex = map_strd_keyindex(spm, key, hash);

	if ((iindex >= 0) && ((!spm->alwaysdefault) || (!((spm->contents[iindex].value) == (spm->defaultvalue)))))
		return &(spm->contents[iindex].value);
	else
		return NULL;
}

mapkit_error map_strd_remove_s(map_strd * spm, char *key, mapkit_hash_t hash)
{
	mapkit_size_t iindex;

	iindex = map_strd_keyindex(spm, key, hash);

	if (iindex < 0) {
		if (spm->alwaysdefault)
			return MAPKIT_OK;
		else
			return MAPKIT_EKEYNOTFOUND;
	}

	spm->contents[iindex].state = MAPKIT_DELETEDSLOT;
	spm->used--;
	if (spm->used < spm->minused) {
#ifdef MAPKIT_DEBUG
		fprintf(stderr, "MAPKIT: used < minused\n");
#endif
		return map_strd_reallocate(spm, map_strd_shrinksize(spm, spm->used));
	}

	return MAPKIT_OK;
}

mapkit_size_t map_strd_keyindex(map_strd * spm, char *key, mapkit_hash_t hash)
{
	mapkit_size_t iindex, decrement;

	signed char state;

	iindex = hash % spm->size;
	decrement = (hash % (spm->size - 2));
	decrement += (decrement == 0);

#ifdef MAPKIT_COLLISIONS
	spm->keyindexs++;
#endif

	while ((state = spm->contents[iindex].state) != MAPKIT_FREESLOT
	       && (state == MAPKIT_DELETEDSLOT || (!(strcmp(spm->contents[iindex].key, key) == 0)))) {
#ifdef MAPKIT_COLLISIONS
		spm->keyindex_collisions++;
#endif
		iindex -= decrement;
		if (iindex < 0)
			iindex += spm->size;
	}

	if (state == MAPKIT_FREESLOT)
		return MAPKIT_KEYNOTFOUND;
	return iindex;
}

mapkit_size_t map_strd_insertionindex(map_strd * spm, char *key, mapkit_hash_t hash)
{
	mapkit_size_t iindex, decrement;
	signed char state;

#ifdef MAPKIT_COLLISIONS
	spm->insertionindexs++;
#endif

	iindex = hash % spm->size;

	/*
	 * Fast path (largely superfluous) 
	 */
	if ((state = spm->contents[iindex].state) == MAPKIT_FREESLOT)
		return iindex;
	if ((state == MAPKIT_FULLSLOT)
	    && (strcmp(spm->contents[iindex].key, key) == 0))
		return -iindex - 1;

	decrement = (hash % (spm->size - 2));
	decrement += (decrement == 0);

	while ((state == MAPKIT_FULLSLOT)
	       && (!(strcmp(spm->contents[iindex].key, key) == 0))) {
#ifdef MAPKIT_COLLISIONS
		spm->insertionindex_collisions++;
#endif
		iindex -= decrement;
		if (iindex < 0)
			iindex += spm->size;
		state = spm->contents[iindex].state;
	}

	if (state == MAPKIT_DELETEDSLOT) {
		mapkit_size_t index2 = iindex;

		while ((state = spm->contents[iindex].state) != MAPKIT_FREESLOT && ((state == MAPKIT_DELETEDSLOT)
										    || (!(strcmp(spm->contents[iindex].key, key) == 0)))) {
			iindex -= decrement;
			if (iindex < 0)
				iindex += spm->size;
			state = spm->contents[iindex].state;
		}
		if (state == MAPKIT_FREESLOT)
			return index2;
	}

	if (state == MAPKIT_FULLSLOT)
		return -iindex - 1;
	return iindex;
}

mapkit_size_t map_strd_next(map_strd * spm, mapkit_size_t iindex)
{
	mapkit_size_t size = spm->size;
	map_strd_storage *pos_contents;
	double defaultvalue = spm->defaultvalue;
	int notalwaysdefault = !spm->alwaysdefault;

	pos_contents = &(spm->contents[++iindex]);

	for (; iindex < size; iindex++, pos_contents++)
		if ((pos_contents->state == MAPKIT_FULLSLOT)
		    && (notalwaysdefault || (!((pos_contents->value) == (defaultvalue)))))
			return iindex;

	return -1;
}

map_strd_storage *map_strd_nextptr(map_strd * spm, map_strd_storage * pos_contents)
{
	map_strd_storage *end = &(spm->contents[spm->size]);
	double defaultvalue = spm->defaultvalue;
	int notalwaysdefault = !spm->alwaysdefault;

	if (pos_contents == NULL)
		pos_contents = spm->contents;
	else {
		pos_contents++;
		if (pos_contents <= spm->contents)
			return NULL;
	}

	for (; pos_contents < end; pos_contents++)
		if ((pos_contents->state == MAPKIT_FULLSLOT)
		    && (notalwaysdefault || (!((pos_contents->value) == (defaultvalue)))))
			return pos_contents;

	return NULL;
}

mapkit_error map_strd_getall(map_strd * spm, map_strd_element ** array, mapkit_size_t * count)
{
	mapkit_size_t iindex;
	mapkit_size_t size = spm->size, vcount = 0;
	map_strd_element *pos_array;
	map_strd_storage *pos_contents;
	double defaultvalue = spm->defaultvalue;
	int notalwaysdefault = !spm->alwaysdefault;

	pos_array = *array = (map_strd_element *) malloc(sizeof(**array) * spm->used);
	if (*array == NULL)
		MAPKIT_ERROR(MAPKIT_ENOMEM);

	pos_contents = spm->contents;

	for (iindex = 0; iindex < size; iindex++, pos_contents++)
		if ((pos_contents->state == MAPKIT_FULLSLOT)
		    && (notalwaysdefault || (!((pos_contents->value) == (defaultvalue))))) {
			pos_array->key = pos_contents->key;
			pos_array->value = pos_contents->value;
			pos_array++;
			vcount++;
		}
	*count = vcount;

	return MAPKIT_OK;
}

mapkit_error map_strd_clean(map_strd * spm)
{
	mapkit_size_t iindex, count = 0;
	mapkit_size_t size = spm->size;
	map_strd_storage *pos_contents;
	double defaultvalue = spm->defaultvalue;

	pos_contents = spm->contents;
	for (iindex = 0; iindex < size; iindex++, pos_contents++)
		if ((pos_contents->state == MAPKIT_FULLSLOT)
		    && ((pos_contents->value) == (defaultvalue))) {
			pos_contents->state = MAPKIT_DELETEDSLOT;
			count++;
		}

	spm->used -= count;
	if (spm->used < spm->minused) {
#ifdef MAPKIT_DEBUG
		fprintf(stderr, "MAPKIT: used < minused\n");
#endif
		return map_strd_reallocate(spm, map_strd_meansize(spm, spm->used));
	}

	return MAPKIT_OK;
}

int map_strd_compare(const void *e1, const void *e2)
{
	char *key1 = ((const map_strd_element *) e1)->key;
	char *key2 = ((const map_strd_element *) e2)->key;

	return strcmp(key1, key2);
}

mapkit_error map_strd_getall_sorted(map_strd * spm, map_strd_element ** array, mapkit_size_t * count)
{
	mapkit_error err;

	err = map_strd_getall(spm, array, count);
	if (err)
		MAPKIT_ERROR(err);

	qsort(*array, *count, sizeof(**array), map_strd_compare);

	return MAPKIT_OK;
}

mapkit_error map_strd_setall(map_strd * spm, map_strd_element * array, mapkit_size_t count)
{
	mapkit_size_t array_iindex;
	mapkit_error err;

	err = map_strd_ensurecapacity(spm, spm->used + count);
	if (err)
		MAPKIT_ERROR(err);

	if (spm->alwaysdefault)
		/*
		 * Prevent shrinking 
		 */
		spm->minused = 0;

	for (array_iindex = 0; array_iindex < count; array_iindex++) {
		err = map_strd_set(spm, array[array_iindex].key, array[array_iindex].value);
		if (err)
			MAPKIT_ERROR(err);
	}

	if (spm->alwaysdefault)
		map_strd_adjustcapacity(spm);

	return MAPKIT_OK;
}

mapkit_error map_strd_removeall(map_strd * spm, char **array, mapkit_size_t count)
{
	mapkit_size_t array_iindex;
	mapkit_error err;

	/*
	 * Prevent shrinking 
	 */
	spm->minused = 0;

	for (array_iindex = 0; array_iindex < count; array_iindex++) {
		err = map_strd_remove(spm, array[array_iindex]);
		if (err)
			MAPKIT_ERROR(err);
	}

	map_strd_adjustcapacity(spm);

	return MAPKIT_OK;
}

void map_strd_printstats(map_strd * spm)
{
	fprintf(stderr, "MAPKIT: map_strd statistics\n");
	fprintf(stderr, "MAPKIT: alwaysdefault = %d\n", spm->alwaysdefault);
	fprintf(stderr, "MAPKIT: minused = %ld, maxfill = %ld\n", (long) spm->minused, (long) spm->maxfill);
	fprintf(stderr, "MAPKIT: minusedfactor = %g, maxfillfactor = %g\n", spm->minusedfactor, spm->maxfillfactor);
#ifdef MAPKIT_COLLISIONS
	fprintf(stderr, "MAPKIT: insertionindexs = %lu, collisions = %lu\n", (unsigned long) spm->insertionindexs,
		(unsigned long) spm->insertionindex_collisions);
	fprintf(stderr, "MAPKIT: keyindexs = %lu, collisions = %lu\n", (unsigned long) spm->keyindexs, (unsigned long) spm->keyindex_collisions);
#endif
}

#endif							       /* MAPKIT_map_strd */

#ifdef MAPKIT_map_strvp

/*
  Implementation for map_strvp (char* -> void*)
  Default value : NULL
  Uses a state field.
*/

/* Static prototypes */

/* Return the iindex of key, or MAPKIT_KEYNOTFOUND */
static mapkit_size_t map_strvp_keyindex(map_strvp * spm, char *key, mapkit_hash_t hash);

/* Return the iindex of key or -(insertion iindex)-1 if key not found */
static mapkit_size_t map_strvp_insertionindex(map_strvp * spm, char *key, mapkit_hash_t hash);

/* Implementation */

mapkit_error map_strvp_init(map_strvp * spm)
{
	return map_strvp_init_hint(spm, MAPKIT_DEFAULT_EXPECTEDUSED);
}

mapkit_error map_strvp_init_hint(map_strvp * spm, mapkit_size_t used)
{
#ifdef MAPKIT_DEBUG
	fprintf(stderr, "MAPKIT: init\n");
#endif

	spm->size = 0;
	spm->fill = 0;
	spm->used = 0;
	spm->maxfillfactor = 0.5;
	spm->minusedfactor = 0.2;
	spm->contents = NULL;
#ifdef MAPKIT_COLLISIONS
	spm->insertionindexs = spm->insertionindex_collisions = 0;
	spm->keyindexs = spm->keyindex_collisions = 0;
#endif
	spm->defaultvalue = NULL;
	spm->alwaysdefault = 0;

	return map_strvp_reallocate(spm, map_strvp_meansize(spm, used));
}

mapkit_error map_strvp_ensurecapacity(map_strvp * spm, mapkit_size_t used)
{
	if (used > (spm->used + spm->maxfill - spm->fill)) {
#ifdef MAPKIT_DEBUG
		fprintf(stderr, "MAPKIT: need more capacity\n");
#endif
		return map_strvp_reallocate(spm, map_strvp_meansize(spm, used));
	} else
		return MAPKIT_OK;
}

mapkit_error map_strvp_adjustcapacity(map_strvp * spm)
{
	spm->minused = (mapkit_size_t) (spm->size * spm->minusedfactor);
	spm->maxfill = (mapkit_size_t) (spm->size * spm->maxfillfactor);

	if (spm->used < spm->minused) {
#ifdef MAPKIT_DEBUG
		fprintf(stderr, "MAPKIT: used < minused\n");
#endif
		return map_strvp_reallocate(spm, map_strvp_meansize(spm, spm->used));
	} else if (spm->fill > spm->maxfill) {
#ifdef MAPKIT_DEBUG
		fprintf(stderr, "MAPKIT: fill > maxfill\n");
#endif
		return map_strvp_reallocate(spm, map_strvp_meansize(spm, spm->used));
	} else
		return MAPKIT_OK;
}

void map_strvp_free(map_strvp * spm)
{
	free(spm->contents);
	spm->contents = NULL;
	spm->size = 0;
	spm->fill = 0;
	spm->used = 0;
	spm->maxfill = 0;
	spm->minused = 0;
#ifdef MAPKIT_COLLISIONS
	spm->insertionindexs = spm->insertionindex_collisions = 0;
	spm->keyindexs = spm->keyindex_collisions = 0;
#endif
}

mapkit_error map_strvp_copy(map_strvp * to, map_strvp * from)
{
	map_strvp_storage *contentscopy;
	size_t size = from->size * sizeof(*from->contents);

	contentscopy = (map_strvp_storage *) malloc(size);
	if (contentscopy == NULL)
		MAPKIT_ERROR(MAPKIT_ENOMEM);

	memcpy(to, from, sizeof(map_strvp));
	to->contents = contentscopy;
	memcpy(to->contents, from->contents, size);

	return MAPKIT_OK;
}

mapkit_size_t map_strvp_growsize(map_strvp * spm, mapkit_size_t used)
{
	return (mapkit_size_t) (4.0 * used / (3 * spm->minusedfactor + spm->maxfillfactor));
}

mapkit_size_t map_strvp_shrinksize(map_strvp * spm, mapkit_size_t used)
{
	return (mapkit_size_t) (4.0 * used / (spm->minusedfactor + 3 * spm->maxfillfactor));
}

mapkit_size_t map_strvp_meansize(map_strvp * spm, mapkit_size_t used)
{
	return (mapkit_size_t) (2.0 * used / (spm->minusedfactor + spm->maxfillfactor));
}

mapkit_error map_strvp_reallocate(map_strvp * spm, mapkit_size_t newsize)
{
	mapkit_size_t iindex;
	mapkit_size_t oldsize;
	map_strvp_storage *newcontents, *oldcontents;

	/*
	 * At least one free entry 
	 */
	if (newsize <= spm->used)
		newsize = spm->used + 1;
	newsize = mapkit_nextprime(newsize);
	if (newsize <= spm->used)
		MAPKIT_ERROR(MAPKIT_ETOOBIG);

	newcontents = (map_strvp_storage *) malloc(newsize * sizeof(*spm->contents));
	if (newcontents == NULL)
		MAPKIT_ERROR(MAPKIT_ENOMEM);

	/*
	 * Initialize all entries to "free" 
	 */
	for (iindex = 0; iindex < newsize; iindex++)
		newcontents[iindex].state = MAPKIT_FREESLOT;

	oldcontents = spm->contents;
	oldsize = spm->size;
	spm->contents = newcontents;
	spm->size = newsize;

#ifdef MAPKIT_DEBUG
	fprintf(stderr, "MAPKIT: reallocate %ld -> %ld\n", (long) oldsize, (long) newsize);
#endif

	spm->maxfill = (mapkit_size_t) (newsize * spm->maxfillfactor);
	/*
	 * At least one free entry 
	 */
	if (spm->maxfill >= newsize)
		spm->maxfill = newsize - 1;
	spm->minused = (mapkit_size_t) (newsize * spm->minusedfactor);
	spm->used = 0;

	if (oldcontents != NULL) {
		int used = 0;
		char *key;
		void *defaultvalue = spm->defaultvalue;
		int notalwaysdefault = !spm->alwaysdefault;

		/*
		 * Copy all entries from old to new 
		 */
		for (iindex = 0; iindex < oldsize; iindex++)
			if (oldcontents[iindex].state == MAPKIT_FULLSLOT) {
				mapkit_size_t ins_iindex;
				mapkit_hash_t hash;
				map_strvp_storage *contents;

				key = oldcontents[iindex].key;

				/*
				 * Fast path 
				 */
				ins_iindex = (hash = mapkit_strhash(key)) % spm->size;
				contents = &(newcontents[ins_iindex]);

				if (contents->state != MAPKIT_FREESLOT) {
					ins_iindex = map_strvp_insertionindex(spm, key, hash);
					contents = &(newcontents[ins_iindex]);
				}
#ifdef MAPKIT_COLLISIONS
				else
					spm->insertionindexs++;
#endif
				if (notalwaysdefault || (!((oldcontents[iindex].value) == (defaultvalue)))) {
					contents->value = oldcontents[iindex].value;
					contents->state = MAPKIT_FULLSLOT;
					contents->key = key;
					used++;
				}
			}
		free(oldcontents);
		spm->used = used;
	}
	spm->fill = spm->used;

	return MAPKIT_OK;
}

void *map_strvp_value_s(map_strvp * spm, char *key, mapkit_hash_t hash)
{
	mapkit_size_t iindex;

	iindex = map_strvp_keyindex(spm, key, hash);

	if (iindex < 0) {
		if (spm->alwaysdefault)
			return spm->defaultvalue;
		else
			MAPKIT_FATAL_ERROR(MAPKIT_EKEYNOTFOUND);
	}

	return spm->contents[iindex].value;
}

mapkit_error map_strvp_get_s(map_strvp * spm, char *key, void **value, mapkit_hash_t hash)
{
	mapkit_size_t iindex;

	iindex = map_strvp_keyindex(spm, key, hash);

	if (iindex < 0) {
		if (spm->alwaysdefault) {
			*value = spm->defaultvalue;
			return MAPKIT_OK;
		} else
			MAPKIT_ERROR(MAPKIT_EKEYNOTFOUND);
	}

	*value = spm->contents[iindex].value;
	return MAPKIT_OK;
}

mapkit_error map_strvp_set_s(map_strvp * spm, char *key, void *value, mapkit_hash_t hash)
{
	mapkit_size_t iindex;

	iindex = map_strvp_insertionindex(spm, key, hash);

	if (iindex < 0)
		/*
		 * FULLSLOT 
		 */
		spm->contents[-iindex - 1].value = value;
	else {
		map_strvp_storage *element = &(spm->contents[iindex]);
		int ffree = element->state == MAPKIT_FREESLOT;

		element->state = MAPKIT_FULLSLOT;
		element->value = value;
		element->key = key;
		spm->used++;

		if (ffree && ((++spm->fill) > spm->maxfill)) {
#ifdef MAPKIT_DEBUG
			fprintf(stderr, "MAPKIT: fill > maxfill\n");
#endif
			return map_strvp_reallocate(spm, map_strvp_growsize(spm, spm->used));
		}
	}
	return MAPKIT_OK;
}

void **map_strvp_insertptr_s(map_strvp * spm, char *key, mapkit_hash_t hash)
{
	mapkit_size_t iindex;

	iindex = map_strvp_insertionindex(spm, key, hash);

	if (iindex < 0)
		return &(spm->contents[-iindex - 1].value);
	else {
		map_strvp_storage *element = &(spm->contents[iindex]);

		if (element->state == MAPKIT_FREESLOT) {
			/*
			 * FREESLOT 
			 */
			if (spm->fill >= spm->maxfill) {
				mapkit_error err;

#ifdef MAPKIT_DEBUG
				fprintf(stderr, "MAPKIT: fill => maxfill before insert\n");
#endif
				/*
				 * Must reallocate -before- inserting defaultvalue 
				 */
				err = map_strvp_reallocate(spm, map_strvp_growsize(spm, spm->used + 1));
				if (err) {
					MAPKIT_ERROR_NORET(err);
					return NULL;
				}

				iindex = map_strvp_insertionindex(spm, key, hash);
				/*
				 * FREESLOT 
				 */
				spm->contents[iindex].state = MAPKIT_FULLSLOT;
				spm->contents[iindex].key = key;
				spm->contents[iindex].value = spm->defaultvalue;
				spm->used++;
				spm->fill++;
				return &(spm->contents[iindex].value);
			} else
				spm->fill++;
		}

		element->state = MAPKIT_FULLSLOT;
		element->key = key;
		element->value = spm->defaultvalue;
		spm->used++;
		return &(element->value);
	}
}

void **map_strvp_ptr_s(map_strvp * spm, char *key, mapkit_hash_t hash)
{
	mapkit_size_t iindex;

	iindex = map_strvp_keyindex(spm, key, hash);

	if ((iindex >= 0) && ((!spm->alwaysdefault) || (!((spm->contents[iindex].value) == (spm->defaultvalue)))))
		return &(spm->contents[iindex].value);
	else
		return NULL;
}

mapkit_error map_strvp_remove_s(map_strvp * spm, char *key, mapkit_hash_t hash)
{
	mapkit_size_t iindex;

	iindex = map_strvp_keyindex(spm, key, hash);

	if (iindex < 0) {
		if (spm->alwaysdefault)
			return MAPKIT_OK;
		else
			return MAPKIT_EKEYNOTFOUND;
	}

	spm->contents[iindex].state = MAPKIT_DELETEDSLOT;
	spm->used--;
	if (spm->used < spm->minused) {
#ifdef MAPKIT_DEBUG
		fprintf(stderr, "MAPKIT: used < minused\n");
#endif
		return map_strvp_reallocate(spm, map_strvp_shrinksize(spm, spm->used));
	}

	return MAPKIT_OK;
}

mapkit_size_t map_strvp_keyindex(map_strvp * spm, char *key, mapkit_hash_t hash)
{
	mapkit_size_t iindex, decrement;

	signed char state;

	iindex = hash % spm->size;
	decrement = (hash % (spm->size - 2));
	decrement += (decrement == 0);

#ifdef MAPKIT_COLLISIONS
	spm->keyindexs++;
#endif

	while ((state = spm->contents[iindex].state) != MAPKIT_FREESLOT
	       && (state == MAPKIT_DELETEDSLOT || (!(strcmp(spm->contents[iindex].key, key) == 0)))) {
#ifdef MAPKIT_COLLISIONS
		spm->keyindex_collisions++;
#endif
		iindex -= decrement;
		if (iindex < 0)
			iindex += spm->size;
	}

	if (state == MAPKIT_FREESLOT)
		return MAPKIT_KEYNOTFOUND;
	return iindex;
}

mapkit_size_t map_strvp_insertionindex(map_strvp * spm, char *key, mapkit_hash_t hash)
{
	mapkit_size_t iindex, decrement;
	signed char state;

#ifdef MAPKIT_COLLISIONS
	spm->insertionindexs++;
#endif

	iindex = hash % spm->size;

	/*
	 * Fast path (largely superfluous) 
	 */
	if ((state = spm->contents[iindex].state) == MAPKIT_FREESLOT)
		return iindex;
	if ((state == MAPKIT_FULLSLOT)
	    && (strcmp(spm->contents[iindex].key, key) == 0))
		return -iindex - 1;

	decrement = (hash % (spm->size - 2));
	decrement += (decrement == 0);

	while ((state == MAPKIT_FULLSLOT)
	       && (!(strcmp(spm->contents[iindex].key, key) == 0))) {
#ifdef MAPKIT_COLLISIONS
		spm->insertionindex_collisions++;
#endif
		iindex -= decrement;
		if (iindex < 0)
			iindex += spm->size;
		state = spm->contents[iindex].state;
	}

	if (state == MAPKIT_DELETEDSLOT) {
		mapkit_size_t index2 = iindex;

		while ((state = spm->contents[iindex].state) != MAPKIT_FREESLOT && ((state == MAPKIT_DELETEDSLOT)
										    || (!(strcmp(spm->contents[iindex].key, key) == 0)))) {
			iindex -= decrement;
			if (iindex < 0)
				iindex += spm->size;
			state = spm->contents[iindex].state;
		}
		if (state == MAPKIT_FREESLOT)
			return index2;
	}

	if (state == MAPKIT_FULLSLOT)
		return -iindex - 1;
	return iindex;
}

mapkit_size_t map_strvp_next(map_strvp * spm, mapkit_size_t iindex)
{
	mapkit_size_t size = spm->size;
	map_strvp_storage *pos_contents;
	void *defaultvalue = spm->defaultvalue;
	int notalwaysdefault = !spm->alwaysdefault;

	pos_contents = &(spm->contents[++iindex]);

	for (; iindex < size; iindex++, pos_contents++)
		if ((pos_contents->state == MAPKIT_FULLSLOT)
		    && (notalwaysdefault || (!((pos_contents->value) == (defaultvalue)))))
			return iindex;

	return -1;
}

map_strvp_storage *map_strvp_nextptr(map_strvp * spm, map_strvp_storage * pos_contents)
{
	map_strvp_storage *end = &(spm->contents[spm->size]);
	void *defaultvalue = spm->defaultvalue;
	int notalwaysdefault = !spm->alwaysdefault;

	if (pos_contents == NULL)
		pos_contents = spm->contents;
	else {
		pos_contents++;
		if (pos_contents <= spm->contents)
			return NULL;
	}

	for (; pos_contents < end; pos_contents++)
		if ((pos_contents->state == MAPKIT_FULLSLOT)
		    && (notalwaysdefault || (!((pos_contents->value) == (defaultvalue)))))
			return pos_contents;

	return NULL;
}

mapkit_error map_strvp_getall(map_strvp * spm, map_strvp_element ** array, mapkit_size_t * count)
{
	mapkit_size_t iindex;
	mapkit_size_t size = spm->size, vcount = 0;
	map_strvp_element *pos_array;
	map_strvp_storage *pos_contents;
	void *defaultvalue = spm->defaultvalue;
	int notalwaysdefault = !spm->alwaysdefault;

	pos_array = *array = (map_strvp_element *) malloc(sizeof(**array) * spm->used);
	if (*array == NULL)
		MAPKIT_ERROR(MAPKIT_ENOMEM);

	pos_contents = spm->contents;

	for (iindex = 0; iindex < size; iindex++, pos_contents++)
		if ((pos_contents->state == MAPKIT_FULLSLOT)
		    && (notalwaysdefault || (!((pos_contents->value) == (defaultvalue))))) {
			pos_array->key = pos_contents->key;
			pos_array->value = pos_contents->value;
			pos_array++;
			vcount++;
		}
	*count = vcount;

	return MAPKIT_OK;
}

mapkit_error map_strvp_clean(map_strvp * spm)
{
	mapkit_size_t iindex, count = 0;
	mapkit_size_t size = spm->size;
	map_strvp_storage *pos_contents;
	void *defaultvalue = spm->defaultvalue;

	pos_contents = spm->contents;
	for (iindex = 0; iindex < size; iindex++, pos_contents++)
		if ((pos_contents->state == MAPKIT_FULLSLOT)
		    && ((pos_contents->value) == (defaultvalue))) {
			pos_contents->state = MAPKIT_DELETEDSLOT;
			count++;
		}

	spm->used -= count;
	if (spm->used < spm->minused) {
#ifdef MAPKIT_DEBUG
		fprintf(stderr, "MAPKIT: used < minused\n");
#endif
		return map_strvp_reallocate(spm, map_strvp_meansize(spm, spm->used));
	}

	return MAPKIT_OK;
}

int map_strvp_compare(const void *e1, const void *e2)
{
	char *key1 = ((const map_strvp_element *) e1)->key;
	char *key2 = ((const map_strvp_element *) e2)->key;

	return strcmp(key1, key2);
}

mapkit_error map_strvp_getall_sorted(map_strvp * spm, map_strvp_element ** array, mapkit_size_t * count)
{
	mapkit_error err;

	err = map_strvp_getall(spm, array, count);
	if (err)
		MAPKIT_ERROR(err);

	qsort(*array, *count, sizeof(**array), map_strvp_compare);

	return MAPKIT_OK;
}

mapkit_error map_strvp_setall(map_strvp * spm, map_strvp_element * array, mapkit_size_t count)
{
	mapkit_size_t array_iindex;
	mapkit_error err;

	err = map_strvp_ensurecapacity(spm, spm->used + count);
	if (err)
		MAPKIT_ERROR(err);

	if (spm->alwaysdefault)
		/*
		 * Prevent shrinking 
		 */
		spm->minused = 0;

	for (array_iindex = 0; array_iindex < count; array_iindex++) {
		err = map_strvp_set(spm, array[array_iindex].key, array[array_iindex].value);
		if (err)
			MAPKIT_ERROR(err);
	}

	if (spm->alwaysdefault)
		map_strvp_adjustcapacity(spm);

	return MAPKIT_OK;
}

mapkit_error map_strvp_removeall(map_strvp * spm, char **array, mapkit_size_t count)
{
	mapkit_size_t array_iindex;
	mapkit_error err;

	/*
	 * Prevent shrinking 
	 */
	spm->minused = 0;

	for (array_iindex = 0; array_iindex < count; array_iindex++) {
		err = map_strvp_remove(spm, array[array_iindex]);
		if (err)
			MAPKIT_ERROR(err);
	}

	map_strvp_adjustcapacity(spm);

	return MAPKIT_OK;
}

void map_strvp_printstats(map_strvp * spm)
{
	fprintf(stderr, "MAPKIT: map_strvp statistics\n");
	fprintf(stderr, "MAPKIT: alwaysdefault = %d\n", spm->alwaysdefault);
	fprintf(stderr, "MAPKIT: minused = %ld, maxfill = %ld\n", (long) spm->minused, (long) spm->maxfill);
	fprintf(stderr, "MAPKIT: minusedfactor = %g, maxfillfactor = %g\n", spm->minusedfactor, spm->maxfillfactor);
#ifdef MAPKIT_COLLISIONS
	fprintf(stderr, "MAPKIT: insertionindexs = %lu, collisions = %lu\n", (unsigned long) spm->insertionindexs,
		(unsigned long) spm->insertionindex_collisions);
	fprintf(stderr, "MAPKIT: keyindexs = %lu, collisions = %lu\n", (unsigned long) spm->keyindexs, (unsigned long) spm->keyindex_collisions);
#endif
}

#endif							       /* MAPKIT_map_strvp */

#ifdef MAPKIT_map_strstr

/*
  Implementation for map_strstr (char* -> char*)
  Default value : ""
  Uses a state field.
*/

/* Static prototypes */

/* Return the iindex of key, or MAPKIT_KEYNOTFOUND */
static mapkit_size_t map_strstr_keyindex(map_strstr * spm, char *key, mapkit_hash_t hash);

/* Return the iindex of key or -(insertion iindex)-1 if key not found */
static mapkit_size_t map_strstr_insertionindex(map_strstr * spm, char *key, mapkit_hash_t hash);

/* Implementation */

mapkit_error map_strstr_init(map_strstr * spm)
{
	return map_strstr_init_hint(spm, MAPKIT_DEFAULT_EXPECTEDUSED);
}

mapkit_error map_strstr_init_hint(map_strstr * spm, mapkit_size_t used)
{
#ifdef MAPKIT_DEBUG
	fprintf(stderr, "MAPKIT: init\n");
#endif

	spm->size = 0;
	spm->fill = 0;
	spm->used = 0;
	spm->maxfillfactor = 0.5;
	spm->minusedfactor = 0.2;
	spm->contents = NULL;
#ifdef MAPKIT_COLLISIONS
	spm->insertionindexs = spm->insertionindex_collisions = 0;
	spm->keyindexs = spm->keyindex_collisions = 0;
#endif
	spm->defaultvalue = GMRFLib_strdup("");
	spm->alwaysdefault = 0;

	return map_strstr_reallocate(spm, map_strstr_meansize(spm, used));
}

mapkit_error map_strstr_ensurecapacity(map_strstr * spm, mapkit_size_t used)
{
	if (used > (spm->used + spm->maxfill - spm->fill)) {
#ifdef MAPKIT_DEBUG
		fprintf(stderr, "MAPKIT: need more capacity\n");
#endif
		return map_strstr_reallocate(spm, map_strstr_meansize(spm, used));
	} else
		return MAPKIT_OK;
}

mapkit_error map_strstr_adjustcapacity(map_strstr * spm)
{
	spm->minused = (mapkit_size_t) (spm->size * spm->minusedfactor);
	spm->maxfill = (mapkit_size_t) (spm->size * spm->maxfillfactor);

	if (spm->used < spm->minused) {
#ifdef MAPKIT_DEBUG
		fprintf(stderr, "MAPKIT: used < minused\n");
#endif
		return map_strstr_reallocate(spm, map_strstr_meansize(spm, spm->used));
	} else if (spm->fill > spm->maxfill) {
#ifdef MAPKIT_DEBUG
		fprintf(stderr, "MAPKIT: fill > maxfill\n");
#endif
		return map_strstr_reallocate(spm, map_strstr_meansize(spm, spm->used));
	} else
		return MAPKIT_OK;
}

void map_strstr_free(map_strstr * spm)
{
	free(spm->contents);
	spm->contents = NULL;
	spm->size = 0;
	spm->fill = 0;
	spm->used = 0;
	spm->maxfill = 0;
	spm->minused = 0;
#ifdef MAPKIT_COLLISIONS
	spm->insertionindexs = spm->insertionindex_collisions = 0;
	spm->keyindexs = spm->keyindex_collisions = 0;
#endif
}

mapkit_error map_strstr_copy(map_strstr * to, map_strstr * from)
{
	map_strstr_storage *contentscopy;
	size_t size = from->size * sizeof(*from->contents);

	contentscopy = (map_strstr_storage *) malloc(size);
	if (contentscopy == NULL)
		MAPKIT_ERROR(MAPKIT_ENOMEM);

	memcpy(to, from, sizeof(map_strstr));
	to->contents = contentscopy;
	memcpy(to->contents, from->contents, size);

	return MAPKIT_OK;
}

mapkit_size_t map_strstr_growsize(map_strstr * spm, mapkit_size_t used)
{
	return (mapkit_size_t) (4.0 * used / (3 * spm->minusedfactor + spm->maxfillfactor));
}

mapkit_size_t map_strstr_shrinksize(map_strstr * spm, mapkit_size_t used)
{
	return (mapkit_size_t) (4.0 * used / (spm->minusedfactor + 3 * spm->maxfillfactor));
}

mapkit_size_t map_strstr_meansize(map_strstr * spm, mapkit_size_t used)
{
	return (mapkit_size_t) (2.0 * used / (spm->minusedfactor + spm->maxfillfactor));
}

mapkit_error map_strstr_reallocate(map_strstr * spm, mapkit_size_t newsize)
{
	mapkit_size_t iindex;
	mapkit_size_t oldsize;
	map_strstr_storage *newcontents, *oldcontents;

	/*
	 * At least one free entry 
	 */
	if (newsize <= spm->used)
		newsize = spm->used + 1;
	newsize = mapkit_nextprime(newsize);
	if (newsize <= spm->used)
		MAPKIT_ERROR(MAPKIT_ETOOBIG);

	newcontents = (map_strstr_storage *) malloc(newsize * sizeof(*spm->contents));
	if (newcontents == NULL)
		MAPKIT_ERROR(MAPKIT_ENOMEM);

	/*
	 * Initialize all entries to "free" 
	 */
	for (iindex = 0; iindex < newsize; iindex++)
		newcontents[iindex].state = MAPKIT_FREESLOT;

	oldcontents = spm->contents;
	oldsize = spm->size;
	spm->contents = newcontents;
	spm->size = newsize;

#ifdef MAPKIT_DEBUG
	fprintf(stderr, "MAPKIT: reallocate %ld -> %ld\n", (long) oldsize, (long) newsize);
#endif

	spm->maxfill = (mapkit_size_t) (newsize * spm->maxfillfactor);
	/*
	 * At least one free entry 
	 */
	if (spm->maxfill >= newsize)
		spm->maxfill = newsize - 1;
	spm->minused = (mapkit_size_t) (newsize * spm->minusedfactor);
	spm->used = 0;

	if (oldcontents != NULL) {
		int used = 0;
		char *key;
		char *defaultvalue = spm->defaultvalue;
		int notalwaysdefault = !spm->alwaysdefault;

		/*
		 * Copy all entries from old to new 
		 */
		for (iindex = 0; iindex < oldsize; iindex++)
			if (oldcontents[iindex].state == MAPKIT_FULLSLOT) {
				mapkit_size_t ins_iindex;
				mapkit_hash_t hash;
				map_strstr_storage *contents;

				key = oldcontents[iindex].key;

				/*
				 * Fast path 
				 */
				ins_iindex = (hash = mapkit_strhash(key)) % spm->size;
				contents = &(newcontents[ins_iindex]);

				if (contents->state != MAPKIT_FREESLOT) {
					ins_iindex = map_strstr_insertionindex(spm, key, hash);
					contents = &(newcontents[ins_iindex]);
				}
#ifdef MAPKIT_COLLISIONS
				else
					spm->insertionindexs++;
#endif
				if (notalwaysdefault || (!(strcmp(oldcontents[iindex].value, defaultvalue) == 0))) {
					contents->value = oldcontents[iindex].value;
					contents->state = MAPKIT_FULLSLOT;
					contents->key = key;
					used++;
				}
			}
		free(oldcontents);
		spm->used = used;
	}
	spm->fill = spm->used;

	return MAPKIT_OK;
}

char *map_strstr_value_s(map_strstr * spm, char *key, mapkit_hash_t hash)
{
	mapkit_size_t iindex;

	iindex = map_strstr_keyindex(spm, key, hash);

	if (iindex < 0) {
		if (spm->alwaysdefault)
			return spm->defaultvalue;
		else
			MAPKIT_FATAL_ERROR(MAPKIT_EKEYNOTFOUND);
	}

	return spm->contents[iindex].value;
}

mapkit_error map_strstr_get_s(map_strstr * spm, char *key, char **value, mapkit_hash_t hash)
{
	mapkit_size_t iindex;

	iindex = map_strstr_keyindex(spm, key, hash);

	if (iindex < 0) {
		if (spm->alwaysdefault) {
			*value = spm->defaultvalue;
			return MAPKIT_OK;
		} else
			MAPKIT_ERROR(MAPKIT_EKEYNOTFOUND);
	}

	*value = spm->contents[iindex].value;
	return MAPKIT_OK;
}

mapkit_error map_strstr_set_s(map_strstr * spm, char *key, char *value, mapkit_hash_t hash)
{
	mapkit_size_t iindex;

	iindex = map_strstr_insertionindex(spm, key, hash);

	if (iindex < 0)
		/*
		 * FULLSLOT 
		 */
		spm->contents[-iindex - 1].value = value;
	else {
		map_strstr_storage *element = &(spm->contents[iindex]);
		int ffree = element->state == MAPKIT_FREESLOT;

		element->state = MAPKIT_FULLSLOT;
		element->value = value;
		element->key = key;
		spm->used++;

		if (ffree && ((++spm->fill) > spm->maxfill)) {
#ifdef MAPKIT_DEBUG
			fprintf(stderr, "MAPKIT: fill > maxfill\n");
#endif
			return map_strstr_reallocate(spm, map_strstr_growsize(spm, spm->used));
		}
	}
	return MAPKIT_OK;
}

char **map_strstr_insertptr_s(map_strstr * spm, char *key, mapkit_hash_t hash)
{
	mapkit_size_t iindex;

	iindex = map_strstr_insertionindex(spm, key, hash);

	if (iindex < 0)
		return &(spm->contents[-iindex - 1].value);
	else {
		map_strstr_storage *element = &(spm->contents[iindex]);

		if (element->state == MAPKIT_FREESLOT) {
			/*
			 * FREESLOT 
			 */
			if (spm->fill >= spm->maxfill) {
				mapkit_error err;

#ifdef MAPKIT_DEBUG
				fprintf(stderr, "MAPKIT: fill => maxfill before insert\n");
#endif
				/*
				 * Must reallocate -before- inserting defaultvalue 
				 */
				err = map_strstr_reallocate(spm, map_strstr_growsize(spm, spm->used + 1));
				if (err) {
					MAPKIT_ERROR_NORET(err);
					return NULL;
				}

				iindex = map_strstr_insertionindex(spm, key, hash);
				/*
				 * FREESLOT 
				 */
				spm->contents[iindex].state = MAPKIT_FULLSLOT;
				spm->contents[iindex].key = key;
				spm->contents[iindex].value = spm->defaultvalue;
				spm->used++;
				spm->fill++;
				return &(spm->contents[iindex].value);
			} else
				spm->fill++;
		}

		element->state = MAPKIT_FULLSLOT;
		element->key = key;
		element->value = spm->defaultvalue;
		spm->used++;
		return &(element->value);
	}
}

char **map_strstr_ptr_s(map_strstr * spm, char *key, mapkit_hash_t hash)
{
	mapkit_size_t iindex;

	iindex = map_strstr_keyindex(spm, key, hash);

	if ((iindex >= 0) && ((!spm->alwaysdefault) || (!(strcmp(spm->contents[iindex].value, spm->defaultvalue) == 0))))
		return &(spm->contents[iindex].value);
	else
		return NULL;
}

mapkit_error map_strstr_remove_s(map_strstr * spm, char *key, mapkit_hash_t hash)
{
	mapkit_size_t iindex;

	iindex = map_strstr_keyindex(spm, key, hash);

	if (iindex < 0) {
		if (spm->alwaysdefault)
			return MAPKIT_OK;
		else
			return MAPKIT_EKEYNOTFOUND;
	}

	spm->contents[iindex].state = MAPKIT_DELETEDSLOT;
	spm->used--;
	if (spm->used < spm->minused) {
#ifdef MAPKIT_DEBUG
		fprintf(stderr, "MAPKIT: used < minused\n");
#endif
		return map_strstr_reallocate(spm, map_strstr_shrinksize(spm, spm->used));
	}

	return MAPKIT_OK;
}

mapkit_size_t map_strstr_keyindex(map_strstr * spm, char *key, mapkit_hash_t hash)
{
	mapkit_size_t iindex, decrement;

	signed char state;

	iindex = hash % spm->size;
	decrement = (hash % (spm->size - 2));
	decrement += (decrement == 0);

#ifdef MAPKIT_COLLISIONS
	spm->keyindexs++;
#endif

	while ((state = spm->contents[iindex].state) != MAPKIT_FREESLOT
	       && (state == MAPKIT_DELETEDSLOT || (!(strcmp(spm->contents[iindex].key, key) == 0)))) {
#ifdef MAPKIT_COLLISIONS
		spm->keyindex_collisions++;
#endif
		iindex -= decrement;
		if (iindex < 0)
			iindex += spm->size;
	}

	if (state == MAPKIT_FREESLOT)
		return MAPKIT_KEYNOTFOUND;
	return iindex;
}

mapkit_size_t map_strstr_insertionindex(map_strstr * spm, char *key, mapkit_hash_t hash)
{
	mapkit_size_t iindex, decrement;
	signed char state;

#ifdef MAPKIT_COLLISIONS
	spm->insertionindexs++;
#endif

	iindex = hash % spm->size;

	/*
	 * Fast path (largely superfluous) 
	 */
	if ((state = spm->contents[iindex].state) == MAPKIT_FREESLOT)
		return iindex;
	if ((state == MAPKIT_FULLSLOT)
	    && (strcmp(spm->contents[iindex].key, key) == 0))
		return -iindex - 1;

	decrement = (hash % (spm->size - 2));
	decrement += (decrement == 0);

	while ((state == MAPKIT_FULLSLOT)
	       && (!(strcmp(spm->contents[iindex].key, key) == 0))) {
#ifdef MAPKIT_COLLISIONS
		spm->insertionindex_collisions++;
#endif
		iindex -= decrement;
		if (iindex < 0)
			iindex += spm->size;
		state = spm->contents[iindex].state;
	}

	if (state == MAPKIT_DELETEDSLOT) {
		mapkit_size_t index2 = iindex;

		while ((state = spm->contents[iindex].state) != MAPKIT_FREESLOT && ((state == MAPKIT_DELETEDSLOT)
										    || (!(strcmp(spm->contents[iindex].key, key) == 0)))) {
			iindex -= decrement;
			if (iindex < 0)
				iindex += spm->size;
			state = spm->contents[iindex].state;
		}
		if (state == MAPKIT_FREESLOT)
			return index2;
	}

	if (state == MAPKIT_FULLSLOT)
		return -iindex - 1;
	return iindex;
}

mapkit_size_t map_strstr_next(map_strstr * spm, mapkit_size_t iindex)
{
	mapkit_size_t size = spm->size;
	map_strstr_storage *pos_contents;
	char *defaultvalue = spm->defaultvalue;
	int notalwaysdefault = !spm->alwaysdefault;

	pos_contents = &(spm->contents[++iindex]);

	for (; iindex < size; iindex++, pos_contents++)
		if ((pos_contents->state == MAPKIT_FULLSLOT)
		    && (notalwaysdefault || (!(strcmp(pos_contents->value, defaultvalue) == 0))))
			return iindex;

	return -1;
}

map_strstr_storage *map_strstr_nextptr(map_strstr * spm, map_strstr_storage * pos_contents)
{
	map_strstr_storage *end = &(spm->contents[spm->size]);
	char *defaultvalue = spm->defaultvalue;
	int notalwaysdefault = !spm->alwaysdefault;

	if (pos_contents == NULL)
		pos_contents = spm->contents;
	else {
		pos_contents++;
		if (pos_contents <= spm->contents)
			return NULL;
	}

	for (; pos_contents < end; pos_contents++)
		if ((pos_contents->state == MAPKIT_FULLSLOT)
		    && (notalwaysdefault || (!(strcmp(pos_contents->value, defaultvalue) == 0))))
			return pos_contents;

	return NULL;
}

mapkit_error map_strstr_getall(map_strstr * spm, map_strstr_element ** array, mapkit_size_t * count)
{
	mapkit_size_t iindex;
	mapkit_size_t size = spm->size, vcount = 0;
	map_strstr_element *pos_array;
	map_strstr_storage *pos_contents;
	char *defaultvalue = spm->defaultvalue;
	int notalwaysdefault = !spm->alwaysdefault;

	pos_array = *array = (map_strstr_element *) malloc(sizeof(**array) * spm->used);
	if (*array == NULL)
		MAPKIT_ERROR(MAPKIT_ENOMEM);

	pos_contents = spm->contents;

	for (iindex = 0; iindex < size; iindex++, pos_contents++)
		if ((pos_contents->state == MAPKIT_FULLSLOT)
		    && (notalwaysdefault || (!(strcmp(pos_contents->value, defaultvalue) == 0)))) {
			pos_array->key = pos_contents->key;
			pos_array->value = pos_contents->value;
			pos_array++;
			vcount++;
		}
	*count = vcount;

	return MAPKIT_OK;
}

mapkit_error map_strstr_clean(map_strstr * spm)
{
	mapkit_size_t iindex, count = 0;
	mapkit_size_t size = spm->size;
	map_strstr_storage *pos_contents;
	char *defaultvalue = spm->defaultvalue;

	pos_contents = spm->contents;
	for (iindex = 0; iindex < size; iindex++, pos_contents++)
		if ((pos_contents->state == MAPKIT_FULLSLOT)
		    && (strcmp(pos_contents->value, defaultvalue) == 0)) {
			pos_contents->state = MAPKIT_DELETEDSLOT;
			count++;
		}

	spm->used -= count;
	if (spm->used < spm->minused) {
#ifdef MAPKIT_DEBUG
		fprintf(stderr, "MAPKIT: used < minused\n");
#endif
		return map_strstr_reallocate(spm, map_strstr_meansize(spm, spm->used));
	}

	return MAPKIT_OK;
}

int map_strstr_compare(const void *e1, const void *e2)
{
	char *key1 = ((const map_strstr_element *) e1)->key;
	char *key2 = ((const map_strstr_element *) e2)->key;

	return strcmp(key1, key2);
}

mapkit_error map_strstr_getall_sorted(map_strstr * spm, map_strstr_element ** array, mapkit_size_t * count)
{
	mapkit_error err;

	err = map_strstr_getall(spm, array, count);
	if (err)
		MAPKIT_ERROR(err);

	qsort(*array, *count, sizeof(**array), map_strstr_compare);

	return MAPKIT_OK;
}

mapkit_error map_strstr_setall(map_strstr * spm, map_strstr_element * array, mapkit_size_t count)
{
	mapkit_size_t array_iindex;
	mapkit_error err;

	err = map_strstr_ensurecapacity(spm, spm->used + count);
	if (err)
		MAPKIT_ERROR(err);

	if (spm->alwaysdefault)
		/*
		 * Prevent shrinking 
		 */
		spm->minused = 0;

	for (array_iindex = 0; array_iindex < count; array_iindex++) {
		err = map_strstr_set(spm, array[array_iindex].key, array[array_iindex].value);
		if (err)
			MAPKIT_ERROR(err);
	}

	if (spm->alwaysdefault)
		map_strstr_adjustcapacity(spm);

	return MAPKIT_OK;
}

mapkit_error map_strstr_removeall(map_strstr * spm, char **array, mapkit_size_t count)
{
	mapkit_size_t array_iindex;
	mapkit_error err;

	/*
	 * Prevent shrinking 
	 */
	spm->minused = 0;

	for (array_iindex = 0; array_iindex < count; array_iindex++) {
		err = map_strstr_remove(spm, array[array_iindex]);
		if (err)
			MAPKIT_ERROR(err);
	}

	map_strstr_adjustcapacity(spm);

	return MAPKIT_OK;
}

void map_strstr_printstats(map_strstr * spm)
{
	fprintf(stderr, "MAPKIT: map_strstr statistics\n");
	fprintf(stderr, "MAPKIT: alwaysdefault = %d\n", spm->alwaysdefault);
	fprintf(stderr, "MAPKIT: minused = %ld, maxfill = %ld\n", (long) spm->minused, (long) spm->maxfill);
	fprintf(stderr, "MAPKIT: minusedfactor = %g, maxfillfactor = %g\n", spm->minusedfactor, spm->maxfillfactor);
#ifdef MAPKIT_COLLISIONS
	fprintf(stderr, "MAPKIT: insertionindexs = %lu, collisions = %lu\n", (unsigned long) spm->insertionindexs,
		(unsigned long) spm->insertionindex_collisions);
	fprintf(stderr, "MAPKIT: keyindexs = %lu, collisions = %lu\n", (unsigned long) spm->keyindexs, (unsigned long) spm->keyindex_collisions);
#endif
}

#endif							       /* MAPKIT_map_strstr */

/* Sparse structures */

#ifdef MAPKIT_spvector

/*
  Implementation for spvector (int -> double)
  Default value : 0.0
  Unassigned values always default to the default value.
*/

/* Static prototypes */

/* Return the iindex of key, or MAPKIT_KEYNOTFOUND */
static mapkit_size_t spvector_keyindex(spvector * spm, int key);

/* Return the iindex of key or -(insertion iindex)-1 if key not found */
static mapkit_size_t spvector_insertionindex(spvector * spm, int key);

/* Implementation */

mapkit_error spvector_init(spvector * spm)
{
	return spvector_init_hint(spm, MAPKIT_DEFAULT_EXPECTEDUSED);
}

mapkit_error spvector_init_hint(spvector * spm, mapkit_size_t used)
{
#ifdef MAPKIT_DEBUG
	fprintf(stderr, "MAPKIT: init\n");
#endif

	spm->size = 0;
	spm->fill = 0;
	spm->used = 0;
	spm->maxfillfactor = 0.5;
	spm->minusedfactor = 0.2;
	spm->contents = NULL;
#ifdef MAPKIT_COLLISIONS
	spm->insertionindexs = spm->insertionindex_collisions = 0;
	spm->keyindexs = spm->keyindex_collisions = 0;
#endif
	spm->defaultvalue = 0.0;
	spm->alwaysdefault = 1;

	return spvector_reallocate(spm, spvector_meansize(spm, used));
}

mapkit_error spvector_ensurecapacity(spvector * spm, mapkit_size_t used)
{
	if (used > (spm->used + spm->maxfill - spm->fill)) {
#ifdef MAPKIT_DEBUG
		fprintf(stderr, "MAPKIT: need more capacity\n");
#endif
		return spvector_reallocate(spm, spvector_meansize(spm, used));
	} else
		return MAPKIT_OK;
}

mapkit_error spvector_adjustcapacity(spvector * spm)
{
	spm->minused = (mapkit_size_t) (spm->size * spm->minusedfactor);
	spm->maxfill = (mapkit_size_t) (spm->size * spm->maxfillfactor);

	if (spm->used < spm->minused) {
#ifdef MAPKIT_DEBUG
		fprintf(stderr, "MAPKIT: used < minused\n");
#endif
		return spvector_reallocate(spm, spvector_meansize(spm, spm->used));
	} else if (spm->fill > spm->maxfill) {
#ifdef MAPKIT_DEBUG
		fprintf(stderr, "MAPKIT: fill > maxfill\n");
#endif
		return spvector_reallocate(spm, spvector_meansize(spm, spm->used));
	} else
		return MAPKIT_OK;
}

void spvector_free(spvector * spm)
{
	free(spm->contents);
	spm->contents = NULL;
	spm->size = 0;
	spm->fill = 0;
	spm->used = 0;
	spm->maxfill = 0;
	spm->minused = 0;
#ifdef MAPKIT_COLLISIONS
	spm->insertionindexs = spm->insertionindex_collisions = 0;
	spm->keyindexs = spm->keyindex_collisions = 0;
#endif
}

mapkit_error spvector_copy(spvector * to, spvector * from)
{
	spvector_storage *contentscopy;
	size_t size = from->size * sizeof(*from->contents);

	contentscopy = (spvector_storage *) malloc(size);
	if (contentscopy == NULL)
		MAPKIT_ERROR(MAPKIT_ENOMEM);

	memcpy(to, from, sizeof(spvector));
	to->contents = contentscopy;
	memcpy(to->contents, from->contents, size);

	return MAPKIT_OK;
}

mapkit_size_t spvector_growsize(spvector * spm, mapkit_size_t used)
{
	return (mapkit_size_t) (4.0 * used / (3 * spm->minusedfactor + spm->maxfillfactor));
}

mapkit_size_t spvector_shrinksize(spvector * spm, mapkit_size_t used)
{
	return (mapkit_size_t) (4.0 * used / (spm->minusedfactor + 3 * spm->maxfillfactor));
}

mapkit_size_t spvector_meansize(spvector * spm, mapkit_size_t used)
{
	return (mapkit_size_t) (2.0 * used / (spm->minusedfactor + spm->maxfillfactor));
}

mapkit_error spvector_reallocate(spvector * spm, mapkit_size_t newsize)
{
	mapkit_size_t iindex;
	mapkit_size_t oldsize;
	spvector_storage *newcontents, *oldcontents;

	/*
	 * At least one free entry 
	 */
	if (newsize <= spm->used)
		newsize = spm->used + 1;
	newsize = mapkit_nextprime(newsize);
	if (newsize <= spm->used)
		MAPKIT_ERROR(MAPKIT_ETOOBIG);

	newcontents = (spvector_storage *) malloc(newsize * sizeof(*spm->contents));
	if (newcontents == NULL)
		MAPKIT_ERROR(MAPKIT_ENOMEM);

	/*
	 * Initialize all entries to "free" 
	 */
	for (iindex = 0; iindex < newsize; iindex++)
		newcontents[iindex].key = MAPKIT_FREESLOT;

	oldcontents = spm->contents;
	oldsize = spm->size;
	spm->contents = newcontents;
	spm->size = newsize;

#ifdef MAPKIT_DEBUG
	fprintf(stderr, "MAPKIT: reallocate %ld -> %ld\n", (long) oldsize, (long) newsize);
#endif

	spm->maxfill = (mapkit_size_t) (newsize * spm->maxfillfactor);
	/*
	 * At least one free entry 
	 */
	if (spm->maxfill >= newsize)
		spm->maxfill = newsize - 1;
	spm->minused = (mapkit_size_t) (newsize * spm->minusedfactor);
	spm->used = 0;

	if (oldcontents != NULL) {
		int used = 0;
		int key;
		double defaultvalue = spm->defaultvalue;
		int notalwaysdefault = !spm->alwaysdefault;

		/*
		 * Copy all entries from old to new 
		 */
		for (iindex = 0; iindex < oldsize; iindex++)
			if ((key = oldcontents[iindex].key) >= MAPKIT_FULLSLOT) {
				mapkit_size_t ins_iindex;
				spvector_storage *contents;

				/*
				 * Fast path 
				 */
				ins_iindex = ((mapkit_hash_t) key) % spm->size;
				contents = &(newcontents[ins_iindex]);

				if (contents->key != MAPKIT_FREESLOT) {
					ins_iindex = spvector_insertionindex(spm, key);
					contents = &(newcontents[ins_iindex]);
				}
#ifdef MAPKIT_COLLISIONS
				else
					spm->insertionindexs++;
#endif
				if (notalwaysdefault || (!((oldcontents[iindex].value) == (defaultvalue)))) {
					contents->value = oldcontents[iindex].value;
					contents->key = key;
					used++;
				}
			}
		free(oldcontents);
		spm->used = used;
	}
	spm->fill = spm->used;

	return MAPKIT_OK;
}

double spvector_value_s(spvector * spm, int key)
{
	mapkit_size_t iindex;

	iindex = spvector_keyindex(spm, key);

	if (iindex < 0) {
		if (spm->alwaysdefault)
			return spm->defaultvalue;
		else
			MAPKIT_FATAL_ERROR(MAPKIT_EKEYNOTFOUND);
	}

	return spm->contents[iindex].value;
}

mapkit_error spvector_get_s(spvector * spm, int key, double *value)
{
	mapkit_size_t iindex;

	iindex = spvector_keyindex(spm, key);

	if (iindex < 0) {
		if (spm->alwaysdefault) {
			*value = spm->defaultvalue;
			return MAPKIT_OK;
		} else
			MAPKIT_ERROR(MAPKIT_EKEYNOTFOUND);
	}

	*value = spm->contents[iindex].value;
	return MAPKIT_OK;
}

mapkit_error spvector_set_s(spvector * spm, int key, double value)
{
	mapkit_size_t iindex;

	iindex = spvector_insertionindex(spm, key);

	if (iindex < 0)
		/*
		 * FULLSLOT 
		 */
		spm->contents[-iindex - 1].value = value;
	else {
		spvector_storage *element = &(spm->contents[iindex]);
		int ffree = element->key == MAPKIT_FREESLOT;

		element->value = value;
		element->key = key;
		spm->used++;

		if (ffree && ((++spm->fill) > spm->maxfill)) {
#ifdef MAPKIT_DEBUG
			fprintf(stderr, "MAPKIT: fill > maxfill\n");
#endif
			return spvector_reallocate(spm, spvector_growsize(spm, spm->used));
		}
	}
	return MAPKIT_OK;
}

double *spvector_insertptr_s(spvector * spm, int key)
{
	mapkit_size_t iindex;

	iindex = spvector_insertionindex(spm, key);

	if (iindex < 0)
		return &(spm->contents[-iindex - 1].value);
	else {
		spvector_storage *element = &(spm->contents[iindex]);

		if (element->key == MAPKIT_FREESLOT) {
			/*
			 * FREESLOT 
			 */
			if (spm->fill >= spm->maxfill) {
				mapkit_error err;

#ifdef MAPKIT_DEBUG
				fprintf(stderr, "MAPKIT: fill => maxfill before insert\n");
#endif
				/*
				 * Must reallocate -before- inserting defaultvalue 
				 */
				err = spvector_reallocate(spm, spvector_growsize(spm, spm->used + 1));
				if (err) {
					MAPKIT_ERROR_NORET(err);
					return NULL;
				}

				iindex = spvector_insertionindex(spm, key);
				/*
				 * FREESLOT 
				 */
				spm->contents[iindex].key = key;
				spm->contents[iindex].value = spm->defaultvalue;
				spm->used++;
				spm->fill++;
				return &(spm->contents[iindex].value);
			} else
				spm->fill++;
		}

		element->key = key;
		element->value = spm->defaultvalue;
		spm->used++;
		return &(element->value);
	}
}

double *spvector_ptr_s(spvector * spm, int key)
{
	mapkit_size_t iindex;

	iindex = spvector_keyindex(spm, key);

	if ((iindex >= 0) && ((!spm->alwaysdefault) || (!((spm->contents[iindex].value) == (spm->defaultvalue)))))
		return &(spm->contents[iindex].value);
	else
		return NULL;
}

mapkit_error spvector_remove_s(spvector * spm, int key)
{
	mapkit_size_t iindex;

	iindex = spvector_keyindex(spm, key);

	if (iindex < 0) {
		if (spm->alwaysdefault)
			return MAPKIT_OK;
		else
			return MAPKIT_EKEYNOTFOUND;
	}

	spm->contents[iindex].key = MAPKIT_DELETEDSLOT;
	spm->used--;
	if (spm->used < spm->minused) {
#ifdef MAPKIT_DEBUG
		fprintf(stderr, "MAPKIT: used < minused\n");
#endif
		return spvector_reallocate(spm, spvector_shrinksize(spm, spm->used));
	}

	return MAPKIT_OK;
}

mapkit_size_t spvector_keyindex(spvector * spm, int key)
{
	mapkit_size_t iindex, decrement;

	int ckey;

#ifdef MAPKIT_DEBUG
	if (key < MAPKIT_FULLSLOT) {
		/*
		 * Not user-called : should never happen 
		 */
		MAPKIT_ERROR_NORET(MAPKIT_EBADKEY);
		return MAPKIT_KEYNOTFOUND;
	}
#endif

	iindex = ((mapkit_hash_t) key) % spm->size;
	decrement = (((mapkit_hash_t) key) % (spm->size - 2));
	decrement += (decrement == 0);

#ifdef MAPKIT_COLLISIONS
	spm->keyindexs++;
#endif

	while ((ckey = spm->contents[iindex].key) != MAPKIT_FREESLOT && (ckey == MAPKIT_DELETEDSLOT || ckey != key)) {
#ifdef MAPKIT_COLLISIONS
		spm->keyindex_collisions++;
#endif
		iindex -= decrement;
		if (iindex < 0)
			iindex += spm->size;
	}

	if (ckey == MAPKIT_FREESLOT)
		return MAPKIT_KEYNOTFOUND;
	return iindex;
}

mapkit_size_t spvector_insertionindex(spvector * spm, int key)
{
	mapkit_size_t iindex, decrement;
	int ckey;

#ifdef MAPKIT_DEBUG
	if (key < MAPKIT_FULLSLOT)
		/*
		 * Not user-called : should never happen 
		 */
		MAPKIT_FATAL_ERROR(MAPKIT_EBADKEY);
#endif

#ifdef MAPKIT_COLLISIONS
	spm->insertionindexs++;
#endif

	iindex = ((mapkit_hash_t) key) % spm->size;

	/*
	 * Fast path (largely superfluous) 
	 */
	if ((ckey = spm->contents[iindex].key) == MAPKIT_FREESLOT)
		return iindex;
	if (ckey == key)
		return -iindex - 1;

	decrement = (((mapkit_hash_t) key) % (spm->size - 2));
	decrement += (decrement == 0);

	while ((ckey >= 0) && (ckey != key)) {
#ifdef MAPKIT_COLLISIONS
		spm->insertionindex_collisions++;
#endif
		iindex -= decrement;
		if (iindex < 0)
			iindex += spm->size;
		ckey = spm->contents[iindex].key;
	}

	if (ckey == MAPKIT_DELETEDSLOT) {
		mapkit_size_t index2 = iindex;

		while ((ckey != MAPKIT_FREESLOT)
		       && (ckey == MAPKIT_DELETEDSLOT || ckey != key)) {
			iindex -= decrement;
			if (iindex < 0)
				iindex += spm->size;
			ckey = spm->contents[iindex].key;
		}
		if (ckey == MAPKIT_FREESLOT)
			return index2;
	}

	if (ckey >= 0)
		return -iindex - 1;
	return iindex;
}

mapkit_size_t spvector_next(spvector * spm, mapkit_size_t iindex)
{
	mapkit_size_t size = spm->size;
	spvector_storage *pos_contents;
	double defaultvalue = spm->defaultvalue;
	int notalwaysdefault = !spm->alwaysdefault;

	pos_contents = &(spm->contents[++iindex]);

	for (; iindex < size; iindex++, pos_contents++)
		if ((pos_contents->key >= MAPKIT_FULLSLOT)
		    && (notalwaysdefault || (!((pos_contents->value) == (defaultvalue)))))
			return iindex;

	return -1;
}

spvector_storage *spvector_nextptr(spvector * spm, spvector_storage * pos_contents)
{
	spvector_storage *end = &(spm->contents[spm->size]);
	double defaultvalue = spm->defaultvalue;
	int notalwaysdefault = !spm->alwaysdefault;

	if (pos_contents == NULL)
		pos_contents = spm->contents;
	else {
		pos_contents++;
		if (pos_contents <= spm->contents)
			return NULL;
	}

	for (; pos_contents < end; pos_contents++)
		if ((pos_contents->key >= MAPKIT_FULLSLOT)
		    && (notalwaysdefault || (!((pos_contents->value) == (defaultvalue)))))
			return pos_contents;

	return NULL;
}

mapkit_error spvector_getall(spvector * spm, spvector_element ** array, mapkit_size_t * count)
{
	mapkit_size_t iindex;
	mapkit_size_t size = spm->size, vcount = 0;
	int key;
	spvector_element *pos_array;
	spvector_storage *pos_contents;
	double defaultvalue = spm->defaultvalue;
	int notalwaysdefault = !spm->alwaysdefault;

	pos_array = *array = (spvector_element *) malloc(sizeof(**array) * spm->used);
	if (*array == NULL)
		MAPKIT_ERROR(MAPKIT_ENOMEM);

	pos_contents = spm->contents;

	for (iindex = 0; iindex < size; iindex++, pos_contents++)
		if (((key = pos_contents->key) >= MAPKIT_FULLSLOT)
		    && (notalwaysdefault || (!((pos_contents->value) == (defaultvalue))))) {
			pos_array->key = key;
			pos_array->value = pos_contents->value;
			pos_array++;
			vcount++;
		}
	*count = vcount;

	return MAPKIT_OK;
}

mapkit_error spvector_clean(spvector * spm)
{
	mapkit_size_t iindex, count = 0;
	mapkit_size_t size = spm->size;
	spvector_storage *pos_contents;
	double defaultvalue = spm->defaultvalue;

	pos_contents = spm->contents;
	for (iindex = 0; iindex < size; iindex++, pos_contents++)
		if ((pos_contents->key >= MAPKIT_FULLSLOT)
		    && ((pos_contents->value) == (defaultvalue))) {
			pos_contents->key = MAPKIT_DELETEDSLOT;
			count++;
		}

	spm->used -= count;
	if (spm->used < spm->minused) {
#ifdef MAPKIT_DEBUG
		fprintf(stderr, "MAPKIT: used < minused\n");
#endif
		return spvector_reallocate(spm, spvector_meansize(spm, spm->used));
	}

	return MAPKIT_OK;
}

int spvector_compare(const void *e1, const void *e2)
{
	int key1 = ((const spvector_element *) e1)->key;
	int key2 = ((const spvector_element *) e2)->key;

	return ((key1) < (key2) ? -1 : ((key1) == (key2) ? 0 : 1));
}

mapkit_error spvector_getall_sorted(spvector * spm, spvector_element ** array, mapkit_size_t * count)
{
	mapkit_error err;

	err = spvector_getall(spm, array, count);
	if (err)
		MAPKIT_ERROR(err);

	qsort(*array, *count, sizeof(**array), spvector_compare);

	return MAPKIT_OK;
}

mapkit_error spvector_setall(spvector * spm, spvector_element * array, mapkit_size_t count)
{
	mapkit_size_t array_iindex;
	mapkit_error err;

	err = spvector_ensurecapacity(spm, spm->used + count);
	if (err)
		MAPKIT_ERROR(err);

	if (spm->alwaysdefault)
		/*
		 * Prevent shrinking 
		 */
		spm->minused = 0;

	for (array_iindex = 0; array_iindex < count; array_iindex++) {
		err = spvector_set(spm, array[array_iindex].key, array[array_iindex].value);
		if (err)
			MAPKIT_ERROR(err);
	}

	if (spm->alwaysdefault)
		spvector_adjustcapacity(spm);

	return MAPKIT_OK;
}

mapkit_error spvector_removeall(spvector * spm, int *array, mapkit_size_t count)
{
	mapkit_size_t array_iindex;
	mapkit_error err;

	/*
	 * Prevent shrinking 
	 */
	spm->minused = 0;

	for (array_iindex = 0; array_iindex < count; array_iindex++) {
		err = spvector_remove(spm, array[array_iindex]);
		if (err)
			MAPKIT_ERROR(err);
	}

	spvector_adjustcapacity(spm);

	return MAPKIT_OK;
}

void spvector_printstats(spvector * spm)
{
	fprintf(stderr, "MAPKIT: spvector statistics\n");
	fprintf(stderr, "MAPKIT: alwaysdefault = %d\n", spm->alwaysdefault);
	fprintf(stderr, "MAPKIT: minused = %ld, maxfill = %ld\n", (long) spm->minused, (long) spm->maxfill);
	fprintf(stderr, "MAPKIT: minusedfactor = %g, maxfillfactor = %g\n", spm->minusedfactor, spm->maxfillfactor);
#ifdef MAPKIT_COLLISIONS
	fprintf(stderr, "MAPKIT: insertionindexs = %lu, collisions = %lu\n", (unsigned long) spm->insertionindexs,
		(unsigned long) spm->insertionindex_collisions);
	fprintf(stderr, "MAPKIT: keyindexs = %lu, collisions = %lu\n", (unsigned long) spm->keyindexs, (unsigned long) spm->keyindex_collisions);
#endif
}

#endif							       /* MAPKIT_spvector */

#ifdef MAPKIT__spmatrix

/*
  Implementation for _spmatrix (spmatrix_key_pair -> double)
  Default value : 0.0
  Uses a state field.
  Unassigned values always default to the default value.
*/

/* Static prototypes */

/* Return the iindex of key, or MAPKIT_KEYNOTFOUND */
static mapkit_size_t _spmatrix_keyindex(_spmatrix * spm, spmatrix_key_pair key, mapkit_hash_t hash);

/* Return the iindex of key or -(insertion iindex)-1 if key not found */
static mapkit_size_t _spmatrix_insertionindex(_spmatrix * spm, spmatrix_key_pair key, mapkit_hash_t hash);

/* Implementation */

mapkit_error _spmatrix_init(_spmatrix * spm)
{
	return _spmatrix_init_hint(spm, MAPKIT_DEFAULT_EXPECTEDUSED);
}

mapkit_error _spmatrix_init_hint(_spmatrix * spm, mapkit_size_t used)
{
#ifdef MAPKIT_DEBUG
	fprintf(stderr, "MAPKIT: init\n");
#endif

	spm->size = 0;
	spm->fill = 0;
	spm->used = 0;
	spm->maxfillfactor = 0.5;
	spm->minusedfactor = 0.2;
	spm->contents = NULL;
#ifdef MAPKIT_COLLISIONS
	spm->insertionindexs = spm->insertionindex_collisions = 0;
	spm->keyindexs = spm->keyindex_collisions = 0;
#endif
	spm->defaultvalue = 0.0;
	spm->alwaysdefault = 1;

	return _spmatrix_reallocate(spm, _spmatrix_meansize(spm, used));
}

mapkit_error _spmatrix_ensurecapacity(_spmatrix * spm, mapkit_size_t used)
{
	if (used > (spm->used + spm->maxfill - spm->fill)) {
#ifdef MAPKIT_DEBUG
		fprintf(stderr, "MAPKIT: need more capacity\n");
#endif
		return _spmatrix_reallocate(spm, _spmatrix_meansize(spm, used));
	} else
		return MAPKIT_OK;
}

mapkit_error _spmatrix_adjustcapacity(_spmatrix * spm)
{
	spm->minused = (mapkit_size_t) (spm->size * spm->minusedfactor);
	spm->maxfill = (mapkit_size_t) (spm->size * spm->maxfillfactor);

	if (spm->used < spm->minused) {
#ifdef MAPKIT_DEBUG
		fprintf(stderr, "MAPKIT: used < minused\n");
#endif
		return _spmatrix_reallocate(spm, _spmatrix_meansize(spm, spm->used));
	} else if (spm->fill > spm->maxfill) {
#ifdef MAPKIT_DEBUG
		fprintf(stderr, "MAPKIT: fill > maxfill\n");
#endif
		return _spmatrix_reallocate(spm, _spmatrix_meansize(spm, spm->used));
	} else
		return MAPKIT_OK;
}

void _spmatrix_free(_spmatrix * spm)
{
	free(spm->contents);
	spm->contents = NULL;
	spm->size = 0;
	spm->fill = 0;
	spm->used = 0;
	spm->maxfill = 0;
	spm->minused = 0;
#ifdef MAPKIT_COLLISIONS
	spm->insertionindexs = spm->insertionindex_collisions = 0;
	spm->keyindexs = spm->keyindex_collisions = 0;
#endif
}

mapkit_error _spmatrix_copy(_spmatrix * to, _spmatrix * from)
{
	_spmatrix_storage *contentscopy;
	size_t size = from->size * sizeof(*from->contents);

	contentscopy = (_spmatrix_storage *) malloc(size);
	if (contentscopy == NULL)
		MAPKIT_ERROR(MAPKIT_ENOMEM);

	memcpy(to, from, sizeof(_spmatrix));
	to->contents = contentscopy;
	memcpy(to->contents, from->contents, size);

	return MAPKIT_OK;
}

mapkit_size_t _spmatrix_growsize(_spmatrix * spm, mapkit_size_t used)
{
	return (mapkit_size_t) (4.0 * used / (3 * spm->minusedfactor + spm->maxfillfactor));
}

mapkit_size_t _spmatrix_shrinksize(_spmatrix * spm, mapkit_size_t used)
{
	return (mapkit_size_t) (4.0 * used / (spm->minusedfactor + 3 * spm->maxfillfactor));
}

mapkit_size_t _spmatrix_meansize(_spmatrix * spm, mapkit_size_t used)
{
	return (mapkit_size_t) (2.0 * used / (spm->minusedfactor + spm->maxfillfactor));
}

mapkit_error _spmatrix_reallocate(_spmatrix * spm, mapkit_size_t newsize)
{
	mapkit_size_t iindex;
	mapkit_size_t oldsize;
	_spmatrix_storage *newcontents, *oldcontents;

	/*
	 * At least one free entry 
	 */
	if (newsize <= spm->used)
		newsize = spm->used + 1;
	newsize = mapkit_nextprime(newsize);
	if (newsize <= spm->used)
		MAPKIT_ERROR(MAPKIT_ETOOBIG);

	newcontents = (_spmatrix_storage *) malloc(newsize * sizeof(*spm->contents));
	if (newcontents == NULL)
		MAPKIT_ERROR(MAPKIT_ENOMEM);

	/*
	 * Initialize all entries to "free" 
	 */
	for (iindex = 0; iindex < newsize; iindex++)
		newcontents[iindex].state = MAPKIT_FREESLOT;

	oldcontents = spm->contents;
	oldsize = spm->size;
	spm->contents = newcontents;
	spm->size = newsize;

#ifdef MAPKIT_DEBUG
	fprintf(stderr, "MAPKIT: reallocate %ld -> %ld\n", (long) oldsize, (long) newsize);
#endif

	spm->maxfill = (mapkit_size_t) (newsize * spm->maxfillfactor);
	/*
	 * At least one free entry 
	 */
	if (spm->maxfill >= newsize)
		spm->maxfill = newsize - 1;
	spm->minused = (mapkit_size_t) (newsize * spm->minusedfactor);
	spm->used = 0;

	if (oldcontents != NULL) {
		int used = 0;
		spmatrix_key_pair key;
		double defaultvalue = spm->defaultvalue;
		int notalwaysdefault = !spm->alwaysdefault;

		/*
		 * Copy all entries from old to new 
		 */
		for (iindex = 0; iindex < oldsize; iindex++)
			if (oldcontents[iindex].state == MAPKIT_FULLSLOT) {
				mapkit_size_t ins_iindex;
				mapkit_hash_t hash;
				_spmatrix_storage *contents;

				key = oldcontents[iindex].key;

				/*
				 * Fast path 
				 */
				ins_iindex = (hash = (mapkit_hash_t) ((key).key1 + (key).key2 * 65519)) % spm->size;
				contents = &(newcontents[ins_iindex]);

				if (contents->state != MAPKIT_FREESLOT) {
					ins_iindex = _spmatrix_insertionindex(spm, key, hash);
					contents = &(newcontents[ins_iindex]);
				}
#ifdef MAPKIT_COLLISIONS
				else
					spm->insertionindexs++;
#endif
				if (notalwaysdefault || (!((oldcontents[iindex].value) == (defaultvalue)))) {
					contents->value = oldcontents[iindex].value;
					contents->state = MAPKIT_FULLSLOT;
					contents->key = key;
					used++;
				}
			}
		free(oldcontents);
		spm->used = used;
	}
	spm->fill = spm->used;

	return MAPKIT_OK;
}

double _spmatrix_value_s(_spmatrix * spm, spmatrix_key_pair key, mapkit_hash_t hash)
{
	mapkit_size_t iindex;

	iindex = _spmatrix_keyindex(spm, key, hash);

	if (iindex < 0) {
		if (spm->alwaysdefault)
			return spm->defaultvalue;
		else
			MAPKIT_FATAL_ERROR(MAPKIT_EKEYNOTFOUND);
	}

	return spm->contents[iindex].value;
}

mapkit_error _spmatrix_get_s(_spmatrix * spm, spmatrix_key_pair key, double *value, mapkit_hash_t hash)
{
	mapkit_size_t iindex;

	iindex = _spmatrix_keyindex(spm, key, hash);

	if (iindex < 0) {
		if (spm->alwaysdefault) {
			*value = spm->defaultvalue;
			return MAPKIT_OK;
		} else {
			MAPKIT_ERROR(MAPKIT_EKEYNOTFOUND);
		}
	}
	*value = spm->contents[iindex].value;
	return MAPKIT_OK;
}

mapkit_error _spmatrix_set_s(_spmatrix * spm, spmatrix_key_pair key, double value, mapkit_hash_t hash)
{
	mapkit_size_t iindex;

	iindex = _spmatrix_insertionindex(spm, key, hash);

	if (iindex < 0)
		/*
		 * FULLSLOT 
		 */
		spm->contents[-iindex - 1].value = value;
	else {
		_spmatrix_storage *element = &(spm->contents[iindex]);
		int ffree = element->state == MAPKIT_FREESLOT;

		element->state = MAPKIT_FULLSLOT;
		element->value = value;
		element->key = key;
		spm->used++;

		if (ffree && ((++spm->fill) > spm->maxfill)) {
#ifdef MAPKIT_DEBUG
			fprintf(stderr, "MAPKIT: fill > maxfill\n");
#endif
			return _spmatrix_reallocate(spm, _spmatrix_growsize(spm, spm->used));
		}
	}
	return MAPKIT_OK;
}

double *_spmatrix_insertptr_s(_spmatrix * spm, spmatrix_key_pair key, mapkit_hash_t hash)
{
	mapkit_size_t iindex;

	iindex = _spmatrix_insertionindex(spm, key, hash);

	if (iindex < 0)
		return &(spm->contents[-iindex - 1].value);
	else {
		_spmatrix_storage *element = &(spm->contents[iindex]);

		if (element->state == MAPKIT_FREESLOT) {
			/*
			 * FREESLOT 
			 */
			if (spm->fill >= spm->maxfill) {
				mapkit_error err;

#ifdef MAPKIT_DEBUG
				fprintf(stderr, "MAPKIT: fill => maxfill before insert\n");
#endif
				/*
				 * Must reallocate -before- inserting defaultvalue 
				 */
				err = _spmatrix_reallocate(spm, _spmatrix_growsize(spm, spm->used + 1));
				if (err) {
					MAPKIT_ERROR_NORET(err);
					return NULL;
				}

				iindex = _spmatrix_insertionindex(spm, key, hash);
				/*
				 * FREESLOT 
				 */
				spm->contents[iindex].state = MAPKIT_FULLSLOT;
				spm->contents[iindex].key = key;
				spm->contents[iindex].value = spm->defaultvalue;
				spm->used++;
				spm->fill++;
				return &(spm->contents[iindex].value);
			} else
				spm->fill++;
		}

		element->state = MAPKIT_FULLSLOT;
		element->key = key;
		element->value = spm->defaultvalue;
		spm->used++;
		return &(element->value);
	}
}

double *_spmatrix_ptr_s(_spmatrix * spm, spmatrix_key_pair key, mapkit_hash_t hash)
{
	mapkit_size_t iindex;

	iindex = _spmatrix_keyindex(spm, key, hash);

	if ((iindex >= 0) && ((!spm->alwaysdefault) || (!((spm->contents[iindex].value) == (spm->defaultvalue)))))
		return &(spm->contents[iindex].value);
	else
		return NULL;
}

mapkit_error _spmatrix_remove_s(_spmatrix * spm, spmatrix_key_pair key, mapkit_hash_t hash)
{
	mapkit_size_t iindex;

	iindex = _spmatrix_keyindex(spm, key, hash);

	if (iindex < 0) {
		if (spm->alwaysdefault)
			return MAPKIT_OK;
		else
			return MAPKIT_EKEYNOTFOUND;
	}

	spm->contents[iindex].state = MAPKIT_DELETEDSLOT;
	spm->used--;
	if (spm->used < spm->minused) {
#ifdef MAPKIT_DEBUG
		fprintf(stderr, "MAPKIT: used < minused\n");
#endif
		return _spmatrix_reallocate(spm, _spmatrix_shrinksize(spm, spm->used));
	}

	return MAPKIT_OK;
}

mapkit_size_t _spmatrix_keyindex(_spmatrix * spm, spmatrix_key_pair key, mapkit_hash_t hash)
{
	mapkit_size_t iindex, decrement;

	signed char state;

	iindex = hash % spm->size;
	decrement = (hash % (spm->size - 2));
	decrement += (decrement == 0);

#ifdef MAPKIT_COLLISIONS
	spm->keyindexs++;
#endif

	while ((state = spm->contents[iindex].state) != MAPKIT_FREESLOT
	       && (state == MAPKIT_DELETEDSLOT
		   || (!(((spm->contents[iindex].key).key1 == (key).key1) && ((spm->contents[iindex].key).key2 == (key).key2))))) {
#ifdef MAPKIT_COLLISIONS
		spm->keyindex_collisions++;
#endif
		iindex -= decrement;
		if (iindex < 0)
			iindex += spm->size;
	}

	if (state == MAPKIT_FREESLOT)
		return MAPKIT_KEYNOTFOUND;
	return iindex;
}

mapkit_size_t _spmatrix_insertionindex(_spmatrix * spm, spmatrix_key_pair key, mapkit_hash_t hash)
{
	mapkit_size_t iindex, decrement;
	signed char state;

#ifdef MAPKIT_COLLISIONS
	spm->insertionindexs++;
#endif

	iindex = hash % spm->size;

	/*
	 * Fast path (largely superfluous) 
	 */
	if ((state = spm->contents[iindex].state) == MAPKIT_FREESLOT)
		return iindex;
	if ((state == MAPKIT_FULLSLOT)
	    && (((spm->contents[iindex].key).key1 == (key).key1) && ((spm->contents[iindex].key).key2 == (key).key2)))
		return -iindex - 1;

	decrement = (hash % (spm->size - 2));
	decrement += (decrement == 0);

	while ((state == MAPKIT_FULLSLOT)
	       && (!(((spm->contents[iindex].key).key1 == (key).key1) && ((spm->contents[iindex].key).key2 == (key).key2)))) {
#ifdef MAPKIT_COLLISIONS
		spm->insertionindex_collisions++;
#endif
		iindex -= decrement;
		if (iindex < 0)
			iindex += spm->size;
		state = spm->contents[iindex].state;
	}

	if (state == MAPKIT_DELETEDSLOT) {
		mapkit_size_t index2 = iindex;

		while ((state = spm->contents[iindex].state) != MAPKIT_FREESLOT && ((state == MAPKIT_DELETEDSLOT)
										    || (!(((spm->contents[iindex].key).key1 == (key).key1)
											  && ((spm->contents[iindex].key).key2 == (key).key2))))) {
			iindex -= decrement;
			if (iindex < 0)
				iindex += spm->size;
			state = spm->contents[iindex].state;
		}
		if (state == MAPKIT_FREESLOT)
			return index2;
	}

	if (state == MAPKIT_FULLSLOT)
		return -iindex - 1;
	return iindex;
}

mapkit_size_t _spmatrix_next(_spmatrix * spm, mapkit_size_t iindex)
{
	mapkit_size_t size = spm->size;
	_spmatrix_storage *pos_contents;
	double defaultvalue = spm->defaultvalue;
	int notalwaysdefault = !spm->alwaysdefault;

	pos_contents = &(spm->contents[++iindex]);

	for (; iindex < size; iindex++, pos_contents++)
		if ((pos_contents->state == MAPKIT_FULLSLOT)
		    && (notalwaysdefault || (!((pos_contents->value) == (defaultvalue)))))
			return iindex;

	return -1;
}

_spmatrix_storage *_spmatrix_nextptr(_spmatrix * spm, _spmatrix_storage * pos_contents)
{
	_spmatrix_storage *end = &(spm->contents[spm->size]);
	double defaultvalue = spm->defaultvalue;
	int notalwaysdefault = !spm->alwaysdefault;

	if (pos_contents == NULL)
		pos_contents = spm->contents;
	else {
		pos_contents++;
		if (pos_contents <= spm->contents)
			return NULL;
	}

	for (; pos_contents < end; pos_contents++)
		if ((pos_contents->state == MAPKIT_FULLSLOT)
		    && (notalwaysdefault || (!((pos_contents->value) == (defaultvalue)))))
			return pos_contents;

	return NULL;
}

mapkit_error _spmatrix_getall(_spmatrix * spm, _spmatrix_element ** array, mapkit_size_t * count)
{
	mapkit_size_t iindex;
	mapkit_size_t size = spm->size, vcount = 0;
	_spmatrix_element *pos_array;
	_spmatrix_storage *pos_contents;
	double defaultvalue = spm->defaultvalue;
	int notalwaysdefault = !spm->alwaysdefault;

	pos_array = *array = (_spmatrix_element *) malloc(sizeof(**array) * spm->used);
	if (*array == NULL)
		MAPKIT_ERROR(MAPKIT_ENOMEM);

	pos_contents = spm->contents;

	for (iindex = 0; iindex < size; iindex++, pos_contents++)
		if ((pos_contents->state == MAPKIT_FULLSLOT)
		    && (notalwaysdefault || (!((pos_contents->value) == (defaultvalue))))) {
			pos_array->key = pos_contents->key;
			pos_array->value = pos_contents->value;
			pos_array++;
			vcount++;
		}
	*count = vcount;

	return MAPKIT_OK;
}

mapkit_error _spmatrix_clean(_spmatrix * spm)
{
	mapkit_size_t iindex, count = 0;
	mapkit_size_t size = spm->size;
	_spmatrix_storage *pos_contents;
	double defaultvalue = spm->defaultvalue;

	pos_contents = spm->contents;
	for (iindex = 0; iindex < size; iindex++, pos_contents++)
		if ((pos_contents->state == MAPKIT_FULLSLOT)
		    && ((pos_contents->value) == (defaultvalue))) {
			pos_contents->state = MAPKIT_DELETEDSLOT;
			count++;
		}

	spm->used -= count;
	if (spm->used < spm->minused) {
#ifdef MAPKIT_DEBUG
		fprintf(stderr, "MAPKIT: used < minused\n");
#endif
		return _spmatrix_reallocate(spm, _spmatrix_meansize(spm, spm->used));
	}

	return MAPKIT_OK;
}

int _spmatrix_compare(const void *e1, const void *e2)
{
	spmatrix_key_pair key1 = ((const _spmatrix_element *) e1)->key;
	spmatrix_key_pair key2 = ((const _spmatrix_element *) e2)->key;

	return (((key1).key1 <
		 (key2).key1) ? -1 : (((key1).key1 >
				       (key2).key1) ? 1 : (((key1).key2 < (key2).key2) ? -1 : (((key1).key2 == (key2).key2) ? 0 : 1))));
}

mapkit_error _spmatrix_getall_sorted(_spmatrix * spm, _spmatrix_element ** array, mapkit_size_t * count)
{
	mapkit_error err;

	err = _spmatrix_getall(spm, array, count);
	if (err)
		MAPKIT_ERROR(err);

	qsort(*array, *count, sizeof(**array), _spmatrix_compare);

	return MAPKIT_OK;
}

mapkit_error _spmatrix_setall(_spmatrix * spm, _spmatrix_element * array, mapkit_size_t count)
{
	mapkit_size_t array_iindex;
	mapkit_error err;

	err = _spmatrix_ensurecapacity(spm, spm->used + count);
	if (err)
		MAPKIT_ERROR(err);

	if (spm->alwaysdefault)
		/*
		 * Prevent shrinking 
		 */
		spm->minused = 0;

	for (array_iindex = 0; array_iindex < count; array_iindex++) {
		err = _spmatrix_set(spm, array[array_iindex].key, array[array_iindex].value);
		if (err)
			MAPKIT_ERROR(err);
	}

	if (spm->alwaysdefault)
		_spmatrix_adjustcapacity(spm);

	return MAPKIT_OK;
}

mapkit_error _spmatrix_removeall(_spmatrix * spm, spmatrix_key_pair * array, mapkit_size_t count)
{
	mapkit_size_t array_iindex;
	mapkit_error err;

	/*
	 * Prevent shrinking 
	 */
	spm->minused = 0;

	for (array_iindex = 0; array_iindex < count; array_iindex++) {
		err = _spmatrix_remove(spm, array[array_iindex]);
		if (err)
			MAPKIT_ERROR(err);
	}

	_spmatrix_adjustcapacity(spm);

	return MAPKIT_OK;
}

void _spmatrix_printstats(_spmatrix * spm)
{
	fprintf(stderr, "MAPKIT: _spmatrix statistics\n");
	fprintf(stderr, "MAPKIT: alwaysdefault = %d\n", spm->alwaysdefault);
	fprintf(stderr, "MAPKIT: minused = %ld, maxfill = %ld\n", (long) spm->minused, (long) spm->maxfill);
	fprintf(stderr, "MAPKIT: minusedfactor = %g, maxfillfactor = %g\n", spm->minusedfactor, spm->maxfillfactor);
#ifdef MAPKIT_COLLISIONS
	fprintf(stderr, "MAPKIT: insertionindexs = %lu, collisions = %lu\n", (unsigned long) spm->insertionindexs,
		(unsigned long) spm->insertionindex_collisions);
	fprintf(stderr, "MAPKIT: keyindexs = %lu, collisions = %lu\n", (unsigned long) spm->keyindexs, (unsigned long) spm->keyindex_collisions);
#endif
}

#endif							       /* MAPKIT__spmatrix */

#ifdef MAPKIT_spmatrix

/* Prototypes */

/* compare the key of two spmatrix_element (for qsort) key1 then key2 */
INLINE int spmatrix_compare1(const void *e1, const void *e2);

/* compare the key of two spmatrix_element (for qsort) key2 then key1 */
INLINE int spmatrix_compare2(const void *e1, const void *e2);

/* Implementation */

mapkit_error spmatrix_init(spmatrix * spm)
{
	return _spmatrix_init(spm);
}

void spmatrix_free(spmatrix * spm)
{
	_spmatrix_free(spm);
}

mapkit_error spmatrix_copy(spmatrix * to, spmatrix * from)
{
	return _spmatrix_copy(to, from);
}

mapkit_error spmatrix_init_hint(spmatrix * spm, mapkit_size_t newsize)
{
	return _spmatrix_init_hint(spm, newsize);
}

mapkit_error spmatrix_ensurecapacity(spmatrix * spm, mapkit_size_t used)
{
	return _spmatrix_ensurecapacity(spm, used);
}

mapkit_error spmatrix_adjustcapacity(spmatrix * spm)
{
	return _spmatrix_adjustcapacity(spm);
}

mapkit_error spmatrix_getall(spmatrix * spm, spmatrix_element ** array, mapkit_size_t * count)
{
	mapkit_size_t iindex;
	mapkit_size_t size = spm->size, vcount = 0;
	spmatrix_key_pair key;
	spmatrix_element *pos_array;
	spmatrix_storage *pos_contents;
	double defaultvalue = spm->defaultvalue;
	int notalwaysdefault = !spm->alwaysdefault;

	pos_array = *array = (spmatrix_element *) malloc(sizeof(**array) * spm->used);
	if (*array == NULL)
		MAPKIT_ERROR(MAPKIT_ENOMEM);

	pos_contents = spm->contents;

	for (iindex = 0; iindex < size; iindex++, pos_contents++)
		if ((pos_contents->state == MAPKIT_FULLSLOT)
		    && (notalwaysdefault || (!((pos_contents->value) == (defaultvalue))))) {
			key = spm->contents[iindex].key;
			pos_array->key1 = ((key).key1);
			pos_array->key2 = ((key).key2);
			pos_array->value = spm->contents[iindex].value;
			pos_array++;
			vcount++;
		}
	*count = vcount;

	return MAPKIT_OK;
}

mapkit_error spmatrix_clean(spmatrix * spm)
{
	return _spmatrix_clean(spm);
}

int spmatrix_compare1(const void *e1, const void *e2)
{
	int key1 = ((const spmatrix_element *) e1)->key1, key2 = ((const spmatrix_element *) e2)->key1;
	int comp = ((key1) < (key2) ? -1 : ((key1) == (key2) ? 0 : 1));

	if (comp)
		return comp;

	key1 = ((const spmatrix_element *) e1)->key2;
	key2 = ((const spmatrix_element *) e2)->key2;

	return ((key1) < (key2) ? -1 : ((key1) == (key2) ? 0 : 1));
}

int spmatrix_compare2(const void *e1, const void *e2)
{
	int key1 = ((const spmatrix_element *) e1)->key2, key2 = ((const spmatrix_element *) e2)->key2;
	int comp = ((key1) < (key2) ? -1 : ((key1) == (key2) ? 0 : 1));

	if (comp)
		return comp;

	key1 = ((const spmatrix_element *) e1)->key1;
	key2 = ((const spmatrix_element *) e2)->key1;

	return ((key1) < (key2) ? -1 : ((key1) == (key2) ? 0 : 1));
}

mapkit_error spmatrix_getall_sorted1(spmatrix * spm, spmatrix_element ** array, mapkit_size_t * count)
{
	mapkit_error err;

	err = spmatrix_getall(spm, array, count);
	if (err)
		MAPKIT_ERROR(err);

	qsort(*array, *count, sizeof(**array), spmatrix_compare1);

	return MAPKIT_OK;
}

mapkit_error spmatrix_getall_sorted2(spmatrix * spm, spmatrix_element ** array, mapkit_size_t * count)
{
	mapkit_error err;

	err = spmatrix_getall(spm, array, count);
	if (err)
		MAPKIT_ERROR(err);

	qsort(*array, *count, sizeof(**array), spmatrix_compare2);

	return MAPKIT_OK;
}

mapkit_error spmatrix_setall(spmatrix * spm, spmatrix_element * array, mapkit_size_t count)
{
	mapkit_size_t array_iindex;
	mapkit_error err;

	err = spmatrix_ensurecapacity(spm, spm->used + count);
	if (err)
		MAPKIT_ERROR(err);

	if (spm->alwaysdefault)
		/*
		 * Prevent shrinking 
		 */
		spm->minused = 0;

	for (array_iindex = 0; array_iindex < count; array_iindex++) {
		err = spmatrix_set(spm, array[array_iindex].key1, array[array_iindex].key2, array[array_iindex].value);
		if (err)
			MAPKIT_ERROR(err);
	}

	if (spm->alwaysdefault)
		spmatrix_adjustcapacity(spm);

	return MAPKIT_OK;
}

mapkit_error spmatrix_removeall(spmatrix * spm, spmatrix_key * array, mapkit_size_t count)
{
	mapkit_size_t array_iindex;
	mapkit_error err;

	/*
	 * Prevent shrinking 
	 */
	spm->minused = 0;

	for (array_iindex = 0; array_iindex < count; array_iindex++) {
		err = spmatrix_remove(spm, array[array_iindex].key1, array[array_iindex].key2);
		if (err)
			MAPKIT_ERROR(err);
	}

	spmatrix_adjustcapacity(spm);

	return MAPKIT_OK;
}

void spmatrix_printstats(spmatrix * spm)
{
	_spmatrix_printstats(spm);
}

#endif							       /* MAPKIT_spmatrix */
