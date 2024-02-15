/*
https://github.com/scandum/fluxsort

This is free and unencumbered software released into the public domain.

Anyone is free to copy, modify, publish, use, compile, sell, or
distribute this software, either in source code form or as a compiled
binary, for any purpose, commercial or non-commercial, and by any
means.

In jurisdictions that recognize copyright laws, the author or authors
of this software dedicate any and all copyright interest in the
software to the public domain. We make this dedication for the benefit
of the public at large and to the detriment of our heirs and
successors. We intend this dedication to be an overt act of
relinquishment in perpetuity of all present and future rights to this
software under copyright law.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT SHALL THE AUTHORS BE LIABLE FOR ANY CLAIM, DAMAGES OR
OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE,
ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
OTHER DEALINGS IN THE SOFTWARE.

For more information, please refer to <http://unlicense.org>
*/

// fluxsort 1.2.1.3 - Igor van den Hoven ivdhoven@gmail.com

#ifndef FLUXSORT_H
#define FLUXSORT_H

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <errno.h>
#include <float.h>
#include <string.h>

typedef int CMPFUNC(const void *a, const void *b);

#ifndef QUADSORT_H
#include "GMRFLib/fsort/quadsort.h"
#endif

// When sorting an array of 32/64 bit pointers, like a string array, QUAD_CACHE
// needs to be adjusted in quadsort.h and here for proper performance when
// sorting large arrays.

#ifdef cmp
#define QUAD_CACHE 4294967295
#else
//#define QUAD_CACHE 131072
#define QUAD_CACHE 262144
//#define QUAD_CACHE 524288
//#define QUAD_CACHE 4294967295
#endif

//////////////////////////////////////////////////////////
// ┌───────────────────────────────────────────────────┐//
// │       ██████┐ ██████┐    ██████┐ ██████┐████████┐ │//
// │       └────██┐└────██┐   ██┌──██┐└─██┌─┘└──██┌──┘ │//
// │        █████┌┘ █████┌┘   ██████┌┘  ██│     ██│    │//
// │        └───██┐██┌───┘    ██┌──██┐  ██│     ██│    │//
// │       ██████┌┘███████┐   ██████┌┘██████┐   ██│    │//
// │       └─────┘ └──────┘   └─────┘ └─────┘   └─┘    │//
// └───────────────────────────────────────────────────┘//
//////////////////////////////////////////////////////////

#define VAR int
#define FUNC(NAME) NAME##32

#include "GMRFLib/fsort/fluxsort.c"

#undef VAR
#undef FUNC

// fluxsort_prim

#define VAR int
#define FUNC(NAME) NAME##_int32
#ifndef cmp
#define cmp(a,b) (*(a) > *(b))
#include "GMRFLib/fsort/fluxsort.c"
#undef cmp
#else
#include "GMRFLib/fsort/fluxsort.c"
#endif
#undef VAR
#undef FUNC

#define VAR unsigned int
#define FUNC(NAME) NAME##_uint32
#ifndef cmp
#define cmp(a,b) (*(a) > *(b))
#include "GMRFLib/fsort/fluxsort.c"
#undef cmp
#else
#include "GMRFLib/fsort/fluxsort.c"
#endif
#undef VAR
#undef FUNC

//////////////////////////////////////////////////////////
// ┌───────────────────────────────────────────────────┐//
// │        █████┐ ██┐  ██┐   ██████┐ ██████┐████████┐ │//
// │       ██┌───┘ ██│  ██│   ██┌──██┐└─██┌─┘└──██┌──┘ │//
// │       ██████┐ ███████│   ██████┌┘  ██│     ██│    │//
// │       ██┌──██┐└────██│   ██┌──██┐  ██│     ██│    │//
// │       └█████┌┘     ██│   ██████┌┘██████┐   ██│    │//
// │        └────┘      └─┘   └─────┘ └─────┘   └─┘    │//
// └───────────────────────────────────────────────────┘//
//////////////////////////////////////////////////////////

#define VAR long long
#define FUNC(NAME) NAME##64

#include "GMRFLib/fsort/fluxsort.c"

#undef VAR
#undef FUNC

// fluxsort_prim

#define VAR long long
#define FUNC(NAME) NAME##_int64
#ifndef cmp
#define cmp(a,b) (*(a) > *(b))
#include "GMRFLib/fsort/fluxsort.c"
#undef cmp
#else
#include "GMRFLib/fsort/fluxsort.c"
#endif
#undef VAR
#undef FUNC

#define VAR unsigned long long
#define FUNC(NAME) NAME##_uint64
#ifndef cmp
#define cmp(a,b) (*(a) > *(b))
#include "GMRFLib/fsort/fluxsort.c"
#undef cmp
#else
#include "GMRFLib/fsort/fluxsort.c"
#endif
#undef VAR
#undef FUNC

// This section is outside of 32/64 bit pointer territory, so no cache checks
// necessary, unless sorting 32+ byte structures.

#undef QUAD_CACHE
#define QUAD_CACHE 4294967295

//////////////////////////////////////////////////////////
//┌────────────────────────────────────────────────────┐//
//│                █████┐    ██████┐ ██████┐████████┐  │//
//│               ██┌──██┐   ██┌──██┐└─██┌─┘└──██┌──┘  │//
//│               └█████┌┘   ██████┌┘  ██│     ██│     │//
//│               ██┌──██┐   ██┌──██┐  ██│     ██│     │//
//│               └█████┌┘   ██████┌┘██████┐   ██│     │//
//│                └────┘    └─────┘ └─────┘   └─┘     │//
//└────────────────────────────────────────────────────┘//
//////////////////////////////////////////////////////////

#define VAR char
#define FUNC(NAME) NAME##8

#include "GMRFLib/fsort/fluxsort.c"

#undef VAR
#undef FUNC

//////////////////////////////////////////////////////////
//┌────────────────────────────────────────────────────┐//
//│           ▄██┐   █████┐    ██████┐ ██████┐████████┐│//
//│          ████│  ██┌───┘    ██┌──██┐└─██┌─┘└──██┌──┘│//
//│          └─██│  ██████┐    ██████┌┘  ██│     ██│   │//
//│            ██│  ██┌──██┐   ██┌──██┐  ██│     ██│   │//
//│          ██████┐└█████┌┘   ██████┌┘██████┐   ██│   │//
//│          └─────┘ └────┘    └─────┘ └─────┘   └─┘   │//
//└────────────────────────────────────────────────────┘//
//////////////////////////////////////////////////////////

#define VAR short
#define FUNC(NAME) NAME##16

#include "GMRFLib/fsort/fluxsort.c"

#undef VAR
#undef FUNC



//////////////////////////////////////////////////////////
//┌────────────────────────────────────────────────────┐//
//│  ▄██┐  ██████┐  █████┐    ██████┐ ██████┐████████┐ │//
//│ ████│  └────██┐██┌──██┐   ██┌──██┐└─██┌─┘└──██┌──┘ │//
//│ └─██│   █████┌┘└█████┌┘   ██████┌┘  ██│     ██│    │//
//│   ██│  ██┌───┘ ██┌──██┐   ██┌──██┐  ██│     ██│    │//
//│ ██████┐███████┐└█████┌┘   ██████┌┘██████┐   ██│    │//
//│ └─────┘└──────┘ └────┘    └─────┘ └─────┘   └─┘    │//
//└────────────────────────────────────────────────────┘//
//////////////////////////////////////////////////////////

#if (DBL_MANT_DIG < LDBL_MANT_DIG)
#define VAR long double
#define FUNC(NAME) NAME##128
#include "GMRFLib/fsort/fluxsort.c"
#undef VAR
#undef FUNC
#endif

//////////////////////////////////////////////////////////////////////////
//┌────────────────────────────────────────────────────────────────────┐//
//│███████┐██┐     ██┐   ██┐██┐  ██┐███████┐ ██████┐ ██████┐ ████████┐ │//
//│██┌────┘██│     ██│   ██│└██┐██┌┘██┌────┘██┌───██┐██┌──██┐└──██┌──┘ │//
//│█████┐  ██│     ██│   ██│ └███┌┘ ███████┐██│   ██│██████┌┘   ██│    │//
//│██┌──┘  ██│     ██│   ██│ ██┌██┐ └────██│██│   ██│██┌──██┐   ██│    │//
//│██│     ███████┐└██████┌┘██┌┘ ██┐███████│└██████┌┘██│  ██│   ██│    │//
//│└─┘     └──────┘ └─────┘ └─┘  └─┘└──────┘ └─────┘ └─┘  └─┘   └─┘    │//
//└────────────────────────────────────────────────────────────────────┘//
//////////////////////////////////////////////////////////////////////////

void fluxsort(void *array, size_t nmemb, size_t size, CMPFUNC *cmp)
{
	if (nmemb < 2) {
		return;
	}

	switch (size) {
	case sizeof(char):
		fluxsort8(array, nmemb, cmp);
		break;

	case sizeof(short):
		fluxsort16(array, nmemb, cmp);
		break;

	case sizeof(int):
		fluxsort32(array, nmemb, cmp);
		break;

	case sizeof(long long):
		fluxsort64(array, nmemb, cmp);
		break;
#if (DBL_MANT_DIG < LDBL_MANT_DIG)
	case sizeof(long double):
		fluxsort128(array, nmemb, cmp);
		break;
#endif

	default:
		qsort(array, nmemb, size, cmp);
	}
	return;
}

// This must match quadsort_prim()

void fluxsort_prim(void *array, size_t nmemb, size_t size)
{
	if (nmemb < 2) {
		return;
	}

	switch (size) {
	case 4:
		fluxsort_int32(array, nmemb, NULL);
		return;
	case 5:
		fluxsort_uint32(array, nmemb, NULL);
		return;
	case 8:
		fluxsort_int64(array, nmemb, NULL);
		return;
	case 9:
		fluxsort_uint64(array, nmemb, NULL);
		return;
	default:
		assert(size == sizeof(int) || size == sizeof(int) + 1 || size == sizeof(long long) || size == sizeof(long long) + 1);
		return;
	}
}

// Sort arrays of structures, the comparison function must be by reference.

void fluxsort_size(void *array, size_t nmemb, size_t size, CMPFUNC *cmp)
{
	char **pti, *pta, *pts;
	size_t index, offset;

	pta = (char *) array;
	pti = (char **) malloc(nmemb * sizeof(char *));

	assert(pti != NULL);

	for (index = offset = 0; index < nmemb; index++) {
		pti[index] = pta + offset;
		offset += size;
	}

	switch (sizeof(size_t)) {
	case 4:
		fluxsort32(pti, nmemb, cmp);
		break;
	case 8:
		fluxsort64(pti, nmemb, cmp);
		break;
	}

	pts = (char *) malloc(nmemb * size);

	assert(pts != NULL);

	for (index = 0; index < nmemb; index++) {
		memcpy(pts, pti[index], size);
		pts += size;
	}
	pts -= nmemb * size;
	memcpy(array, pts, nmemb * size);

	free(pti);
	free(pts);
}

#undef QUAD_CACHE
#endif // #ifndef FLUXSORT_H
