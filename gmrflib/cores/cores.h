
/*
 * Copyright (c) Przemyslaw Skibinski, Yann Collet, Facebook, Inc.
 * All rights reserved.
 *
 * This source code is licensed under both the BSD-style license (found in the
 * LICENSE file in the root directory of this source tree) and the GPLv2 (found
 * in the COPYING file in the root directory of this source tree).
 * You may select, at your option, one of the above-listed licenses.
 */

#ifndef CORES_H_MODULE
#define CORES_H_MODULE

#undef __BEGIN_DECLS
#undef __END_DECLS
#ifdef __cplusplus
#define __BEGIN_DECLS extern "C" {
#define __END_DECLS }
#else
#define __BEGIN_DECLS					       /* empty */
#define __END_DECLS					       /* empty */
#endif

__BEGIN_DECLS

/*-****************************************
*  Dependencies
******************************************/
#include <stddef.h>					       /* size_t, ptrdiff_t */
#include <sys/types.h>					       /* stat, utime */
#include <sys/stat.h>					       /* stat, chmod */

/***************************************************************
*  Basic Types
*****************************************************************/
#if  !defined (__VMS) && (defined (__cplusplus) || (defined (__STDC_VERSION__) && (__STDC_VERSION__ >= 199901L) /* C99 */) )
#if defined(_AIX)
#include <inttypes.h>
#else
#include <stdint.h>					       /* intptr_t */
#endif
typedef uint8_t BYTE;
typedef uint8_t U8;
typedef int8_t S8;
typedef uint16_t U16;
typedef int16_t S16;
typedef uint32_t U32;
typedef int32_t S32;
typedef uint64_t U64;
typedef int64_t S64;
#else
#include <limits.h>
#if CHAR_BIT != 8
#error "this implementation requires char to be exactly 8-bit type"
#endif
typedef unsigned char BYTE;
typedef unsigned char U8;
typedef signed char S8;
#if USHRT_MAX != 65535
#error "this implementation requires short to be exactly 16-bit type"
#endif
typedef unsigned short U16;
typedef signed short S16;
#if UINT_MAX != 4294967295
#error "this implementation requires int to be exactly 32-bit type"
#endif
typedef unsigned int U32;
typedef signed int S32;

/* note : there are no limits defined for long long type in C90.
 * limits exist in C99, however, in such case, <stdint.h> is preferred */
typedef unsigned long long U64;
typedef signed long long S64;
#endif

int UTIL_countCores(int logical);
int UTIL_countPhysicalCores(void);
int UTIL_countLogicalCores(void);

__END_DECLS
#endif
