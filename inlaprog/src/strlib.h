
/*-------------------------------------------------------------------------*/

/**
  @file     strlib.h
  @author   N. Devillard
  @date     Jan 2001
  @version  $Revision: 1.4 $
  @brief    Various string handling routines to complement the C lib.

  This modules adds a few complementary string routines usually missing
  in the standard C library.
*/

/*--------------------------------------------------------------------------*/

/*
	$Id: strlib.h,v 1.4 2008/08/16 20:13:43 hrue Exp $
	$Author: hrue $
	$Date: 2008/08/16 20:13:43 $
	$Revision: 1.4 $
*/

#ifndef _STRLIB_H_
#define _STRLIB_H_
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

/*---------------------------------------------------------------------------
   								Includes
 ---------------------------------------------------------------------------*/
#include <stdio.h>
#include <stdlib.h>

/*---------------------------------------------------------------------------
  							Function codes
 ---------------------------------------------------------------------------*/

/*-------------------------------------------------------------------------*/

/**
  @brief    Convert a string to lowercase.
  @param    s   String to convert.
  @return   ptr to statically allocated string.

  This function returns a pointer to a statically allocated string
  containing a lowercased version of the input string. Do not free
  or modify the returned string! Since the returned string is statically
  allocated, it will be modified at each function call (not re-entrant).
 */

/*--------------------------------------------------------------------------*/
char *strlwc(const char *s);

/*-------------------------------------------------------------------------*/

/**
  @brief    Convert a string to uppercase.
  @param    s   String to convert.
  @return   ptr to statically allocated string.

  This function returns a pointer to a statically allocated string
  containing an uppercased version of the input string. Do not free
  or modify the returned string! Since the returned string is statically
  allocated, it will be modified at each function call (not re-entrant).
 */

/*--------------------------------------------------------------------------*/
char *strupc(char *s);

/*-------------------------------------------------------------------------*/

/**
  @brief    Skip blanks until the first non-blank character.
  @param    s   String to parse.
  @return   Pointer to char inside given string.

  This function returns a pointer to the first non-blank character in the
  given string.
 */

/*--------------------------------------------------------------------------*/
char *strskp(char *s);

/*-------------------------------------------------------------------------*/

/**
  @brief    Remove blanks at the end of a string.
  @param    s   String to parse.
  @return   ptr to statically allocated string.

  This function returns a pointer to a statically allocated string,
  which is identical to the input string, except that all blank
  characters at the end of the string have been removed.
  Do not free or modify the returned string! Since the returned string
  is statically allocated, it will be modified at each function call
  (not re-entrant).
 */

/*--------------------------------------------------------------------------*/
char *strcrop(char *s);

/*-------------------------------------------------------------------------*/

/**
  @brief    Remove blanks at the beginning and the end of a string.
  @param    s   String to parse.
  @return   ptr to statically allocated string.

  This function returns a pointer to a statically allocated string,
  which is identical to the input string, except that all blank
  characters at the end and the beg. of the string have been removed.
  Do not free or modify the returned string! Since the returned string
  is statically allocated, it will be modified at each function call
  (not re-entrant).
 */

/*--------------------------------------------------------------------------*/
char *strstrip(char *s);


__END_DECLS
#endif
