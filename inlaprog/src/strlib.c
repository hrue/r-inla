
/*-------------------------------------------------------------------------*/

/**
  @file		strlib.c
  @author	N. Devillard
  @date		Jan 2001
  @version	$Revision: 1.5 $
  @brief	Various string handling routines to complement the C lib.

  This modules adds a few complementary string routines usually missing
  in the standard C library.
*/

/*--------------------------------------------------------------------------*/

/*
	$Id: strlib.c,v 1.5 2008/07/17 10:11:36 hrue Exp $
	$Author: hrue $
	$Date: 2008/07/17 10:11:36 $
	$Revision: 1.5 $
*/

/*---------------------------------------------------------------------------
   								Includes
 ---------------------------------------------------------------------------*/

#include <assert.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include "strlib.h"

/*---------------------------------------------------------------------------
   							    Defines	
 ---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------
  							Function codes
 ---------------------------------------------------------------------------*/

/*-------------------------------------------------------------------------*/

/**
  @brief	Convert a string to lowercase.
  @param	s	String to convert.
  @return	ptr to statically allocated string.

  This function returns a pointer to a statically allocated string
  containing a lowercased version of the input string. Do not free
  or modify the returned string! Since the returned string is statically
  allocated, it will be modified at each function call (not re-entrant).
 */

/*--------------------------------------------------------------------------*/

char *strlwc(const char *s)
{
	if (s == NULL) {
		return NULL;
	}

	static char *l = NULL;
	l = (char *) realloc(l, (size_t) (strlen(s) + 1) * sizeof(char));
	assert(l);

	size_t i;
	i = 0;
	while (i < strlen(s) && s[i]) {
		l[i] = (char) tolower((int) s[i]);
		i++;
	}
	l[strlen(s)] = (char) 0;
	return l;

}

/*-------------------------------------------------------------------------*/

/**
  @brief	Convert a string to uppercase.
  @param	s	String to convert.
  @return	ptr to statically allocated string.

  This function returns a pointer to a statically allocated string
  containing an uppercased version of the input string. Do not free
  or modify the returned string! Since the returned string is statically
  allocated, it will be modified at each function call (not re-entrant).
 */

/*--------------------------------------------------------------------------*/

char *strupc(char *s)
{
	if (s == NULL) {
		return NULL;
	}

	assert(s);
	static char *l = NULL;
	l = (char *) realloc(l, (size_t) (strlen(s) + 1) * sizeof(char));
	assert(l);
	size_t i;
	i = 0;
	while (i < strlen(s) && s[i]) {
		l[i] = (char) toupper((int) s[i]);
		i++;
	}
	l[strlen(s)] = (char) 0;
	return l;
}

/*-------------------------------------------------------------------------*/

/**
  @brief	Skip blanks until the first non-blank character.
  @param	s	String to parse.
  @return	Pointer to char inside given string.

  This function returns a pointer to the first non-blank character in the
  given string.
 */

/*--------------------------------------------------------------------------*/

char *strskp(char *s)
{
	char *skip = s;

	if (s == NULL)
		return NULL;
	while (isspace((int) *skip) && *skip)
		skip++;
	return skip;
}

/*-------------------------------------------------------------------------*/

/**
  @brief	Remove blanks at the end of a string.
  @param	s	String to parse.
  @return	ptr to statically allocated string.

  This function returns a pointer to a statically allocated string,
  which is identical to the input string, except that all blank
  characters at the end of the string have been removed.
  Do not free or modify the returned string! Since the returned string
  is statically allocated, it will be modified at each function call
  (not re-entrant).
 */

/*--------------------------------------------------------------------------*/

char *strcrop(char *s)
{
	if (s == NULL) {
		return NULL;
	}

	assert(s);
	static char *l = NULL;
	l = (char *) realloc(l, (size_t) (strlen(s) + 1) * sizeof(char));
	assert(l);

	char *last;
	strcpy(l, s);
	last = l + strlen(l);
	while (last > l) {
		if (!isspace((int) *(last - 1)))
			break;
		last--;
	}
	*last = (char) 0;
	return l;
}

/*-------------------------------------------------------------------------*/

/* Test code */
#ifdef TEST
int main(int argc, char *argv[])
{
	char *str;

	str = "\t\tI'm a lumberkack and I'm OK      ";
	printf("lowercase: [%s]\n", strlwc(str));
	printf("uppercase: [%s]\n", strupc(str));
	printf("skipped  : [%s]\n", strskp(str));
	printf("cropped  : [%s]\n", strcrop(str));

	return 0;
}
#endif

/* vim: set ts=4 et sw=4 tw=75 */
