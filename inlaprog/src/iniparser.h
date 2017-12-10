
/**
   @file    iniparser.h
   @author  N. Devillard
   @date    Mar 2000
   @version $Revision: 1.7 $
   @brief   Parser for ini files.
*/

/*--------------------------------------------------------------------------*/

/*
	$Id: iniparser.h,v 1.7 2009/03/30 16:47:05 hrue Exp $
	$Author: hrue $
	$Date: 2009/03/30 16:47:05 $
	$Revision: 1.7 $
*/

#ifndef _INIPARSER_H_
#define _INIPARSER_H_

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
#include <string.h>

/* #include <unistd.h> */
#include "dictionary.h"
#define INIPARSER_SECSEP "!"  // the character used to define section names
#define INIPARSER_SEP '!'     // the internal separator
#define LEN_INIPARSER_SEP ((size_t)1)
int iniparser_getnsec(dictionary * d);
char *iniparser_getsecname(dictionary * d, int n);
void iniparser_dump_ini(dictionary * d, FILE * f);
void iniparser_dump(dictionary * d, FILE * f);
char *iniparser_getstr(dictionary * d, const char *key);
char *iniparser_getstring(dictionary * d, const char *key, char *def);
int iniparser_getint(dictionary * d, const char *key, int notfound);
double iniparser_getdouble(dictionary * d, const char *key, double notfound);
int iniparser_getboolean(dictionary * d, const char *key, int notfound);
int iniparser_setstr(dictionary * ini, char *entry, char *val);
void iniparser_unset(dictionary * ini, char *entry);
int iniparser_find_entry(dictionary * ini, char *entry);
char *iniparser_getline(FILE * f);
dictionary *iniparser_load(const char *ininame);
void iniparser_freedict(dictionary * d);

__END_DECLS
#endif
