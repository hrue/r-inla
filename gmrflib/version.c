
/* version.c
 * 
 * Copyright (C) 2007 Havard Rue
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
  \file version.c
  \brief Print out the RCSId's for the files in this build
*/

#ifndef HGVERSION
#define HGVERSION
#endif
static const char RCSId[] = "file: " __FILE__ "  " HGVERSION;

/* Pre-hg-Id: $Id: version.c,v 1.5 2007/03/01 23:46:00 hrue Exp $ */

#if !defined(__FreeBSD__)
#include <stddef.h>
#include <string.h>
#include <malloc.h>
#endif
#include <stdlib.h>

#include "GMRFLib/GMRFLib.h"
#include "GMRFLib/GMRFLibP.h"

int GMRFLib_version(FILE * fp)
{
	const char *RCSIds[] = {
#include "RCSId.all"
		NULL
	};
	int i = 0;

	fp = (fp ? fp : stdout);

	fprintf(fp, "GMRFLib %s\n", GMRFLib_VERSION);
	while (RCSIds[i]) {
		fprintf(fp, "\t%s\n", RCSIds[i]);
		i++;
	}

	return GMRFLib_SUCCESS;
}
