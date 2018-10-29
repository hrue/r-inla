
/* io.c
 * 
 * Copyright (C) 2005-2006 Havard Rue
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
  \file io.c
  \brief Functions for input and output
*/

#include <stdarg.h>
#include <string.h>
#include <stdio.h>
#if !defined(__FreeBSD__)
#include <malloc.h>
#endif
#include <stdlib.h>
#include "GMRFLib/GMRFLib.h"
#include "GMRFLib/GMRFLibP.h"

#ifndef HGVERSION
#define HGVERSION
#endif
static const char RCSId[] = "file: " __FILE__ "  " HGVERSION;

/* Pre-hg-Id: $Id: io.c,v 1.53 2009/05/23 06:16:16 hrue Exp $ */

int GMRFLib_sprintf(char **ptr, const char *fmt, ...)
{
	/*
	 * parts of this code is copied from the manual page of snprintf. 
	 */

	int n, size = 128 + 1;
	char *p;
	va_list ap;

	GMRFLib_ASSERT(ptr, GMRFLib_EINVARG);
	GMRFLib_ASSERT(fmt, GMRFLib_EINVARG);

	p = Calloc(size, char);

	while (1) {
		/*
		 * Try to print in the allocated space. 
		 */
		va_start(ap, fmt);
		n = vsnprintf(p, (unsigned int) size, fmt, ap);
		va_end(ap);

		/*
		 * if that worked, return the string, 
		 */
		if (n > -1 && n < size) {
			*ptr = p;
			return GMRFLib_SUCCESS;
		}

		/*
		 * ...else try again with more space 
		 */
		if (n > -1) {
			size = n + 1;
		} else {
			size *= 2;
		}
		p = Realloc(p, size, char);
	}

	return GMRFLib_SUCCESS;
}
int GMRFLib_io_find_file_in_path(char **ptr, const char *filename, int must_find)
{
	/*
	 * find readable file looking for files in GMRFLib_path. if must_find, then error if the file is not found (unless
	 * filename == NULL). this function return a malloc'ed string int *ptr with the filename (path included) to the file. 
	 */

	int found = 0;
	char *path = NULL, *p = NULL, *pp = NULL, *fnm = NULL, *strtok_ptr = NULL;
	const char *delim = ":";
	FILE *fp = NULL;

	GMRFLib_ASSERT(ptr, GMRFLib_EINVARG);
	*ptr = NULL;					       /* not found */

	if (filename) {
		p = getenv("GMRFLib_PATH");		       /* look for this, or */
		if (p) {
			GMRFLib_EWRAP0(GMRFLib_sprintf(&pp, "%s:%s", p, GMRFLib_path));
			p = pp;
		} else {
			p = GMRFLib_strdup(GMRFLib_path);      /* use default if it does not exists */
		}

		path = p;
		while ((pp = GMRFLib_strtok_r(path, delim, &strtok_ptr)) && !found) {
			path = NULL;
			GMRFLib_EWRAP0(GMRFLib_sprintf(&fnm, "%s/%s", pp, filename));
			if ((fp = fopen(fnm, "r"))) {
				found = 1;
				fclose(fp);
				*ptr = fnm;
			} else {
				Free(fnm);
			}
		}
		Free(p);

		if (must_find && !found) {
			GMRFLib_EWRAP0(GMRFLib_sprintf(&pp, "File [%s]", filename));
			GMRFLib_ERROR_MSG(GMRFLib_EOPENFILE, pp);
			Free(pp);
			return GMRFLib_EOPENFILE;
		}
	}

	return (found ? GMRFLib_SUCCESS : GMRFLib_EOPENFILE);
}
int GMRFLib_io_open(GMRFLib_io_tp ** io, const char *filename, const char *mode)
{
	GMRFLib_ASSERT(filename, GMRFLib_EPARAMETER);
	GMRFLib_ASSERT(io, GMRFLib_EPARAMETER);

	*io = Calloc(1, GMRFLib_io_tp);
	(*io)->filename = GMRFLib_strdup(filename);
	(*io)->mode = GMRFLib_strdup(mode);
	(*io)->strtok_ptr = NULL;
	(*io)->lines_read = 0;
	(*io)->tokens_read = 0;
	(*io)->bytes_read = 0;
	(*io)->bytes_written = 0;
	(*io)->fp = gzopen(filename, mode);

	if (!((*io)->fp)) {
		return GMRFLib_io_error(*io, GMRFLib_IO_ERR_OPEN);
	} else {
		return GMRFLib_SUCCESS;
	}
}
int GMRFLib_io_seek(GMRFLib_io_tp * io, size_t offset, int whence)
{
	/*
	 * whence is one of SEEK_SET and SEEK_CURRENT. SEEK_END is not supported. 
	 */

	return ((int) gzseek(io->fp, (z_off_t) offset, whence));
}
int GMRFLib_io_close(GMRFLib_io_tp * io)
{
	if (io) {
		if (io->fp) {
			gzclose(io->fp);
		}
		Free(io->filename);
		Free(io->mode);
		Free(io);
		GMRFLib_io_next_token(NULL, NULL);	       /* special: reset */
	}
	return GMRFLib_SUCCESS;
}
int GMRFLib_io_error(GMRFLib_io_tp * io, int error)
{
	char *msg = NULL;
	int ecode;

	switch (error) {
	case GMRFLib_IO_ERR_OPEN:
		GMRFLib_EWRAP0(GMRFLib_sprintf(&msg, "Fail to open file[%s] with mode[%s]", io->filename, io->mode));
		ecode = GMRFLib_EOPENFILE;
		break;
	case GMRFLib_IO_ERR_NOLINE:
		GMRFLib_EWRAP0(GMRFLib_sprintf(&msg, "Fail to read line[%1d] in file[%s]", io->lines_read + 1, io->filename));
		ecode = GMRFLib_EREADFILE;
		break;
	case GMRFLib_IO_ERR_READLINE:
		if (io->tokens_read) {
			GMRFLib_EWRAP0(GMRFLib_sprintf
				       (&msg, "Fail to read from or get, line[%1d] token[%1d] in file[%s]", io->lines_read, io->tokens_read,
					io->filename));
		} else {
			GMRFLib_EWRAP0(GMRFLib_sprintf(&msg, "Fail to read from or get, line[%1d] in file[%s]", io->lines_read, io->filename));
		}
		ecode = GMRFLib_EREADFILE;
		break;
	case GMRFLib_IO_ERR_READBYTES:
		GMRFLib_EWRAP0(GMRFLib_sprintf(&msg, "Fail to read more after [%1d] bytes are read in file[%s]", io->bytes_read, io->filename));
		ecode = GMRFLib_EREADFILE;
		break;
	case GMRFLib_IO_ERR_WRITEBYTES:
		GMRFLib_EWRAP0(GMRFLib_sprintf
			       (&msg, "Fail to write more after [%1d] bytes are written to file[%s]", io->bytes_written, io->filename));
		ecode = GMRFLib_EWRITE;
		break;

	default:
		GMRFLib_ASSERT(0 == 1, GMRFLib_ESNH);
	}

	GMRFLib_ERROR_MSG(ecode, msg);
	Free(msg);

	return ecode;
}
int GMRFLib_io_strip_blanks(char *line)
{
	char *ptr = NULL;
	int i;

	if (!line)
		return GMRFLib_SUCCESS;

	ptr = line;					       /* strip leading white-spaces (including CR) */
	while ((!strncmp(ptr, " ", 1) || !strncmp(ptr, "\t", 1) || !strncmp(ptr, "\r", 1)) && *ptr) {
		ptr++;
	}

	for (i = 0; i < (int) strlen(ptr) + 1; i++) {
		line[i] = ptr[i];			       /* +1 : yes */
	}

	if (!strlen(ptr)) {
		line = ptr;
		return GMRFLib_SUCCESS;
	}

	ptr = line + strlen(line) - 1;			       /* strip trailing white-spaces (including CR) */
	while ((!strncmp(ptr, " ", 1) || !strncmp(ptr, "\t", 1) || !strncmp(ptr, "\r", 1)) && *ptr && ptr > line) {
		ptr--;
	}

	line[ptr - line + 1] = '\0';

	return GMRFLib_SUCCESS;
}
int GMRFLib_io_nextline(char **ptr, GMRFLib_io_tp * io)
{
	char *line = NULL;
	int maxlen = 4096, read_ok;
	z_off_t tell;

	*ptr = NULL;					       /* default */

	line = Calloc(maxlen, char);

	do {
		/*
		 * read until a non-empty line and non-comment line is found 
		 */
		read_ok = 0;
		do {
			/*
			 * read one line into `line', we need to be careful about long lines... 
			 */
			tell = gztell(io->fp);
			if (gzgets(io->fp, line, maxlen) == Z_NULL) {
				/*
				 * nothing to read 
				 */
				Free(line);
				return GMRFLib_SUCCESS;
			}

			if ((int) strlen(line) == maxlen - 1) {
				/*
				 * line is to short, increase length and reread line 
				 */
				maxlen *= 4;
				line = Realloc(line, maxlen, char);

				gzseek(io->fp, tell, SEEK_SET);
			} else {
				read_ok = 1;
			}
		}
		while (!read_ok);

		io->lines_read++;
	}
	while (!strncmp(GMRFLib_IO_COMMENT_CHAR, line, 1) || strlen(line) <= 1);

	line[strlen(line) - 1] = '\0';			       /* ignore the newline */
	GMRFLib_io_strip_blanks(line);			       /* no need to wrap */

	if (!strlen(line)) {
		/*
		 * this is a line with spaces/tabs, return then the next one. 
		 */
		Free(line);
		return GMRFLib_io_nextline(ptr, io);
	} else {
		*ptr = line;
	}

	io->tokens_read = 0;

	return GMRFLib_SUCCESS;
}
int GMRFLib_io_next_token(char **ptr, GMRFLib_io_tp * io)
{
	char *tok = NULL;
	static char *line = NULL;
#pragma omp threadprivate(line)

	if (io == NULL) {				       /* special: reset strtok */
		Free(line);
		return GMRFLib_SUCCESS;
	}

	*ptr = NULL;					       /* default */

	if (io->strtok_ptr) {
		tok = GMRFLib_strtok_r(NULL, GMRFLib_IO_SEP, &(io->strtok_ptr));
		if (tok) {
			io->tokens_read++;
			*ptr = tok;
			return GMRFLib_SUCCESS;
		}
	} else {
		tok = NULL;
	}

	/*
	 * no token. read next line and extract the first token 
	 */
	Free(line);

	GMRFLib_EWRAP0(GMRFLib_io_nextline(&line, io));
	if (line) {
		tok = GMRFLib_strtok_r(line, GMRFLib_IO_SEP, &(io->strtok_ptr));
		io->tokens_read++;
		*ptr = tok;
		return GMRFLib_SUCCESS;
	} else {
		*ptr = NULL;
		return GMRFLib_EREADFILE;
	}

	return GMRFLib_SUCCESS;
}
int GMRFLib_io_read_next(GMRFLib_io_tp * io, void *ptr, const char *fmt)
{
	/*
	 * read next int/double etc, signalling an error if not successful 
	 */
	char *token = NULL;

	GMRFLib_ASSERT(io, GMRFLib_EINVARG);

	GMRFLib_io_next_token(&token, io);		       /* do NOT check error here */
	if (token) {
		if (sscanf(token, fmt, ptr) != 1) {
			return GMRFLib_io_error(io, GMRFLib_IO_ERR_READLINE);
		} else {
			return GMRFLib_SUCCESS;
		}
	} else {
		return GMRFLib_io_error(io, GMRFLib_IO_ERR_NOLINE);
	}

	GMRFLib_ASSERT(1 == 0, GMRFLib_ESNH);

	return GMRFLib_SUCCESS;
}
int GMRFLib_io_read(GMRFLib_io_tp * io, void *buf, size_t len)
{
	/*
	 * binary: read the next LEN bytes into BUF 
	 */
	unsigned int nr;

	if (io && len) {
		nr = gzread(io->fp, buf, (unsigned int) len);
		if (nr != len) {
			return GMRFLib_io_error(io, GMRFLib_IO_ERR_READBYTES);
		} else {
			io->bytes_read += (size_t) len;
		}
	}
	return GMRFLib_SUCCESS;
}
int GMRFLib_io_write(GMRFLib_io_tp * io, const void *buf, size_t len)
{
	/*
	 * binary: write LEN bytes in BUF to file 
	 */
	unsigned int nw;

	if (io && len) {
		nw = gzwrite(io->fp, buf, (unsigned int) len);
		if (nw != len) {
			return GMRFLib_io_error(io, GMRFLib_IO_ERR_WRITEBYTES);
		} else {
			io->bytes_written += (size_t) len;
		}
	}
	return GMRFLib_SUCCESS;
}
