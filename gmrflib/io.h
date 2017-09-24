
/* io.h
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
 *
 */

/*!
  \file io.h
  \brief Typedefs for input/output-routines
*/

#ifndef __GMRFLib_IO_H__
#define __GMRFLib_IO_H__

#include <stdlib.h>
#include <zlib.h>
#if !defined(__FreeBSD__)
#include <malloc.h>
#endif

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

/* 
   to define strtok_r()
*/
#if defined(__sun)
#define __EXTENSIONS__
#undef _STRING_H
#include <string.h>
#endif
#define GMRFLib_IO_COMMENT_CHAR "#"
#define GMRFLib_IO_SEP          " \t"

/* 
   internal error codes
*/
#define GMRFLib_IO_ERR_OPEN       1
#define GMRFLib_IO_ERR_NOLINE     2
#define GMRFLib_IO_ERR_READLINE   3
#define GMRFLib_IO_ERR_READBYTES  4
#define GMRFLib_IO_ERR_WRITEBYTES 5
    typedef struct {
	char *filename;
	char *mode;
	char *strtok_ptr;

	gzFile fp;
	size_t lines_read;
	size_t tokens_read;
	size_t bytes_read;
	size_t bytes_written;
} GMRFLib_io_tp;

int GMRFLib_io_close(GMRFLib_io_tp * io);
int GMRFLib_io_error(GMRFLib_io_tp * io, int error);
int GMRFLib_io_find_file_in_path(char **ptr, const char *filename, int must_find);
int GMRFLib_io_next_token(char **ptr, GMRFLib_io_tp * io);
int GMRFLib_io_nextline(char **ptr, GMRFLib_io_tp * io);
int GMRFLib_io_open(GMRFLib_io_tp ** io, const char *filename, const char *mode);
int GMRFLib_io_seek(GMRFLib_io_tp * io, size_t offset, int whence);
int GMRFLib_io_read_next(GMRFLib_io_tp * io, void *ptr, const char *fmt);
int GMRFLib_io_strip_blanks(char *line);
int GMRFLib_sprintf(char **ptr, const char *fmt, ...);
int GMRFLib_io_read(GMRFLib_io_tp * io, void *buf, size_t len);
int GMRFLib_io_write(GMRFLib_io_tp * io, const void *buf, size_t len);


__END_DECLS
#endif
