#ifndef __GMRFLib_IO_H__
#define __GMRFLib_IO_H__

#include <stdlib.h>
#include <zlib.h>

/* 
   to define strtok_r()
*/
#if defined(__sun)
#define __EXTENSIONS__
#undef _STRING_H
#include <string.h>
#endif

#include "GMRFLib/GMRFLib.h"
#include "GMRFLib/GMRFLibP.h"

#undef __BEGIN_DECLS
#undef __END_DECLS
#ifdef __cplusplus
#define __BEGIN_DECLS extern "C" {
#define __END_DECLS }
#else
#define __BEGIN_DECLS					       /* empty */
#define __END_DECLS					       /* empty */
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

__BEGIN_DECLS typedef struct {
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
int GMRFLib_io_read(GMRFLib_io_tp * io, void *buf, size_t len);
int GMRFLib_io_write(GMRFLib_io_tp * io, const void *buf, size_t len);

intmax_t GMRFLib_io_file_size(const char *filename);

__END_DECLS
#endif
