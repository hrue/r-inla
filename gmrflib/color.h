#ifndef __GMRFLib_COLOR_H__
#       define __GMRFLib_COLOR_H__

#       include <time.h>
#       include <math.h>
#       include <stdlib.h>
#       include <assert.h>

#       undef __BEGIN_DECLS
#       undef __END_DECLS
#       ifdef __cplusplus
#              define __BEGIN_DECLS extern "C" {
#              define __END_DECLS }
#       else
#              define __BEGIN_DECLS			       /* empty */
#              define __END_DECLS			       /* empty */
#       endif

__BEGIN_DECLS
//
    typedef enum {
	COLOR_RESET,
	COLOR_BG_RED_FG_WHITE,
	COLOR_BG_BLUE_FG_WHITE,
	COLOR_BG_GREEN_FG_WHITE
} GMRFLib_color_tp;

// not yet convinced that this is the way to do this
#if defined(INLA_WITH_COLOR)
#       if defined(NDEBUG)
#              define _assert(condition) ((void)0)
#       else
#              define _assert(condition) ((void)((condition) ? 0 : \
						 (GMRFLib_color_start_error(), \
						  assert(condition), \
						  0)))
#       endif
#else
#define _assert(...) assert(__VA_ARGS__)
#endif

void GMRFLib_color_start_error(void);
void GMRFLib_color_start_warning(void);
void GMRFLib_color_start_info(void);
void GMRFLib_color_reset(void);

#       define LOG_ERROR(fmt, ...)					\
	do {								\
		GMRFLib_color_start_error();				\
		fprintf(stderr, "*** ERROR *** [%1d] %s:%1d (%s)\n",	\
			omp_get_thread_num(), __FILE__,  __LINE__, __GMRFLib_FuncName); \
		fprintf(stderr, "*** ERROR *** " fmt "", ##__VA_ARGS__); \
		fprintf(stderr, "\n");					\
		GMRFLib_color_reset();					\
	} while(0)

#       define LOG_WARNING(fmt, ...)					\
	do {								\
		GMRFLib_color_start_warning();				\
		fprintf(stderr, "*** WARNING *** [%1d] %s:%1d (%s)\n", omp_get_thread_num(), \
			__FILE__,  __LINE__, __GMRFLib_FuncName);	\
		fprintf(stderr, "*** WARNING *** " fmt "", ##__VA_ARGS__); \
		fprintf(stderr, "\n");					\
		GMRFLib_color_reset();					\
	} while(0)

#       define LOG_INFO(fmt, ...)					\
	do {								\
		GMRFLib_color_start_info();				\
		fprintf(stderr, "*** INFO *** [%1d] %s:%1d (%s)\n", omp_get_thread_num(), \
			__FILE__,  __LINE__, __GMRFLib_FuncName);	\
		fprintf(stderr, "*** INFO *** " fmt "", ##__VA_ARGS__);	\
		fprintf(stderr, "\n");					\
		GMRFLib_color_reset();					\
	} while(0)

#       define LOG_ERROR_WRAPPER(...) \
	do {				     \
		GMRFLib_color_start_error(); \
		__VA_ARGS__;		     \
		GMRFLib_color_reset();	     \
	} while (0)

#       define LOG_INFO_WRAPPER(...) \
	do {				     \
		GMRFLib_color_start_info(); \
		__VA_ARGS__;		     \
		GMRFLib_color_reset();	     \
	} while (0)

#       define LOG_WARNING_WRAPPER(...) \
	do {				     \
		GMRFLib_color_start_warning();	\
		__VA_ARGS__;			\
		GMRFLib_color_reset();	     \
	} while (0)

void GMRFLib_color_init(void);
void GMRFLib_color_set(GMRFLib_color_tp color);

__END_DECLS
#endif
