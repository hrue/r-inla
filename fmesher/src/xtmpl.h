/* Time-stamp: <98/01/06 10:30:04 hrue> */

#define XTMPL_OK 0
#define XTMPL_ERR (-1)
#define XTMPL_MAXWINDOWS 64

#undef __BEGIN_DECLS
#undef __END_DECLS
#ifdef __cplusplus
# define __BEGIN_DECLS extern "C" {
# define __END_DECLS }
#else
# define __BEGIN_DECLS                                         /* empty */
# define __END_DECLS                                           /* empty */
#endif

__BEGIN_DECLS


#ifndef XTMPL_INCLUDE_FROM_SOURCE_FILE
extern int xtmpl_window;
extern int xtmpl_view;
#endif

extern int xtmpl_clear(void);
extern int xtmpl_close(void);
extern int xtmpl_draw_line(int  , int  , int  , int );
extern int xtmpl_erase_line(int  , int  , int  , int );
extern int xtmpl_line(int  , int  , int  , int  , int );
extern int xtmpl_draw_text(int  , int  , char *, int );
extern int xtmpl_erase_text(int  , int  , char *, int );
extern int xtmpl_text(int, int  , int  , char *, int );
extern int xtmpl_line_width(int );
extern int xtmpl_open(int  , int  , char * );
extern int xtmpl_record(char *, char *);
extern int xtmpl_record_end(void);
extern int xtmpl_record_event(int, char *, char *, char *, char *);
extern int xtmpl_playback(char *, char *, double );
extern int xtmpl_dot(int, int, int, int);
extern int xtmpl_box(int, int, int, double, int);
extern int xtmpl_press_ret(char *do_what);

__END_DECLS
