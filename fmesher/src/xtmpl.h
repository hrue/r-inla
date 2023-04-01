#ifdef FMESHER_WITH_X

/* Time-stamp: <98/01/06 10:30:04 hrue> */

#define XTMPL_OK 0
#define XTMPL_ERR (-1)
#define XTMPL_MAXWINDOWS 64

#undef __BEGIN_DECLS
#undef __END_DECLS
#ifdef __cplusplus
#define __BEGIN_DECLS extern "C" {
#define __END_DECLS }
#else
#define __BEGIN_DECLS /* empty */
#define __END_DECLS   /* empty */
#endif

__BEGIN_DECLS

#ifndef XTMPL_INCLUDE_FROM_SOURCE_FILE
extern int xtmpl_window;
extern int xtmpl_view;
#endif

extern int xtmpl_open(int nx, int ny, const char *title);
extern int xtmpl_line_width(int lw);
extern int xtmpl_clear();
extern int xtmpl_draw_text(int x0, int y0, char *str, int nstr);
extern int xtmpl_erase_text(int x0, int y0, char *str, int nstr);
extern int xtmpl_text(int fg, int x0, int y0, char *str, int nstr);
extern int xtmpl_line(int fg, int x0, int y0, int x1, int y1);
extern int xtmpl_draw_line(int x0, int y0, int x1, int y1);
extern int xtmpl_erase_line(int x0, int y0, int x1, int y1);
extern int xtmpl_dot(int x0, int y0, int width, int fg);
extern int xtmpl_box(int x0, int y0, int width, double yoffset, int fg);
extern int xtmpl_close();
extern int xtmpl_record(const char *filename, const char *channels);
extern int xtmpl_record_end();
extern int xtmpl_record_event(int type, char *a0, char *a1, const char *a2,
                              char *a3);
extern int xtmpl_press_ret(const char *do_what);
extern int xtmpl_playback(const char *filename, const char *channels,
                          double magnify);

__END_DECLS

#endif
