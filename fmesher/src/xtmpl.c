#ifdef FMESHER_WITH_X

/* Time-stamp: <98/01/06 10:36:34 hrue> */
#include <X11/Xlib.h>
#include <malloc.h>
#include <math.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define XTMPL_INCLUDE_FROM_SOURCE_FILE
#include "xtmpl.h"
#undef XTMPL_INCLUDE_FROM_SOURCE_FILE

int xtmpl_window = 0;                                   /* current X--window */
#define WIN(what) what[xtmpl_window % XTMPL_MAXWINDOWS] /* helpful macro */

#define PERR                                                                   \
  fprintf(stderr, "\nxtmpl_playback: file[%s] cmd[%d] byte[%ld] in error\n",   \
          filename, eno, ftell(fp));

static Display *dpy[XTMPL_MAXWINDOWS];
static int dim_x[XTMPL_MAXWINDOWS], dim_y[XTMPL_MAXWINDOWS];
static Window window[XTMPL_MAXWINDOWS];
static int screen[XTMPL_MAXWINDOWS];
static XEvent xevent;
static GC draw_gc[XTMPL_MAXWINDOWS], erase_gc[XTMPL_MAXWINDOWS];
static XGCValues gcv[XTMPL_MAXWINDOWS];
static unsigned int line_width[XTMPL_MAXWINDOWS];

static FILE *rec_fp = NULL;
static char *rec_channels = NULL;
int xtmpl_view = 1; /* if FALSE, don't draw or open windows */

#define DO_BYTESWAP                                                            \
  (!big_endian_query()) /* write file in big-endian format                     \
                         */

static int forever = 0; /* used by the `playback' program */
#define OPEN 0
#define NEW_WIN 1
#define TITLE 2
#define LINE_W 3
#define CLEAR 4
#define DRAW 5
#define ERASE 6
#define CLOSE 7

int xtmpl_open(int nx, int ny, const char *title) {
  static int first = 1;
  if (first) {
    int i;
    for (i = 0; i < XTMPL_MAXWINDOWS; i++)
      dpy[i] = NULL;
    first = 0;
  }

  xtmpl_record_event(OPEN, (char *)&nx, (char *)&ny, title, NULL);

  if (xtmpl_view) {
    if (!(WIN(dpy) = XOpenDisplay(NULL))) {
      fprintf(stderr, "xtmpl_open: unable to connect with X server[%s]\n",
              XDisplayName(NULL));
      return XTMPL_ERR;
    }
    WIN(dim_x) = nx;
    WIN(dim_y) = ny;
    WIN(screen) = DefaultScreen(WIN(dpy));
    WIN(window) = XCreateSimpleWindow(
        WIN(dpy), RootWindow(WIN(dpy), WIN(screen)), 0, 0, WIN(dim_x),
        WIN(dim_y), 0, BlackPixel(WIN(dpy), WIN(screen)),
        WhitePixel(WIN(dpy), WIN(screen)));
    XStoreName(WIN(dpy), WIN(window), (title ? title : "Untiteled window"));
    XSelectInput(WIN(dpy), WIN(window), ExposureMask);
    XMapWindow(WIN(dpy), WIN(window));
    XWindowEvent(WIN(dpy), WIN(window), ExposureMask, &xevent);

    WIN(gcv).foreground = BlackPixel(WIN(dpy), DefaultScreen(WIN(dpy)));
    WIN(draw_gc) = XCreateGC(WIN(dpy), WIN(window), GCForeground, &WIN(gcv));
    WIN(gcv).foreground = WhitePixel(WIN(dpy), DefaultScreen(WIN(dpy)));
    WIN(erase_gc) = XCreateGC(WIN(dpy), WIN(window), GCForeground, &WIN(gcv));

    XSetWindowBackground(WIN(dpy), WIN(window),
                         WhitePixel(WIN(dpy), DefaultScreen(WIN(dpy))));
  } else {
    /*
    this is needed for the record-function!!!
    */

    WIN(dim_x) = nx;
    WIN(dim_y) = ny;
  }

  xtmpl_clear();
  xtmpl_line_width(0); /* default line-width */

  return XTMPL_OK;
}
int xtmpl_line_width(int lw) {
  xtmpl_record_event(LINE_W, (char *)&lw, NULL, NULL, NULL);

  if (WIN(dpy) && xtmpl_view) {
    WIN(line_width) = WIN(gcv).line_width = lw;
    XChangeGC(WIN(dpy), WIN(erase_gc), GCLineWidth, &WIN(gcv));
    XChangeGC(WIN(dpy), WIN(draw_gc), GCLineWidth, &WIN(gcv));
  }
  return XTMPL_OK;
}
int xtmpl_clear() {
  xtmpl_record_event(CLEAR, NULL, NULL, NULL, NULL);

  if (WIN(dpy) && xtmpl_view) {
    XClearWindow(WIN(dpy), WIN(window));
    XFlush(WIN(dpy));
    return XTMPL_OK;
  } else
    return XTMPL_ERR;
}
int xtmpl_draw_a_line(int x0, int y0, int x1, int y1, GC gc) {
  if (WIN(dpy) && xtmpl_view) {
    XDrawLine(WIN(dpy), WIN(window), gc, x0, WIN(dim_y) - 1 - y0, x1,
              WIN(dim_y) - 1 - y1);
    XFlush(WIN(dpy));
    return XTMPL_OK;
  } else
    return XTMPL_ERR;
}
int xtmpl_draw_a_text(int x0, int y0, char *str, int nstr, GC gc) {
  XTextItem str_item[1];
  str_item[0].chars = str;
  str_item[0].nchars = nstr;
  str_item[0].delta = 0;
  str_item[0].font = None;

  if (WIN(dpy) && xtmpl_view) {
    XDrawText(WIN(dpy), WIN(window), gc, x0, WIN(dim_y) - 1 - y0, str_item, 1);
    XFlush(WIN(dpy));
    return XTMPL_OK;
  } else
    return XTMPL_ERR;
}
int xtmpl_draw_text(int x0, int y0, char *str, int nstr) {
  return xtmpl_draw_a_text(x0, y0, str, nstr, WIN(draw_gc));
}
int xtmpl_erase_text(int x0, int y0, char *str, int nstr) {
  return xtmpl_draw_a_text(x0, y0, str, nstr, WIN(erase_gc));
}
int xtmpl_text(int fg, int x0, int y0, char *str, int nstr) {
  return xtmpl_draw_a_text(x0, y0, str, nstr,
                           (fg ? WIN(draw_gc) : WIN(erase_gc)));
}
int xtmpl_line(int fg, int x0, int y0, int x1, int y1) {
  /*
  draw a line in Foreground if (fg), otherwise draw with Background.
  (0,0) is lower-left
  (nx-1,ny-1) is upper-right
  */
  xtmpl_record_event((fg ? DRAW : ERASE), (char *)&x0, (char *)&y0, (char *)&x1,
                     (char *)&y1);

  return xtmpl_draw_a_line(x0, y0, x1, y1, (fg ? WIN(draw_gc) : WIN(erase_gc)));
}
int xtmpl_draw_line(int x0, int y0, int x1, int y1) {
  xtmpl_record_event(DRAW, (char *)&x0, (char *)&y0, (char *)&x1, (char *)&y1);

  return xtmpl_draw_a_line(x0, y0, x1, y1, WIN(draw_gc));
}
int xtmpl_erase_line(int x0, int y0, int x1, int y1) {
  xtmpl_record_event(ERASE, (char *)&x0, (char *)&y0, (char *)&x1, (char *)&y1);

  return xtmpl_draw_a_line(x0, y0, x1, y1, WIN(erase_gc));
}
int xtmpl_dot(int x0, int y0, int width, int fg) {
  int w, lw;

  w = (width > 0 ? width / 2 : 0);
  lw = WIN(line_width);

  xtmpl_line_width(2 * w + 1);
  if (w)
    xtmpl_line(fg, x0 - w, y0 - w, x0 + w, y0 + w);
  else
    xtmpl_line(fg, x0, y0, x0 + 1, y0 + 1);
  xtmpl_line_width(lw);

  return 0;
}
int xtmpl_box(int x0, int y0, int width, double yoffset, int fg) {
  int w, lw, yoff;

  w = (width > 0 ? width / 2 : 0);
  lw = WIN(line_width);

  yoff = yoffset * width;

  xtmpl_line_width(2 * w + 1);
  if (w)
    xtmpl_line(fg, x0 - w, y0 + yoff, x0 + w, y0 + yoff);
  else
    xtmpl_line(fg, x0, y0 + yoff, x0 + 1, y0 + yoff);
  xtmpl_line_width(lw);

  return 0;
}
int xtmpl_close() {
  xtmpl_record_event(CLOSE, NULL, NULL, NULL, NULL);

  if (WIN(dpy) && xtmpl_view) {
    XCloseDisplay(WIN(dpy));
    WIN(dpy) = NULL;
    return XTMPL_OK;
  } else
    return XTMPL_ERR;
}
int xtmpl_record(const char *filename, const char *channels) {
  if ((rec_fp = fopen(filename, "w+b"))) {
    if (channels) {
      rec_channels = (char *)malloc(strlen(channels) + 1);
      strcpy(rec_channels, channels);
    }
    return 0;
  } else {
    fprintf(stderr, "xtmpl_record_start: failed to open file[%s]\n", filename);
    return 1;
  }
}
int xtmpl_record_end() {
  if (rec_fp) {
    fclose(rec_fp);
    rec_fp = NULL;
    if (rec_channels) {
      free(rec_channels);
      rec_channels = NULL;
    }
  }
  return 0;
}
int record_channel_q() {
  char w[2];
  static int prev_w = -1, prev_result;

  if (xtmpl_window == prev_w)
    return prev_result;
  else
    prev_w = xtmpl_window;

  if (!rec_fp)
    return (prev_result = 0);

  if (!rec_channels)
    return (prev_result = 1);

  sprintf(w, "%1d", xtmpl_window);
  w[1] = '\0';
  return (prev_result = (strpbrk(w, rec_channels) ? 1 : 0));
}
int big_endian_query() {
  short int i = 1;
  char *ic = (char *)&i;
  return (int)ic[sizeof(i) - 1];
}
int byte_swap(char *what, int size, int n) {
  static char hold[256];
  int i, j, k, pos;
  for (i = pos = 0; i < n; i++, pos += size) {
    memcpy(hold, &what[pos], size);
    for (j = 0, k = size - 1; j < size; j++, k--)
      what[pos + j] = hold[k];
  }
  return XTMPL_OK;
}
int xtmpl_record_event(int type, char *a0, char *a1, const char *a2, char *a3) {
  static int win = 0;
  static int first = 1;
  static int do_byteswap;
  double x[4];
  int i;
  char c[16];
  short int si[6], len, ilen;
  int si_siz = sizeof(short int);

  if (!record_channel_q())
    return XTMPL_OK;

  if (xtmpl_window != win || first) {
    c[0] = NEW_WIN;
    si[0] = win = xtmpl_window;
    fwrite(c, 1, 1, rec_fp);
    fwrite(si, si_siz, 1, rec_fp);
    do_byteswap = DO_BYTESWAP;
    first = 0;
  }

  switch (type) {
  case OPEN:
    c[0] = TITLE;
    len = strlen(a2) + 1;
    fwrite(c, 1, 1, rec_fp);
    ilen = (int)len;
    if (do_byteswap)
      byte_swap((char *)&len, si_siz, 1);
    fwrite(&len, si_siz, 1, rec_fp);
    fwrite(a2, ilen, 1, rec_fp);

    c[0] = OPEN;
    si[0] = *((int *)a0);
    si[1] = *((int *)a1);
    if (do_byteswap)
      byte_swap((char *)si, si_siz, 2);
    fwrite(c, 1, 1, rec_fp);
    fwrite(si, si_siz, 2, rec_fp);

    break;

  case LINE_W:
    c[0] = LINE_W;
    si[0] = *((int *)a0);
    if (do_byteswap)
      byte_swap((char *)si, si_siz, 1);
    fwrite(c, 1, 1, rec_fp);
    fwrite(si, si_siz, 1, rec_fp);

    break;

  case CLEAR:
    c[0] = CLEAR;
    fwrite(c, 1, 1, rec_fp);

    break;

  case DRAW:
  case ERASE:

    x[0] = *((int *)a0) / (double)WIN(dim_x);
    x[1] = *((int *)a1) / (double)WIN(dim_y);
    x[2] = *((int *)a2) / (double)WIN(dim_x);
    x[3] = *((int *)a3) / (double)WIN(dim_y);

    for (i = 0; i < 4; i++)
      si[i] = (short int)(1000.0 * x[i]);

    if (do_byteswap)
      byte_swap((char *)si, si_siz, 4);

    c[0] = (type == DRAW ? DRAW : ERASE);
    fwrite(c, 1, 1, rec_fp);
    fwrite(si, si_siz, 4, rec_fp);

    break;

  case CLOSE:
    c[0] = CLOSE;
    fwrite(c, 1, 1, rec_fp);

    break;

  default:
    break;
  }

  return XTMPL_OK;
}
int press_ret(const char *do_what) {
  char str[4];
  int how_long;
  if (do_what)
    printf("press return to [%s] ", do_what);
  fgets(str, sizeof(str), stdin);
  if (1 == sscanf(str, "%d", &how_long)) {
    char s[10];
    sprintf(s, "sleep %d", how_long);
    system(s);
  }
  return XTMPL_OK;
}
int xtmpl_press_ret(const char *do_what) { return press_ret(do_what); }
int my_fread(void *ptr, size_t size, size_t nitems, FILE *stream) {
  long int position;
  int i;

  if (!forever)
    return fread(ptr, size, nitems, stream);

  position = ftell(stream);
  while ((i = fread(ptr, size, nitems, stream)) != (int)nitems) {
    fseek(stream, position, SEEK_SET);
    system("sleep 1");
  }
  return i;
}
int xtmpl_playback(const char *filename, const char *channels, double magnify) {
  FILE *fp;
  int type;
  int ix[4];
  int eno = 0;
  int sti = 0;
  char *title = NULL, w[2];
  int view = 1;
  char c[16];
  short int si[6], len;
  int si_siz = sizeof(short int);
  int do_byteswap = DO_BYTESWAP;

  if (!strcmp(filename, "-")) {
    fp = stdin;
    sti = 1;
  } else if (!(fp = fopen(filename, "rb"))) {
    fprintf(stderr, "xtmpl_playback: cannot open file [%s]\n", filename);
    return XTMPL_ERR;
  }

  while (my_fread(c, 1, 1, fp)) {
    eno++;
    type = (int)c[0];
    switch (type) {
    case OPEN:
      if (my_fread(si, si_siz, 2, fp) != 2) {
        PERR;
        break;
      }
      if (!view)
        break;
      if (do_byteswap)
        byte_swap((char *)si, si_siz, 2);
      xtmpl_open((int)(magnify * si[0]), (int)(magnify * si[1]), title);
      if (!sti)
        press_ret("continue");
      if (title)
        free(title);
      title = NULL;
      break;

    case TITLE:
      if (my_fread(&len, si_siz, 1, fp) != 1) {
        PERR;
        break;
      }
      if (do_byteswap)
        byte_swap((char *)&len, si_siz, 1);

      title = (char *)malloc((size_t)len);
      if (my_fread((view ? title : (char *)c), (int)len, 1, fp) != 1) {
        PERR;
        break;
      }
      break;

    case NEW_WIN:
      if (my_fread(si, si_siz, 1, fp) != 1) {
        PERR;
        break;
      }
      if (do_byteswap)
        byte_swap((char *)si, si_siz, 1);

      xtmpl_window = (int)si[0];
      if (channels) {
        sprintf(w, "%1d", xtmpl_window);
        w[1] = '\0';
        view = (strpbrk(w, channels) ? 1 : 0);
      } else
        view = 1;

      break;

    case LINE_W:
      if (my_fread(si, si_siz, 1, fp) != 1) {
        PERR;
        break;
      }
      if (!view)
        break;
      if (do_byteswap)
        byte_swap((char *)si, si_siz, 1);

      xtmpl_line_width((int)(magnify * si[0]));
      break;

    case CLEAR:
      if (!view)
        break;
      xtmpl_clear();
      break;

    case CLOSE:
      if (!view)
        break;
      xtmpl_close();
      break;

    case DRAW:
    case ERASE:
      if (my_fread(si, si_siz, 4, fp) != 4) {
        PERR;
        break;
      }
      if (!view)
        break;

      if (do_byteswap)
        byte_swap((char *)si, si_siz, 4);

      ix[0] = si[0] * WIN(dim_x) / 1000.;
      ix[1] = si[1] * WIN(dim_y) / 1000.;
      ix[2] = si[2] * WIN(dim_x) / 1000.;
      ix[3] = si[3] * WIN(dim_y) / 1000.;

      if (type == DRAW)
        xtmpl_draw_line(ix[0], ix[1], ix[2], ix[3]);
      else
        xtmpl_erase_line(ix[0], ix[1], ix[2], ix[3]);

      break;

    default:
      PERR;
      break;
    }
  }

  fclose(fp);
  return XTMPL_OK;
}
#ifdef PLAYBACK
int main(int argc, char **argv) {
  char *channels = NULL;
  double magnify = 1;

  extern char *optarg;
  extern int optind;
  int c;

  while ((c = getopt(argc, argv, "m:c:f")) != EOF)
    switch (c) {
    case 'm':
      magnify = atof(optarg);
      break;

    case 'c':
      channels = (char *)malloc(strlen(optarg) + 1);
      strcpy(channels, optarg);
      break;

    case 'f':
      forever = 1;
      break;

    case '?':
      fprintf(stderr, "Usage: %s [-m mag] [-c channels] [-f] FILE\n", argv[0]);
      exit(1);
      break;
    }

  xtmpl_view = 1;
  return xtmpl_playback((argv[optind] ? argv[optind] : "-"), channels, magnify);
}
#endif
#ifdef TEST
int doit(int op, int cl) {
#define DIM 100
  if (op)
    xtmpl_open(DIM, DIM, "title");
  xtmpl_draw_line(0, 0, DIM - 1, DIM - 1);
  system("sleep 1");
  xtmpl_erase_line(0, 0, DIM - 1, DIM - 1);

  xtmpl_line_width(10);
  xtmpl_draw_line(0, 0, DIM - 1, DIM - 1);
  system("sleep 1");
  xtmpl_erase_line(0, 0, DIM - 1, DIM - 1);

  xtmpl_line_width(0);
  xtmpl_draw_line(0, 0, DIM - 1, DIM - 1);
  system("sleep 1");
  xtmpl_erase_line(0, 0, DIM - 1, DIM - 1);

  if (cl)
    xtmpl_close();
}
main() {
  char a[] = "1234ABCDabcd";
  char b[] = "123ABCabc";
  char c[] = "12ABab";

  if (0) {
    printf("%s %s %s\n", a, b, c);
    byte_swap(a, 4, 3);
    byte_swap(b, 3, 3);
    byte_swap(c, 2, 3);
    printf("%s %s %s\n", a, b, c);
  }

  xtmpl_record("test.movie", NULL);
  doit(1, 0); /* xtmpl_window is default 0 */
  xtmpl_window = 1;
  doit(1, 0);
  xtmpl_window = 0;
  doit(0, 1);
  xtmpl_window = 1;
  doit(0, 1);
}
#endif
#ifdef HILDE
main() {
  int i, j = 10;

  /* xtmpl_view = 0; */
  /* xtmpl_record("test.movie", NULL); */

  xtmpl_open(100, 100, "hilde");
  for (i = 0; i < 10; i++) {
    xtmpl_box(50 + i * 4, 50, 4, (i % 2 ? 0.0 : -0.5), 1);
    system("sleep 1");
  }
}
#endif

#endif
