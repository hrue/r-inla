#ifdef FMESHER_WITH_X

#include <cerrno>
#include <cmath>
#include <cstddef>
#include <cstring>
#include <ctime>
#include <map>
#include <set>
#include <sstream>

#include "mesh.h"
#include "x11utils.h"
#include "xtmpl.h"

#include "fmesher_debuglog.h"

using std::endl;

namespace fmesh {

int Xtmpl::window_next_ = 0;

Xtmpl::Xtmpl(const Xtmpl &X)
    : window_(X.window_ + 1), name_char_(NULL), sx_(X.sx_), sy_(X.sy_),
      minx_(X.minx_), maxx_(X.maxx_), miny_(X.miny_), maxy_(X.maxy_),
      draw_text_(true), delay_(0.0) {
  open(std::string(X.name_char_), X.sx_, X.sy_);
  setAxis(X.minx_, X.maxx_, X.miny_, X.maxy_);
};
Xtmpl::Xtmpl(bool draw_text, int sx, int sy, double minx, double maxx,
             double miny, double maxy, std::string name)
    : window_(-1), name_char_(NULL), sx_(sx), sy_(sy), minx_(minx), maxx_(maxx),
      miny_(miny), maxy_(maxy), draw_text_(draw_text), delay_(0.0) {
  open(name, sx_, sy_);
  setAxis(minx, maxx, miny, maxy);
};
void Xtmpl::reopen(int sx, int sy) {
  if (!(window_ < 0))
    close();
  else
    window_ = window_next_++;
  sx_ = sx;
  sy_ = sy;
  xtmpl_window = window_;
  xtmpl_open(sx_, sy_, name_char_);
};
void Xtmpl::reopen(int sx, int sy, bool draw_text) {
  reopen(sx, sy);
  draw_text_ = draw_text;
};
void Xtmpl::open(std::string name, int sx, int sy) {
  if (!(window_ < 0))
    close();
  else
    window_ = window_next_++;
  sx_ = sx;
  sy_ = sy;
  if (name_char_)
    delete[] name_char_;
  name_char_ = new char[name.length() + 1];
  name.copy(name_char_, name.length(), 0);
  name_char_[name.length()] = '\0';
  xtmpl_window = window_;
  xtmpl_open(sx, sy, name_char_);
  setAxis(-0.05, 1.05, -0.05, 1.05);
};
void Xtmpl::close() {
  if (window_ < 0)
    return;
  xtmpl_window = window_;
  xtmpl_close();
  window_ = -1;
};
Xtmpl::~Xtmpl() {
  close();
  if (name_char_)
    delete[] name_char_;
};

void Xtmpl::clear() {
  xtmpl_window = window_;
  xtmpl_clear();
}

void Xtmpl::delay(double set_delay) { delay_ = set_delay; };
void Xtmpl::delay() const {
  if (delay_ > 0.0) {
    struct timespec req;
    struct timespec rem;
    req.tv_sec = time_t(delay_);
    req.tv_nsec = long((delay_ - double(req.tv_sec)) * 1.e9);
    while ((::nanosleep(&req, &rem) != -1) && (errno == EINTR)) {
      req = rem;
    }
  }
};

void Xtmpl::setSize(int sx, int sy) { reopen(sx, sy); };
void Xtmpl::setAxis(double minx, double maxx, double miny, double maxy) {
  clear();
  minx_ = minx;
  maxx_ = maxx;
  miny_ = miny;
  maxy_ = maxy;
};
double Xtmpl::width() const { return (maxx_ - minx_); };

void Xtmpl::dot(bool fg, const Point &s0, int sz) {
  xtmpl_window = window_;
  double x0 = (s0[0] - minx_) / (maxx_ - minx_);
  double y0 = (s0[1] - miny_) / (maxy_ - miny_);
  if ((0.0 <= x0) && (x0 <= 1.0) && (0.0 <= y0) && (y0 <= 1.0))
    xtmpl_dot((int)(sx_ * x0), (int)(sy_ * y0), sz, (int)fg);
};

void Xtmpl::dot_on_sphere(bool fg, const Point &s0, int sz, double xoffset) {
  xtmpl_window = window_;
  double x0 = (s0[0] + xoffset - minx_) / (maxx_ - minx_);
  double y0 = (s0[1] - miny_) / (maxy_ - miny_);
  if ((0.0 <= x0) && (x0 <= 1.0) && (0.0 <= y0) && (y0 <= 1.0))
    xtmpl_dot((int)(sx_ * x0), (int)(sy_ * y0), sz, (int)fg);
};

void Xtmpl::draw_line(bool fg, double x0, double y0, double x1, double y1) {
  if ((((x0 <= 0.0) && (x1 <= 0.0)) || ((x0 >= 1.0) && (x1 >= 1.0))) ||
      (((y0 <= 0.0) && (y1 <= 0.0)) || ((y0 >= 1.0) && (y1 >= 1.0))))
    return;
  if (x1 < x0) {
    /* Swap */
    double x = x1;
    double y = y1;
    x1 = x0;
    y1 = y0;
    x0 = x;
    y0 = y;
  }
  if (x0 < x1) {
    /* Truncate at x=0 and 1 */
    if (x0 < 0.0) {
      /* x = (1-a)*x0+a*x1 = 0
         y = (1-a)*y0+a*y1
         a = -x0/(x1-x0), 1-a = x1/(x1-x0)
         y = (x1*y0-x0*y1)/(x1-x0)
      */
      double x = 0.0;
      double y = (x1 * y0 - x0 * y1) / (x1 - x0);
      x0 = x;
      y0 = y;
    }
    if (x1 > 1.0) {
      /* x = (1-a)*x0+a*x1 = 1
         y = (1-a)*y0+a*y1
         a = (1-x0)/(x1-x0), 1-a = (x1-1)/(x1-x0)
         y = ((x1-1)*y0-(x0-1)*y1)/(x1-x0)
      */
      double x = 1.0;
      double y = ((x1 - 1) * y0 - (x0 - 1) * y1) / (x1 - x0);
      x1 = x;
      y1 = y;
    }
  }

  if (((y0 <= 0.0) && (y1 <= 0.0)) || ((y0 >= 1.0) && (y1 >= 1.0)))
    return;
  if (y1 < y0) {
    /* Swap */
    double x = x1;
    double y = y1;
    x1 = x0;
    y1 = y0;
    x0 = x;
    y0 = y;
  }
  if (y0 < y1) {
    /* Truncate at y=0 and 1 */
    if (y0 < 0.0) {
      double x = (y1 * x0 - y0 * x1) / (y1 - y0);
      double y = 0.0;
      x0 = x;
      y0 = y;
    }
    if (y1 > 1.0) {
      double x = ((y1 - 1) * x0 - (y0 - 1) * x1) / (y1 - y0);
      double y = 1.0;
      x1 = x;
      y1 = y;
    }
  }

  xtmpl_window = window_;
  if (fg)
    xtmpl_draw_line((int)(sx_ * x0), (int)(sy_ * y0), (int)(sx_ * x1),
                    (int)(sy_ * y1));
  else
    xtmpl_erase_line((int)(sx_ * x0), (int)(sy_ * y0), (int)(sx_ * x1),
                     (int)(sy_ * y1));
}

void Xtmpl::arc(bool fg, const Point &s0, const Point &s1, double xoffset) {
  int n = 8;
  xtmpl_window = window_;
  double p0[2];
  double p1[2];
  double s[3];
  double l;
  int dim;
  p1[0] = s0[0];
  p1[1] = s0[1];
  for (int i = 1; i <= n; i++) {
    l = 0.0;
    p0[0] = p1[0];
    p0[1] = p1[1];
    for (dim = 0; dim < 3; dim++) {
      s[dim] = ((n - i) * s0[dim] + i * s1[dim]) / n;
      l += s[dim] * s[dim];
    }
    l = std::sqrt(l);
    for (dim = 0; dim < 2; dim++)
      p1[dim] = s[dim] / l;

    draw_line(fg, (p0[0] + xoffset - minx_) / (maxx_ - minx_),
              (p0[1] - miny_) / (maxy_ - miny_),
              (p1[0] + xoffset - minx_) / (maxx_ - minx_),
              (p1[1] - miny_) / (maxy_ - miny_));
  }
}

void Xtmpl::line(bool fg, const Point &s0, const Point &s1) {
  xtmpl_window = window_;
  draw_line(
      fg, (s0[0] - minx_) / (maxx_ - minx_), (s0[1] - miny_) / (maxy_ - miny_),
      (s1[0] - minx_) / (maxx_ - minx_), (s1[1] - miny_) / (maxy_ - miny_));
};

void Xtmpl::text(bool fg, const Point &s0, std::string str) {
  if (!draw_text_)
    return;
  char *str_ = new char[str.length() + 1];
  str.copy(str_, str.length(), 0);
  str_[str.length()] = '\0';
  xtmpl_window = window_;
  double x0 = (s0[0] - minx_) / (maxx_ - minx_);
  double y0 = (s0[1] - miny_) / (maxy_ - miny_);
  if ((0.0 <= x0) && (x0 <= 1.0) && (0.0 <= y0) && (y0 <= 1.0)) {
    if (fg)
      xtmpl_draw_text((int)(sx_ * x0), (int)(sy_ * y0), str_, str.length());
    else
      xtmpl_erase_text((int)(sx_ * x0), (int)(sy_ * y0), str_, str.length());
  }
  delete[] str_;
};

} /* namespace fmesh */

#endif
