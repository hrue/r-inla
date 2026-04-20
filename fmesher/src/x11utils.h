#ifndef _FMESH_X11UTILS_
#define _FMESH_X11UTILS_ 1
#ifdef FMESHER_WITH_X

#include <cstddef>
//#include <cstring>
#include <iomanip>
#include <iostream>
#include <list>
#include <map>
#include <set>
#include <string>
#include <vector>

#include "mesh.h"
#include "vector.h"

namespace fmesh {

class Mesh;

class Xtmpl {
private:
  static int window_next_;
  int window_;
  char *name_char_;
  int sx_, sy_;
  double minx_, maxx_, miny_, maxy_;
  bool draw_text_;
  double delay_;

public:
  Xtmpl(const Xtmpl &X);
  Xtmpl(bool draw_text, int sx, int sy, double minx, double maxx, double miny,
        double maxy, std::string name = "fmesher::Mesh");
  void reopen(int sx, int sy);
  void reopen(int sx, int sy, bool draw_text);
  void open(std::string name, int sx, int sy);
  void close();
  ~Xtmpl();

  void clear();

  void delay(double set_delay);
  void delay() const;

  void setSize(int sx, int sy);
  void setAxis(double minx, double maxx, double miny, double maxy);
  double width() const;

  void draw_line(bool fg, double x0, double y0, double x1, double y1);

  void dot(bool fg, const Point &s0, int sz);
  void dot_on_sphere(bool fg, const Point &s0, int sz, double xoffset);
  void arc(bool fg, const Point &s0, const Point &s1, double xoffset);
  void line(bool fg, const Point &s0, const Point &s1);
  void text(bool fg, const Point &s0, std::string str);
};

} /* namespace fmesh */

#endif
#endif
