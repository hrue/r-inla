#ifndef FMESHER_NO_X

#include <cstddef>
#include <cstring>
#include <set>
#include <map>
#include <sstream>
#include <cmath>
#include <ctime>
#include <cerrno>

#include "xtmpl.h"
#include "mesh.hh"
#include "x11utils.hh"

#define WHEREAMI __FILE__ << "(" << __LINE__ << ")\t"

#ifdef DEBUG
#define MESH_LOG(msg) std::cout << WHEREAMI << msg;
#else
#define MESH_LOG(msg)
#endif


using std::cout;
using std::endl;

namespace fmesh {


  int Xtmpl::window_next_ = 0;


  Xtmpl::Xtmpl(const Xtmpl& X)
    : window_(X.window_+1), name_char_(NULL),
      sx_(X.sx_), sy_(X.sy_),
      minx_(X.minx_), maxx_(X.maxx_),
      miny_(X.miny_), maxy_(X.maxy_), draw_text_(true),
      delay_(0.0) {
    open(std::string(X.name_char_),X.sx_,X.sy_);
    setAxis(X.minx_, X.maxx_, X.miny_, X.maxy_);
  };
  Xtmpl::Xtmpl(bool draw_text, int sx, int sy,
	       double minx,
	       double maxx,
	       double miny,
	       double maxy,
	       std::string name)
    : window_(-1), name_char_(NULL),
      sx_(sx), sy_(sy),
      minx_(minx), maxx_(maxx),
      miny_(miny), maxy_(maxy), draw_text_(draw_text),
      delay_(0.0) {
    open(name,sx_,sy_);
    setAxis(minx, maxx, miny, maxy);
  };
  void Xtmpl::reopen(int sx, int sy) {
    if (!(window_<0))
      close();
    else
      window_ = window_next_++;
    sx_ = sx;
    sy_ = sy;
    xtmpl_window = window_;
    xtmpl_open(sx_,sy_,name_char_);
  };
  void Xtmpl::reopen(int sx, int sy, bool draw_text) {
    reopen(sx,sy);
    draw_text_ = draw_text;
  };
  void Xtmpl::open(std::string name,
		   int sx, int sy) {
    if (!(window_<0))
      close();
    else
      window_ = window_next_++;
    sx_ = sx;
    sy_ = sy;
    if (name_char_) delete[] name_char_;
    name_char_ = new char[name.length()+1];
    name.copy(name_char_,name.length(),0);
    name_char_[name.length()] = '\0';
    xtmpl_window = window_;
    xtmpl_open(sx,sy,name_char_);
    setAxis(-0.05,1.05,-0.05,1.05);
  };
  void Xtmpl::close() {
    if (window_<0)
      return;
    xtmpl_window = window_;
    xtmpl_close();
    window_ = -1;
  };
  Xtmpl::~Xtmpl() {
    close();
    if (name_char_) delete[] name_char_;
  };
  
  void Xtmpl::clear() {
    xtmpl_window = window_;
    xtmpl_clear();
  }
  
  void Xtmpl::delay(double set_delay) { delay_ = set_delay; };
  void Xtmpl::delay() const {
    if (delay_>0.0) {
      struct timespec req;
      struct timespec rem;
      req.tv_sec = time_t(delay_);
      req.tv_nsec = long((delay_-double(req.tv_sec))*1.e9);
      while ((::nanosleep(&req,&rem) != -1) &&
	     (errno == EINTR)) {
	req = rem;
      }
    }
  };
  

  void Xtmpl::setSize(int sx, int sy) {
    reopen(sx,sy);
  };
  void Xtmpl::setAxis(double minx, double maxx,
		      double miny, double maxy) {
    clear();
    minx_ = minx;
    maxx_ = maxx;
    miny_ = miny;
    maxy_ = maxy;
  };
  double Xtmpl::width() const { return (maxx_-minx_); };

  
  
  
  void Xtmpl::dot(bool fg, const Point& s0, int sz)
  {
    xtmpl_window = window_;
    xtmpl_dot((int)(sx_*(s0[0]-minx_)/(maxx_-minx_)),
	      (int)(sy_*(s0[1]-miny_)/(maxy_-miny_)),
	      sz,
	      (int)fg);
  };

  void Xtmpl::dot_on_sphere(bool fg, const Point& s0, int sz, double xoffset)
  {
    xtmpl_window = window_;
    xtmpl_dot((int)(sx_*(s0[0]+xoffset-minx_)/(maxx_-minx_)),
	      (int)(sy_*(s0[1]-miny_)/(maxy_-miny_)),
	      sz,
	      (int)fg);
  };

  void Xtmpl::arc(bool fg, const Point& s0, const Point& s1, double xoffset)
  {
    int n = 8;
    xtmpl_window = window_;
    double p0[2];
    double p1[2];
    double s[3];
    double l;
    int dim;
    p1[0] = s0[0];
    p1[1] = s0[1];
    for (int i=1;i<=n;i++) {
      l = 0.0;
      p0[0] = p1[0]; p0[1] = p1[1];
      for (dim=0;dim<3;dim++) {
	s[dim] = ((n-i)*s0[dim]+i*s1[dim])/n;
	l += s[dim]*s[dim];
      }
      l = std::sqrt(l);
      for (dim=0;dim<2;dim++)
	p1[dim] = s[dim]/l;
      
      if (fg)
	xtmpl_draw_line((int)(sx_*(p0[0]+xoffset-minx_)/(maxx_-minx_)),
			(int)(sy_*(p0[1]-miny_)/(maxy_-miny_)),
			(int)(sx_*(p1[0]+xoffset-minx_)/(maxx_-minx_)),
			(int)(sy_*(p1[1]-miny_)/(maxy_-miny_)));
      else
	xtmpl_erase_line((int)(sx_*(p0[0]+xoffset-minx_)/(maxx_-minx_)),
			 (int)(sy_*(p0[1]-miny_)/(maxy_-miny_)),
			 (int)(sx_*(p1[0]+xoffset-minx_)/(maxx_-minx_)),
			 (int)(sy_*(p1[1]-miny_)/(maxy_-miny_)));
    }
  };

  void Xtmpl::line(bool fg, const Point& s0, const Point& s1)
  {
    xtmpl_window = window_;
    if (fg)
      xtmpl_draw_line((int)(sx_*(s0[0]-minx_)/(maxx_-minx_)),
		      (int)(sy_*(s0[1]-miny_)/(maxy_-miny_)),
		      (int)(sx_*(s1[0]-minx_)/(maxx_-minx_)),
		      (int)(sy_*(s1[1]-miny_)/(maxy_-miny_)));
    else
      xtmpl_erase_line((int)(sx_*(s0[0]-minx_)/(maxx_-minx_)),
		       (int)(sy_*(s0[1]-miny_)/(maxy_-miny_)),
		       (int)(sx_*(s1[0]-minx_)/(maxx_-minx_)),
		       (int)(sy_*(s1[1]-miny_)/(maxy_-miny_)));
  };

  void Xtmpl::text(bool fg, const Point& s0, std::string str)
  {
    if (!draw_text_) return;
    char* str_ = new char[str.length()+1];
    str.copy(str_,str.length(),0);
    str_[str.length()] = '\0';
    xtmpl_window = window_;
    if (fg)
      xtmpl_draw_text((int)(sx_*(s0[0]-minx_)/(maxx_-minx_)),
		      (int)(sy_*(s0[1]-miny_)/(maxy_-miny_)),
		      str_,str.length());
    else
      xtmpl_erase_text((int)(sx_*(s0[0]-minx_)/(maxx_-minx_)),
		       (int)(sy_*(s0[1]-miny_)/(maxy_-miny_)),
		       str_,str.length());
    delete[] str_;
  };





} /* namespace fmesh */

#endif
