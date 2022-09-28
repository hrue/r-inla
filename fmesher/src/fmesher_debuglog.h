#ifndef FMESHER_DEBUGLOG_HH
#define FMESHER_DEBUGLOG_HH

#include <iostream>

// Define NDEBUG to disable assert
#include <cassert>

#ifndef FM_CIN
#define FM_CIN std::cin
#endif
#ifndef FM_COUT
#define FM_COUT std::cout
#endif

#ifndef WHEREAMI
#define WHEREAMI __FILE__ << "(" << __LINE__ << ")\t"
#endif

#ifndef FMLOG_
#define FMLOG_(msg) FM_COUT << WHEREAMI << msg;
#endif

#ifndef FMLOG
#ifdef FMESHER_DEBUG
#define FMLOG(msg) FMLOG_(msg)
#else
#define FMLOG(msg)
#endif
#endif

#ifndef NOT_IMPLEMENTED
#define NOT_IMPLEMENTED (FM_COUT					 \
<< __FILE__ << "(" << __LINE__ << ")\t"	\
<< "NOT IMPLEMENTED: "				              \
<< __PRETTY_FUNCTION__ << std::endl);
#endif


#endif
