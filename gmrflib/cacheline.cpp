#include <new>

using namespace std;
extern "C" int GMRFLib_get_cachelinesize(void);

int GMRFLib_get_cachelinesize(void)
{
#ifdef __cpp_lib_hardware_interference_size
	return (std::hardware_destructive_interference_size);
#else
	return 64; // best guess
#endif	
}
