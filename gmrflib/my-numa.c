#include "GMRFLib/GMRFLib.h"
#include "GMRFLib/GMRFLibP.h"

#if defined(__linux__)
#include <sched.h>
void GMRFLib_numa_get(int *cpu, int *numa) 
{
	int c = sched_getcpu();
	int n = numa_node_of_cpu(c);
	if (cpu) *cpu = c;
	if (numa) *numa = n;
}
#else
void GMRFLib_numa_get(int *cpu, int *numa) 
{
	if (cpu) *cpu = 0;
	if (numa) *numa = 0;
}
#endif
