#include "GMRFLib/GMRFLib.h"
#include "GMRFLib/GMRFLibP.h"

#if defined(__linux__)
#include <sched.h>
void GMRFLib_numa_get(int *cpu, int *numa)
{
	int c = sched_getcpu();
	int n = numa_node_of_cpu(c);
	if (cpu)
		*cpu = c;
	if (numa)
		*numa = n;
}

int GMRFLib_numa(void) 
{
	return ((numa_available() > -1) && (numa_num_configured_nodes() > 1)
		? 1 : 0);
}

int GMRFLib_numa_nodes(void) 
{
	if (GMRFLib_numa()) {
		return (numa_num_configured_nodes());
	} else {
		return 0;
	}
}

#else

void GMRFLib_numa_get(int *cpu, int *numa)
{
	if (cpu)
		*cpu = 0;
	if (numa)
		*numa = 0;
}

int GMRFLib_numa(void) 
{
	return 0;
}

int GMRFLib_numa_nodes(void) 
{
	return 0;
}

#endif
