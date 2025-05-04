#include "GMRFLib/GMRFLib.h"
#include "GMRFLib/GMRFLibP.h"

#if defined(INLA_WITH_NUMA)

#if !defined(INLA_WITH_HWLOC)
#define INLA_WITH_HWLOC
#endif

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

#if defined(INLA_WITH_HWLOC)
#include <hwloc.h>

size_t GMRFLib_get_L3_cache(void)
{
	// this returns the first one found, must be checked that this is what we want
	hwloc_topology_t topology;
	size_t l3 = 0;

	// Initialize and load the topology
	hwloc_topology_init(&topology);
	hwloc_topology_load(topology);

	int depth = hwloc_topology_get_depth(topology);
	for (int i = 0; i < depth; i++) {
		int num_objs = hwloc_get_nbobjs_by_depth(topology, i);
		for (int j = 0; j < num_objs && l3 == 0; j++) {
			hwloc_obj_t obj = hwloc_get_obj_by_depth(topology, i, j);
			if (obj->type == HWLOC_OBJ_L3CACHE &&
			    (obj->attr->cache.type == HWLOC_OBJ_CACHE_UNIFIED || obj->attr->cache.type == HWLOC_OBJ_CACHE_DATA)) {
				l3 = (size_t) (obj->attr->cache.size / ISQR(1024));
			}
		}
	}

	hwloc_topology_destroy(topology);
	return l3;
}

size_t GMRFLib_numa_get_L3_cache(int nnode)
{
	// this returns the first one found, must be checked that this is what we want
	hwloc_topology_t topology;
	hwloc_obj_t numa_node;
	hwloc_obj_t obj = NULL;
	size_t l3 = 0;

	hwloc_topology_init(&topology);
	hwloc_topology_load(topology);

	if ((numa_node = hwloc_get_obj_by_type(topology, HWLOC_OBJ_NUMANODE, nnode))) {
		while (l3 == 0 && (obj = hwloc_get_next_obj_by_type(topology, HWLOC_OBJ_L3CACHE, obj)) != NULL) {
			if (obj->type == HWLOC_OBJ_L3CACHE &&
			    (obj->attr->cache.type == HWLOC_OBJ_CACHE_UNIFIED || obj->attr->cache.type == HWLOC_OBJ_CACHE_DATA)) {
				// Check if the cache is in the subtree of NUMA node 0
				if (hwloc_obj_is_in_subtree(topology, obj, numa_node)) {
					l3 = (size_t) (obj->attr->cache.size / ISQR(1024));
				}
			}
		}
	}
	hwloc_topology_destroy(topology);
	return l3;
}

#else

size_t GMRFLib_get_L3_cache(void)
{
	return 0;
}

size_t GMRFLib_numa_L3_cache(int nnode)
{
	return 0;
}

#endif
