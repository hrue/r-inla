#include "GMRFLib/GMRFLib.h"
#include "GMRFLib/GMRFLibP.h"

#include <sched.h>
#include "my-numa.h"

#if defined(INLA_WITH_NUMA)

#if !defined(INLA_WITH_HWLOC)
#define INLA_WITH_HWLOC
#endif

#include <numa.h>
#include <numaif.h>

// oops: need to call GMRFLib_numa_init() before use (this is done from main()...)
static int numa_have = -1;
static int numa_nodes = -1;

void GMRFLib_numa_init(void)
{
	if (numa_have < 0)
		GMRFLib_numa_have();
}

int GMRFLib_numa_have(void)
{
	if (numa_have >= 0) {
		return numa_have;
	} else {
		numa_have = ((numa_available() > -1) && (numa_num_configured_nodes() >= 1) ? 1 : 0);
		numa_nodes = IMAX(1, numa_num_configured_nodes());
		return numa_have;
	}
}

void GMRFLib_numa_get(int *cpu, int *numa_node)
{
	unsigned int ucpu, unode; 

	getcpu(&ucpu, &unode); 
	if (cpu) {
		*cpu = (int) ucpu;
	}
	if (numa_node) {
		*numa_node = (int) unode;
	}
}

int GMRFLib_numa_nodes(void)
{
	return numa_nodes;
}

int GMRFLib_numa_node_of_ptr(void *ptr)
{
	int numa_node = -1;
	if (numa_have == 1) {
		get_mempolicy(&numa_node, NULL, 0, (void *) ptr, MPOL_F_NODE | MPOL_F_ADDR);
	}
	return numa_node;
}

int GMRFLib_numa_cache_hitmiss_core(void *ptr, const char *filename, int lineno) 
{
	// return -1 if not in use, 0=hit, 1=miss
	int numa_ptr = GMRFLib_numa_node_of_ptr(ptr);
	if (numa_ptr >= 0) {
		char *nm = NULL;
		GMRFLib_sprintf(&nm, "%s:%1d", filename, lineno);

		if (GMRFLib_trace_cache_hitmiss((const char *) nm)) {
			int numa_cpu;
			GMRFLib_numa_get(NULL, &numa_cpu);
			return (numa_cpu == numa_ptr ? 0 : 1);
		} else {
			return -1;
		}
	}
	return -1;
}

void *GMRFLib_numa_alloc_onnode(size_t size, int node) 
{
	return numa_alloc_onnode(size, node);
}

#else

void *GMRFLib_numa_alloc_onnode(size_t size, int node) 
{
	return malloc(size);
}

void GMRFLib_numa_init(void)
{
	// nothing to do
}

void GMRFLib_numa_get(int *cpu, int *numa_node)
{
	if (cpu) {
		*cpu = sched_getcpu();
	}
	if (numa_node) {
		*numa_node = 0;
	}
}

int GMRFLib_numa_have(void)
{
	return 0;
}

int GMRFLib_numa_nodes(void)
{
	return 1;
}

int GMRFLib_numa_node_of_ptr(void *UNUSED(ptr))
{
	return -1;
}

#endif


#if defined(INLA_WITH_HWLOC)
#include <hwloc.h>

size_t GMRFLib_get_L3_cache(void)
{
	// this returns the first one found, must be checked that this is what we want!!!!!!!!!
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
	// this returns the first one found, must be checked that this is what we want!!!!!!!!
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

size_t GMRFLib_numa_get_L3_cache(int UNUSED(nnode))
{
	return 0;
}

#endif
