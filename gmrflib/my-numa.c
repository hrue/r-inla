#include "GMRFLib/GMRFLib.h"
#include "GMRFLib/GMRFLibP.h"

// >>>>>>>> LINUX ONLY <<<<<<<<

// oops: need to call GMRFLib_numa_init() before use (this is done from main()...)
static int POSSIBLY_UNUSED(NUMA_have) = -1;		       // we have (=1) NUMA support or not (=0)
static int POSSIBLY_UNUSED(NUMA_nodes) = -1;		       // number of NUMA nodes. =1 if no NUMA */
static int NUMA_enable = 1;				       // if not enabled, then all NUMA support is disabled (and we return
							       // to the behaviour as if INLA_WITH_NUMA was not defined)
#include "my-numa.h"

#if defined(INLA_WITH_NUMA)

#if !defined(__linux__)
#error NUMA is for Linux only, please compile without 'INLA_WITH_NUMA'
#endif

#include <sched.h>

// older linux (Rocky8 and older)
#if defined(GMRFLib_GETCPU_FIX)
int getcpu(unsigned int *cpu, unsigned int *node)
{
	unsigned int ucpu = sched_getcpu();
	unsigned int unode = (unsigned int) numa_node_of_cpu((int) ucpu);
	if (cpu) {
		*cpu = ucpu;
	}
	if (node) {
		*node = unode;
	}
	return 0;
}
#endif



#if !defined(INLA_WITH_HWLOC)
#define INLA_WITH_HWLOC
#endif

#include <numa.h>
#include <numaif.h>

void GMRFLib_numa_init(void)
{
	if (NUMA_have < 0)
		GMRFLib_numa_have();
}

void GMRFLib_numa_set_ctl(int enable)
{
	NUMA_enable = enable;
}

int GMRFLib_numa_have(void)
{
	if (NUMA_enable) {
		if (NUMA_have >= 0) {
			return NUMA_have;
		} else {
			NUMA_have = ((numa_available() > -1) && (numa_num_configured_nodes() >= 1) ? 1 : 0);
			NUMA_nodes = IMAX(1, numa_num_configured_nodes());
			return NUMA_have;
		}
	} else {
		return 0;
	}
}

int GMRFLib_numa_get_node(void)
{
	int nnode = 0;
	GMRFLib_numa_get(NULL, &nnode);
	return nnode;
}

void GMRFLib_numa_get(int *cpu, int *numa_node)
{
	unsigned int unode;
	if (cpu) {
		unsigned int ucpu;
		getcpu(&ucpu, &unode);
		*cpu = (int) ucpu;
	} else {
		getcpu(NULL, &unode);
	}
		
	if (numa_node) {
		if (NUMA_enable) {
			*numa_node = (int) unode;
		} else {
			*numa_node = 0;
		}
	}
}

int GMRFLib_numa_nodes(void)
{
	return (NUMA_enable ? NUMA_nodes : 1);
}

int GMRFLib_numa_node_of_ptr(void *ptr)
{
	if (NUMA_enable) {
		int numa_node = -1;
		if (NUMA_have == 1) {
			get_mempolicy(&numa_node, NULL, 0, (void *) ptr, MPOL_F_NODE | MPOL_F_ADDR);
		}
		return numa_node;
	} else {
		return 0;
	}
}

int GMRFLib_numa_cache_hitmiss_core(void *ptr, int numa, const char *filename, int lineno)
{
	// return -1 if not in use, 0=hit, 1=miss
	if (NUMA_enable) {
		int numa_ptr = GMRFLib_numa_node_of_ptr(ptr);
		if (numa_ptr >= 0) {
			char *nm = NULL;
			GMRFLib_sprintf(&nm, "%s:%1d", filename, lineno);
			if (GMRFLib_trace_cache_hitmiss((const char *) nm)) {
				return (numa == numa_ptr ? 0 : 1);
			} else {
				return -1;
			}
		}
		return -1;
	} else {
		return -1;
	}
}

void *GMRFLib_numa_alloc_onnode(size_t size, int node)
{
	if (NUMA_enable) {
		return (size > 0 ? numa_alloc_onnode(size, node) : NULL);
	} else {
		return (size > 0 ? malloc(size) : NULL);
	}
}

void GMRFLib_numa_free(void *start, size_t size)
{
	if (NUMA_enable) {
		if (size > 0 && start) {
			numa_free(start, size);
		}
	} else {
		if (size > 0) {
			Free(start);
		}
	}
}

#else

void GMRFLib_numa_set_ctl(int UNUSED(enable))
{
	NUMA_enable = 0;
}


int GMRFLib_numa_get_node(void)
{
	int nnode = 0;
	GMRFLib_numa_get(NULL, &nnode);
	return nnode;
}

void GMRFLib_numa_free(void *start, size_t size)
{
	if (size > 0) {
		Free(start);
	}
}

void *GMRFLib_numa_alloc_onnode(size_t size, int UNUSED(node))
{
	return (size > 0 ? malloc(size) : NULL);
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

int GMRFLib_numa_cache_hitmiss_core(void *UNUSED(ptr), int UNUSED(numa), const char *UNUSED(filename), int UNUSED(lineno))
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
