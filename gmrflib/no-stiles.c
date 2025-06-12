#include <omp.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <strings.h>

#if !defined(INLA_WITH_STILES)

#include "no-stiles.h"

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-parameter"

#define EMPTY_FUNCTION								\
	{								\
		fprintf(stderr, "\n\n\n *** sTiles is not available in this build, abort()\n\n\n"); \
		abort();						\
	}								\

int sTiles_assign_graph(int group_index, sTiles_object *stile, int N, int NNZ, int *row_indices, int *col_indices) EMPTY_FUNCTION;
int sTiles_init(sTiles_object **obj) EMPTY_FUNCTION;
int sTiles_init_group(int group_index, sTiles_object **obj) EMPTY_FUNCTION;
int sTiles_assign_values(int group_index, int call_index, sTiles_object **obj, double *x) EMPTY_FUNCTION;
int sTiles_chol(int group_index, int call_index, sTiles_object **obj) EMPTY_FUNCTION;
int sTiles_selinv(int group_index, int call_index, sTiles_object **obj) EMPTY_FUNCTION;
double sTiles_get_selinv_elm(int group_index, int call_index, int irow, int icol, sTiles_object **obj) EMPTY_FUNCTION;
double *sTiles_get_selinv_row(int group_index, int call_index, int node, int *node_neighbors, int size, sTiles_object **obj) EMPTY_FUNCTION;
double sTiles_get_logdet(int group_index, int call_index, sTiles_object **obj) EMPTY_FUNCTION;
int sTiles_solve_LLT(int group_index, int call_index, sTiles_object **obj, double *B, int sindex) EMPTY_FUNCTION;
int sTiles_solve_L(int group_index, int call_index, sTiles_object **obj, double *B, int sindex) EMPTY_FUNCTION;
int sTiles_solve_LT(int group_index, int call_index, sTiles_object **obj, double *B, int sindex) EMPTY_FUNCTION;
int *sTiles_return_perm_vec(int group_index, sTiles_object **obj) EMPTY_FUNCTION;
int *sTiles_return_iperm_vec(int group_index, sTiles_object **obj) EMPTY_FUNCTION;
int sTiles_clear_selinv(int group_index, int call_index, sTiles_object **obj) EMPTY_FUNCTION;
double sTiles_GetGroupMemoryUsage(int group_ID) EMPTY_FUNCTION;
double sTiles_GetGroupsMemoryUsage() EMPTY_FUNCTION;
void sTiles_freeGroup(int group_ID) EMPTY_FUNCTION;
void sTiles_quit() EMPTY_FUNCTION;
int sTiles_create(sTiles_object **obj, int num_call_groups, const int *calls_per_group, const int *cores_per_group,
		  const int *factor_type_per_group, const bool *get_inverse, const int *rhs) EMPTY_FUNCTION;
int sTiles_create_expert(sTiles_object **obj, int num_call_groups, const int *calls_per_group,
			 const int *cores_per_group, const int *factor_type_per_group,
			 const bool *get_inverse, const int *rhs, const int *arrowhead_size,
			 const int *arrowhead_size_per_group, const int *user_params) EMPTY_FUNCTION;
void sTiles_map_group_call_to_group_call(sTiles_object **obj, int group_index1, int call_index1, int group_index2, int call_index2) EMPTY_FUNCTION;
int sTiles_bind(int group_index, int call_index, sTiles_object **obj) EMPTY_FUNCTION;
int sTiles_unbind(int group_index, int call_index, sTiles_object **obj) EMPTY_FUNCTION;

// exceptions
void sTiles_set_tile_size(int tile_size)
{
}
int sTiles_return_tile_size()
{
	return 40;
}

#undef EMPTY_FUNCTION
#pragma GCC diagnostic pop

#endif
