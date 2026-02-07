#include <strings.h>
#include "GMRFLib/GMRFLib.h"

#if defined(INLA_WITH_STILES)
#       include "stiles.h"
#else

//...

double *sTiles_get_selinv_row(int group_index, int call_index, int node, int *node_neighbors, int size, void **obj);
double sTiles_GetGroupMemoryUsage(int group_ID);
double sTiles_GetGroupsMemoryUsage();
double sTiles_get_logdet(int group_index, int call_index, void **obj);
double sTiles_get_selinv_elm(int group_index, int call_index, int irow, int icol, void **obj);
int *sTiles_return_iperm_vec(int group_index, void **obj);
int *sTiles_return_perm_vec(int group_index, void **obj);
int sTiles_assign_graph(int group_index, void **stile, int N, int NNZ, int *row_indices, int *col_indices);
int sTiles_assign_values(int group_index, int call_index, void **obj, double *x);
int sTiles_bind(int group_index, int call_index, void **obj);
int sTiles_chol(int group_index, int call_index, void **obj);
int sTiles_clear_selinv(int group_index, int call_index, void **obj);
int sTiles_create(void **obj, int num_call_groups, const int *calls_per_group, const int *cores_per_group,
		  const int *factor_type_per_group, const bool *get_inverse);
int sTiles_create_expert(void **obj, int num_call_groups, const int *calls_per_group, const int *cores_per_group,
			 const int *factor_type_per_group, const bool *get_inverse, const int *rhs,
			 const int *arrowhead_size, const int *arrowhead_size_per_group, const int *user_params);
int sTiles_init(void **obj);
int sTiles_init_group(int group_index, void **obj);
int sTiles_return_tile_size();
int sTiles_selinv(int group_index, int call_index, void **obj);
int sTiles_solve_L(int group_index, int call_index, void **obj, double *B, int nrhs);
int sTiles_solve_LLT(int group_index, int call_index, void **obj, double *B, int nrhs);
int sTiles_solve_LT(int group_index, int call_index, void **obj, double *B, int nrhs);
int sTiles_unbind(int group_index, int call_index, void **obj);
void sTiles_expert_user(void);
void sTiles_freeGroup(int group_ID);
void sTiles_map_group_call_to_group_call(void **obj, int group_index1, int call_index1, int group_index2, int call_index2);
void sTiles_quit();
void sTiles_set_correction_mode(int correction_mode);
void sTiles_set_ordering_mode(int reordering);
void sTiles_set_tile_ordering_min_dim(int min_dim);
void sTiles_set_tile_ordering_mode(int reordering);
void sTiles_set_tile_ordering_size(int tile_size);
void sTiles_set_tile_size(int tile_size);
void sTiles_set_tile_type_mode(int tile_type);

#endif
