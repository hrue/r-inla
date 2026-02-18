#include <omp.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <strings.h>

#if !defined(INLA_WITH_STILES)
#       include "no-stiles.h"

#       pragma GCC diagnostic push
#       pragma GCC diagnostic ignored "-Wunused-parameter"

#       ifdef __cplusplus
extern "C" {
#       endif

#       define EMPTY_FUNCTION {						\
		fprintf(stderr, "\n\n\n *** sTiles is not available in this build, abort()\n\n\n"); \
		abort();						\
	}

	void sTiles_set_tile_size(int) {
		// empty
	}
	// ... int sTiles_return_tile_size(void) { return 40;

	const char *sTiles_get_version(void) {
		return (const char *) NULL;
	}
	double sTiles_GetGroupMemoryUsage(int) EMPTY_FUNCTION;
	double sTiles_GetGroupsMemoryUsage(void) EMPTY_FUNCTION;
	double sTiles_debug_matrix(int, int, void **) EMPTY_FUNCTION;
	double sTiles_get_chol_timing(int, int, void **) EMPTY_FUNCTION;
	double sTiles_get_logdet(int, int, void **) EMPTY_FUNCTION;
	double sTiles_get_selinv_elm(int, int, int, int, void **) EMPTY_FUNCTION;
	double sTiles_get_selinv_timing(int, int, void **) EMPTY_FUNCTION;
	int *sTiles_return_iperm_vec(int, void **) EMPTY_FUNCTION;
	int *sTiles_return_perm_vec(int, void **) EMPTY_FUNCTION;
	int sTiles_assign_graph(int, void **, int, int, int *, int *) EMPTY_FUNCTION;
	int sTiles_assign_graph_one_call(int, int, void **, int, int, int *, int *) EMPTY_FUNCTION;
	int sTiles_assign_values(int, int, void **, double *) EMPTY_FUNCTION;
	int sTiles_bind(int, int, void **) EMPTY_FUNCTION;
	int sTiles_chol(int, int, void **) EMPTY_FUNCTION;
	int sTiles_clear_selinv(int, int, void **) EMPTY_FUNCTION;
	int sTiles_create(void **, int, const int *, const int *, const int *, const bool *) EMPTY_FUNCTION;
	int sTiles_create_expert(void **, int, const int *, const int *, const int *, const bool *, const int *, const int *, const int *, const int *) EMPTY_FUNCTION;
	int sTiles_get_auto_tile_size(int) EMPTY_FUNCTION;
	int sTiles_get_num_calls(void *, int) EMPTY_FUNCTION;
	int sTiles_init(void **) EMPTY_FUNCTION;
	int sTiles_init_group(int, void **) EMPTY_FUNCTION;
	int sTiles_selinv(int, int, void **) EMPTY_FUNCTION;
	int sTiles_solve_L(int, int, void **, double *, int) EMPTY_FUNCTION;
	int sTiles_solve_LLT(int, int, void **, double *, int) EMPTY_FUNCTION;
	int sTiles_solve_LT(int, int, void **, double *, int) EMPTY_FUNCTION;
	int sTiles_unbind(int, int, void **) EMPTY_FUNCTION;
	void sTiles_expert_user() EMPTY_FUNCTION;
	void sTiles_freeGroup(int) EMPTY_FUNCTION;
	void sTiles_map_group_call_to_group_call(void **, int, int, int, int) EMPTY_FUNCTION;
	void sTiles_print_chol_timings(int, void **) EMPTY_FUNCTION;
	void sTiles_print_logdets(int, void **) EMPTY_FUNCTION;
	void sTiles_print_selinv_timings(int, void **) EMPTY_FUNCTION;
	void sTiles_print_version(void) EMPTY_FUNCTION;
	void sTiles_quit(void) EMPTY_FUNCTION;
	void sTiles_set_log_level(int) EMPTY_FUNCTION;

	void sTiles_set_ordering_mode(int UNUSED(reordering))
	{
	}
	void sTiles_set_correction_mode(int UNUSED(correction_mode))
	{
	}
	void sTiles_set_tile_type_mode(int UNUSED(tile_type))
	{
	}
	void sTiles_set_tile_ordering_mode(int UNUSED(reordering))
	{
	}
	void sTiles_set_tile_ordering_size(int UNUSED(tile_size))
	{
	}
	void sTiles_set_tile_ordering_min_dim(int UNUSED(min_dim))
	{
	}
	void sTiles_turn_off_rescale(int, void *) EMPTY_FUNCTION;
	void sTiles_turn_on_rescale(int, void *) EMPTY_FUNCTION;

#       ifdef __cplusplus
}
#       endif

#       pragma GCC diagnostic pop
#endif
