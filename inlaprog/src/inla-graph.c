
/* inla-graph.c
 * 
 * Copyright (C) 2007-2024 Havard Rue
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or (at
 * your option) any later version.
 *  * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 *
 * The author's contact information:
 *
 *        Haavard Rue
 *        CEMSE Division
 *        King Abdullah University of Science and Technology
 *        Thuwal 23955-6900, Saudi Arabia
 *        Email: haavard.rue@kaust.edu.sa
 *        Office: +966 (0)12 808 0640
 *
 */

int inla_make_ar1c_graph(GMRFLib_graph_tp **graph_out, inla_ar1c_arg_tp *def)
{
	int i, j;
	GMRFLib_ged_tp *ged = NULL;

	GMRFLib_ged_init(&ged, NULL);
	for (i = 0; i < def->n - 1; i++) {		       /* this is the AR1 part */
		GMRFLib_ged_add(ged, i, i + 1);		       /* diag is added in _ged_add() */
	}
	for (i = 0; i < def->m; i++) {			       /* fill the dense beta-block */
		for (j = i + 1; j < def->m; j++) {
			GMRFLib_ged_add(ged, def->n + i, def->n + j);
		}
	}
	for (i = 0; i < def->n; i++) {			       /* the interaction */
		for (j = 0; j < def->m; j++) {
			GMRFLib_ged_add(ged, i, def->n + j);
		}
	}
	assert(ged->n == def->N);
	GMRFLib_ged_build(graph_out, ged);
	GMRFLib_ged_free(ged);

	return GMRFLib_SUCCESS;
}

int inla_make_besag2_graph(GMRFLib_graph_tp **graph_out, GMRFLib_graph_tp *graph)
{
	int i;
	GMRFLib_ged_tp *ged = NULL;

	GMRFLib_ged_init(&ged, graph);
	for (i = 0; i < graph->n; i++) {
		GMRFLib_ged_add(ged, i, i + graph->n);
	}
	GMRFLib_ged_build(graph_out, ged);
	GMRFLib_ged_free(ged);

	return GMRFLib_SUCCESS;
}

int inla_make_2diid_graph(GMRFLib_graph_tp **graph, inla_2diid_arg_tp *arg)
{
	int i;
	GMRFLib_ged_tp *ged = NULL;

	GMRFLib_ged_init(&ged, NULL);
	for (i = 0; i < 2 * arg->n; i += 2) {
		GMRFLib_ged_add(ged, i, i + 1);
	}
	GMRFLib_ged_build(graph, ged);
	GMRFLib_ged_free(ged);

	return GMRFLib_SUCCESS;
}

int inla_make_2diid_wishart_graph(GMRFLib_graph_tp **graph, inla_2diid_arg_tp *arg)
{
	int i;
	GMRFLib_ged_tp *ged = NULL;

	GMRFLib_ged_init(&ged, NULL);
	for (i = 0; i < arg->n; i++) {
		GMRFLib_ged_add(ged, i, i + arg->n);
	}
	GMRFLib_ged_build(graph, ged);
	GMRFLib_ged_free(ged);

	return GMRFLib_SUCCESS;
}

int inla_make_iid_wishart_graph(GMRFLib_graph_tp **graph, inla_iid_wishart_arg_tp *arg)
{
	int i, j, k, n = arg->n, dim = arg->dim;

	GMRFLib_ged_tp *ged = NULL;

	GMRFLib_ged_init(&ged, NULL);
	for (i = 0; i < n; i++) {
		GMRFLib_ged_add(ged, i, i);
		for (j = 0; j < dim; j++) {
			for (k = j + 1; k < dim; k++) {
				GMRFLib_ged_add(ged, i + j * n, i + k * n);
			}
		}
	}
	GMRFLib_ged_build(graph, ged);
	GMRFLib_ged_free(ged);

	return GMRFLib_SUCCESS;
}

int inla_make_iid_wishartk_graph(GMRFLib_graph_tp **graph, inla_iid_wishartk_arg_tp *arg)
{
	int i, j, k, n = arg->n, dim = arg->dim;
	GMRFLib_ged_tp *ged = NULL;

	GMRFLib_ged_init(&ged, NULL);
	for (i = 0; i < n; i++) {
		GMRFLib_ged_add(ged, i, i);
		for (j = 0; j < dim; j++) {
			for (k = j + 1; k < dim; k++) {
				GMRFLib_ged_add(ged, i + j * n, i + k * n);
			}
		}
	}
	GMRFLib_ged_build(graph, ged);
	GMRFLib_ged_free(ged);

	return GMRFLib_SUCCESS;
}

int inla_make_rw2diid_graph(GMRFLib_graph_tp **graph, GMRFLib_rw2ddef_tp *def)
{
	int i, n;
	GMRFLib_ged_tp *ged = NULL;
	GMRFLib_graph_tp *g;

	n = def->nrow * def->ncol;
	GMRFLib_make_rw2d_graph(&g, def);
	GMRFLib_ged_init2(&ged, n);
	GMRFLib_ged_append_graph(ged, g);
	for (i = 0; i < n; i++) {
		GMRFLib_ged_add(ged, i, i + n);
	}
	GMRFLib_ged_build(graph, ged);
	GMRFLib_ged_free(ged);

	return GMRFLib_SUCCESS;
}

int inla_make_bym_graph(GMRFLib_graph_tp **new_graph, GMRFLib_graph_tp *graph)
{
	/*
	 * for layout: see Qfunc_bym() 
	 */

	int i;
	int n = graph->n;
	GMRFLib_ged_tp *ged = NULL;

	GMRFLib_ged_init2(&ged, n);
	GMRFLib_ged_append_graph(ged, graph);
	for (i = 0; i < n; i++) {
		GMRFLib_ged_add(ged, i, i + n);
	}
	GMRFLib_ged_build(new_graph, ged);
	GMRFLib_ged_free(ged);

	return GMRFLib_SUCCESS;
}

int inla_make_group_graph(GMRFLib_graph_tp **new_graph, GMRFLib_graph_tp *graph, int ngroup, int type, int cyclic, int order,
			  GMRFLib_graph_tp *group_graph)
{
	int i, j, n = graph->n;
	const int debug = 0;
	GMRFLib_ged_tp *ged = NULL;

	GMRFLib_ged_init(&ged, NULL);
	for (i = 0; i < ngroup; i++) {
		GMRFLib_ged_append_graph(ged, graph);
	}

	if (debug) {
		P(n);
		P(ngroup);
		P(graph->n);
		P(group_graph->n);
		P(ged->n);
		P(n * ngroup - 1);
	}

	switch (type) {
	case G_EXCHANGEABLE:
	case G_EXCHANGEABLE_POS:
	{
		assert(cyclic == 0);
		for (i = 0; i < ngroup; i++) {
			for (j = i + 1; j < ngroup; j++) {
				GMRFLib_ged_insert_graph2(ged, graph, i * n, j * n);
			}
		}
	}
		break;

	case G_AR1:
	{
		assert(ngroup >= 2);
		for (i = 0; i < ngroup - 1; i++) {
			GMRFLib_ged_insert_graph2(ged, graph, i * n, (i + 1) * n);
		}
		if (cyclic) {
			GMRFLib_ged_insert_graph2(ged, graph, 0 * n, (ngroup - 1) * n);
		}
	}
		break;

	case G_AR:
	{
		assert(ngroup >= 2);
		for (i = 0; i < ngroup - 1; i++) {
			for (j = 1; j <= order; j++) {
				if (i + j < ngroup) {
					GMRFLib_ged_insert_graph2(ged, graph, i * n, (i + j) * n);
				}
			}
		}
	}
		break;

	case G_RW1:
	{
		assert(ngroup >= 2);
		for (i = 0; i < ngroup - 1; i++) {
			GMRFLib_ged_insert_graph2(ged, graph, i * n, (i + 1) * n);
		}
		if (cyclic) {
			GMRFLib_ged_insert_graph2(ged, graph, 0 * n, (ngroup - 1) * n);
		}
	}
		break;

	case G_RW2:
	{
		assert(ngroup >= 3);
		for (i = 0; i < ngroup - 2; i++) {
			GMRFLib_ged_insert_graph2(ged, graph, i * n, (i + 1) * n);
			GMRFLib_ged_insert_graph2(ged, graph, i * n, (i + 2) * n);
		}
		GMRFLib_ged_insert_graph2(ged, graph, (ngroup - 2) * n, (ngroup - 1) * n);
		if (cyclic) {
			GMRFLib_ged_insert_graph2(ged, graph, 0 * n, (ngroup - 1) * n);
			GMRFLib_ged_insert_graph2(ged, graph, 0 * n, (ngroup - 2) * n);
			GMRFLib_ged_insert_graph2(ged, graph, 1 * n, (ngroup - 1) * n);
		}
	}
		break;

	case G_BESAG:
	{
		assert(group_graph);
		for (i = 0; i < group_graph->n; i++) {
			int jj;
			for (jj = 0; jj < group_graph->lnnbs[i]; jj++) {
				j = group_graph->lnbs[i][jj];
				GMRFLib_ged_insert_graph2(ged, graph, i * n, j * n);
			}
		}
	}
		break;

	case G_IID:
	{
		assert(ngroup >= 1);
		for (i = 0; i < ngroup; i++) {
			GMRFLib_ged_insert_graph2(ged, graph, i * n, i * n);
		}
	}
		break;

	default:
		inla_error_general("This should not happen");
		abort();
	}

	if (0) {
		FILE *fp = fopen("g.dat", "w");
		GMRFLib_printf_graph(fp, new_graph[0]);
		fclose(fp);
	}

	assert(ged->n == n * ngroup);
	GMRFLib_ged_build(new_graph, ged);
	GMRFLib_ged_free(ged);

	return GMRFLib_SUCCESS;
}

int inla_replicate_graph(GMRFLib_graph_tp **g, int replicate)
{
	/*
	 * replace the graph G, with on that is replicated REPLICATE times 
	 */
	int i;
	GMRFLib_ged_tp *ged;

	if (!g || !*g || replicate <= 1) {
		return GMRFLib_SUCCESS;
	}
	GMRFLib_ged_init(&ged, NULL);
	for (i = 0; i < replicate; i++) {
		GMRFLib_ged_append_graph(ged, *g);
	}
	GMRFLib_graph_free(*g);
	GMRFLib_ged_build(g, ged);
	GMRFLib_ged_free(ged);

	return GMRFLib_SUCCESS;
}

int inla_make_ar1_graph(GMRFLib_graph_tp **graph, inla_ar1_arg_tp *arg)
{
	return GMRFLib_graph_mk_linear(graph, arg->n, 1, arg->cyclic);
}

int inla_make_ou_graph(GMRFLib_graph_tp **graph, inla_ou_arg_tp *arg)
{
	return GMRFLib_graph_mk_linear(graph, arg->n, 1, 0);
}

int inla_make_intslope_graph(GMRFLib_graph_tp **graph, inla_intslope_arg_tp *arg)
{
	int idx, subject;
	GMRFLib_ged_tp *ged = NULL;
	GMRFLib_graph_tp *g = NULL;

	// the iid2d part of the model is stored at the end
	GMRFLib_ged_init(&ged, NULL);
	inla_make_iid_wishart_graph(&g, arg->warg);
	GMRFLib_ged_insert_graph(ged, g, arg->n);

	for (idx = 0; idx < arg->n; idx++) {
		subject = (int) GMRFLib_matrix_get(idx, INTSLOPE_SUBJECT, arg->def);
		GMRFLib_ged_add(ged, idx, arg->n + subject);
		GMRFLib_ged_add(ged, idx, arg->n + arg->nsubject + subject);
	}
	GMRFLib_ged_build(graph, ged);
	GMRFLib_ged_free(ged);
	GMRFLib_graph_free(g);

	if (0) {
		GMRFLib_printf_graph(stdout, *graph);
		GMRFLib_graph_write("GRAPH", *graph);
	}

	return GMRFLib_SUCCESS;
}
