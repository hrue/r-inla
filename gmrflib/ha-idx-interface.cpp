/* ha-idx-interface.c
 * 
 * Copyright (C) 2022-2022 Havard Rue
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
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

#ifndef GITCOMMIT
#define GITCOMMIT
#endif

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-const-variable"
static const char GitID[] = "file: " __FILE__ "  " GITCOMMIT;
#pragma GCC diagnostic pop

#include "GMRFLib/HArray/stdafx.h"
#include <assert.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include "GMRFLib/HArray/HArrayInt.h"

/// these functions check if 'i' is in the set of 'keys' or not. PS: note that all keys > 0

extern "C" int ha_idx_q(void *ha, int key);
extern "C" void *ha_idx_init(void);
extern "C" void *ha_idx_init_hint(int);
extern "C" void ha_idx_free(void *);
extern "C" void ha_idx_set(void *ha, int key);
extern "C" void ha_idx_sets(void *ha, int n, int *keys);
extern "C" void ha_idx_stats(void *, int print, double *mb);

void *ha_idx_init_hint(int siz)
{
	HArrayInt *ha = NULL;

	if (siz <= 0) {
		siz = 256;
	}
	ha = (HArrayInt *) calloc(1, sizeof(HArrayInt));

       // I am not sure if +1 is needed
	(*ha).init((uint32_t) (log2((double) siz) + 0));

	return (void *) ha;
}

void *ha_idx_init(void)
{
	return ha_idx_init_hint(0);
}

void ha_idx_free(void *ha)
{
	if (ha)
		free(ha);
}

void ha_idx_set(void *ha, int key)
{
	ha_idx_sets(ha, 1, &key);
}

void ha_idx_sets(void *ha, int n, int *keys)
{
	HArrayInt *h = (HArrayInt *) ha;
	for (int i = 0; i < n; i++) {
		assert(keys[i] > 0);
		(*h).insert((uint32_t) keys[i], (uint32_t) 1);
	}
}

int ha_idx_q(void *ha, int key)
{
	HArrayInt *h = (HArrayInt *) ha;
	return ((*h).getValueByKey((uint32_t) key) != 0);
}

void ha_idx_stats(void *ha, int print, double *mb)
{
	HArrayInt *h = (HArrayInt *) ha;
	double mem = (*h).getTotalMemory() / 1024.0 / 1024.0;
	if (print) {
		printf("ha_idx[%p]: Total memory %.2fMb\n", h, mem);
	}
	if (mb) {
		*mb = mem;
	}
}
