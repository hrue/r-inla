#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <stdio.h>

#include "cgeneric.h"
#include "cgeneric-defs.h"

typedef struct {
	const char *name;
	inla_cgeneric_func_tp *func;
} inla_func_mapper_elm_tp;

void inla_cgeneric_mapper_list(FILE *fp)
{
	inla_func_mapper_elm_tp table[] = {
#include "cgeneric-table.h"
		{ (const char *) NULL, (inla_cgeneric_func_tp *) NULL }
	};

	if (!fp)
		fp = stdout;
	int i = 0;
	while (table[i].name && table[i].func) {
		fprintf(fp, "\ttable[%1d] = { name = %s, func.ptr = %p }\n", i, table[i].name, (void *) table[i].func);
		i++;
	}
}

inla_cgeneric_func_tp *inla_cgeneric_mapper(char *name)
{
	inla_func_mapper_elm_tp table[] = {
#include "cgeneric-table.h"
		{ (const char *) NULL, (inla_cgeneric_func_tp *) NULL }
	};

	int i = 0;
	while (name && table[i].name) {
		if (!strcmp(name, table[i].name)) {
			return table[i].func;
		}
		i++;
	}
	return NULL;
}
