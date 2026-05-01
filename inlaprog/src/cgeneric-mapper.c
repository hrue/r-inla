#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <stdio.h>

#include "cgeneric.h"
#include "../../external-packages/cgeneric-defs.h"

typedef struct {
	char *name;
	inla_cgeneric_func_tp *func;
} inla_func_mapper_elm_tp;

void inla_cgeneric_mapper_list(FILE *fp)
{
	inla_func_mapper_elm_tp table[] = {
#if __has_include("../../external-packages/cgeneric-table.h")
#       include "../../external-packages/cgeneric-table.h"
#elif __has_include("../../../external-packages/cgeneric-table.h")
#       include "../../../external-packages/cgeneric-table.h"
#else
#       error "Cannot find file 'cgeneric-table.h'"
#endif
		{ (char *) NULL, (inla_cgeneric_func_tp *) NULL }
	};

	if (!fp)
		fp = stdout;
	int i = 0;
	while (table[i].name && table[i].func) {
		fprintf(fp, "\ttable[%1d] = { name = %s, cgeneric_func.ptr = %p }\n", i, table[i].name, (void *) table[i].func);
		i++;
	}
}

inla_cgeneric_func_tp *inla_cgeneric_mapper(char *name)
{
	inla_func_mapper_elm_tp table[] = {
#if __has_include("../../external-packages/cgeneric-table.h")
#       include "../../external-packages/cgeneric-table.h"
#elif __has_include("../../../external-packages/cgeneric-table.h")
#       include "../../../external-packages/cgeneric-table.h"
#else
#       error "Cannot find file 'cgeneric-table.h'"
#endif
		{ (char *) NULL, (inla_cgeneric_func_tp *) NULL }
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
