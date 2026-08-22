#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <stdio.h>

#include "cgeneric.h"
#if __has_include("cgeneric-defs.h")
#       include "cgeneric-defs.h"
#endif
#if __has_include("cloglike-defs.h")
#       include "cloglike-defs.h"
#endif

typedef struct {
	const char *name;
	inla_cgeneric_func_tp *func;
} inla_cgeneric_mapper_elm_tp;

static inla_cgeneric_mapper_elm_tp table_cgeneric[] = {
#if __has_include("cgeneric-table.h")
#       include "cgeneric-table.h"
#endif
	{ (const char *) NULL, (inla_cgeneric_func_tp *) NULL }
};

void inla_cgeneric_mapper_list(FILE *fp)
{
	fp = (fp ? fp : stdout);
	int i = 0;
	while (table_cgeneric[i].name && table_cgeneric[i].func) {
		fprintf(fp, "\ttable_cgeneric[%1d] = { name = %s, func.ptr = %p }\n", i, table_cgeneric[i].name, (void *) table_cgeneric[i].func);
		i++;
	}
}

inla_cgeneric_func_tp *inla_cgeneric_mapper(char *name)
{
	int i = 0;
	while (name && table_cgeneric[i].name) {
		if (!strcmp(name, table_cgeneric[i].name)) {
			return table_cgeneric[i].func;
		}
		i++;
	}
	return NULL;
}

typedef struct {
	const char *name;
	inla_cloglike_func_tp *func;
} inla_cloglike_mapper_elm_tp;

static inla_cloglike_mapper_elm_tp table_cloglike[] = {
#if __has_include("cloglike-table.h")
#       include "cloglike-table.h"
#endif
	{ (const char *) NULL, (inla_cloglike_func_tp *) NULL }
};

void inla_cloglike_mapper_list(FILE *fp)
{
	fp = (fp ? fp : stdout);
	int i = 0;
	while (table_cloglike[i].name && table_cloglike[i].func) {
		fprintf(fp, "\ttable_cloglike[%1d] = { name = %s, func.ptr = %p }\n", i, table_cloglike[i].name, (void *) table_cloglike[i].func);
		i++;
	}
}

inla_cloglike_func_tp *inla_cloglike_mapper(char *name)
{
	int i = 0;
	while (name && table_cloglike[i].name) {
		if (!strcmp(name, table_cloglike[i].name)) {
			return table_cloglike[i].func;
		}
		i++;
	}
	return NULL;
}
