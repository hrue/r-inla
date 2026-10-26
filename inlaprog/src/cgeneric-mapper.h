#include <stdio.h>

#ifndef __INLA_CGENERIC_MAPPER_H__
#       define __INLA_CGENERIC_MAPPER_H__
#       undef __BEGIN_DECLS
#       undef __END_DECLS
#       ifdef __cplusplus
#              define __BEGIN_DECLS extern "C" {
#              define __END_DECLS }
#       else
#              define __BEGIN_DECLS			       /* empty */
#              define __END_DECLS			       /* empty */
#       endif

#       include "cgeneric.h"

__BEGIN_DECLS typedef struct {
	const char *name;
	inla_cgeneric_func_tp *func;
} inla_cgeneric_mapper_elm_tp;

typedef struct {
	const char *name;
	inla_cloglike_func_tp *func;
} inla_cloglike_mapper_elm_tp;

void inla_cgeneric_mapper_list(FILE * fp);
inla_cgeneric_func_tp *inla_cgeneric_mapper(char *name);
void inla_cloglike_mapper_list(FILE * fp);
inla_cloglike_func_tp *inla_cloglike_mapper(char *name);

#       if __has_include("cgeneric-defs.h")
#              include "cgeneric-defs.h"
#       endif
#       if __has_include("cloglike-defs.h")
#              include "cloglike-defs.h"
#       endif

__END_DECLS
#endif
