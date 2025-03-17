#ifndef __GMRFLib_SMTP_STILES_H__
#define __GMRFLib_SMTP_STILES_H__

#include <stdlib.h>

#undef __BEGIN_DECLS
#undef __END_DECLS
#ifdef __cplusplus
#define __BEGIN_DECLS extern "C" {
#define __END_DECLS }
#else
#define __BEGIN_DECLS					       /* empty */
#define __END_DECLS					       /* empty */
#endif

__BEGIN_DECLS
//
//

typedef struct 
{
	int verbose;
	int debug;
}
	GMRFLib_stiles_ctl_tp;

int GMRFLib_stiles_set_param(int verbose, int debug);

__END_DECLS
#endif
