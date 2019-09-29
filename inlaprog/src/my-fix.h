#ifndef _MY_FIX
#define _MY_FIX

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
// ...

#if defined(INLA_WINDOWS32_FIX)
void _mm_pause(void);
#endif

int my_is_double(char *str);
int my_is_int(char *str);
char *my_strlwc(const char *str);

__END_DECLS
#endif
