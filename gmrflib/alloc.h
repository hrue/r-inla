#ifndef __GMRFLib_ALLOC_H__
#define __GMRFLib_ALLOC_H__

#undef __BEGIN_DECLS
#undef __END_DECLS

#ifdef __cplusplus
#define __BEGIN_DECLS extern "C" {
#define __END_DECLS }
#else
#define __BEGIN_DECLS                                          /* empty */
#define __END_DECLS                                            /* empty */
#endif

__BEGIN_DECLS

// disable 32u for the moment
#if 0 && defined(__linux) && defined(__AVX2__)
#define GMRFLib_MEM_ALIGN 32u
#else
#define GMRFLib_MEM_ALIGN 16u
#endif

size_t GMRFLib_align_len(size_t n, size_t size);
void *calloc_intern(size_t nmemb, size_t size);
void *malloc_intern(size_t size);

#define GMRFLib_MEM_ALIGN_TEST_(p_) ((uintptr_t) (p_) % GMRFLib_MEM_ALIGN == 0)
#define GMRFLib_is_aligned(p_) GMRFLib_MEM_ALIGN_TEST_(p_)
#define GMRFLib_is_aligned2(p_, pp_) (GMRFLib_MEM_ALIGN_TEST_(p_) && GMRFLib_MEM_ALIGN_TEST_(pp_))
#define GMRFLib_is_aligned3(p_, pp_, ppp_) (GMRFLib_MEM_ALIGN_TEST_(p_) && GMRFLib_MEM_ALIGN_TEST_(pp_) && GMRFLib_MEM_ALIGN_TEST_(ppp_))

__END_DECLS
#endif
