#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

#define MSG(fun_) fprintf(stderr, "\n\n*** Warning *** Fake function [%s] called. This may, or may not, go well.\n\n\n", fun_)
#define mkfun(fun_) void fun_(void) { static char first = 1;  if (first) MSG(# fun_); first = 0; return; }

#if defined(WINDOWS) && defined(INLA_WITH_MKL) && (defined(__MINGW32__) || defined(__MINGW64__))

mkfun(__GSHandlerCheck);

//mkfun(__chkstk);

#       if 0
// Force the compiler to recognize this as a naked function (no prologue/epilogue)
__declspec(naked)
void __chkstk()
{
	__asm {
		push rcx				       // Save registers we use
		 push rax cmp rax, 0x1000		       // If size < 4096, only one probe is needed
		 lea rcx,[rsp + 24]			       // Calculate the original RSP
		jb last_probe probing_loop: sub rcx, 0x1000    // Move down one page
		 test[rcx], rcx				       // "Touch" memory to trigger guard page
		 sub rax, 0x1000			       // Decrement remaining size
		 cmp rax, 0x1000 ja probing_loop last_probe: sub rcx, rax	// Subtract final remainder
		 test[rcx], rcx				       // Final touch
		 pop rax				       // Restore original registers
 pop rcx ret}}
#       endif
// does not work void __chkstk(size_t size)
{
	if (size == 0)
		return;
	size_t pages = (size + 4095) / 4096;		       // Round up to page boundary
	volatile char *probe_ptr = (volatile char *) &probe_ptr - size;
	for (size_t i = 0; i < pages; ++i) {
		*probe_ptr = 0;
		probe_ptr += 4096;
	}
}

uint64_t __security_cookie;
void __security_init_cookie()
{
	__security_cookie = 0;
}
void __security_check_cookie(uint64_t retrieved)
{
	if (__security_cookie != retrieved) {
		// abort();
	}
}

#endif

#if defined(WINDOWS) && defined(INLA_WITH_OPENBLAS) && (defined(__MINGW32__) || defined(__MINGW64__))
mkfun(__imp__cprintf);
#endif
