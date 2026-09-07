## TEMPORARY PIN, remove when upstream compiles again.
##
## The tip of the inla branch does not build, and because every lane compiles
## the external packages, it fails ALL of them on every platform:
##
##   src/cgeneric_kronecker.c:201,227,423,426
##       lt_dlclose(d12cache>-handle1);      <- the arrow is transposed, so it
##                                              parses as d12cache > (-handle1)
##                                              and handle1 reads as undeclared
##   src/cgeneric_kronecker.c:182 vs 421
##       static int ltdl_cgkron is declared inside the INITIAL block but read
##       in the QUIT block, so it is out of scope there
##
## Introduced by 563df4368 (2026-09-07). e2369e260 (2026-09-02) is the newest
## commit that compiles, so pin there rather than pinning a whole branch or
## carrying a patch. Reported upstream; delete the checkout line below once it
## is fixed, so we go back to tracking the branch.
INLATOOLS_PIN = e2369e260

download :; git clone --single-branch --branch inla git@github.com:eliaskrainski/INLAtools && git -C INLAtools checkout -q $(INLATOOLS_PIN)
build :;  make -C INLAtools/src -f ../../../Makefile.generic -B CC="$(CC)" CXX="$(CXX)" FC="$(FC)" FLAGS="$(FLAGS)" INC="$(INC)"
clean :; rm -rf INLAtools
.PHONY: download build clean
