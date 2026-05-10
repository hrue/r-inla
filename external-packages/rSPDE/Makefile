download :; git clone --branch inla --single-branch git@github.com:davidbolin/rSPDE.git
build :;  make -C rSPDE/src -f ../../../Makefile.generic -B CC="$(CC)" CXX="$(CXX)" FC="$(FC)" FLAGS="$(FLAGS)" INC="$(INC)"
clean :; rm -rf rSPDE
.PHONY: download build clean
