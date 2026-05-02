download :; git clone --single-branch --branch inla git@github.com:eliaskrainski/graphpcor
build :;  make -C graphpcor/src -f ../../../Makefile.generic CC="$(CC)" CXX="$(CXX)" FC="$(FC)" FLAGS="$(FLAGS)" INC="$(INC)"
clean :; rm -rf graphpcor
.PHONY: download build clean
