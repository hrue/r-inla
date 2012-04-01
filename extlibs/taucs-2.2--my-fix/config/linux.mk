#########################################################
# Linux                                                 #
#########################################################
OBJEXT=.o
LIBEXT=.a
EXEEXT= 
F2CEXT=.f
PATHSEP=/
DEFFLG=-D

CC = gcc
FC = gfortran

FFLAGS    = -O3 -march=native -mfpmath=sse -msse4.1 -fomit-frame-pointer -funroll-loops -ftracer -floop-block -floop-strip-mine
##FFLAGS    = -m32 -mfpmath=sse -msse2 -O3 -mtune=generic -fomit-frame-pointer -funroll-all-loops -ftracer -floop-block -floop-strip-mine

CFLAGS= -O3 -march=native -mtune=native -mfpmath=sse -msse4 -funroll-loops -ftracer -pipe -floop-block -floop-strip-mine
##CFLAGS= -m32 -mfpmath=sse -msse2 -O3 -mtune=generic -fomit-frame-pointer -funroll-all-loops -ftracer -floop-block -floop-strip-mine


##CFLAGS= -m32 -mfpmath=sse -msse2 -O3 -mtune=generic -ftracer

FOUTFLG   =-o 
COUTFLG   = -o
#CC = icc
#FC = ifort
#FFLAGS = -O3 -no-prec-div -xHost
#CFLAGS= -O3 -no-prec-div -xHost 




LD        = $(CC) 
LDFLAGS   = 
LOUTFLG   = $(COUTFLG)

AR        = ar cr
AOUTFLG   =

RANLIB    = ranlib
RM        = rm -rf

# These are for a Pentium4 version of ATLAS (not bundled with taucs)
#LIBBLAS   = -L /home/stoledo/Public/Linux_P4SSE2/lib -lf77blas -lcblas -latlas \
#            -L /usr/lib/gcc-lib/i386-redhat-linux/2.96 -lg2c
#LIBLAPACK = -L /home/stoledo/Public/Linux_P4SSE2/lib -llapack

LIBBLAS   = -L external/lib/linux -lf77blas -lcblas -latlas
LIBLAPACK = -L external/lib/linux -llapack

LIBMETIS  = -L external/lib/linux -lmetis 

LIBF77 = -lgfortran
LIBC   = -lm 

#########################################################







