CC = gcc
FC = gfortran
FFLAGS = -O3 
CFLAGS = -O3 -fopenmp
AR = ar
RANLIB = ranlib
LIBNAME = libtaucs.a

LIBOBJ = obj/amdatr.o obj/amdbarNEW.o obj/amdbar.o obj/amdexa.o \
	obj/amdhaf.o obj/amdhat.o obj/amdNEW.o obj/amdpre.o obj/amdtru.o \
	obj/colamd.o obj/genmmd.o obj/readhb.o obj/taucs_ccs_base_D.o \
	obj/taucs_ccs_base.o obj/taucs_ccs_factor_llt_D.o obj/taucs_ccs_factor_llt.o \
	obj/taucs_ccs_ops_D.o obj/taucs_ccs_ops.o obj/taucs_ccs_order.o \
	obj/taucs_ccs_solve_llt_D.o obj/taucs_ccs_solve_llt.o obj/taucs_malloc.o \
	obj/taucs_sn_llt_D.o obj/taucs_sn_llt.o obj/taucs_timer.o \
	amd/src/amd_1.o	amd/src/amd_2.o amd/src/amd_aat.o amd/src/amdbar.o \
	amd/src/amd_control.o amd/src/amd_defaults.o amd/src/amd_dump.o \
	amd/src/amd_global.o amd/src/amd_info.o amd/src/amd.o amd/src/amd_order.o \
	amd/src/amd_postorder.o amd/src/amd_post_tree.o amd/src/amd_preprocess.o \
	amd/src/amd_valid.o

$(LIBNAME) :
	$(MAKE) -f Makefile.taucs CC="$(CC)" FC="$(FC)" FFLAGS="$(FFLAGS)" CFLAGS="$(CFLAGS)"
	$(MAKE) -f Makefile.amd   CC="$(CC)" FC="$(FC)" FFLAGS="$(FFLAGS)" CFLAGS="$(CFLAGS)"
	$(RM) -f $(LIBNAME)
	$(AR) crv $(LIBNAME) $(LIBOBJ)
	$(RANLIB) $(LIBNAME)

clean :; $(RM) -f $(LIBOBJ) $(LIBNAME)
.PHONY: clean all

