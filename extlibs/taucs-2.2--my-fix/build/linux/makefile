OSTYPE = linux
CONFIGURATION = anonymous
INCS = \
  build/linux/taucs_config_build.h\
  build/linux/taucs_config_tests.h\
  $(DIRSRC)taucs.h\
  $(DIRSRC)taucs_private.h\
  $(DIRBLD)taucs_config_tests.h
CILKC=$(CC)
CILKFLAGS=$(CFLAGS)
CILKOUTFLG=$(COUTFLG)


include config/$(OSTYPE).mk

LIBS = $(LIBLAPACK) $(LIBBLAS) $(LIBMETIS) $(LIBF77) $(LIBC)
default: all

include config/std.mk

build/linux/taucs_config_tests.h: $(DIROBJ)exists.log 
#	- $(CC) -c $(CFLAGS) $(STDDEFS) $(STDINCS) \
#	  $(DIRPROGS)taucs_blas_underscore_test.c \
#	  $(COUTFLG)$(DIROBJ)taucs_blas_underscore_test$(OBJEXT)
#	- $(LD) $(LDFLAGS) \
#	  $(LOUTFLG)$(DIROBJ)taucs_blas_underscore_test$(EXEEXT) \
#	  $(DIROBJ)taucs_blas_underscore_test$(OBJEXT) $(LIBS)
#	- $(DIROBJ)taucs_blas_underscore_test$(EXEEXT) build/linux/taucs_config_tests.h
#	- $(CC) -c $(CFLAGS) $(STDDEFS) $(STDINCS) \
#	  $(DIRPROGS)taucs_blas_nounderscore_test.c \
#	  $(COUTFLG)$(DIROBJ)taucs_blas_nounderscore_test$(OBJEXT)
#	- $(LD) $(LDFLAGS) \
#	  $(LOUTFLG)$(DIROBJ)taucs_blas_nounderscore_test$(EXEEXT) \
#	  $(DIROBJ)taucs_blas_nounderscore_test$(OBJEXT) $(LIBS)
#	- $(DIROBJ)taucs_blas_nounderscore_test$(EXEEXT) build/linux/taucs_config_tests.h
#	- $(CC) -c $(CFLAGS) $(STDDEFS) $(STDINCS) \
#	  $(DIRPROGS)taucs_c99_complex_test.c \
#	  $(COUTFLG)$(DIROBJ)taucs_c99_complex_test$(OBJEXT)
#	- $(LD) $(LDFLAGS) \
#	  $(LOUTFLG)$(DIROBJ)taucs_c99_complex_test$(EXEEXT) \
#	  $(DIROBJ)taucs_c99_complex_test$(OBJEXT) $(LIBS)
#	- $(DIROBJ)taucs_c99_complex_test$(EXEEXT) build/linux/taucs_config_tests.h
#	- $(CILKC) -c $(CILKFLAGS) $(STDDEFS) $(STDINCS) \
#	  $(DIRPROGS)taucs_cilk_test.c \
#	  $(CILKOUTFLG)$(DIROBJ)taucs_cilk_test$(OBJEXT)
#	- $(LD) $(LDFLAGS) \
#	  $(LOUTFLG)$(DIROBJ)taucs_cilk_test$(EXEEXT) \
#	  $(DIROBJ)taucs_cilk_test$(OBJEXT) $(LIBS)
#	- $(DIROBJ)taucs_cilk_test$(EXEEXT) build/linux/taucs_config_tests.h
$(DIROBJ)taucs_sn_llt$(OBJEXT): $(DIRSRC)taucs_sn_llt.c $(INCS) $(STDDEPS)
	$(CILKC) -c $(CILKFLAGS) $(STDDEFS) $(STDINCS) \
	-DTAUCS_CORE_GENERAL \
	$(DIRSRC)taucs_sn_llt.c \
	$(CILKOUTFLG)$(DIROBJ)taucs_sn_llt$(OBJEXT)
$(DIROBJ)taucs_sn_llt_D$(OBJEXT): $(DIRSRC)taucs_sn_llt.c $(INCS) $(STDDEPS)
	$(CILKC) -c $(CILKFLAGS) $(STDDEFS) $(STDINCS) \
	-DTAUCS_CORE_DOUBLE \
	$(DIRSRC)taucs_sn_llt.c \
	$(CILKOUTFLG)$(DIROBJ)taucs_sn_llt_D$(OBJEXT)
$(DIROBJ)taucs_sn_llt_S$(OBJEXT): $(DIRSRC)taucs_sn_llt.c $(INCS) $(STDDEPS)
	$(CILKC) -c $(CILKFLAGS) $(STDDEFS) $(STDINCS) \
	-DTAUCS_CORE_SINGLE \
	$(DIRSRC)taucs_sn_llt.c \
	$(CILKOUTFLG)$(DIROBJ)taucs_sn_llt_S$(OBJEXT)
$(DIROBJ)taucs_sn_llt_Z$(OBJEXT): $(DIRSRC)taucs_sn_llt.c $(INCS) $(STDDEPS)
	$(CILKC) -c $(CILKFLAGS) $(STDDEFS) $(STDINCS) \
	-DTAUCS_CORE_DCOMPLEX \
	$(DIRSRC)taucs_sn_llt.c \
	$(CILKOUTFLG)$(DIROBJ)taucs_sn_llt_Z$(OBJEXT)
$(DIROBJ)taucs_sn_llt_C$(OBJEXT): $(DIRSRC)taucs_sn_llt.c $(INCS) $(STDDEPS)
	$(CILKC) -c $(CILKFLAGS) $(STDDEFS) $(STDINCS) \
	-DTAUCS_CORE_SCOMPLEX \
	$(DIRSRC)taucs_sn_llt.c \
	$(CILKOUTFLG)$(DIROBJ)taucs_sn_llt_C$(OBJEXT)
$(DIROBJ)taucs_linsolve$(OBJEXT): $(DIRSRC)taucs_linsolve.c $(INCS) $(STDDEPS)
	$(CILKC) -c $(CILKFLAGS) $(STDDEFS) $(STDINCS) \
	-DTAUCS_CORE_GENERAL \
	$(DIRSRC)taucs_linsolve.c \
	$(CILKOUTFLG)$(DIROBJ)taucs_linsolve$(OBJEXT)
$(DIROBJ)taucs_logging$(OBJEXT): $(DIRSRC)taucs_logging.c $(INCS) $(STDDEPS)
	$(CC) -c $(CFLAGS) $(STDDEFS) $(STDINCS) \
	-DTAUCS_CORE_GENERAL \
	$(DIRSRC)taucs_logging.c \
	$(COUTFLG)$(DIROBJ)taucs_logging$(OBJEXT)
$(DIROBJ)taucs_memory$(OBJEXT): $(DIRSRC)taucs_memory.c $(INCS) $(STDDEPS)
	$(CC) -c $(CFLAGS) $(STDDEFS) $(STDINCS) \
	-DTAUCS_CORE_GENERAL \
	$(DIRSRC)taucs_memory.c \
	$(COUTFLG)$(DIROBJ)taucs_memory$(OBJEXT)
$(DIROBJ)taucs_timer$(OBJEXT): $(DIRSRC)taucs_timer.c $(INCS) $(STDDEPS)
	$(CC) -c $(CFLAGS) $(STDDEFS) $(STDINCS) \
	-DTAUCS_CORE_GENERAL \
	$(DIRSRC)taucs_timer.c \
	$(COUTFLG)$(DIROBJ)taucs_timer$(OBJEXT)
$(DIROBJ)taucs_ccs_base$(OBJEXT): $(DIRSRC)taucs_ccs_base.c $(INCS) $(STDDEPS)
	$(CC) -c $(CFLAGS) $(STDDEFS) $(STDINCS) \
	-DTAUCS_CORE_GENERAL \
	$(DIRSRC)taucs_ccs_base.c \
	$(COUTFLG)$(DIROBJ)taucs_ccs_base$(OBJEXT)
$(DIROBJ)taucs_ccs_base_D$(OBJEXT): $(DIRSRC)taucs_ccs_base.c $(INCS) $(STDDEPS)
	$(CC) -c $(CFLAGS) $(STDDEFS) $(STDINCS) \
	-DTAUCS_CORE_DOUBLE \
	$(DIRSRC)taucs_ccs_base.c \
	$(COUTFLG)$(DIROBJ)taucs_ccs_base_D$(OBJEXT)
$(DIROBJ)taucs_ccs_base_S$(OBJEXT): $(DIRSRC)taucs_ccs_base.c $(INCS) $(STDDEPS)
	$(CC) -c $(CFLAGS) $(STDDEFS) $(STDINCS) \
	-DTAUCS_CORE_SINGLE \
	$(DIRSRC)taucs_ccs_base.c \
	$(COUTFLG)$(DIROBJ)taucs_ccs_base_S$(OBJEXT)
$(DIROBJ)taucs_ccs_base_Z$(OBJEXT): $(DIRSRC)taucs_ccs_base.c $(INCS) $(STDDEPS)
	$(CC) -c $(CFLAGS) $(STDDEFS) $(STDINCS) \
	-DTAUCS_CORE_DCOMPLEX \
	$(DIRSRC)taucs_ccs_base.c \
	$(COUTFLG)$(DIROBJ)taucs_ccs_base_Z$(OBJEXT)
$(DIROBJ)taucs_ccs_base_C$(OBJEXT): $(DIRSRC)taucs_ccs_base.c $(INCS) $(STDDEPS)
	$(CC) -c $(CFLAGS) $(STDDEFS) $(STDINCS) \
	-DTAUCS_CORE_SCOMPLEX \
	$(DIRSRC)taucs_ccs_base.c \
	$(COUTFLG)$(DIROBJ)taucs_ccs_base_C$(OBJEXT)
$(DIROBJ)taucs_vec_base$(OBJEXT): $(DIRSRC)taucs_vec_base.c $(INCS) $(STDDEPS)
	$(CC) -c $(CFLAGS) $(STDDEFS) $(STDINCS) \
	-DTAUCS_CORE_GENERAL \
	$(DIRSRC)taucs_vec_base.c \
	$(COUTFLG)$(DIROBJ)taucs_vec_base$(OBJEXT)
$(DIROBJ)taucs_vec_base_D$(OBJEXT): $(DIRSRC)taucs_vec_base.c $(INCS) $(STDDEPS)
	$(CC) -c $(CFLAGS) $(STDDEFS) $(STDINCS) \
	-DTAUCS_CORE_DOUBLE \
	$(DIRSRC)taucs_vec_base.c \
	$(COUTFLG)$(DIROBJ)taucs_vec_base_D$(OBJEXT)
$(DIROBJ)taucs_vec_base_S$(OBJEXT): $(DIRSRC)taucs_vec_base.c $(INCS) $(STDDEPS)
	$(CC) -c $(CFLAGS) $(STDDEFS) $(STDINCS) \
	-DTAUCS_CORE_SINGLE \
	$(DIRSRC)taucs_vec_base.c \
	$(COUTFLG)$(DIROBJ)taucs_vec_base_S$(OBJEXT)
$(DIROBJ)taucs_vec_base_Z$(OBJEXT): $(DIRSRC)taucs_vec_base.c $(INCS) $(STDDEPS)
	$(CC) -c $(CFLAGS) $(STDDEFS) $(STDINCS) \
	-DTAUCS_CORE_DCOMPLEX \
	$(DIRSRC)taucs_vec_base.c \
	$(COUTFLG)$(DIROBJ)taucs_vec_base_Z$(OBJEXT)
$(DIROBJ)taucs_vec_base_C$(OBJEXT): $(DIRSRC)taucs_vec_base.c $(INCS) $(STDDEPS)
	$(CC) -c $(CFLAGS) $(STDDEFS) $(STDINCS) \
	-DTAUCS_CORE_SCOMPLEX \
	$(DIRSRC)taucs_vec_base.c \
	$(COUTFLG)$(DIROBJ)taucs_vec_base_C$(OBJEXT)
$(DIROBJ)taucs_ccs_ops$(OBJEXT): $(DIRSRC)taucs_ccs_ops.c $(INCS) $(STDDEPS)
	$(CC) -c $(CFLAGS) $(STDDEFS) $(STDINCS) \
	-DTAUCS_CORE_GENERAL \
	$(DIRSRC)taucs_ccs_ops.c \
	$(COUTFLG)$(DIROBJ)taucs_ccs_ops$(OBJEXT)
$(DIROBJ)taucs_ccs_ops_D$(OBJEXT): $(DIRSRC)taucs_ccs_ops.c $(INCS) $(STDDEPS)
	$(CC) -c $(CFLAGS) $(STDDEFS) $(STDINCS) \
	-DTAUCS_CORE_DOUBLE \
	$(DIRSRC)taucs_ccs_ops.c \
	$(COUTFLG)$(DIROBJ)taucs_ccs_ops_D$(OBJEXT)
$(DIROBJ)taucs_ccs_ops_S$(OBJEXT): $(DIRSRC)taucs_ccs_ops.c $(INCS) $(STDDEPS)
	$(CC) -c $(CFLAGS) $(STDDEFS) $(STDINCS) \
	-DTAUCS_CORE_SINGLE \
	$(DIRSRC)taucs_ccs_ops.c \
	$(COUTFLG)$(DIROBJ)taucs_ccs_ops_S$(OBJEXT)
$(DIROBJ)taucs_ccs_ops_Z$(OBJEXT): $(DIRSRC)taucs_ccs_ops.c $(INCS) $(STDDEPS)
	$(CC) -c $(CFLAGS) $(STDDEFS) $(STDINCS) \
	-DTAUCS_CORE_DCOMPLEX \
	$(DIRSRC)taucs_ccs_ops.c \
	$(COUTFLG)$(DIROBJ)taucs_ccs_ops_Z$(OBJEXT)
$(DIROBJ)taucs_ccs_ops_C$(OBJEXT): $(DIRSRC)taucs_ccs_ops.c $(INCS) $(STDDEPS)
	$(CC) -c $(CFLAGS) $(STDDEFS) $(STDINCS) \
	-DTAUCS_CORE_SCOMPLEX \
	$(DIRSRC)taucs_ccs_ops.c \
	$(COUTFLG)$(DIROBJ)taucs_ccs_ops_C$(OBJEXT)

$(DIROBJ)taucs_ccs_io$(OBJEXT): $(DIRSRC)taucs_ccs_io.c $(INCS) $(STDDEPS)
	$(CC) -c $(CFLAGS) $(STDDEFS) $(STDINCS) \
	-DTAUCS_CORE_GENERAL \
	$(DIRSRC)taucs_ccs_io.c \
	$(COUTFLG)$(DIROBJ)taucs_ccs_io$(OBJEXT)
$(DIROBJ)taucs_ccs_io_D$(OBJEXT): $(DIRSRC)taucs_ccs_io.c $(INCS) $(STDDEPS)
	$(CC) -c $(CFLAGS) $(STDDEFS) $(STDINCS) \
	-DTAUCS_CORE_DOUBLE \
	$(DIRSRC)taucs_ccs_io.c \
	$(COUTFLG)$(DIROBJ)taucs_ccs_io_D$(OBJEXT)
$(DIROBJ)taucs_ccs_io_S$(OBJEXT): $(DIRSRC)taucs_ccs_io.c $(INCS) $(STDDEPS)
	$(CC) -c $(CFLAGS) $(STDDEFS) $(STDINCS) \
	-DTAUCS_CORE_SINGLE \
	$(DIRSRC)taucs_ccs_io.c \
	$(COUTFLG)$(DIROBJ)taucs_ccs_io_S$(OBJEXT)
$(DIROBJ)taucs_ccs_io_Z$(OBJEXT): $(DIRSRC)taucs_ccs_io.c $(INCS) $(STDDEPS)
	$(CC) -c $(CFLAGS) $(STDDEFS) $(STDINCS) \
	-DTAUCS_CORE_DCOMPLEX \
	$(DIRSRC)taucs_ccs_io.c \
	$(COUTFLG)$(DIROBJ)taucs_ccs_io_Z$(OBJEXT)
$(DIROBJ)taucs_ccs_io_C$(OBJEXT): $(DIRSRC)taucs_ccs_io.c $(INCS) $(STDDEPS)
	$(CC) -c $(CFLAGS) $(STDDEFS) $(STDINCS) \
	-DTAUCS_CORE_SCOMPLEX \
	$(DIRSRC)taucs_ccs_io.c \
	$(COUTFLG)$(DIROBJ)taucs_ccs_io_C$(OBJEXT)

$(DIROBJ)taucs_ccs_order$(OBJEXT): $(DIRSRC)taucs_ccs_order.c $(INCS) $(STDDEPS)
	$(CC) -c $(CFLAGS) $(STDDEFS) $(STDINCS) \
	-DTAUCS_CORE_GENERAL \
	$(DIRSRC)taucs_ccs_order.c \
	$(COUTFLG)$(DIROBJ)taucs_ccs_order$(OBJEXT)
$(DIROBJ)taucs_ccs_factor_llt$(OBJEXT): $(DIRSRC)taucs_ccs_factor_llt.c $(INCS) $(STDDEPS)
	$(CC) -c $(CFLAGS) $(STDDEFS) $(STDINCS) \
	-DTAUCS_CORE_GENERAL \
	$(DIRSRC)taucs_ccs_factor_llt.c \
	$(COUTFLG)$(DIROBJ)taucs_ccs_factor_llt$(OBJEXT)
$(DIROBJ)taucs_ccs_factor_llt_D$(OBJEXT): $(DIRSRC)taucs_ccs_factor_llt.c $(INCS) $(STDDEPS)
	$(CC) -c $(CFLAGS) $(STDDEFS) $(STDINCS) \
	-DTAUCS_CORE_DOUBLE \
	$(DIRSRC)taucs_ccs_factor_llt.c \
	$(COUTFLG)$(DIROBJ)taucs_ccs_factor_llt_D$(OBJEXT)
$(DIROBJ)taucs_ccs_factor_llt_S$(OBJEXT): $(DIRSRC)taucs_ccs_factor_llt.c $(INCS) $(STDDEPS)
	$(CC) -c $(CFLAGS) $(STDDEFS) $(STDINCS) \
	-DTAUCS_CORE_SINGLE \
	$(DIRSRC)taucs_ccs_factor_llt.c \
	$(COUTFLG)$(DIROBJ)taucs_ccs_factor_llt_S$(OBJEXT)
$(DIROBJ)taucs_ccs_factor_llt_Z$(OBJEXT): $(DIRSRC)taucs_ccs_factor_llt.c $(INCS) $(STDDEPS)
	$(CC) -c $(CFLAGS) $(STDDEFS) $(STDINCS) \
	-DTAUCS_CORE_DCOMPLEX \
	$(DIRSRC)taucs_ccs_factor_llt.c \
	$(COUTFLG)$(DIROBJ)taucs_ccs_factor_llt_Z$(OBJEXT)
$(DIROBJ)taucs_ccs_factor_llt_C$(OBJEXT): $(DIRSRC)taucs_ccs_factor_llt.c $(INCS) $(STDDEPS)
	$(CC) -c $(CFLAGS) $(STDDEFS) $(STDINCS) \
	-DTAUCS_CORE_SCOMPLEX \
	$(DIRSRC)taucs_ccs_factor_llt.c \
	$(COUTFLG)$(DIROBJ)taucs_ccs_factor_llt_C$(OBJEXT)
$(DIROBJ)taucs_ccs_solve_llt$(OBJEXT): $(DIRSRC)taucs_ccs_solve_llt.c $(INCS) $(STDDEPS)
	$(CC) -c $(CFLAGS) $(STDDEFS) $(STDINCS) \
	-DTAUCS_CORE_GENERAL \
	$(DIRSRC)taucs_ccs_solve_llt.c \
	$(COUTFLG)$(DIROBJ)taucs_ccs_solve_llt$(OBJEXT)
$(DIROBJ)taucs_ccs_solve_llt_D$(OBJEXT): $(DIRSRC)taucs_ccs_solve_llt.c $(INCS) $(STDDEPS)
	$(CC) -c $(CFLAGS) $(STDDEFS) $(STDINCS) \
	-DTAUCS_CORE_DOUBLE \
	$(DIRSRC)taucs_ccs_solve_llt.c \
	$(COUTFLG)$(DIROBJ)taucs_ccs_solve_llt_D$(OBJEXT)
$(DIROBJ)taucs_ccs_solve_llt_S$(OBJEXT): $(DIRSRC)taucs_ccs_solve_llt.c $(INCS) $(STDDEPS)
	$(CC) -c $(CFLAGS) $(STDDEFS) $(STDINCS) \
	-DTAUCS_CORE_SINGLE \
	$(DIRSRC)taucs_ccs_solve_llt.c \
	$(COUTFLG)$(DIROBJ)taucs_ccs_solve_llt_S$(OBJEXT)
$(DIROBJ)taucs_ccs_solve_llt_Z$(OBJEXT): $(DIRSRC)taucs_ccs_solve_llt.c $(INCS) $(STDDEPS)
	$(CC) -c $(CFLAGS) $(STDDEFS) $(STDINCS) \
	-DTAUCS_CORE_DCOMPLEX \
	$(DIRSRC)taucs_ccs_solve_llt.c \
	$(COUTFLG)$(DIROBJ)taucs_ccs_solve_llt_Z$(OBJEXT)
$(DIROBJ)taucs_ccs_solve_llt_C$(OBJEXT): $(DIRSRC)taucs_ccs_solve_llt.c $(INCS) $(STDDEPS)
	$(CC) -c $(CFLAGS) $(STDDEFS) $(STDINCS) \
	-DTAUCS_CORE_SCOMPLEX \
	$(DIRSRC)taucs_ccs_solve_llt.c \
	$(COUTFLG)$(DIROBJ)taucs_ccs_solve_llt_C$(OBJEXT)
$(DIROBJ)taucs_complex$(OBJEXT): $(DIRSRC)taucs_complex.c $(INCS) $(STDDEPS)
	$(CC) -c $(CFLAGS) $(STDDEFS) $(STDINCS) \
	-DTAUCS_CORE_GENERAL \
	$(DIRSRC)taucs_complex.c \
	$(COUTFLG)$(DIROBJ)taucs_complex$(OBJEXT)
$(DIROBJ)taucs_complex_D$(OBJEXT): $(DIRSRC)taucs_complex.c $(INCS) $(STDDEPS)
	$(CC) -c $(CFLAGS) $(STDDEFS) $(STDINCS) \
	-DTAUCS_CORE_DOUBLE \
	$(DIRSRC)taucs_complex.c \
	$(COUTFLG)$(DIROBJ)taucs_complex_D$(OBJEXT)
$(DIROBJ)taucs_complex_S$(OBJEXT): $(DIRSRC)taucs_complex.c $(INCS) $(STDDEPS)
	$(CC) -c $(CFLAGS) $(STDDEFS) $(STDINCS) \
	-DTAUCS_CORE_SINGLE \
	$(DIRSRC)taucs_complex.c \
	$(COUTFLG)$(DIROBJ)taucs_complex_S$(OBJEXT)
$(DIROBJ)taucs_complex_Z$(OBJEXT): $(DIRSRC)taucs_complex.c $(INCS) $(STDDEPS)
	$(CC) -c $(CFLAGS) $(STDDEFS) $(STDINCS) \
	-DTAUCS_CORE_DCOMPLEX \
	$(DIRSRC)taucs_complex.c \
	$(COUTFLG)$(DIROBJ)taucs_complex_Z$(OBJEXT)
$(DIROBJ)taucs_complex_C$(OBJEXT): $(DIRSRC)taucs_complex.c $(INCS) $(STDDEPS)
	$(CC) -c $(CFLAGS) $(STDDEFS) $(STDINCS) \
	-DTAUCS_CORE_SCOMPLEX \
	$(DIRSRC)taucs_complex.c \
	$(COUTFLG)$(DIROBJ)taucs_complex_C$(OBJEXT)

$(DIROBJ)taucs_ccs_ooc_llt$(OBJEXT): $(DIRSRC)taucs_ccs_ooc_llt.c $(INCS) $(STDDEPS)
	$(CC) -c $(CFLAGS) $(STDDEFS) $(STDINCS) \
	-DTAUCS_CORE_GENERAL \
	$(DIRSRC)taucs_ccs_ooc_llt.c \
	$(COUTFLG)$(DIROBJ)taucs_ccs_ooc_llt$(OBJEXT)
$(DIROBJ)taucs_ccs_ooc_llt_D$(OBJEXT): $(DIRSRC)taucs_ccs_ooc_llt.c $(INCS) $(STDDEPS)
	$(CC) -c $(CFLAGS) $(STDDEFS) $(STDINCS) \
	-DTAUCS_CORE_DOUBLE \
	$(DIRSRC)taucs_ccs_ooc_llt.c \
	$(COUTFLG)$(DIROBJ)taucs_ccs_ooc_llt_D$(OBJEXT)
$(DIROBJ)taucs_ccs_ooc_llt_S$(OBJEXT): $(DIRSRC)taucs_ccs_ooc_llt.c $(INCS) $(STDDEPS)
	$(CC) -c $(CFLAGS) $(STDDEFS) $(STDINCS) \
	-DTAUCS_CORE_SINGLE \
	$(DIRSRC)taucs_ccs_ooc_llt.c \
	$(COUTFLG)$(DIROBJ)taucs_ccs_ooc_llt_S$(OBJEXT)
$(DIROBJ)taucs_ccs_ooc_llt_Z$(OBJEXT): $(DIRSRC)taucs_ccs_ooc_llt.c $(INCS) $(STDDEPS)
	$(CC) -c $(CFLAGS) $(STDDEFS) $(STDINCS) \
	-DTAUCS_CORE_DCOMPLEX \
	$(DIRSRC)taucs_ccs_ooc_llt.c \
	$(COUTFLG)$(DIROBJ)taucs_ccs_ooc_llt_Z$(OBJEXT)
$(DIROBJ)taucs_ccs_ooc_llt_C$(OBJEXT): $(DIRSRC)taucs_ccs_ooc_llt.c $(INCS) $(STDDEPS)
	$(CC) -c $(CFLAGS) $(STDDEFS) $(STDINCS) \
	-DTAUCS_CORE_SCOMPLEX \
	$(DIRSRC)taucs_ccs_ooc_llt.c \
	$(COUTFLG)$(DIROBJ)taucs_ccs_ooc_llt_C$(OBJEXT)
$(DIROBJ)taucs_ccs_ooc_lu$(OBJEXT): $(DIRSRC)taucs_ccs_ooc_lu.c $(INCS) $(STDDEPS)
	$(CC) -c $(CFLAGS) $(STDDEFS) $(STDINCS) \
	-DTAUCS_CORE_GENERAL \
	$(DIRSRC)taucs_ccs_ooc_lu.c \
	$(COUTFLG)$(DIROBJ)taucs_ccs_ooc_lu$(OBJEXT)
$(DIROBJ)taucs_ccs_ooc_lu_D$(OBJEXT): $(DIRSRC)taucs_ccs_ooc_lu.c $(INCS) $(STDDEPS)
	$(CC) -c $(CFLAGS) $(STDDEFS) $(STDINCS) \
	-DTAUCS_CORE_DOUBLE \
	$(DIRSRC)taucs_ccs_ooc_lu.c \
	$(COUTFLG)$(DIROBJ)taucs_ccs_ooc_lu_D$(OBJEXT)
$(DIROBJ)taucs_ccs_ooc_lu_S$(OBJEXT): $(DIRSRC)taucs_ccs_ooc_lu.c $(INCS) $(STDDEPS)
	$(CC) -c $(CFLAGS) $(STDDEFS) $(STDINCS) \
	-DTAUCS_CORE_SINGLE \
	$(DIRSRC)taucs_ccs_ooc_lu.c \
	$(COUTFLG)$(DIROBJ)taucs_ccs_ooc_lu_S$(OBJEXT)
$(DIROBJ)taucs_ccs_ooc_lu_Z$(OBJEXT): $(DIRSRC)taucs_ccs_ooc_lu.c $(INCS) $(STDDEPS)
	$(CC) -c $(CFLAGS) $(STDDEFS) $(STDINCS) \
	-DTAUCS_CORE_DCOMPLEX \
	$(DIRSRC)taucs_ccs_ooc_lu.c \
	$(COUTFLG)$(DIROBJ)taucs_ccs_ooc_lu_Z$(OBJEXT)
$(DIROBJ)taucs_ccs_ooc_lu_C$(OBJEXT): $(DIRSRC)taucs_ccs_ooc_lu.c $(INCS) $(STDDEPS)
	$(CC) -c $(CFLAGS) $(STDDEFS) $(STDINCS) \
	-DTAUCS_CORE_SCOMPLEX \
	$(DIRSRC)taucs_ccs_ooc_lu.c \
	$(COUTFLG)$(DIROBJ)taucs_ccs_ooc_lu_C$(OBJEXT)
$(DIROBJ)taucs_ooc_io$(OBJEXT): $(DIRSRC)taucs_ooc_io.c $(INCS) $(STDDEPS)
	$(CC) -c $(CFLAGS) $(STDDEFS) $(STDINCS) \
	-DTAUCS_CORE_GENERAL \
	$(DIRSRC)taucs_ooc_io.c \
	$(COUTFLG)$(DIROBJ)taucs_ooc_io$(OBJEXT)

$(DIROBJ)taucs_iter_D$(OBJEXT): $(DIRSRC)taucs_iter.c $(INCS) $(STDDEPS)
	$(CC) -c $(CFLAGS) $(STDDEFS) $(STDINCS) \
	-DTAUCS_CORE_DOUBLE \
	$(DIRSRC)taucs_iter.c \
	$(COUTFLG)$(DIROBJ)taucs_iter_D$(OBJEXT)

$(DIROBJ)taucs_vaidya_D$(OBJEXT): $(DIRSRC)taucs_vaidya.c $(INCS) $(STDDEPS)
	$(CC) -c $(CFLAGS) $(STDDEFS) $(STDINCS) \
	-DTAUCS_CORE_DOUBLE \
	$(DIRSRC)taucs_vaidya.c \
	$(COUTFLG)$(DIROBJ)taucs_vaidya_D$(OBJEXT)

$(DIROBJ)taucs_recvaidya_D$(OBJEXT): $(DIRSRC)taucs_recvaidya.c $(INCS) $(STDDEPS)
	$(CC) -c $(CFLAGS) $(STDDEFS) $(STDINCS) \
	-DTAUCS_CORE_DOUBLE \
	$(DIRSRC)taucs_recvaidya.c \
	$(COUTFLG)$(DIROBJ)taucs_recvaidya_D$(OBJEXT)

$(DIROBJ)taucs_gremban_D$(OBJEXT): $(DIRSRC)taucs_gremban.c $(INCS) $(STDDEPS)
	$(CC) -c $(CFLAGS) $(STDDEFS) $(STDINCS) \
	-DTAUCS_CORE_DOUBLE \
	$(DIRSRC)taucs_gremban.c \
	$(COUTFLG)$(DIROBJ)taucs_gremban_D$(OBJEXT)

$(DIROBJ)taucs_ccs_xxt_D$(OBJEXT): $(DIRSRC)taucs_ccs_xxt.c $(INCS) $(STDDEPS)
	$(CC) -c $(CFLAGS) $(STDDEFS) $(STDINCS) \
	-DTAUCS_CORE_DOUBLE \
	$(DIRSRC)taucs_ccs_xxt.c \
	$(COUTFLG)$(DIROBJ)taucs_ccs_xxt_D$(OBJEXT)

$(DIROBJ)taucs_ccs_generators_D$(OBJEXT): $(DIRSRC)taucs_ccs_generators.c $(INCS) $(STDDEPS)
	$(CC) -c $(CFLAGS) $(STDDEFS) $(STDINCS) \
	-DTAUCS_CORE_DOUBLE \
	$(DIRSRC)taucs_ccs_generators.c \
	$(COUTFLG)$(DIROBJ)taucs_ccs_generators_D$(OBJEXT)

$(DIROBJ)taucs_malloc$(OBJEXT): $(DIRSRC)taucs_malloc.c $(INCS) $(STDDEPS)
	$(CC) -c $(CFLAGS) $(STDDEFS) $(STDINCS) \
	-DTAUCS_CORE_GENERAL \
	$(DIRSRC)taucs_malloc.c \
	$(COUTFLG)$(DIROBJ)taucs_malloc$(OBJEXT)

$(DIROBJ)readhb$(OBJEXT): $(DIREXTSRC)readhb$(F2CEXT) $(STDDEPS)
	$(FC) -c $(FFLAGS) \
	$(DIREXTSRC)readhb$(F2CEXT) \
	$(FOUTFLG)$(DIROBJ)readhb$(OBJEXT)

$(DIROBJ)amdatr$(OBJEXT): $(DIREXTSRC)amdatr$(F2CEXT) $(STDDEPS)
	$(FC) -c $(FFLAGS) \
	$(DIREXTSRC)amdatr$(F2CEXT) \
	$(FOUTFLG)$(DIROBJ)amdatr$(OBJEXT)

$(DIROBJ)amdbar$(OBJEXT): $(DIREXTSRC)amdbar$(F2CEXT) $(STDDEPS)
	$(FC) -c $(FFLAGS) \
	$(DIREXTSRC)amdbar$(F2CEXT) \
	$(FOUTFLG)$(DIROBJ)amdbar$(OBJEXT)

$(DIROBJ)amdNEW$(OBJEXT): $(DIREXTSRC)amdNEW$(F2CEXT) $(STDDEPS)
	$(FC) -c $(FFLAGS) -O0 \
	$(DIREXTSRC)amdNEW$(F2CEXT) \
	$(FOUTFLG)$(DIROBJ)amdNEW$(OBJEXT)

$(DIROBJ)amdbarNEW$(OBJEXT): $(DIREXTSRC)amdbarNEW$(F2CEXT) $(STDDEPS)
	$(FC) -c $(FFLAGS) \
	$(DIREXTSRC)amdbarNEW$(F2CEXT) \
	$(FOUTFLG)$(DIROBJ)amdbarNEW$(OBJEXT)

$(DIROBJ)amdexa$(OBJEXT): $(DIREXTSRC)amdexa$(F2CEXT) $(STDDEPS)
	$(FC) -c $(FFLAGS) \
	$(DIREXTSRC)amdexa$(F2CEXT) \
	$(FOUTFLG)$(DIROBJ)amdexa$(OBJEXT)

$(DIROBJ)amdhaf$(OBJEXT): $(DIREXTSRC)amdhaf$(F2CEXT) $(STDDEPS)
	$(FC) -c $(FFLAGS) \
	$(DIREXTSRC)amdhaf$(F2CEXT) \
	$(FOUTFLG)$(DIROBJ)amdhaf$(OBJEXT)

$(DIROBJ)amdhat$(OBJEXT): $(DIREXTSRC)amdhat$(F2CEXT) $(STDDEPS)
	$(FC) -c $(FFLAGS) \
	$(DIREXTSRC)amdhat$(F2CEXT) \
	$(FOUTFLG)$(DIROBJ)amdhat$(OBJEXT)

$(DIROBJ)amdpre$(OBJEXT): $(DIREXTSRC)amdpre$(F2CEXT) $(STDDEPS)
	$(FC) -c $(FFLAGS) \
	$(DIREXTSRC)amdpre$(F2CEXT) \
	$(FOUTFLG)$(DIROBJ)amdpre$(OBJEXT)

$(DIROBJ)amdtru$(OBJEXT): $(DIREXTSRC)amdtru$(F2CEXT) $(STDDEPS)
	$(FC) -c $(FFLAGS) \
	$(DIREXTSRC)amdtru$(F2CEXT) \
	$(FOUTFLG)$(DIROBJ)amdtru$(OBJEXT)

$(DIROBJ)genmmd$(OBJEXT): $(DIREXTSRC)genmmd$(F2CEXT) $(STDDEPS)
	$(FC) -c $(FFLAGS) \
	$(DIREXTSRC)genmmd$(F2CEXT) \
	$(FOUTFLG)$(DIROBJ)genmmd$(OBJEXT)

$(DIROBJ)colamd$(OBJEXT): $(DIREXTSRC)colamd.c $(INCS) $(STDDEPS)
	$(CC) -c $(CFLAGS) $(STDDEFS) $(STDINCS) \
	 $(DIREXTSRC)colamd.c \
	$(COUTFLG)$(DIROBJ)colamd$(OBJEXT)

$(DIROBJ)direct$(OBJEXT): $(DIRPROGS)direct.c $(INCS) $(STDDEPS)
	$(CC) -c $(CFLAGS) $(STDDEFS) $(STDINCS) \
	 $(DIRPROGS)direct.c \
	$(COUTFLG)$(DIROBJ)direct$(OBJEXT)

$(DIROBJ)taucs_run$(OBJEXT): $(DIRPROGS)taucs_run.c $(INCS) $(STDDEPS)
	$(CC) -c $(CFLAGS) $(STDDEFS) $(STDINCS) \
	 $(DIRPROGS)taucs_run.c \
	$(COUTFLG)$(DIROBJ)taucs_run$(OBJEXT)

$(DIROBJ)iter$(OBJEXT): $(DIRPROGS)iter.c $(INCS) $(STDDEPS)
	$(CC) -c $(CFLAGS) $(STDDEFS) $(STDINCS) \
	 $(DIRPROGS)iter.c \
	$(COUTFLG)$(DIROBJ)iter$(OBJEXT)

libtaucs_content = \
  $(DIROBJ)taucs_sn_llt$(OBJEXT) \
  $(DIROBJ)taucs_sn_llt_D$(OBJEXT) \
  $(DIROBJ)taucs_linsolve$(OBJEXT) \
  $(DIROBJ)taucs_logging$(OBJEXT) \
  $(DIROBJ)taucs_timer$(OBJEXT) \
  $(DIROBJ)taucs_ccs_base$(OBJEXT) \
  $(DIROBJ)taucs_ccs_base_D$(OBJEXT) \
  $(DIROBJ)taucs_vec_base$(OBJEXT) \
  $(DIROBJ)taucs_vec_base_D$(OBJEXT) \
  $(DIROBJ)taucs_ccs_ops$(OBJEXT) \
  $(DIROBJ)taucs_ccs_ops_D$(OBJEXT) \
  $(DIROBJ)taucs_ccs_order$(OBJEXT) \
  $(DIROBJ)taucs_ccs_factor_llt$(OBJEXT) \
  $(DIROBJ)taucs_ccs_factor_llt_D$(OBJEXT) \
  $(DIROBJ)taucs_ccs_solve_llt$(OBJEXT) \
  $(DIROBJ)taucs_ccs_solve_llt_D$(OBJEXT) \
  $(DIROBJ)taucs_complex$(OBJEXT) \
  $(DIROBJ)taucs_complex_D$(OBJEXT) \
  $(DIROBJ)taucs_ccs_xxt_D$(OBJEXT) \
  $(DIROBJ)taucs_ccs_generators_D$(OBJEXT) \
  $(DIROBJ)taucs_malloc$(OBJEXT) \
  $(DIROBJ)readhb$(OBJEXT) \
  $(DIROBJ)amdatr$(OBJEXT) \
  $(DIROBJ)amdbar$(OBJEXT) \
  $(DIROBJ)amdbarNEW$(OBJEXT) \
  $(DIROBJ)amdNEW$(OBJEXT) \
  $(DIROBJ)amdexa$(OBJEXT) \
  $(DIROBJ)amdhaf$(OBJEXT) \
  $(DIROBJ)amdhat$(OBJEXT) \
  $(DIROBJ)amdpre$(OBJEXT) \
  $(DIROBJ)amdtru$(OBJEXT) \
  $(DIROBJ)genmmd$(OBJEXT) \
  $(DIROBJ)colamd$(OBJEXT) \
#  $(DIROBJ)taucs_memory$(OBJEXT) \
#  $(DIROBJ)taucs_ccs_io$(OBJEXT) \
#  $(DIROBJ)taucs_ccs_io_D$(OBJEXT) \
#  $(DIROBJ)taucs_ccs_io_S$(OBJEXT) \
#  $(DIROBJ)taucs_ccs_io_Z$(OBJEXT) \
#  $(DIROBJ)taucs_ccs_io_C$(OBJEXT) \
#  $(DIROBJ)taucs_ccs_ooc_llt$(OBJEXT) \
#  $(DIROBJ)taucs_ccs_ooc_llt_D$(OBJEXT) \
#  $(DIROBJ)taucs_ccs_ooc_llt_S$(OBJEXT) \
#  $(DIROBJ)taucs_ccs_ooc_llt_Z$(OBJEXT) \
#  $(DIROBJ)taucs_ccs_ooc_llt_C$(OBJEXT) \
#  $(DIROBJ)taucs_ccs_ooc_lu$(OBJEXT) \
#  $(DIROBJ)taucs_ccs_ooc_lu_D$(OBJEXT) \
#  $(DIROBJ)taucs_ccs_ooc_lu_S$(OBJEXT) \
#  $(DIROBJ)taucs_ccs_ooc_lu_Z$(OBJEXT) \
#  $(DIROBJ)taucs_ccs_ooc_lu_C$(OBJEXT) \
#  $(DIROBJ)taucs_ooc_io$(OBJEXT) \
#  $(DIROBJ)taucs_iter_D$(OBJEXT) \
#  $(DIROBJ)taucs_vaidya_D$(OBJEXT) \
#  $(DIROBJ)taucs_recvaidya_D$(OBJEXT) \
#  $(DIROBJ)taucs_gremban_D$(OBJEXT) \

$(DIRLIB)libtaucs$(LIBEXT): $(libtaucs_content) $(STDDEPS)
	- $(RM) $(DIRLIB)libtaucs$(LIBEXT)
	$(AR) $(AOUTFLG)$(DIRLIB)libtaucs$(LIBEXT) $(libtaucs_content)
	$(RANLIB) $(DIRLIB)libtaucs$(LIBEXT)
#$(DIREXE)direct$(EXEEXT): $(STDDEPS) $(DIROBJ)direct$(OBJEXT) $(DIRLIB)libtaucs$(LIBEXT)
#	$(LD) $(LDFLAGS) \
#	$(LOUTFLG)$(DIREXE)direct$(EXEEXT) \
#	$(DIROBJ)direct$(OBJEXT) \
#	-L$(DIRLIB) -ltaucs \
#	$(LIBS)
#$(DIREXE)taucs_run$(EXEEXT): $(STDDEPS) $(DIROBJ)taucs_run$(OBJEXT) $(DIRLIB)libtaucs$(LIBEXT)
#	$(LD) $(LDFLAGS) \
#	$(LOUTFLG)$(DIREXE)taucs_run$(EXEEXT) \
#	$(DIROBJ)taucs_run$(OBJEXT) \
#	-L$(DIRLIB) -ltaucs \
#	$(LIBS)
#$(DIREXE)iter$(EXEEXT): $(STDDEPS) $(DIROBJ)iter$(OBJEXT) $(DIRLIB)libtaucs$(LIBEXT)
#	$(LD) $(LDFLAGS) \
#	$(LOUTFLG)$(DIREXE)iter$(EXEEXT) \
#	$(DIROBJ)iter$(OBJEXT) \
#	-L$(DIRLIB) -ltaucs \
#	$(LIBS)
#all: $(STD_PRE_TARGETS) $(DIRLIB)libtaucs$(LIBEXT) $(DIREXE)direct$(EXEEXT) $(DIREXE)taucs_run$(EXEEXT) $(DIREXE)iter$(EXEEXT)
all:$(DIRLIB)libtaucs$(LIBEXT)
