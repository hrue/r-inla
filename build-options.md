This document discuss various options available to compile and build
the inla-program using the gcc-suite.

The compile-options has normally the format

```gcc -DINLA_WITH_FEATURE ...```

allowing ```FEATURE``` to be enabled this ```flag``` is set, and if
this ```flag``` is not set then this feature is off. In one instance, a
feature is default on and ```flag``` can then turn it off,

```gcc -DINLA_WITHOUT_FAST_GAMMA ...```

The main features are as follows (some developer options are not
shown):

- INLA_WITH_ARMPL Add this ```flag``` if you link with [ARM
  Performance Libraries](https://developer.arm.com/tools-and-software/arm-performance-libraries):
  Link with **-L$(ARMPL_PATH)/lib -larmpl -lamath -lastring**

- INLA_WITH_CLONE_TARGETS This ```flag`` will enable multi-target
  builds, like building for ```AVX2```-support and ```AVX512F```-support in the
  CPU, which is determined at run-time. Extentions are defined for ```Mac
  Intel``` and ```M-series```, and Linux on ```AARCH64```.

- INLA_WITH_COLOR This ```flag``` will enable some error-messages,
  warnings-messages and some other messages, to be shown in colors.
  (This feature is experimental.)

- INLA_WITH_EXTERNAL_PACKAGES This ```flag``` will enable the
  inclusion of code from some selected external repositories, like
  **rSPDE**, **INLAtools**, **INLAspacetime**, etc, which provides
  extended features, new models and likelihoods. These repositories
  needs to be compiled spesifically, and the example **Makefile** will
  show how to set this up.

- INLA_WITH_FRAMEWORK_ACCELERATE Add this ```flag``` to enable and
  link with the **Accelerate Framework** in **MacOS**

- INLA_WITH_HWLOC Add this ```flag``` to enable and link with the
  **libhwloc**-library, a `Portable Hardware Locality (hwloc) is an
  open-source, BSD-licensed software project that provides a portable
  abstraction of the hierarchical topology of modern computing
  architectures'.

- INLA_WITH_LIBR Add this ```flag``` to link with **libR** library,
  which enables some **R**-like-features: **rgeneric** and **rprior**.
  The **libR** library, can be linked wither explicite or implicite,
  depending on if the **INLA_WITH_LIBR_DLOPEN**-flag is set or not

- INLA_WITH_LIBR_DLOPEN Enable this ```flag``` to link with **libR**
  implicite, meaning that **libR** is read and loaded dynamically
  entirely from within the inla-program. This has the advantage that
  the build is not tied up to a spesific version of **libR**, but we
  can use the version that is already in use. This require linking
  with the **libltdl**-library (``libltdl (Libtool Dynamic Module
  Loader) is a small, portable library part of GNU Libtool that hides
  the differences in dynamic library loading across different
  operating systems'').

- INLA_WITH_MKL Enable this ```flag`` when linking with the **Intel
  MKL library** (``The Intel oneAPI Math Kernel Library (oneMKL) is a
  computing math library of highly optimized and extensively
  parallelized routines for applications that require maximum
  performance''). We recommend linking with the sequential version (as
  the parallelisation in the inla-program is done on a higher level):
  **-L$(MKL_PATH)/lib -lmkl_gf_lp64 -lmkl_sequential -lmkl_core**

- INLA_WITH_MUPARSER Enable this ```flag`` to link with the
  **libmuparser**-library, which is required for **expression**-priors
  and some other minor features in the inla-program.

- INLA_WITH_MVEC Enable this ```flag``` to use the vectorized version
  of basic numerical functions (**exp**, **log**, etc) from the
  **libmvec** library (in the gcc-suite) instead of using the
  corresponding versions in **MKL**, **ARMPL** and **ACCELERATE
  FRAMEWORK**

- INLA_WITH_NUMA Enable this ```flag``` to link with the **libnuma**
  library (Linux only): ``The NUMA library (libnuma) provides a
  high-level programming interface that allows developers to manage
  Non-Uniform Memory Access (NUMA) policies''. This feature is
  experimental.

- INLA_WITH_OPENBLAS Add this ```flag``` when linking with the
  **libopenblas**-library: ``OpenBLAS is an optimized BLAS (Basic
  Linear Algebra Subprograms) library based on GotoBLAS2 1.13 BSD
  version''.

- INLA_WITH_OPENSSL Add to ```flag``` to link with the SHA256 hash
  algorithm implementation in **libcrypt** (in the **openssl**
  package). This implementation is faster than the built-in version
  that will be used if this ```flag``` is not set.

- INLA_WITH_SIMDE Enable this ```flag``` to enable SIMDE (header-only
  library) which provides a faster implementation of some critial
  low-level functions: **SIMDE** ``provides fast, portable
  implementations of SIMD intrinsics on hardware which doesn't
  natively support them, such as calling SSE functions on ARM. There
  is no performance penalty if the hardware supports the native
  implementation (e.g., SSE/AVX runs at full speed on x86, NEON on
  ARM, etc.)''.

- INLA_WITH_STILES Enable this ```flag``` to link with the
  **libsTiles** library: ``sTiles targets the full Cholesky → solve →
  selected inverse pipeline on a single shared-memory node, using
  tile-based parallelism.''

- INLA_WITHOUT_FAST_LGAMMA This feature is default enabled, but can be
  turned off if this ```flag``` is set. By default, the inla-program
  use a faster implementation (that also improve the ``SIMD``-speedup)
  of the **lgamma** function than is provided in **libm**, at the
  slight cost of accuracy.


