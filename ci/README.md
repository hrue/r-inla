# How these builds drive the makefiles

The makefiles under `gmrflib/` and `inlaprog/` hold the build *rules* only.
Everything machine-specific (compilers, paths, library choices, modes) is a
*value*, and values are supplied by whoever runs make. GNU make lets a
command-line assignment override anything the makefile sets, so one set of
makefiles builds every platform and every configuration without editing:

    make -C gmrflib CC=$TRIPLET-gcc CXX=$TRIPLET-g++ FC=$TRIPLET-gfortran

That single idea is why the CI never patches a makefile: `ci/build.sh`,
`ci/build-macos.sh` and `ci/build-mingw.sh` compute the values per lane and
pass them in. The same makefiles produce the Linux bundles, the macOS
bundles and the Windows cross build.

## Where each kind of value lives

| kind | lives in | example |
|---|---|---|
| versions and pins | `ci/toolchain.env` (the single source of truth) | `GCC_PREFER`, `MINGW_TRIPLET`, `OPENBLAS_VERSION` |
| per-lane wiring | the `ci/build-*.sh` script | `RLIB_INC`, `RLIB_LIB`, `EXTLIBS2` overrides |
| build rules | the makefiles, unchanged | compile and link recipes |
| mode switches | compile-time defines in the source | `-DINLA_WITH_LIBR_DLOPEN`, `-DINLA_WITH_STILES` |

The last row matters: when a *behavior* changes (how R is embedded, whether
sTiles is present), the switch is a define handled in the C source, not a
different link recipe per makefile. The makefile never learns which mode it
is building.

## Two worked examples

**External packages.** Instead of naming each archive, glob the directory:

    EPATH = $(PREFIXROOT)/external-packages-compile
    ELIBS = $(shell echo $(EPATH)/lib*.a)

A new package only has to build into that directory; no makefile changes,
ever again. (`$(shell echo ...)` rather than `$(wildcard ...)` on purpose:
with no archives present, the literal `lib*.a` stays on the link line and the
link fails loudly, where `wildcard` would silently link an inla with no
external packages.)

**libR.** Deriving the values instead of writing them makes R changes free
the same way:

    R_HOME   = $(shell R RHOME)
    RLIB_INC = -DINLA_WITH_LIBR -I$(shell "$(R_HOME)/bin/Rscript" -e 'cat(R.home("include"))')
    RLIB_LIB = -lRmath -L$(R_HOME)/lib -lR

After that one edit, a new R version or layout needs nothing; switching the
embedding mode is a define (`ci/build.sh` builds all three libR modes,
linked, dlopen and none, from the unchanged makefiles this way).

## The rules that keep it working

- Never write a version, path or triplet into a makefile or a CI script when
  it can live in `ci/toolchain.env` or be derived at run time.
- Never branch on platform inside a makefile when the caller can pass the
  right value instead.
- When a change would mean "edit every makefile", it is in the wrong layer:
  move it to a variable the caller supplies, a derived value, or a define.
