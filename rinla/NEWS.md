# INLA 24.10.29
* Fixed a bug when using family 'egp'
* Update vignette 'SPDE2mdatasurv'
* Code cleanup and improvement

# INLA 24.10.13
* Changes in family binomialmix
* Code optimization for Mac Intel/ARM
* Code cleanup and optimization

# INLA 24.10.07
* Some minor changes in family=egp and binomialmix
* Fix a constant in barrier model (see diff for details)
* Code cleanup

# INLA 24.10.02
* New family: 'egp' (experimental)
* Minor fixes and improvements

# INLA 24.09.29
* Code improvement and optimisation
* Fixed a regression for Intel Mac run-script

# INLA 24.09.26
* Link with armpl-library for Mac-ARM64
* New experimental likelihood 'binomialmix'
* Various code improvements

# INLA 24.09.14
* Increased number of parameters for model 'intslope'
* Various code improvements

# INLA 24.06.29
* Fixed a regression from Version 24.06.09
* Code optimization and development

# INLA 24.06.19
* Code improvement 
* Fixed regression using `inla.group.cv`
* Better warning for utf8 characters in home-directory

# INLA 24.06.14
* More work on `control.mode`
* Fixed an issue with `A.local`

# INLA 24.06.08
* Documentation fixes
* More robust code for ``survival` models
* Remove inla.mode="twostage"
* Fixed regression in `control.mode`

# INLA 24.06.04
* Support for `inla.mdata`-objects in `inla.stack` (experimental)
* New likelihood `exppower` (experimental)

# INLA 24.06.02
* Fixed a regressin in model `iidkd`
* Fixed missprints in `group cv` vignette
* Updated `plot`
* Added new experimental likelihood

# INLA 24.05.27
* Change: do not load INLA-package by default using `rgeneric` or `rprior`
* New (experimental) likelihood: `occupancy`
* Optimized family `betabinomial`
* Minor code work

# INLA 24.05.22
* Increased max.dim for model `iidkd` from 20 to 24
* Optimized family `betabinomial`

# INLA 24.05.18
* Fixed (I hope) PARDISO issue for Mac (Intel and M...) due to OpenMP
  issues
* New and smaller default value for `tolerance.step`
* Fixed a regression in `cgeneric.h` in cache-macros
* Revise `plot(..., plot.opt.trace=TRUE)` plot

# INLA 24.05.13
* Rewrote 'read.graph' to avoid recursion

# INLA 24.05.10
* First build with R-4.4

# INLA 24.05.01
* New likelihood `mgamma`
* Fixes in `inla.knmodels`

# INLA 24.04.25
* Fixed a regression in VB correction for the mean
* Updated vignette (barrier global)
* Minor code changes

# INLA 24.04.22
* New likelihood 'npoisson'
* Robustify VB mean correction
* Various updates

# INLA 24.04.01
* New likelihood `stochvolln`
* Updated vignette `barrier_global`
* More code-cleanup

# INLA 24.03.29
* Limit dimension of the low-rank VB correction when replications
  and/or group is very high
* Code optimization and cleanup

# INLA 24.03.24
* Code optimization

# INLA 24.03.20
* Fixed issue in `rprior`
* More development for `gcpo`
* Update vignette about CODA

# INLA 24.03.18
* gcpo development (new features) and code improvement
* likelihood `fl` changes
* `MatrixModels` moved to Imports

# INLA 24.03.09
* Code optimization and improvement

# INLA 24.03.02
* Updated PARDISO library to version 8.2
* Various code improvements

# INLA 24.02.27
* Fixed an regression since cf5714584c6ea59911526fb3952193bfb9f550c6
  for the marginal of one hyperparameter only
* Improved function `inla.spy`
* Minor code edits

# INLA 24.02.23
* Code optimization
* Remove `experimental` status of many likelihoods and model
  components
* Fixed an issue with `scopy`

# INLA 24.02.21
* New vignette added (barrier model on the world's oceans)
* Small fix for 'malloc.lib' for Mac's

# INLA 24.02.20
* Code cleanup and improvements
* New option 'malloc.lib' in 'inla.setOption' (experimental)

# INLA 24.02.09
* Fixed a bug: prior means was not accounted for, doing the VB mean correction
* Added new (experimental) likelihood: `tpoisson`
* Some code cleanup

# INLA 24.02.07
* Fixed a bug: prior means for fixed effects was not included for
  compact mode
* New experimental likelihoods: `rcpoisson` and `fl`
* New arguments to `control.vb$f.enable.limit` to limit
  maximum dimension of vb-corrections (which is multiplied with
  `replicates` and `group`)
* Code cleanup

# INLA 24.01.29
* Minor bug fix
* Improved control of rounding errors for numerical derivatives
* New experimental link functions, `gev` and `cgec`

# INLA 24.01.20
* Code cleanup
* Remove library jemalloc in builds

# INLA 24.01.18
* Cleanup code for sparse dot-products
* Intel MKL no longer supports MacOS and the Accelerate framework is now
  used instead in the build
* Some code optimization

# INLA 24.01.14
* Optimization of sparse dot-products
* Code cleanup

# INLA 24.01.09
* Fixe issue with 'r$mode$x'

# INLA 24.01.07
* Some code improvement and minor fixes
* Fixed false warning for 'scopy'
* Fixed issue for 'safe'-option

# INLA 23.12.17
* New version of `scopy` (experimental)
* Added support for defining priors in `R` (experimental, see `rprior`)
* Some code improvements
* Updated vignette about `group-cv`

# INLA 23.12.06
* Fixed issue (hopefully) with saturated DIC
* New option, 'internal.opt' in 'inla.setOption'
* Fixed regression error: 'result$misc$nfunc'
* Made link 'cloglog' and 'ccloglog' more numerical stable
* Link with PARDISO8-library (not yet Windows), get license at 'https://panua.ch'
* Switch to static MKL libraries for Mac Intel


# INLA 23.11.26
* New (experimental) likelihood `ggaussian` and `ggaussianS`
* Fix for `readlink` on older MacOS

# INLA 23.11.15
* Development build

# INLA 23.11.12
* Development build

# INLA 23.11.01
* Fix for sparse-matrix in data-argument
* Refactor internal functions

# INLA 23.10.28
* Cleanup after a merge done by error

# INLA 23.10.27
* Improve storage when `selection` is used
* Various minor development

# INLA 23.10.25
* Development bug-fix
* Fixed memory leaks


# INLA 23.10.19
* Bux-fix to work around regression bug in package 'excursions'
* Internal development 

# INLA 23.10.17
* Fixed issue with `inla.call="submit"`
* Increased precision for `copy`

# INLA 23.10.14
* New family `stdgaussian` (where the precision is fixed to be 1)
* New family `nzPoisson` (for Poisson without zero's)
* Internal code-cleanup and some minor fixes

# INLA 23.09.09
* A new approach to recover from failure (option `safe=TRUE`)
* Option `control.inla=list(b.strategy="keep")` is a new default.
* Various code improvement

# INLA 23.08.26
* Bug-fix release (see commit log)
* New vignette added

# INLA 23.08.18
* Some improvement in the optimiser
* Fix option `control.compute=list(q=TRUE)` for the default mode
* Depend on the new `fmesher` package. By default, the new methods are used
  silently instead of the old fmesher standalone binary, and are meant to produce
  the same meshes etc as before. During a transition period, one can switch between
  the two code bases, as well as turn on informative deprecation warnings
  that point to which fmesher R function replaces the existing ones.
  See https://inlabru-org.github.io/fmesher/articles/inla_conversion.html
  for more details.

# INLA 23.08.08
* Fixed regression with argument `cdf=`. 
* Various code improvement and minor fixes.

# INLA 23.06.29
* Some code rewrite
* Fixed an issue if joint-priors (`control.expert`) for R-4.3.

# INLA 23.06.25
* Fix vignettes

# INLA 23.06.24
* Some improvement in TAUCS linear solve(s)
* Fixed an issue with Gaussian likelihood and non-identity link
* Fixed issue with priors with no parameters and `inla.rerun`
* Fixed an issue with computing moments from a marginal
* Update deprecation information (sp/rgdal)
* Remove CYGWIN code
* Various namespace and inla/inlabru connection updates
* Documentation updates to protect the user's filespace
* Various code-cleanup work


# INLA 23.06.15
* Tolerance was set too low
* Robustify deviance residuals

# INLA 23.06.14
* Improved version of `inla.rerun`
* Issues reported in `23.06.12` fixed

# INLA 23.06.12
* Documentation and dependency updates (roxygen transition is now complete)
* Removed AMDC from default reordering (TAUCS only)
* Internal code reorganisation 
* Optimisation of cpo/pit/po/dic
* Optimisation of spde2-models
* Optimisation of poisson/binomial/nbinomial (with default link and no offset)
* Various optimisation work
* New toolchain for building on Windows (MinGW)
* Removed openblas-related code (and specific builds)
* Removed `inla.setOption`-options: `mkl, blas.num.threads, vecLib, vecLibPath` and `CYGWIN`
* Added new external package: `fbesag`

# INLA 23.05.30
* First build with R-4.3 and MacOS 12.6
* Code improvement and optimisation
* Updated some documentation

# INLA 23.05.21
* Reduced memory usage for data-rich models
* Added new option `control.compute=list(save.memory=FALSE|TRUE)` for
  more aggressive savings
* Change how equal correlations are defined in `gcpo`
* Code optimization

# INLA 23.05.13
* Reduce memory consumption
* New diagnostic: `result$misc$warnings`

# INLA 23.05.10
* Added documentation for `copy`
* New `scopy` model added (experimental)
* Added option `hessian.correct.skewness.only` in `control.inla`
* Minor code changes

# INLA 23.04.24
* Fixed inconsistent treatment of `offset()`
  and `offset=` in the output, compared with 
  `mode='classic'`

# INLA 23.04.20
* Small fix in `inla.dryrun`
* Change to dynamic libs for Mac

# INLA 23.04.19
* Update config-output with vb-variance corrections
* Added `INLAspacetime` functions in `barrier.R`
* Added `inla.dryrun` (experimental)

# INLA 23.04.11
* Fixed an issue with `cgeneric`
* Dot-product optimization for MKL
* Improved Qinv for TAUCS
* Added global-constraint(s) to the expert option

# INLA 23.04.02
* Fixed an issue with group-cv vigette

# INLA 23.04.01
* Fixed issue with `model="slm"` (thx to RB and VGB)
* New vigette about CV and group-CV
* New option to turn off online-optimisations
* New option to measure gain in dot-product optimisations

# INLA 23.03.26
* Various code-cleanup and optimisation
* Improved code for `family="nbinom"` and `"bell"`
* New vignette (SPDE2mdatasurv)
* Small change in how saturated deviance is computed

# INLA 23.03.19
* Updated family=bell
* Fixed an issue with selection of nodes for vb.correct

# INLA 23.03.17
* Added family=bell (experimental)
* Workaround for PARDISO issue with diagonal matrix
* Various code-cleanup

# INLA 23.02.27
* New option to `inla.binary.install()`
* Fix a broken build

# INLA 23.02.22
* Add linear predictor summaries to config
* Vignette about `rgeneric` is updated
* Updated doc about `inla.posterior.sample`
* Fixed issue with missing covariance in A-matrix

# INLA 23.02.04
* Fix for `control.mode=list(..., fixed=TRUE)`
* Various code-cleanup and improvement

# INLA 23.01.31
* Remove `rgeos` dependency
* Remove `rgdal` dependency
* Add new link-function: `ccloglog`
* New argument to `inla.group.cv`
* `safe`-mode improved in the `inla()`-call
* Various code-cleanups

# INLA 23.01.12
* Fixed an issue with ptweedie
* Code cleanup

# INLA 23.01.02
* Sync C++ code from the fmesher package
* Added log-likelihood output to configs
* Fixed a regression for `model="fgn"`
* Improved code for `model="binomial"` and logit-link
* Increased kmax from 10 to 20 for `model="iidkd"`

# INLA 22.12.16
* Fixed an issue with lincombs and PARDISO library

# INLA 22.12.15
* Minor changes to integrate cgeneric-code better
* New stable build

# INLA 22.12.14
* Minor changes to integrate cgeneric-code better

# INLA 22.12.12
* Some minor code fixed
* Development build to test cgeneric update

# INLA 22.12.03
* Fixed an issue with cache and nested parallelism (mostly
  `cgeneric`-related) 

# INLA 22.11.28
* Added friends-list to gcpo
* Expanded support for non-unit-radius spheres in fmesher
* Experimental support for external-packages
* Experimental support for new 0inflated models
* Various code cleanup

# INLA 22.11.22
* Fix class inheritance checks to use inherits(), to allow e.g. inlabru
  output to work with `inla.cpo`, etc.
* Fixed issue with DIC in experimental-mode
* rename experimental-mode into compact-mode (experimental-mode will
  still work)
* Make compact-mode the default.   

# INLA 22.11.08
* bug-fix build

# INLA 22.11.07
* Added new vignette
* Added default values of `safe`, `verbose`' and `keep` to options
* Added new check for valid arguments (default off for the moment)
* `control.compute=list(residuals=TRUE)` provide expected deviance
  residuals, but it is experimental for the moment

# INLA 22.10.23
* Optimisation of gcpo and inner-products
* Port for upcoming Matrix 1.5-2 change in `kronecker()`
* Various minor code improvements

# INLA 22.10.15
* Bug fix

# INLA 22.10.14
* Internal changes and bug fix.

# INLA 22.10.13
* Bug fix

# INLA 22.10.12
* Added new experimental cure-model for parametric survival
  (documentation will come later).
* The `weibullcure`-model is removed as its covered by the new
  cure-model.
* Optimization: Models with many constraints 
* Optimization: Models with 'betabinomial' likelihood

# INLA 22.10.07
* Fixed regression in initial value for `model="copy"`
* Code cleanup

# INLA 22.10.06
* Rebuild package to due broken '22.10.05' version

# INLA 22.10.05
* Code cleanup and optimization of dot-products
* Fixed an issue with 'control.fixed=list(remove.names=...)'

# INLA 22.10.01
* MKL code for nested indexed inner products
* Computation of initial values for experimental mode
* Various more minor code improvements

# INLA 22.09.28
* Updated version of `fmesher`, the program which do the mesh-stuff
* Various code-improvement
* Some internal (work-in-progress) updates

# INLA 22.09.26
* Various code-improvement and cleanup
* Code-optimization of spde2-models
* New likelihood model `cennbinomial2` (experimental)
* New likelihood model `gaussianjw` (experimental)
* New function `inla.group.cv` (experimental)

# INLA 22.09.15
* Added support for `Predictor` and `APredictor` in lincomb for
  experimental mode
* Code improvement for `spde2` models
* Add support for ELLIPSOID radius get/set for epsg:4326
* Improved `gcpo` with singular covariance matrices
* Various minor tweaks

# INLA 22.09.02
* Make argument `f.enable.limit` work the replications and groups
* Fixed an issue with `vb` variance correction (experimental code)
* Added new vignette about `inla.posterior.sample.eval`
* Improve the openmp code for `gcpo`

# INLA 22.08.24
* Fixed an issue with the deprecated class coercions in `inla.as.dgTMatrix()`

# INLA 22.08.23
* Avoid deprecated class coercions in `inla.as.dgTMatrix()`
* Fixed an issue with `rgeneric` 
* New vignette added
* Various code improvements

# INLA 22.07.23
*  Fixed regression from `Version_22.07.15`

# INLA 22.07.21
*  Code-improvement and minor fixes
*  Fixed regression from `Version_22.07.15`

# INLA 22.07.15
*  Code-improvement and minor fixes

# INLA 22.07.01
*  Code-improvement/fix

# INLA 22.06.30
*  Code-improvement
*  Changed stopping criteria for VB corrections

# INLA 22.06.25
*  Code-improvement
*  More work on the new experimental feature (internal)

# INLA 22.06.20
*  Code-improvement
*  Adding new experimental feature (internal)

# INLA 22.06.11
*  Added more to `gcpo`-output
*  Code-improvement

# INLA 22.06.03
*  Bug-fix release

# INLA 22.06.02
*  The `mode` in the output summary is back.
*  Family `gev` is disabled, using `bgev` instead
*  Code improvement, cleanup and bug-fixes

# INLA 22.05.18
*  Changes in work-in-progress features

# INLA 22.05.07
*  Minor changes to move to R-4.2

# INLA 22.05.03
*  Code improvement, cleanup and bug-fixes
*  More work on the `gcpo`-feature
*  Removed threadprivate from code (almost), which was a huge rewrite.
*  Made linesearch in the optimization more robust and stable
*  Improved `stupid`-search
*  Internal changes in the `smart`-optimization
*  (This is likely the last build with R-4.1 before we switch to R-4.2)

# INLA 22.04.16
*  Added `theta`-correction for `gcpo` computations
*  Code improvement and cleanup

# INLA 22.04.14
*  Improve `gcpo`-handling of singular cases. `group.size=-1` now
    gives the CPO.

# INLA 22.04.13
*  Bug-fix release

# INLA 22.04.09
*  Bug-fix and code improvement release

# INLA 22.04.06
*  Bug-fix and code improvement release

# INLA 22.03.28
*  Bug-fix release

# INLA 22.03.27
*  Bug-fix release

# INLA 22.03.26
*  Bug-fix release

# INLA 22.03.25
*  Bug-fix release

# INLA 22.03.24
*  Improved the code for output which is now way faster for
      data-rich models
*  Internal changes in how integrals are computed and how
      marginals densities are represented. This which will
      make the results tiny tiny different from earlier versions.
*  Code-improvement for data-rich models
*  Minor bug-fixes
*  Fixed regression issue with `strategy="laplace"`
      which has been broken
      from version `22.02.16` (ok in `22.01.25`).
*  The value of `mode` in summary output is not longer computed
      for some densities (like for mixtures and non-linear
      transformation), as its pretty expensive to do so accurately.
      It is then shown as `NA`. This might be fixed later.

# INLA 22.03.16
*  Code cleanup and a bug-fix
*  Dense matrices are now store column-wise in `cgeneric`
*  New example added to `cgeneric`, see vignette
*  Changed default integration stratey for `experimental` mode

# INLA 22.03.12
*  Option `safe=TRUE` is now default

*  Code cleanup
*  Code optimization for `TACUS`
*  Nested parallelism is now enabled on Windows

# INLA 22.03.09
*  Issue with `gcpo` and `safe=T` fixed
*  Fixed issue with `A`-matrix in `gcpo`

# INLA 22.03.08
*  Better TAUCS performance: improved storage scheme and
      new code tp allow for linear solves with many rhs.
*  `gcpo`-improvement. The default is now `group.size=1`,
      and it acts like cpo.

*  Fixed an issue for propagation an error
      from the cpo-calculations
*  Various smaller changes

# INLA 22.03.06
*  Improved the TAUCS-interface performance.
*  Various code improvements.

# INLA 22.03.03
*  New options for the group-cpo feature.

*  Fixed CRAN-check issues.

# INLA 22.02.28
*  VB correction for experimental mode do now include nodes
      from all model components, not only the short ones and fixed
      effects.
*  New group-cpo feature, very experimental at the moment
      (internal use only).

*  Parts of the VB code has been improved for
      simplicity and speed.
*  `inla.call="remote"` now use `zstd` for
      parallel compression, meaning that `zstd` must be
      available on both sides.
*  PARDISO license help and startup msg has been disabled.
*  Minor fixes, code improvements and cleanup

# INLA 22.02.16
*  New feature `inla.cgeneric.q`
*  New and updated PARDISO library

*  Code improvements and cleanup

# INLA 22.01.25
*  Fixed regression bug in `knmodels`
*  Code improvements and cleanup

# INLA 22.01.19
*  Fixed bug in the variable expansion in some arguments
      that could give different results than before

# INLA 22.01.16
*  Added new variant for nbinomial likelihood

# INLA 22.01.12
*  Added `graph$cc$mean` to help defining intercepts for
      disconnected regions.

*  Can now define the license key for PARDISO directly as
      `inla.setOption(pardiso.license="<KEY>")`

*  Print total RAM in Gb in the inla-program output, for easier debugging
      of these problems.

*  Using `inla.call="remote"` will now return error code that can be
      controlled by `try`, so that `r=try(inla(..., inla.call="remote"))`
      should work as intended.

*  Add a refine step for `parallel.linesearch`, so if this is
      enabled, then it will do a restart at the end with this option
      disabled.

*  Minor changes when option `safe` is enabled.

*  Fixed regression error with `plot(result, plot.prior=TRUE)`
*  Fixed a regression error from version 21.12.21 which expanded some
      variables incorrect.

# INLA 21.12.21
*  Option `safe` in function `inla`

*  Minor bug fixes

# INLA 21.12.18
*  Improved handling of `NAN` and `INF`

# INLA 21.12.16
*  A bug-fix for `matern2d` model
*  Various code improvements

# INLA 21.12.11
*  Minor bug fixes for `cgeneric` for Windows

# INLA 21.12.10
*  A `cgeneric` interface
      like `rgeneric` but with C-code.
      This is for internal testing only.

*  Various code improvements and minor fixes

# INLA 21.12.01
*  Allow for probit cdf in likelihood model `pom`

*  Various code improvements

# INLA 21.11.29
*  Small fixes

# INLA 21.11.28
*  added new link `powerprobit`. Not completely done yet.

# INLA 21.11.22
*  Remove jemalloc-option that creates warning message for Mac

# INLA 21.11.21
*  Minor internal changes
*  Add jemalloc options for Linux/Mac

# INLA 21.11.16
*  Minor internal changes for robustness
*  Upgraded compilers to gcc-11 on Windows

# INLA 21.11.15
*  Option for parallel linesearch has been renamed into
      `parallel.linesearch` in `control.inla`.

*  This is fix for a build-error in the Windows binary for
      version `21.11.14`
*  Various minor internal code changes and improvements

# INLA 21.11.14
*  Parallel linesearch implemented, see option
       `bfgs.version` in `control.inla`. This option is highly
       experimental at the moment and work in progress, and should
       not be used.

*  Change default settings for `family="tweedie"`
*  Implemented an adaptive parallel/serial version of `Qx`
*  Improved perforance for parallel linear solve with `PARDISO`
*  Use stable Legendre polynomial evaluation for spherical covariances
*  Add support for `Fedora Linux 35`
*  Fixed `dic.sat` when all response are `NA`
*  Various internal code changes and improvements

# INLA 21.11.01
*  New experimental feature `lp.scale`. No good documentation yet, but
      that will come soon.

*  Robustify `inla.iidkd.sample` function for numerical singular matrices
*  Corrected `inla.stack.join` handling of factor variables that in some
      cases converted variables to/from factors unintentionally.
*  `inla.mode="experimental"` got its fitted values back.

# INLA 21.10.03
*  Internal code improvement

# INLA 21.09.29-1
*  Better special number

*  New option to `inla.iidkd.sample`

# INLA 21.09.29
*  New latent model `iidkd` which use a different
      parameterisation than `iid3d` and similar ones.

# INLA 21.09.13
*  Fixed crash when using `int.strategy="user"` with
      `inla.mode="experimental"`
*  Minor fixes

# INLA 21.09.11
*  Fixed crash when using linear combinations with
      `inla.mode="experimental"`
*  Minor fixes

# INLA 21.09.01
*  Added option `scale` to `family="xbinomial"`

# INLA 21.08.31
*  Various bug-fixes

# INLA 21.07.10
      
*  If `control.predictor=list(compute=TRUE)` then the
      marginal densities for the linear predictor and fitted values are
      computed, but the new default is that only the summary statistics
      are returned not the marginal densities. This is because the
      storage requirement for the marginals can be substantial. The
      marginal densities can be returned setting
      `control.compute=list(return.marginals.predictor=TRUE)`
      
*  The argument `control.results` are removed as it was
      not really in use.
      
*  The new default value for `b.strategy` is
      `b.strategy="skip"`. See `?control.inla` for details

*  New option `inla.timeout`, similar to
      `fmesher.timeout`; see `?inla.setOption`
      
*  `CPO-marginals` are returned in when
      `control.compute=list(config=TRUE)`
      
*  There is a change how a (numerical) singular constraint is
      treated, hopefully this new way is more robust.

  

*  New option `inla.mode`, which define how to arrange the
      internal model formulation. One of `"classic"`,
      `"twostage"` and `"experimental"`. The default is
      `"classic"`, which is unchanged behaviour compared for
      earlier versions. The other two are highly experimental for the
      moment. See `?inla`

  

*  Massive code clean-up and some minor fixes.

# INLA 21.06.11
*  Added support for Mac M1 (native build with R-4.1)
*  New `twostage` options (highly experimental)

*  Massive code overhaul
*  Various minor fixes

# INLA 21.05.02
*  Some minor bug-fixes and code improvements

# INLA 21.04.22
*  Mac only: Added optional path to `vecLib`
      BLAS and LAPACK libraries.
*  Adding argument `.special` to `inla.surv`-object.

*  Some code improvement 

# INLA 21.04.16
*  Mac only: Link with `vecLib` BLAS and LAPACK libraries
      by default. Turn off with `inla.setOption(vecLib=FALSE)`.
*  Adding an experimental improvement
      for `strategy=simplified.laplace` (not enabled by default)

*  Some code improvement 

# INLA 21.03.31
*  General code improvement 

# INLA 21.03.27
*  New family `stochvolsn`
*  New family `cenpoisson2`

*  Improve robustness for survival models
*  General code improvement

# INLA 21.03.21
*  Changes in `inla.jp`; see help-page for details.

# INLA 21.03.20
*  New family `gompertz`
*  `inla.binary.install()` now do md5-checksum check

*  Code optimization for `barrier` models

# INLA 21.03.17
*  Internal code changes for `family=tweedie`

# INLA 21.03.16
*  `inla.binary.install()` can now also run
            non-interactive

*  Internal code changes

# INLA 21.03.14
*  Internal code changes

# INLA 21.03.13
*  Fixed some missprints in doc
*  Trying to get vigettes to work again

# INLA 21.03.08
*  Added family `logperiodogram` back in
*  Added new family `agaussian`

*  Improved initial values for `x` 
*  Improved computations for saturated likelihood
*  Various minor code cleanup

# INLA 21.02.23
*  Default `strategy` is now `simplified.laplace`
      for smaller models (< 5000 nodes) and `adaptive` for larger ones

*  Minor bug-fixes

# INLA 21.01.26
*  Improve parallel rerodering for PARDISO
*  Improve `coxph` for large data
*  Improve internal CPU-timing output

# INLA 21.01.18
*  Use METIS5 with PARDISO

# INLA 21.01.13
*  Family `fmri`

*  Speed improvements for family `tweedie`

# INLA 21.01.08
*  New family `tweedie`

*  Optimization work, code improvement and maintenance

# INLA 20.12.10
*  Added a new experimental feature for running PARDISO in parallel

*  Various minor fixes and improvements

# INLA 20.12.04
*  Bug-fix release only

# INLA 20.11.29
*  Added new section to the vignette about `rgeneric`

*  Option `optimize` in `inla.rgeneric.define`

*  Cleanup in some output format.
*  Various internal changes. 

# INLA 20.11.22
*  Code improvement and optimization for `rgeneric`.
*  Some OpenMP improvements

# INLA 20.11.18
*  New family `gammajw`

*  Fixed an OpenMP deadlock case, using `rgeneric`

# INLA 20.11.16
*  PARDISO version 7 is included for Mac and Linux.
*  Improved paralellism and nested paralellism
*  Speedup improvements, especially for models with many
      constraints.
*  Improved default `plot(result)`. Argument
      `cex=..` will now work

*  New shorthand feature in `inla.posterior.sample.eval`
      for extracing model components from samples
*  Added `fmesher.timeout` option, see `?inla.setOption`

*  Using a non-zero seed in `inla.posterior.sample` will
      now force serial computations.
*  A lot of internal code cleanup and improvements

# INLA 20.10.11
*  New likelihood `poisson.special1`

*  Minor code cleanup
*  Minor changes in code/doc to adapt to R-4.0

# INLA 20.09.25
*  Bug fixes and some optimisation work

# INLA 20.09.16
*  Likelihood model `sn` is redone.

*  Various minor changes.

# INLA 20.09.13
*  Model `iid` added to possible models
      for the baseline hazard in the Cox-ph model

*  Likelihood models `sn` and `sn2` have been
      replaced with a new version `sn`. The parameterisation is
      different and now done correct.

# INLA 20.09.09
*  The likelihood models `sn` and `sn2`
      are now disabled. They need a rewrite to be done right.

*  For Linux, then build so-libraries are loaded before the
      system ones. This behaviour can be reverted by setting the
      environment variable `INLA_NATIVE_LD_LIBRARY_PATH`.

# INLA 20.09.07
*  Dimension for CCD integration have been
      increased from 38 to 52.

*  Default prior for the `intercept` and the `skewness` in the
      `skew-normal link` have changed. 

*  Various minor changes

# INLA 20.08.31
*  A small bug-fix with argument `-t`
*  Some improvements in the install scripts

# INLA 20.08.30
*  More work on the nested parallelism stuff

# INLA 20.08.29
*  Argument `num.threads="A"` now means
      `num.threads="A:1"`
*  Some changes for nested paralellism

# INLA 20.08.28
*  There is a minor change in `inla.qsample` and `inla.posterior.sample`.
      If argument `seed!=0` then serial mode is now forced.
      Earlier, an error would be raised if parallel mode was also requested.

*  Argument `remove.names` in `control.fixed`

*  Bug-fixes for nested parallism
*  Some optimization improvement

# INLA 20.08.24
*  Some minor fixes in the install

*  Linux builds now link with PARDISO version 7 (beta)

# INLA 20.08.22
*  MKL is now default for both Mac and Linux.
      If this cause any issues, revert back by setting
      `inla.setOption(mkl=FALSE)`.
*  As of today, the TAUCS library does not play 
      well with MKL on Mac for some reason, and this is
      silent accounted for.
      It is possible to bypass the internal check,
      and if that becomes in issue, just define
      `mkl=FALSE` as above.

*  Updated some install scripts

# INLA 20.08.21
*  Various minor fixes
*  MKL binaries are now linked with libiomp5

# INLA 20.08.19
*  Changed to the gcc-10 compiler suite for Ubuntu 1804
      and 2004.
*  Fixed an issue if the optimiser does not move
      and the use of directions
*  Minor fixes

# INLA 20.08.18
*  Code is prepared for nested parallism, which most gain
      obtained when using the PARDISO library
*  Argument `num.threads` is now in the format `A:B`,
      where A threads for the outer layer and B threads for the inner
      layer. And this applies for several functions.
*  The speed have improved, mostly when using the PARDISO
      library. Hopefully, nothing is broken.

*  likelihood `beta` now allow for cencoring near 0 and 1,
      see `inla.doc("beta", section="likelihood")`

*  Default value for `control.inla$h` is now 0.005, which
      replace the previous value of 0.01. We might change this later,
      but we have enough experience to know that 0.01 is slightly to
      large.
*  `control.inla$optimise.strategy="smart"` is now
      default. By construction, this should be both faster and more
      robust. Using `"plain"` reverts back the old behaviour. 
*  `control.inla$use.directions=TRUE` is now default. This option
      estimate numerical gradients/Hessian in directions with
      more change than for coordinate-wise directions.

# INLA 20.08.11
*  Add new experimental feature
      `control.inla=list(use.directions=TRUE/FALSE)`

# INLA 20.08.09
*  New experimental optimise strategy
      `control.inla=list(optmise.strategy="smart"`
      that is hopefully faster and as safe as the default one. Maybe
      even more safe and robust.

*  More work on the nested paralellism 

# INLA 20.08.06
*  More work on nested parallelism and some code changes
*  Code changes to reduce memory usage
*  Linux/Mac binaries are now linked with `jemalloc`

# INLA 20.08.04
*  Improved nested parallelism and some code changes

# INLA 20.08.03
*  More work on nested parallelism

# INLA 20.08.02
*  Option `control.inla$lincomb.derived.only` is now disabled.

*  Testing nested parallelism `openmp.nested` with
      `num.threads="A,B"`. Work in progress

# INLA 20.07.27
*  `plot(result)` will now produce a plot of the
      `CPO/PIT` for each likelihood family (if available), instead
      of a joint plot as earlier.
*  New vignette about `jmarginal`

*  Added new likelihood models `zeroinflatedcenpoisson0`
      and `zeroinflatedcenpoisson1`
*  Link-model `sn` is updated, as well as the
      PC-prior for the skewness therein, and the added intercept model.

*  Revised the PIT calculations for family `cenpoisson`
*  Code rewrite to (try to) prevent `Inf` for `DIC` calculations
*  Minor fixed in `inla.binary.install`

# INLA 20.07.18
*  Build-script changes and misc fixes

# INLA 20.07.16
*  Packages `mpoly` and `symmoments` are
      added to the Suggests-list.

*  Added info about `inla.prune()` to the
      startup message

*  More features for `jmarginal` added

*  Build-scripts fixes
*  Fix for a rare `fmesher` issue
*  Improved the code for the DIC calculations to make
      them more stable
*  Some improvment in the PC-prior for the SN-link

# INLA 20.07.12
*  Fixed `inla.link.sn` to vectorise over argument
      `a` as `sn`-package do not do that properly itself. 

# INLA 20.07.09
*  Package built with R-4.0

# INLA 20.07.04
*  Vignette added for `family=bGEV`
*  Some internal changes due to the migration to `git`

# INLA 20.06.29
*  Added `Qprior.diag` to the output when `config=TRUE`.
      The off-diagonals of this matrix are the same as `Q` in the
      same configuration, so only the diagonal of `Qprior`
      is stored.
*  Added some internal experimental code

*  PARDISO interface: internal check added
*  Fixed an wrong assert with family=bgev

# INLA 20.06.22
*  Improve some code in the PARDISO interface
*  Improved the computation of the third derivative
      in the log likelihood.

# INLA 20.06.18
*  Improved WKT support for PROJ6/GDAL3

# INLA 20.05.16
*  Support for PROJ6/RGDAL3 for handling CRS information for
        spatial objects.

# INLA 20.06.15
*  For the `intslope`-model: made all `gamma`'s
      default fixed to 1, so its similar in style the copy-feature.
*  Added argument `constr` it `inla.rjmarginal`
*  Added argument `ask` to function `inla.prune`

*  Argument `cyclic=TRUE` in `f()` should not set
            `constr=FALSE` when default is `constr=TRUE`
*  Change the `scale.model=TRUE` code for `RW1/RW2` so the
      scaling for the continous case is the same as for the discrete
      case when the locations are eqvidistant.
*  Disable link `sslogit`

# INLA 20.05.12
*  Fixed an issue with model `besag2`
*  Fixed an issue with `plot(r,plot.prior=TRUE)` for some priors

# INLA 20.05.04
*  Remove the experimental status of `inla.posterior.sample.eval`

*  Added function `inla.prune` which will remove binaries
      not supported by the running OS, to reduce the size of the
      package.
*  Added method `summary` and `print` to class
      `inla.jmarginal`

*  Add check for `NA/NaN/Inf` in mesh creation input
      ocations
*  Make sure that skewness is not to high in `inla.posterior.sample`

# INLA 20.04.18
*  Added new argument `tag` to `inla.coxph`

*  `inla.rjmarginal.eval`, to evaluate samples from a join
      approximations

*  Names of samples are now "sample:1", "sample:2", and should
      be coherent over all functions. Similar, their contents, its like
      "x:1", "x:2", etc.
*  Fixed a bug setting prior for the log baseline hazard in `inla.coxph`

# INLA 20.04.14
*  Small fix so that `result$mode$x` is written out in the
      case where `nhyper=0` and `num_threads>1`
*  Minor internal changes.

# INLA 20.04.06
*  Added link `loga`. Not yet documented.
*  First try on a new feature to more easily approximate
      the joint marginal for a subset of the latent field. This is a new
      option `selection` and corresponding `inla.rjmarginal()`
      to sample from it. 

*  Added check that `model="linear"` is not used with
      replicate or group, which is not intention.

# INLA 20.03.29
*  `MCMC` mode is now disabled

*  Skewness correction is now back as default, in 
      `inla.posterior.sample()`

*  Added family `xbinomial` that allow non-integer
      response.
*  Likelihood model `bgev` add (not yet complete), and was
      renamed from the experimental likelihood model `gev2`.

*  If `inla.call="remote"` is set,
      then `INLA:::inla.call.builtin()` is used
      if `inla.qinv()` and/or `inla.qsolve()` are
      used while constructing the model.

# INLA 20.03.17
*  Updated file `jointdataCD4.rds` in `exampledata/`

# INLA 20.03.09
*  Fixed a bug in the PIT calculations for the zeroinflated,
      type 0, of poisson, binomial and nbinomial.

# INLA 20.03.08
*  Added option `b.strategy` in `control.inla` to
      control what to do with the linear term when the `cmin` option is
      in effect
*  Added in-interval observed event in `inla.surv`

*  Added `dplyr` as suggested package as
      `dplyr::bind_rows` can replace
      `INLA::inla.rbind.data.frames`

# INLA 20.02.19
*  Added argument `E`, or `log(offset)`, to
      likelihood `gammacount`, so its equal to family `poisson`
      for `alpha=1`.

*  Minor changes

# INLA 20.01.25
*  Added a check that discrete observations are indeed
      integers, like for Poisson, Binomial, etc

*  The function `inla.binary.install` is now exported.
*  Added new likelihood family, `xpoisson`, which allows
      continous response: see the documentation for details (and note
      the error-check now done for discrete observations)

*  Added new likelihood `dgp` (discrete generalized Pareto)

*  Code clean-up (`contpoisson` and `qcontpoisson`)
*  Made `inla.pardiso.check()` a bit more informative if
      there is an error.

# INLA 19.12.10
*  Improved documentation of `inla.posterior.sample` and
      `inla.coxph`
*  Fixed an issue with `NA` data in the family `gev2`

# INLA 19.12.03
*  Updated some documentation about the `pc.gevtail` prior.
*  Reverted `inla.posterior.sample` back to the old
      version, the new experimental version is available as
      `INLA:::inla.posterior.sample.new`
*  Error in `Epil` data-set, `y[31]` should be
      23 not 21. 

# INLA 19.11.17
*  Updated the vignette about the multinomial distribution

*  New experimental windows binary built with
  `x86_64-w64-mingw32-gcc`, version 7.3, and linked with the
  pardiso library. Its stored in `bin/windows/experimental`

*  Updated `inla.qreordering` and updated `leuk-demo.R`
      example file (and the corresponding zip-file).

# INLA 19.11.10
*  Cache values of `qgamma` to speedup
      Gamma quantile regression

# INLA 19.10.30
*  Added a scaling constant for the precision parameter in the
      `qkumar` likelihood (to avoid instabilities). See updated
      documentation for details.

*  `inla.posterior.sample` now correct for possible skewness
      by default: see `?inla.posterior.sample` for details.

# INLA 19.10.16
*  Likelihoodmodel `betabinomialna`

# INLA 19.10.15
*  Default prior for the tail parameter in likelihood model
      `gp`, have changed to `pc.gevtail`, and the name change
      from `shape` to `tail`. It is now required to define a
      interval for the tail parameter, similar to `pc.gevtail`.

# INLA 19.10.06
*  Code-improvement for the `loggamma`-function
*  `barrier.R` updated (minor fix and code edits)

# INLA 19.10.02
*  Disable some debug output

# INLA 19.10.01
*  Fixed a bug in the `nmixnb` likelihood. 
*  Preserve names in `inla.posterior.sample.eval`
      if present.

# INLA 19.09.18
*  More work on the skew-normal link model

# INLA 19.09.15
*  `INLA:::inla.binary.install()` is a new interactive tool
      to install alternative Linux builds. 

# INLA 19.09.10
*  Added skew-normal link-model `sn` for binary data,
            with its PC-prior

# INLA 19.09.03
*  Added `robit` link model.

*  Improved the stability of the saturated deviance
      calculations
*  Fixed `INLA:::inla.is.list.of.lists` to cover the
      case where the arguments are a list of named lists

# INLA 19.07.27
*  New (experimental) likelihood: gev2

*  Fixed, again, an issue with (parallel) PARDISO
      and many linear combinations.
*  Minor code changes in `doc.R`

# INLA 19.07.21
*  Removed must-be-enabled warnings in some
      surival models, from Oct  25 2017

*  Added PC-prior for the Weibull likelihood models. The prior
      is derived
      for `variant = 1`, which is the good parameterisation.

*  Added missing `to.theta` and `from.theta`
      functions in likelihoods `sn` and `sn2`
*  Fix some documentation in `marginal.R` (refering to the
      obsolete function `inla.marginal.transform`)
*  Fixed an issue with (parallel) PARDISO and many linear combinations.

# INLA 19.05.19
*  Set `StagedInstall:no` to work around
      installation problems for MacOS and R-3.6

# INLA 19.05.17
*  The internal parameterisation of the alpha-parameter for the
  Weibull likelihood familes, has been redefined/scaled, to fix some
  optimisation issues. This means that the default prior has changed (a
  little) and user-defined priors has to change to account for this new
  internal parameterisation (sorry about that). See the documentation
  for details.

# INLA 19.05.16
*  Option `short.summary` will use a version of
      `summary` with less output, maybe more suitable for
      Markdown documents. 

# INLA 19.05.13
*  Added exampledata directory for various example datasets

*  Code cleanup and improved some input-error checking.

# INLA 19.04.16
*  Fixed an error in the cache-system for
      `model="rgeneric"` and `model="dmatern"`.
      Most notably with option `openmp.strategy="pardiso.parallel"`.

# INLA 19.04.14
*  Removed the weight correction for the computation of the cpo
            for `int.design="user.expert"`

# INLA 19.04.09
*  Option `int.strategy="user.expert"`, see the vignette
      about user-defined integration points.
*  Merge also `cpo` and `po` results in `inla.merge()`

# INLA 19.04.01
*  Fixed an issue with AR-model and group

# INLA 19.03.16
*  Small fix in model `dmatern`

# INLA 19.03.04
*  Redirect error output of some warning messages in the remote-feture
      section from MacOSX to Linux.
*  Faster return when `mu` is zero for `rgeneric`

# INLA 19.03.02
*  Changed from PARDISO.PARALLEL to PARDISO.SERIAL in `inla.qsample`
*  Optimize the `nhrs` for `inla.qsolve` for PARDISO

# INLA 19.02.28
*  Several fgn-models are now fine
*  Fixed CPU timing with the PARDISO library

# INLA 19.02.26
*  Do not need to optimize reordering when PARDISO is used. 

# INLA 19.02.17
*  Fixed input-test using `inla.qsample` with
      `selection`-argument.
*  Added back `family = "normal"` which is now
      translated to `family = "gaussian"` internally. 

# INLA 19.02.14
*  More work and fixes in  `inla.merge`

# INLA 19.02.12
*  Simplied `print.inla` output

*  New method `merge` and function `inla.merge`,
      for merging `inla`-objects

*  Store `control.family` after processing, in the
      `result$.args` argument, not just the calling value.

# INLA 19.02.09
*  New parameter for Gaussian likelihood: Fixed offset in the
      variance. 

*  Updated `envir` definition in the `rgeneric`
      documentation and examples. 

# INLA 19.02.06
*  Removed testing code for likelihood model `testbinomial1`

*  Added new likelihood `gamma.surv`

*  Cleaned up the use of temporary dir and files
*  General code clean-up

# INLA 19.01.29
*  Increased maximum number of covariates in likelihood models
      `nmix` and `nmixnb` from 10 to 15

# INLA 19.01.24
*  Added a new test-script

# INLA 18.12.12
*  New models, `loggamma` and `mloggamma`
      in `mix`.

*  Minor changes in some build scripts.

# INLA 18.12.01
*  New option `mkl` in `inla.setOption()` to chose
      MKL-buildt binaries.
*  Linux binaries now buildt with Ubuntu1804.
*  MKL-versions are included for MacOSX, and Linux (both dynamic
      and static).

# INLA 18.11.28
*  New latent model `intslope`

# INLA 18.11.22
*  Improved `control.mix` interface and code

# INLA 18.10.29
*  Likelihood model `nbinomial2`

# INLA 18.10.28
*  New function `inla.priors.used`

# INLA 18.10.17
*  Export class `inla`

# INLA 18.10.16
*  New latent model: `dmatern`

*  Improved the numerics for computing the scaling of the RW1 and RW2 models.

# INLA 18.10.09
*  New option `control.inla=list(tolerance.step=)`, to
      control the RMS of the step-size for the inner optimization.
*  Changed, slightly, the initial values for the exponent in
      the Weibull likelihood models, to a value close to zero instead of
      zero.
*  New vignette about how to deal with multinomial data.

*  Added option `verbose` to
      `inla.qsample()` and `inla.posterior.sample()`

# INLA 18.09.24
*  Performance improvement using the PARDISO library

# INLA 18.09.21
*  Argument `selection` in `inla.posterior.sample`
      and `inla.qsample`.

# INLA 18.09.19
*  Fix for `num.threads` in `inla.qinv()`

# INLA 18.09.18
*  Allow better user control of sparse matrix library in
      `inla.qinv()`, `inla.qsample()` and
      `inla.posterior.sample()`

# INLA 18.09.14
*  New example added to `inla.posterior.sample()`
*  Slight changes in the default `print`, and
      `summary` for an `inla`-object

*  Fixed the issue when `lincomb.derived.only=FALSE` and
      then using `inla.posterior.sample()`

# INLA 18.08.26
*  Added 32-bit builds for windows (upon request)

*  Added function `inla.posterior.sample.eval()`

# INLA 18.08.09
*  Added new function `inla.pardiso.check()`
*  Added COPYRIGHTS file

*  Separated the quantile link for the binomial response,
      into individual (`model="quantile"`) and population
	(`model="pquantile"`)

*  Added new strategy
      `control.inla=list(strategy="adaptive")` which use the
      `simplified.laplace` approximations for fixed effects and
      low-dimensional model components, and the `gaussian`
      approximation
      otherwise. The argument 
      `adaptive.max` in `control.inla`
      determines what is low-dimensional in this
      context (default 10).

*  Removed some code not used anymore

# INLA 18.07.27
*  NEWS page created (see `news(package="INLA")`)
*  Added vignette about the conditional logit model (thanks to
      Stefani Muff)
*  Fixed missprints in the documentation for model `ar1c`
      (Thanks to Virgilio Gomez Rubio)
*  Fixed documentation about argument `blas.num.threads` in `inla()`

# INLA 18.07.12
*  Package built with `R-3.5`, both stable and testing

# INLA 18.07.11
*  Package built for `R-3.4`, both stable and testing. 

# INLA 18.07.12
*  Likelihood model `pom` (proportional odds model)

