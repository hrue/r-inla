# INLA 23.02.27
* New option to \code{inla.binary.install()}
* Fix a broken build

# INLA 23.02.22
* Add linear predictor summaries to config
* Vignette about \code{rgeneric} is updated
* Updated doc about \code{inla.posterior.sample}
* Fixed issue with missing covariance in A-matrix

# INLA 23.02.04
* Fix for \code{control.mode=list(..., fixed=TRUE)}
* Various code-cleanup and improvement

# INLA 23.01.31
* Remove `rgeos` dependency
* Remove `rgdal` dependency
* Add new link-function: \code{ccloglog}
* New argument to \code{inla.group.cv}
* \code{safe}-mode improved in the \code{inla()}-call
* Various code-cleanups

# INLA 23.01.12
* Fixed an issue with ptweedie
* Code cleanup

# INLA 23.01.02
* Sync C++ code from the fmesher package
* Added log-likelihood output to configs
* Fixed a regression for \code{model="fgn"}
* Improved code for \code{model="binomial"} and logit-link
* Increased kmax from 10 to 20 for \code{model="iidkd"}

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
* Fixed an issue with cache and nested parallelism (mostly \code{cgeneric}-related)

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
* Added default values of 'safe', 'verbose' and 'keep' to options
* Added new check for valid arguments (default off for the moment)
* 'control.compute=list(residuals=TRUE)' provide expected deviance
  residuals, but it is experimental for the moment

# INLA 22.10.23
* Optimisation of gcpo and inner-products
* Port for upcoming Matrix 1.5-2 change in 'kronecker()'
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
* The 'weibullcure'-model is removed as its covered by the new
  cure-model.
* Optimization: Models with many constraints 
* Optimization: Models with 'betabinomial' likelihood

# INLA 22.10.07
* Fixed regression in initial value for model="copy"
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
## BUG FIXES
* Fixed an issue with the deprecated class coercions in `inla.as.dgTMatrix()`


# INLA 22.08.23
## BUG FIXES
* Avoid deprecated class coercions in `inla.as.dgTMatrix()`
* Fixed an issue with `rgeneric` 
## CHANGES
* New vignette added
* Various code improvements


# INLA 22.07.23
## CHANGES
*  Fixed regression from `Version_22.07.15`



# INLA 22.07.21
## CHANGES
*  Code-improvement and minor fixes
*  Fixed regression from `Version_22.07.15`



# INLA 22.07.15
## CHANGES
*  Code-improvement and minor fixes



# INLA 22.07.01
## CHANGES
*  Code-improvement/fix



# INLA 22.06.30
## CHANGES
*  Code-improvement
*  Changed stopping criteria for VB corrections



# INLA 22.06.25
## CHANGES
*  Code-improvement
*  More work on the new experimental feature (internal)



# INLA 22.06.20
## CHANGES
*  Code-improvement
*  Adding new experimental feature (internal)



# INLA 22.06.11
## CHANGES
*  Added more to `gcpo`-output
*  Code-improvement



# INLA 22.06.03
## CHANGES
*  Bug-fix release



# INLA 22.06.02
## CHANGES
*  The `mode` in the output summary is back.
*  Family `gev` is disabled, using `bgev` instead
*  Code improvement, cleanup and bug-fixes



# INLA 22.05.18
## CHANGES
*  Changes in work-in-progress features



# INLA 22.05.07
## CHANGES
*  Minor changes to move to R-4.2



# INLA 22.05.03
## CHANGES
*  Code improvement, cleanup and bug-fixes
*  More work on the `gcpo`-feature
*  Removed threadprivate from code (almost), which was a huge rewrite.
*  Made linesearch in the optimization more robust and stable
*  Improved `stupid`-search
*  Internal changes in the `smart`-optimization
*  (This is likely the last build with R-4.1 before we switch to R-4.2)



# INLA 22.04.16
## CHANGES
*  Added `theta`-correction for `gcpo` computations
*  Code improvement and cleanup



# INLA 22.04.14
## CHANGES
*  Improve `gcpo`-handling of singular cases. `group.size=-1` now
    gives the CPO.



# INLA 22.04.13
## CHANGES
*  Bug-fix release



# INLA 22.04.09
## CHANGES
*  Bug-fix and code improvement release



# INLA 22.04.06
## CHANGES
*  Bug-fix and code improvement release



# INLA 22.03.28
## CHANGES
*  Bug-fix release



# INLA 22.03.27
## CHANGES
*  Bug-fix release



# INLA 22.03.26
## CHANGES
*  Bug-fix release



# INLA 22.03.25
## CHANGES
*  Bug-fix release



# INLA 22.03.24
## CHANGES

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
## CHANGES

*  Code cleanup and a bug-fix
*  Dense matrices are now store column-wise in `cgeneric`
*  New example added to `cgeneric`, see vignette
*  Changed default integration stratey for `experimental` mode




# INLA 22.03.12
## NEW EXPERIMENTAL FEATURES

*  Option `safe=TRUE` is now default


## CHANGES

*  Code cleanup
*  Code optimization for `TACUS`
*  Nested parallelism is now enabled on Windows




# INLA 22.03.09
## BUG FIXES

*  Issue with `gcpo` and `safe=T` fixed
*  Fixed issue with `A`-matrix in `gcpo`




# INLA 22.03.08
## CHANGES

*  Better TAUCS performance: improved storage scheme and
      new code tp allow for linear solves with many rhs.
*  `gcpo`-improvement. The default is now `group.size=1`,
      and it acts like cpo.


## BUG FIXES

*  Fixed an issue for propagation an error
      from the cpo-calculations
*  Various smaller changes




# INLA 22.03.06
## CHANGES

*  Improved the TAUCS-interface performance.
*  Various code improvements.




# INLA 22.03.03
## NEW EXPERIMENTAL FEATURES

*  New options for the group-cpo feature.


## BUG FIXES

*  Fixed CRAN-check issues.




# INLA 22.02.28
## NEW EXPERIMENTAL FEATURES

*  VB correction for experimental mode do now include nodes
      from all model components, not only the short ones and fixed
      effects.
*  New group-cpo feature, very experimental at the moment
      (internal use only).


## VARIOUS CHANGES

*  Parts of the VB code has been improved for
      simplicity and speed.
*  `inla.call="remote"` now use `zstd` for
      parallel compression, meaning that `zstd` must be
      available on both sides.
*  PARDISO license help and startup msg has been disabled.
*  Minor fixes, code improvements and cleanup




# INLA 22.02.16
## NEW EXPERIMENTAL FEATURES

*  New feature `inla.cgeneric.q`
*  New and updated PARDISO library


## BUG FIXES

*  Code improvements and cleanup




# INLA 22.01.25
## BUG FIXES

*  Fixed regression bug in `knmodels`
*  Code improvements and cleanup




# INLA 22.01.19
## BUG FIXES

*  Fixed bug in the variable expansion in some arguments
      that could give different results than before




# INLA 22.01.16
## NEW EXPERIMENTAL FEATURES

*  Added new variant for nbinomial likelihood




# INLA 22.01.12
## IMPROVEMENTS

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


## BUG FIXES

*  Fixed regression error with `plot(result, plot.prior=TRUE)`
*  Fixed a regression error from version 21.12.21 which expanded some
      variables incorrect.




# INLA 21.12.21
## NEW EXPERIMENTAL FEATURES

*  Option `safe` in function `inla`


## BUG FIXES

*  Minor bug fixes




# INLA 21.12.18
## BUG FIXES

*  Improved handling of `NAN` and `INF`




# INLA 21.12.16
## BUG FIXES

*  A bug-fix for `matern2d` model
*  Various code improvements




# INLA 21.12.11
## BUG FIXES

*  Minor bug fixes for `cgeneric` for Windows




# INLA 21.12.10
## NEW EXPERIMENTAL FEATURES

*  A `cgeneric` interface
      like `rgeneric` but with C-code.
      This is for internal testing only.


## BUG FIXES

*  Various code improvements and minor fixes





# INLA 21.12.01
## NEW FEATURES

*  Allow for probit cdf in likelihood model `pom`


## BUG FIXES

*  Various code improvements




# INLA 21.11.29
## NEW EXPERIMENTAL FEATURES

*  Small fixes




# INLA 21.11.28
## NEW EXPERIMENTAL FEATURES

*  added new link `powerprobit`. Not completely done yet.




# INLA 21.11.22
## BUG FIXES AND IMPROVEMENTS

*  Remove jemalloc-option that creates warning message for Mac




# INLA 21.11.21
## BUG FIXES AND IMPROVEMENTS

*  Minor internal changes
*  Add jemalloc options for Linux/Mac




# INLA 21.11.16
## BUG FIXES AND IMPROVEMENTS

*  Minor internal changes for robustness
*  Upgraded compilers to gcc-11 on Windows




# INLA 21.11.15
## NEW EXPERIMENTAL FEATURES

*  Option for parallel linesearch has been renamed into
      `parallel.linesearch` in `control.inla`.


## BUG FIXES AND IMPROVEMENTS

*  This is fix for a build-error in the Windows binary for
      version `21.11.14`
*  Various minor internal code changes and improvements




# INLA 21.11.14
## NEW EXPERIMENTAL FEATURES

*  Parallel linesearch implemented, see option
       `bfgs.version` in `control.inla`. This option is highly
       experimental at the moment and work in progress, and should
       not be used.


## BUG FIXES AND IMPROVEMENTS

*  Change default settings for `family="tweedie"`
*  Implemented an adaptive parallel/serial version of `Qx`
*  Improved perforance for parallel linear solve with `PARDISO`
*  Use stable Legendre polynomial evaluation for spherical covariances
*  Add support for `Fedora Linux 35`
*  Fixed `dic.sat` when all response are `NA`
*  Various internal code changes and improvements





# INLA 21.11.01
## NEW EXPERIMENTAL FEATURES

*  New experimental feature `lp.scale`. No good documentation yet, but
      that will come soon.


## BUG FIXES

*  Robustify `inla.iidkd.sample` function for numerical singular matrices
*  Corrected `inla.stack.join` handling of factor variables that in some
      cases converted variables to/from factors unintentionally.
*  `inla.mode="experimental"` got its fitted values back.




# INLA 21.10.03
## BUG FIXES

*  Internal code improvement




# INLA 21.09.29-1
## BUG FIXES

*  Better special number


## NEW EXPERIMENTAL FEATURES

*  New option to `inla.iidkd.sample`




# INLA 21.09.29
## NEW EXPERIMENTAL FEATURES

*  New latent model `iidkd` which use a different
      parameterisation than `iid3d` and similar ones.




# INLA 21.09.13
## BUG FIXES

*  Fixed crash when using `int.strategy="user"` with
      `inla.mode="experimental"`
*  Minor fixes




# INLA 21.09.11
## BUG FIXES

*  Fixed crash when using linear combinations with
      `inla.mode="experimental"`
*  Minor fixes




# INLA 21.09.01
## NEW EXPERIMENTAL FEATURES

*  Added option `scale` to `family="xbinomial"`




# INLA 21.08.31
## BUG FIXES

*  Various bug-fixes





# INLA 21.07.10
## USER-VISIBLE CHANGES

      
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


## NEW FEATURES

*  New option `inla.timeout`, similar to
      `fmesher.timeout`; see `?inla.setOption`
      
*  `CPO-marginals` are returned in when
      `control.compute=list(config=TRUE)`
      
*  There is a change how a (numerical) singular constraint is
      treated, hopefully this new way is more robust.


  
## NEW EXPERIMENTAL FEATURES

*  New option `inla.mode`, which define how to arrange the
      internal model formulation. One of `"classic"`,
      `"twostage"` and `"experimental"`. The default is
      `"classic"`, which is unchanged behaviour compared for
      earlier versions. The other two are highly experimental for the
      moment. See `?inla`


  
## BUG FIXES

*  Massive code clean-up and some minor fixes.





# INLA 21.06.11
## NEW EXPERIMENTAL FEATURES

*  Added support for Mac M1 (native build with R-4.1)
*  New `twostage` options (highly experimental)


## BUG FIXES

*  Massive code overhaul
*  Various minor fixes




# INLA 21.05.02
## BUG FIXES

*  Some minor bug-fixes and code improvements




# INLA 21.04.22
## NEW EXPERIMENTAL FEATURES

*  Mac only: Added optional path to `vecLib`
      BLAS and LAPACK libraries.
*  Adding argument `.special` to `inla.surv`-object.


## BUG FIXES

*  Some code improvement 





# INLA 21.04.16
## NEW EXPERIMENTAL FEATURES

*  Mac only: Link with `vecLib` BLAS and LAPACK libraries
      by default. Turn off with `inla.setOption(vecLib=FALSE)`.
*  Adding an experimental improvement
      for `strategy=simplified.laplace` (not enabled by default)


## BUG FIXES

*  Some code improvement 




# INLA 21.03.31
## BUG FIXES

*  General code improvement 




# INLA 21.03.27
## NEW EXPERIMENTAL FEATURES

*  New family `stochvolsn`
*  New family `cenpoisson2`


## BUG FIXES

*  Improve robustness for survival models
*  General code improvement




# INLA 21.03.21
## NEW EXPERIMENTAL FEATURES

*  Changes in `inla.jp`; see help-page for details.




# INLA 21.03.20
## NEW EXPERIMENTAL FEATURES

*  New family `gompertz`
*  `inla.binary.install()` now do md5-checksum check


## BUG FIXES

*  Code optimization for `barrier` models





# INLA 21.03.17
## BUG FIXES

*  Internal code changes for `family=tweedie`




# INLA 21.03.16
## NEW EXPERIMENTAL FEATURES

*  `inla.binary.install()` can now also run
            non-interactive


## BUG FIXES

*  Internal code changes




# INLA 21.03.14
## BUG FIXES

*  Internal code changes




# INLA 21.03.13
## BUG FIXES

*  Fixed some missprints in doc
*  Trying to get vigettes to work again




# INLA 21.03.08
## NEW EXPERIMENTAL FEATURES

*  Added family `logperiodogram` back in
*  Added new family `agaussian`


## BUG FIXES

*  Improved initial values for `x` 
*  Improved computations for saturated likelihood
*  Various minor code cleanup




# INLA 21.02.23
## NEW FEATURES

*  Default `strategy` is now `simplified.laplace`
      for smaller models (< 5000 nodes) and `adaptive` for larger ones


## BUG FIXES

*  Minor bug-fixes




# INLA 21.01.26
## BUG FIXES

*  Improve parallel rerodering for PARDISO
*  Improve `coxph` for large data
*  Improve internal CPU-timing output




# INLA 99.99.99
## USER-VISIBLE CHANGES

*  


## NEW FEATURES

*  


## NEW EXPERIMENTAL FEATURES

* 


## BUG FIXES

*  




# INLA 21.01.18
## BUG FIXES

*  Use METIS5 with PARDISO




# INLA 21.01.13
## NEW EXPERIMENTAL FEATURES

*  Family `fmri`


## BUG FIXES

*  Speed improvements for family `tweedie`




# INLA 21.01.08
## NEW EXPERIMENTAL FEATURES

*  New family `tweedie`


## BUG FIXES

*  Optimization work, code improvement and maintenance




# INLA 20.12.10
## NEW EXPERIMENTAL FEATURES

*  Added a new experimental feature for running PARDISO in parallel


## BUG FIXES

*  Various minor fixes and improvements




# INLA 20.12.04
## BUG FIXES

*  Bug-fix release only




# INLA 20.11.29
## USER-VISIBLE CHANGES

*  Added new section to the vignette about `rgeneric`


## NEW FEATURES

*  Option `optimize` in `inla.rgeneric.define`


## BUG FIXES

*  Cleanup in some output format.
*  Various internal changes. 




# INLA 20.11.22
## BUG FIXES

*  Code improvement and optimization for `rgeneric`.
*  Some OpenMP improvements




# INLA 20.11.18
## NEW FEATURES

*  New family `gammajw`


## BUG FIXES

*  Fixed an OpenMP deadlock case, using `rgeneric`




# INLA 20.11.16
## USER-VISIBLE CHANGES

*  PARDISO version 7 is included for Mac and Linux.
*  Improved paralellism and nested paralellism
*  Speedup improvements, especially for models with many
      constraints.
*  Improved default `plot(result)`. Argument
      `cex=..` will now work


## NEW FEATURES

*  New shorthand feature in `inla.posterior.sample.eval`
      for extracing model components from samples
*  Added `fmesher.timeout` option, see `?inla.setOption`


## BUG FIXES

*  Using a non-zero seed in `inla.posterior.sample` will
      now force serial computations.
*  A lot of internal code cleanup and improvements




# INLA 20.10.11
## NEW EXPERIMENTAL FEATURES

*  New likelihood `poisson.special1`


## BUG FIXES

*  Minor code cleanup
*  Minor changes in code/doc to adapt to R-4.0




# INLA 20.09.25
## BUG FIXES

*  Bug fixes and some optimisation work




# INLA 20.09.16
## NEW EXPERIMENTAL FEATURES

*  Likelihood model `sn` is redone.


## BUG FIXES

*  Various minor changes.




# INLA 20.09.13
## NEW FEATURES

*  Model `iid` added to possible models
      for the baseline hazard in the Cox-ph model


## NEW EXPERIMENTAL FEATURES

*  Likelihood models `sn` and `sn2` have been
      replaced with a new version `sn`. The parameterisation is
      different and now done correct.




# INLA 20.09.09
## USER-VISIBLE CHANGES

*  The likelihood models `sn` and `sn2`
      are now disabled. They need a rewrite to be done right.


## NEW EXPERIMENTAL FEATURES

*  For Linux, then build so-libraries are loaded before the
      system ones. This behaviour can be reverted by setting the
      environment variable `INLA_NATIVE_LD_LIBRARY_PATH`.




# INLA 20.09.07
## USER-VISIBLE CHANGES

*  Dimension for CCD integration have been
      increased from 38 to 52.


## NEW EXPERIMENTAL FEATURES

*  Default prior for the `intercept` and the `skewness` in the
      `skew-normal link` have changed. 


## BUG FIXES

*  Various minor changes




# INLA 20.08.31
## BUG FIXES

*  A small bug-fix with argument `-t`
*  Some improvements in the install scripts




# INLA 20.08.30
## BUG FIXES

*  More work on the nested parallelism stuff




# INLA 20.08.29
## BUG FIXES

*  Argument `num.threads="A"` now means
      `num.threads="A:1"`
*  Some changes for nested paralellism




# INLA 20.08.28
## USER-VISIBLE CHANGES

*  There is a minor change in `inla.qsample` and `inla.posterior.sample`.
      If argument `seed!=0` then serial mode is now forced.
      Earlier, an error would be raised if parallel mode was also requested.


## NEW EXPERIMENTAL FEATURES

*  Argument `remove.names` in `control.fixed`


## BUG FIXES

*  Bug-fixes for nested parallism
*  Some optimization improvement




# INLA 20.08.24
## BUG FIXES

*  Some minor fixes in the install


## MISC

*  Linux builds now link with PARDISO version 7 (beta)




# INLA 20.08.22
## NEW EXPERIMENTAL FEATURES

*  MKL is now default for both Mac and Linux.
      If this cause any issues, revert back by setting
      `inla.setOption(mkl=FALSE)`.
*  As of today, the TAUCS library does not play 
      well with MKL on Mac for some reason, and this is
      silent accounted for.
      It is possible to bypass the internal check,
      and if that becomes in issue, just define
      `mkl=FALSE` as above.


## BUG FIXES

*  Updated some install scripts




# INLA 20.08.21
## BUG FIXES

*  Various minor fixes
*  MKL binaries are now linked with libiomp5




# INLA 20.08.19
## BUG FIXES

*  Changed to the gcc-10 compiler suite for Ubuntu 1804
      and 2004.
*  Fixed an issue if the optimiser does not move
      and the use of directions
*  Minor fixes




# INLA 20.08.18
## USER-VISIBLE CHANGES

*  Code is prepared for nested parallism, which most gain
      obtained when using the PARDISO library
*  Argument `num.threads` is now in the format `A:B`,
      where A threads for the outer layer and B threads for the inner
      layer. And this applies for several functions.
*  The speed have improved, mostly when using the PARDISO
      library. Hopefully, nothing is broken.


## NEW FEATURES

*  likelihood `beta` now allow for cencoring near 0 and 1,
      see `inla.doc("beta", section="likelihood")`


## NEW EXPERIMENTAL FEATURES

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
## NEW EXPERIMENTAL FEATURES

*  Add new experimental feature
      `control.inla=list(use.directions=TRUE/FALSE)`




# INLA 20.08.09
## NEW EXPERIMENTAL FEATURES

*  New experimental optimise strategy
      `control.inla=list(optmise.strategy="smart"`
      that is hopefully faster and as safe as the default one. Maybe
      even more safe and robust.


## BUG FIXES

*  More work on the nested paralellism 




# INLA 20.08.06
## NEW EXPERIMENTAL FEATURES

*  More work on nested parallelism and some code changes
*  Code changes to reduce memory usage
*  Linux/Mac binaries are now linked with `jemalloc`




# INLA 20.08.04
## NEW EXPERIMENTAL FEATURES

*  Improved nested parallelism and some code changes




# INLA 20.08.03
## NEW EXPERIMENTAL FEATURES

*  More work on nested parallelism




# INLA 20.08.02
## USER-VISIBLE CHANGES

*  Option `control.inla$lincomb.derived.only` is now disabled.


## NEW EXPERIMENTAL FEATURES

*  Testing nested parallelism `openmp.nested` with
      `num.threads="A,B"`. Work in progress




# INLA 20.07.27
## NEW FEATURES

*  `plot(result)` will now produce a plot of the
      `CPO/PIT` for each likelihood family (if available), instead
      of a joint plot as earlier.
*  New vignette about `jmarginal`


## NEW EXPERIMENTAL FEATURES

*  Added new likelihood models `zeroinflatedcenpoisson0`
      and `zeroinflatedcenpoisson1`
*  Link-model `sn` is updated, as well as the
      PC-prior for the skewness therein, and the added intercept model.


## BUG FIXES

*  Revised the PIT calculations for family `cenpoisson`
*  Code rewrite to (try to) prevent `Inf` for `DIC` calculations
*  Minor fixed in `inla.binary.install`




# INLA 20.07.18
## BUG FIXES

*  Build-script changes and misc fixes




# INLA 20.07.16
## USER-VISIBLE CHANGES

*  Packages `mpoly` and `symmoments` are
      added to the Suggests-list.


## NEW FEATURES

*  Added info about `inla.prune()` to the
      startup message


## NEW EXPERIMENTAL FEATURES

*  More features for `jmarginal` added


## BUG FIXES

*  Build-scripts fixes
*  Fix for a rare `fmesher` issue
*  Improved the code for the DIC calculations to make
      them more stable
*  Some improvment in the PC-prior for the SN-link




# INLA 20.07.12
## BUG FIXES

*  Fixed `inla.link.sn` to vectorise over argument
      `a` as `sn`-package do not do that properly itself. 




# INLA 20.07.09
## NEW FEATURES

*  Package built with R-4.0




# INLA 20.07.04
## USER-VISIBLE CHANGES

*  Vignette added for `family=bGEV`
*  Some internal changes due to the migration to `git`



# INLA 20.06.29
## NEW EXPERIMENTAL FEATURES

*  Added `Qprior.diag` to the output when `config=TRUE`.
      The off-diagonals of this matrix are the same as `Q` in the
      same configuration, so only the diagonal of `Qprior`
      is stored.
*  Added some internal experimental code


## BUG FIXES

*  PARDISO interface: internal check added
*  Fixed an wrong assert with family=bgev




# INLA 20.06.22
## BUG FIXES

*  Improve some code in the PARDISO interface
*  Improved the computation of the third derivative
      in the log likelihood.






# INLA 20.06.18
## BUG FIXES

*  Improved WKT support for PROJ6/GDAL3




# INLA 20.05.16
## NEW FEATURES

*  Support for PROJ6/RGDAL3 for handling CRS information for
        spatial objects.




# INLA 20.06.15
## NEW EXPERIMENTAL FEATURES

*  For the `intslope`-model: made all `gamma`'s
      default fixed to 1, so its similar in style the copy-feature.
*  Added argument `constr` it `inla.rjmarginal`
*  Added argument `ask` to function `inla.prune`


## BUG FIXES

*  Argument `cyclic=TRUE` in `f()` should not set
            `constr=FALSE` when default is `constr=TRUE`
*  Change the `scale.model=TRUE` code for `RW1/RW2` so the
      scaling for the continous case is the same as for the discrete
      case when the locations are eqvidistant.
*  Disable link `sslogit`




# INLA 20.05.12
## BUG FIXES

*  Fixed an issue with model `besag2`
*  Fixed an issue with `plot(r,plot.prior=TRUE)` for some priors




# INLA 20.05.04
## NEW FEATURES

*  Remove the experimental status of `inla.posterior.sample.eval`


## NEW EXPERIMENTAL FEATURES

*  Added function `inla.prune` which will remove binaries
      not supported by the running OS, to reduce the size of the
      package.
*  Added method `summary` and `print` to class
      `inla.jmarginal`


## BUG FIXES

*  Add check for `NA/NaN/Inf` in mesh creation input
      ocations
*  Make sure that skewness is not to high in `inla.posterior.sample`




# INLA 20.04.18
## NEW FEATURES

*  Added new argument `tag` to `inla.coxph`


## NEW EXPERIMENTAL FEATURES

*  `inla.rjmarginal.eval`, to evaluate samples from a join
      approximations


## BUG FIXES

*  Names of samples are now "sample:1", "sample:2", and should
      be coherent over all functions. Similar, their contents, its like
      "x:1", "x:2", etc.
*  Fixed a bug setting prior for the log baseline hazard in `inla.coxph`




# INLA 20.04.14
## BUG FIXES

*  Small fix so that `result$mode$x` is written out in the
      case where `nhyper=0` and `num_threads>1`
*  Minor internal changes.




# INLA 20.04.06
## NEW EXPERIMENTAL FEATURES

*  Added link `loga`. Not yet documented.
*  First try on a new feature to more easily approximate
      the joint marginal for a subset of the latent field. This is a new
      option `selection` and corresponding `inla.rjmarginal()`
      to sample from it. 


## BUG FIXES

*  Added check that `model="linear"` is not used with
      replicate or group, which is not intention.




# INLA 20.03.29
## USER-VISIBLE CHANGES

*  `MCMC` mode is now disabled


## NEW FEATURES

*  Skewness correction is now back as default, in 
      `inla.posterior.sample()`


## NEW EXPERIMENTAL FEATURES

*  Added family `xbinomial` that allow non-integer
      response.
*  Likelihood model `bgev` add (not yet complete), and was
      renamed from the experimental likelihood model `gev2`.


## BUG FIXES

*  If `inla.call="remote"` is set,
      then `INLA:::inla.call.builtin()` is used
      if `inla.qinv()` and/or `inla.qsolve()` are
      used while constructing the model.




# INLA 20.03.17
## BUG FIXES

*  Updated file `jointdataCD4.rds` in `exampledata/`





# INLA 20.03.09
## BUG FIXES

*  Fixed a bug in the PIT calculations for the zeroinflated,
      type 0, of poisson, binomial and nbinomial.




# INLA 20.03.08
## NEW EXPERIMENTAL FEATURES

*  Added option `b.strategy` in `control.inla` to
      control what to do with the linear term when the `cmin` option is
      in effect
*  Added in-interval observed event in `inla.surv`


## BUG FIXES

*  Added `dplyr` as suggested package as
      `dplyr::bind_rows` can replace
      `INLA::inla.rbind.data.frames`




# INLA 20.02.19
## USER-VISIBLE CHANGES

*  Added argument `E`, or `log(offset)`, to
      likelihood `gammacount`, so its equal to family `poisson`
      for `alpha=1`.


## BUG FIXES

*  Minor changes




# INLA 20.01.25
## USER-VISIBLE CHANGES

*  Added a check that discrete observations are indeed
      integers, like for Poisson, Binomial, etc


## NEW FEATURES

*  The function `inla.binary.install` is now exported.
*  Added new likelihood family, `xpoisson`, which allows
      continous response: see the documentation for details (and note
      the error-check now done for discrete observations)


## NEW EXPERIMENTAL FEATURES

*  Added new likelihood `dgp` (discrete generalized Pareto)


## BUG FIXES

*  Code clean-up (`contpoisson` and `qcontpoisson`)
*  Made `inla.pardiso.check()` a bit more informative if
      there is an error.




# INLA 19.12.10
## BUG FIXES

*  Improved documentation of `inla.posterior.sample` and
      `inla.coxph`
*  Fixed an issue with `NA` data in the family `gev2`




# INLA 19.12.03
## BUG FIXES

*  Updated some documentation about the `pc.gevtail` prior.
*  Reverted `inla.posterior.sample` back to the old
      version, the new experimental version is available as
      `INLA:::inla.posterior.sample.new`
*  Error in `Epil` data-set, `y[31]` should be
      23 not 21. 




# INLA 19.11.17
## USER-VISIBLE CHANGES

*  Updated the vignette about the multinomial distribution


## NEW EXPERIMENTAL FEATURES

*  New experimental windows binary built with
  `x86_64-w64-mingw32-gcc`, version 7.3, and linked with the
  pardiso library. Its stored in `bin/windows/experimental`


## BUG FIXES

*  Updated `inla.qreordering` and updated `leuk-demo.R`
      example file (and the corresponding zip-file).




# INLA 19.11.10
## NEW EXPERIMENTAL FEATURES

*  Cache values of `qgamma` to speedup
      Gamma quantile regression




# INLA 19.10.30
## USER-VISIBLE CHANGES

*  Added a scaling constant for the precision parameter in the
      `qkumar` likelihood (to avoid instabilities). See updated
      documentation for details.


## NEW FEATURES

*  `inla.posterior.sample` now correct for possible skewness
      by default: see `?inla.posterior.sample` for details.





# INLA 19.10.16
## NEW EXPERIMENTAL FEATURES

*  Likelihoodmodel `betabinomialna`





# INLA 19.10.15
## USER-VISIBLE CHANGES

*  Default prior for the tail parameter in likelihood model
      `gp`, have changed to `pc.gevtail`, and the name change
      from `shape` to `tail`. It is now required to define a
      interval for the tail parameter, similar to `pc.gevtail`.




# INLA 19.10.06
## BUG FIXES

*  Code-improvement for the `loggamma`-function
*  `barrier.R` updated (minor fix and code edits)




# INLA 19.10.02
## BUG FIXES

*  Disable some debug output




# INLA 19.10.01
## BUG FIXES

*  Fixed a bug in the `nmixnb` likelihood. 
*  Preserve names in `inla.posterior.sample.eval`
      if present.




# INLA 19.09.18
## BUG FIXES

*  More work on the skew-normal link model




# INLA 19.09.15
## NEW EXPERIMENTAL FEATURES

*  `INLA:::inla.binary.install()` is a new interactive tool
      to install alternative Linux builds. 




# INLA 19.09.10
## NEW EXPERIMENTAL FEATURES

*  Added skew-normal link-model `sn` for binary data,
            with its PC-prior





# INLA 19.09.03
## NEW EXPERIMENTAL FEATURES

*  Added `robit` link model.


## BUG FIXES

*  Improved the stability of the saturated deviance
      calculations
*  Fixed `INLA:::inla.is.list.of.lists` to cover the
      case where the arguments are a list of named lists





# INLA 19.07.27
## NEW FEATURES

*  New (experimental) likelihood: gev2


## BUG FIXES

*  Fixed, again, an issue with (parallel) PARDISO
      and many linear combinations.
*  Minor code changes in `doc.R`




# INLA 19.07.21
## USER-VISIBLE CHANGES

*  Removed must-be-enabled warnings in some
      surival models, from Oct  25 2017


## NEW FEATURES

*  Added PC-prior for the Weibull likelihood models. The prior
      is derived
      for `variant = 1`, which is the good parameterisation.


## BUG FIXES

*  Added missing `to.theta` and `from.theta`
      functions in likelihoods `sn` and `sn2`
*  Fix some documentation in `marginal.R` (refering to the
      obsolete function `inla.marginal.transform`)
*  Fixed an issue with (parallel) PARDISO and many linear combinations.




# INLA 19.05.19
## BUG FIXES

*  Set `StagedInstall:no` to work around
      installation problems for MacOS and R-3.6




# INLA 19.05.17
## USER-VISIBLE CHANGES

*  The internal parameterisation of the alpha-parameter for the
  Weibull likelihood familes, has been redefined/scaled, to fix some
  optimisation issues. This means that the default prior has changed (a
  little) and user-defined priors has to change to account for this new
  internal parameterisation (sorry about that). See the documentation
  for details.




# INLA 19.05.16
## NEW FEATURES

*  Option `short.summary` will use a version of
      `summary` with less output, maybe more suitable for
      Markdown documents. 




# INLA 19.05.13
## NEW EXPERIMENTAL FEATURES

*  Added exampledata directory for various example datasets


## BUG FIXES

*  Code cleanup and improved some input-error checking.




# INLA 19.04.16
## BUG FIXES

*  Fixed an error in the cache-system for
      `model="rgeneric"` and `model="dmatern"`.
      Most notably with option `openmp.strategy="pardiso.parallel"`.




# INLA 19.04.14
## BUG FIXES

*  Removed the weight correction for the computation of the cpo
            for `int.design="user.expert"`




# INLA 19.04.09
## NEW EXPERIMENTAL FEATURES

*  Option `int.strategy="user.expert"`, see the vignette
      about user-defined integration points.
*  Merge also `cpo` and `po` results in `inla.merge()`




# INLA 19.04.01
## BUG FIXES

*  Fixed an issue with AR-model and group




# INLA 19.03.16
## BUG FIXES

*  Small fix in model `dmatern`





# INLA 19.03.04
## BUG FIXES

*  Redirect error output of some warning messages in the remote-feture
      section from MacOSX to Linux.
*  Faster return when `mu` is zero for `rgeneric`




# INLA 19.03.02
## BUG FIXES

*  Changed from PARDISO.PARALLEL to PARDISO.SERIAL in `inla.qsample`
*  Optimize the `nhrs` for `inla.qsolve` for PARDISO




# INLA 19.02.28
## BUG FIXES

*  Several fgn-models are now fine
*  Fixed CPU timing with the PARDISO library




# INLA 19.02.26
## BUG FIXES

*  Do not need to optimize reordering when PARDISO is used. 




# INLA 19.02.17
## BUG FIXES

*  Fixed input-test using `inla.qsample` with
      `selection`-argument.
*  Added back `family = "normal"` which is now
      translated to `family = "gaussian"` internally. 




# INLA 19.02.14
## BUG FIXES

*  More work and fixes in  `inla.merge`




# INLA 19.02.12
## USER-VISIBLE CHANGES

*  Simplied `print.inla` output


## NEW EXPERIMENTAL FEATURES

*  New method `merge` and function `inla.merge`,
      for merging `inla`-objects


## BUG FIXES

*  Store `control.family` after processing, in the
      `result$.args` argument, not just the calling value.





# INLA 19.02.09
## NEW FEATURES

*  New parameter for Gaussian likelihood: Fixed offset in the
      variance. 


## BUG FIXES

*  Updated `envir` definition in the `rgeneric`
      documentation and examples. 




# INLA 19.02.06
## USER-VISIBLE CHANGES

*  Removed testing code for likelihood model `testbinomial1`


## NEW EXPERIMENTAL FEATURES

*  Added new likelihood `gamma.surv`


## BUG FIXES

*  Cleaned up the use of temporary dir and files
*  General code clean-up




# INLA 19.01.29
## USER-VISIBLE CHANGES

*  Increased maximum number of covariates in likelihood models
      `nmix` and `nmixnb` from 10 to 15




# INLA 19.01.24
## BUG FIXES

*  Added a new test-script




# INLA 18.12.12
## NEW EXPERIMENTAL FEATURES

*  New models, `loggamma` and `mloggamma`
      in `mix`.


## BUG FIXES

*  Minor changes in some build scripts.




# INLA 18.12.01
## NEW EXPERIMENTAL FEATURES

*  New option `mkl` in `inla.setOption()` to chose
      MKL-buildt binaries.
*  Linux binaries now buildt with Ubuntu1804.
*  MKL-versions are included for MacOSX, and Linux (both dynamic
      and static).




# INLA 18.11.28
## NEW EXPERIMENTAL FEATURES

*  New latent model `intslope`




# INLA 18.11.22
## BUG FIXES

*  Improved `control.mix` interface and code




# INLA 18.10.29
## NEW EXPERIMENTAL FEATURES

*  Likelihood model `nbinomial2`




# INLA 18.10.28
## NEW EXPERIMENTAL FEATURES

*  New function `inla.priors.used`




# INLA 18.10.17
## NEW FEATURES

*  Export class `inla`




# INLA 18.10.16
## NEW EXPERIMENTAL FEATURES

*  New latent model: `dmatern`


## BUG FIXES

*  Improved the numerics for computing the scaling of the RW1 and RW2 models.





# INLA 18.10.09
## USER-VISIBLE CHANGES

*  New option `control.inla=list(tolerance.step=)`, to
      control the RMS of the step-size for the inner optimization.
*  Changed, slightly, the initial values for the exponent in
      the Weibull likelihood models, to a value close to zero instead of
      zero.
*  New vignette about how to deal with multinomial data.


## NEW EXPERIMENTAL FEATURES

*  Added option `verbose` to
      `inla.qsample()` and `inla.posterior.sample()`






# INLA 18.09.24
## BUG FIXES

*  Performance improvement using the PARDISO library




# INLA 18.09.21
## NEW EXPERIMENTAL FEATURES

*  Argument `selection` in `inla.posterior.sample`
      and `inla.qsample`.





# INLA 18.09.19
## BUG FIXES

*  Fix for `num.threads` in `inla.qinv()`




# INLA 18.09.18
## BUG FIXES

*  Allow better user control of sparse matrix library in
      `inla.qinv()`, `inla.qsample()` and
      `inla.posterior.sample()`




# INLA 18.09.14
## USER-VISIBLE CHANGES

*  New example added to `inla.posterior.sample()`
*  Slight changes in the default `print`, and
      `summary` for an `inla`-object


## BUG FIXES

*  Fixed the issue when `lincomb.derived.only=FALSE` and
      then using `inla.posterior.sample()`




# INLA 18.08.26
## USER-VISIBLE CHANGES

*  Added 32-bit builds for windows (upon request)


## NEW EXPERIMENTAL FEATURES

*  Added function `inla.posterior.sample.eval()`




# INLA 18.08.09
## USER-VISIBLE CHANGES

*  Added new function `inla.pardiso.check()`
*  Added COPYRIGHTS file


## NEW EXPERIMENTAL FEATURES

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


## BUG FIXES

*  Removed some code not used anymore





# INLA 18.07.27
## USER-VISIBLE CHANGES

*  NEWS page created (see `news(package="INLA")`)
*  Added vignette about the conditional logit model (thanks to
      Stefani Muff)
*  Fixed missprints in the documentation for model `ar1c`
      (Thanks to Virgilio Gomez Rubio)
*  Fixed documentation about argument `blas.num.threads` in `inla()`





# INLA 18.07.12
## USER-VISIBLE CHANGES

*  Package built with `R-3.5`, both stable and testing




# INLA 18.07.11
## USER-VISIBLE CHANGES

*  Package built for `R-3.4`, both stable and testing. 




# INLA 18.07.12
## NEW EXPERIMENTAL FEATURES

*  Likelihood model `pom` (proportional odds model)



