\newcommand{\PR}{\Sexpr[results=rd]{tools:::Rd_expr_PR(#1)}}

\name{NEWS}
\title{R-INLA News}
\encoding{UTF-8}

\section{CHANGES IN VERSION 22.02.16}{
  \subsection{NEW EXPERIMENTAL FEATURES}{
    \itemize{
      \item New feature \code{inla.cgeneric.q}
      \item New and updated PARDISO library
    }
  }
  \subsection{BUG FIXES}{
    \itemize{
      \item Code improvements and cleanup
    }
  }
}

\section{CHANGES IN VERSION 22.01.25}{
  \subsection{BUG FIXES}{
    \itemize{
      \item Fixed regression bug in \code{knmodels}
      \item Code improvements and cleanup
    }
  }
}

\section{CHANGES IN VERSION 22.01.19}{
  \subsection{BUG FIXES}{
    \itemize{
      \item Fixed bug in the variable expansion in some arguments
      that could give different results than before
    }
  }
}

\section{CHANGES IN VERSION 22.01.16}{
  \subsection{NEW EXPERIMENTAL FEATURES}{
    \itemize{
      \item Added new variant for nbinomial likelihood
    }
  }
}

\section{CHANGES IN VERSION 22.01.12}{
  \subsection{IMPROVEMENTS}{
    \itemize{
      \item Added \code{graph$cc$mean} to help defining intercepts for
      disconnected regions.

      \item Can now define the license key for PARDISO directly as
      \code{inla.setOption(pardiso.license="<KEY>")}

      \item Print total RAM in Gb in the inla-program output, for easier debugging
      of these problems.

      \item Using \code{inla.call="remote"} will now return error code that can be
      controlled by \code{try}, so that \code{r=try(inla(..., inla.call="remote"))}
      should work as intended.

      \item Add a refine step for \code{parallel.linesearch}, so if this is
      enabled, then it will do a restart at the end with this option
      disabled.

      \item Minor changes when option \code{safe} is enabled.
    }
  }
  \subsection{BUG FIXES}{
    \itemize{
      \item Fixed regression error with \code{plot(result, plot.prior=TRUE)}
      \item Fixed a regression error from version 21.12.21 which expanded some
      variables incorrect.
    }
  }
}

\section{CHANGES IN VERSION 21.12.21}{
  \subsection{NEW EXPERIMENTAL FEATURES}{
    \itemize{
      \item Option \code{safe} in function \code{inla}
    }
  }
  \subsection{BUG FIXES}{
    \itemize{
      \item Minor bug fixes
    }
  }
}

\section{CHANGES IN VERSION 21.12.18}{
  \subsection{BUG FIXES}{
    \itemize{
      \item Improved handling of \code{NAN} and \code{INF}
    }
  }
}

\section{CHANGES IN VERSION 21.12.16}{
  \subsection{BUG FIXES}{
    \itemize{
      \item A bug-fix for \code{matern2d} model
      \item Various code improvements
    }
  }
}

\section{CHANGES IN VERSION 21.12.11}{
  \subsection{BUG FIXES}{
    \itemize{
      \item Minor bug fixes for \code{cgeneric} for Windows
    }
  }
}

\section{CHANGES IN VERSION 21.12.10}{
  \subsection{NEW EXPERIMENTAL FEATURES}{
    \itemize{
      \item A \code{cgeneric} interface
      like \code{rgeneric} but with C-code.
      This is for internal testing only.
    }
  }
  \subsection{BUG FIXES}{
    \itemize{
      \item Various code improvements and minor fixes
    }
  }
}


\section{CHANGES IN VERSION 21.12.01}{
  \subsection{NEW FEATURES}{
    \itemize{
      \item Allow for probit cdf in likelihood model \code{pom}
    }
  }
  \subsection{BUG FIXES}{
    \itemize{
      \item Various code improvements
    }
  }
}

\section{CHANGES IN VERSION 21.11.29}{
  \subsection{NEW EXPERIMENTAL FEATURES}{
    \itemize{
      \item Small fixes
    }
  }
}

\section{CHANGES IN VERSION 21.11.28}{
  \subsection{NEW EXPERIMENTAL FEATURES}{
    \itemize{
      \item added new link \code{powerprobit}. Not completely done yet.
    }
  }
}

\section{CHANGES IN VERSION 21.11.22}{
  \subsection{BUG FIXES AND IMPROVEMENTS}{
    \itemize{
      \item Remove jemalloc-option that creates warning message for Mac
    }
  }
}

\section{CHANGES IN VERSION 21.11.21}{
  \subsection{BUG FIXES AND IMPROVEMENTS}{
    \itemize{
      \item Minor internal changes
      \item Add jemalloc options for Linux/Mac
    }
  }
}

\section{CHANGES IN VERSION 21.11.16}{
  \subsection{BUG FIXES AND IMPROVEMENTS}{
    \itemize{
      \item Minor internal changes for robustness
      \item Upgraded compilers to gcc-11 on Windows
    }
  }
}

\section{CHANGES IN VERSION 21.11.15}{
   \subsection{NEW EXPERIMENTAL FEATURES}{
     \itemize{
       \item Option for parallel linesearch has been renamed into
       \code{parallel.linesearch} in \code{control.inla}.
     }
   }
  \subsection{BUG FIXES AND IMPROVEMENTS}{
    \itemize{
      \item This is fix for a build-error in the Windows binary for
      version \code{21.11.14}
      \item Various minor internal code changes and improvements
    }
  }
}

\section{CHANGES IN VERSION 21.11.14}{
   \subsection{NEW EXPERIMENTAL FEATURES}{
     \itemize{
       \item Parallel linesearch implemented, see option
       \code{bfgs.version} in \code{control.inla}. This option is highly
       experimental at the moment and work in progress, and should
       not be used.
     }
   }
  \subsection{BUG FIXES AND IMPROVEMENTS}{
    \itemize{
      \item Change default settings for \code{family="tweedie"}
      \item Implemented an adaptive parallel/serial version of
        \code{Q%*%x}
      \item Improved perforance for parallel linear solve with \code{PARDISO}
      \item Use stable Legendre polynomial evaluation for spherical covariances
      \item Add support for \code{Fedora Linux 35}
      \item Fixed \code{dic.sat} when all response are \code{NA}
      \item Various internal code changes and improvements
    }
  }
}


\section{CHANGES IN VERSION 21.11.01}{
   \subsection{NEW EXPERIMENTAL FEATURES}{
     \itemize{
       \item New experimental feature \code{lp.scale}. No good documentation yet, but
       that will come soon.
     }
   }
  \subsection{BUG FIXES}{
    \itemize{
      \item Robustify \code{inla.iidkd.sample} function for numerical singular matrices
      \item Corrected \code{inla.stack.join} handling of factor variables that in some
            cases converted variables to/from factors unintentionally.
      \item \code{inla.mode="experimental"} got its fitted values back.
    }
  }
}

\section{CHANGES IN VERSION 21.10.03}{
  \subsection{BUG FIXES}{
    \itemize{
      \item Internal code improvement
    }
  }
}

\section{CHANGES IN VERSION 21.09.29-1}{
   \subsection{BUG FIXES}{
     \itemize{
       \item Better special number
     }
   }
   \subsection{NEW EXPERIMENTAL FEATURES}{
     \itemize{
       \item New option to \code{inla.iidkd.sample}
     }
   }
}

\section{CHANGES IN VERSION 21.09.29}{
   \subsection{NEW EXPERIMENTAL FEATURES}{
     \itemize{
       \item New latent model \code{iidkd} which use a different
       parameterisation than \code{iid3d} and similar ones.
     }
   }
}

\section{CHANGES IN VERSION 21.09.13}{
   \subsection{BUG FIXES}{
     \itemize{
       \item Fixed crash when using \code{int.strategy="user"} with
       \code{inla.mode="experimental"}
       \item Minor fixes
     }
   }
 }

\section{CHANGES IN VERSION 21.09.11}{
   \subsection{BUG FIXES}{
     \itemize{
       \item Fixed crash when using linear combinations with
       \code{inla.mode="experimental"}
       \item Minor fixes
     }
   }
 }

\section{CHANGES IN VERSION 21.09.01}{
   \subsection{NEW EXPERIMENTAL FEATURES}{
     \itemize{
       \item Added option \code{scale} to \code{family="xbinomial"}
     }
   }
 }

\section{CHANGES IN VERSION 21.08.31}{
   \subsection{BUG FIXES}{
     \itemize{
       \item Various bug-fixes
     }
   }
 }


\section{CHANGES IN VERSION 21.07.10}{
  \subsection{USER-VISIBLE CHANGES}{
    \itemize{
      
      \item If \code{control.predictor=list(compute=TRUE)} then the
      marginal densities for the linear predictor and fitted values are
      computed, but the new default is that only the summary statistics
      are returned not the marginal densities. This is because the
      storage requirement for the marginals can be substantial. The
      marginal densities can be returned setting
      \code{control.compute=list(return.marginals.predictor=TRUE)}

      \item The argument \code{control.results} are removed as it was
       not really in use.

     \item The new default value for \code{b.strategy} is
       \code{b.strategy="skip"}. See \code{?control.inla} for details
       
     }
   }

   \subsection{NEW FEATURES}{
     \itemize{
       \item New option \code{inla.timeout}, similar to
       \code{fmesher.timeout}; see \code{?inla.setOption}

       \item \code{CPO-marginals} are returned in when
       \code{control.compute=list(config=TRUE)}

       \item There is a change how a (numerical) singular constraint is
       treated, hopefully this new way is more robust.
     }
   }

   \subsection{NEW EXPERIMENTAL FEATURES}{
     \itemize{
       \item New option \code{inla.mode}, which define how to arrange the
       internal model formulation. One of \code{"classic"},
       \code{"twostage"} and \code{"experimental"}. The default is
       \code{"classic"}, which is unchanged behaviour compared for
       earlier versions. The other two are highly experimental for the
       moment. See \code{?inla}
     }
   }

   \subsection{BUG FIXES}{
     \itemize{
       \item Massive code clean-up and some minor fixes.
     }
   }
 }


\section{CHANGES IN VERSION 21.06.11}{
  \subsection{NEW EXPERIMENTAL FEATURES}{
    \itemize{
      \item Added support for Mac M1 (native build with R-4.1)
      \item New \code{twostage} options (highly experimental)
    }
  }
  \subsection{BUG FIXES}{
    \itemize{
      \item Massive code overhaul
      \item Various minor fixes
    }
  }
}

\section{CHANGES IN VERSION 21.05.02}{
  \subsection{BUG FIXES}{
    \itemize{
      \item Some minor bug-fixes and code improvements
    }
  }
}

\section{CHANGES IN VERSION 21.04.22}{
  \subsection{NEW EXPERIMENTAL FEATURES}{
    \itemize{
      \item Mac only: Added optional path to \code{vecLib}
      BLAS and LAPACK libraries.
      \item Adding argument \code{.special} to \code{inla.surv}-object.
    }
  }
  \subsection{BUG FIXES}{
    \itemize{
      \item Some code improvement 
    }
  }
}


\section{CHANGES IN VERSION 21.04.16}{
  \subsection{NEW EXPERIMENTAL FEATURES}{
    \itemize{
      \item Mac only: Link with \code{vecLib} BLAS and LAPACK libraries
            by default. Turn off with \code{inla.setOption(vecLib=FALSE)}.
      \item Adding an experimental improvement
            for \code{strategy=simplified.laplace} (not enabled by default)
    }
  }
  \subsection{BUG FIXES}{
    \itemize{
      \item Some code improvement 
    }
  }
}

\section{CHANGES IN VERSION 21.03.31}{
  \subsection{BUG FIXES}{
    \itemize{
      \item General code improvement 
    }
  }
}

\section{CHANGES IN VERSION 21.03.27}{
  \subsection{NEW EXPERIMENTAL FEATURES}{
    \itemize{
      \item New family \code{stochvolsn}
      \item New family \code{cenpoisson2}
    }
  }
  \subsection{BUG FIXES}{
    \itemize{
      \item Improve robustness for survival models
      \item General code improvement
    }
  }
}

\section{CHANGES IN VERSION 21.03.21}{
  \subsection{NEW EXPERIMENTAL FEATURES}{
    \itemize{
      \item Changes in \code{inla.jp}; see help-page for details.
    }
  }
}

\section{CHANGES IN VERSION 21.03.20}{
  \subsection{NEW EXPERIMENTAL FEATURES}{
    \itemize{
      \item New family \code{gompertz}
      \item \code{inla.binary.install()} now do md5-checksum check
    }
  }
  \subsection{BUG FIXES}{
    \itemize{
      \item Code optimization for \code{barrier} models
    }
  }
}


\section{CHANGES IN VERSION 21.03.17}{
  \subsection{BUG FIXES}{
    \itemize{
      \item Internal code changes for \code{family=tweedie}
    }
  }
}

\section{CHANGES IN VERSION 21.03.16}{
  \subsection{NEW EXPERIMENTAL FEATURES}{
    \itemize{
      \item \code{inla.binary.install()} can now also run
            non-interactive
    }
  }
  \subsection{BUG FIXES}{
    \itemize{
      \item Internal code changes
    }
  }
}

\section{CHANGES IN VERSION 21.03.14}{
  \subsection{BUG FIXES}{
    \itemize{
      \item Internal code changes
    }
  }
}

\section{CHANGES IN VERSION 21.03.13}{
  \subsection{BUG FIXES}{
    \itemize{
      \item Fixed some missprints in doc
      \item Trying to get vigettes to work again
    }
  }
}

\section{CHANGES IN VERSION 21.03.08}{
  \subsection{NEW EXPERIMENTAL FEATURES}{
    \itemize{
      \item Added family \code{logperiodogram} back in
      \item Added new family \code{agaussian}
    }
  }
  \subsection{BUG FIXES}{
    \itemize{
      \item Improved initial values for \code{x} 
      \item Improved computations for saturated likelihood
      \item Various minor code cleanup
    }
  }
}

\section{CHANGES IN VERSION 21.02.23}{
  \subsection{NEW FEATURES}{
    \itemize{
      \item Default \code{strategy} is now \code{simplified.laplace}
      for smaller models (< 5000 nodes) and \code{adaptive} for larger ones
    }
  }
  \subsection{BUG FIXES}{
    \itemize{
      \item Minor bug-fixes
    }
  }
}

\section{CHANGES IN VERSION 21.01.26}{
  \subsection{BUG FIXES}{
    \itemize{
      \item Improve parallel rerodering for PARDISO
      \item Improve \code{coxph} for large data
      \item Improve internal CPU-timing output
    }
  }
}

\section{CHANGES IN VERSION xx.xx.xx}{
  \subsection{USER-VISIBLE CHANGES}{
    \itemize{
      \item 
    }
  }
  \subsection{NEW FEATURES}{
    \itemize{
      \item 
    }
  }
  \subsection{NEW EXPERIMENTAL FEATURES}{
    \itemize{
      \item
    }
  }
  \subsection{BUG FIXES}{
    \itemize{
      \item 
    }
  }
}

\section{CHANGES IN VERSION 21.01.18}{
  \subsection{BUG FIXES}{
    \itemize{
      \item Use METIS5 with PARDISO
    }
  }
}

\section{CHANGES IN VERSION 21.01.13}{
  \subsection{NEW EXPERIMENTAL FEATURES}{
    \itemize{
      \item Family \code{fmri}
    }
  }
  \subsection{BUG FIXES}{
    \itemize{
      \item Speed improvements for family \code{tweedie}
    }
  }
}

\section{CHANGES IN VERSION 21.01.08}{
  \subsection{NEW EXPERIMENTAL FEATURES}{
    \itemize{
      \item New family \code{tweedie}
    }
  }
  \subsection{BUG FIXES}{
    \itemize{
      \item Optimization work, code improvement and maintenance
    }
  }
}

\section{CHANGES IN VERSION 20.12.10}{
  \subsection{NEW EXPERIMENTAL FEATURES}{
    \itemize{
      \item Added a new experimental feature for running PARDISO in parallel
    }
  }
  \subsection{BUG FIXES}{
    \itemize{
      \item Various minor fixes and improvements
    }
  }
}

\section{CHANGES IN VERSION 20.12.04}{
  \subsection{BUG FIXES}{
    \itemize{
      \item Bug-fix release only
    }
  }
}

\section{CHANGES IN VERSION 20.11.29}{
  \subsection{USER-VISIBLE CHANGES}{
    \itemize{
      \item Added new section to the vignette about \code{rgeneric}
    }
  }
  \subsection{NEW FEATURES}{
    \itemize{
      \item Option \code{optimize} in \code{inla.rgeneric.define}
    }
  }
  \subsection{BUG FIXES}{
    \itemize{
      \item Cleanup in some output format.
      \item Various internal changes. 
    }
  }
}

\section{CHANGES IN VERSION 20.11.22}{
  \subsection{BUG FIXES}{
    \itemize{
      \item Code improvement and optimization for \code{rgeneric}.
      \item Some OpenMP improvements
    }
  }
}

\section{CHANGES IN VERSION 20.11.18}{
  \subsection{NEW FEATURES}{
    \itemize{
      \item New family \code{gammajw}
    }
  }
  \subsection{BUG FIXES}{
    \itemize{
      \item Fixed an OpenMP deadlock case, using \code{rgeneric}
    }
  }
}

\section{CHANGES IN VERSION 20.11.16}{
  \subsection{USER-VISIBLE CHANGES}{
    \itemize{
      \item PARDISO version 7 is included for Mac and Linux.
      \item Improved paralellism and nested paralellism
      \item Speedup improvements, especially for models with many
      constraints.
      \item Improved default \code{plot(result)}. Argument
      \code{cex=..} will now work
    }
  }
  \subsection{NEW FEATURES}{
    \itemize{
      \item New shorthand feature in \code{inla.posterior.sample.eval}
      for extracing model components from samples
      \item Added \code{fmesher.timeout} option, see \code{?inla.setOption}
    }
  }
  \subsection{BUG FIXES}{
    \itemize{
      \item Using a non-zero seed in \code{inla.posterior.sample} will
      now force serial computations.
      \item A lot of internal code cleanup and improvements
    }
  }
}

\section{CHANGES IN VERSION 20.10.11}{
  \subsection{NEW EXPERIMENTAL FEATURES}{
    \itemize{
      \item New likelihood \code{poisson.special1}
    }
  }
  \subsection{BUG FIXES}{
    \itemize{
      \item Minor code cleanup
      \item Minor changes in code/doc to adapt to R-4.0
    }
  }
}

\section{CHANGES IN VERSION 20.09.25}{
  \subsection{BUG FIXES}{
    \itemize{
      \item Bug fixes and some optimisation work
    }
  }
}

\section{CHANGES IN VERSION 20.09.16}{
  \subsection{NEW EXPERIMENTAL FEATURES}{
    \itemize{
      \item Likelihood model \code{sn} is redone.
    }
  }
  \subsection{BUG FIXES}{
    \itemize{
      \item Various minor changes.
    }
  }
}

\section{CHANGES IN VERSION 20.09.13}{
  \subsection{NEW FEATURES}{
    \itemize{
      \item Model \code{iid} added to possible models
      for the baseline hazard in the Cox-ph model
    }
  }
  \subsection{NEW EXPERIMENTAL FEATURES}{
    \itemize{
      \item Likelihood models \code{sn} and \code{sn2} have been
      replaced with a new version \code{sn}. The parameterisation is
      different and now done correct.
    }
  }
}

\section{CHANGES IN VERSION 20.09.09}{
  \subsection{USER-VISIBLE CHANGES}{
    \itemize{
      \item The likelihood models \code{sn} and \code{sn2}
      are now disabled. They need a rewrite to be done right.
    }
  }
  \subsection{NEW EXPERIMENTAL FEATURES}{
    \itemize{
      \item For Linux, then build so-libraries are loaded before the
      system ones. This behaviour can be reverted by setting the
      environment variable \code{INLA_NATIVE_LD_LIBRARY_PATH}.
    }
  }
}

\section{CHANGES IN VERSION 20.09.07}{
  \subsection{USER-VISIBLE CHANGES}{
    \itemize{
      \item Dimension for CCD integration have been
      increased from 38 to 52.
    }
  }
  \subsection{NEW EXPERIMENTAL FEATURES}{
    \itemize{
      \item Default prior for the \code{intercept} and the \code{skewness} in the
      \code{skew-normal link} have changed. 
    }
  }
  \subsection{BUG FIXES}{
    \itemize{
      \item Various minor changes
    }
  }
}

\section{CHANGES IN VERSION 20.08.31}{
  \subsection{BUG FIXES}{
    \itemize{
      \item A small bug-fix with argument \code{-t}
      \item Some improvements in the install scripts
    }
  }
}

\section{CHANGES IN VERSION 20.08.30}{
  \subsection{BUG FIXES}{
    \itemize{
      \item More work on the nested parallelism stuff
    }
  }
}

\section{CHANGES IN VERSION 20.08.29}{
  \subsection{BUG FIXES}{
    \itemize{
      \item Argument \code{num.threads="A"} now means
      \code{num.threads="A:1"}
      \item Some changes for nested paralellism
    }
  }
}

\section{CHANGES IN VERSION 20.08.28}{
  \subsection{USER-VISIBLE CHANGES}{
    \itemize{
      \item There is a minor change in \code{inla.qsample} and \code{inla.posterior.sample}.
      If argument \code{seed!=0} then serial mode is now forced.
      Earlier, an error would be raised if parallel mode was also requested.
    }
  }
  \subsection{NEW EXPERIMENTAL FEATURES}{
    \itemize{
      \item Argument \code{remove.names} in \code{control.fixed}
    }
  }
  \subsection{BUG FIXES}{
    \itemize{
      \item Bug-fixes for nested parallism
      \item Some optimization improvement
    }
  }
}

\section{CHANGES IN VERSION 20.08.24}{
  \subsection{BUG FIXES}{
    \itemize{
      \item Some minor fixes in the install
    }
  }
  \subsection{MISC}{
    \itemize{
      \item Linux builds now link with PARDISO version 7 (beta)
    }
  }
}

\section{CHANGES IN VERSION 20.08.22}{
  \subsection{NEW EXPERIMENTAL FEATURES}{
    \itemize{
      \item MKL is now default for both Mac and Linux.
      If this cause any issues, revert back by setting
      \code{inla.setOption(mkl=FALSE)}.
      \item As of today, the TAUCS library does not play 
      well with MKL on Mac for some reason, and this is
      silent accounted for.
      It is possible to bypass the internal check,
      and if that becomes in issue, just define
      \code{mkl=FALSE} as above.
    }
  }
  \subsection{BUG FIXES}{
    \itemize{
      \item Updated some install scripts
    }
  }
}

\section{CHANGES IN VERSION 20.08.21}{
  \subsection{BUG FIXES}{
    \itemize{
      \item Various minor fixes
      \item MKL binaries are now linked with libiomp5
    }
  }
}

\section{CHANGES IN VERSION 20.08.19}{
  \subsection{BUG FIXES}{
    \itemize{
      \item Changed to the gcc-10 compiler suite for Ubuntu 1804
      and 2004.
      \item Fixed an issue if the optimiser does not move
      and the use of directions
      \item Minor fixes
    }
  }
}

\section{CHANGES IN VERSION 20.08.18}{
  \subsection{USER-VISIBLE CHANGES}{
    \itemize{
      \item Code is prepared for nested parallism, which most gain
      obtained when using the PARDISO library
      \item Argument \code{num.threads} is now in the format \code{A:B},
      where A threads for the outer layer and B threads for the inner
      layer. And this applies for several functions.
      \item The speed have improved, mostly when using the PARDISO
      library. Hopefully, nothing is broken.
    }
  }
  \subsection{NEW FEATURES}{
    \itemize{
      \item likelihood \code{beta} now allow for cencoring near 0 and 1,
      see \code{inla.doc("beta", section="likelihood")}
    }
  }
  \subsection{NEW EXPERIMENTAL FEATURES}{
    \itemize{
      \item Default value for \code{control.inla$h} is now 0.005, which
      replace the previous value of 0.01. We might change this later,
      but we have enough experience to know that 0.01 is slightly to
      large.
      \item \code{control.inla$optimise.strategy="smart"} is now
      default. By construction, this should be both faster and more
      robust. Using \code{"plain"} reverts back the old behaviour. 
      \item \code{control.inla$use.directions=TRUE} is now default. This option
      estimate numerical gradients/Hessian in directions with
      more change than for coordinate-wise directions.
    }
  }
}

\section{CHANGES IN VERSION 20.08.11}{
  \subsection{NEW EXPERIMENTAL FEATURES}{
    \itemize{
      \item Add new experimental feature
      \code{control.inla=list(use.directions=TRUE/FALSE)}
    }
  }
}

\section{CHANGES IN VERSION 20.08.09}{
  \subsection{NEW EXPERIMENTAL FEATURES}{
    \itemize{
      \item New experimental optimise strategy
      \code{control.inla=list(optmise.strategy="smart"}
      that is hopefully faster and as safe as the default one. Maybe
      even more safe and robust.
    }
  }
  \subsection{BUG FIXES}{
    \itemize{
      \item More work on the nested paralellism 
    }
  }
}

\section{CHANGES IN VERSION 20.08.06}{
  \subsection{NEW EXPERIMENTAL FEATURES}{
    \itemize{
      \item More work on nested parallelism and some code changes
      \item Code changes to reduce memory usage
      \item Linux/Mac binaries are now linked with \code{jemalloc}
    }
  }
}

\section{CHANGES IN VERSION 20.08.04}{
  \subsection{NEW EXPERIMENTAL FEATURES}{
    \itemize{
      \item Improved nested parallelism and some code changes
    }
  }
}

\section{CHANGES IN VERSION 20.08.03}{
  \subsection{NEW EXPERIMENTAL FEATURES}{
    \itemize{
      \item More work on nested parallelism
    }
  }
}

\section{CHANGES IN VERSION 20.08.02}{
  \subsection{USER-VISIBLE CHANGES}{
    \itemize{
      \item Option \code{control.inla$lincomb.derived.only} is now disabled.
    }
  }
  \subsection{NEW EXPERIMENTAL FEATURES}{
    \itemize{
      \item Testing nested parallelism \code{openmp.nested} with
      \code{num.threads="A,B"}. Work in progress
    }
  }
}

\section{CHANGES IN VERSION 20.07.27}{
  \subsection{NEW FEATURES}{
    \itemize{
      \item \code{plot(result)} will now produce a plot of the
      \code{CPO/PIT} for each likelihood family (if available), instead
      of a joint plot as earlier.
      \item New vignette about \code{jmarginal}
    }
  }
  \subsection{NEW EXPERIMENTAL FEATURES}{
    \itemize{
      \item Added new likelihood models \code{zeroinflatedcenpoisson0}
      and \code{zeroinflatedcenpoisson1}
      \item Link-model \code{sn} is updated, as well as the
      PC-prior for the skewness therein, and the added intercept model.
    }
  }
  \subsection{BUG FIXES}{
    \itemize{
      \item Revised the PIT calculations for family \code{cenpoisson}
      \item Code rewrite to (try to) prevent \code{Inf} for \code{DIC} calculations
      \item Minor fixed in \code{inla.binary.install}
    }
  }
}

\section{CHANGES IN VERSION 20.07.18}{
  \subsection{BUG FIXES}{
    \itemize{
      \item Build-script changes and misc fixes
    }
  }
}

\section{CHANGES IN VERSION 20.07.16}{
  \subsection{USER-VISIBLE CHANGES}{
    \itemize{
      \item Packages \code{mpoly} and \code{symmoments} are
      added to the Suggests-list.
    }
  }
  \subsection{NEW FEATURES}{
    \itemize{
      \item Added info about \code{inla.prune()} to the
      startup message
    }
  }
  \subsection{NEW EXPERIMENTAL FEATURES}{
    \itemize{
      \item More features for \code{jmarginal} added
    }
  }
  \subsection{BUG FIXES}{
    \itemize{
      \item Build-scripts fixes
      \item Fix for a rare \code{fmesher} issue
      \item Improved the code for the DIC calculations to make
      them more stable
      \item Some improvment in the PC-prior for the SN-link
    }
  }
}

\section{CHANGES IN VERSION 20.07.12}{
  \subsection{BUG FIXES}{
    \itemize{
      \item Fixed \code{inla.link.sn} to vectorise over argument
      \code{a} as \code{sn}-package do not do that properly itself. 
    }
  }
}

\section{CHANGES IN VERSION 20.07.09}{
  \subsection{NEW FEATURES}{
    \itemize{
      \item Package built with R-4.0
    }
  }
}

\section{CHANGES IN VERSION 20.07.04}{
  \subsection{USER-VISIBLE CHANGES}{
    \itemize{
      \item Vignette added for \code{family=bGEV}
      \item Some internal changes due to the migration to \code{git}
    }
  }
}
\section{CHANGES IN VERSION 20.06.29}{
  \subsection{NEW EXPERIMENTAL FEATURES}{
    \itemize{
      \item Added \code{Qprior.diag} to the output when \code{config=TRUE}.
      The off-diagonals of this matrix are the same as \code{Q} in the
      same configuration, so only the diagonal of \code{Qprior}
      is stored.
      \item Added some internal experimental code
    }
  }
  \subsection{BUG FIXES}{
    \itemize{
      \item PARDISO interface: internal check added
      \item Fixed an wrong assert with family=bgev
    }
  }
}

\section{CHANGES IN VERSION 20.06.22}{
  \subsection{BUG FIXES}{
    \itemize{
      \item Improve some code in the PARDISO interface
      \item Improved the computation of the third derivative
      in the log likelihood.
    }
  }
}



\section{CHANGES IN VERSION 20.06.18}{
  \subsection{BUG FIXES}{
    \itemize{
      \item Improved WKT support for PROJ6/GDAL3
    }
  }
}

\section{CHANGES IN VERSION 20.05.16}{
  \subsection{NEW FEATURES}{
    \itemize{
      \item Support for PROJ6/RGDAL3 for handling CRS information for
        spatial objects.
    }
  }
}

\section{CHANGES IN VERSION 20.06.15}{
  \subsection{NEW EXPERIMENTAL FEATURES}{
    \itemize{
      \item For the \code{intslope}-model: made all \code{gamma}'s
      default fixed to 1, so its similar in style the copy-feature.
      \item Added argument \code{constr} it \code{inla.rjmarginal}
      \item Added argument \code{ask} to function \code{inla.prune}
    }
  }
  \subsection{BUG FIXES}{
    \itemize{
      \item Argument \code{cyclic=TRUE} in \code{f()} should not set
            \code{constr=FALSE} when default is \code{constr=TRUE}
      \item Change the \code{scale.model=TRUE} code for \code{RW1/RW2} so the
      scaling for the continous case is the same as for the discrete
      case when the locations are eqvidistant.
      \item Disable link \code{sslogit}
    }
  }
}

\section{CHANGES IN VERSION 20.05.12}{
  \subsection{BUG FIXES}{
    \itemize{
      \item Fixed an issue with model \code{besag2}
      \item Fixed an issue with \code{plot(r,plot.prior=TRUE)} for some priors
    }
  }
}

\section{CHANGES IN VERSION 20.05.04}{
  \subsection{NEW FEATURES}{
    \itemize{
      \item Remove the experimental status of \code{inla.posterior.sample.eval}
    }
  }
  \subsection{NEW EXPERIMENTAL FEATURES}{
    \itemize{
      \item Added function \code{inla.prune} which will remove binaries
      not supported by the running OS, to reduce the size of the
      package.
      \item Added method \code{summary} and \code{print} to class
      \code{inla.jmarginal}
    }
  }
  \subsection{BUG FIXES}{
    \itemize{
      \item Add check for \code{NA/NaN/Inf} in mesh creation input
      ocations
      \item Make sure that skewness is not to high in \code{inla.posterior.sample}
    }
  }
}

\section{CHANGES IN VERSION 20.04.18}{
  \subsection{NEW FEATURES}{
    \itemize{
      \item Added new argument \code{tag} to \code{inla.coxph}
    }
  }
  \subsection{NEW EXPERIMENTAL FEATURES}{
    \itemize{
      \item \code{inla.rjmarginal.eval}, to evaluate samples from a join
      approximations
    }
  }
  \subsection{BUG FIXES}{
    \itemize{
      \item Names of samples are now "sample:1", "sample:2", and should
      be coherent over all functions. Similar, their contents, its like
      "x:1", "x:2", etc.
      \item Fixed a bug setting prior for the log baseline hazard in \code{inla.coxph}
    }
  }
}

\section{CHANGES IN VERSION 20.04.14}{
  \subsection{BUG FIXES}{
    \itemize{
      \item Small fix so that \code{result$mode$x} is written out in the
      case where \code{nhyper=0} and \code{num_threads>1}
      \item Minor internal changes.
    }
  }
}

\section{CHANGES IN VERSION 20.04.06}{
  \subsection{NEW EXPERIMENTAL FEATURES}{
    \itemize{
      \item Added link \code{loga}. Not yet documented.
      \item First try on a new feature to more easily approximate
      the joint marginal for a subset of the latent field. This is a new
      option \code{selection} and corresponding \code{inla.rjmarginal()}
      to sample from it. 
    }
  }
  \subsection{BUG FIXES}{
    \itemize{
      \item Added check that \code{model="linear"} is not used with
      replicate or group, which is not intention.
    }
  }
}

\section{CHANGES IN VERSION 20.03.29}{
  \subsection{USER-VISIBLE CHANGES}{
    \itemize{
      \item \code{MCMC} mode is now disabled
    }
  }
  \subsection{NEW FEATURES}{
    \itemize{
      \item Skewness correction is now back as default, in 
      \code{inla.posterior.sample()}
    }
  }
  \subsection{NEW EXPERIMENTAL FEATURES}{
    \itemize{
      \item Added family \code{xbinomial} that allow non-integer
      response.
      \item Likelihood model \code{bgev} add (not yet complete), and was
      renamed from the experimental likelihood model \code{gev2}.
    }
  }
  \subsection{BUG FIXES}{
    \itemize{
      \item If \code{inla.call="remote"} is set,
      then \code{INLA:::inla.call.builtin()} is used
      if \code{inla.qinv()} and/or \code{inla.qsolve()} are
      used while constructing the model.
    }
  }
}

\section{CHANGES IN VERSION 20.03.17}{
  \subsection{BUG FIXES}{
    \itemize{
      \item Updated file \code{jointdataCD4.rds} in \code{exampledata/}
    }
  }
}


\section{CHANGES IN VERSION 20.03.09}{
  \subsection{BUG FIXES}{
    \itemize{
      \item Fixed a bug in the PIT calculations for the zeroinflated,
      type 0, of poisson, binomial and nbinomial.
    }
  }
}

\section{CHANGES IN VERSION 20.03.08}{
  \subsection{NEW EXPERIMENTAL FEATURES}{
    \itemize{
      \item Added option \code{b.strategy} in \code{control.inla} to
      control what to do with the linear term when the \code{cmin} option is
      in effect
      \item Added in-interval observed event in \code{inla.surv}
    }
  }
  \subsection{BUG FIXES}{
    \itemize{
      \item Added \code{dplyr} as suggested package as
      \code{dplyr::bind_rows} can replace
      \code{INLA::inla.rbind.data.frames}
    }
  }
}

\section{CHANGES IN VERSION 20.02.19}{
  \subsection{USER-VISIBLE CHANGES}{
    \itemize{
      \item Added argument \code{E}, or \code{log(offset)}, to
      likelihood \code{gammacount}, so its equal to family \code{poisson}
      for \code{alpha=1}.
    }
  }
  \subsection{BUG FIXES}{
    \itemize{
      \item Minor changes
    }
  }
}

\section{CHANGES IN VERSION 20.01.25}{
  \subsection{USER-VISIBLE CHANGES}{
    \itemize{
      \item Added a check that discrete observations are indeed
      integers, like for Poisson, Binomial, etc
    }
  }
  \subsection{NEW FEATURES}{
    \itemize{
      \item The function \code{inla.binary.install} is now exported.
      \item Added new likelihood family, \code{xpoisson}, which allows
      continous response: see the documentation for details (and note
      the error-check now done for discrete observations)
    }
  }
  \subsection{NEW EXPERIMENTAL FEATURES}{
    \itemize{
      \item Added new likelihood \code{dgp} (discrete generalized Pareto)
    }
  }
  \subsection{BUG FIXES}{
    \itemize{
      \item Code clean-up (\code{contpoisson} and \code{qcontpoisson})
      \item Made \code{inla.pardiso.check()} a bit more informative if
      there is an error.
    }
  }
}

\section{CHANGES IN VERSION 19.12.10}{
  \subsection{BUG FIXES}{
    \itemize{
      \item Improved documentation of \code{inla.posterior.sample} and
      \code{inla.coxph}
      \item Fixed an issue with \code{NA} data in the family \code{gev2}
    }
  }
}

\section{CHANGES IN VERSION 19.12.03}{
  \subsection{BUG FIXES}{
    \itemize{
      \item Updated some documentation about the \code{pc.gevtail} prior.
      \item Reverted \code{inla.posterior.sample} back to the old
      version, the new experimental version is available as
      \code{INLA:::inla.posterior.sample.new}
      \item Error in \code{Epil} data-set, \code{y[31]} should be
      23 not 21. 
    }
  }
}

\section{CHANGES IN VERSION 19.11.17}{
  \subsection{USER-VISIBLE CHANGES}{
    \itemize{
      \item Updated the vignette about the multinomial distribution
    }
  }
  \subsection{NEW EXPERIMENTAL FEATURES}{
    \itemize{
      \item New experimental windows binary built with
  \code{x86_64-w64-mingw32-gcc}, version 7.3, and linked with the
  pardiso library. Its stored in \code{bin/windows/experimental}
    }
  }
  \subsection{BUG FIXES}{
    \itemize{
      \item Updated \code{inla.qreordering} and updated \code{leuk-demo.R}
      example file (and the corresponding zip-file).
    }
  }
}

\section{CHANGES IN VERSION 19.11.10}{
  \subsection{NEW EXPERIMENTAL FEATURES}{
    \itemize{
      \item Cache values of \code{qgamma} to speedup
      Gamma quantile regression
    }
  }
}

\section{CHANGES IN VERSION 19.10.30}{
  \subsection{USER-VISIBLE CHANGES}{
    \itemize{
      \item Added a scaling constant for the precision parameter in the
      \code{qkumar} likelihood (to avoid instabilities). See updated
      documentation for details.
    }
  }
  \subsection{NEW FEATURES}{
    \itemize{
      \item \code{inla.posterior.sample} now correct for possible skewness
      by default: see \code{?inla.posterior.sample} for details.
    }
  }
}


\section{CHANGES IN VERSION 19.10.16}{
  \subsection{NEW EXPERIMENTAL FEATURES}{
    \itemize{
      \item Likelihoodmodel \code{betabinomialna}
    }
  }
}


\section{CHANGES IN VERSION 19.10.15}{
  \subsection{USER-VISIBLE CHANGES}{
    \itemize{
      \item Default prior for the tail parameter in likelihood model
      \code{gp}, have changed to \code{pc.gevtail}, and the name change
      from \code{shape} to \code{tail}. It is now required to define a
      interval for the tail parameter, similar to \code{pc.gevtail}.
    }
  }
}

\section{CHANGES IN VERSION 19.10.06}{
  \subsection{BUG FIXES}{
    \itemize{
      \item Code-improvement for the \code{loggamma}-function
      \item \code{barrier.R} updated (minor fix and code edits)
    }
  }
}

\section{CHANGES IN VERSION 19.10.02}{
  \subsection{BUG FIXES}{
    \itemize{
      \item Disable some debug output
    }
  }
}

\section{CHANGES IN VERSION 19.10.01}{
  \subsection{BUG FIXES}{
    \itemize{
      \item Fixed a bug in the \code{nmixnb} likelihood. 
      \item Preserve names in \code{inla.posterior.sample.eval}
      if present.
    }
  }
}

\section{CHANGES IN VERSION 19.09.18}{
  \subsection{BUG FIXES}{
    \itemize{
      \item More work on the skew-normal link model
    }
  }
}

\section{CHANGES IN VERSION 19.09.15}{
  \subsection{NEW EXPERIMENTAL FEATURES}{
    \itemize{
      \item \code{INLA:::inla.binary.install()} is a new interactive tool
      to install alternative Linux builds. 
    }
  }
}

\section{CHANGES IN VERSION 19.09.10}{
  \subsection{NEW EXPERIMENTAL FEATURES}{
    \itemize{
      \item Added skew-normal link-model \code{sn} for binary data,
            with its PC-prior
    }
  }
}


\section{CHANGES IN VERSION 19.09.03}{
  \subsection{NEW EXPERIMENTAL FEATURES}{
    \itemize{
      \item Added \code{robit} link model.
    }
  }
  \subsection{BUG FIXES}{
    \itemize{
      \item Improved the stability of the saturated deviance
      calculations
      \item Fixed \code{INLA:::inla.is.list.of.lists} to cover the
      case where the arguments are a list of named lists
    }
  }
}


\section{CHANGES IN VERSION 19.07.27}{
  \subsection{NEW FEATURES}{
    \itemize{
      \item New (experimental) likelihood: gev2
    }
  }
  \subsection{BUG FIXES}{
    \itemize{
      \item Fixed, again, an issue with (parallel) PARDISO
      and many linear combinations.
      \item Minor code changes in \code{doc.R}
    }
  }
}

\section{CHANGES IN VERSION 19.07.21}{
  \subsection{USER-VISIBLE CHANGES}{
    \itemize{
      \item Removed must-be-enabled warnings in some
      surival models, from Oct  25 2017
    }
  }
  \subsection{NEW FEATURES}{
    \itemize{
      \item Added PC-prior for the Weibull likelihood models. The prior
      is derived
      for \code{variant = 1}, which is the good parameterisation.
    }
  }
  \subsection{BUG FIXES}{
    \itemize{
      \item Added missing \code{to.theta} and \code{from.theta}
      functions in likelihoods \code{sn} and \code{sn2}
      \item Fix some documentation in \code{marginal.R} (refering to the
      obsolete function \code{inla.marginal.transform})
      \item Fixed an issue with (parallel) PARDISO and many linear combinations.
    }
  }
}

\section{CHANGES IN VERSION 19.05.19}{
  \subsection{BUG FIXES}{
    \itemize{
      \item Set \code{StagedInstall:no} to work around
      installation problems for MacOS and R-3.6
    }
  }
}

\section{CHANGES IN VERSION 19.05.17}{
  \subsection{USER-VISIBLE CHANGES}{
    \itemize{
      \item The internal parameterisation of the alpha-parameter for the
  Weibull likelihood familes, has been redefined/scaled, to fix some
  optimisation issues. This means that the default prior has changed (a
  little) and user-defined priors has to change to account for this new
  internal parameterisation (sorry about that). See the documentation
  for details.
    }
  }
}

\section{CHANGES IN VERSION 19.05.16}{
  \subsection{NEW FEATURES}{
    \itemize{
      \item Option \code{short.summary} will use a version of
      \code{summary} with less output, maybe more suitable for
      Markdown documents. 
    }
  }
}

\section{CHANGES IN VERSION 19.05.13}{
  \subsection{NEW EXPERIMENTAL FEATURES}{
    \itemize{
      \item Added exampledata directory for various example datasets
    }
  }
  \subsection{BUG FIXES}{
    \itemize{
      \item Code cleanup and improved some input-error checking.
    }
  }
}

\section{CHANGES IN VERSION 19.04.16}{
  \subsection{BUG FIXES}{
    \itemize{
      \item Fixed an error in the cache-system for
      \code{model="rgeneric"} and \code{model="dmatern"}.
      Most notably with option \code{openmp.strategy="pardiso.parallel"}.
    }
  }
}

\section{CHANGES IN VERSION 19.04.14}{
  \subsection{BUG FIXES}{
    \itemize{
      \item Removed the weight correction for the computation of the cpo
            for \code{int.design="user.expert"}
    }
  }
}

\section{CHANGES IN VERSION 19.04.09}{
  \subsection{NEW EXPERIMENTAL FEATURES}{
    \itemize{
      \item Option \code{int.strategy="user.expert"}, see the vignette
      about user-defined integration points.
      \item Merge also \code{cpo} and \code{po} results in \code{inla.merge()}
    }
  }
}

\section{CHANGES IN VERSION 19.04.01}{
  \subsection{BUG FIXES}{
    \itemize{
      \item Fixed an issue with AR-model and group
    }
  }
}

\section{CHANGES IN VERSION 19.03.16}{
  \subsection{BUG FIXES}{
    \itemize{
      \item Small fix in model \code{dmatern}
    }
  }
}


\section{CHANGES IN VERSION 19.03.04}{
  \subsection{BUG FIXES}{
    \itemize{
      \item Redirect error output of some warning messages in the remote-feture
      section from MacOSX to Linux.
      \item Faster return when \code{mu} is zero for \code{rgeneric}
    }
  }
}

\section{CHANGES IN VERSION 19.03.02}{
  \subsection{BUG FIXES}{
    \itemize{
      \item Changed from PARDISO.PARALLEL to PARDISO.SERIAL in \code{inla.qsample}
      \item Optimize the \code{nhrs} for \code{inla.qsolve} for PARDISO
    }
  }
}

\section{CHANGES IN VERSION 19.02.28}{
  \subsection{BUG FIXES}{
    \itemize{
      \item Several fgn-models are now fine
      \item Fixed CPU timing with the PARDISO library
    }
  }
}

\section{CHANGES IN VERSION 19.02.26}{
  \subsection{BUG FIXES}{
    \itemize{
      \item Do not need to optimize reordering when PARDISO is used. 
    }
  }
}

\section{CHANGES IN VERSION 19.02.17}{
  \subsection{BUG FIXES}{
    \itemize{
      \item Fixed input-test using \code{inla.qsample} with
      \code{selection}-argument.
      \item Added back \code{family = "normal"} which is now
      translated to \code{family = "gaussian"} internally. 
    }
  }
}

\section{CHANGES IN VERSION 19.02.14}{
  \subsection{BUG FIXES}{
    \itemize{
      \item More work and fixes in  \code{inla.merge}
    }
  }
}

\section{CHANGES IN VERSION 19.02.12}{
  \subsection{USER-VISIBLE CHANGES}{
    \itemize{
      \item Simplied \code{print.inla} output
    }
  }
  \subsection{NEW EXPERIMENTAL FEATURES}{
    \itemize{
      \item New method \code{merge} and function \code{inla.merge},
      for merging \code{inla}-objects
    }
  }
  \subsection{BUG FIXES}{
    \itemize{
      \item Store \code{control.family} after processing, in the
      \code{result$.args} argument, not just the calling value.
    }
  }
}


\section{CHANGES IN VERSION 19.02.09}{
  \subsection{NEW FEATURES}{
    \itemize{
      \item New parameter for Gaussian likelihood: Fixed offset in the
      variance. 
    }
  }
  \subsection{BUG FIXES}{
    \itemize{
      \item Updated \code{envir} definition in the \code{rgeneric}
      documentation and examples. 
    }
  }
}

\section{CHANGES IN VERSION 19.02.06}{
  \subsection{USER-VISIBLE CHANGES}{
    \itemize{
      \item Removed testing code for likelihood model \code{testbinomial1}
    }
  }
  \subsection{NEW EXPERIMENTAL FEATURES}{
    \itemize{
      \item Added new likelihood \code{gamma.surv}
    }
  }
  \subsection{BUG FIXES}{
    \itemize{
      \item Cleaned up the use of temporary dir and files
      \item General code clean-up
    }
  }
}

\section{CHANGES IN VERSION 19.01.29}{
  \subsection{USER-VISIBLE CHANGES}{
    \itemize{
      \item Increased maximum number of covariates in likelihood models
      \code{nmix} and \code{nmixnb} from 10 to 15
    }
  }
}

\section{CHANGES IN VERSION 19.01.24}{
  \subsection{BUG FIXES}{
    \itemize{
      \item Added a new test-script
    }
  }
}

\section{CHANGES IN VERSION 18.12.12}{
  \subsection{NEW EXPERIMENTAL FEATURES}{
    \itemize{
      \item New models, \code{loggamma} and \code{mloggamma}
      in \code{mix}.
    }
  }
  \subsection{BUG FIXES}{
    \itemize{
      \item Minor changes in some build scripts.
    }
  }
}

\section{CHANGES IN VERSION 18.12.01}{
  \subsection{NEW EXPERIMENTAL FEATURES}{
    \itemize{
      \item New option \code{mkl} in \code{inla.setOption()} to chose
      MKL-buildt binaries.
      \item Linux binaries now buildt with Ubuntu1804.
      \item MKL-versions are included for MacOSX, and Linux (both dynamic
      and static).
    }
  }
}

\section{CHANGES IN VERSION 18.11.28}{
  \subsection{NEW EXPERIMENTAL FEATURES}{
    \itemize{
      \item New latent model \code{intslope}
    }
  }
}

\section{CHANGES IN VERSION 18.11.22}{
  \subsection{BUG FIXES}{
    \itemize{
      \item Improved \code{control.mix} interface and code
    }
  }
}

\section{CHANGES IN VERSION 18.10.29}{
  \subsection{NEW EXPERIMENTAL FEATURES}{
    \itemize{
      \item Likelihood model \code{nbinomial2}
    }
  }
}

\section{CHANGES IN VERSION 18.10.28}{
  \subsection{NEW EXPERIMENTAL FEATURES}{
    \itemize{
      \item New function \code{inla.priors.used}
    }
  }
}

\section{CHANGES IN VERSION 18.10.17}{
  \subsection{NEW FEATURES}{
    \itemize{
      \item Export class \code{inla}
    }
  }
}

\section{CHANGES IN VERSION 18.10.16}{
  \subsection{NEW EXPERIMENTAL FEATURES}{
    \itemize{
      \item New latent model: \code{dmatern}
    }
  }
  \subsection{BUG FIXES}{
    \itemize{
      \item Improved the numerics for computing the scaling of the RW1 and RW2 models.
    }
  }
}


\section{CHANGES IN VERSION 18.10.09}{
  \subsection{USER-VISIBLE CHANGES}{
    \itemize{
      \item New option \code{control.inla=list(tolerance.step=)}, to
      control the RMS of the step-size for the inner optimization.
      \item Changed, slightly, the initial values for the exponent in
      the Weibull likelihood models, to a value close to zero instead of
      zero.
      \item New vignette about how to deal with multinomial data.
    }
  }
  \subsection{NEW EXPERIMENTAL FEATURES}{
    \itemize{
      \item Added option \code{verbose} to
      \code{inla.qsample()} and \code{inla.posterior.sample()}
    }
  }
}



\section{CHANGES IN VERSION 18.09.24}{
  \subsection{BUG FIXES}{
    \itemize{
      \item Performance improvement using the PARDISO library
    }
  }
}

\section{CHANGES IN VERSION 18.09.21}{
  \subsection{NEW EXPERIMENTAL FEATURES}{
    \itemize{
      \item Argument \code{selection} in \code{inla.posterior.sample}
      and \code{inla.qsample}.
    }
  }
}


\section{CHANGES IN VERSION 18.09.19}{
  \subsection{BUG FIXES}{
    \itemize{
      \item Fix for \code{num.threads} in \code{inla.qinv()}
    }
  }
}

\section{CHANGES IN VERSION 18.09.18}{
  \subsection{BUG FIXES}{
    \itemize{
      \item Allow better user control of sparse matrix library in
      \code{inla.qinv()}, \code{inla.qsample()} and
      \code{inla.posterior.sample()}
    }
  }
}

\section{CHANGES IN VERSION 18.09.14}{
  \subsection{USER-VISIBLE CHANGES}{
    \itemize{
      \item New example added to \code{inla.posterior.sample()}
      \item Slight changes in the default \code{print}, and
      \code{summary} for an \code{inla}-object
    }
  }
  \subsection{BUG FIXES}{
    \itemize{
      \item Fixed the issue when \code{lincomb.derived.only=FALSE} and
      then using \code{inla.posterior.sample()}
    }
  }
}

\section{CHANGES IN VERSION 18.08.26}{
  \subsection{USER-VISIBLE CHANGES}{
    \itemize{
      \item Added 32-bit builds for windows (upon request)
    }
  }
  \subsection{NEW EXPERIMENTAL FEATURES}{
    \itemize{
      \item Added function \code{inla.posterior.sample.eval()}
    }
  }
}

\section{CHANGES IN VERSION 18.08.09}{
  \subsection{USER-VISIBLE CHANGES}{
    \itemize{
      \item Added new function \code{inla.pardiso.check()}
      \item Added COPYRIGHTS file
    }
  }
  \subsection{NEW EXPERIMENTAL FEATURES}{
    \itemize{
      \item Separated the quantile link for the binomial response,
      into individual (\code{model="quantile"}) and population
	(\code{model="pquantile"})

      \item Added new strategy
      \code{control.inla=list(strategy="adaptive")} which use the
      \code{simplified.laplace} approximations for fixed effects and
      low-dimensional model components, and the \code{gaussian}
      approximation
      otherwise. The argument 
      \code{adaptive.max} in \code{control.inla}
      determines what is low-dimensional in this
      context (default 10).
      }
  }
  \subsection{BUG FIXES}{
    \itemize{
      \item Removed some code not used anymore
    }
  }
}


\section{CHANGES IN VERSION 18.07.27}{
  \subsection{USER-VISIBLE CHANGES}{
    \itemize{
      \item NEWS page created (see \code{news(package="INLA")})
      \item Added vignette about the conditional logit model (thanks to
      Stefani Muff)
      \item Fixed missprints in the documentation for model \code{ar1c}
      (Thanks to Virgilio Gomez Rubio)
      \item Fixed documentation about argument \code{blas.num.threads} in \code{inla()}
    }
  }
}


\section{CHANGES IN VERSION 18.07.12}{
  \subsection{USER-VISIBLE CHANGES}{
    \itemize{
      \item Package built with \code{R-3.5}, both stable and testing
    }
  }
}

\section{CHANGES IN VERSION 18.07.11}{
  \subsection{USER-VISIBLE CHANGES}{
    \itemize{
      \item Package built for \code{R-3.4}, both stable and testing. 
    }
  }
}

\section{CHANGES IN VERSION 18.07.12}{
  \subsection{NEW EXPERIMENTAL FEATURES}{
    \itemize{
      \item Likelihood model \code{pom} (proportional odds model)
    }
  }
}
