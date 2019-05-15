\newcommand{\PR}{\Sexpr[results=rd]{tools:::Rd_expr_PR(#1)}}

\name{NEWS}
\title{R-INLA News}
\encoding{UTF-8}



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
