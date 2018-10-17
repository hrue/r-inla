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
