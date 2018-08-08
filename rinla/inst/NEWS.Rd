\newcommand{\PR}{\Sexpr[results=rd]{tools:::Rd_expr_PR(#1)}}

\name{NEWS}
\title{R-INLA News}
\encoding{UTF-8}


\section{CHANGES IN VERSION XX.XX.XX}{
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
      \item Added new strategy
      \code{control.inla=list(strategy="adaptive")} which use the
      simplified.laplace approximations for fixed effects and low-dimensional model
      components, and the gaussian approximation otherwise. 
      }
  }
  \subsection{BUG FIXES}{
    \itemize{
      \item 
    }
  }
}

\section{CHANGES IN VERSION XX.XX.XX}{
  \subsection{USER-VISIBLE CHANGES}{
    \itemize{
      \item Added new function \code{inla.pardiso.check()}
      \item Added COPYRIGHTS file
    }
  }
  \subsection{NEW FEATURES}{
    \itemize{
      \item
    }
  }
  \subsection{NEW EXPERIMENTAL FEATURES}{
    \itemize{
      \item Separated the quantile link for the binomial response,
      into individual (\texttt{model="quantile"}) and population
	(\texttt{model="pquantile"})
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
