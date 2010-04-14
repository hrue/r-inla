\name{INLA-package}
\alias{INLA-package}
\alias{INLA}
\docType{package}
\title{
Integrated Nested Laplace Approximation 
}

\description{ Package to perform full Bayesian analysis on generalised
additive mixed models using Integrated Nested Laplace Approximations.
}


\details{
\tabular{ll}{
Package: \tab INLA\cr
Type: \tab Package\cr
Version: \tab 0.1\cr
Date: \tab 2009-01-01 \cr
License: \tab GPL\cr
NOTE: \tab This package has no version number yet; it is to heavily developed
at the moment...
}

The INLA package is a R interface to the GMRFLib library and the
\emph{inla} program. These are written in (ANSI) C and both need to be
downloaded and installed for the INLA package to work.  The \emph{inla}
program can be downloaded from
\url{http://www.math.ntnu.no/~hrue/GMRFLib/} in a binary format.

The INLA package provides with full Bayesian analysis for generalized
additive models using Integrated Nested Laplace Approximation. For
details about Integrated Nested Laplace Approximation see \cite{Rue et
al. (2008)}, for a wide series of worked out examples solved using INLA
see \cite{ Martino and Rue (2008)}.%%
} 
%%

\author{Sara Martino and Havard Rue.  Maintainer: Havard Rue
<hrue@math.ntnu.no> (On leave 2008: Sara Martino <martino@math.ntnu.no>)
}

\references{Rue, H. and Martino, S. and Chopin, N. (2009)
  \emph{Approximate Bayesian Inference for latent Gaussian models using
  Integrated Nested Laplace Approximations}, Journal of Royal
  Statistical Society, Series B, with discussion. (To appear.)

  Rue, H and Held, L. (2005) \emph{Gaussian Markov Random Fields -
  Theory and Applications}, Monographs on Statistics and Applied
  Probability, no 104, Chapman \& Hall, London.
  
  Martino, S. and Rue, H. (2008)
  \emph{Implementing Approximate Bayesian Inference using Integrated
    Nested Laplace Approximation: a manual for the inla program}
  Preprint,  No.2, Dep. of Mathematical Sciences, NTNU,  Norway
}

\keyword{Approximate Bayesian inference, Laplace approximations,
  Integrated nested Laplace approximations}
