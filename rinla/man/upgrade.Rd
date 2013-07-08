\name{inla.upgrade}
\alias{inla.upgrade}

\title{Upgrade the INLA-package}

\description{This function will upgrade the \code{INLA}-package to the
  current version.}


\usage{
inla.upgrade(lib = NULL, testing=FALSE, force=FALSE)
}

\arguments{

  \item{lib}{Location to install the library.}

  \item{testing}{If \code{TRUE}, then look for a test-version if the
    INLA-package.}

  \item{force}{If \code{TRUE}, then force a new install}
}

\value{\code{inla.upgrade} will download and install the current version
  of \code{INLA}.}

\author{Havard Rue \email{hrue@math.ntnu.no}}

\seealso{\code{install.packages}}

