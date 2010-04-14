\name{interpret.formula}
\alias{interpret.formula}

\title{ Interpret a INLA formula }
\description{
  This is an internal function of package \code{INLA}.
  It is a service routine for gam which splits off the strictly
  parametric part of the model formula,
  returning it as a formula, and interprets the different parts of the model formula.

Not normally called directly.

}
\usage{
interpret.formula(gf)
}

\arguments{
  \item{gf}{ An formula supplied to \code{inla} }
}
\author{ Sara Martino \email{martino@math.ntnu.no} }
\seealso{ \code{\link{inla.models}},  \code{\link{inla}}}

\keyword{ internal }

