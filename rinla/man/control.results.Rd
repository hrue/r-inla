
\name{control.results}
\alias{control.results}
\alias{inla.control.results}
\title{Control variables in control.results}
\description{Control variables in \code{control.results} for use in \code{inla}}
\usage{
inla(formula, ..., control.results=list( <arg>=<value>[, <arg>=<value>]))
}
\arguments{
\item{return.marginals.random}{  A boolean variable; read the marginals for the fterms? (Default TRUE)}
\item{return.marginals.predictor}{  A boolean variable; read the marginals for the linear predictor? (Default TRUE)}
\item{image.dim}{  The dimension to reduce the matrix-images returned by the argument 'control.compute(q=TRUE)' when displaying the images}
}
\seealso{
\code{\link{control.lincomb}}, \code{\link{control.group}}, \code{\link{control.expert}}, \code{\link{control.compute}}, \code{\link{control.data}}, \code{\link{control.fixed}}, \code{\link{control.inla}}, \code{\link{control.predictor}}, \code{\link{control.results}}, \code{\link{control.mode}}, \code{\link{control.hazard}}, 
\code{\link{inla}}
}
