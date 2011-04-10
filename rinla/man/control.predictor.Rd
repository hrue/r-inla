
\name{control.predictor}
\alias{control.predictor}
\alias{inla.control.predictor}
\title{Control variables in control.predictor}
\description{Control variables in \code{control.predictor} for use in \code{inla}}
\usage{
inla(formula, ..., control.predictor=list( <arg>=<value>[, <arg>=<value>]))
}
\arguments{
\item{hyper}{  Definition of the hyperparameters.}
\item{fixed}{  (OBSOLETE!) If the precision for the artificial noise is fixed or not (defualt TRUE)}
\item{prior}{  (OBSOLETE!) The prior for the artificial noise}
\item{param}{  (OBSOLETE!) Prior parameters for the artificial noise}
\item{initial}{  (OBSOLETE!) The value of the log precision of the artificial noise}
\item{compute}{  A boolean variable; should the marginals for the linear predictor be computed? (Default FALSE.)}
\item{cdf}{  A list of values to compute the CDF for the linear predictor}
\item{quantiles}{  A list of quantiles to compute for the linear predictor}
\item{cross}{  Cross-sum-to-zero constraints}
\item{A}{  The observation matrix A = list(i=..., j=..., values=...) or a filename.}
}
\seealso{
\code{\link{control.lincomb}}, \code{\link{control.group}}, \code{\link{control.expert}}, \code{\link{control.compute}}, \code{\link{control.data}}, \code{\link{control.fixed}}, \code{\link{control.inla}}, \code{\link{control.predictor}}, \code{\link{control.results}}, \code{\link{control.mode}}, \code{\link{control.hazard}}, 
\code{\link{inla}}
}
