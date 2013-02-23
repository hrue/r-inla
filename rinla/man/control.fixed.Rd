
\name{control.fixed}
\alias{control.fixed}
\alias{inla.control.fixed}
\title{Control variables in control.fixed}
\description{Control variables in \code{control.fixed} for use in \code{inla}}
\usage{
inla(formula, ..., control.fixed=list( <arg>=<value>[, <arg>=<value>]))
}
\arguments{
\item{cdf}{  A list of values to compute the CDF for, for all fixed effects}
\item{quantiles}{  A list of quantiles to compute for all fixed effects}
\item{mean}{  Prior mean for all fixed effects except the intercept. Alternatively, a named list with specific means where name=default applies to unmatched names. For example \code{control.fixed=list(mean=list(a=1, b=2, default=0))} assign 'mean=1' to fixed effect 'a' , 'mean=2' to effect 'b' and 'mean=0' to all others.}
\item{mean.intercept}{  Prior mean for the intercept}
\item{prec}{  Default precision for all fixed effects except the intercept. Alternatively, a named list with specific means where name=default applies to unmatched names. For example \code{control.fixed=list(prec=list(a=1, b=2, default=0.01))} assign 'prec=1' to fixed effect 'a' , 'prec=2' to effect 'b' and 'prec=0.01' to all others.}
\item{prec.intercept}{  Default precision the intercept (default 0.0)}
\item{compute}{  Compute marginals for the fixed effects ? (default TRUE)}
\item{correlation.matrix}{  Compute the posterior correlation matrix for all fixed effects? (default FALSE) OOPS: This option will set up appropriate linear combinations and the results are shown as the posterior correlation matrix of the linear combinations. This option will imply \code{control.inla=list(lincomb.derived.correlation.matrix=TRUE)}.}
}
\seealso{
\code{\link{control.lincomb}}, \code{\link{control.group}}, \code{\link{control.mix}}, \code{\link{control.link}}, \code{\link{control.expert}}, \code{\link{control.compute}}, \code{\link{control.family}}, \code{\link{control.data}}, \code{\link{control.fixed}}, \code{\link{control.inla}}, \code{\link{control.predictor}}, \code{\link{control.results}}, \code{\link{control.mode}}, \code{\link{control.hazard}}, 
\code{\link{inla}}
}
