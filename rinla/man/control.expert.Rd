
\name{control.expert}
\alias{control.expert}
\alias{inla.control.expert}
\title{Control variables in control.expert}
\description{Control variables in \code{control.expert} for use in \code{inla}}
\usage{
inla(formula, ..., control.expert=list( <arg>=<value>[, <arg>=<value>]))
inla.set.control.expert.default()
}
\arguments{
The function \code{inla.set.control.expert.default} returns the default value of all parameters within this control group. These are as follows.
\item{cpo.manual}{  A boolean variable to decide if the inla-program is to be runned in a manual-cpo-mode. (EXPERT OPTION.)}
\item{cpo.idx}{  The index of the data point to remove. (EXPERT OPTION.)}
}
\seealso{
\code{\link{control.lincomb}}, \code{\link{control.group}}, \code{\link{control.mix}}, \code{\link{control.link}}, \code{\link{control.expert}}, \code{\link{control.compute}}, \code{\link{control.family}}, \code{\link{control.data}}, \code{\link{control.fixed}}, \code{\link{control.inla}}, \code{\link{control.predictor}}, \code{\link{control.results}}, \code{\link{control.mode}}, \code{\link{control.hazard}}, 
\code{\link{inla}}
}
