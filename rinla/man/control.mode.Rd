
\name{control.mode}
\alias{control.mode}
\alias{inla.control.mode}
\alias{inla.set.control.mode.default}
\alias{set.control.mode.default}
\alias{control.mode.default}
\title{Control variables in control.mode}
\description{Control variables in \code{control.mode} for use in \code{inla}}
\usage{
inla.set.control.mode.default(...)
}
\arguments{
\item{...}{Additional arguments; none for the moment}
}
\value{
The function \code{inla.set.control.mode.default} returns a list with all the default values of all parameters within this control statement. These are as follows.
\item{result}{  Prevous result from inla(). Use the theta- and x-mode from this run.}
\item{theta}{  The theta-mode/initial values for theta. This option has preference over result$mode$theta.}
\item{x}{  The x-mode/intitial values for x. This option has preference over result$mode$x.}
\item{restart}{  A boolean variable; should we restart the optimisation from this configuration or fix the mode at this configuration? (Default FALSE.)}
\item{fixed}{  A boolean variable. If TRUE then treat all thetas as known and fixed, and if FALSE then treat all thetas as unknown and random (default).}
}
\seealso{
\code{\link{control.lincomb}}, \code{\link{control.group}}, \code{\link{control.mix}}, \code{\link{control.link}}, \code{\link{control.expert}}, \code{\link{control.compute}}, \code{\link{control.family}}, \code{\link{control.data}}, \code{\link{control.fixed}}, \code{\link{control.inla}}, \code{\link{control.predictor}}, \code{\link{control.results}}, \code{\link{control.mode}}, \code{\link{control.hazard}}, 
\code{\link{inla}}
}
