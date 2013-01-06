
\name{control.group}
\alias{control.group}
\alias{inla.control.group}
\title{Control variables in control.group}
\description{Control variables in \code{control.group} for use in \code{inla}}
\usage{
inla(formula, ..., control.group=list( <arg>=<value>[, <arg>=<value>]))
}
\arguments{
\item{model}{  Group model (one of 'exchangable', 'ar1', 'ar', 'rw1' or 'rw2')}
\item{order}{  Defines the \code{order} of the model: for model \code{ar} this defines the order p, in AR(p). Not used for other models at the time being.}
\item{cyclic}{  Make the group model cyclic? (Only applies to models 'ar1', 'rw1' and 'rw2')}
\item{hyper}{  Definition of the hyperparameter(s)}
\item{initial}{  (OBSOLETE!) The initial value for the group correlation or precision in the internal scale.}
\item{fixed}{  (OBSOLETE!) A boolean variable if the group correction or precision is assumed to be fixed or random.}
\item{prior}{  (OBSOLETE!) The name of the prior distribution for the group correlation or precision in the internal scale}
\item{param}{  (OBSOLETE!) Prior parameters}
}
\seealso{
\code{\link{control.lincomb}}, \code{\link{control.group}}, \code{\link{control.mix}}, \code{\link{control.expert}}, \code{\link{control.compute}}, \code{\link{control.family}}, \code{\link{control.data}}, \code{\link{control.fixed}}, \code{\link{control.inla}}, \code{\link{control.predictor}}, \code{\link{control.results}}, \code{\link{control.mode}}, \code{\link{control.hazard}}, 
\code{\link{inla}}
}
