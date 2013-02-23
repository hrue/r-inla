
\name{control.hazard}
\alias{control.hazard}
\alias{inla.control.hazard}
\title{Control variables in control.hazard}
\description{Control variables in \code{control.hazard} for use in \code{inla}}
\usage{
inla(formula, ..., control.hazard=list( <arg>=<value>[, <arg>=<value>]))
}
\arguments{
\item{model}{  The model for the baseline hazard model. One of 'rw1' or 'rw2'. (Default 'rw1'.)}
\item{hyper}{  The definition of the hyperparameters.}
\item{fixed}{  (OBSOLETE!) A boolean variable; is the precision for 'model' fixed? (Default FALSE.)}
\item{initial}{  (OBSOLETE!) The initial value for the precision.}
\item{prior}{  (OBSOLETE!) The prior distribution for the precision for 'model'}
\item{param}{  (OBSOLETE!) The parameters in the prior distribution}
\item{constr}{  A boolean variable; shall the 'model' be constrained to sum to zero?}
\item{n.intervals}{  Number of intervals in the baseline hazard. (Default 15)}
\item{cutpoints}{  The cutpoints to use. If not specified the they are compute from 'n.intervals' and the maximum length of the interval. (Default NULL)}
\item{strata.name}{  The name of the stratefication variable for the baseline hazard in the data.frame}
\item{si}{  A boolean variable; should all the Gaussian approximations configurations be written to files?}
}
\seealso{
\code{\link{control.lincomb}}, \code{\link{control.group}}, \code{\link{control.mix}}, \code{\link{control.link}}, \code{\link{control.expert}}, \code{\link{control.compute}}, \code{\link{control.family}}, \code{\link{control.data}}, \code{\link{control.fixed}}, \code{\link{control.inla}}, \code{\link{control.predictor}}, \code{\link{control.results}}, \code{\link{control.mode}}, \code{\link{control.hazard}}, 
\code{\link{inla}}
}
