
\name{control.mix}
\alias{control.mix}
\alias{inla.control.mix}
\title{Control variables in control.mix}
\description{Control variables in \code{control.mix} for use in \code{inla}}
\usage{
inla(formula, ..., control.mix=list( <arg>=<value>[, <arg>=<value>]))
}
\arguments{
\item{use}{  Integrate the likelihood with a random-effect added to the linear predictor? Default FALSE}
\item{model}{  The model for the random effect. Currently, only \code{model='gaussian'} is implemented}
\item{hyper}{  Definition of the hyperparameter(s) for the random effect model chosen}
\item{initial}{  (OBSOLETE!) The initial value(s) for the hyperparameter(s)}
\item{fixed}{  (OBSOLETE!) A boolean variable if hyperparmater(s) is/are fixed or random}
\item{prior}{  (OBSOLETE!) The name of the prior distribution(s) for the hyperparmater(s)}
\item{param}{  (OBSOLETE!) The parameters for the prior distribution(s) for the hyperparmater(s)}
}
\seealso{
\code{\link{control.lincomb}}, \code{\link{control.group}}, \code{\link{control.mix}}, \code{\link{control.expert}}, \code{\link{control.compute}}, \code{\link{control.family}}, \code{\link{control.data}}, \code{\link{control.fixed}}, \code{\link{control.inla}}, \code{\link{control.predictor}}, \code{\link{control.results}}, \code{\link{control.mode}}, \code{\link{control.hazard}}, 
\code{\link{inla}}
}
