
\name{control.family}
\alias{control.family}
\alias{inla.control.family}
\title{Control variables in control.family}
\description{Control variables in \code{control.family} for use in \code{inla}}
\usage{
inla(formula, ..., control.family=list( <arg>=<value>[, <arg>=<value>]))
}
\arguments{
\item{hyper}{  Definition of the hyperparameters}
\item{initial}{  (OBSOLETE!) Initial value for the hyperparameter(s) of the likelihood in the internal scale.}
\item{prior}{  (OBSOLETE!) The name of the prior distribution(s) for othe hyperparameter(s).}
\item{param}{  (OBSOLETE!) The parameters for the prior distribution}
\item{fixed}{  (OBSOLETE!) Boolean variable(s) to say if the hyperparameter(s) is fixed or random.}
\item{link}{  The link function to use.}
\item{alpha}{  The parameter 'alpha' for the asymmetric Laplace likelihood (default 0.5)}
\item{epsilon}{  The parameter 'epsilon' for the asymmetric Laplace likelihood (default 0.01)}
\item{gamma}{  The parameter 'gamma' for the asymmetric Laplace likelihood (default 1.0)}
\item{sh.shape.max}{  Maximum value for the shape-parameter for Skew Normal observations}
\item{gev.scale.xi}{  The scaling of the shape-parameter for the GEV distribution. (default 0.01)}
\item{variant}{  This variable is used to give options for various variants of the likelihood, like chosing different parameterisations for example. See the relevant likelihood documentations for options (does only apply to some likelihoods).}
\item{control.mix}{  See \code{?control.mix}}
}
\seealso{
\code{\link{control.lincomb}}, \code{\link{control.group}}, \code{\link{control.mix}}, \code{\link{control.expert}}, \code{\link{control.compute}}, \code{\link{control.family}}, \code{\link{control.data}}, \code{\link{control.fixed}}, \code{\link{control.inla}}, \code{\link{control.predictor}}, \code{\link{control.results}}, \code{\link{control.mode}}, \code{\link{control.hazard}}, 
\code{\link{inla}}
}
