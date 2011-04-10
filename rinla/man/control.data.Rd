
\name{control.data}
\alias{control.data}
\alias{inla.control.data}
\title{Control variables in control.data}
\description{Control variables in \code{control.data} for use in \code{inla}}
\usage{
inla(formula, ..., control.data=list( <arg>=<value>[, <arg>=<value>]))
}
\arguments{
\item{hyper}{  Definition of the hyperparameters}
\item{initial}{  (OBSOLETE!) Initial value for the hyperparameter(s) of the likelihood in the internal scale.}
\item{prior}{  (OBSOLETE!) The name of the prior distribution(s) for othe hyperparameter(s).}
\item{param}{  (OBSOLETE!) The parameters for the prior distribution}
\item{fixed}{  (OBSOLETE!) Boolean variable(s) to say if the hyperparameter(s) is fixed or random.}
\item{link}{  The link function to use.}
\item{dof.max}{  Maximum degrees of freedom for the student-t distribution before its is considered as Normal (default 100)}
\item{alpha}{  The parameter 'alpha' for the asymmetric Laplace likelihood (default 0.5)}
\item{epsilon}{  The parameter 'epsilon' for the asymmetric Laplace likelihood (default 0.01)}
\item{gamma}{  The parameter 'gamma' for the asymmetric Laplace likelihood (default 1.0)}
\item{sh.shape.max}{  Maximum value for the shape-parameter for Skew Normal observations}
\item{gev.scale.xi}{  The scaling of the shape-parameter for the GEV distribution. (default 0.01)}
}
\seealso{
\code{\link{control.lincomb}}, \code{\link{control.group}}, \code{\link{control.expert}}, \code{\link{control.compute}}, \code{\link{control.data}}, \code{\link{control.fixed}}, \code{\link{control.inla}}, \code{\link{control.predictor}}, \code{\link{control.results}}, \code{\link{control.mode}}, \code{\link{control.hazard}}, 
\code{\link{inla}}
}
