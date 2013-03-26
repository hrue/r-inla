
\name{control.link}
\alias{control.link}
\alias{inla.control.link}
\title{Control variables in control.link}
\description{Control variables in \code{control.link} for use in \code{inla}}
\usage{
inla(formula, ..., control.link=list( <arg>=<value>[, <arg>=<value>]))
}
\arguments{
\item{link}{  The name of the link function/model}
\item{order}{  The \code{order} of the link function, where the interpretation of \code{order} is model-dependent.}
\item{nq}{  Number of quadrature-points used to do the numerical integration}
\item{hyper}{  Definition of the hyperparameter(s) for the link model chosen}
\item{initial}{  (OBSOLETE!) The initial value(s) for the hyperparameter(s)}
\item{fixed}{  (OBSOLETE!) A boolean variable if hyperparmater(s) is/are fixed or random}
\item{prior}{  (OBSOLETE!) The name of the prior distribution(s) for the hyperparmater(s)}
\item{param}{  (OBSOLETE!) The parameters for the prior distribution(s) for the hyperparmater(s)}
}
\seealso{
\code{\link{control.lincomb}}, \code{\link{control.group}}, \code{\link{control.mix}}, \code{\link{control.link}}, \code{\link{control.expert}}, \code{\link{control.compute}}, \code{\link{control.family}}, \code{\link{control.data}}, \code{\link{control.fixed}}, \code{\link{control.inla}}, \code{\link{control.predictor}}, \code{\link{control.results}}, \code{\link{control.mode}}, \code{\link{control.hazard}}, 
\code{\link{inla}}
}
