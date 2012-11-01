
\name{control.compute}
\alias{control.compute}
\alias{inla.control.compute}
\title{Control variables in control.compute}
\description{Control variables in \code{control.compute} for use in \code{inla}}
\usage{
inla(formula, ..., control.compute=list( <arg>=<value>[, <arg>=<value>]))
}
\arguments{
\item{strategy}{  The computational strategy to use: 'small', 'medium', 'large', 'huge' and 'default'. The difference is how the parallelisation is done, and is tuned for 'small'-sized models, 'medium'-sized models, etc. The default option tries to make an educated guess, but this allows to overide this selection. Default is 'default'}
\item{hyperpar}{  A boolean variable if the marginal for the hyperparameters should be computed. Default TRUE.}
\item{return.marginals}{  A boolean variable if the marginals for the latent field should be returned (although it is computed). Default TRUE}
\item{dic}{  A boolean variable if the DIC-value should be computed. Default FALSE.}
\item{mlik}{  A boolean variable if the marginal likelihood should be computed. Default FALSE.}
\item{cpo}{  A boolean variable if the cross-validated predictive measures (cpo, pit) should be computed}
\item{q}{  A boolean variable if binary images of the precision matrix, the reordered precision matrix and the Cholesky triangle should be generated. (Default FALSE.)}
\item{config}{  A boolean variable if the internal GMRF approximations be stored. (Default FALSE. EXPERIMENTAL)}
\item{smtp}{  The sparse-matrix solver, one of 'smtp' (default) or 'band'}
\item{graph}{  A boolean variable if the graph itself should be returned. (Default FALSE.)}
}
\seealso{
\code{\link{control.lincomb}}, \code{\link{control.group}}, \code{\link{control.expert}}, \code{\link{control.compute}}, \code{\link{control.family}}, \code{\link{control.data}}, \code{\link{control.fixed}}, \code{\link{control.inla}}, \code{\link{control.predictor}}, \code{\link{control.results}}, \code{\link{control.mode}}, \code{\link{control.hazard}}, 
\code{\link{inla}}
}
