\name{plot.inla}
\alias{plot.inla}
\alias{plot.surv.inla}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Default INLA plotting }
\description{
  Takes a fitted \code{inla} or \code{surv.inla} object produced by
  \code{inla} or \code{surv.inla} and plots
}
\usage{
\method{plot}{inla}(x, 
             plot.fixed.effects = TRUE,
             plot.lincomb = TRUE,
             plot.random.effects = TRUE,
             plot.hyperparameters = TRUE,
             plot.predictor = TRUE,
             plot.q = TRUE,
             plot.cpo = TRUE,
             ...)

\method{plot}{surv.inla}(x, 
             plot.fixed.effects = TRUE,
             plot.lincomb = TRUE,
             plot.random.effects = TRUE,
             plot.hyperparameters = TRUE,
             plot.predictor = TRUE,
             plot.q = TRUE,
             ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{A fitted  \code{inla} object produced by \code{inla} }
  \item{plot.fixed.effects}{Boolean indicating if posterior marginals
    for the fixed effects in the model should be plotted }
  \item{plot.lincomb}{Boolean indicating if posterior marginals
    for the linear combinations should be plotted }
  \item{plot.random.effects}{Boolean indicating if posterior mean and quantiles
    for the random effects in the model should be plotted  }
  \item{plot.hyperparameters}{Boolean indicating if posterior marginals
    for the hyperparameters in the model should be plotted }
  \item{plot.predictor}{Boolean indicating if posterior mean and quantiles
    for the linear predictor in the model should be plotted }
  \item{plot.q}{Boolean indicating if precision matrix should be displayed}
  \item{plot.cpo}{Boolean indicating if CPO/PIT valuesshould be plotted}
  \item{...}{Additional arguments to plot()}
}

\author{Havard Rue \email{hrue@math.ntnu.no} }
\seealso{\code{\link{inla}}, \code{\link{inla.surv}}}
\keyword{plot}
