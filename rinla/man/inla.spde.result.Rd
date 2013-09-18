\name{inla.spde.result}
\alias{inla.spde.result}
\alias{inla.spde.result.inla.spde1}
\alias{inla.spde.result.inla.spde2}
\alias{inla.spde1.result}
\alias{inla.spde2.result}
\title{
  SPDE result extraction from INLA estimation results
}
\description{
  Exctract field and parameter values and distributions for an
  \code{inla.spde} SPDE effect from an inla result
  object.
}
\usage{
inla.spde.result(...)

\method{inla.spde.result}{inla.spde1}(inla, name, spde, do.transform = TRUE, ...)
\method{inla.spde.result}{inla.spde2}(inla, name, spde, do.transform = TRUE, ...)

## Direct function call for class 'inla.spde1':
inla.spde1.result(inla, name, spde, do.transform = TRUE, ...)
## Direct function call for class 'inla.spde2':
inla.spde2.result(inla, name, spde, do.transform = TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{inla}{
    An \code{inla} object obtained from a call to \code{\link{inla}}
  }
  \item{name}{
    A character string with the name of the SPDE effect in the inla formula.
  }
  \item{spde}{
    The \code{inla.spde} object used for the effect in the inla
    formula. (Note: this could have been stored in the inla output, but
    isn't.)
    Usually the result of a call to \code{\link{inla.spde2.matern}}.
  }
  \item{do.transform}{
    If \code{TRUE}, also calculate marginals transformed to user-scale.
    Setting to \code{FALSE} is useful for large non-stationary models,
    as transforming many marginal densities is time-consuming.
  }
  \item{\dots}{
    Further arguments passed to and from other methods.
  }
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
  For \code{inla.spde2} models, a list, where the nominal range and
  variance are defined as the values that would have been obtained with
  a stationary model and no boundary effects:
  \item{marginals.kappa }{Marginal densities for kappa}
  \item{marginals.log.kappa }{Marginal densities for log(kappa)}
  \item{marginals.log.range.nominal }{Marginal densities for log(range)}
  \item{marginals.log.tau }{Marginal densities for log(tau)}
  \item{marginals.log.variance.nominal }{Marginal densities for log(variance)}
  \item{marginals.range.nominal }{Marginal densities for range}
  \item{marginals.tau }{Marginal densities for tau}
  \item{marginals.theta }{Marginal densities for the theta parameters}
  \item{marginals.values }{Marginal densities for the field values}
  \item{marginals.variance.nominal }{Marginal densities for variance}
  \item{summary.hyperpar }{The SPDE related part of the inla hyperpar
    output summary}
  \item{summary.log.kappa }{Summary statistics for log(kappa)}
  \item{summary.log.range.nominal }{Summary statistics for log(range)}
  \item{summary.log.tau }{Summary statistics for log(tau)}
  \item{summary.log.variance.nominal }{Summary statistics for log(kappa)}
  \item{summary.theta }{Summary statistics for the theta parameters}
  \item{summary.values }{Summary statistics for the field values}
}
%\references{
%% ~put references to the literature/web site here ~
%}
\author{
  Finn Lindgren \email{finn.lindgren@gmail.com}
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{inla.spde.models}},
\code{\link{inla.spde2.matern}}
}
\examples{
loc = matrix(runif(100*2),100,2)
mesh = inla.mesh.create.helper(points.domain=loc, max.edge=c(0.1,0.5))
spde = inla.spde2.matern(mesh)
index = inla.spde.make.index("spatial", mesh$n, n.repl=2)
spatial.A = inla.spde.make.A(mesh, loc,
                             index=rep(1:nrow(loc), 2),
                             repl=rep(1:2, each=nrow(loc)))
## Toy example with no spatial correlation (range=zero)
y = 10+rnorm(100*2)
stack = inla.stack(data=list(y=y),
                   A=list(spatial.A),
                   effects=list(c(index, list(intercept=1))),
                   tag="tag")
data = inla.stack.data(stack, spde=spde)
formula = y ~ -1 + intercept + f(spatial, model=spde,
                                 replicate=spatial.repl)
result = inla(formula, family="gaussian", data=data, 
              control.predictor=list(A=inla.stack.A(stack)))
spde.result = inla.spde.result(result, "spatial", spde)
plot(spde.result$marginals.range.nominal[[1]], type="l")
}
%% Note: the same example is used by inla.spde.make.index.Rd
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
