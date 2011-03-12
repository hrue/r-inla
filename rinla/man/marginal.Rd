\name{marginal}
\alias{marginal}
\alias{pmarginal}
\alias{inla.pmarginal}
\alias{qmarginal}
\alias{inla.qmarginal}
\alias{dmarginal}
\alias{inla.dmarginal}
\alias{rmarginal}
\alias{inla.rmarginal}
\alias{inla.expectation}
\alias{inla.marginal.expectation}
\alias{marginal.expectation}
\alias{inla.spline}
\alias{inla.marginal}
\alias{inla.marginal.transform}
\alias{marginal.transform}


\title{Functions which operates on marginals}

\description{Density, distribution function, quantile function, random
     generation, interpolation and computation of expectations for
     marginals obtained by \code{inla} or \code{inla.hyperpar()}.}

\usage{
inla.dmarginal = function(x, marginal, log = FALSE)
inla.pmarginal = function(q, marginal, normalize = TRUE)
inla.qmarginal = function(p, marginal, len = 1024)
inla.rmarginal = function(n, marginal)
inla.marginal.expectation = function(fun, marginal, ...)
inla.marginal.transform = function(fun, marginal, n, h.diff, ...)
inla.spline = function(marginal, log = FALSE, extrapolate = 0.0)
}

\arguments{

  \item{marginal}{A marginal object from either \code{inla} or
    \code{inla.hyperpar()}, which is either \code{list(x=c(), y=c())}
    with density values \code{y} at locations \code{x}, or a
    \code{matrix(,n,2)} for which the density values are the second
    column and the locations in the first column.}

  \item{fun}{A (vectorised) function like \code{function(x) exp(x)} to
    compute the expectation against, or which define the transformation
    new = fun(old)}
  
  \item{x}{Evaluation points}

  \item{q}{Quantiles}

  \item{p}{Probabilities}

  \item{n}{The number of observations. If \code{length(n) > 1}, the
    length is taken to be the number required. For
    \code{inla.marginal.transform}, its the number of points to use
    in the new density.}

  \item{h.diff}{The step-length for the numerical differeniation inside \code{inla.marginal.transform}}
    
  \item{...}{Further arguments to be passed to function which
    expectation is to be computed.}

  \item{log}{Return density or interpolated density in log-scale?}

  \item{normalize}{Renormalise the density after interpolation?}
  \item{len}{Number of locations used to interpolate the distribution
  function.}

}

\value{%%
  \code{inla.spline} returns \code{list=c(x=c(), y=c())} of
  interpolated values do extrapolation using the factor given, whereas
  the remaining function returns what they say they should do.  }
%%

\author{Havard Rue \email{hrue@math.ntnu.no}}

\seealso{\code{\link{inla}}, \code{\link{inla.hyperpar}}}

\examples{
## a simple linear regression example
n = 10
x = rnorm(n)
sd = 0.1
y = 1+x + rnorm(n,sd=sd)
res = inla(y ~ 1 + x, data = data.frame(x,y),
           control.data=list(initial = log(1/sd^2),fixed=TRUE))

## chose a marginal and compare the with the results computed by the
## inla-program
r = res$summary.fixed["x",]
m = res$marginals.fixed$x

## compute the the density for exp(r), version 1
r.exp = inla.marginal.transform(exp, m)
## or version 2
r.exp = inla.marginal.transform(function(x) exp(x), m)

## to plot the marginal, we use the inla.spline, which interpolates (in
## log-scale). Compare with some samples.
plot(inla.spline(m), type="l")
s = inla.rmarginal(1000, m)
hist(inla.rmarginal(1000, m), add=TRUE, prob=TRUE)
lines(density(s), lty=2)

m1 = inla.marginal.expectation(function(x) x^1, m)
m2 = inla.marginal.expectation(function(x) x^2, m)
stdev = sqrt(m2 - m1^2)
q = inla.qmarginal(c(0.025,0.975), m)

## inla-program results
print(r)

## inla.marginal-results (they shouldn't be perfect!)
print(c(mean=m1, sd=stdev, "0.025quant" = q[1], "0.975quant" = q[2]))
}

