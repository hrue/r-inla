## Export: inla.pc.rgammacount inla.pc.dgammacount inla.pc.qgammacount inla.pc.pgammacount

##! \name{pc.gammacount}
##! \alias{inla.pc.gammacount}
##! \alias{pc.gammacount}
##! \alias{pc.rgammacount}
##! \alias{inla.pc.rgammacount}
##! \alias{pc.dgammacount}
##! \alias{inla.pc.dgammacount}
##! \alias{pc.pgammacount}
##! \alias{inla.pc.pgammacount}
##! \alias{pc.qgammacount}
##! \alias{inla.pc.qgammacount}
##! 
##! \title{Utility functions for the PC prior for the \code{gammacount} likelihood}
##! 
##! \description{Functions to evaluate, sample, compute quantiles and
##!              percentiles of the PC prior for the \code{gammacount} likelihood}
##! \usage{
##! inla.pc.rgammacount(n, lambda = 1)
##! inla.pc.dgammacount(x, lambda = 1, log = FALSE)
##! inla.pc.qgammacount(p, lambda = 1)
##! inla.pc.pgammacount(q, lambda = 1)
##! }
##! \arguments{
##!   \item{n}{Number of observations}
##!   \item{lambda}{The rate parameter (see Details)}
##!   \item{x}{Evaluation points}
##!   \item{log}{Logical. Return the density in natural or log-scale.}
##!   \item{p}{Vector of probabilities}
##!   \item{q}{Vector of quantiles}
##! }
##! \details{
##! This gives the PC prior for the \code{gammacount} likelihood,  which is the PC prior for
##! \code{a} in \code{Gamma(a, 1)} where \code{Gamma(1, 1)} is the base model.
##! }
##!\value{%%
##!  \code{inla.pc.dgammacount} gives the density,
##!  \code{inla.pc.pgammacount} gives the distribution function,
##!  \code{inla.pc.qgammacount} gives the quantile function, and
##!  \code{inla.pc.rgammacount} generates random deviates.
##! }
##! \seealso{inla.doc("pc.gammacount")}
##! \author{Havard Rue \email{hrue@r-inla.org}}
##! \examples{
##! x = inla.pc.rgammacount(100,  lambda = 1)
##! d = inla.pc.dgammacount(x, lambda = 1)
##! x = inla.pc.qgammacount(0.5, lambda = 1)
##! inla.pc.pgammacount(x, lambda = 1)
##! }

## I think this could have be done better than using the generic routines, by tabulating the d()
## function.

inla.pc.gammacount.intern = function(lambda = 1)
{
    ## return a marginal-object in the log(x)-scale
    f = if (lambda >= 1.0) 1.0 else 0.25
    log.x = (seq(-10, 10, len=4096))/lambda^f
    idx = which(log.x == 0.0)
    if (length(idx) > 0) {
        log.x = log.x[-idx]
    }
    x = exp(log.x)
    marg = list(x = log.x, y = inla.pc.dgammacount(x, lambda = lambda, log=FALSE) * x)
    return (marg)
}

inla.pc.rgammacount = function(n, lambda = 1)
{
    log.x = inla.rmarginal(n, inla.pc.gammacount.intern(lambda = lambda))
    return (exp(log.x))
}

inla.pc.dgammacount = function(x, lambda = 1, log = FALSE)
{
    ## > d := a -> (-2*ln(GAMMA(a)) + 2*(a - 1)*Psi(a))^(1/2);
    ## > R(log(lambda) -lambda*d(x) + log(abs(diff(d(x),x))), optimize);

    t1 <- log(lambda)
    t3 <- lgamma(x)
    t4 <- x - 1
    t5 <- digamma(x)
    t8 <- sqrt(2) * sqrt(t5 * t4 - t3)
    t12 <- trigamma(x)
    t14 <- abs(t12 * t4 / t8)
    t15 <- log(t14)
    t16 <- -t8 * lambda + t1 + t15
    ld <- t16 - log(2.0)
    
    return (if (log) ld else exp(ld))
}

inla.pc.qgammacount = function(p, lambda = 1)
{
    log.x = inla.pmarginal(p, inla.pc.gammacount.intern(lambda = lambda))
    return (exp(log.x))
}

inla.pc.pgammacount = function(q, lambda = 1)
{
    p = inla.qmarginal(log(q), inla.pc.gammacount.intern(lambda = lambda))
    return (p)
}
