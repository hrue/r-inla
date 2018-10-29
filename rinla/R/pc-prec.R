## Export: inla.pc.rprec inla.pc.dprec inla.pc.qprec inla.pc.pprec

##! \name{pc.prec}
##! \alias{inla.pc.prec}
##! \alias{pc.prec}
##! \alias{pc.rprec}
##! \alias{inla.pc.rprec}
##! \alias{pc.dprec}
##! \alias{inla.pc.dprec}
##! \alias{pc.pprec}
##! \alias{inla.pc.pprec}
##! \alias{pc.qprec}
##! \alias{inla.pc.qprec}
##! 
##! \title{Utility functions for the PC prior for the precision}
##! 
##! \description{Functions to evaluate, sample, compute quantiles and
##!              percentiles of the PC prior for the precision
##!              in the Gaussian distribution.}
##! \usage{
##! inla.pc.rprec(n, u, alpha, lambda)
##! inla.pc.dprec(prec, u, alpha, lambda, log = FALSE)
##! inla.pc.qprec(p, u, alpha, lambda)
##! inla.pc.pprec(q, u, alpha, lambda)
##! }
##! \arguments{
##!   \item{n}{Number of observations}
##!   \item{u}{The upper limit (see Details)}
##!   \item{alpha}{The probability going above the upper limit (see Details)}
##!   \item{lambda}{The rate parameter (see Details)}
##!   \item{prec}{Vector of precisions}
##!   \item{log}{Logical. Return the density in natural or log-scale.}
##!   \item{p}{Vector of probabilities}
##!   \item{q}{Vector of quantiles}
##! }
##! \details{
##!     The statement \code{Prob(1/sqrt(prec) > u) = alpha} is used to
##!     determine \code{lambda} unless \code{lambda} is given.
##!     Either \code{lambda}  must be given,  or
##!     \code{u} AND \code{alpha}.
##! }
##!\value{%%
##!  \code{inla.pc.dprec} gives the density,
##!  \code{inla.pc.pprec} gives the distribution function,
##!  \code{inla.pc.qprec} gives the quantile function, and
##!  \code{inla.pc.rprec} generates random deviates.
##! }
##! \seealso{inla.doc("pc.prec")}
##! \author{Havard Rue \email{hrue@r-inla.org}}
##! \examples{
##! prec = inla.pc.rprec(100,  lambda = 1)
##! d = inla.pc.dprec(prec, lambda = 1)
##! prec = inla.pc.qprec(0.5, u = 1, alpha=0.01)
##! inla.pc.pprec(prec, u = 1, alpha=0.01)
##! }

inla.pc.prec.lambda = function(u, alpha, lambda)
{
    if (missing(lambda)) {
        stopifnot(!missing(u) && !missing(alpha))
        lambda = -log(alpha)/u
    }
    return (lambda)
}

inla.pc.rprec = function(n, u, alpha, lambda)
{
    lambda = inla.pc.prec.lambda(u, alpha, lambda)
    s = rexp(n, rate = lambda)
    prec = 1/s^2
    return(prec)
}

inla.pc.dprec = function(prec, u, alpha, lambda, log = FALSE)
{
    lambda = inla.pc.prec.lambda(u, alpha, lambda)
    s = 1/sqrt(prec)
    d = dexp(s, rate = lambda, log=log)
    log.jac = -log(2) - 3/2 * log(prec)
    if (log) {
        return (d + log.jac)
    } else {
        return (d*exp(log.jac))
    }
}

inla.pc.qprec = function(p, u, alpha, lambda)
{
    lambda = inla.pc.prec.lambda(u, alpha, lambda)
    q = qexp(1-p, rate = lambda)
    qprec = 1/q^2
    return (qprec)
}

inla.pc.pprec = function(q, u, alpha, lambda)
{
    lambda = inla.pc.prec.lambda(u, alpha, lambda)
    p = pexp(1/sqrt(q), rate = lambda)
    return (1-p)
}
