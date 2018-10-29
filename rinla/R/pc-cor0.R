## Export: inla.pc.rcor0 inla.pc.dcor0 inla.pc.qcor0 inla.pc.pcor0

##! \name{pc.cor0}
##! \alias{inla.pc.cor0}
##! \alias{pc.cor0}
##! \alias{pc.rcor0}
##! \alias{inla.pc.rcor0}
##! \alias{pc.dcor0}
##! \alias{inla.pc.dcor0}
##! \alias{pc.pcor0}
##! \alias{inla.pc.pcor0}
##! \alias{pc.qcor0}
##! \alias{inla.pc.qcor0}
##! 
##! \title{Utility functions for the PC prior for correlation in AR(1)}
##! 
##! \description{Functions to evaluate, sample, compute quantiles and
##!              percentiles of the PC prior for the correlation
##!              in the Gaussian AR(1) model where the base-model
##!              is zero correlation.}
##! \usage{
##! inla.pc.rcor0(n, u, alpha, lambda)
##! inla.pc.dcor0(cor, u, alpha, lambda, log = FALSE)
##! inla.pc.qcor0(p, u, alpha, lambda)
##! inla.pc.pcor0(q, u, alpha, lambda)
##! }
##! \arguments{
##!   \item{n}{Number of observations}
##!   \item{u}{The upper limit (see Details)}
##!   \item{alpha}{The probability going above the upper limit (see Details)}
##!   \item{lambda}{The rate parameter (see Details)}
##!   \item{cor}{Vector of correlations}
##!   \item{log}{Logical. Return the density in natural or log-scale.}
##!   \item{p}{Vector of probabilities}
##!   \item{q}{Vector of quantiles}
##! }
##! \details{
##!     The statement \code{Prob(|cor| > u) = alpha} is used to
##!     determine \code{lambda} unless \code{lambda} is given.
##!     Either \code{lambda}  must be given,  or
##!     \code{u} AND \code{alpha}. The density is symmetric around zero.
##! }
##!\value{%%
##!  \code{inla.pc.dcor0} gives the density,
##!  \code{inla.pc.pcor0} gives the distribution function,
##!  \code{inla.pc.qcor0} gives the quantile function, and
##!  \code{inla.pc.rcor0} generates random deviates.
##! }
##! \seealso{inla.doc("pc.rho0")}
##! \author{Havard Rue \email{hrue@r-inla.org}}
##! \examples{
##! cor = inla.pc.rcor0(100,  lambda = 1)
##! d = inla.pc.dcor0(cor, lambda = 1)
##! cor = inla.pc.qcor0(c(0.3, 0.7), u = 0.5, alpha=0.01)
##! inla.pc.pcor0(cor, u = 0.5, alpha=0.01)
##! }

inla.pc.cor0.lambda = function(u, alpha, lambda)
{
    if (missing(lambda)) {
        stopifnot(!missing(u) && !missing(alpha))
        lambda = -log(alpha)/sqrt(-log(1-u^2))
    }
    return (lambda)
}

inla.pc.rcor0 = function(n, u, alpha, lambda)
{
    lambda = inla.pc.cor0.lambda(u, alpha, lambda)
    d = rexp(n, rate = lambda)
    sign = sample(c(-1, 1), size = n, replace = TRUE)
    cor = sqrt(1-exp(-d^2)) * sign
    return (cor)
}

inla.pc.dcor0 = function(cor, u, alpha, lambda, log = FALSE)
{
    lambda = inla.pc.cor0.lambda(u, alpha, lambda)
    mu = sqrt(-log(1-cor^2))
    jac = abs(cor)/mu/(1-cor^2)
    if (log) {
        return (dexp(mu, rate = lambda, log=TRUE) + log(jac/2.0))
    } else {
        return (dexp(mu, rate = lambda) * jac/2.0)
    }
}

inla.pc.qcor0 = function(p, u, alpha, lambda)
{
    lambda = inla.pc.cor0.lambda(u, alpha, lambda)
    sign = (p >= 0.5)
    pp = sign * p + (1-sign)*(1-p)
    pp = 2*pp - 1
    qq = qexp(pp, rate = lambda)
    q = sqrt(1-exp(-qq^2))
    q = sign*q + (1-sign)*(-q)
    return (q)
}

inla.pc.pcor0 = function(q, u, alpha, lambda)
{
    lambda = inla.pc.cor0.lambda(u, alpha, lambda)
    sign = (q >= 0)
    qq = sign*q + (1-sign)*(-q)
    pp = sqrt(-log(1-qq^2))
    pp = pexp(pp, rate = lambda)
    p = sign*(1/2 + 1/2*pp) + (1-sign)*(1/2-1/2*pp)
    return (p)
}
