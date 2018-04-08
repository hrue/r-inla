## Export: inla.pc.rgamma inla.pc.dgamma inla.pc.qgamma inla.pc.pgamma

##! \name{pc.gamma}
##! \alias{inla.pc.gamma}
##! \alias{pc.gamma}
##! \alias{pc.rgamma}
##! \alias{inla.pc.rgamma}
##! \alias{pc.dgamma}
##! \alias{inla.pc.dgamma}
##! \alias{pc.pgamma}
##! \alias{inla.pc.pgamma}
##! \alias{pc.qgamma}
##! \alias{inla.pc.qgamma}
##! 
##! \title{Utility functions for the PC prior for \code{Gamma(1/a, 1/a)}}
##! 
##! \description{Functions to evaluate, sample, compute quantiles and
##!              percentiles of the PC prior for \code{Gamma(1/a, 1/a)}}
##! \usage{
##! inla.pc.rgamma(n, lambda = 1)
##! inla.pc.dgamma(x, lambda = 1, log = FALSE)
##! inla.pc.qgamma(p, lambda = 1)
##! inla.pc.pgamma(q, lambda = 1)
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
##! This gives the PC prior for the \code{Gamma(1/a, 1/a)} case,  where \code{a=0} is the
##! base model. 
##! }
##!\value{%%
##!  \code{inla.pc.dgamma} gives the density,
##!  \code{inla.pc.pgamma} gives the distribution function,
##!  \code{inla.pc.qgamma} gives the quantile function, and
##!  \code{inla.pc.rgamma} generates random deviates.
##! }
##! \seealso{inla.doc("pc.gamma")}
##! \author{Havard Rue \email{hrue@r-inla.org}}
##! \examples{
##! x = inla.pc.rgamma(100,  lambda = 1)
##! d = inla.pc.dgamma(x, lambda = 1)
##! x = inla.pc.qgamma(0.5, lambda = 1)
##! inla.pc.pgamma(x, lambda = 1)
##! }

## I think this could have be done better than using the generic routines, by tabulating the d()
## function. 

inla.pc.gamma.intern = function(lambda = 1)
{
    ## return a marginal-object in the log(x)-scale. the important variable is x*lambda, 
    log.x = seq(-20, 10, len=4096) - log(lambda)
    x = exp(log.x)
    marg = list(x = log.x, y = inla.pc.dgamma(x, lambda = lambda, log=FALSE) * x)
    return (marg)
}

inla.pc.rgamma = function(n, lambda = 1)
{
    log.x = inla.rmarginal(n, inla.pc.gamma.intern(lambda = lambda))
    return (exp(log.x))
}

inla.pc.dgamma = function(x, lambda = 1, log = FALSE)
{
    inv.x = 1/x
    d = sqrt(2*(log(inv.x) - psigamma(inv.x)))
    ld = -lambda * d - log(d) + log(lambda) - 2*log(x) +
        log(psigamma(inv.x, deriv=1) - x)
    return (if (log) ld else exp(ld))
}

inla.pc.qgamma = function(p, lambda = 1)
{
    log.x = inla.pmarginal(p, inla.pc.gamma.intern(lambda = lambda))
    return (exp(log.x))
}

inla.pc.pgamma = function(q, lambda = 1)
{
    p = inla.qmarginal(log(q), inla.pc.gamma.intern(lambda = lambda))
    return (p)
}
