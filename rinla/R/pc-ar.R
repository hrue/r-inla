## Export: inla.pc.ar.rpacf inla.pc.ar.dpacf

##! \name{pc.ar}
##! \alias{inla.pc.ar}
##! \alias{inla.pc.ar.rpacf}
##! \alias{ar.rpacf}
##! \alias{inla.pc.ar.dpacf}
##! \alias{ar.dpacf}
##! 
##! \title{Utility functions for the PC prior for a an AR(p) model}
##! 
##! \description{Functions to evaluate and sample from the 
##!              PC prior for an AR(p) model}
##! \usage{
##!    inla.pc.ar.rpacf(n=1, p, lambda = 1)
##!    inla.pc.ar.dpacf(pac, lambda = 1, log = TRUE)
##! }
##! \arguments{
##!   \item{p}{The order of the AR-model}
##!   \item{pac}{A vector of partial autocorrelation coefficients}
##!   \item{n}{Number of observations}
##!   \item{lambda}{The rate parameter in the prior}
##!   \item{log}{Logical. Return the density in natural or log-scale.}
##! }
##! \details{
##! }
##! \value{%%
##!     \code{inla.pc.ar.rpac} generate samples from the prior, returning a matrix
##!     where each row is a sample of \code{theta}.
##!     \code{inla.pc.ar.dpac} evaluates the density of \code{pac}.
##!     Use \code{inla.ar.pacf2phi}, \code{inla.ar.phi2pacf},
##!     \code{inla.ar.pacf2acf} and \code{inla.ar.acf2pacf} to convert between various
##!    parameterisations.
##! }
##! \author{Havard Rue \email{hrue@math.ntnu.no}}
##! \examples{
##! }

inla.pc.ar.rpacf = function(n=1, p, lambda = 1)
{
    stopifnot(!missing(p) && p > 1)
    stopifnot(lambda > 0)
    stopifnot(n >= 1)
    
    rpac = function(n, p, lambda)
    {
        r = rexp(1, rate = lambda)
        gamma = rexp(p)
        gamma = (gamma/sum(gamma)) * r^2
        pac = sqrt(1-exp(-gamma)) * sample(c(-1, 1), size = p, replace=TRUE)
        return (pac)
    }
    
    x = unlist(lapply(rep(1, n), rpac, p=p, lambda = lambda))
    x = matrix(x, n, p, byrow=TRUE)
    return (x)
}

inla.pc.ar.dpacf = function(pac, lambda = 1, log = TRUE)
{
    p = length(pac)
    stopifnot(p > 1)
    d = sqrt(sum(-log(1-pac^2)))
    ld = (log(lambda) 
          + lfactorial(p) 
          + (5*p/2 + 1) * log(2.0)
          - 0.5 * log(p+1) 
          - (2*p+1) * log(d)
          - lambda * d
          + sum(log(abs(pac/(1-pac^2)))))
    
    return (if (log) ld else exp(ld))
}    
