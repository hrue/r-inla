## Export: inla.pc.rcor1 inla.pc.dcor1 inla.pc.qcor1 inla.pc.pcor1

##! \name{pc.cor1}
##! \alias{inla.pc.cor1}
##! \alias{pc.cor1}
##! \alias{pc.rcor1}
##! \alias{inla.pc.rcor1}
##! \alias{pc.dcor1}
##! \alias{inla.pc.dcor1}
##! \alias{pc.pcor1}
##! \alias{inla.pc.pcor1}
##! \alias{pc.qcor1}
##! \alias{inla.pc.qcor1}
##! 
##! \title{Utility functions for the PC prior for correlation in AR(1)}
##! 
##! \description{Functions to evaluate, sample, compute quantiles and
##!              percentiles of the PC prior for the correlation
##!              in the Gaussian AR(1) model where the base-model
##!              is correlation one.}
##! \usage{
##! inla.pc.rcor1(n, u, alpha, lambda)
##! inla.pc.dcor1(cor, u, alpha, lambda, log = FALSE)
##! inla.pc.qcor1(p, u, alpha, lambda)
##! inla.pc.pcor1(q, u, alpha, lambda)
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
##!     The statement \code{Prob(cor > u) = alpha} is used to
##!     determine \code{lambda} unless \code{lambda} is given.
##!     Either \code{lambda}  must be given,  or
##!     \code{u} AND \code{alpha}. 
##! }
##!\value{%%
##!  \code{inla.pc.dcor1} gives the density,
##!  \code{inla.pc.pcor1} gives the distribution function,
##!  \code{inla.pc.qcor1} gives the quantile function, and
##!  \code{inla.pc.rcor1} generates random deviates.
##! }
##! \seealso{inla.doc("pc.rho1")}
##! \author{Havard Rue \email{hrue@r-inla.org}}
##! \examples{
##! cor = inla.pc.rcor1(100,  lambda = 1)
##! d = inla.pc.dcor1(cor, lambda = 1)
##! cor = inla.pc.qcor1(c(0.3, 0.7), u = 0.5, alpha=0.75)
##! inla.pc.pcor1(cor, u = 0.5, alpha=0.75)
##! }

inla.pc.cor1.lambda = function(u, alpha, lambda)
{
    if (missing(lambda)) {
        stopifnot(!missing(u) && !missing(alpha))
        alpha.min = sqrt((1-u)/2)
        if (!(alpha > alpha.min))
            stop(paste("inla.pc.cor1.lambda: alpha >", alpha.min))
        fun = function(lam, u, alpha) {
            F = (1-exp(-lam*sqrt(1-u)))/(1-exp(-lam*sqrt(2)))
            return ((F-alpha)^2)
        }
        lambdas = unique(c(seq(0.0000001, 10, len=100),  seq(10, 100, len=10)))
        idx = which.min(fun(lambdas, u, alpha))
        stopifnot(idx>1 && idx < length(lambdas))
        lambda = optimize(fun, interval = lambdas[c(idx-1, idx+1)], 
                maximum = FALSE, u=u, alpha=alpha)$minimum
        stopifnot(fun(lambda, u, alpha) < 1e-8)
    }
    return (lambda)
}

inla.pc.rcor1 = function(n, u, alpha, lambda)
{
    lambda = inla.pc.cor1.lambda(u, alpha, lambda)
    u = runif(n)
    cor = -1/lambda * log(1 - u*(1-exp(-sqrt(2)*lambda)))
    cor = 1-cor^2
    return(cor)
}

inla.pc.dcor1 = function(cor, u, alpha, lambda, log = FALSE)
{
    lambda = inla.pc.cor1.lambda(u, alpha, lambda)
    log.dens = log(lambda) - lambda * sqrt(1-cor) -
        log(1 - exp(-sqrt(2)*lambda)) - log(2*sqrt(1-cor))
    return (inla.ifelse(log,  log.dens,  exp(log.dens)))
}

inla.pc.qcor1 = function(p, u, alpha, lambda)
{
    lambda = inla.pc.cor1.lambda(u, alpha, lambda)
    pp = -1/lambda * log(1 - (1-p) * (1-exp(-lambda*sqrt(2))))
    q = 1-pp^2
    return (q)
}

inla.pc.pcor1 = function(q, u, alpha, lambda)
{
    lambda = inla.pc.cor1.lambda(u, alpha, lambda)
    qq = (1-exp(-lambda*sqrt(1-q)))/(1-exp(-lambda*sqrt(2)))
    return (1-qq)
}
