#' Utility functions for the PC prior for \code{Gamma(1/a, 1/a)}
#' 
#' Functions to evaluate, sample, compute quantiles and percentiles of the PC
#' prior for \code{Gamma(1/a, 1/a)}
#' 
#' This gives the PC prior for the \code{Gamma(1/a, 1/a)} case, where
#' \code{a=0} is the base model.
#' 
#' @aliases inla.pc.gamma pc.gamma pc.rgamma inla.pc.rgamma pc.dgamma
#' inla.pc.dgamma pc.pgamma inla.pc.pgamma pc.qgamma inla.pc.qgamma
#' @param n Number of observations
#' @param lambda The rate parameter (see Details)
#' @param x Evaluation points
#' @param log Logical. Return the density in natural or log-scale.
#' @param p Vector of probabilities
#' @param q Vector of quantiles
#' @returns \code{inla.pc.dgamma} gives the density, \code{inla.pc.pgamma}
#' gives the distribution function, \code{inla.pc.qgamma} gives the quantile
#' function, and \code{inla.pc.rgamma} generates random deviates.
#' @author Havard Rue \email{hrue@@r-inla.org}
#' @seealso inla.doc("pc.gamma")
#' @examples
#' 
#'  x = inla.pc.rgamma(100,  lambda = 1)
#'  d = inla.pc.dgamma(x, lambda = 1)
#'  x = inla.pc.qgamma(0.5, lambda = 1)
#'  inla.pc.pgamma(x, lambda = 1)
#'  
#' @name pc.gamma
#' @rdname pc-gamma
NULL


## I think this could have be done better than using the generic routines, by tabulating the d()
## function.

inla.pc.gamma.intern <- function(lambda = 1) {
    ## return a marginal-object in the log(x)-scale. the important variable is x*lambda,
    log.x <- seq(-20, 10, len = 4096) - log(lambda)
    x <- exp(log.x)
    marg <- list(x = log.x, y = inla.pc.dgamma(x, lambda = lambda, log = FALSE) * x)
    return(marg)
}

#' @rdname pc-gamma
#' @export
inla.pc.rgamma <- function(n, lambda = 1) {
    log.x <- inla.rmarginal(n, inla.pc.gamma.intern(lambda = lambda))
    return(exp(log.x))
}

#' @rdname pc-gamma
#' @export
inla.pc.dgamma <- function(x, lambda = 1, log = FALSE) {
    inv.x <- 1 / x
    d <- sqrt(2 * (log(inv.x) - psigamma(inv.x)))
    ld <- -lambda * d - log(d) + log(lambda) - 2 * log(x) +
        log(psigamma(inv.x, deriv = 1) - x)
    return(if (log) ld else exp(ld))
}

#' @rdname pc-gamma
#' @export
inla.pc.qgamma <- function(p, lambda = 1) {
    log.x <- inla.qmarginal(p, inla.pc.gamma.intern(lambda = lambda))
    return(exp(log.x))
}

#' @rdname pc-gamma
#' @export
inla.pc.pgamma <- function(q, lambda = 1) {
    p <- inla.pmarginal(log(q), inla.pc.gamma.intern(lambda = lambda))
    return(p)
}
