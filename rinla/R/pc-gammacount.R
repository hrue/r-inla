#' Utility functions for the PC prior for the \code{gammacount} likelihood
#' 
#' Functions to evaluate, sample, compute quantiles and percentiles of the PC
#' prior for the \code{gammacount} likelihood
#' 
#' This gives the PC prior for the \code{gammacount} likelihood, which is the
#' PC prior for \code{a} in \code{Gamma(a, 1)} where \code{Gamma(1, 1)} is the
#' base model.
#' 
#' @aliases inla.pc.gammacount pc.gammacount pc.rgammacount inla.pc.rgammacount
#' pc.dgammacount inla.pc.dgammacount pc.pgammacount inla.pc.pgammacount
#' pc.qgammacount inla.pc.qgammacount
#' @param n Number of observations
#' @param lambda The rate parameter (see Details)
#' @param x Evaluation points
#' @param log Logical. Return the density in natural or log-scale.
#' @param p Vector of probabilities
#' @param q Vector of quantiles
#' @returns \code{inla.pc.dgammacount} gives the density,
#' \code{inla.pc.pgammacount} gives the distribution function,
#' \code{inla.pc.qgammacount} gives the quantile function, and
#' \code{inla.pc.rgammacount} generates random deviates.
#' @author Havard Rue \email{hrue@@r-inla.org}
#' @seealso inla.doc("pc.gammacount")
#' @examples
#' 
#'  x = inla.pc.rgammacount(100,  lambda = 1)
#'  d = inla.pc.dgammacount(x, lambda = 1)
#'  x = inla.pc.qgammacount(0.5, lambda = 1)
#'  inla.pc.pgammacount(x, lambda = 1)
#'  
#' @name pc.gammacount
#' @rdname pc-gammacount
NULL


## I think this could have be done better than using the generic routines, by tabulating the d()
## function.

inla.pc.gammacount.intern <- function(lambda = 1) {
    ## return a marginal-object in the log(x)-scale
    f <- if (lambda >= 1.0) 1.0 else 0.25
    log.x <- (seq(-10, 10, len = 4096)) / lambda^f
    idx <- which(log.x == 0.0)
    if (length(idx) > 0) {
        log.x <- log.x[-idx]
    }
    x <- exp(log.x)
    marg <- list(x = log.x, y = inla.pc.dgammacount(x, lambda = lambda, log = FALSE) * x)
    return(marg)
}

#' @rdname pc-gammacount
#' @export
inla.pc.rgammacount <- function(n, lambda = 1) {
    log.x <- inla.rmarginal(n, inla.pc.gammacount.intern(lambda = lambda))
    return(exp(log.x))
}

#' @rdname pc-gammacount
#' @export
inla.pc.dgammacount <- function(x, lambda = 1, log = FALSE) {
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

    return(if (log) ld else exp(ld))
}

#' @rdname pc-gammacount
#' @export
inla.pc.qgammacount <- function(p, lambda = 1) {
    log.x <- inla.qmarginal(p, inla.pc.gammacount.intern(lambda = lambda))
    return(exp(log.x))
}

#' @rdname pc-gammacount
#' @export
inla.pc.pgammacount <- function(q, lambda = 1) {
    p <- inla.pmarginal(log(q), inla.pc.gammacount.intern(lambda = lambda))
    return(p)
}
