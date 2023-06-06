#' Utility functions for the PC prior for the precision
#' 
#' Functions to evaluate, sample, compute quantiles and percentiles of the PC
#' prior for the precision in the Gaussian distribution.
#' 
#' The statement `Prob(1/sqrt(prec) > u) = alpha` is used to determine
#' `lambda` unless `lambda` is given.  Either `lambda` must be
#' given, or `u` AND `alpha`.
#' 
#' @aliases inla.pc.prec pc.prec pc.rprec inla.pc.rprec pc.dprec inla.pc.dprec
#' pc.pprec inla.pc.pprec pc.qprec inla.pc.qprec
#' @param n Number of observations
#' @param u The upper limit (see Details)
#' @param alpha The probability going above the upper limit (see Details)
#' @param lambda The rate parameter (see Details)
#' @param prec Vector of precisions
#' @param log Logical. Return the density in natural or log-scale.
#' @param p Vector of probabilities
#' @param q Vector of quantiles
#' @return `inla.pc.dprec` gives the density, `inla.pc.pprec`
#' gives the distribution function, `inla.pc.qprec` gives the quantile
#' function, and `inla.pc.rprec` generates random deviates.
#' @author Havard Rue \email{hrue@@r-inla.org}
#' @seealso inla.doc("pc.prec")
#' @examples
#' 
#'  prec = inla.pc.rprec(100,  lambda = 1)
#'  d = inla.pc.dprec(prec, lambda = 1)
#'  prec = inla.pc.qprec(0.5, u = 1, alpha=0.01)
#'  inla.pc.pprec(prec, u = 1, alpha=0.01)
#' 
#' @name pc.prec
#' @rdname pc-prec 
NULL

inla.pc.prec.lambda <- function(u, alpha, lambda) {
    if (missing(lambda)) {
        stopifnot(!missing(u) && !missing(alpha))
        lambda <- -log(alpha) / u
    }
    return(lambda)
}

#' @rdname pc-prec
#' @export
inla.pc.rprec <- function(n, u, alpha, lambda) {
    lambda <- inla.pc.prec.lambda(u, alpha, lambda)
    s <- rexp(n, rate = lambda)
    prec <- 1 / s^2
    return(prec)
}

#' @rdname pc-prec
#' @export
inla.pc.dprec <- function(prec, u, alpha, lambda, log = FALSE) {
    lambda <- inla.pc.prec.lambda(u, alpha, lambda)
    s <- 1 / sqrt(prec)
    d <- dexp(s, rate = lambda, log = log)
    log.jac <- -log(2) - 3 / 2 * log(prec)
    if (log) {
        return(d + log.jac)
    } else {
        return(d * exp(log.jac))
    }
}

#' @rdname pc-prec
#' @export
inla.pc.qprec <- function(p, u, alpha, lambda) {
    lambda <- inla.pc.prec.lambda(u, alpha, lambda)
    q <- qexp(1 - p, rate = lambda)
    qprec <- 1 / q^2
    return(qprec)
}

#' @rdname pc-prec
#' @export
inla.pc.pprec <- function(q, u, alpha, lambda) {
    lambda <- inla.pc.prec.lambda(u, alpha, lambda)
    p <- pexp(1 / sqrt(q), rate = lambda)
    return(1 - p)
}
