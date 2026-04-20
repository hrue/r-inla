#' Utility functions for the PC prior for a an AR(p) model
#' 
#' Functions to evaluate and sample from the PC prior for an AR(p) model
#' 
#' 
#' @aliases inla.pc.ar inla.pc.ar.rpacf ar.rpacf inla.pc.ar.dpacf ar.dpacf
#' @param p The order of the AR-model
#' @param pac A vector of partial autocorrelation coefficients
#' @param n Number of observations
#' @param lambda The rate parameter in the prior
#' @param log Logical. Return the density in natural or log-scale.
#' @return `inla.pc.ar.rpac` generate samples from the prior, returning
#' a matrix where each row is a sample of `theta`.  `inla.pc.ar.dpac`
#' evaluates the density of `pac`.  Use `inla.ar.pacf2phi`,
#' `inla.ar.phi2pacf`, `inla.ar.pacf2acf` and `inla.ar.acf2pacf`
#' to convert between various parameterisations.
#' @author Havard Rue \email{hrue@@r-inla.org}
#' 
#' @name pc.ar
#' @rdname pc-ar
#' @export

inla.pc.ar.rpacf <- function(n = 1, p, lambda = 1) {
    stopifnot(!missing(p) && p >= 1)
    stopifnot(lambda > 0)
    stopifnot(n >= 1)

    rpac <- function(n, p, lambda) {
        gamma <- inla.pc.multvar.simplex.r(n = n, b = rep(1 / 2, p), lambda = lambda)
        pac <- sqrt(1 - exp(-gamma)) * sample(c(-1, 1), size = p, replace = TRUE)
        return(pac)
    }

    x <- unlist(lapply(rep(1, n), rpac, p = p, lambda = lambda))
    x <- matrix(x, n, p, byrow = TRUE)
    return(x)
}

#' @rdname pc-ar
#' @export
inla.pc.ar.dpacf <- function(pac, lambda = 1, log = TRUE) {
    p <- length(pac)
    stopifnot(p >= 1)
    gamma <- -log(1 - pac^2)
    ld <- (inla.pc.multvar.simplex.d(x = gamma, lambda = lambda, b = rep(1 / 2, p), log = TRUE)
    +
        sum(log(abs(pac / (1 - pac^2)))))
    return(if (log) ld else exp(ld))
}
