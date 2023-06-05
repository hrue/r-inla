#' Utility functions for the PC prior for correlation in AR(1)
#' 
#' Functions to evaluate, sample, compute quantiles and percentiles of the PC
#' prior for the correlation in the Gaussian AR(1) model where the base-model
#' is zero correlation.
#' 
#' The statement `Prob(|cor| > u) = alpha` is used to determine
#' `lambda` unless `lambda` is given.  Either `lambda` must be
#' given, or `u` AND `alpha`. The density is symmetric around zero.
#' 
#' @aliases inla.pc.cor0 pc.cor0 pc.rcor0 inla.pc.rcor0 pc.dcor0 inla.pc.dcor0
#' pc.pcor0 inla.pc.pcor0 pc.qcor0 inla.pc.qcor0
#' @param n Number of observations
#' @param u The upper limit (see Details)
#' @param alpha The probability going above the upper limit (see Details)
#' @param lambda The rate parameter (see Details)
#' @param cor Vector of correlations
#' @param log Logical. Return the density in natural or log-scale.
#' @param p Vector of probabilities
#' @param q Vector of quantiles
#' @returns `inla.pc.dcor0` gives the density, `inla.pc.pcor0`
#' gives the distribution function, `inla.pc.qcor0` gives the quantile
#' function, and `inla.pc.rcor0` generates random deviates.
#' @author Havard Rue \email{hrue@@r-inla.org}
#' @seealso inla.doc("pc.rho0")
#' @examples
#' 
#'  cor = inla.pc.rcor0(100,  lambda = 1)
#'  d = inla.pc.dcor0(cor, lambda = 1)
#'  cor = inla.pc.qcor0(c(0.3, 0.7), u = 0.5, alpha=0.01)
#'  inla.pc.pcor0(cor, u = 0.5, alpha=0.01)
#'  
#' @name pc.cor0
#' @rdname pc-cor0
NULL

inla.pc.cor0.lambda <- function(u, alpha, lambda) {
    if (missing(lambda)) {
        stopifnot(!missing(u) && !missing(alpha))
        lambda <- -log(alpha) / sqrt(-log(1 - u^2))
    }
    return(lambda)
}

#' @rdname pc-cor0
#' @export
inla.pc.rcor0 <- function(n, u, alpha, lambda) {
    lambda <- inla.pc.cor0.lambda(u, alpha, lambda)
    d <- rexp(n, rate = lambda)
    sign <- sample(c(-1, 1), size = n, replace = TRUE)
    cor <- sqrt(1 - exp(-d^2)) * sign
    return(cor)
}

#' @rdname pc-cor0
#' @export
inla.pc.dcor0 <- function(cor, u, alpha, lambda, log = FALSE) {
    lambda <- inla.pc.cor0.lambda(u, alpha, lambda)
    mu <- sqrt(-log(1 - cor^2))
    jac <- abs(cor) / mu / (1 - cor^2)
    if (log) {
        return(dexp(mu, rate = lambda, log = TRUE) + log(jac / 2.0))
    } else {
        return(dexp(mu, rate = lambda) * jac / 2.0)
    }
}

#' @rdname pc-cor0
#' @export
inla.pc.qcor0 <- function(p, u, alpha, lambda) {
    lambda <- inla.pc.cor0.lambda(u, alpha, lambda)
    sign <- (p >= 0.5)
    pp <- sign * p + (1 - sign) * (1 - p)
    pp <- 2 * pp - 1
    qq <- qexp(pp, rate = lambda)
    q <- sqrt(1 - exp(-qq^2))
    q <- sign * q + (1 - sign) * (-q)
    return(q)
}

#' @rdname pc-cor0
#' @export
inla.pc.pcor0 <- function(q, u, alpha, lambda) {
    lambda <- inla.pc.cor0.lambda(u, alpha, lambda)
    sign <- (q >= 0)
    qq <- sign * q + (1 - sign) * (-q)
    pp <- sqrt(-log(1 - qq^2))
    pp <- pexp(pp, rate = lambda)
    p <- sign * (1 / 2 + 1 / 2 * pp) + (1 - sign) * (1 / 2 - 1 / 2 * pp)
    return(p)
}
