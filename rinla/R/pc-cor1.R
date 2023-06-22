#' Utility functions for the PC prior for correlation in AR(1)
#' 
#' Functions to evaluate, sample, compute quantiles and percentiles of the PC
#' prior for the correlation in the Gaussian AR(1) model where the base-model
#' is correlation one.
#' 
#' The statement `Prob(cor > u) = alpha` is used to determine
#' `lambda` unless `lambda` is given.  Either `lambda` must be
#' given, or `u` AND `alpha`.
#' 
#' @aliases inla.pc.cor1 pc.cor1 pc.rcor1 inla.pc.rcor1 pc.dcor1 inla.pc.dcor1
#' pc.pcor1 inla.pc.pcor1 pc.qcor1 inla.pc.qcor1
#' @param n Number of observations
#' @param u The upper limit (see Details)
#' @param alpha The probability going above the upper limit (see Details)
#' @param lambda The rate parameter (see Details)
#' @param cor Vector of correlations
#' @param log Logical. Return the density in natural or log-scale.
#' @param p Vector of probabilities
#' @param q Vector of quantiles
#' @returns `inla.pc.dcor1` gives the density, `inla.pc.pcor1`
#' gives the distribution function, `inla.pc.qcor1` gives the quantile
#' function, and `inla.pc.rcor1` generates random deviates.
#' @author Havard Rue \email{hrue@@r-inla.org}
#' @seealso inla.doc("pc.rho1")
#' @examples
#' 
#'  cor = inla.pc.rcor1(100,  lambda = 1)
#'  d = inla.pc.dcor1(cor, lambda = 1)
#'  cor = inla.pc.qcor1(c(0.3, 0.7), u = 0.5, alpha=0.75)
#'  inla.pc.pcor1(cor, u = 0.5, alpha=0.75)
#'  
#' @name pc.cor1
#' @rdname pc-cor1
NULL

inla.pc.cor1.lambda <- function(u, alpha, lambda) {
    if (missing(lambda)) {
        stopifnot(!missing(u) && !missing(alpha))
        alpha.min <- sqrt((1 - u) / 2)
        if (!(alpha > alpha.min)) {
              stop(paste("inla.pc.cor1.lambda: alpha >", alpha.min))
          }
        fun <- function(lam, u, alpha) {
            F <- (1 - exp(-lam * sqrt(1 - u))) / (1 - exp(-lam * sqrt(2)))
            return((F - alpha)^2)
        }
        lambdas <- unique(c(seq(0.0000001, 10, length.out = 100),
                            seq(10, 100, length.out = 10)))
        idx <- which.min(fun(lambdas, u, alpha))
        stopifnot(idx > 1 && idx < length(lambdas))
        lambda <- optimize(fun,
            interval = lambdas[c(idx - 1, idx + 1)],
            maximum = FALSE, u = u, alpha = alpha
        )$minimum
        stopifnot(fun(lambda, u, alpha) < 1e-8)
    }
    return(lambda)
}

#' @rdname pc-cor1
#' @export
inla.pc.rcor1 <- function(n, u, alpha, lambda) {
    lambda <- inla.pc.cor1.lambda(u, alpha, lambda)
    u <- runif(n)
    cor <- -1 / lambda * log(1 - u * (1 - exp(-sqrt(2) * lambda)))
    cor <- 1 - cor^2
    return(cor)
}

#' @rdname pc-cor1
#' @export
inla.pc.dcor1 <- function(cor, u, alpha, lambda, log = FALSE) {
    lambda <- inla.pc.cor1.lambda(u, alpha, lambda)
    log.dens <- log(lambda) - lambda * sqrt(1 - cor) -
        log(1 - exp(-sqrt(2) * lambda)) - log(2 * sqrt(1 - cor))
    return(inla.ifelse(log, log.dens, exp(log.dens)))
}

#' @rdname pc-cor1
#' @export
inla.pc.qcor1 <- function(p, u, alpha, lambda) {
    lambda <- inla.pc.cor1.lambda(u, alpha, lambda)
    pp <- -1 / lambda * log(1 - (1 - p) * (1 - exp(-lambda * sqrt(2))))
    q <- 1 - pp^2
    return(q)
}

#' @rdname pc-cor1
#' @export
inla.pc.pcor1 <- function(q, u, alpha, lambda) {
    lambda <- inla.pc.cor1.lambda(u, alpha, lambda)
    qq <- (1 - exp(-lambda * sqrt(1 - q))) / (1 - exp(-lambda * sqrt(2)))
    return(1 - qq)
}
