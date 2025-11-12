#' Utility functions for the PC prior for range in PRW2
#' 
#' Functions to evaluate, sample, compute quantiles and percentiles of the PC prior for the
#' range parameter in the PRW2 model.
#' 
#' The statement `Prob(range > r0) = alpha` is used to determine
#' `lambda` unless `lambda` is given.  
#' 
#' @aliases inla.pc.prw2.drange inla.pc.prw2.rrange inla.pc.prw2.qrange inla.pc.prw2.prange
#' @param n Number of observations
#' @param param Vector of parameters (see 'inla.doc("pc.prw2.range")')
#' @param p Vector of probabilities
#' @param q Vector of quantiles
#' @param x Vector of quantiles
#' @param p Vector of probabilities
#' @returns `inla.prw2.drange` gives the density, `inla.prw2.prange` gives the distribution
#'     function, `inla.prw2.qrange` gives the quantile function, and `inla.prw2.rrange`
#'     generates random deviates.
#' @note The 'r' and 'q' functions are quite slow as the quantile-function needs to be computed.
#' @author Havard Rue \email{hrue@@r-inla.org}
#' @seealso inla.doc("pc.prw2.range")
#' @examples
#'
#'     param <- c(10, 0.5, 1, 0)
#'     p <- inla.prw2.prange(10, param)
#'     q <- inla.prw2.qrange(p, param)
#'     print(c(10, q))
#'  
#' @name pc.prw2.range
#' @rdname prw2.range

prw2.d <- function(rho) {
    lrho <- log1p(rho-1)
    return (sqrt(-(3+rho)* expm1(lrho)^3 / (1+rho^2)))
}

prw2.range2rho <- function(range, hsize) {
    return (exp(-sqrt(12) * hsize / range))
}

prw2.F <- function(range, lambda, hsize) {
    v <- exp(-sqrt(3) * lambda)
    return ((exp(-lambda * prw2.d(prw2.range2rho(range, hsize))) - v) / (1.0 - v))
}
  
prw2.check.param <- function(param) {
    if (param[1] <= 0 || param[2] <= 0 || param[3] <= 0) {
        stop("inla.prw2.{p, d, r, q}range does not support undefined values in 'param'")
    }
}

prw2.calibrate.lambda <- function(param) {
    prw2.check.param(param)
    if (is.nan(param[4]) || param[4] <= 0) {
        ## find lambda so that Prob(range > range0) = alpha
        fun.opt2 <- function(log.lambda,  range0, alpha, hsize) {
            F <- prw2.F(range0, exp(log.lambda), hsize)
            return ((log(F/(1.0 - F)) - log(alpha / (1.0 - alpha)))^2)
        }
        param[4] <- exp(optim(0, gr = NULL, fun.opt2, method = "BFGS",
                              range0 = param[1], alpha = param[2], hsize = param[3])$par)
    }
    return (param)
}

prw2.get.values <- function(param,  eps = 1e-6) {
    param <- prw2.calibrate.lambda(param)
    iter <- 1
    iter.max <- 10000

    fac <- 1.5
    mid <- 1
    pmid <- inla.prw2.prange(mid, param)

    low <- mid
    plow <- pmid
    while (plow > eps) {
        low <- low / fac
        ploww <- inla.prw2.prange(low, param)
        if (ploww == plow) break
        plow <- ploww
        iter <- iter + 1
        stopifnot(iter < iter.max)
    }

    high <- mid
    phigh <- pmid
    while (phigh < 1-eps) {
        high <- high * fac
        phighh <- inla.prw2.prange(high, param)
        if (phighh == phigh) break
        phigh <- phighh
        iter <- iter + 1
        stopifnot(iter < iter.max)
    }
    return (exp(seq(log(low), log(high), len = 2048)))
}

prw2.rF <- function(param) {
    r <- prw2.get.values(param)
    F <- inla.prw2.prange(r, param)

    idx <- which(is.finite(F))
    r <- r[idx]
    F <- F[idx]

    small <- (1.0 + F == F)
    r <- r[!small]
    F <- F[!small]

    not.dup <- !duplicated(F)
    r <- r[not.dup]
    F <- F[not.dup]
    
    eps <- 1e-12
    F <- (F - F[1] + eps) / (F[length(F)] - F[1] + 2*eps)

    return (list(r = r, F = F))
}

#' @export
inla.prw2.prange <- function(q, param) {
    param <- prw2.calibrate.lambda(param)
    return(prw2.F(q, param[4], param[3]))
}

#' @export
inla.prw2.drange <- function(r, param) {
    param <- prw2.calibrate.lambda(param)
    h <- 1e-4
    dd <- (prw2.F(r*exp(h), param[4], param[3]) -
           prw2.F(r*exp(-h), param[4], param[3])) / (2*h) / r
    return (dd)
}

#' @export
inla.prw2.qrange <- function(p, param) {
    d <- prw2.rF(param)
    fun <- splinefun(log(d$F), log(d$r), method = "monoH.FC")
    return (exp(fun(log(p))))
}

#' @export
inla.prw2.rrange <- function(n, param) {
    return (inla.prw2.qrange(runif(n), param))
}
