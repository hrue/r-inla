#' Utility functions for the PC prior for the concentration of von Mises distribution with point
#' mass base model
#' 
#' Functions to evaluate, sample, compute quantiles and percentiles of the PC prior for the
#' concentration in the von Mises distribution.
#' 
#' The statement `Prob(2*pi/(1+k) > u) = alpha` is used to determine `lambda` unless `lambda` is
#' given. Either `lambda` must be given, or `u` AND `alpha`.
#' 
#' @details Due to limitations in handling extreme values for special functions, the output of
#'     these functions may exhibit bias when the input parameter values are either excessively
#'     large or very close to zero.
#' 
#' @aliases inla.pc.vminf pc.vminf pc.rvminf inla.pc.rvminf pc.dvminf 
#' inla.pc.dvminf pc.pvminf inla.pc.pvminf pc.qvminf inla.pc.qvminf
#' @param n Number of observations
#' @param u The upper limit (0 < u < 2*pi). The small values of u indicate a 
#' high concentration to a point mass, whilst large values of u mean that the 
#' user believes the data spread widely.
#' @param alpha The probability going above the upper limit (the probability 
#' assigned to the event Prob(2*pi/(1+k) > u)).
#' @param lambda The rate parameter.
#' @param k The concentration of von Mises distribution
#' @param log Logical. Return the density in natural or log-scale.
#' @param p Vector of probabilities.
#' @param q Vector of quantiles.
#' @return `inla.pc.dvminf` gives the density, `inla.pc.pvminf`
#' gives the distribution function, `inla.pc.qvminf` gives the quantile
#' function, and `inla.pc.rvminf` generates random deviates.
#' @author Xiang Ye \email{xiang.ye@@kaust.edu.sa}
#' @seealso inla.doc("pc.vminf")
#' @examples
#' 
#'  k = inla.pc.rvminf(100,  lambda = 1)
#'  d = inla.pc.dvminf(1, lambda = 1)
#'  k = inla.pc.qvminf(0.5, u = 1, alpha=0.01)
#'  inla.pc.pvminf(5, u = 1, alpha=0.01)
#' 
#' @name pc.vminf
#' @rdname pc-vminf 
NULL

inla.pc.vminf.lambda <- function(u, alpha, lambda) {
    if (missing(lambda)) {
        stopifnot(!missing(u) && !missing(alpha))
        lambda <- -log(1 - alpha) / sqrt(1 - u)
    }
    return(lambda)
}

#' @rdname pc-vminf
#' @export
inla.pc.rvminf <- function(n, u, alpha, lambda) {
    return(inla.pc.qvminf(runif(n),u,alpha,lambda))
}

#' @rdname pc-vminf
#' @export
inla.pc.dvminf <- function(k, u, alpha, lambda, log = FALSE) {
    if (any(log(k) > .Machine$double.xmax)) {
        warning("log(k) exceeds log(.Machine$double.xmax); k is too large.")
    }
    lambda <- inla.pc.vminf.lambda(u, alpha, lambda)
    I0.exp <- besselI(x=k,nu=0,expon.scaled = TRUE)
    I1.exp <- besselI(x=k,nu=1,expon.scaled = TRUE)
    I2.exp <- besselI(x=k,nu=2,expon.scaled = TRUE)
    log.distance.asympt <- -log(2)/2 - log(k)/2 + 1/(8*k) + 7/(64*k^2) + 
        1/(6*k^3) + 715/(2048*k^4)+293/(320*k^5)
    distance <- ifelse(k<=1e+5,
                       sqrt(1-I1.exp/I0.exp),
                       exp(log.distance.asympt))
    log.distance <- ifelse(k<=1e+5,
                           log(distance),
                           log.distance.asympt)
    log.dens.exp <- dexp(distance, rate = lambda, log = TRUE)
    log.jac.partial <- ifelse(k<=1e+5,
                              log((I0.exp+I2.exp)/(2*I0.exp)-I1.exp^2/I0.exp^2),
                              -log(2) - 2*log(k) + 1/(2*k) + 5/(8*k^2) + 
                              59/(48*k^3) + 203/(64*k^4) + 12743/(1280*k^5))
    log.jac <- -log(2) - log.distance + log.jac.partial
    if (log) {
        return(log.dens.exp + log.jac)
    } else {
        return(exp(log.dens.exp) * exp(log.jac))
    }
}

#' @rdname pc-vminf
#' @export
inla.pc.qvminf <- function(p, u, alpha, lambda, len = 2048L) {
    lambda <- inla.pc.vminf.lambda(u, alpha, lambda)
    q.grid <- exp(seq(-6, 66, length.out = len))
    Fq <- inla.pc.pvminf(q.grid, u, alpha, lambda)
    Fq.log <- inla.pc.pvminf(q.grid, u, alpha, lambda,log=TRUE)
    x <- Fq.log - log(1 - Fq)
    y <- log(q.grid)
    
    keep <- !duplicated(round(x, 8))
    x <- x[keep]; y <- y[keep]
    
    spline_inv <- splinefun(x, y, method = "monoH.FC")
    p <- pmin(pmax(p, 0), 1)
    q.out <- sapply(p, function(prob) {
        if (prob == 0) return(0)
        if (prob == 1) return(Inf)
        logit_p <- log(prob / (1 - prob))
        exp(spline_inv(logit_p))
    })
    return(q.out)
}

#' @rdname pc-vminf
#' @export
inla.pc.pvminf <- function(q, u, alpha, lambda, log=FALSE) {
    lambda <- inla.pc.vminf.lambda(u, alpha, lambda)
    I0.exp <- besselI(x=q,nu=0,expon.scaled = TRUE)
    I1.exp <- besselI(x=q,nu=1,expon.scaled = TRUE)
    distance <- ifelse(q<=1e+5,
                       sqrt(1-I1.exp/I0.exp),
                       sqrt(1/(2*q) + 1/(8*q^2) + 1/(8*q^3) + 
                            25/(128*q^4) + 13/(32*q^5)))
    if (log) {
        return(-lambda*distance)
    } else {
        return(1-pexp(distance, rate = lambda))
    }
}
