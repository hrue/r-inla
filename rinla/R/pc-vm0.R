#' Utility functions for the PC prior for the precision
#' 
#' Functions to evaluate, sample, compute quantiles and percentiles of the PC prior for the
#' precision in the Gaussian distribution.
#' 
#' The statement `Prob(2*pi/(1+k) > u) = alpha` is used to determine `lambda` unless `lambda` is
#' given. Either `lambda` must be given, or `u` AND `alpha`.
#' 
#' @details Due to limitations in handling extreme values for special functions, the output of
#'     these functions may exhibit bias when the input parameter values are either excessively
#'     large or very close to zero.
#' 
#' @aliases inla.pc.vm0 pc.vm0 pc.rvm0 inla.pc.rvm0 pc.dvm0 
#' inla.pc.dvm0 pc.pvm0 inla.pc.pvm0 pc.qvm0 inla.pc.qvm0
#' @param n Number of observations
#' @param u The upper limit (0 < u < 2*pi). The small values of u indicate a 
#' high concentration to a point mass, whilst large values of u mean that the 
#' user believes the data spread widely.
#' @param alpha The probability going above the upper limit (the probability 
#' assigned to the event Prob(2*pi/(1+k) > u)).
#' @param lambda The rate parameter (see Details)
#' @param k The concentration of von Mises distribution
#' @param log Logical. Return the density in natural or log-scale.
#' @param p Vector of probabilities
#' @param q Vector of quantiles
#' @return `inla.pc.dvm0` gives the density, `inla.pc.pvm0`
#' gives the distribution function, `inla.pc.qvm0` gives the quantile
#' function, and `inla.pc.rvm0` generates random deviates.
#' @author Xiang Ye \email{xiang.ye@@kaust.edu.sa}
#' @seealso inla.doc("pc.vm0")
#' @examples
#' 
#'  k = inla.pc.rvm0(100,  lambda = 1)
#'  d = inla.pc.dvm0(1, lambda = 1)
#'  k = inla.pc.qvm0(0.5, u = 1, alpha=0.01)
#'  inla.pc.pvm0(5, u = 1, alpha=0.01)
#' 
#' @name pc.vm0
#' @rdname pc-vm0 
NULL

inla.pc.vm0.lambda <- function(u, alpha, lambda) {
    if (missing(lambda)) {
        stopifnot(!missing(u) && !missing(alpha))
        k <- 2*pi/u - 1
        I0 <- besselI(x=k,nu=0)
        I1 <- besselI(x=k,nu=1)
        lambda <- -log(1-alpha) / sqrt(k*I1/I0-log(I0))
    }
    return(lambda)
}

#' @rdname pc-vm0
#' @export
inla.pc.rvm0 <- function(n, u, alpha, lambda) {
    return(inla.pc.qvm0(runif(n),u,alpha,lambda))
}

#' @rdname pc-vm0
#' @export
inla.pc.dvm0 <- function(k, u, alpha, lambda, log = FALSE) {
    if (any(log(k) > .Machine$double.xmax)) {
        warning("log(k) exceeds log(.Machine$double.xmax); k is too large.")
    }
    lambda <- inla.pc.vm0.lambda(u, alpha, lambda)
    f_single <- function(kk) {
        if (kk < 1e-4) {
            dens <- lambda/2 - lambda^2 * kk/4 +
                (-9 * lambda/64 + lambda^3/16) * kk^2 +
                (12 * lambda^2/128 - lambda^4/96) * kk^3 +
                (1195 * lambda/36864 + 3 * lambda^3/512 + lambda^5/768) * kk^4 +
                (-79 * lambda^2/6144 - lambda^6/7680) * kk^5 # Here I remain to 5 order first for you to check
            if (log) {
                return(log(dens))
            } else {
                return(dens)
            }
        } else {
                                        # For k >= 1e-4
            I0 <- besselI(x = kk, nu = 0)
            I1 <- besselI(x = kk, nu = 1)
            I2 <- besselI(x = kk, nu = 2)
            
            dis.square <- if (kk <= 700) {
                              kk * I1 / I0 - log(I0)
                          } else {
                              -1/2 + kk - (kk - 1/2 * log(2*pi)) +
                                  log(kk)/2 - 1/(4*kk) - 3/(16*kk^2) -
                                      25/(96*kk^3) - 65/(128*kk^4) - 3219/(2560*kk^5)
                          }
            
            log.jac.partial <- if (kk <= 350) {
                                   log(kk * (I0 + I2)/(2 * I0) - kk * I1^2/I0^2)
                               } else {
                                   -log(2) - log(kk) + 1/(2*kk) + 5/(8*kk^2) +
                                       59/(48*kk^3) + 203/(64*kk^4) + 12743/(1280*kk^5)
                               }
            
            distance <- sqrt(dis.square)
            log.dens.exp <- dexp(distance, rate = lambda, log = TRUE)
            log.jac <- -log(2) - log(distance) + log.jac.partial
            
            ldens <- log.dens.exp + log.jac
            return (if (log) ldens else exp(ldens))
        }
    }
    f_vec <- Vectorize(f_single)
    return(f_vec(k))
}

#' @rdname pc-vm0
#' @export
inla.pc.qvm0 <- function(p, u, alpha, lambda, len = 2048L) {
    lambda <- inla.pc.vm0.lambda(u, alpha, lambda)
    q.grid <- exp(seq(-16, 22, length.out = len))
    Fq <- inla.pc.pvm0(q.grid, u, alpha, lambda)
    Fq.log <- inla.pc.pvm0(q.grid, u, alpha, lambda,log=TRUE)
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

#' @rdname pc-vm0
#' @export
inla.pc.pvm0 <- function(q, u, alpha, lambda, log = FALSE) {
    if (any(log(q) > .Machine$double.xmax)) {
        warning("log(q) exceeds log(.Machine$double.xmax); q is too large.")
    }
    lambda <- inla.pc.vm0.lambda(u, alpha, lambda)
    I0 <- besselI(x=q,nu=0)
    I1 <- besselI(x=q,nu=1)
    dis.square <- ifelse(q <= 700,
                         q*I1/I0-log(I0),
                         -1/2 + q - (q - 1/2*log(2*pi)) + log(q)/2 - 1/(4*q) - 3/(16*q^2) - 25/(96*q^3) - 65/(128*q^4) - 3219/(2560*q^5) )  # Here I remain to 5 order first for you to check
    distance <- sqrt(dis.square) # It is very hard to get the asympt behavior of distance and log-distance, so I tried for square of distance
    p <- pexp(distance, rate = lambda, log.p = log)
    return(p)
}
