#' Utility functions for the PC prior for the precision
#' 
#' Functions to evaluate, sample, compute quantiles and percentiles of the PC prior for the
#' precision in the Gaussian distribution.
#' 
#' The statement `Prob(2*pi/(1+k) > u) = alpha` is used to determine
#' `lambda` unless `lambda` is given.  Either `lambda` must be
#' given, or `u` AND `alpha`.
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
        obj <- function(k) {
            ratio <- besselI(k, 1L) / besselI(k, 0L)
            (ratio - u)^2
        }
        fit <- optim(par = 1, fn = obj, method = "L-BFGS-B",
                     lower = .Machine$double.eps)
        kappa_u <- fit$par
        T_val <- sqrt(kappa_u * u - log(besselI(kappa_u, 0L)))
        lambda <- -log(alpha) / T_val
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
    dens <- numeric(length(k))
    
    ## Case 1: k < 1e-4
    idx_small <- k < 1e-4
    if (any(idx_small)) {
        k_small <- k[idx_small]
        dens_small <- lambda/2 - lambda^2 * k_small/4 +
            (-9 * lambda/64 + lambda^3/16) * k_small^2 +
            (12 * lambda^2/128 - lambda^4/96) * k_small^3 +
            (1195 * lambda/36864 + 3 * lambda^3/512 + lambda^5/768) * k_small^4 +
            (-79 * lambda^2/6144 - lambda^6/7680) * k_small^5
        dens[idx_small] <- if (log) {
                               log(dens_small)
                           } else {
                               dens_small
                           }
    }
    
    ## Case 2: k >= 1e-4
    idx_large <- !idx_small
    if (any(idx_large)) {
        kk <- k[idx_large]
        I0.exp <- besselI(x = kk, nu = 0, expon.scaled = TRUE)
        I1.exp <- besselI(x = kk, nu = 1, expon.scaled = TRUE)
        I2.exp <- besselI(x = kk, nu = 2, expon.scaled = TRUE)

        idx_mid <- kk <= 1e+5
        idx_high <- !idx_mid
        
        dis.square <- numeric(length(kk))
        log.jac.partial <- numeric(length(kk))
        
        if (any(idx_mid)) {
            kk_mid <- kk[idx_mid]
            dis.square[idx_mid] <- kk_mid * I1.exp[idx_mid] / I0.exp[idx_mid] - kk_mid - log(I0.exp[idx_mid])
            log.jac.partial[idx_mid] <- log(kk_mid) + 
                log(( (I0.exp[idx_mid] + I2.exp[idx_mid]) / (2 * I0.exp[idx_mid]) - I1.exp[idx_mid]^2 / I0.exp[idx_mid]^2 ))
        }
        if (any(idx_high)) {
            kk_high <- kk[idx_high]
            dis.square[idx_high] <- -1/2 + kk_high - (kk_high - 1/2 * log(2*pi)) +
                log(kk_high)/2 - 1/(4 * kk_high) - 3/(16 * kk_high^2) -
                25/(96 * kk_high^3) - 65/(128 * kk_high^4) - 3219/(2560 * kk_high^5)
            log.jac.partial[idx_high] <- -log(2) - log(kk_high) + 1/(2 * kk_high) + 
                5/(8 * kk_high^2) + 59/(48 * kk_high^3) + 203/(64 * kk_high^4) + 12743/(1280 * kk_high^5)
        }
        
        distance <- sqrt(dis.square)
        log.dens.exp <- dexp(distance, rate = lambda, log = TRUE)
        log.jac <- -log(2) - log(distance) + log.jac.partial
        
        dens_large <- if (log) {
                          log.dens.exp + log.jac
                      } else {
                          exp(log.dens.exp) * exp(log.jac)
                      }
        dens[idx_large] <- dens_large
    }
    return(dens)
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
    I0.exp <- besselI(x=q,nu=0,expon.scaled = TRUE)
    I1.exp <- besselI(x=q,nu=1,expon.scaled = TRUE)
    dis.square <- ifelse(q <= 1e+5,
                         q*I1.exp/I0.exp-q-log(I0.exp),
                         -1/2 + q - (q - 1/2*log(2*pi)) + log(q)/2 - 1/(4*q) - 
                         3/(16*q^2) - 25/(96*q^3) - 65/(128*q^4) - 
                         3219/(2560*q^5))
    p <- pexp(sqrt(dis.square), rate = lambda, log.p = log)
    return(p)
}

