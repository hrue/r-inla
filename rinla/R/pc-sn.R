#' Utility functions for the PC prior for `skewness` in the skew-normal
#' linkfunction and likelihood
#' 
#' Functions to evaluate, sample, compute quantiles and percentiles of the PC
#' prior for the `skewness` in the skew-normal link-function and
#' likelihood
#' 
#' Defines the PC prior for the `skewness` for the skew-normal
#' linkfunction and likelihood, where `skew=0` is the base model. The
#' skewness range from -0.99527... to 0.99527....  ca.
#' 
#' @aliases inla.pc.sn pc.sn pc.rsn inla.pc.rsn pc.dsn inla.pc.dsn pc.psn
#' inla.pc.psn pc.qsn inla.pc.qsn
#' @param n number of observations
#' @param lambda the rate parameter in the PC prior
#' @param skew vector of evaluation points
#' @param log logical. return the density in natural or log-scale.
#' @param p vector of probabilities
#' @param q vector of quantiles
#' @returns `inla.pc.dsn` gives the density, `inla.pc.psn` gives
#' the distribution function, `inla.pc.qsn` gives the quantile function,
#' and `inla.pc.rsn` generates random deviates.
#' @author Havard Rue \email{hrue@@r-inla.org}
#' @seealso inla.doc("pc.sn")
#' @examples
#' 
#'  x = inla.pc.rsn(100,  lambda = 40)
#'  d = inla.pc.dsn(x, lambda = 40)
#'  x = inla.pc.qsn(0.5, lambda = 40)
#'  inla.pc.psn(x, lambda = 40)
#'  
#' @name pc.sn
#' @rdname pc-sn
NULL

inla.pc.sn.cache <- function(force = FALSE, write.files = FALSE) {
    sn.skew <- function(alpha) {
        delta <- alpha / sqrt(1 + alpha^2)
        skew <- (4 - pi) / 2 * (delta * sqrt(2 / pi))^3 / (1 - 2 * delta^2 / pi)^(3 / 2)
        return(skew)
    }

    sn <- function(x, alpha) {
        delta <- alpha / sqrt(1 + alpha^2)
        omega <- 1 / sqrt(1 - 2 * delta^2 / pi)
        xi <- -omega * delta * sqrt(2 / pi)
        xx <- (x - xi) / omega
        return(2 / omega * dnorm(xx) * pnorm(alpha * xx))
    }

    kld.value <- function(x, alpha) {
        dsn <- sn(x, alpha)
        idx <- which(dsn > 0.0)
        dsn <- dsn[idx]
        x <- x[idx]
        if (length(x) %% 2 == 0) {
            x <- x[-length(x)]
            dsn <- dsn[-length(dsn)]
        }
        w <- c(1, rep(c(4, 2), (length(x) - 3) %/% 2), 4, 1) / 3
        dn <- dnorm(x)
        dsn <- dsn / sum(w * dsn)
        dn <- dn / sum(w * dn)
        return(sum(w * dsn * log(dsn / dn)))
    }

    kld.t <- function(alpha) {
        t1 <- alpha * alpha
        t8 <- t1 * t1
        t10 <- t8 * t1 * (0.3960827996046746E-2 + t1 * (-0.4043436622715706E-2 + t1 * (
            0.205737163227263E-2 + 0.2400653839003813E-4 * t1)))
        return(t10)
    }

    kld <- function(alpha) {
        ## 0.4 is decided empirically, so I do not see any discontinuities
        if (alpha < 0.4) {
            return(kld.t(alpha))
        } else {
            dx <- 0.005
            ran <- 8
            x <- seq(-ran, ran, by = dx)
            return(kld.value(x, alpha))
        }
    }

    dist <- function(alpha) {
        d <- numeric(length(alpha))
        for (i in seq_along(alpha)) {
            d[i] <- sqrt(2.0 * max(0, kld(alpha[i])))
        }
        idx <- which(!(is.na(d) | is.infinite(d)))
        return(data.frame(alpha = alpha[idx], dist = d[idx]))
    }

    tag <- "cache.pc.sn"
    if (force || !exists(tag, envir = inla.get.inlaEnv())) {
        ran <- 25
        alphas <- seq(0, ran, length.out = 500)[-1]
        alphas <- c(-rev(alphas), 0, alphas)
        alphas.pos <- alphas[which(alphas >= 0)]
        alphas.neg <- alphas[which(alphas <= 0)]
        dist.pos <- dist(alphas.pos)$dist
        dist.neg <- rev(dist.pos)

        skews <- sn.skew(alphas)
        skews.pos <- sn.skew(alphas.pos)
        skews.neg <- sn.skew(alphas.neg)

        if (write.files) {
            print("write file [s-sn-pc-prior.dat]")
            write(skews.pos, ncolumns = 1, file = "s-sn-pc-prior.dat")
            print("write file [d-sn-pc-prior.dat]")
            write(dist.pos, ncolumns = 1, file = "d-sn-pc-prior.dat")
        }

        ## we could improve here, as d'(s) = d*'(s*) d(s) / s, where d* = log(d), and s* =
        ## log(s), but since we're using the Taylor expansion near zero, it is ok not to.

        m <- "monoH.FC"
        assign(tag,
            list(
                skews.range = range(skews),
                dist.max = max(dist.pos),
                pos = list(
                    x = skews.pos,
                    a = splinefun(skews.pos, alphas.pos, method = m),
                    dist = splinefun(skews.pos, dist.pos, method = m),
                    idist = splinefun(dist.pos, skews.pos, method = m)
                ),
                neg = list(
                    x = skews.neg,
                    a = splinefun(skews.neg, alphas.neg, method = m),
                    dist = splinefun(skews.neg, dist.neg, method = m),
                    idist = splinefun(dist.neg, skews.neg, method = m)
                )
            ),
            envir = inla.get.inlaEnv()
        )
    }
    return(get(tag, envir = inla.get.inlaEnv()))
}

#' @rdname pc-sn
#' @export
inla.pc.rsn <- function(n, lambda = 40) {
    cache <- inla.pc.sn.cache()
    fac <- 1.0 - exp(-lambda * cache$dist.max)
    x <- numeric(n)
    ind <- sample(c(-1, 1), n, replace = TRUE)
    idx.pos <- which(ind > 0)
    idx.neg <- which(ind < 0)
    z <- -1.0 / lambda * log(1.0 - runif(n) * fac)
    if (length(idx.pos) > 0) {
        x[idx.pos] <- cache$pos$idist(z[idx.pos])
    }
    if (length(idx.neg) > 0) {
        x[idx.neg] <- cache$neg$idist(z[idx.neg])
    }
    return(x)
}

#' @rdname pc-sn
#' @export
inla.pc.dsn <- function(skew, lambda = 40, log = FALSE) {
    cache <- inla.pc.sn.cache()
    fac <- 1.0 - exp(-lambda * cache$dist.max)
    d <- numeric(length(skew))
    idx.pos <- which(skew >= 0)
    idx.neg <- which(skew < 0)
    if (length(idx.pos) > 0) {
        la <- skew[idx.pos]
        dist <- cache$pos$dist(la)
        deriv <- cache$pos$dist(la, deriv = 1)
        d[idx.pos] <- log(0.5) + log(lambda) - lambda * dist + log(abs(deriv)) - log(fac)
    }
    if (length(idx.neg) > 0) {
        la <- skew[idx.neg]
        dist <- cache$neg$dist(la)
        deriv <- cache$neg$dist(la, deriv = 1)
        d[idx.neg] <- log(0.5) + log(lambda) - lambda * dist + log(abs(deriv)) - log(fac)
    }
    return(if (log) d else exp(d))
}

#' @rdname pc-sn
#' @export
inla.pc.qsn <- function(p, lambda = 40) {
    cache <- inla.pc.sn.cache()
    fac <- 1.0 - exp(-lambda * cache$dist.max)
    n <- length(p)
    q <- numeric(n)
    idx.pos <- which(p >= 0.5)
    idx.neg <- which(p < 0.5)
    if (length(idx.pos) > 0) {
        pp <- 2.0 * (p[idx.pos] - 0.5) * fac
        qe <- qexp(pp, rate = lambda)
        q[idx.pos] <- cache$pos$idist(qe)
    }
    if (length(idx.neg) > 0) {
        pp <- 2.0 * (0.5 - p[idx.neg]) * fac
        qe <- qexp(pp, rate = lambda)
        q[idx.neg] <- cache$neg$idist(qe)
    }
    return(q)
}

#' @rdname pc-sn
#' @export
inla.pc.psn <- function(q, lambda = 40) {
    cache <- inla.pc.sn.cache()
    fac <- 1.0 - exp(-lambda * cache$dist.max)
    n <- length(q)
    p <- numeric(n)
    idx.pos <- which(q >= 0.0)
    idx.neg <- which(q < 0.0)
    if (length(idx.pos) > 0) {
        qq <- cache$pos$dist(q[idx.pos])
        pp <- pexp(qq, rate = lambda) / fac
        p[idx.pos] <- (1.0 + pp) / 2.0
    }
    if (length(idx.neg) > 0) {
        qq <- cache$neg$dist(q[idx.neg])
        pp <- pexp(qq, rate = lambda) / fac
        p[idx.neg] <- (1.0 - pp) / 2.0
    }
    return(p)
}

inla.pc.sn.test <- function(lambda = 40, n = 10^6) {
    ## this is for testing only
    f <- function(...) format(..., digits = 3)
    x <- inla.pc.rsn(n, lambda = lambda)
    skew <- seq(min(x), max(x), by = 0.01)
    hist(x, n = 3000, prob = TRUE, xlim = quantile(x, c(0.001, 0.999)))
    lines(skew, inla.pc.dsn(skew, lambda), lwd = 3, col = "blue")
    p <- c(0.001, 0.01, seq(0.05, 0.95, by = 0.1), 0.99, 0.999)
    print(cbind(
        p = f(p), q.emp = f(as.numeric(quantile(x, prob = p))),
        q.true = f(inla.pc.qsn(p, lambda))
    ))
    q <- as.numeric(quantile(x, prob = p))
    cdf <- ecdf(x)
    print(cbind(q = f(q), p.emp = f(cdf(q)), p.true = f(inla.pc.psn(q, lambda))))
}

inla.pc.sn.test2 <- function() {
    fun <- function(intern.intercept, skew) {
        cache <- inla.pc.sn.cache()
        alpha <- sign(skew) * cache$pos$a(abs(skew))
        delta <- alpha / sqrt(1 + alpha^2)
        omega <- 1 / sqrt(1 - 2 * delta^2 / pi)
        xi <- -omega * delta * sqrt(2 / pi)
        print(c(skew = skew, xi = xi, omega = omega, alpha = alpha))
        return(sn::qsn(intern.intercept, xi = xi, omega = omega, alpha = alpha))
    }

    print(fun(0.43, 0.123))
    print(fun(0.823, -0.123))
}

if (FALSE) {
    inla.get.inlaEnv <- function(...) inla.get.inlaEnv(...)
}
