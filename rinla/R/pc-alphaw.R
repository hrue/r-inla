#' Utility functions for the PC prior for the `alpha` parameter in the
#' Weibull likelihood
#' 
#' Functions to evaluate, sample, compute quantiles and percentiles of the PC
#' prior for the `alpha` parameter in the Weibull likelihood
#' 
#' This gives the PC prior for the `alpha` parameter for the Weibull
#' likelihood, where `alpha=1` is the base model.
#' 
#' @aliases inla.pc.alphaw pc.alphaw pc.ralphaw inla.pc.ralphaw pc.dalphaw
#' inla.pc.dalphaw pc.palphaw inla.pc.palphaw pc.qalphaw inla.pc.qalphaw
#' @param n Number of observations
#' @param lambda The rate parameter in the PC-prior
#' @param alpha Vector of evaluation points, where `alpha>0`.
#' @param log Logical. Return the density in natural or log-scale.
#' @param p Vector of probabilities
#' @param q Vector of quantiles
#' @returns `inla.pc.dalphaw` gives the density, `inla.pc.palphaw`
#' gives the distribution function, `inla.pc.qalphaw` gives the quantile
#' function, and `inla.pc.ralphaw` generates random deviates.
#' @author Havard Rue \email{hrue@@r-inla.org}
#' @seealso inla.doc("pc.alphaw")
#' @examples
#' 
#'  x = inla.pc.ralphaw(100,  lambda = 5)
#'  d = inla.pc.dalphaw(x, lambda = 5)
#'  x = inla.pc.qalphaw(0.5, lambda = 5)
#'  inla.pc.palphaw(x, lambda = 5)
#'  
#' @name pc.alphaw
#' @rdname pc-alphaw
NULL

inla.pc.alphaw.cache <- function() {
    ## return the cache for these functions

    dist <- function(lalpha) {
        alpha <- exp(lalpha)
        gam <- 0.5772156649
        kld <- (gamma((1.0 + alpha) / alpha) * alpha + alpha * lalpha - alpha * gam + gam - alpha) / alpha
        return(sqrt(2.0 * kld))
    }

    tag <- "cache.pc.alphaw"
    if (!exists(tag, envir = inla.get.inlaEnv())) {
        lalphas <- seq(-4, 5, by = 0.005)
        lalphas[which.min(abs(lalphas))] <- 0 ## yes, make it exactly 0
        lalphas.pos <- lalphas[which(lalphas >= 0)]
        lalphas.neg <- lalphas[which(lalphas <= 0)]
        dist.pos <- dist(lalphas.pos)
        dist.neg <- dist(lalphas.neg)
        assign(tag,
            list(
                pos = list(
                    x = lalphas.pos,
                    dist = splinefun(lalphas.pos, dist.pos),
                    idist = splinefun(dist.pos, lalphas.pos)
                ),
                neg = list(
                    x = lalphas.neg,
                    dist = splinefun(lalphas.neg, dist.neg),
                    idist = splinefun(dist.neg, lalphas.neg)
                )
            ),
            envir = inla.get.inlaEnv()
        )
    }
    return(get(tag, envir = inla.get.inlaEnv()))
}

#' @rdname pc-alphaw
#' @aliases pc.ralphaw
#' @export
inla.pc.ralphaw <- function(n, lambda = 5) {
    cache <- inla.pc.alphaw.cache()
    x <- numeric(n)
    ind <- sample(c(-1, 1), n, replace = TRUE)
    idx.pos <- which(ind > 0)
    idx.neg <- which(ind < 0)
    z <- rexp(n, rate = lambda)
    if (length(idx.pos) > 0) {
        x[idx.pos] <- cache$pos$idist(z[idx.pos])
    }
    if (length(idx.neg) > 0) {
        x[idx.neg] <- cache$neg$idist(z[idx.neg])
    }
    return(exp(x))
}

#' @rdname pc-alphaw
#' @aliases pc.dalphaw
#' @export
inla.pc.dalphaw <- function(alpha, lambda = 5, log = FALSE) {
    cache <- inla.pc.alphaw.cache()
    d <- numeric(length(alpha))
    lalpha <- log(alpha)
    idx.pos <- which(lalpha >= 0)
    idx.neg <- which(lalpha < 0)
    if (length(idx.pos) > 0) {
        la <- lalpha[idx.pos]
        dist <- cache$pos$dist(la)
        deriv <- cache$pos$dist(la, deriv = 1)
        d[idx.pos] <- log(0.5) + log(lambda) - lambda * dist + log(abs(deriv)) - la
    }
    if (length(idx.neg) > 0) {
        la <- lalpha[idx.neg]
        dist <- cache$neg$dist(la)
        deriv <- cache$neg$dist(la, deriv = 1)
        d[idx.neg] <- log(0.5) + log(lambda) - lambda * dist + log(abs(deriv)) - la
    }
    return(if (log) d else exp(d))
}

#' @rdname pc-alphaw
#' @aliases pc.qalphaw
#' @export
inla.pc.qalphaw <- function(p, lambda = 5) {
    cache <- inla.pc.alphaw.cache()
    n <- length(p)
    q <- numeric(n)
    idx.pos <- which(p >= 0.5)
    idx.neg <- which(p < 0.5)
    if (length(idx.pos) > 0) {
        pp <- 2.0 * (p[idx.pos] - 0.5)
        qe <- qexp(pp, rate = lambda)
        q[idx.pos] <- cache$pos$idist(qe)
    }
    if (length(idx.neg) > 0) {
        pp <- 2.0 * (0.5 - p[idx.neg])
        qe <- qexp(pp, rate = lambda)
        q[idx.neg] <- cache$neg$idist(qe)
    }
    return(exp(q))
}

#' @rdname pc-alphaw
#' @aliases pc.palphaw
#' @export
inla.pc.palphaw <- function(q, lambda = 5) {
    cache <- inla.pc.alphaw.cache()
    n <- length(q)
    p <- numeric(n)
    idx.pos <- which(q >= 1.0)
    idx.neg <- which(q < 1.0)
    if (length(idx.pos) > 0) {
        qq <- cache$pos$dist(log(q[idx.pos]))
        pp <- pexp(qq, rate = lambda)
        p[idx.pos] <- 0.5 * (1 + pp)
    }
    if (length(idx.neg) > 0) {
        qq <- cache$neg$dist(log(q[idx.neg]))
        pp <- pexp(qq, rate = lambda)
        p[idx.neg] <- 0.5 * (1 - pp)
    }
    return(p)
}

inla.pc.alphaw.test <- function(lambda = 5) {
    ## this is just an internal test, and is not exported
    n <- 10^6
    x <- inla.pc.ralphaw(n, lambda = lambda)
    x <- x[x < quantile(x, prob = 0.999)]
    alpha <- seq(min(x), max(x), by = 0.001)
    hist(x, n = 300, prob = TRUE)
    lines(alpha, inla.pc.dalphaw(alpha, lambda), lwd = 3, col = "blue")
    print(sum(alpha * inla.pc.dalphaw(alpha, lambda) * diff(alpha)[1]))
    p <- seq(0.1, 0.9, by = 0.05)
    print(cbind(p = p, q.emp = as.numeric(quantile(x, prob = p)), q.true = inla.pc.qalphaw(p, lambda)))
    q <- as.numeric(quantile(x, prob = p))
    cdf <- ecdf(x)
    print(cbind(q = q, p.emp = cdf(q), p.true = inla.pc.palphaw(q, lambda)))
}

## inla.pc.alphaw.test()
