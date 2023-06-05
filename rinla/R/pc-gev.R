#' Utility functions for the PC prior for the \code{tail} parameter in the GEV
#' likelihood
#' 
#' Functions to evaluate, sample, compute quantiles and percentiles of the PC
#' prior for the \code{tail} parameter in the GEV likelihood
#' 
#' This gives the PC prior for the \code{tail} parameter for the GEV
#' likelihood, where \code{xi=0} is the base model.
#' 
#' @aliases inla.pc.gevtail pc.gevtail pc.rgevtail inla.pc.rgevtail pc.dgevtail
#' inla.pc.dgevtail pc.pgevtail inla.pc.pgevtail pc.qgevtail inla.pc.qgevtail
#' @param n Number of observations
#' @param lambda The rate parameter in the PC-prior
#' @param xi Vector of evaluation points, where \code{1>xi>0}.
#' @param log Logical. Return the density in natural or log-scale.
#' @param p Vector of probabilities
#' @param q Vector of quantiles
#' @returns \code{inla.pc.dgevtail} gives the density,
#' \code{inla.pc.pgevtail} gives the distribution function,
#' \code{inla.pc.qgevtail} gives the quantile function, and
#' \code{inla.pc.rgevtail} generates random deviates.
#' @author Havard Rue \email{hrue@@r-inla.org}
#' @seealso inla.doc("pc.gevtail")
#' @examples
#' 
#'  xi = inla.pc.rgevtail(100,  lambda = 7)
#'  d = inla.pc.dgevtail(xi, lambda = 7)
#'  xi = inla.pc.qgevtail(0.5, lambda = 7)
#'  inla.pc.pgevtail(xi, lambda = 7)
#'  
#' @name pc.gevtail
#' @rdname pc-gev
NULL

inla.pc.gevtail.map <- function(x, inverse = FALSE, deriv = 0) {
    ## map from R to (0, 1)
    if (!inverse) {
        ex <- exp(x)
        if (deriv == 0) {
            return(ex / (1 + ex))
        } else if (deriv == 1) {
            return(ex / (1 + ex)^2)
        } else {
            stop("error")
        }
    } else {
        if (deriv == 0) {
            return(log(x / (1.0 - x)))
        } else {
            stop("error")
        }
    }
}
inla.pc.gevtail.cache <- function() {
    ## return the cache for these functions. even though we could do all this analytic, the
    ## simplicity of just reusing the code was to tempting...

    dist <- function(lxi) {
        xi <- inla.pc.gevtail.map(lxi)
        kld <- xi^2 / (1 - xi)
        return(sqrt(2.0 * kld))
    }

    tag <- "cache.pc.gevtail"
    if (!exists(tag, envir = inla.get.inlaEnv())) {
        lxis <- seq(-15, 10, by = 0.01)
        dist <- dist(lxis)
        assign(tag,
            list(
                x = lxis,
                dist = splinefun(lxis, dist),
                idist = splinefun(dist, lxis)
            ),
            envir = inla.get.inlaEnv()
        )
    }
    return(get(tag, envir = inla.get.inlaEnv()))
}

#' @rdname pc-gev
#' @export
inla.pc.rgevtail <- function(n, lambda = 7) {
    cache <- inla.pc.gevtail.cache()
    x <- numeric(n)
    z <- rexp(n, rate = lambda)
    x <- cache$idist(z)
    return(inla.pc.gevtail.map(x))
}

#' @rdname pc-gev
#' @export
inla.pc.dgevtail <- function(xi, lambda = 7, log = FALSE) {
    cache <- inla.pc.gevtail.cache()
    d <- numeric(length(xi))
    lxi <- inla.pc.gevtail.map(xi, inverse = TRUE)
    dist <- cache$dist(lxi)
    deriv <- cache$dist(lxi, deriv = 1)
    d <- log(lambda) - lambda * dist + log(abs(deriv)) - log(abs(inla.pc.gevtail.map(lxi, deriv = 1)))
    return(if (log) d else exp(d))
}

#' @rdname pc-gev
#' @export
inla.pc.qgevtail <- function(p, lambda = 7) {
    cache <- inla.pc.gevtail.cache()
    n <- length(p)
    q <- numeric(n)
    qe <- qexp(p, rate = lambda)
    q <- cache$idist(qe)
    return(inla.pc.gevtail.map(q))
}

#' @rdname pc-gev
#' @export
inla.pc.pgevtail <- function(q, lambda = 7) {
    cache <- inla.pc.gevtail.cache()
    n <- length(q)
    p <- numeric(n)
    qq <- cache$dist(inla.pc.gevtail.map(q, inverse = TRUE))
    p <- pexp(qq, rate = lambda)
    return(p)
}

inla.pc.gevtail.test <- function(lambda = 7) {
    ## this is just an internal test, and is not exported
    n <- 10^6
    x <- inla.pc.rgevtail(n, lambda = lambda)
    x <- x[x < quantile(x, prob = 0.999)]
    xi <- seq(min(x), max(x), by = 0.001)
    hist(x, n = 300, prob = TRUE)
    lines(xi, inla.pc.dgevtail(xi, lambda), lwd = 3, col = "blue")
    p <- seq(0.1, 0.9, by = 0.05)
    print(cbind(p = p, q.emp = as.numeric(quantile(x, prob = p)), q.true = inla.pc.qgevtail(p, lambda)))
    q <- as.numeric(quantile(x, prob = p))
    cdf <- ecdf(x)
    print(cbind(q = q, p.emp = cdf(q), p.true = inla.pc.pgevtail(q, lambda)))
}

## How to run test outside the library:
## inla.get.inlaEnv = function(...) inla.get.inlaEnv(...)
## inla.pc.gevtail.test(lambda = .01)
