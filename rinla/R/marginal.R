#' @title Functions which operates on marginals
#' 
#' @description
#' Density, distribution function, quantile function, random generation,
#' hpd-interval, interpolation, expectations, mode and transformations of
#' marginals obtained by `inla` or `inla.hyperpar()`.  These
#' functions computes the density (inla.dmarginal), the distribution function
#' (inla.pmarginal), the quantile function (inla.qmarginal), random generation
#' (inla.rmarginal), spline smoothing (inla.smarginal), computes expected
#' values (inla.emarginal), computes the mode (inla.mmarginal), transforms the
#' marginal (inla.tmarginal), and provide summary statistics (inla.zmarginal).
#' 
#' @rdname marginal
#' @name marginal
#' @aliases inla.marginal marginal pmarginal inla.pmarginal qmarginal
#' inla.qmarginal dmarginal inla.dmarginal rmarginal inla.rmarginal
#' inla.hpdmarginal hpdmarginal inla.emarginal emarginal inla.smarginal
#' smarginal inla.tmarginal inla.mmarginal mmarginal inla.zmarginal zmarginal
#' @param marginal A marginal object from either `inla` or
#' `inla.hyperpar()`, which is either `list(x=c(), y=c())` with
#' density values `y` at locations `x`, or a `matrix(,n,2)` for
#' which the density values are the second column and the locations in the
#' first column.  The`inla.hpdmarginal()`-function assumes a unimodal
#' density.
#' @param fun A (vectorised) function like `function(x) exp(x)` to compute
#' the expectation against, or which define the transformation new = fun(old)
#' @param x Evaluation points
#' @param q Quantiles
#' @param p Probabilities
#' @param n The number of observations. If `length(n) > 1`, the length is
#' taken to be the number required.
#' @param h.diff The step-length for the numerical differeniation inside
#' `inla.tmarginal`
#' @param ... Further arguments to be passed to function which expectation is
#' to be computed.
#' @param log Return density or interpolated density in log-scale?
#' @param normalize Renormalise the density after interpolation?
#' @param len Number of locations used to interpolate the distribution
#' function.
#' @param keep.type If `FALSE` then return a `list(x=, y=)`,
#' otherwise if `TRUE`, then return a matrix if the input is a matrix
#' @param extrapolate How much to extrapolate on each side when computing the
#' interpolation. In fraction of the range.
#' @param factor The number of points after interpolation is `factor`
#' times the original number of points; which is argument `n` in
#' `spline`
#' @param method Which method should be used to layout points for where the
#' transformation is computed.
#' @param silent Output the result visually (TRUE) or just through the call.
#' @returns `inla.smarginal` returns `list=c(x=c(), y=c())` of
#' interpolated values do extrapolation using the factor given, and the
#' remaining function returns what they say they should do.
#' @author Havard Rue \email{hrue@@r-inla.org}
#' @seealso [inla()], [inla.hyperpar()]
#' @examples
#' 
#'  ## a simple linear regression example
#'  n = 10
#'  x = rnorm(n)
#'  sd = 0.1
#'  y = 1+x + rnorm(n,sd=sd)
#'  res = inla(y ~ 1 + x, data = data.frame(x,y),
#'             control.family=list(initial = log(1/sd^2L),fixed=TRUE))
#' 
#'  ## chose a marginal and compare the with the results computed by the
#'  ## inla-program
#'  r = res$summary.fixed["x",]
#'  m = res$marginals.fixed$x
#' 
#'  ## compute the 95% HPD interval
#'  inla.hpdmarginal(0.95, m)
#' 
#'  x = seq(-6, 6, length.out = 1000)
#'  y = dnorm(x)
#'  inla.hpdmarginal(0.95, list(x=x, y=y))
#' 
#'  ## compute the the density for exp(r), version 1
#'  r.exp = inla.tmarginal(exp, m)
#'  ## or version 2
#'  r.exp = inla.tmarginal(function(x) exp(x), m)
#' 
#'  ## to plot the marginal, we use the inla.smarginal, which interpolates (in
#'  ## log-scale). Compare with some samples.
#'  plot(inla.smarginal(m), type="l")
#'  s = inla.rmarginal(1000, m)
#'  hist(inla.rmarginal(1000, m), add=TRUE, prob=TRUE)
#'  lines(density(s), lty=2)
#' 
#'  m1 = inla.emarginal(function(x) x, m)
#'  m2 = inla.emarginal(function(x) x^2L, m)
#'  stdev = sqrt(m2 - m1^2L)
#'  q = inla.qmarginal(c(0.025,0.975), m)
#' 
#'  ## inla-program results
#'  print(r)
#' 
#'  ## inla.marginal-results (they shouldn't be perfect!)
#'  print(c(mean=m1, sd=stdev, "0.025quant" = q[1], "0.975quant" = q[2L]))
#'  ## using the buildt-in function
#'  inla.zmarginal(m)
#'  
NULL



### functions to work with the marginal, either defined as a matrix
### x[2, n], or a list(x=, y=).  NOTE:: there are NO EXTRAPOLATION, so
### the marginal is defined FROM...TO!

`inla.marginal.fix` <- function(marginal) {
    ## just remove points where the density is <= 0 or where the
    ## density is to small compared to the maximum density. (othewise
    ## we can get trouble with the spline interpolation). same with
    ## 'x'? No...
    eps <- .Machine[["double.eps"]] * 1000.0
    ## marginal = spline(marginal)
    if (is.matrix(marginal)) {
        if (any(is.na(marginal[, 2L]))) {
            idx <- which(is.na(marginal[, 2L]))


            marginal <- marginal[-idx, ]
        }
        i <- (marginal[, 2L] > 0.0) & (abs(marginal[, 2L] / max(marginal[, 2L])) > eps)
        m <- list(x = marginal[i, 1L], y = marginal[i, 2L])
        ## i = c(diff(marginal[, 1L]) > eps, TRUE)
        ## m = list(x=marginal[i, 1L], y=marginal[i, 2L])
    } else {
        if (any(is.na(marginal[["y"]]))) {
            idx <- which(is.na(marginal[["y"]]))
            marginal[["x"]] <- marginal[["x"]][-idx]
            marginal[["y"]] <- marginal[["y"]][-idx]
        }
        i <- (marginal[["y"]] > 0.0) & (abs(marginal[["y"]] / max(marginal[["y"]])) > eps)
        m <- list(x = marginal[["x"]][i], y = marginal[["y"]][i])
        ## i = c(diff(marginal[["x"]]) > eps, TRUE)
        ## m = list(x = marginal[["x"]][i], y = marginal[["y"]][i])
    }

    return(m)
}

`inla.spline` <- function(x, ...) {
    s <- spline(x, ...)
    if (is.matrix(x)) {
        m <- cbind(x = s[["x"]], y = s[["y"]])
    } else if (is.list(x)) {
        m <- list(x = s[["x"]], y = s[["y"]])
    } else {
        m <- s
    }
    return(m)
}

#' @rdname marginal
#' @export
`inla.smarginal` <- function(marginal, log = FALSE, extrapolate = 0.0, keep.type = FALSE, factor = 15L) {
    ## for marginal in matrix MARGINAL, which is a marginal density,
    ## return the nice interpolated (x, y) where the interpolation is
    ## done in log(y)

    is.mat <- is.matrix(marginal)
    m <- inla.marginal.fix(marginal)
    r <- diff(range(m[["x"]]))
    xmin <- min(m[["x"]]) - extrapolate * r
    xmax <- max(m[["x"]]) + extrapolate * r
    n <- factor * length(m[["x"]])
    xx <- seq(xmin, xmax, length.out = n)
    if (extrapolate) {
        xx <- c(xmin, m[["x"]], xmax)
    } else {
        xx <- m[["x"]]
    }
    nx <- length(xx)
    dx <- nx * diff(xx) / mean(diff(xx))
    xnew <- c(0.0, cumsum(sqrt(dx)))
    xnew <- xmin + (xmax - xmin) * ((xnew - min(xnew)) / (max(xnew) - min(xnew)))
    fun <- splinefun(xnew, xx, method = "hyman")
    fun.inv <- splinefun(xx, xnew, method = "hyman")
    ans <- spline(fun.inv(m[["x"]]), log(m[["y"]]), xmin = fun.inv(xmin), xmax = fun.inv(xmax), n = n, method = "fmm")
    ans[["x"]] <- fun(ans[["x"]])
    if (!log) {
        ans[["y"]] <- exp(ans[["y"]])
        ans <- inla.marginal.fix(ans)
    }
    if (is.mat && keep.type) {
        return(cbind(x = ans[["x"]], y = ans[["y"]]))
    } else {
        return(ans)
    }
}

`inla.sfmarginal` <- function(marginal) {
    ## for marginal in matrix MARGINAL, which is a marginal density,
    ## return the spline-function which returns the log-density

    m <- inla.marginal.fix(marginal)
    r <- range(m[["x"]])
    return(list(range = r, fun = splinefun(m[["x"]], log(m[["y"]]))))
}

#' @rdname marginal
#' @export
`inla.emarginal` <- function(fun, marginal, ...) {
    ## compute E(FUN(x)), where the marginal of x is given in
    ## `marginal'; see inla.smarginal(). Also work for FUN(x)
    ## returning a vector.

    xx <- inla.smarginal(marginal)
    n <- length(xx[["x"]])
    if (n %% 2L == 0L) {
        n <- n - 1L
        xx[["x"]] <- xx[["x"]][1L:n]
        xx[["y"]] <- xx[["y"]][1L:n]
    }

    ## use Simpsons integration rule
    i.0 <- c(1L, n)
    i.4 <- seq(2L, n - 1L, by = 2L)
    i.2 <- seq(3L, n - 2L, by = 2L)

    dx <- diff(xx[["x"]])
    dx <- 0.5 * (c(dx, 0.0) + c(0.0, dx))
    fun <- match.fun(fun)
    ff <- fun(xx[["x"]][1L:n], ...) * xx[["y"]][1L:n] * dx
    nf <- length(ff) %/% n
    e <- numeric(nf)
    off <- 0L
    for (i in 1L:nf) {
        e[i] <- sum(sum(ff[i.0 + off]) + 4.0 * sum(ff[i.4 + off]) + 2.0 * sum(ff[i.2 + off]))
        off <- off + n
    }

    ## normalise, so that E(1) = 1
    ff <- dx * xx[["y"]][1L:n]
    e.1 <- sum(sum(ff[i.0]) + 4.0 * sum(ff[i.4]) + 2.0 * sum(ff[i.2]))

    return(e / e.1)
}

#' @rdname marginal
#' @export
`inla.dmarginal` <- function(x, marginal, log = FALSE) {
    ## return the density of the marginal. if LOG, in log-scale.  this
    ## density is not renormalized (but should be already normalized
    ## in some scale).

    ## return list(range=, fun=..)
    f <- inla.sfmarginal(inla.smarginal(marginal))

    n <- length(x)
    d <- numeric(n)

    for (i in 1L:n) {
        if (x[i] >= f[["range"]][[1L]] && x[i] <= f[["range"]][[2L]]) {
            if (log) {
                d[i] <- f[["fun"]](x[i])
            } else {
                d[i] <- exp(f[["fun"]](x[i]))
            }
        } else {
            if (log) {
                d[i] <- -Inf
            } else {
                d[i] <- 0.0
            }
        }
    }
    return(d)
}

#' @rdname marginal
#' @export
`inla.pmarginal` <- function(q, marginal, normalize = TRUE, len = 2048L) {
    f <- inla.sfmarginal(inla.smarginal(marginal))
    xx <- seq(f[["range"]][[1L]], f[["range"]][[2L]], length.out = len)
    d <- cumsum(exp(f[["fun"]](xx)))
    d <- d / d[length(d)]

    ## just spline-interpolate the mapping
    fq <- splinefun(xx, d, method = "hyman")

    ## just make sure the p's are in [0, 1]
    n <- length(q)
    xx <- pmin(pmax(q, rep(f[["range"]][[1L]], n)), rep(f[["range"]][[2L]], n))

    return(fq(xx))
}

#' @rdname marginal
#' @export
`inla.qmarginal` <- function(p, marginal, len = 2048L) {
    f <- inla.sfmarginal(inla.smarginal(marginal))
    xx <- seq(f[["range"]][[1L]], f[["range"]][[2L]], length.out = len)
    d <- cumsum(exp(f[["fun"]](xx)))
    d <- d / d[length(d)]

    ## for the moment, we remove only duplicated zero's and one's.
    eps <- .Machine[["double.eps"]] * 1000.0
    for (val in c(0.0, 1.0)) {
        is.val <- which(abs(d - val) <= eps)
        if (length(is.val) > 1L) {
            ## keep the first, ie remove it from the list to the removed
            is.val <- is.val[-1L]
            ## and remove the rest
            d <- d[-is.val]
            xx <- xx[-is.val]
        }
    }

    ## just spline-interpolate the inverse mapping
    fq <- splinefun(d, xx, method = "hyman")

    ## just make sure the p's are in [0, 1]
    n <- length(p)
    pp <- pmin(pmax(p, rep(0.0, n)), rep(1.0, n))

    return(fq(pp))
}
#' @rdname marginal
#' @export
`inla.hpdmarginal` <- function(p, marginal, len = 2048L) {
    sm <- inla.smarginal(marginal, keep.type = FALSE)
    f <- inla.sfmarginal(sm)
    x.range <- f[["range"]]
    xx <- seq(x.range[[1L]], x.range[[2L]], length.out = len)
    ## simply add 0 density outside the range.
    dx <- diff(xx)[[1L]]
    ef <- exp(f[["fun"]](xx))
    xx <- c(min(xx) - dx, xx, max(xx) + dx)
    d <- c(0.0, ef, 0.0)
    d <- cumsum(d)
    d <- d / d[length(d)]

    ## for the moment, we remove only duplicated zero's and one's.
    eps <- .Machine[["double.eps"]] * 1000.0
    for (val in c(0.0, 1.0)) {
        is.val <- which(abs(d - val) <= eps)
        if (length(is.val) > 1L) {
            ## keep the first, ie remove it from the list to the removed
            is.val <- is.val[-1L]
            ## and remove the rest
            d <- d[-is.val]
            xx <- xx[-is.val]
        }
    }
    ## just spline-interpolate the inverse mapping
    fq <- splinefun(d, xx, method = "hyman")

    ## just make sure the p's are in [0, 1]
    np <- length(p)
    pp <- 1.0 - pmin(pmax(p, rep(0.0, np)), rep(1.0, np))

    ## parts of this code below is taken from TeachingDemo::hpd
    f <- function(x, posterior.icdf, conf) {
        return(posterior.icdf(1.0 - conf + x) - posterior.icdf(x))
    }

    tol <- sqrt(.Machine[["double.eps"]])
    tol <- 1E-6
    result <- matrix(NA, np, 2L)
    for (i in 1L:np) {
        out <- optimize(f, c(0.0, pp[i]), posterior.icdf = fq, conf = pp[i], tol = tol)
        result[i, ] <- c(fq(out[["minimum"]]), fq(1.0 - pp[i] + out[["minimum"]]))
        ## make sure its in the range of the marginal itself as we expanded the domain
        result[i, ] <- c(
            min(x.range[[2L]], max(result[i, 1L], x.range[[1L]])),
            min(x.range[[2L]], max(result[i, 2L], x.range[[1L]]))
        )
    }
    colnames(result) <- c("low", "high")
    rownames(result) <- paste("level:", format(1.0 - pp, digits = 6L, justify = "left", trim = TRUE), sep = "")

    return(result)
}
#' @rdname marginal
#' @export
`inla.rmarginal` <- function(n, marginal) {
    return(inla.qmarginal(runif(n), marginal))
}

`inla.marginal.transform` <- function(fun, marginal, n = 2048L, h.diff =
                                          .Machine[["double.eps"]]^(1.0 / 3.0),
                                      method = c("quantile", "linear")) {
    return(inla.tmarginal(fun, marginal, n, h.diff, method = method))
}

`inla.deriv.func` <- function(fun,
                              step.size = .Machine[["double.eps"]]^(1.0 / 4.0)) {
    ## return a function computing the derivatives to 'fun = function(x)' which must
    ## vectorize
    func <- match.fun(fun)
    if (inla.require("Deriv")) {
        fd <- try(Deriv::Deriv(func, names(formals(func))[[1L]]), silent = TRUE)
        if (!inherits(fd, "try-error")) {
            return(fd)
        }
    }
    if (inla.require("numDeriv")) {
        fd <- function(x) unlist(lapply(x, function(xi) numDeriv::grad(func, xi)))
    } else {
        fd <- function(x) {
            ((-func(x + 2.0 * step.size) + 8.0 * func(x + step.size) -
                8.0 * func(x - step.size) + func(x - 2.0 * step.size)) / (12.0 * step.size))
        }
    }

    return(fd)
}

#' @rdname marginal
#' @export
`inla.tmarginal` <- function(fun, marginal, n = 2048L,
                             h.diff = .Machine[["double.eps"]]^(1.0 / 3.0),
                             method = c("quantile", "linear")) {
    ff <- match.fun(fun)

    is.mat <- is.matrix(marginal)
    m <- inla.smarginal(marginal)
    r <- range(m[["x"]])

    method <- match.arg(method)
    if (inla.strcasecmp(method, "quantile")) {
        x <- inla.qmarginal((1L:n) / (n + 1.0), marginal)
    } else if (inla.strcasecmp(method, "linear")) {
        x <- seq(r[[1L]], r[[2L]], length.out = n)
    } else {
        stop("unknown method")
    }
    xx <- ff(x)

    fd <- inla.deriv.func(ff)
    log.dens <- inla.dmarginal(x, marginal, log = FALSE) / abs(fd(x))

    ## if decreasing function, then reverse
    if (xx[[1L]] > xx[n]) {
        xx[1L:n] <- xx[n:1L]
        log.dens[1L:n] <- log.dens[n:1L]
    }

    if (is.mat) {
        ret <- cbind(x = xx, y = log.dens)
    } else {
        ret <- list(x = xx, y = log.dens)
    }

    if (FALSE) {
        class(ret) <- "inla.marginal"
        attr(ret, "inla.tag") <- paste(attr(marginal, "inla.tag"), "transformed")
    }

    return(ret)
}

#' @rdname marginal
#' @export
`inla.mmarginal` <- function(marginal) {
    p <- seq(0.01, 0.99, by = 0.01)
    n <- length(p)
    x <- inla.qmarginal(p, marginal)
    d <- inla.dmarginal(x, marginal, log = TRUE)
    idx <- which.max(d)
    res <- optimize(
        f = inla.dmarginal,
        interval = c(x[max(1L, idx - 1L)], x[min(n, idx + 1L)]),
        maximum = TRUE,
        ## arguments to inla.dmarginal
        marginal = marginal, log = TRUE
    )
    return(res[["maximum"]])
}

#' @rdname marginal
#' @export
`inla.zmarginal` <- function(marginal, silent = FALSE) {
    stopifnot(inla.is.marginal(marginal))

    m <- inla.emarginal(function(xx) c(xx, xx^2L), marginal)
    q <- inla.qmarginal(c(0.025, 0.25, 0.5, 0.75, 0.975), marginal)
    s <- sqrt(max(0.0, m[[2L]] - m[[1L]]^2L))

    if (!silent) {
        ## cat(as.character(match.call()[-1L]), "\n")
        ## cat("+--------+-----+------------\n")
        cat("Mean           ", format(m[[1L]], digits = 6L), "\n")
        cat("Stdev          ", format(s, digits = 6L), "\n")
        cat("Quantile  0.025", format(q[[1L]], digits = 6L), "\n")
        cat("Quantile  0.25 ", format(q[[2L]], digits = 6L), "\n")
        cat("Quantile  0.5  ", format(q[[3L]], digits = 6L), "\n")
        cat("Quantile  0.75 ", format(q[[4L]], digits = 6L), "\n")
        cat("Quantile  0.975", format(q[[5L]], digits = 6L), "\n")
    }
    return(invisible(list(
        mean = m[[1L]], sd = s, "quant0.025" = q[[1L]],
        "quant0.25" = q[[2L]], "quant0.5" = q[[3L]], "quant0.75" = q[[4L]], "quant0.975" = q[[5L]]
    )))
}

#' @rdname marginal
#' @export
`inla.is.marginal` <- function(marginal) {
    return((is.matrix(marginal) && ncol(marginal) == 2L && nrow(marginal) > 2L) ||
        (is.list(marginal) && all(names(marginal) == c("x", "y"))))
}


## NOT-YET-IN-USE
## 'plot' and 'summary'-methods for marginals
#' Not yet in use
#' 
#' Method not yet in use
#' 
#' @param x,y,... Not yet in use
#' 
#' @keywords internal
#' @exportS3Method plot inla.marginal
plot.inla.marginal <- function(x, y, ...) {
    m <- inla.emarginal(function(xx) c(xx, xx^2L), x)
    xlab <- paste(
        "x (mean", format(m[[1L]], digits = 4L),
        ", stdev", format(sqrt(max(0.0, m[[2L]] - m[[1L]]^2L)), digits = 4L), ")"
    )
    plot(inla.smarginal(x, extrapolate = 0.0), type = "l", xlab = xlab, ylab = "marginal density", ...)

    return(invisible())
}
## NOT-YET-IN-USE
#' Not yet in use
#' 
#' Method not yet in use
#' 
#' @usage
#' \method{summary}{inla.marginal}(object, ...)
#' 
#' @param object,... Not yet in use
#' 
#' @keywords internal
#' @exportS3Method summary inla.marginal
summary.inla.marginal <- function(object, ...) {
    ## option
    silent <- FALSE
    inla.eval.dots(...)

    m <- inla.emarginal(function(xx) c(xx, xx^2L), object)
    q <- inla.qmarginal(c(0.025, 0.25, 0.5, 0.75, 0.975), object)
    s <- sqrt(max(0.0, m[[2L]] - m[[1L]]^2L))

    if (!silent) {
        cat("Properties of: ", "FIXME", "\n")
        cat("+--------------+------------------\n")
        cat("Mean           ", format(m[[1L]], digits = 6L), "\n")
        cat("Stdev          ", format(s, digits = 6L), "\n")
        cat("Quantile  0.025", format(q[[1L]], digits = 6L), "\n")
        cat("Quantile  0.25 ", format(q[[2L]], digits = 6L), "\n")
        cat("Quantile  0.5  ", format(q[[3L]], digits = 6L), "\n")
        cat("Quantile  0.75 ", format(q[[4L]], digits = 6L), "\n")
        cat("Quantile  0.975", format(q[[5L]], digits = 6L), "\n")
    }
    return(invisible(list(
        mean = m[[1L]], sd = s, "quant0.025" = q[[1L]],
        "quant0.25" = q[[2L]], "quant0.5" = q[[3L]], "quant0.75" = q[[4L]], "quant0.975" = q[[5L]]
    )))
}

`inla.plot.inla.marginals` <- function(x, ...) {
    stop("NOT IN USE FOR THE MOMENT")
}
