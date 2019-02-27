## Export: inla.dmarginal inla.pmarginal inla.qmarginal inla.rmarginal inla.hpdmarginal
## Export: inla.smarginal inla.emarginal inla.tmarginal inla.mmarginal inla.zmarginal
## NOT-YET-IN-USE-Export: summary!inla.marginal plot!inla.marginal

##! \name{marginal}
##! \alias{inla.marginal}
##! \alias{marginal}
##! \alias{pmarginal}
##! \alias{inla.pmarginal}
##! \alias{qmarginal}
##! \alias{inla.qmarginal}
##! \alias{dmarginal}
##! \alias{inla.dmarginal}
##! \alias{rmarginal}
##! \alias{inla.rmarginal}
##! \alias{inla.hpdmarginal}
##! \alias{hpdmarginal}
##! \alias{inla.emarginal}
##! \alias{emarginal}
##! \alias{inla.marginal.expectation}
##! \alias{marginal.expectation}
##! \alias{inla.smarginal}
##! \alias{smarginal}
##! \alias{inla.marginal.transform}
##! \alias{marginal.transform}
##! \alias{inla.tmarginal}
##! \alias{inla.mmarginal}
##! \alias{mmarginal}
##! \alias{inla.zmarginal}
##! \alias{zmarginal}
##! 
##! \title{Functions which operates on marginals}
##! 
##! \description{Density, distribution function, quantile function, random
##!      generation, hpd-interval, interpolation, expectations, mode and transformations of
##!      marginals obtained by \code{inla} or \code{inla.hyperpar()}.
##! These functions computes the density (inla.dmarginal), 
##! the distribution function (inla.pmarginal), 
##! the quantile function (inla.qmarginal), 
##! random generation (inla.rmarginal), 
##! spline smoothing (inla.smarginal), 
##! computes expected values (inla.emarginal), 
##! computes the mode (inla.mmarginal), 
##! transforms the marginal (inla.tmarginal), and provide summary statistics (inla.zmarginal).
##! }
##! 
##! \usage{
##! inla.dmarginal(x, marginal, log = FALSE)
##! inla.pmarginal(q, marginal, normalize = TRUE, len = 2048L)
##! inla.qmarginal(p, marginal, len = 2048L)
##! inla.rmarginal(n, marginal)
##! inla.hpdmarginal(p, marginal, len = 2048L)
##! inla.smarginal(marginal, log = FALSE, extrapolate = 0.0, keep.type = FALSE, factor=15L)
##! inla.emarginal(fun, marginal, ...)
##! inla.mmarginal(marginal)
##! inla.tmarginal(fun, marginal, n=2048L, h.diff = .Machine$double.eps^(1/3),
##!                method = c("quantile", "linear")) 
##! inla.zmarginal(marginal, silent = FALSE)
##! }
##! \arguments{
##! 
##!   \item{marginal}{A marginal object from either \code{inla} or
##!     \code{inla.hyperpar()}, which is either \code{list(x=c(), y=c())}
##!     with density values \code{y} at locations \code{x}, or a
##!     \code{matrix(,n,2)} for which the density values are the second
##!     column and the locations in the first column.
##!     The\code{inla.hpdmarginal()}-function
##!     assumes a unimodal density.}
##! 
##!   \item{fun}{A (vectorised) function like \code{function(x) exp(x)} to
##!     compute the expectation against, or which define the transformation
##!     new = fun(old)}
##!   
##!   \item{x}{Evaluation points}
##! 
##!   \item{q}{Quantiles}
##! 
##!   \item{p}{Probabilities}
##! 
##!   \item{n}{The number of observations. If \code{length(n) > 1}, the
##!     length is taken to be the number required. For
##!     \code{inla.marginal.transform}, its the number of points to use
##!     in the new density.}
##! 
##!   \item{h.diff}{The step-length for the numerical differeniation inside \code{inla.marginal.transform}}
##!     
##!   \item{...}{Further arguments to be passed to function which
##!     expectation is to be computed.}
##! 
##!   \item{log}{Return density or interpolated density in log-scale?}
##! 
##!   \item{normalize}{Renormalise the density after interpolation?}
##!   \item{len}{Number of locations used to interpolate the distribution
##!   function.}
##!
##!   \item{keep.type}{If \code{FALSE} then return a \code{list(x=, y=)},  otherwise if \code{TRUE},
##!                   then return a matrix if the input is a matrix}
##! 
##!   \item{extrapolate}{How much to extrapolate on each side when computing the
##!   interpolation. In fraction of the range.}
##!
##!   \item{factor}{The number of points after interpolation is \code{factor} times the original number of points;
##!                 which is argument \code{n} in \code{spline}}
##!   
##!    \item{method}{Which method should be used to layout points for where the transformation is computed.}
##!
##!    \item{silent}{Output the result visually (TRUE) or just through the call.}
##! }
##! 
##! \value{%%
##!   \code{inla.smarginal} returns \code{list=c(x=c(), y=c())} of
##!   interpolated values do extrapolation using the factor given, 
##!   and the remaining function returns what they say they should do.  }
##! %%
##! 
##! \author{Havard Rue \email{hrue@r-inla.org}}
##! 
##! \seealso{\code{\link{inla}}, \code{\link{inla.hyperpar}}}
##! 
##! \examples{
##! ## a simple linear regression example
##! n = 10
##! x = rnorm(n)
##! sd = 0.1
##! y = 1+x + rnorm(n,sd=sd)
##! res = inla(y ~ 1 + x, data = data.frame(x,y),
##!            control.family=list(initial = log(1/sd^2),fixed=TRUE))
##! 
##! ## chose a marginal and compare the with the results computed by the
##! ## inla-program
##! r = res$summary.fixed["x",]
##! m = res$marginals.fixed$x
##!
##! ## compute the 95% HPD interval
##! inla.hpdmarginal(0.95, m)
##!
##! x = seq(-6, 6, len = 1000)
##! y = dnorm(x)
##! inla.hpdmarginal(0.95, list(x=x, y=y))
##!
##! ## compute the the density for exp(r), version 1
##! r.exp = inla.tmarginal(exp, m)
##! ## or version 2
##! r.exp = inla.tmarginal(function(x) exp(x), m)
##! 
##! ## to plot the marginal, we use the inla.smarginal, which interpolates (in
##! ## log-scale). Compare with some samples.
##! plot(inla.smarginal(m), type="l")
##! s = inla.rmarginal(1000, m)
##! hist(inla.rmarginal(1000, m), add=TRUE, prob=TRUE)
##! lines(density(s), lty=2)
##! 
##! m1 = inla.emarginal(function(x) x^1, m)
##! m2 = inla.emarginal(function(x) x^2, m)
##! stdev = sqrt(m2 - m1^2)
##! q = inla.qmarginal(c(0.025,0.975), m)
##! 
##! ## inla-program results
##! print(r)
##! 
##! ## inla.marginal-results (they shouldn't be perfect!)
##! print(c(mean=m1, sd=stdev, "0.025quant" = q[1], "0.975quant" = q[2]))
##! ## using the buildt-in function
##! inla.zmarginal(m)
##! }
##! 

### functions to work with the marginal, either defined as a matrix
### x[2, n], or a list(x=, y=).  NOTE:: there are NO EXTRAPOLATION, so
### the marginal is defined FROM...TO!

`inla.marginal.fix` = function(marginal)
{
    ## just remove points where the density is <= 0 or where the
    ## density is to small compared to the maximum density. (othewise
    ## we can get trouble with the spline interpolation). same with
    ## 'x'? No...
    eps = .Machine$double.eps * 1000
    ##marginal = spline(marginal)
    if (is.matrix(marginal)) {
        if (any(is.na(marginal[, 2]))) {
            idx = which(is.na(marginal[, 2]))
            marginal = marginal[-idx, ]
        }
        i = (marginal[, 2] > 0) & (abs(marginal[, 2]/max(marginal[, 2])) > eps)
        m = list(x=marginal[i, 1], y=marginal[i, 2])
        ## i = c(diff(marginal[, 1]) > eps, TRUE)
        ## m = list(x=marginal[i, 1], y=marginal[i, 2])
    } else {
        if (any(is.na(marginal$y))) {
            idx = which(is.na(marginal$y))
            marginal$x = marginal$x[-idx]
            marginal$y = marginal$y[-idx]
        }
        i = (marginal$y > 0) & (abs(marginal$y/max(marginal$y)) > eps)
        m = list(x = marginal$x[i], y = marginal$y[i])
        ## i = c(diff(marginal$x) > eps, TRUE)
        ## m = list(x = marginal$x[i], y = marginal$y[i])
    }

    return (m)
}

`inla.spline` = function(x, ...)
{
    s = spline(x, ...)
    if (is.matrix(x)) {
        m = cbind(x = s$x, y = s$y)
    } else if (is.list(x)) {
        m = list(x = s$x, y = s$y)
    } else {
        m = s
    }
    return (m)
}

`inla.smarginal` = function(marginal, log = FALSE, extrapolate = 0.0, keep.type = FALSE, factor=15L)
{
    ## for marginal in matrix MARGINAL, which is a marginal density,
    ## return the nice interpolated (x, y) where the interpolation is
    ## done in log(y)

    is.mat = is.matrix(marginal)
    m = inla.marginal.fix(marginal)
    r = diff(range(m$x))
    xmin = min(m$x) - extrapolate * r
    xmax = max(m$x) + extrapolate * r
    n = factor * length(m$x)
    xx = seq(xmin, xmax, len = n)
    if (extrapolate) {
        xx = c(xmin, m$x, xmax)
    } else {
        xx = m$x
    }
    nx = length(xx)
    dx = nx * diff(xx) / mean(diff(xx))
    xnew = c(0, cumsum(sqrt(dx)))
    xnew = xmin + (xmax - xmin) * ((xnew - min(xnew))/(max(xnew) - min(xnew)))
    fun = splinefun(xnew, xx, method = "hyman")
    fun.inv = splinefun(xx, xnew, method = "hyman")
    ans = spline(fun.inv(m$x), log(m$y), xmin = fun.inv(xmin), xmax = fun.inv(xmax), n = n, method = "fmm")
    ans$x = fun(ans$x)
    if (!log) {
        ans$y = exp(ans$y)
        ans = inla.marginal.fix(ans)
    }
    if (is.mat && keep.type) {
        return (cbind(x=ans$x, y=ans$y))
    } else {
        return (ans)
    }
}

`inla.sfmarginal` = function(marginal)
{
    ## for marginal in matrix MARGINAL, which is a marginal density,
    ## return the spline-function which returns the log-density

    m = inla.marginal.fix(marginal)
    r = range(m$x)
    return (list(range = r, fun = splinefun(m$x, log(m$y))))
}

`inla.emarginal` = function(fun, marginal, ...)
{
    ## compute E(FUN(x)), where the marginal of x is given in
    ## `marginal'; see inla.smarginal(). Also work for FUN(x)
    ## returning a vector.

    xx = inla.smarginal(marginal)
    n = length(xx$x)
    if (n%%2 == 0) {
        n = n -1
        xx$x = xx$x[1:n]
        xx$y = xx$y[1:n]
    }
    
    ## use Simpsons integration rule
    i.0 = c(1, n)
    i.4 = seq(2, n-1, by=2)
    i.2 = seq(3, n-2, by=2)

    dx = diff(xx$x)
    dx = 0.5 * (c(dx, 0) + c(0, dx))
    fun = match.fun(fun)
    ff = fun(xx$x[1:n], ...) * xx$y[1:n] * dx
    nf = length(ff) %/% n
    e = numeric(nf)
    off = 0L
    for(i in 1:nf) {
        e[i] = sum(sum(ff[i.0 + off]) + 4*sum(ff[i.4 + off]) + 2*sum(ff[i.2 + off]))
        off = off + n
    }

    ## normalise, so that E(1) = 1
    ff = dx * xx$y[1:n] 
    e.1 = sum(sum(ff[i.0]) + 4*sum(ff[i.4]) + 2*sum(ff[i.2]))

    return (e/e.1)
}

`inla.dmarginal` = function(x, marginal, log = FALSE)
{
    ## return the density of the marginal. if LOG, in log-scale.  this
    ## density is not renormalized (but should be already normalized
    ## in some scale).

    ## return list(range=, fun=..)
    f = inla.sfmarginal(inla.smarginal(marginal))

    n = length(x)
    d = numeric(n)

    for(i in 1:n) {
        if (x[i] >= f$range[1] && x[i] <= f$range[2]) {
            if (log)
                d[i] = f$fun(x[i])
            else
                d[i] = exp(f$fun(x[i]))
        } else {
            if (log)
                d[i] = -Inf
            else
                d[i] = 0.0
        }
    }
    return (d)
}

`inla.pmarginal` = function(q, marginal, normalize = TRUE, len = 2048L)
{
    f = inla.sfmarginal(inla.smarginal(marginal))
    xx = seq(f$range[1], f$range[2], length = len)
    d = cumsum(exp(f$fun(xx)))
    d = d/d[length(d)]

    ## just spline-interpolate the mapping
    fq = splinefun(xx, d, method = "hyman")

    ## just make sure the p's are in [0, 1]
    n = length(q)
    xx = pmin(pmax(q, rep(f$range[1], n)), rep(f$range[2], n))

    return (fq(xx))
}

`inla.qmarginal` = function(p, marginal, len = 2048L)
{
    f = inla.sfmarginal(inla.smarginal(marginal))
    xx = seq(f$range[1], f$range[2], length = len)
    d = cumsum(exp(f$fun(xx)))
    d = d/d[length(d)]

    ## for the moment, we remove only duplicated zero's and one's.
    eps = .Machine$double.eps * 1000.0
    for(val in c(0.0, 1.0)) {
        is.val = which(abs(d - val) <= eps)
        if (length(is.val) > 1) {
            ## keep the first, ie remove it from the list to the removed
            is.val = is.val[-1]
            ## and remove the rest
            d = d[-is.val]
            xx = xx[-is.val]
        }
    }

    ## just spline-interpolate the inverse mapping
    fq = splinefun(d, xx, method = "hyman")

    ## just make sure the p's are in [0, 1]
    n = length(p)
    pp = pmin(pmax(p, rep(0, n)), rep(1, n))

    return (fq(pp))
}
`inla.hpdmarginal` = function(p, marginal, len = 2048L)
{
    sm = inla.smarginal(marginal, keep.type = FALSE)
    f = inla.sfmarginal(sm)
    x.range = f$range
    xx = seq(x.range[1], x.range[2], length = len)
    ## simply add 0 density outside the range.
    dx = diff(xx)[1]
    ef = exp(f$fun(xx))
    xx = c(min(xx) - dx,  xx,  max(xx) + dx)
    d = c(0, ef, 0)
    d = cumsum(d)
    d = d/d[length(d)]

    ## for the moment, we remove only duplicated zero's and one's.
    eps = .Machine$double.eps * 1000.0
    for(val in c(0.0, 1.0)) {
        is.val = which(abs(d - val) <= eps)
        if (length(is.val) > 1) {
            ## keep the first, ie remove it from the list to the removed
            is.val = is.val[-1]
            ## and remove the rest
            d = d[-is.val]
            xx = xx[-is.val]
        }
    }
    ## just spline-interpolate the inverse mapping
    fq = splinefun(d, xx, method = "hyman")

    ## just make sure the p's are in [0, 1]
    np = length(p)
    pp = 1-pmin(pmax(p, rep(0, np)), rep(1, np))
    
    ## parts of this code below is taken from TeachingDemo::hpd
    f = function(x, posterior.icdf, conf) {
        return (posterior.icdf(1 - conf + x) - posterior.icdf(x))
    }

    tol= sqrt(.Machine$double.eps)
    tol = 1E-6
    result = matrix(NA, np, 2)
    for(i in 1:np) {
        out = optimize(f, c(0, pp[i]), posterior.icdf = fq, conf = pp[i], tol = tol)
        result[i, ] = c(fq(out$minimum), fq(1 - pp[i] + out$minimum))
        ## make sure its in the range of the marginal itself as we expanded the domain
        result[i, ] = c(min(x.range[2], max(result[i, 1], x.range[1])),
                        min(x.range[2], max(result[i, 2], x.range[1])))
    }
    colnames(result) = c("low", "high")
    rownames(result) = paste("level:", format(1-pp, digits=6, justify="left", trim=TRUE), sep="")

    return (result)
}
`inla.rmarginal` = function(n, marginal)
{
    return (inla.qmarginal(runif(n), marginal))
}

`inla.marginal.transform` = function(fun, marginal, n=2048L, h.diff = .Machine$double.eps^(1/3),
        method = c("quantile", "linear"))
{
    return (inla.tmarginal(fun, marginal, n, h.diff, method = method))
}

`inla.deriv.func` = function(fun, step.size = .Machine$double.eps^(1/4))
{
    ## return a function computing the derivatives to 'fun = function(x)' which must
    ## vectorize
    func = match.fun(fun)
    if (inla.require("Deriv")) {
        fd = try(Deriv::Deriv(func, names(formals(func))[1]), silent=TRUE)
        if (!inherits(fd, "try-error"))
            return (fd)
    }
    if (inla.require("numDeriv")) {
        fd = function(x) unlist(lapply(x, function(xi) numDeriv::grad(func, xi)))
    } else {
        fd = function(x) ((-func(x + 2*step.size) + 8*func(x + step.size) -
                           8*func(x - step.size) + func(x - 2*step.size))/(12*step.size))
    }
    
    return(fd)
}

`inla.tmarginal` = function(fun, marginal, n=2048L, h.diff = .Machine$double.eps^(1/3),
        method = c("quantile", "linear")) 
{
    ff = match.fun(fun)

    is.mat = is.matrix(marginal)
    m = inla.smarginal(marginal)
    r = range(m$x)

    method = match.arg(method)
    if (inla.strcasecmp(method, "quantile")) {
        x = inla.qmarginal((1:n)/(n+1), marginal)
    } else if (inla.strcasecmp(method, "linear")) {
        x = seq(r[1], r[2], length = n)
    } else {
        stop("unknown method")
    }
    xx = ff(x)

    fd = inla.deriv.func(ff)
    log.dens = inla.dmarginal(x, marginal, log=FALSE)/abs(fd(x))

    ## if decreasing function, then reverse
    if (xx[1] > xx[n]) {
        xx[1:n] = xx[n:1]
        log.dens[1:n] = log.dens[n:1]
    }

    if (is.mat) {
        ret = cbind(x = xx, y = log.dens)
    } else {
        ret = list(x = xx, y = log.dens)
    }

    if (FALSE) {
        class(ret) = "inla.marginal"
        attr(ret, "inla.tag") = paste(attr(marginal, "inla.tag"), "transformed")
    }

    return (ret)
}

`inla.mmarginal` = function(marginal)
{
    p = seq(0.01, 0.99, by=0.01)
    n = length(p)
    x = inla.qmarginal(p, marginal)
    d = inla.dmarginal(x, marginal, log=TRUE)
    idx = which.max(d)
    res = optimize(
            f = inla.dmarginal,
            interval = c(x[max(1L, idx-1L)], x[min(n, idx+1L)]),
            maximum = TRUE,
            ## arguments to inla.dmarginal
            marginal = marginal, log=TRUE)
    return(res$maximum)
}

`inla.zmarginal` = function(marginal, silent = FALSE)
{
    stopifnot(inla.is.marginal(marginal))

    m = inla.emarginal(function(xx) c(xx, xx^2), marginal)
    q = inla.qmarginal(c(0.025, 0.25, 0.5, 0.75, 0.975), marginal)
    s = sqrt(max(0, m[2]-m[1]^2))

    if (!silent) {
        ##cat(as.character(match.call()[-1]), "\n")
        ##cat("+--------+-----+------------\n")
        cat("Mean           ", format(m[1], digits=6), "\n")
        cat("Stdev          ", format(s, digits=6), "\n")
        cat("Quantile  0.025", format(q[1], digits=6), "\n")
        cat("Quantile  0.25 ", format(q[2], digits=6), "\n")
        cat("Quantile  0.5  ", format(q[3], digits=6), "\n")
        cat("Quantile  0.75 ", format(q[4], digits=6), "\n")
        cat("Quantile  0.975", format(q[5], digits=6), "\n")
    }
    return (invisible(list(mean = m[1],  sd = s,  "quant0.025" = q[1],
                           "quant0.25" = q[2], "quant0.5" = q[3], "quant0.75"=q[4], "quant0.975" = q[5])))
}

`inla.is.marginal` = function(m)
{
    return ((is.matrix(m) && ncol(m) == 2L && nrow(m) > 2L) ||
            (is.list(m) && all(names(m) == c("x", "y"))))
}
    

## 'plot' and 'summary'-methods for marginals
plot.inla.marginal = function(x, y, ...)
{
    m = inla.emarginal(function(xx) c(xx, xx^2), x)
    xlab = paste("x (mean", format(m[1], digits=4),
            ", stdev", format(sqrt(max(0, m[2]-m[1]^2)), digits=4), ")")
    plot(inla.smarginal(x, extrapolate=0), type = "l", xlab = xlab, ylab = "marginal density", ...)

    return (invisible())
}
summary.inla.marginal = function(object, ...)
{
    ## option
    silent=FALSE
    inla.eval.dots(...)

    m = inla.emarginal(function(xx) c(xx, xx^2), object)
    q = inla.qmarginal(c(0.025, 0.25, 0.5, 0.75, 0.975), object)
    s = sqrt(max(0, m[2]-m[1]^2))

    if (!silent) {
        cat("Properties of: ", "FIXME", "\n")
        cat("+--------------+------------------\n")
        cat("Mean           ", format(m[1], digits=6), "\n")
        cat("Stdev          ", format(s, digits=6), "\n")
        cat("Quantile  0.025", format(q[1], digits=6), "\n")
        cat("Quantile  0.25 ", format(q[2], digits=6), "\n")
        cat("Quantile  0.5  ", format(q[3], digits=6), "\n")
        cat("Quantile  0.75 ", format(q[4], digits=6), "\n")
        cat("Quantile  0.975", format(q[5], digits=6), "\n")
    }
    return (invisible(list(mean = m[1],  sd = s,  "quant0.025" = q[1],
                           "quant0.25" = q[2], "quant0.5" = q[3], "quant0.75"=q[4], "quant0.975" = q[5])))
}

`inla.plot.inla.marginals` = function(x, ...)
{
    stop("NOT IN USE FOR THE MOMENT")

    ## input here is a list of marginals

    n = length(x)
    if (n == 1) {
        ## length is 1,  so here we plot the marginal itself
        plot(x[[1]], ...)
    } else {
        xx = sapply(x, summary, silent=TRUE)

        stopifnot(dim(xx)[1] == 7)
        stopifnot(dim(xx)[2] == n)

        xx.min = min(as.numeric(xx))
        xx.max = max(as.numeric(xx))

        ylim = range(pretty(c(xx.min, xx.max)))
        lab = c("mean",  "quant0.5", "quant0.025", "quant0.975", "quant0.25", "quant0.75")
        xval = 1:n

        plot(xval, rep(ylim[2]*1000, n),
             ylim = ylim, main = "FIXME", 
             xlab = "x", ylab = inla.paste(lab, sep=", "), ...)
    
        polygon(c(xval, rev(xval)), c(xx[lab[3], ], rev(xx[lab[4], ])), col="lightgray", border=NA, ...)
        polygon(c(xval, rev(xval)), c(xx[lab[5], ], rev(xx[lab[6], ])), col="gray", border=NA, ...)

        ## mean
        j=1; lines(xval, xx[lab[j], ], lwd=1, lty=j, ...)
        j=1; points(xval, xx[lab[j], ], pch=20, cex=0.4)
        ## median
        j=2; lines(xval, xx[lab[j], ], lwd=1, lty=j, ...)
    
        return (invisible())
    }
}
