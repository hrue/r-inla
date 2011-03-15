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
##! \alias{inla.expectation}
##! \alias{inla.emarginal}
##! \alias{emarginal}
##! \alias{inla.marginal.expectation}
##! \alias{marginal.expectation}
##! \alias{inla.spline}
##! \alias{inla.smarginal}
##! \alias{smarginal}
##! \alias{inla.marginal.transform}
##! \alias{marginal.transform}
##! \alias{inla.tmarginal}
##! 
##! \title{Functions which operates on marginals}
##! 
##! \description{Density, distribution function, quantile function, random
##!      generation, interpolation, expectations and transformations of
##!      marginals obtained by \code{inla} or \code{inla.hyperpar()}.}
##! 
##! \usage{
##! inla.dmarginal = function(x, marginal, log = FALSE)
##! inla.pmarginal = function(q, marginal, normalize = TRUE)
##! inla.qmarginal = function(p, marginal, len = 1024)
##! inla.rmarginal = function(n, marginal)
##! inla.smarginal = function(marginal, log = FALSE, extrapolate = 0.0, keep.type = FALSE)
##! inla.emarginal = function(fun, marginal, ...)
##! inla.tmarginal = function(fun, marginal, n, h.diff, ...)
##!
##! These functions computes the density (inla.dmarginal), 
##! the distribution function (inla.pmarginal), 
##! the quantile function (inla.qmarginal), 
##! random generation (inla.rmarginal), 
##! spline smoothing (inla.smarginal), 
##! computes expected values (inla.emarginal), 
##! and transforms the marginal (inla.tmarginal).
##! }
##! 
##! \arguments{
##! 
##!   \item{marginal}{A marginal object from either \code{inla} or
##!     \code{inla.hyperpar()}, which is either \code{list(x=c(), y=c())}
##!     with density values \code{y} at locations \code{x}, or a
##!     \code{matrix(,n,2)} for which the density values are the second
##!     column and the locations in the first column.}
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
##! }
##! 
##! \value{%%
##!   \code{inla.smarginal} returns \code{list=c(x=c(), y=c())} of
##!   interpolated values do extrapolation using the factor given, whereas
##!   the remaining function returns what they say they should do.  }
##! %%
##! 
##! \author{Havard Rue \email{hrue@math.ntnu.no}}
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
##!            control.data=list(initial = log(1/sd^2),fixed=TRUE))
##! 
##! ## chose a marginal and compare the with the results computed by the
##! ## inla-program
##! r = res$summary.fixed["x",]
##! m = res$marginals.fixed$x
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
##! }
##! 

### functions to work with the marginal, either defined as a matrix
### x[2, n], or a list(x=, y=).  NOTE:: there are NO EXTRAPOLATION, so
### the marginal is defined FROM...TO!

`inla.marginal.fix` = function(marginal)
{
    ## just remove points where the density is <= 0
    if (is.matrix(marginal)) {
        i = (marginal[, 2] > 0)
        m = list(x=marginal[i, 1], y=marginal[i, 2])
    } else {
        i = (marginal$y > 0)
        m = list(x= marginal$x[i], y = marginal$y[i])
    }

    return (m)
}
`inla.spline` = function(marginal, log = FALSE, extrapolate = 0.0) {
    return (inla.smarginal(marginal, log, extrapolate))
}
`inla.smarginal` = function(marginal, log = FALSE, extrapolate = 0.0, keep.type = FALSE)
{
    ## for marginal in matrix MARGINAL, which is a marginal density,
    ## return the nice interpolated (x, y) where the interpolation is
    ## done in log(y)

    is.mat = is.matrix(marginal)
    m = inla.marginal.fix(marginal)
    r = range(m$x)
    r = r[2] - r[1]
    ans = spline(m$x, log(m$y), xmin = min(m$x) - extrapolate * r, xmax = max(m$x) + extrapolate * r)
    if (!log) {
        ans$y = exp(ans$y)
    }

    if (is.mat && keep.type) {
        return (cbind(ans$x, ans$y))
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
`inla.expectation` = function(fun, marginal, ...) {
    return (inla.emarginal(fun, marginal, ...))
}
`inla.emarginal` = function(fun, marginal, ...)
{
    ## compute E(FUN(x)), where the marginal of x is given in
    ## `marginal'; see inla.smarginal().

    xx = inla.smarginal(marginal)
    n = length(xx$x)
    if (n%%2 == 0)
        n = n -1

    ## use Simpsons integration rule
    i.0 = c(1, n)
    i.4 = seq(2, n-1, by=2)
    i.2 = seq(3, n-2, by=2)

    fun = match.fun(fun)
    ff = fun(xx$x[1:n], ...) * xx$y[1:n]
    e = sum(sum(ff[i.0]) + 4*sum(ff[i.4]) + 2*sum(ff[i.2]))

    ## normalise, so that E(1) = 1
    ff = 1 * xx$y[1:n] 
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
`inla.pmarginal` = function(x, marginal, normalize = TRUE, len = 1024)
{
    f = inla.sfmarginal(inla.smarginal(marginal))
    xx = seq(f$range[1], f$range[2], length = len)
    d = cumsum(exp(f$fun(xx)))
    d = d/d[length(d)]

    ## just spline-interpolate the mapping
    fq = splinefun(xx, d)

    ## just make sure the p's are in [0, 1]
    n = length(x)
    xx = pmin(pmax(x, rep(f$range[1], n)), rep(f$range[2], n))

    return (fq(xx))
}
`inla.qmarginal` = function(p, marginal, len = 1024)
{
    f = inla.sfmarginal(inla.smarginal(marginal))
    xx = seq(f$range[1], f$range[2], length = len)
    d = cumsum(exp(f$fun(xx)))
    d = d/d[length(d)]

    ## just spline-interpolate the inverse mapping
    fq = splinefun(d, xx)

    ## just make sure the p's are in [0, 1]
    n = length(p)
    pp = pmin(pmax(p, rep(0, n)), rep(1, n))

    return (fq(pp))
}
`inla.rmarginal` = function(n, marginal)
{
    return (inla.qmarginal(runif(n), marginal))
}

`inla.marginal.transform` = function(fun, marginal, n=256, h.diff = .Machine$double.eps^(1/3), ...) {
    return (inla.tmarginal(fun, marginal, n, h.diff, ...))
}
`inla.tmarginal` = function(fun, marginal, n=256, h.diff = .Machine$double.eps^(1/3), ...)
{
    f = match.fun(fun)
    ff = function(x) f(x, ...)

    is.mat = is.matrix(marginal)
    m = inla.smarginal(marginal)
    r = range(m$x)
    x = seq(r[1], r[2], length = n)
    xx = ff(x)

    ## use the numDeriv library if present
    w = getOption("warn")
    options(warn = -1)
    present = require("numDeriv", quietly = TRUE)
    options(warn = w)

    if (present) {
        ## using numDeriv
        dif = numeric(length(x))
        for(i in 1:length(x))
            dif[i] = numDeriv::grad(ff, x[i])
        log.dens = inla.dmarginal(x, marginal, log=FALSE)/ abs(dif)
    } else {
        ## use a simple algorithm
        log.dens = inla.dmarginal(x, marginal, log=FALSE)/ abs((ff(x + h.diff) - ff(x-h.diff))/(2*h.diff))
    }
    ## if decreasing function, then reverse
    if (xx[1] > xx[n]) {
        xx[1:n] = xx[n:1]
        log.dens[1:n] = log.dens[n:1]
    }

    if (is.mat) {
        return (cbind(x = xx, y = log.dens))
    } else {
        return (list(x = xx, y = log.dens))
    }
}
