### RCSId = "$Id: marginal.R,v 1.17 2010/02/25 21:48:11 hrue Exp $"

### functions to work with the marginal, either defined as a matrix
### x[2,n], or a list(x=, y=).  NOTE:: there are NO EXTRAPOLATION, so
### the marginal is defined FROM...TO!

`inla.marginal.fix` = function(marginal)
{
    ## just remove points where the density is <= 0
    if (is.matrix(marginal)) {
        i = (marginal[,2] > 0)
        m = list(x=marginal[i,1], y=marginal[i,2])
    } else {
        i = (marginal$y > 0)
        m = list(x= marginal$x[i], y = marginal$y[i])
    }

    return (m)
}
`inla.spline` = function(marginal, log = FALSE, extrapolate = 0.0)
{
    ## for marginal in matrix MARGINAL, which is a marginal density,
    ## return the nice interpolated (x,y) where the interpolation is
    ## done in log(y)

    m = inla.marginal.fix(marginal)
    r = range(m$x)
    r = r[2] - r[1]
    ans = spline(m$x, log(m$y), xmin = min(m$x) - extrapolate * r, xmax = max(m$x) + extrapolate * r)
    if (!log)
        ans$y = exp(ans$y)
    return (ans)
}
`inla.splinefun` = function(marginal)
{
    ## for marginal in matrix MARGINAL, which is a marginal density,
    ## return the spline-function which returns the log-density

    m = inla.marginal.fix(marginal)
    r = range(m$x)
    return (list(range = r, fun = splinefun(m$x, log(m$y))))
}
`inla.expectation` = function(fun, marginal, ...)
{
    ## compute E(FUN(x)), where the marginal of x is given in
    ## `marginal'; see inla.spline().

    xx = inla.spline(marginal)
    n = length(xx$x)
    if (n%%2 == 0)
        n = n -1

    ## use Simpsons integration rule
    i.0 = c(1,n)
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
    f = inla.splinefun(inla.spline(marginal))

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
`inla.pmarginal` = function(x, marginal, normalize = TRUE)
{
    n = length(x)
    p = numeric(n)
    for(i in 1:n) {
        xx.xx.xx = x[i]
        p[i] = inla.expectation(function(x) (1.0*(x < xx.xx.xx)), marginal)
    }
    return (p)
}
`inla.qmarginal` = function(p, marginal, len = 1024)
{
    f = inla.splinefun(inla.spline(marginal))
    xx = seq(f$range[1], f$range[2], length = len)
    d = cumsum(exp(f$fun(xx)))
    d = d/d[length(d)]

    ## just spline-interpolate the inverse mapping
    fq = splinefun(d, xx)

    ## just make sure the p's are in [0,1]
    n = length(p)
    pp = pmin(pmax(p, rep(0,n)), rep(1,n))

    return (fq(pp))
}
`inla.rmarginal` = function(n, marginal)
{
    return (inla.qmarginal(runif(n), marginal))
}

`inla.marginal.transform` = function(fun, marginal, n=256, h.diff = .Machine$double.eps^(1/3), ...)
{
    f = match.fun(fun)
    ff = function(x) f(x, ...)

    m = inla.spline(marginal)
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

    return (list(x = xx, y = log.dens))
}
