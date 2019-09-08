## export: inla.pc.rsn inla.pc.dsn inla.pc.qsn inla.pc.psn

##! \name{pc.sn}
##! \alias{inla.pc.sn}
##! \alias{pc.sn}
##! \alias{pc.rsn}
##! \alias{inla.pc.rsn}
##! \alias{pc.dsn}
##! \alias{inla.pc.dsn}
##! \alias{pc.psn}
##! \alias{inla.pc.psn}
##! \alias{pc.qsn}
##! \alias{inla.pc.qsn}
##! 
##! \title{Utility functions for the PC prior for the \code{alpha} parameter in the skew-normal
##! linkfunction}
##! 
##! \description{Functions to evaluate, sample, compute quantiles and
##!              percentiles of the PC prior for the \code{alpha} parameter
##!              in the skew-normal link-function}
##! \usage{
##! inla.pc.rsn(n, lambda = 30)
##! inla.pc.dsn(alpha, lambda = 30, log = FALSE)
##! inla.pc.qsn(p, lambda = 30)
##! inla.pc.psn(q, lambda = 30)
##! }
##! \arguments{
##!   \item{n}{number of observations}
##!   \item{lambda}{the rate parameter in the PC prior}
##!   \item{alpha}{vector of evaluation points}
##!   \item{log}{logical. return the density in natural or log-scale.}
##!   \item{p}{vector of probabilities}
##!   \item{q}{vector of quantiles}
##! }
##! \details{
##! Defines the PC prior for the \code{alpha} parameter for the skew-normal linkfunction
##! where \code{alpha=0} is the base model. 
##! }
##!\value{%%
##!  \code{inla.pc.dsn} gives the density,
##!  \code{inla.pc.psn} gives the distribution function,
##!  \code{inla.pc.qsn} gives the quantile function, and
##!  \code{inla.pc.rsn} generates random deviates.
##! }
##! \seealso{inla.doc("pc.sn")}
##! \author{havard rue \email{hrue@r-inla.org}}
##! \examples{
##! x = inla.pc.rsn(100,  lambda = 30)
##! d = inla.pc.dsn(x, lambda = 30)
##! x = inla.pc.qsn(0.5, lambda = 30)
##! inla.pc.psn(x, lambda = 30)
##! }

map = function(alpha) sign(alpha) * abs(alpha)^(1/3)

inla.pc.sn.cache = function() 
{
    ## return the cache for these functions
    sn = function(x, alpha) {
        a = map(alpha)
        delta = a/sqrt(1 + a^2)
        omega = 1/sqrt(1-2*delta^2/pi)
        xi = -omega * delta * sqrt(2/pi)
        xx = (x - xi)/omega
        return (2/omega * dnorm(xx) * pnorm(a*xx))
    }

    kld.arg = function(x, alpha) {
        dsn = sn(x, alpha)
        idx = which(dsn > 0.0)
        dsn = dsn[idx]
        x = x[idx]

        if (length(x) %%2==0) {
            x = x[-length(x)]
            dsn = dsn[-length(dsn)]
        }
        w = c(1, rep(c(4, 2), (length(x)-3) %/% 2), 4, 1)/3

        dn = dnorm(x)
        dsn = dsn / sum(w*dsn)
        dn = dn / sum(w*dn)
        return (sum(w*dsn*log(dsn/dn)))
    }

    kld = function(alpha) {
        dx = 0.001
        ran = 10
        x = seq(-ran, ran, by = dx)
        return (kld.arg(x, alpha))
    }

    dist = function(alpha) {
        d = numeric(length(alpha))
        for(i in seq_along(alpha)) {
            d[i] = sqrt(2.0 * max(0, kld(alpha[i])))
        }
        idx = which(!(is.na(d) | is.infinite(d) | d < 0.0))
        return (data.frame(alpha[idx], dist=d[idx]))
    }

    tag = "cache.pc.sn"
    if (!exists(tag, envir = inla.get.inlaEnv())) {
        ran = 100
        alphas = exp(seq(log(0.001), log(ran), len=99))
        alphas = c(-rev(alphas), 0, alphas)
        alphas.pos = alphas[which(alphas >= 0)]
        alphas.neg = alphas[which(alphas <= 0)]
        dist.pos = dist(alphas.pos)$dist
        dist.neg = rev(dist.pos)
        ## write(alphas.pos, ncol = 1, file="aa")
        ## write(dist.pos, ncol = 1, file="dd")
        assign(tag,
               list(alpha.range = range(alphas),
                    pos = list(
                        x = alphas.pos, 
                        dist = splinefun(alphas.pos, dist.pos),
                        idist = splinefun(dist.pos, alphas.pos)), 
                    neg = list(
                        x = alphas.neg, 
                        dist = splinefun(alphas.neg, dist.neg), 
                        idist = splinefun(dist.neg, alphas.neg))), 
               envir = inla.get.inlaEnv())
    }
    return (get(tag, envir = inla.get.inlaEnv()))
}

inla.pc.rsn = function(n, lambda = 30)
{
    cache = inla.pc.sn.cache()
    x = numeric(n)
    ind = sample(c(-1, 1), n, replace=TRUE)
    idx.pos = which(ind > 0)
    idx.neg = which(ind < 0)
    z = rexp(n, rate = lambda)
    if (length(idx.pos) > 0) {
        x[idx.pos] = cache$pos$idist(z[idx.pos])
    }
    if (length(idx.neg) > 0) {
        x[idx.neg] = cache$neg$idist(z[idx.neg])
    }
    return (x)
}

inla.pc.dsn = function(alpha, lambda = 30, log = FALSE)
{
    cache = inla.pc.sn.cache()
    d = numeric(length(alpha))
    idx.pos = which(alpha >= 0)
    idx.neg = which(alpha <  0)
    if (length(idx.pos) > 0) {
        la = alpha[idx.pos]
        dist = cache$pos$dist(la)
        deriv = cache$pos$dist(la, deriv = 1)
        d[idx.pos] = log(0.5) + log(lambda) - lambda * dist + log(abs(deriv)) 
    }        
    if (length(idx.neg) > 0) {
        la = alpha[idx.neg]
        dist = cache$neg$dist(la)
        deriv = cache$neg$dist(la, deriv = 1)
        d[idx.neg] = log(0.5) + log(lambda) - lambda * dist + log(abs(deriv))
    }        
    return (if (log) d else exp(d))
}

inla.pc.qsn = function(p, lambda = 30)
{
    cache = inla.pc.sn.cache()
    n = length(p)
    q = numeric(n)
    idx.pos = which(p >= 0.5)
    idx.neg = which(p <  0.5)
    if (length(idx.pos) > 0) {
        pp = 2.0*(p[idx.pos] - 0.5)
        qe = qexp(pp, rate = lambda)
        q[idx.pos] = cache$pos$idist(qe)
    }
    if (length(idx.neg) > 0) {
        pp = 2.0*(0.5 - p[idx.neg])
        qe = qexp(pp, rate = lambda)
        q[idx.neg] = cache$neg$idist(qe)
    }
    return(q)
}

inla.pc.psn = function(q, lambda = 30)
{
    cache = inla.pc.sn.cache()
    n = length(q)
    p = numeric(n)
    idx.pos = which(q >= 0.0)
    idx.neg = which(q <  0.0)
    if (length(idx.pos) > 0) {
        qq = cache$pos$dist(q[idx.pos])
        pp = pexp(qq, rate = lambda)
        p[idx.pos] = (1.0 + pp)/2.0
    }
    if (length(idx.neg) > 0) {
        qq = cache$neg$dist(q[idx.neg])
        pp = pexp(qq, rate = lambda)
        p[idx.neg] = (1.0 - pp)/2.0
    }
    return(p)
}

inla.pc.sn.test = function(lambda = 30) 
{
    ## this is just an internal test, and is not exported
    f = function(...) format(...,  digits = 3)
    n = 10^6
    x = inla.pc.rsn(n, lambda = lambda)
    x = x[x < quantile(x, prob = 0.999)]
    x = x[x > quantile(x, prob = 1-0.999)]
    alpha = seq(min(x), max(x), by = 0.001)
    hist(x, n = 300,  prob = TRUE)
    lines(alpha, inla.pc.dsn(alpha, lambda), lwd=3, col="blue")
    p = seq(0.1, 0.9, by = 0.05)
    print(cbind(p=f(p), q.emp = f(as.numeric(quantile(x,  prob = p))),
                q.true = f(inla.pc.qsn(p, lambda))))
    q = as.numeric(quantile(x, prob = p))
    cdf = ecdf(x)
    print(cbind(q=f(q), p.emp = f(cdf(q)), p.true = f(inla.pc.psn(q, lambda))))
}

## to run the test with a sourced file,  do
##
##     inla.get.inlaEnv = function(...) INLA:::inla.get.inlaEnv(...)
##     inla.pc.sn.test()
##
