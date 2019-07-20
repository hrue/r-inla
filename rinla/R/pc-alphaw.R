## Export: inla.pc.ralphaw inla.pc.dalphaw inla.pc.qalphaw inla.pc.palphaw

##! \name{pc.alphaw}
##! \alias{inla.pc.alphaw}
##! \alias{pc.alphaw}
##! \alias{pc.ralphaw}
##! \alias{inla.pc.ralphaw}
##! \alias{pc.dalphaw}
##! \alias{inla.pc.dalphaw}
##! \alias{pc.palphaw}
##! \alias{inla.pc.palphaw}
##! \alias{pc.qalphaw}
##! \alias{inla.pc.qalphaw}
##! 
##! \title{Utility functions for the PC prior for the \code{alpha} parameter in the Weibull likelihood}
##! 
##! \description{Functions to evaluate, sample, compute quantiles and
##!              percentiles of the PC prior for the \code{alpha} parameter
##!              in the Weibull likelihood}
##! \usage{
##! inla.pc.ralphaw(n, lambda = 1)
##! inla.pc.dalphaw(x, lambda = 1, log = FALSE)
##! inla.pc.qalphaw(p, lambda = 1)
##! inla.pc.palphaw(q, lambda = 1)
##! }
##! \arguments{
##!   \item{n}{Number of observations}
##!   \item{lambda}{The rate parameter (see Details)}
##!   \item{x}{Evaluation points}
##!   \item{log}{Logical. Return the density in natural or log-scale.}
##!   \item{p}{Vector of probabilities}
##!   \item{q}{Vector of quantiles}
##! }
##! \details{
##! This gives the PC prior for the \code{alpha} parameter for the Weibull likelihood,
##! where \code{alpha=1} is the base model. 
##! }
##!\value{%%
##!  \code{inla.pc.dalphaw} gives the density,
##!  \code{inla.pc.palphaw} gives the distribution function,
##!  \code{inla.pc.qalphaw} gives the quantile function, and
##!  \code{inla.pc.ralphaw} generates random deviates.
##! }
##! \seealso{inla.doc("pc.alphaw")}
##! \author{Havard Rue \email{hrue@r-inla.org}}
##! \examples{
##! x = inla.pc.ralphaw(100,  lambda = 1)
##! d = inla.pc.dalphaw(x, lambda = 1)
##! x = inla.pc.qalphaw(0.5, lambda = 1)
##! inla.pc.palphaw(x, lambda = 1)
##! }

inla.pc.alphaw.dist = function(lalpha) 
{
    dist = function(lalpha) {
        alpha = exp(lalpha)
        gam = 0.5772156649
        d = (gamma((1.0 + alpha)/alpha)*alpha + alpha * lalpha - alpha*gam + gam - alpha)/alpha
        return (sqrt(2.0*d))
    }

    tag = "pc.alphaw.dist.cache"
    if (!exists(tag, envir = INLA:::inla.get.inlaEnv())) {
        lalphas = seq(-3, 2, by = 0.01)
        lalphas[which.min(abs(lalphas - 0.0))] = 01.0 ## make it exactly 0
        idx.pos = which(lalphas >= 0)
        idx.neg = which(lalphas <= 0)
        dist.pos = dist(lalphas[idx.pos])
        dist.neg = dist(lalphas[idx.neg])
        assign(tag,
               list(lalpha.pos = lalphas[idx.pos],
                    lalphas.neg = lalphas[idx.neg], 
                    dist.pos = splinefun(lalphas[idx.pos], dist.pos),
                    dist.neg = splinefun(lalphas[idx.neg], dist.neg), 
                    inv.dist.pos = splinefun(dist.pos, lalphas[idx.pos]), 
                    inv.dist.neg = splinefun(dist.neg, lalphas[idx.neg])), 
               envir = INLA:::inla.get.inlaEnv())
    }
    dist.cache = get(tag, envir = INLA:::inla.get.inlaEnv())
    if (missing(lalpha) || length(lalpha) == 0) {
        return(dist.cache)
    } else {
        len = length(lalpha)
        idx.pos = which(lalpha >= 0.0)
        idx.neg = which(lalpha < 0.0)
        d = numeric(len)
        dd = numeric(len)
        d[idx.pos] = dist.cache$dist.pos(lalpha[idx.pos])
        d[idx.neg] = dist.cache$dist.neg(lalpha[idx.neg])
        dd[idx.pos] = dist.cache$dist.pos(lalpha[idx.pos], deriv = 1) 
        dd[idx.neg] = dist.cache$dist.neg(lalpha[idx.neg], deriv = 1) 
    }
    return (list(dist = d, deriv = dd))
}

inla.pc.alphaw.intern = function(lambda = 1)
{
    d = inla.pc.alphaw.dist()
    idx.pos = which(d$alpha >= 1.0)
    idx.neg = which(d$alpha <= 1.0)
    alpha.pos = d$alpha[idx.pos]
    alpha.neg = d$alpha[idx.neg]
    marg = list(marg.neg = list(
                    x = alpha.neg, 
                    y = 0.5 * lambda * exp(-lambda * d$dist.neg(alpha.neg)) *
                        abs(d$dist.neg(alpha.neg, deriv=1))),
                marg.pos = list(
                    x = alpha.pos, 
                    y = 0.5 * lambda * exp(-lambda * d$dist.pos(alpha.pos)) *
                        abs(d$dist.pos(alpha.pos, deriv=1))))
    return (marg)
}

inla.pc.ralphaw = function(n, lambda = 1)
{
    marg = inla.pc.alphaw.intern()
    ind = sample(c(0, 1), n, replace=TRUE)
    idx.pos = which(ind == 1)
    idx.neg = which(ind == 0)
    x = numeric(n)

    m = length(idx.pos)
    if (m > 0) {
        x[idx.pos] = inla.rmarginal(runif(m), marg$marg.pos)
    }
    m = length(idx.neg)
    if (m > 0) {
        x[idx.neg] = inla.rmarginal(runif(m), marg$marg.neg)
    }
    return (x)
}

inla.pc.dgamma = function(x, lambda = 1, log = FALSE)
{
    inv.x = 1/x
    d = sqrt(2*(log(inv.x) - psigamma(inv.x)))
    ld = -lambda * d - log(d) + log(lambda) - 2*log(x) +
        log(psigamma(inv.x, deriv=1) - x)
    return (if (log) ld else exp(ld))
}

inla.pc.qgamma = function(p, lambda = 1)
{
    log.x = inla.pmarginal(p, inla.pc.gamma.intern(lambda = lambda))
    return (exp(log.x))
}

inla.pc.pgamma = function(q, lambda = 1)
{
    p = inla.qmarginal(log(q), inla.pc.gamma.intern(lambda = lambda))
    return (p)
}
