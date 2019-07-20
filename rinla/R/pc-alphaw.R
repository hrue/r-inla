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
##! inla.pc.ralphaw(n, lambda = 7)
##! inla.pc.dalphaw(alpha, lambda = 7, log = FALSE)
##! inla.pc.qalphaw(p, lambda = 7)
##! inla.pc.palphaw(q, lambda = 7)
##! }
##! \arguments{
##!   \item{n}{Number of observations}
##!   \item{lambda}{The rate parameter in the PC-prior}
##!   \item{alpha}{Vector of evaluation points, where \code{alpha>0}.}
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
##! x = inla.pc.ralphaw(100,  lambda = 7)
##! d = inla.pc.dalphaw(x, lambda = 7)
##! x = inla.pc.qalphaw(0.5, lambda = 7)
##! inla.pc.palphaw(x, lambda = 7)
##! }

inla.pc.alphaw.cache = function() 
{
    dist = function(lalpha) {
        alpha = exp(lalpha)
        gam = 0.5772156649
        kld = (gamma((1.0 + alpha)/alpha)*alpha + alpha * lalpha - alpha*gam + gam - alpha)/alpha
        return (sqrt(2.0*kld))
    }

    tag = "pc.alphaw.dist.cache"
    if (!exists(tag, envir = INLA:::inla.get.inlaEnv())) {
        lalphas = seq(-5, 3, by = 0.01)
        lalphas[which.min(abs(lalphas))] = 0 ## yes, make it exactly 0
        lalphas.pos = lalphas[which(lalphas >= 0)]
        lalphas.neg = lalphas[which(lalphas <= 0)]
        dist.pos = dist(lalphas.pos)
        dist.neg = dist(lalphas.neg)
        assign(tag,
               list(pos = list(
                        x = lalphas.pos, 
                        dist = splinefun(lalphas.pos, dist.pos),
                        idist = splinefun(dist.pos, lalphas.pos)), 
                    neg = list(
                        x = lalphas.neg, 
                        dist = splinefun(lalphas.neg, dist.neg), 
                        idist = splinefun(dist.neg, lalphas.neg))), 
               envir = INLA:::inla.get.inlaEnv())
    }
    return (get(tag, envir = INLA:::inla.get.inlaEnv()))
}

inla.pc.ralphaw = function(n, lambda = 7)
{
    cac = inla.pc.alphaw.cache()
    x = numeric(n)
    ind = sample(c(-1, 1), n, replace=TRUE)
    idx.pos = which(ind > 0)
    idx.neg = which(ind < 0)
    z = rexp(n, rate = lambda)
    if (length(idx.pos) > 0) {
        x[idx.pos] = cac$pos$idist(z[idx.pos])
    }
    if (length(idx.neg) > 0) {
        x[idx.neg] = cac$neg$idist(z[idx.neg])
    }
    return (exp(x))
}

inla.pc.dalphaw = function(alpha, lambda = 7, log = FALSE)
{
    cac = inla.pc.alphaw.cache()
    d = numeric(length(alpha))
    lalpha = log(alpha)
    idx.pos = which(lalpha >= 0)
    idx.neg = which(lalpha <  0)
    if (length(idx.pos) > 0) {
        la = lalpha[idx.pos]
        dist = cac$pos$dist(la)
        deriv = cac$pos$dist(la, deriv = 1)
        d[idx.pos] = log(0.5) + log(lambda) - lambda * dist + log(abs(deriv)) - la
    }        
    if (length(idx.neg) > 0) {
        la = lalpha[idx.neg]
        dist = cac$neg$dist(la)
        deriv = cac$neg$dist(la, deriv = 1)
        d[idx.neg] = log(0.5) + log(lambda) - lambda * dist + log(abs(deriv)) - la
    }        
    return (if (log) d else exp(d))
}

inla.pc.qalphaw = function(p, lambda = 7)
{
    cac = inla.pc.alphaw.cache()
    n = length(p)
    q = numeric(n)
    idx.pos = which(p >= 0.5)
    idx.neg = which(p <  0.5)
    if (length(idx.pos) > 0) {
        pp = 2.0*(p[idx.pos] - 0.5)
        qe = qexp(pp, rate = lambda)
        q[idx.pos] = cac$pos$idist(qe)
    }
    if (length(idx.neg) > 0) {
        pp = 2.0*(0.5 - p[idx.neg])
        qe = qexp(pp, rate = lambda)
        q[idx.neg] = cac$neg$idist(qe)
    }
    return(exp(q))
}

inla.pc.palphaw = function(q, lambda = 7)
{
    cac = inla.pc.alphaw.cache()
    n = length(q)
    p = numeric(n)
    idx.pos = which(q >= 1.0)
    idx.neg = which(q <  1.0)
    if (length(idx.pos) > 0) {
        qq = cac$pos$dist(log(q[idx.pos]))
        pp = pexp(qq, rate = lambda)
        p[idx.pos] = 0.5 * (1 + pp)
    }
    if (length(idx.neg) > 0) {
        qq = cac$neg$dist(log(q[idx.neg]))
        pp = pexp(qq, rate = lambda)
        p[idx.neg] = 0.5 * (1 - pp)
    }
    return(p)
}

inla.pc.alphaw.test = function(lambda = 7) 
{
    n = 10^7
    x = inla.pc.ralphaw(n, lambda = lambda)
    x = x[x < quantile(x, prob = 0.999)]
    alpha = seq(min(x), max(x), by = 0.001)
    hist(x, n = 300,  prob = TRUE)
    lines(alpha, inla.pc.dalphaw(alpha, lambda), lwd=3, col="blue")
    print(sum(alpha * inla.pc.dalphaw(alpha, lambda) * diff(alpha)[1]))

    p = seq(0.1, 0.9, by = 0.05)
    print(cbind(p=p, q.emp = as.numeric(quantile(x,  prob = p)), q.true = inla.pc.qalphaw(p, lambda)))

    q = as.numeric(quantile(x, prob = p))
    cdf = ecdf(x)
    print(cbind(q=q, p.emp = cdf(q), p.true = inla.pc.palphaw(q, lambda)))
}

## inla.pc.alphaw.test()
