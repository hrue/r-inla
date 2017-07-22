## utility function for quantile-regression.
##     - continous Poisson and quantile regression for Poisson

inla.incGamma = function(x, lambda=0, log=FALSE) {
    ##library(gsl); return (gamma_inc(x, lambda))
    lres = lgamma(x) + pgamma(lambda, x, lower = FALSE, log.p=TRUE)
    return (if (log) lres else exp(lres))
}

inla.qcontpoisson = function(p, lambda, print.level = 0) {
    fun.min = function(log.x, lambda, prob) {
        p.est = inla.pcontpoisson(exp(log.x), lambda)
        eps = 1E-12
        p.est = max(eps, min(p.est, 1 - eps))
        return ((inla.link.logit(p.est) - inla.link.logit(prob))^2)
    }
    initial.value = log(qpois(p, lambda) + 0.5)
    return (exp(nlm(fun.min, p = initial.value, print.level = print.level,
                    lambda = lambda, prob = p)$estimate))
}

inla.pcontpoisson = function(x, lambda, log=FALSE, deriv=0) {
    ## > F := (x, lambda) -> GAMMA(x, lambda) / GAMMA(x);
    ##                                     GAMMA(x, lambda)
    ##                F := (x, lambda) -> ----------------
    ##                                         GAMMA(x)
    ##
    ## > simplify(diff(F(x, lambda), lambda));           
    ##                             (x - 1)
    ##                       lambda        exp(-lambda)
    ##                     - --------------------------
    ##                                GAMMA(x)
    ##
    ## > simplify(diff(F(x, lambda), lambda$2));         
    ##                   (x - 2)
    ##             lambda        exp(-lambda) (-x + 1 + lambda)
    ##             --------------------------------------------
    ##                               GAMMA(x)
    ##
    if (deriv == 0) {
        ## can use gsl::gamma_inc_Q() instead
        lres = inla.incGamma(x, lambda, log=TRUE) - inla.incGamma(x, log=TRUE)
        return (if (log) lres else exp(lres))
    } else if (deriv == 1) {
        stopifnot(!log)
        return (-exp((x-1)*log(lambda) -lambda -inla.incGamma(x, log=TRUE)))
    } else if (deriv == 2) {
        stopifnot(!log)
        return ((lambda-x+1) * exp((x-2)*log(lambda) -lambda -inla.incGamma(x, log=TRUE)))
    } else {
        stop("deriv != 0, 1, 2")
    }
}

inla.pcontpoisson.eta = function(x, eta, deriv = 0, log=FALSE) {
    ## the cdf of the contpoisson parameterised by the linear predictor and the log-link
    lambda = exp(eta)
    if (deriv == 0) {
        return (inla.pcontpoisson(x, lambda, log=log))
    } else if (deriv == 1) {
        stopifnot(!log)
        return (inla.pcontpoisson(x, lambda, deriv=1) * lambda)
    } else if (deriv == 2) {
        stopifnot(!log)
        return (lambda * (inla.pcontpoisson(x, lambda, deriv=1) +
                          inla.pcontpoisson(x, lambda, deriv=2) * lambda))
    } else {
        stop("deriv != 0, 1, 2")
    }
}

inla.contpoisson.solve.lambda = function(quantile, alpha, iter.max = 1000, max.step = 3,
                                           tol = sqrt(.Machine$double.eps), verbose=FALSE)
{
    ## solve quantile=inla.pcontpoisson(lambda, alpha),  for lambda
    stopifnot(length(quantile) == 1 && length(alpha) == 1)
    return(exp(inla.contpoisson.solve.eta(quantile, alpha, iter.max, max.step, tol, verbose)))
}

inla.contpoisson.solve.eta = function(quantile, alpha, iter.max = 1000, max.step = 3,
                                               tol = sqrt(.Machine$double.eps), verbose=FALSE) 
{
    ## solve quantile=inla.pcontpoisson(lambda=exp(eta), alpha),  for eta
    stopifnot(length(quantile) == 1 && length(alpha) == 1)
    eta.0 = log((sqrt(quantile) - qnorm(alpha,  sd = 0.5))^2)
    for(i in 1:iter.max) {
        f = inla.pcontpoisson(quantile, lambda = exp(eta.0)) - alpha
        fd = inla.pcontpoisson.eta(quantile, eta.0, deriv=1)
        d = -min(max.step, max(-max.step, f/fd))
        eta = eta.0 = eta.0 + d
        if (verbose)
            print(round(c(iter=i, eta = eta, f=f, f.deriv = fd, err = d), digits = 6))
        if (abs(d) < tol)
            return(eta)
    }
    stop(paste("FAILURE", quantile, alpha))
}
