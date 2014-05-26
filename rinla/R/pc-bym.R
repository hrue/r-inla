## NOT-YET-Export: 

inla.pc.bym.Q = function(graph, rankdef = 1, scale.model = FALSE)
{
    Q = -inla.graph2matrix(graph)
    diag(Q) = 0
    diag(Q) = -rowSums(Q)
    if (scale.model) {
        fac = exp(mean(log(diag(INLA:::inla.ginv(
                x=as.matrix(Q), rankdef=rankdef)))))
        Q = fac * Q
    }
    return (Q)
}

inla.pc.bym.phi = function(
        graph,
        Q, 
        eigenvalues = NULL,
        marginal.variances = NULL, 
        rankdef = 1L,
        ## if alpha is not set,  it will be computed from alpha.min
        alpha,
        u = 1/2,
        lambda,
        ## use with care, as this only applied to the case where Q is
        ## already scaled.
        scale.model = TRUE,
        phi.s = 1/(1+exp(-seq(-12, 12,  len = 10000))),
        return.as.table = FALSE, 
        debug = FALSE)
{
    my.debug = function(...) if (debug) cat("*** debug *** inla.pc.bym.phi: ", ... , "\n")
    
    if (missing(eigenvalues) && missing(marginal.variances)) {
        ## computes the prior for the mixing parameter `phi'.
        if (missing(graph) && !missing(Q)) {
            Q = as.matrix(Q)
        } else if (!missing(graph) && missing(Q)) {
            Q = inla.pc.bym.Q(graph, rankdef)
        } else if (missing(graph) && missing(Q)) {
            stop("Either <Q> or <graph> must be given.")
        } else {
            stop("Only one of <Q> and <graph> can be given.")
        }
        
        if (scale.model) {
            fac = exp(mean(log(diag(
                    INLA:::inla.ginv(x=as.matrix(Q),
                                     rankdef=rankdef)))))
            Q = fac * Q
        }
        Q = as.matrix(Q)
        
        n = dim(Q)[1]
        e = eigen(Q)
        gamma.inv = c(1/e$values[1:(n-rankdef)], rep(0, rankdef))
        Qinv.d = c(diag(e$vectors %*% diag(gamma.inv) %*% t(e$vectors)))
        f = mean(Qinv.d)-1
    } else {
        n = length(eigenvalues)
        gamma.inv = c(1/eigenvalues[1:(n-rankdef)], rep(0, rankdef))
        if (length(marginal.variances) == 1) {
            Qinv.d = rep(marginal.variances,  n)
        } else {
            Qinv.d = marginal.variances
        }
        f = mean(Qinv.d)-1
    }

    d = numeric(length(phi.s))
    k = 1
    for(phi in phi.s) {
        aa = n*phi*f
        bb = sum(log(1-phi + phi*gamma.inv))
        d[k] = sqrt(aa - bb)
        k = k + 1
    }
    phi.s = c(0, phi.s)
    d = c(0, d)
    f.d = splinefun(phi.s, d)
    d.max = f.d(1.0)
    
    if (missing(lambda)) {
        ## Prob(phi < u) = alpha
        alpha.min = f.d(u)/d.max
        my.debug("compute alpha.min = ", alpha.min)
        if (missing(alpha) || (alpha <= 0.0 || alpha >= 1.0)) {
            ## set the default value a little over the minimum one.
            fac = 0.975
            alpha = alpha.min * fac + 1 * (1-fac)
            my.debug("missing(alpha) == TRUE. set alpha = ", alpha, "using u=", u)
        }
        ## initial value
        if (!(alpha > alpha.min)) {
            my.debug("This must be true: alpha > ", alpha.min)
            stopifnot(alpha > alpha.min)
        }

        f.target = function(lam, alpha, d.u, d.max) {
            return ((1-exp(-lam*d.u)/(1-exp(-lam*d.max)) - alpha)^2)
        }

        lam = exp(seq(-5, 3, len=1000))
        idx = which.min(f.target(lam, alpha, f.d(u), d.max))
        r = optimise(f.target, interval = c(lam[idx-1], lam[idx+1]),
                maximum = FALSE,
                alpha = alpha, d.u = f.d(u), d.max = d.max)
        stopifnot(abs(r$objective) < 1e-5)
        lambda = r$minimum
        my.debug("found lambda = ", lambda)
    }
    
    h = min(min(1e-5, min(phi.s[phi.s > 0])), 1-max(phi.s[phi.s<1]))
    deriv = (f.d(phi.s + h) - f.d(phi.s - h))/(2*h)
    jac = abs(deriv)
    log.prior = log(lambda) - lambda * d + log(jac) - log(1-exp(-lambda * f.d(1)))

    if (return.as.table) {
        ## return this prior as a "table:" converted to the internal
        ## scale for phi.
        my.debug("return prior as a table on phi.internal")
        legal.idx = (phi.s > 0 & phi.s < 1)
        phi.s = phi.s[legal.idx]
        log.prior = log.prior[legal.idx]
        prior.phi = splinefun(phi.s, log.prior)
        theta = log(phi.s/(1-phi.s))
        log.prior.theta = prior.phi(phi.s) + log(exp(-theta)/(1+exp(-theta))^2)
        theta.prior.table = paste(c("table:", cbind(theta, log.prior.theta)),
                sep = "", collapse = " ")
        return (theta.prior.table)
    } else {
        ## return a function which evaluates the log-prior as a
        ## function of phi.
        my.debug("return prior as a function of phi")
        prior.phi = splinefun(phi.s, log.prior)
        return (prior.phi)
    }
}
