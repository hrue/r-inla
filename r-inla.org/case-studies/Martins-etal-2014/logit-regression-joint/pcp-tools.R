library("numDeriv")

my.logdet = function(A) {
    return (determinant(A, log=TRUE)$modulus)
}

kld = function(S, S.base, solve.S.base, log.det.S.base)
{
    ## only compute these if they are missing
    if (missing(solve.S.base))
        solve.S.base = solve(S.base)
    if (missing(log.det.S.base))
        log.det.S.base = my.logdet(S.base)

    return (0.5*(sum(diag(solve.S.base %*% S)) - dim(S)[1]
                 - (my.logdet(S) - log.det.S.base)))
}

dist = function(S, S.base, solve.S.base, log.det.S.base)
{
    return (sqrt(max(0, 2*kld(S, S.base, solve.S.base, log.det.S.base))))
}

pcp.case1 = function(S.base, S.alt, U = 0.5, alpha = 0.5,
    use.phi.intern = FALSE, debug = FALSE)
{
    ## case 1.
    ## S(phi) = (1-phi) * S.base + phi * S.alt
    ## S(0) = S.base
    ## return a function which computes the log-pcp.
    ## Prob(phi < U) = alpha
    my.debug = function(...) if (debug) 
        cat("*** debug *** pcp.case1: ", ..., "\n")
 
    map = function(intern) 1/(1+exp(-intern))
    map.inv = function(phi) log(phi/(1-phi))
    n = 25
    lim = 10
    phis.intern = seq(-lim, lim, len=n)
    phis = map(phis.intern)
    d = numeric(n)

    solve.S.base = solve(S.base)
    log.det.S.base = my.logdet(S.base)
    for(i in 1:n) {
        phi = phis[i]
        d[i] = dist((1-phi) * S.base + phi * S.alt, S.base,
             solve.S.base, log.det.S.base)
    }

    s1 = spline(phis.intern, d, xmin = -lim, xmax = lim, n = n)
    f1 = splinefun(s1$x, s1$y)
    dfunc = function(phi) {
        intern = map.inv(phi)
        return (f1(intern))
    }

    if (FALSE) {
        ## Prob(phi < U) = alpha
        lambda = - log(1-alpha)/dfunc(U)
        print(paste("pcp.case1: guess lambda", lambda))
    }

    d.max = max(d)
    phis.intern.max = max(phis)
    alpha.min = dfunc(U)/d.max
    my.debug("compute alpha.min = ", alpha.min)
    if (missing(alpha) || (alpha <= 0 || alpha >= 1)) {
        fac = 0.975
        alpha = alpha.min * fac + 1 * (1 - fac)
        my.debug("missing(alpha) == TRUE. set alpha = ", 
                 alpha, "using u=", U)
    }
    if (!(alpha > alpha.min)) {
        my.debug("This must be true: alpha > ", alpha.min)
        print(paste("alpha", alpha, "alpha.min", alpha.min, "must have alpha > alpha.min"))
        stopifnot(alpha > alpha.min)
    }
    f.target = function(lam, alpha, d.u, d.max) {
        value = ((1 - exp(-lam * d.u)/(1 - exp(-lam * d.max)) - 
                  alpha)^2)
        return(value)
    }
    lam = exp(seq(-10, 3, len = 25))
    idx = which.min(f.target(lam, alpha, dfunc(U), d.max))
    lam = exp(seq(log(lam[idx - 1]), log(lam[idx + 1]), len = 25))
    idx = which.min(f.target(lam, alpha, dfunc(U), d.max))
    r = (optimise(f.target,
                  interval = c(lam[idx - 1], lam[idx + 1]),
                  maximum = FALSE, alpha = alpha, d.u = dfunc(U), 
                  d.max = d.max))
    stopifnot(abs(r$objective) < 0.001)
    lambda = r$minimum
    my.debug("found lambda = ", lambda)
    print(paste("pcp.case1: find lambda", lambda))

    ## this function do one numerical diff without touching the boundary
    grad = function(fun, x, h = 1e-4, xlim = c(0, 1))
    {
        n = length(x)
        g = numeric(n)
        for(i in 1:n) {
            if ((x[i] > range(xlim)[1] + h) &&
                (x[i] < range(xlim)[2] - h)) {
                g[i] = (fun(x[i]+h) - fun(x[i]-h))/(2*h)
            } else if (x[i] < mean(xlim)) {
                g[i] = (fun(x[i]+h)-fun(x[i]))/h
            } else {
                g[i] = (fun(x[i])-fun(x[i]-h))/h
            }
        }
        return (g)
    }

    prior = function(phi, lambda, xlim, use.phi.intern = FALSE, d.max) {
        ldens = log(lambda) - lambda * dfunc(phi) - log(1-exp(-lambda * d.max))
        if (use.phi.intern) {
            return (ldens)
        } else {
            return (ldens + log(abs(grad(dfunc, phi, xlim = xlim))))
        }
    }

    if (use.phi.intern) {
        phis = map(seq(-lim, lim, len=n*10))
        prior.fun = splinefun(map.inv(phis), prior(phis, lambda, xlim = range(phis),
            use.phi.intern =TRUE, d.max = d.max))
    } else {
        phis = map(seq(-lim, lim, len=n*10))
        prior.fun = splinefun(phis, prior(phis, lambda, xlim = range(phis), d.max = d.max))
    }

    return (prior.fun)
}

## convert between phi[i] = log(w[i]/w[m]) and vica versa. this is the good parameterisation for
## the weights and is widely used for compositional data.
w2phi = function(w)
{
    m = length(w)
    phi = log(w[1:(m-1)]/w[m])
    return (phi)
}
phi2w = function(phi)
{
    n = length(phi)
    m = n + 1
    w = numeric(m)
    w[m] = 1/(1 + sum(exp(phi)))
    w[1:n] = w[m] * exp(phi)
    return (w)
}

pcp.case2.H = function(S, eps = 1e-4) 
{
    ## case 2
    ## S(w) = \sum_i w[i] S[[i]]
    ## S.base = S(w = 1/length(w)) 
    ## return the Hessian to be used to compute the pcp
    ## using the sphere approximation using the 'phi' parameterisation
    m = length(S)
    n = m-1

    S.base = S[[1]]
    for(i in 2:m) {
        S.base = S.base + S[[i]]
    }
    S.base = S.base/m
    solve.S.base = solve(S.base)
    log.det.S.base = my.logdet(S.base)

    ## create a function computing KLD as a function of the phi-parameters
    kld.case2 = function(phi, SS, S.base, solve.S.base, log.det.S.base)
    {
        n = length(phi)
        m = n + 1
        stopifnot(length(SS) == m)
        ## parameters phi[i] = log(w[i]/w[m]), i=1...m-1
        w = phi2w(phi)
        S = w[1] * SS[[1]]
        for(i in 2:m) {
            S = S + w[i] * SS[[i]]
        }
        return (kld(S, S.base, solve.S.base, log.det.S.base))
    }

    H = hessian(kld.case2, 
        ## the base-model
        x = rep(0, n), 
        ## args to my.kld
        SS = S, 
        S.base = S.base,
        solve.S.base = solve.S.base,
        log.det.S.base = log.det.S.base,
        ## default args to hessian: eps=1e-4. the rest is copied from the
        ## numDeriv:::hessian.default code.
        method.args = list(
            eps = eps, d = 0.1,
            zero.tol = sqrt(.Machine$double.eps/7e-07),
            r = 4, v = 2, show.details = FALSE))

    if (!all(eigen(H, only.values=TRUE)$values > 0)) {
        warning("H is not positive definite.")
    }

    return (H)
}


if (FALSE) {
    
    n = 300
    S.base = matrix(1, n, n) + 1e-6 * diag(n)
    S.alt = INLA:::inla.ginv(
        INLA:::inla.rw(n, order = 2, scale.model=TRUE, sparse=FALSE))

    pcp = pcp.case1(S.base, S.alt, U = 0.5, alpha = 0.5)

    phis = 1/(1+exp(-seq(-5, 5, len=1000)))
    plot(phis,  exp(pcp(phis)),  type = "l", lwd=3)
}

if (FALSE) {
    make.A = function(n, m)
    {
        A = matrix(0, n, m)
        s = sample(1:m, n, replace=TRUE)
        for(i in 1:n) {
            A[i, s[i]]=1
        }
        return (A)
    }
    
    m = 3
    n.rw = 30
    S.rw = INLA:::inla.ginv(
        INLA:::inla.rw(n.rw, order = 2, scale.model=TRUE, sparse=FALSE))
    SS = list()
    n = 50
    
    for(i in 1:m) {
        A = make.A(n, n.rw)
        SS[[i]] = A %*% S.rw %*% t(A) + 1e-4 * diag(n)
    }
    
    H = pcp.case2.H(SS)
    x = inla.pc.multvar.sphere.general.r(10000, lambda = 1, H=H)
    w = t(apply(x, 1, phi2w))
    hist(w[, ], n=1000, prob=TRUE)
}

ensure.symmetry = function(X)
{
    return ((X + t(X))/2.0)
}
