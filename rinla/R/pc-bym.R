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

## also used for the rw2d+iid model (eigenvalues/marginal.variances
## arguments...)
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
        phi.s = 1/(1+exp(-seq(-12, 12,  len = 1000))),
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
                    INLA:::inla.ginv(x=as.matrix(Q), rankdef=rankdef)))))
            Q = fac * Q
        }
        Q = as.matrix(Q)
        
        n = dim(Q)[1]
        e = eigen(Q)
        gamma.inv = c(1/e$values[1:(n-rankdef)], rep(0, rankdef))
        gamma.invm1 = gamma.inv - 1.0
        Qinv.d = c(diag(e$vectors %*% diag(gamma.inv) %*% t(e$vectors)))
        f = mean(Qinv.d)-1
    } else {
        n = length(eigenvalues)
        eigenvalues = pmax(0, sort(eigenvalues, decreasing=TRUE))
        gamma.invm1= c(1/eigenvalues[1:(n-rankdef)], rep(0, rankdef)) - 1.0
        f = mean(marginal.variances) - 1.0
    }

    d = numeric(length(phi.s))
    k = 1
    for(phi in phi.s) {
        aa = n*phi*f
        ##bb = sum(log(1-phi + phi*gamma.inv))
        bb = sum(log1p(phi*gamma.invm1))
        ##bb = sum(log(1+phi*gamma.invm1))
        ## sometimes we *might* experience numerical problems, so we need
        ## to remove (if any) small phi's.
        if (aa >= bb) {
            d[k] = sqrt(aa - bb)
        } else {
            d[k] = NA
        }
        k = k + 1
    }
    ## remove phi's that failed
    remove = is.na(d)
    d = d[!remove]
    phi.s = phi.s[!remove]
    ##
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
            value = ((1-exp(-lam*d.u)/(1-exp(-lam*d.max)) - alpha)^2)
            ##my.debug("lam=", lam, "value=", value)
            return (value)
        }

        ## do the screning in two rounds
        lam = exp(seq(-10, 3, len=25))
        idx = which.min(f.target(lam, alpha, f.d(u), d.max))
        lam = exp(seq(log(lam[idx-1]), log(lam[idx+1]), len=25))
        idx = which.min(f.target(lam, alpha, f.d(u), d.max))
        ## then call the optimiser...
        r = optimise(f.target, interval = c(lam[idx-1], lam[idx+1]),
                maximum = FALSE,
                alpha = alpha, d.u = f.d(u), d.max = d.max)
        stopifnot(abs(r$objective) < 1e-3)
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

inla.pc.rw2diid.phi = function(
        size, 
        alpha,
        u = 1/2,
        lambda,
        phi.s = 1/(1+exp(-seq(-12, 12,  len = 1000))),
        return.as.table = FALSE, 
        debug = FALSE)
{
    ## not all functions here are used, but...
    DFT2 = function(x)
    {
        fft(x,inverse=FALSE)/sqrt(length(x))
    }

    IDFT2 = function(x)
    {
        fft(x,inverse=TRUE)/sqrt(length(x))
    }

    make.base.1 = function(size, delta = 0)
    {
        if (length(size) == 1) size = c(size, size)
        base      = matrix(data=0, nrow=size[1], ncol=size[2])
        base[1,1] = 4 + delta
        base[1,2] = base[2,1] = base[size[1],1] = base[1,size[2]] = -1
        return(base)
    }

    eigenvalues = function(base)
    {
        return (Re(sqrt(length(base))*DFT2(base)))
    }

    eigenvalues.rw2d = function(base)
    {
        n = dim(base)[1]
        outer(0:(n-1), 0:(n-1),
              function(i, j) {
                  two.pi = 2*pi
                  return (base[1, 1]
                          + base[2, 1] * cos(two.pi * i / n) 
                          + base[1, 2] * cos(two.pi * j / n) 
                          + base[n, n] * cos(two.pi * ((n-1) * i / n + (n-1) * j / n)) 
                          + base[n, 1] * cos(two.pi * (n-1) * i / n) 
                          + base[1, n] * cos(two.pi * (n-1) * j / n))
              })^2
    }

    make.base = function(size, delta = 0)
    {
        base = make.base.1(size, delta)
        base2 = Re(IDFT2(DFT2(base)^2)) * sqrt(length(base))
        return(base2)
    }

    inverse = function(base, rankdef, factor = 100)
    {
        eig = Re(DFT2(base))
        eig.max = max(eig)
        if (missing(rankdef)) {
            pos = eig > eig.max * factor * .Machine$double.eps
        } else {
            eig.sort = sort(eig)
            if (rankdef > 0) {
                pos = eig > eig.sort[rankdef]
            } else {
                pos = eig > 0
            }
        }
        eig[pos] = 1/eig[pos]
        eig[!pos] = 0
        base.inv = Re(IDFT2(eig))/length(base)
        return (base.inv)
    }

    tozero = function(base, factor = 100)
    {
        elm.max = max(abs(base))
        base[abs(base) < factor*.Machine$double.eps*elm.max] = 0
        return(base)
    }
    
    mult = function(base, bbase)
    {
        a = Re(DFT2(base))
        b = Re(DFT2(bbase))
        ab = a*b
        xx = Re(IDFT2(ab)) * sqrt(length(base))
        return(xx)
    }
    
    check.inla = function(size)
    {
        if (length(size) == 1) size = c(size, size)
        n = prod(size)
        y = numeric(n)
        y[] = 0
        idx = 1:n
        r = inla(y ~ -1 + f(idx, model="rw2d", nrow = size[1], ncol=size[2], scale.model=TRUE, fixed=TRUE),
                data = data.frame(y, idx),
                family = "gaussian",
                control.family = list(hyper = list(prec = list(initial = 8, fixed=TRUE))))
        return(r$logfile[ grep(" prec_scale", r$logfile) ])
    }

    if (length(size) == 1) {
        size = rep(size, 2)
    }
    size = size * 1.55

    ## numbers of the form 2^i*3^j*5^k
    good.size = c(8, 9, 10, 12, 15, 16, 18, 20, 24, 25, 27, 30, 32,
            36, 40, 45, 48, 50, 54, 60, 64, 72, 75, 80, 81, 90, 96,
            100, 108, 120, 125, 128, 135, 144, 150, 160, 162, 180,
            192, 200, 216, 225, 240, 243, 250, 256, 270, 288, 300,
            320, 324, 360, 375, 384, 400, 405, 432, 450, 480, 486,
            500, 512, 540, 576, 600, 625, 640, 648, 675, 720, 729,
            750, 768, 800, 810, 864, 900, 960, 972, 1000, 1024, 1080,
            1125, 1152, 1200, 1215, 1250, 1280, 1296, 1350, 1440,
            1458, 1500, 1536, 1600, 1620, 1728, 1800, 1875, 1920,
            1944, 2000, 2025, 2048, 2160, 2187, 2250, 2304, 2400,
            2430, 2500, 2560, 2592, 2700, 2880, 2916, 3000, 3072,
            3125, 3200, 3240, 3375, 3456, 3600, 3645, 3750, 3840,
            3888, 4000, 4050, 4096)

    for(k in 1:length(size)) {
        size[k] = good.size[min(which(size[k] <= good.size))]
    }

    base = make.base(size)
    base = base * inverse(base, rankdef = 1L)[1, 1]
    marg.var = 1.0
    e.values = pmax(0, eigenvalues(base))

    return (inla.pc.bym.phi(
        eigenvalues = c(e.values),
        marginal.variances = marg.var,
        phi.s = phi.s,
        return.as.table = return.as.table, 
        debug = debug,
        alpha = alpha,
        u = u,
        lambda = lambda, 
        rankdef = 1L))
}
