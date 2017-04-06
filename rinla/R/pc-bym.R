## NOT-YET-Export: 


inla.bym.constr.internal = function(Q, adjust.for.con.comp = TRUE)
{
    ## return the rankdef and the constr for the BYM model given in Q

    n = dim(Q)[1]
    g = inla.read.graph(Q)
    cc.n = sapply(g$cc$nodes, length)
    cc.n1 = sum(cc.n == 1L)
    cc.n2 = sum(cc.n >= 2L)

    if (adjust.for.con.comp) {
        constr = list(A = matrix(0, cc.n2, n), e = rep(0, cc.n2))
        kk = 0
        for(k in which(cc.n >= 2L)) {
            kk = kk + 1
            constr$A[kk, g$cc$nodes[[k]]] = 1
        }
        stopifnot(kk == cc.n2)
        rankdef = cc.n2
    } else {
        nc = 1
        constr = list(A = matrix(1, nc, n), e = rep(1, nc))
        rankdef = cc.n1 + 1
    }

    return (list(rankdef = nrow(constr$A),
                 constr = constr, 
                 rankdef = rankdef,
                 cc.n = cc.n,
                 cc.n1 = cc.n1,
                 cc.n2 = cc.n2))
}

inla.sparse.det.bym = function(Q,
                               rankdef,
                               adjust.for.con.comp = TRUE,
                               log=TRUE,
                               constr = NULL, 
                               eps = 0.01 * sqrt(.Machine$double.eps))
{
    ## compute the (log-)determinant. Assume a sum to zero constraint on each connected
    ## component >1. if 'constr' is given,  it is supposed to be the 'constr' that is
    ## constructed below.

    stopifnot(adjust.for.con.comp == TRUE)
    Q = inla.as.sparse(Q)
    n = dim(Q)[1]
    res = NULL
    
    if (is.null(constr)) {
        res = inla.bym.constr.internal(Q, adjust.for.con.comp = adjust.for.con.comp)
        constr = res$constr
    }
    if (missing(rankdef)) {
        if (!is.null(constr)) {
            rankdef = nrow(constr$A)
        } else {
            rankdef = 0
        }
    }
    d = diag(Q)
    diag(Q) = d + mean(d) * eps
    res = inla.qsample(n=1, Q=Q, sample = matrix(0, n, 1), constr = constr)
    logdet = 2.0 * (res$logdens + (n-rankdef)/2.0*log(2.0*pi))
    return (if (log) logdet else exp(logdet))
}

inla.scale.model.bym = function(Q, eps = sqrt(.Machine$double.eps), adjust.for.con.comp = TRUE) 
{
    ## the old version
    res = inla.scale.model.bym.internal(Q=Q, adjust.for.con.comp = adjust.for.con.comp, eps = eps)
    return (res$Q)

}

inla.scale.model.bym.internal = function(Q, eps = sqrt(.Machine$double.eps), adjust.for.con.comp = TRUE) 
{
    ## the new one which also return the (scaled) marginal variances
    n = dim(Q)[1]
    constr = list(A = matrix(1, nrow=1, ncol=n), e=0)

    if (adjust.for.con.comp) {
        res = inla.scale.model.internal(Q, constr = constr, eps = eps)
    } else {
        mvar = rep(NA, n)
        idx = which(diag(Q) > 0)
        QQ = Q[idx, idx, drop=FALSE]
        QQ = inla.as.sparse(QQ)
        n = dim(QQ)[1]
        constr = list(A = matrix(1, nrow=1, ncol=n), e=0)
        res = inla.qinv(QQ + Diagonal(n) * mean(diag(QQ)) * eps, constr = constr)
        fac = exp(mean(log(diag(res))))
        QQ = fac * QQ
        Q[idx, idx] = QQ
        mvar[idx] = diag(res)/fac
        res = list(Q=Q, var = mvar)
    }
    return (res)        
}

inla.pc.bym.Q = function(graph)
{
    Q = -inla.graph2matrix(graph)
    diag(Q) = 0
    diag(Q) = -rowSums(Q)
    return (Q)
}

## also used for the rw2d+iid model (eigenvalues/marginal.variances arguments...)
inla.pc.bym.phi = function(graph,
                           Q, 
                           eigenvalues = NULL,
                           marginal.variances = NULL, 
                           rankdef,
                           alpha,
                           u = 1/2,
                           lambda,
                           scale.model = TRUE,
                           return.as.table = FALSE, 
                           adjust.for.con.comp = TRUE, 
                           ## where to switch to alternative strategy
                           eps = sqrt(.Machine$double.eps), 
                           debug = FALSE)
{
    my.debug = function(...) if (debug) cat("*** debug *** inla.pc.bym.phi: ", ... , "\n")

    ## I must assume this!!!
    stopifnot(scale.model == TRUE)
    stopifnot(adjust.for.con.comp == TRUE)
    res = NULL
    
    if (missing(eigenvalues) && missing(marginal.variances)) {
        ## computes the prior for the mixing parameter `phi'.
        if (missing(graph) && !missing(Q)) {
            Q = inla.as.sparse(Q)
        } else if (!missing(graph) && missing(Q)) {
            Q = inla.pc.bym.Q(graph)
        } else if (missing(graph) && missing(Q)) {
            stop("Either <Q> or <graph> must be given.")
        } else {
            stop("Only one of <Q> and <graph> can be given.")
        }
        
        res = inla.bym.constr.internal(Q, adjust.for.con.comp = adjust.for.con.comp)
        if (missing(rankdef)) {
            rankdef = res$rankdef
        }

        n = dim(Q)[1]
        res = inla.scale.model.bym.internal(Q, adjust.for.con.comp = adjust.for.con.comp)
        Q = res$Q
        f = mean(res$var, na.rm=TRUE) - 1.0 
        use.eigenvalues = FALSE
    } else {
        n = length(eigenvalues)
        eigenvalues = pmax(0, sort(eigenvalues, decreasing=TRUE))
        gamma.invm1= c(1/eigenvalues[1:(n-rankdef)], rep(0, rankdef)) - 1.0
        f = mean(marginal.variances) - 1.0
        use.eigenvalues = TRUE
    }
    
    if (use.eigenvalues) {
        ## this is fast for low dimension where we can compute the
        ## eigenvalues
        phi.s = 1/(1+exp(-seq(-15, 12,  len = 1000)))
        d = numeric(length(phi.s))
        k = 1
        for(phi in phi.s) {
            aa = n*phi*f
            ##bb = sum(log(1-phi + phi*gamma.inv))
            bb = sum(log1p(phi*gamma.invm1))
            ##bb = sum(log(1+phi*gamma.invm1))
            ## sometimes we *might* experience numerical problems, so we need
            ## to remove (if any) small phi's.
            d[k] = (if (aa >= bb) sqrt(aa -bb) else NA)
            k = k + 1
        }
    } else {
        ## alternative strategy for larger matrices
        phi.s = 1/(1+exp(-c(seq(-15, 0, len=40), 1:12)))
        d = numeric(length(phi.s))
        log.q1.det = inla.sparse.det.bym(Q, adjust.for.con.comp = adjust.for.con.comp,
                                         constr = res$constr, rankdef = rankdef)
        d = c(unlist(
            inla.mclapply(phi.s, 
                          function(phi) {
                aa = n*phi*f
                bb = (n * log((1.0 - phi)/phi) +
                      inla.sparse.det.bym(Q + phi/(1-phi) * Diagonal(n),  
                                          adjust.for.con.comp = adjust.for.con.comp,
                                          constr = res$constr, rankdef = rankdef) -
                      (log.q1.det - n * log(phi)))
                my.debug("aa=", aa, " bb=", bb, " sqrt(", aa-bb, ")")
                return (if (aa >= bb) sqrt(aa -bb) else NA)
            }
            )
        ))
    }
    names(d) = NULL

    ## remove phi's that failed, if any
    remove = is.na(d)
    my.debug(paste(sum(remove), "out of", length(remove), "removed"))
    d = d[!remove]
    phi.s = phi.s[!remove]
    remove = 1:6
    d = d[-remove]
    phi.s = phi.s[-remove]
    phi.intern = log(phi.s/(1-phi.s))
    ff.d = splinefun(phi.intern, log(d))
    phi.intern = seq(min(phi.intern)-0, max(phi.intern), len = 10000)
    phi.s = 1/(1+exp(-phi.intern))
    d = exp(ff.d(phi.intern))
    
    ## ff.d:  return distance as a function of phi.intern
    ff.d.core = splinefun(phi.intern, log(d))
    ff.d = function(phi.intern, deriv=0) {
        if (deriv == 0) {
            return (exp(ff.d.core(phi.intern)))
        } else if (deriv == 1) {
            return (exp(ff.d.core(phi.intern)) * ff.d.core(phi.intern,  deriv=1))
        } else {
            stop("ERROR")
        }
    }
    ## f.d: return distance as a function of phi
    f.d = function(phi.s) ff.d(log(phi.s/(1-phi.s)))
    d = f.d(phi.s)
    if (missing(lambda)) {
        ## Prob(phi < u) = alpha gives an analytical solution
        stopifnot(alpha > 0.0 && alpha < 1.0 && u > 0)
        lambda = -log(1-alpha)/f.d(u)
    }

    ## this is the log.prior wrt phi.intern
    log.jac = log(abs(ff.d(phi.intern, deriv=1))) 
    log.prior = log(lambda) - lambda * d + log.jac
    ssum = sum(exp(log.prior) * (c(0, diff(phi.intern)) + c(diff(phi.intern), 0)))/2.0
    my.debug("empirical integral of the prior = ", ssum, ". correct it.")
    log.prior = log.prior - log(ssum)

    if (return.as.table) {
        ## return this prior as a "table:" converted to the intern scale for phi.
        my.debug("return prior as a table on phi.intern")
        theta.prior.table = paste(c("table:", cbind(phi.intern, log.prior)),
                                  sep = "", collapse = " ")
        return (theta.prior.table)
    } else {
        ## return a function which evaluates the log-prior as a function of phi.
        my.debug("return prior as a function of phi")
        ## do the interpolation in the intern scale
        prior.theta.1 = splinefun(phi.intern, log.prior +phi.intern +2*log(1+exp(-phi.intern)))
        prior.phi = function(phi) prior.theta.1(log(phi/(1.0-phi)))
        return (prior.phi)
    }
}

inla.pc.rw2diid.phi = function(
        size, 
        alpha,
        u = 1/2,
        lambda,
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
        return.as.table = return.as.table, 
        debug = debug,
        alpha = alpha,
        u = u,
        lambda = lambda, 
        rankdef = 1L))
}


