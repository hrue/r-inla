require(parallel)
require(INLA)

save.and.quit = FALSE
plot.figures = TRUE
write.new.figures = TRUE
new.windows = TRUE

my.grad = function(func, x)
{
    eps = 1e-6
    return ((func(x + eps)-func(x))/eps)
    
    require(numDeriv)
    return (grad(func, x))
}

rt.scaled = function(n, dof)
{
    sigma = sqrt(dof/(dof-2))
    x = rt(n, df=dof)
    return (x/sigma)
}

dt.scaled = function(x, dof, log=FALSE)
{
    sigma = sqrt(dof/(dof-2))
    if (!log) {
        return (dt(x*sigma, df=dof, log=FALSE)*sigma)
    } else {
        return (dt(x*sigma, df=dof, log=TRUE) + log(sigma))
    }
}

kld.approx = function(dof)
{
    stopifnot(length(dof) == 1L)

    ## plain version
    ##return (1/dof^2 * (3/4 + 1/dof * (3/2 + 1/dof * (17/8 + 1/dof * (29/10 + 1/dof * 61/12  etc....)))))
    df = dof
    t1 = df * df;
    t4 = t1 * df;
    t7 = t1 * t1;
    t10 = t7 * df;
    t13 = t7 * t1;
    t16 = t7 * t4;
    t19 = t7 * t7;
    t43 = t19 * t19;
    t60 = 0;
    t61 = 0.7500000000000000e0 / t1 + 0.1500000000000000e1 / t4 + 0.2125000000000000e1 / t7 + 0.2900000000000000e1 / t10 + 0.5083333333333333e1 / t13 + 0.1035714285714286e2 / t16 + 0.1706250000000000e2 / t19 + 0.1983333333333333e2 / t19 / df + 0.4345000000000000e2 / t19 / t1 + 0.1873181818181818e3 / t19 / t4 + 0.2570416666666667e3 / t19 / t7 - 0.1155192307692308e4 / t19 / t10 - 0.7801071428571429e3 / t19 / t13 + 0.3207790000000000e5 / t19 / t16 + 0.3109703125000000e5 / t43 - 0.8438102058823529e6 / t43 / df - 0.7932909722222222e6 / t43 / t1 + 0.2921518928947368e8 / t43 / t4 + 0.2776753702500000e8 / t43 / t7 + t60;
    return(t61)
}

kld = function(dof, lim=30)
{
    stopifnot(length(dof) == 1L)
    if (dof > 9) {
        return (kld.approx(dof))
    }

    func = function(x, dof) {
        dlog = dt.scaled(x, dof=dof, log=TRUE)
        dnlog = dnorm(x, log=TRUE)
        return (exp(dlog)*(dlog - dnlog))
    }

    s = sqrt(dof/(dof-2))
    fac = 2*s
    ret = 2*integrate(func, 0.0, lim*fac, dof=dof,
            subdivisions = 1000000L,
            stop.on.error = FALSE, 
            rel.tol = 10*.Machine$double.eps)$value
    return (ret)
}

## testing
stopifnot(abs(2*integrate(dt.scaled, 0, 20, dof=1000+runif(1))$value-1) < 0.0001) ## integral is 1
stopifnot(abs(2*integrate( ## variance is 1 for high dof
    function(x, dof) x^2*dt.scaled(x,dof=dof),
                        0, 100, 
                        subdivisions = 1000000L,
                        stop.on.error = FALSE, 
                        rel.tol = 100*.Machine$double.eps, 
              dof=100+runif(1))$value-1) < 0.0001) 
stopifnot(abs(2*integrate( ## variance is 1 for high dof
    function(x, dof) x^2*dt.scaled(x,dof=dof),
                          0, 100, 
                        subdivisions = 100000L,
                        stop.on.error = FALSE, 
                        rel.tol = 100*.Machine$double.eps, 
              dof=3.01)$value-1) < 0.02) 
stopifnot(kld(1000) < 0.0001) ## goes to zero

n = 10000
dofs = 2 + exp(seq(-6, 15,  len=n))
##kl = sapply(dofs, kld)
kl = unlist(mclapply(dofs, kld, mc.cores = detectCores()))
stopifnot(all(diff(kl) < 0.0))

bty = "l"
lwd = 2L
cex.lab = 1.4
cex.axis = 1.4
idx = which(dofs < 1000)
if (plot.figures) {
    if (new.windows) inla.dev.new()
    par(cex.lab=cex.lab, cex.axis=cex.axis)
    plot(dofs[idx], kl[idx], log="y", xlab="Degrees of freedom",  ylab = "KLD",
         type = "l",  lwd=lwd, bty = bty)
    if (write.new.figures) {
        dev.print(postscript,  file="tkld.ps")
    }
}

d = sqrt(2*kl)
## internal ones
f..ld = splinefun(log(dofs-2), log(d), method="hyman")
f..ldof = splinefun(log(d), log(dofs-2), method = "hyman")
## user
f.d = function(dof) exp(f..ld(log(dof-2)))
f.dof = function(d) 2+exp(f..ldof(log(d)))

stopifnot(abs(f.dof(f.d(3.123)) - 3.123) < 0.0001)

new.prior = function(dofs, u=NULL, p=NULL, lambda)
{
    stopifnot(!missing(dofs))
    x = dofs
    if (missing(lambda)) {
        ## our new prior
        lambda = -log(p)/f.d(u)
        y = dexp(f.d(x), rate = lambda) * abs(my.grad(f.d, x))
    } else {
        ## the exponential prior
        y = dexp(dofs-2.0, rate = lambda)
    }
    ## the prior is already scaled, but we scale it on dofs
    y = y / sum(y * c(diff(dofs), 0))
    idx = which(y > 0)
    return (list(x=x[idx], y=y[idx], f=splinefun(x[idx], y[idx]), u=u, p=p, lambda = lambda))
}
        
if (plot.figures) {
    if (new.windows) inla.dev.new()
    par(cex.lab=cex.lab, cex.axis=cex.axis)
    df.star = 10
    p.star = c(0.2, 0.5, 0.8)
    k = 1
    for(k in seq_along(p.star)) {
        np = new.prior(dofs, df.star, p.star[k])
        x = np$x
        y = np$y
        a = (which.max(y) > 1)
        if (k == 1) {
            plot(INLA:::inla.ifelse(a, c(2, x), x),
                 INLA:::inla.ifelse(a, c(0, y), y),
                 type="l", lwd=lwd, lty=1,
                 xlim = c(0, 80),
                 ylim = c(0, .2), 
                 xlab = "Degrees of freedom",
                 ylab = "Density",
                 bty = bty)
        } else {
            lines(
                INLA:::inla.ifelse(a, c(2, x), x),
                INLA:::inla.ifelse(a, c(0, y), y),
                lty=k, lwd=lwd)
        }
        k = k + 1
    }
    if (write.new.figures) {
        dev.print(postscript, file="tpriors.ps")
    }
}

if (plot.figures) {
    if (new.windows) inla.dev.new()
    par(cex.lab=cex.lab, cex.axis=cex.axis)
    lambda = 1/(5-2)
    x = f.d(dofs)
    xlim = c(0, 1.5)
    ylim = c(0, 8)
    y = dexp(dofs-2, rate = lambda) * abs(my.grad(f.dof, x))
    plot(x, y, type="l", lwd=lwd, lty=1,
         xlim = xlim, 
         ylim = ylim, , 
         xlab = "Distance",
         ylab = "Density",
         bty = bty)

    idx = which(x > xlim[1] & x < xlim[2])
    xx = x[idx]
    yy = y[idx]
    print(paste("mean 5: DOF at max", f.dof(xx[which.max(yy)])))

    lambda = (1/(10-2))
    y = dexp(dofs-2, rate = lambda) * abs(my.grad(f.dof, x))
    lines(x, y, lwd=lwd, lty=2)

    idx = which(x > xlim[1] & x < xlim[2])
    xx = x[idx]
    yy = y[idx]
    print(paste("mean 10: DOF at max", f.dof(xx[which.max(yy)])))

    lambda = (1/(20-2))
    y = dexp(dofs-2, rate = lambda) * abs(my.grad(f.dof, x))
    lines(x, y, lwd=lwd, lty=3)

    idx = which(x > xlim[1] & x < xlim[2])
    xx = x[idx]
    yy = y[idx]
    print(paste("mean 20: DOF at max", f.dof(xx[which.max(yy)])))

    if (write.new.figures) {
        dev.print(postscript,  file="tpriorsexp.ps")
    }

    if (new.windows) inla.dev.new()
    par(cex.lab=cex.lab, cex.axis=cex.axis)
    xlim = c(0, 0.3)
    ylim = c(0, 80)
    y = dunif(dofs, min = 2, max = 20) * abs(my.grad(f.dof, x))
    plot(x, y, type="l", lwd=lwd, lty=1,
         xlim = xlim, 
         ylim = ylim, , 
         xlab = "Distance",
         ylab = "Density",
         bty = bty)

    y = dunif(dofs, min = 2, max = 50) * abs(my.grad(f.dof, x))
    lines(x, y, lwd=lwd, lty=2)
    y = dunif(dofs, min = 2, max = 100) * abs(my.grad(f.dof, x))
    lines(x, y, lwd=lwd, lty=3)

    if (write.new.figures) {
        dev.print(postscript,  file="tpriorsunif.ps")
    }
}


run.sim = function(
        pcode = 1:11,
        dofs = 2 + exp(seq(-3, 10, len=10000)), 
        n = 100,
        n.sim = 100,
        dof.true = 100,
        llik.matrix)
{
    stopifnot(length(pcode) == 1)
    stopifnot(pcode %in% 1:11)
    ## do the simulation experiments
    post.sum = numeric(length(dofs))
    for(k in 1:n.sim) {
        print(k)
        if (missing(llik.matrix)) {
            y = rt.scaled(n, dof=dof.true)
            llik = sapply(dofs, function(dof, y) sum(dt.scaled(y, dof=dof, log=TRUE)), y=y)
        } else {
            stopifnot(is.matrix(llik.matrix) && nrow(llik.matrix) == n.sim && ncol(llik.matrix) == length(dofs))
            llik = as.vector(llik.matrix[k, ])
        }
        lik = exp(llik - max(llik))
        ## new prior
        if (pcode == 1) {
            prior = new.prior(dofs, u=10, p=0.2)
        } else if (pcode == 2) {
            prior = new.prior(dofs, u=10, p=0.3)
        } else if (pcode == 3) {
            prior = new.prior(dofs, u=10, p=0.4)
        } else if (pcode == 4) {
            prior = new.prior(dofs, u=10, p=0.5)
        } else if (pcode == 5) {
            prior = new.prior(dofs, u=10, p=0.6)
        } else if (pcode == 6) {
            prior = new.prior(dofs, u=10, p=0.7)
        } else if (pcode == 7) {
            prior = new.prior(dofs, u=10, p=0.8)
        } else if (pcode == 8) {
            prior = new.prior(dofs, lambda = 1/3.)
        } else if (pcode == 9) {
            prior = new.prior(dofs, lambda = 1/8.)
        } else if (pcode == 10) {
            prior = new.prior(dofs, lambda = 1/18.)
        } else if (pcode == 11) {
            prior = new.prior(dofs, lambda = 1/98.)
        } else {
            stop(paste("This should not happen, pcode=", pcode))
        }
        post = prior$y * lik
        post = post / max(post)
        z = sum(post*c(diff(dofs), 0))
        post = post / z
        post.sum = post.sum + post / n.sim
    }

    return(list(x = dofs,  y = post.sum, pcode = pcode))
}

run.sims = function(pcodes = 1:11, mc.cores = 4, n = 100, n.sim = 100, dof.true = 100)
{
    require(parallel)
    dofs = 2 + exp(seq(-3, 10, len=1000))

    ## old code
    ## y.matrix = matrix(rt.scaled(n*n.sim, dof=dof.true), nrow = n.sim, ncol = n)
    ##llik.matrix = matrix(0, nrow = n.sim, ncol = length(dofs))
    ##for(k in 1:n.sim) {
    ##    y = y.matrix[k, ]
    ##    llik.matrix[k, ] = sapply(dofs, function(dof, y) sum(dt.scaled(y, dof=dof, log=TRUE)), y=y)
    ##}
    
    llik.matrix = matrix(
            unlist(
                mclapply(1:n.sim,
                         function(k) {
                             sapply(dofs,
                                    function(dof, y) sum(dt.scaled(y, dof=dof, log=TRUE)),
                                    y=rt.scaled(n, dof=dof.true))
                         }, mc.cores = mc.cores)), 
            nrow = n.sim, ncol = length(dofs), byrow = TRUE)
    
    return (mclapply(pcodes, run.sim,
                     mc.cores = mc.cores,
                     dofs = dofs, n=n, n.sim=n.sim, dof.true=dof.true,
                     llik.matrix = llik.matrix))
}

if (FALSE) {
    res.dof5.n10      = run.sims(n=    10, n.sim = 1000, dof.true =   5, mc.cores = 12)
    res.dof5.n100     = run.sims(n=   100, n.sim = 1000, dof.true =   5, mc.cores = 12)
    res.dof5.n1000    = run.sims(n=  1000, n.sim = 1000, dof.true =   5, mc.cores = 12)
    res.dof5.n10000   = run.sims(n= 10000, n.sim = 1000, dof.true =   5, mc.cores = 12)

    res.dof10.n10     = run.sims(n=    10, n.sim = 1000, dof.true =  10, mc.cores = 12)
    res.dof10.n100    = run.sims(n=   100, n.sim = 1000, dof.true =  10, mc.cores = 12)
    res.dof10.n1000   = run.sims(n=  1000, n.sim = 1000, dof.true =  10, mc.cores = 12)
    res.dof10.n10000  = run.sims(n= 10000, n.sim = 1000, dof.true =  10, mc.cores = 12)

    res.dof20.n10     = run.sims(n=    10, n.sim = 1000, dof.true =  20, mc.cores = 12)
    res.dof20.n100    = run.sims(n=   100, n.sim = 1000, dof.true =  20, mc.cores = 12)
    res.dof20.n1000   = run.sims(n=  1000, n.sim = 1000, dof.true =  20, mc.cores = 12)
    res.dof20.n10000  = run.sims(n= 10000, n.sim = 1000, dof.true =  20, mc.cores = 12)

    res.dof100.n10    = run.sims(n=    10, n.sim = 1000, dof.true = 100, mc.cores = 12)
    res.dof100.n100   = run.sims(n=   100, n.sim = 1000, dof.true = 100, mc.cores = 12)
    res.dof100.n1000  = run.sims(n=  1000, n.sim = 1000, dof.true = 100, mc.cores = 12)
    res.dof100.n10000 = run.sims(n= 10000, n.sim = 1000, dof.true = 100, mc.cores = 12)

    res.dof1000000.n10    = run.sims(n=    10, n.sim = 1000, dof.true = 1000000, mc.cores = 12)
    res.dof1000000.n100   = run.sims(n=   100, n.sim = 1000, dof.true = 1000000, mc.cores = 12)
    res.dof1000000.n1000  = run.sims(n=  1000, n.sim = 1000, dof.true = 1000000, mc.cores = 12)
    res.dof1000000.n10000 = run.sims(n= 10000, n.sim = 1000, dof.true = 1000000, mc.cores = 12)
}


if (FALSE) {
    ##
    ## this is the loss of efficiency example
    ##    
    n.sim = 100
    n = 100
    s = 0.5

    if (FALSE) {
        ## no longer needed. this new prior for 'dof' is added to R-INLA
        
        ## our new robust solution
        dofs = 2 + exp(seq(-5, 15, len=10000))
        prior = new.prior(dofs, u=10, p=0.5)
        ##prior = new.prior(dofs, lambda=1/8)
        ## theta log.dens(nu(theta))+log.jacobian(theta)
        prior.table = paste(c("table:",
                cbind(log(prior$x - 2.0), log(prior$y) + log(abs(prior$x - 2.0)))) , sep = "", collapse = " ")
    } else {
        ## set it to something as its passed on later...
        prior.table = NA
    }
    
    r.opt = r.new = NULL

    run.a.sim = function(k, n, s, prior.table)
    {
        print(k)
        x = rnorm(n)
        y = 1 + x + rnorm(n, sd = s)
        
        ## optimal solution
        r.opt = inla(y ~ 1 + x,
                data = data.frame(x, y), 
                family = "gaussian",
                control.fixed = list(prec.intercept = 0, prec = 0), 
                control.compute = list(cpo=TRUE), 
                control.family = list(
                        hyper = list(
                                prec = list(
                                        initial = log(1/s^2),
                                        fixed = TRUE))))
        r.new = inla(y ~ 1 + x, 
                data = data.frame(x, y), 
                family = "t",
                control.fixed = list(prec.intercept = 0, prec = 0), 
                control.compute = list(cpo=TRUE), 
                control.family = list(
                        hyper = list(
                                prec = list(
                                        initial = log(1/s^2),
                                        fixed = TRUE),
                                dof = list(
                                        initial = 4, 
                                        prior = "dof",
                                        param = c(20, 0.1)))), 
                control.inla = list(
                        cmin = -100, 
                        int.strategy = "grid",
                        diff.logdens = 4,
                        dz = 0.25))
        res = list(
                z.intercept = ((r.opt$summary.fixed[1, "mean"] - r.new$summary.fixed[1, "mean"]) /
                               r.opt$summary.fixed[1, "sd"]),
                z.beta = ((r.opt$summary.fixed[2, "mean"] - r.new$summary.fixed[2, "mean"]) /
                          r.opt$summary.fixed[2, "sd"]), 
                sd.ratio.intercept = r.new$summary.fixed[1, "sd"]/r.opt$summary.fixed[1, "sd"],
                sd.ratio.beta = r.new$summary.fixed[2, "sd"]/r.opt$summary.fixed[2, "sd"],
                log.cpo.diff = sum(log(r.opt$cpo$cpo)) - sum(log(r.new$cpo$cpo))
                )
        return (res)
    }

    res = mclapply(1:n.sim, run.a.sim, n=n, s=s, prior.table = prior.table, mc.cores = detectCores())
}

if (plot.figures) {
    if (write.new.figures)
        system("which psfix && psfix t*.ps")
}

if (save.and.quit) {
    q(save="yes")
}
