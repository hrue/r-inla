library(INLA)
inla.pardiso.check()

library(mvtnorm)
n = 1000
A = matrix(rnorm(n^2), n, n)
Q = A %*% t(A)
y = c(rmvnorm(1, sigma = solve(Q)))

rg.test = function (cmd = c("graph", "Q", "mu", "initial",
                            "log.norm.const", "log.prior",
                            "quit"),
                    theta = NULL) 
{
    envir = parent.env(environment())
    interpret.theta = function() {
        return(list(prec = exp(theta[1L])))
    }
    graph = function() {
        return (Q())
    }
    Q = function() {
        prec = interpret.theta()$prec
        return(prec * Q.matrix)
    }
    mu = function() {
        return(numeric(0))
    }
    log.norm.const = function() {
        ## ignore the |Q.matrix| like model=generic0
        prec = interpret.theta()$prec
        return (dim(Q.matrix)[1]/2.0 * (log(prec) - log(2*pi)))
    }
    log.prior = function() {
        prec = interpret.theta()$prec
        val = dgamma(prec, shape = 1, rate = 1, log = TRUE) + 
            theta[1L]
        return(val)
    }
    initial = function() {
        return (1)
    }
    quit = function() {
        return(invisible())
    }
    if (is.null(theta)) {
        theta = initial()
    }

    val = do.call(match.arg(cmd), args = list())
    return(val)
}

rg = inla.rgeneric.define(rg.test,  Q.matrix = Q, debug=FALSE)
models = c('y ~ -1 + f(idx, model = rg)',
           'y ~ -1 + f(idx, model="generic0", Cmatrix=Q, hyper = list(prec = list(prior = "loggamma", param = c(1, 1))))')
           
##
for(model in models) {
    formula = as.formula(model)
    
    inla.setOption(mkl=FALSE)
    r.pp = inla(formula, 
                data = data.frame(y, idx = 1:n),
                control.compute = list(smtp = "pardiso", openmp.strategy = "pardiso.parallel"))

    inla.setOption(mkl=TRUE)
    r.ppm = inla(formula, 
                 data = data.frame(y, idx = 1:n),
                 control.compute = list(smtp = "pardiso", openmp.strategy = "pardiso.parallel"))

    inla.setOption(mkl=FALSE)
    r.p = inla(formula, 
               data = data.frame(y, idx = 1:n),
               control.compute = list(smtp = "pardiso", openmp.strategy = "pardiso.serial"))

    inla.setOption(mkl=TRUE)
    r.pm = inla(formula, 
                data = data.frame(y, idx = 1:n),
                control.inla = list(int.strategy = "eb"), 
                control.compute = list(smtp = "pardiso", openmp.strategy = "pardiso.serial"))

    inla.setOption(mkl=FALSE)
    r.t = inla(formula, 
               data = data.frame(y, idx = 1:n),
               control.compute = list(smtp = "taucs", openmp.strategy = "huge"))

    inla.setOption(mkl=TRUE)
    r.tm = inla(formula, 
                data = data.frame(y, idx = 1:n),
                control.compute = list(smtp = "taucs", openmp.strategy = "huge"))

    inla.setOption(mkl=FALSE)
    r.b = inla(formula, 
               data = data.frame(y, idx = 1:n),
               control.compute = list(smtp = "band", openmp.strategy = "huge"))

    inla.setOption(mkl=TRUE)
    r.bm = inla(formula, 
                data = data.frame(y, idx = 1:n),
                control.compute = list(smtp = "band", openmp.strategy = "huge"))

    cat("\n\n\n")
    print(paste0("Using model = [", model, "]"))
    print(rbind(PARDISO.PARALLEL=c(round(r.pp$cpu, digits = 2), "#func"=r.pp$misc$nfunc), 
                PARDISO.PARALLEL.MKL=c(round(r.ppm$cpu, digits = 2), "#func"=r.ppm$misc$nfunc), 
                PARDISO.SERIAL=c(round(r.p$cpu, digits = 2), "#func"=r.p$misc$nfunc), 
                PARDISO.SERIAL.MKL=c(round(r.pm$cpu, digits = 2), "#func"=r.pm$misc$nfunc), 
                TAUCS=c(round(r.t$cpu, digits = 2), "#func" = r.t$misc$nfunc), 
                TAUCS.MKL=c(round(r.tm$cpu, digits = 2), "#func" = r.tm$misc$nfunc), 
                BAND=c(round(r.b$cpu, digits = 2), "#func" = r.b$misc$nfunc), 
                BAND.MKL=c(round(r.bm$cpu, digits = 2), "#func" = r.bm$misc$nfunc)))
}
