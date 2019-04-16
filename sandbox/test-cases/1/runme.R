library(mvtnorm)
n = 200
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
        return(numeric(0))
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
inla.pardiso.check()

inla.setOption(mkl=FALSE)
r.p = inla(y ~ -1 + f(idx, model = rg),
           data = data.frame(y, idx = 1:n),
           control.compute = list(smtp = "pardiso", openmp.strategy = "huge"))

inla.setOption(mkl=TRUE)
r.pm = inla(y ~ -1 + f(idx, model = rg),
           data = data.frame(y, idx = 1:n),
           control.compute = list(smtp = "pardiso", openmp.strategy = "huge"))

inla.setOption(mkl=FALSE)
r.t = inla(y ~ -1 + f(idx, model = rg),
         data = data.frame(y, idx = 1:n),
         control.compute = list(smtp = "taucs", openmp.strategy = "huge"))

inla.setOption(mkl=TRUE)
r.tm = inla(y ~ -1 + f(idx, model = rg),
         data = data.frame(y, idx = 1:n),
         control.compute = list(smtp = "taucs", openmp.strategy = "huge"))

inla.setOption(mkl=FALSE)
r.b = inla(y ~ -1 + f(idx, model = rg),
         data = data.frame(y, idx = 1:n),
         control.compute = list(smtp = "band", openmp.strategy = "huge"))

inla.setOption(mkl=TRUE)
r.bm = inla(y ~ -1 + f(idx, model = rg),
         data = data.frame(y, idx = 1:n),
         control.compute = list(smtp = "band", openmp.strategy = "huge"))

cat("\n\n\n")
print(rbind(PARDISO=c(round(r.p$cpu, digits = 2), "#func"=r.p$misc$nfunc), 
            PARDISO.MKL=c(round(r.pm$cpu, digits = 2), "#func"=r.pm$misc$nfunc), 
            TAUCS=c(round(r.t$cpu, digits = 2), "#func" = r.t$misc$nfunc), 
            TAUCS.MKL=c(round(r.tm$cpu, digits = 2), "#func" = r.tm$misc$nfunc), 
            BAND=c(round(r.b$cpu, digits = 2), "#func" = r.b$misc$nfunc), 
            BAND.MKL=c(round(r.bm$cpu, digits = 2), "#func" = r.bm$misc$nfunc)))
