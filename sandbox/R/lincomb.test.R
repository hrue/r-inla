remove(list=ls())
require(INLA)
source("~/hg/inla/rinla/R/spde2.R")

n = 10
n.repl = 2
loc = matrix(runif(n*2), n, 2)

mesh = inla.mesh.create.helper(loc, max.edge=c(0.02,1))
spde = inla.spde2.matern(mesh)

Q = inla.spde.precision(spde, theta=c(0,0))

## Simulated data, n.repl replicates:
field = c()
for (k in 1:n.repl)
    field = c(field, inla.spde.sample(Q))

## Observations:
obs.idx0 = rep(mesh$idx$loc, n.repl)
obs.idx = obs.idx0 + rep((0:(n.repl-1))*mesh$n, each=n)
y = field[obs.idx] + rnorm(n*n.repl)*0.01
repl = rep(1:n.repl, each=n)

formula = y ~ -1 + f(field, model=spde, replicate=repl)

A.matrices =
    list(c(),
         c(),
         kronecker(Diagonal(n.repl), inla.mesh.project(mesh, loc)$A),
         rBind(kronecker(Diagonal(n.repl), inla.mesh.project(mesh, loc)$A),
               Diagonal(mesh$n*n.repl)),
         c()
         )

lincombs =
    list(c(),
         c(),
         c(),
         c(),
         inla.make.lincombs(field = (as(sparseMatrix(i=1:(mesh$n*n.repl),
                                                     j=1:(mesh$n*n.repl),
                                                     x=rep(1.0, mesh$n*n.repl)),
                                        "dgTMatrix")))
         )

data =
    list(
         list(spde=spde,
              y = y,
              field = obs.idx0,
              repl = repl)
         ,
         list(spde=spde,
              y = c(y, rep(NA, mesh$n*n.repl)),
              field = c(obs.idx0, rep(1:mesh$n, n.repl)),
              repl = c(repl, rep(c(1:n.repl), each=mesh$n)))
         ,
         list(spde=spde,
              y = y,
              field = rep(1:mesh$n, n.repl),
              repl = rep(c(1:n.repl), each=mesh$n))
         ,
         list(spde=spde,
              y = c(y, rep(NA, mesh$n*n.repl)),
              field = rep(1:mesh$n, n.repl),
              repl = rep(c(1:n.repl), each=mesh$n))
         ,
         list(spde=spde,
              y = y,
              field = obs.idx0,
              repl = repl)
         )

result = list()
for (k in 1:length(data)) {
    print(k)
    if (is.null(A.matrices[[k]])) {
        if (is.null(lincombs[[k]])) {
            result =
                c(result,
                  list(inla(formula, data=data[[k]], family="gaussian",
                            control.predictor=list(compute=TRUE),
                            control.compute=list(cpo=TRUE),
                            verbose=FALSE)))
        } else {
            result =
                c(result,
                  list(inla(formula, data=data[[k]], family="gaussian",
                            lincomb = lincombs[[k]],
                            control.lincomb = list(verbose=FALSE),
                            control.predictor=list(compute=TRUE),
                            control.compute=list(cpo=TRUE),
                            verbose=FALSE)))
            }
    } else {
        result =
            c(result,
              list(inla(formula, data=data[[k]], family="gaussian",
                        control.lincomb = list(verbose=FALSE),
                        control.predictor=list(A=A.matrices[[k]], compute=TRUE),
                        control.compute=list(cpo=TRUE),
                        verbose=FALSE)))
    }
    print(result[[k]]$cpu.used)
}
print("------------")
cpu.used = result[[1]]$cpu.used
if (length(data)>1)
    for (k in 2:length(data)) {
        cpu.used = rbind(cpu.used, result[[k]]$cpu.used)
    }
cpu.used.rel =
    cpu.used / matrix(cpu.used[1,], length(data), 4, byrow=TRUE)
print(cpu.used)
print(cpu.used.rel)
