remove(list=ls())
require(INLA)
##inla.my.update()
##source("~/hg/inla/rinla/R/interpret.formula.R")
##source("~/hg/inla/rinla/R/f.R")
##source("~/hg/inla/rinla/R/hyper.R")
##source("~/hg/inla/rinla/R/models.R")
##source("~/hg/inla/rinla/R/mesh.R")
##source("~/hg/inla/rinla/R/fmesher.R")
##source("~/hg/inla/rinla/R/spde.R")
##source("~/hg/inla/rinla/R/spde2.R")

mesh =
    inla.mesh.create(matrix(runif(8),4,2),
                     refine=list(min.angle=30, max.edge=0.1))
spde0 = inla.spde.create(mesh)
field = inla.spde.query(spde0, sample=list(tau=1, kappa2=8/0.3^2))$sample

spde1 = inla.spde2.matern(mesh)
spde2 = inla.spde2.matern(mesh)

data =
    list(
         y = field+rnorm(length(field))*diff(range(field))*0.1,
         field1=1:mesh$n,
         field2=mesh$n:1,
         spde0=spde0,
         spde1=spde1,
         spde2=spde2
         )

formula0 = as.formula("y ~ 1 + f(field1, model=spde0)")
formula = as.formula("y ~ -1 + f(field1, model=spde1)+ f(field2, model=spde2)")

##ow=options("warn")
##options(warn=2)
##result0 =
##    inla(formula0, data=data, family="gaussian",
##         verbose=TRUE)
result =
    inla(formula, data=data, family="gaussian",
         verbose=FALSE)
print(result)

##options(ow)

result.field1 = inla.spde2.result(result,"field1",spde1)
result.field2 = inla.spde2.result(result,"field2",spde2)
