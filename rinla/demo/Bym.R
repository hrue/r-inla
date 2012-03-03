data(Germany)
g = system.file("demodata/germany.graph", package="INLA")
source(system.file("demodata/Bym-map.R", package="INLA"))
summary(Germany)

## just make a duplicated column
Germany = cbind(Germany,region.struct=Germany$region)

# standard BYM model
formula1 = Y ~ f(region.struct,model="besag",graph=g) +
               f(region,model="iid")

result1  =  inla(formula1,family="poisson",data=Germany,E=E)

# with linear covariate
formula2 = Y ~ f(region.struct,model="besag",graph=g) +
               f(region,model="iid") + x

result2 =  inla(formula2,family="poisson",data=Germany,E=E)

# with smooth covariate
formula3 = Y ~ f(region.struct,model="besag",graph=g) +
               f(region,model="iid") + f(x, model="rw2")

result3 =  inla(formula3,family="poisson",data=Germany,E=E)

dev.new()
par(mfrow=c(2,2))
Bym.map(result1$summary.random$region.struct$mean)
Bym.map(result2$summary.random$region.struct$mean)
Bym.map(result3$summary.random$region.struct$mean)

## An alternative to above is to use the combined model BYM, which
## contains both the "iid" and the "besag" model. This makes it
## possible to get the marginals of "idd" + "besag" which otherwise is
## only possible using linear-combinations. Here, I just repeat the
## first example. I also add some prior paramters and initial values
## as they are in the order ("idd", "besag"). Further, the internal
## representation is ("iid"+"besag", "besag"), which implies that the
## first 'n' elements of the mean vector (of length '2*n') is the sum
## and the remaining 'n' elements is the spatial term. If the
## constr=TRUE, the only the "besag" part satisfy the sum-to-zero
## constraint.
prior.iid = c(1,0.01)
prior.besag = c(1,0.001)
initial.iid = 4
initial.besag = 3

formula1.bym = Y ~ f(region, model = "bym", graph = g,
                     param = c(prior.iid, prior.besag),
                     initial = c(initial.iid, initial.besag))
result1.bym = inla(formula1.bym,family="poisson",data=Germany,E=E)


