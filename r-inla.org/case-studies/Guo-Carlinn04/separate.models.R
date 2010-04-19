## MODEL FOR LONGITUDINAL DATA

##read the data set
data1 <- read.table("longitudinal.txt", header=TRUE)
N <- 467

## the model can be defined both using the 2diidwishard(part1 and part2) model or the
## 2diidwishard using "copy"

##definition using 2diidwishartp1 and 2diidwishartp2
## we need another colum for the random effect
data1$ID2 <- data1$ID

formula = y ~ b12.time + b13.timedrug + b14.sex + b15.prevoi + b16.stratum +
  f(ID , model="2diidwishartp1",param = c(23,100,100,0),initial = c(-2.7,0.9,-0.22)) +
  f(ID2, b12.time, model="2diidwishartp2")

mod = inla(formula, family = "gaussian", verbose = TRUE, data = data1,
  control.fixed = list(prec = 0.01, prec.intercept = 0.01),
  control.data = list(param =c(0.1,0.1), initial = log(0.345)))

## alternative definition using 2diidwishart and copy
data1$U11 = data1$ID
data1$U21 = N + data1$ID

formula = y ~ b12.time + b13.timedrug + b14.sex + b15.prevoi + b16.stratum +
  f(U11 , model="2diidwishart",param = c(23,100,100,0),
    initial = c(-2.7,0.9,-0.22),
    n=N, values = 1:(2*N)) +
  f(U21, b12.time, copy="U11",values = 1:(2*N))

mod1 = inla(formula, family = "gaussian", verbose = TRUE, data = data1,
  control.fixed = list(prec = 0.01, prec.intercept = 0.01),
  control.data = list(param =c(0.1,0.1), initial = log(0.345)))



## MODEL FOR SURVIVAL DATA

data2 <- read.table("survival.data.txt",header=T)

formula <- inla.surv(time,event) ~ b22.drug + b23.sex + b24.prevoi + b25.stratum
##fit two models, one with exponential likelihood and one with
##piecewise constant hazard 
mod.surv <- inla(formula, family = "exponential", data = data2)
mod.surv1 <- inla(formula, family = "piecewiseconstant", data=data2)


