## MODEL FOR LONGITUDINAL DATA

##read the data set
data1 <- read.table("longitudinal.txt", header=TRUE)
N <- 467


data1$U11 = data1$ID
data1$U21 = N + data1$ID

formula = y ~ b12.time + b13.timedrug + b14.sex + b15.prevoi + b16.stratum +
  f(U11 , model="iid2d",param = c(23,100,100,0),
    initial = c(-2.7,0.9,-0.22), n=2*N) +
  f(U21, b12.time, copy="U11",values = 1:(2*N))

mod1 = inla(formula, family = "gaussian", verbose = TRUE, data = data1,
  control.fixed = list(prec = 0.01, prec.intercept = 0.01),
  control.family = list(param =c(0.1,0.1), initial = log(0.345)))



## MODEL FOR SURVIVAL DATA

data2 <- read.table("survival.data.txt",header=T)

formula <- inla.surv(time,event) ~ b22.drug + b23.sex + b24.prevoi + b25.stratum
##fit two models, one with exponential likelihood and one with
##piecewise constant hazard 
mod.surv <- inla(formula, family = "exponential", data = data2)
mod.surv1 <- inla(formula, family = "coxph", data=data2)
