## JOINT MODELS
## Here is the code to run all joint models proposed in the papaer, except for model X and model XII
##read the data set
data1 <- read.table("longitudinal.txt", header=TRUE)
data2 <- read.table("survival.data.txt",header=T)

N <- 467

##prepare the data set
ng = dim(data1)[1]
ns  = N

## prepare the response variable
y.long <- c(data1$y, rep(NA, ns))
y.surv <- inla.surv(time = c(rep(NA, ng), data2$time), event = c(rep(NA, ng),data2$event))
Yjoint <- list(y.long, y.surv)

##prepare the fixed covariate
##NB for fixed effects covariate use 0's where the covariate is not defined!!!
## for randon effects you have to use NA!!!

linear.covariate <- data.frame(mu = as.factor(c(rep(1, ng), rep(2, ns))),
                               b12.time = c(data1$b12.time, rep(0, ns)),
                               b13.timedrug = c(data1$b13.timedrug, rep(0, ns)),
                               b14.sex = c(data1$b14.sex, rep(0, ns)),
                               b15.prevoi =  c(data1$b15.prevoi, rep(0, ns)), 
                               b16.stratum =  c(data1$b16.stratum, rep(0, ns)),
                               b22.drug = c(rep(0, ng), data2$b22.drug),
                               b23.sex = c(rep(0, ng), data2$b23.sex),
                               b24.prevoi = c(rep(0, ng), data2$b24.prevoi),
                               b25.stratum = c(rep(0, ng), data2$b25.stratum))


random.covariate <- list(U11 = c(rep(1:N, each=5),rep(NA, ns)),
                         U21 = c(rep(N+(1:N), each=5),rep(NA, ns)),
                         U12 = c(rep(NA,ng), 1:N),
                         U22 = c(rep(NA,ng), N+(1:N)),
                         U3 = c(rep(NA,ng),1:N))

joint.data <- c(linear.covariate,random.covariate)
joint.data$Y <- Yjoint


##MODEL I
formula1 = Y ~ mu + b12.time + b13.timedrug + b14.sex + b15.prevoi +
  b16.stratum + b22.drug + b23.sex + b24.prevoi + b25.stratum - 1

mod1 = inla(formula1, family = c("gaussian","exponential"),
  data = joint.data, verbose=TRUE, control.compute=list(dic=TRUE))

##MODEL II
formula2 = Y ~ mu + b12.time + b13.timedrug + b14.sex + b15.prevoi +
  b16.stratum + b22.drug + b23.sex + b24.prevoi + b25.stratum - 1 +
  f(U3, model="iid")

mod2 = inla(formula2, family = c("gaussian","exponential"),
  data = joint.data, verbose=TRUE, control.compute=list(dic=TRUE))

##MODEL III
formula3 = Y ~ mu + b12.time + b13.timedrug + b14.sex + b15.prevoi +
  b16.stratum + b22.drug + b23.sex + b24.prevoi + b25.stratum - 1 +
  f(U11 , model="iid")

mod3 = inla(formula3, family = c("gaussian","exponential"),
  data = joint.data, verbose=TRUE, control.compute=list(dic=TRUE))

##MODEL IV
formula4 = Y ~ mu + b12.time + b13.timedrug + b14.sex + b15.prevoi +
  b16.stratum + b22.drug + b23.sex + b24.prevoi + b25.stratum - 1 +
  f(U11 , model="iid") + f(U3,model="iid")

mod4 = inla(formula4, family = c("gaussian","exponential"),
  data = joint.data, verbose=TRUE, control.compute=list(dic=TRUE))

##MODEL V
formula5 = Y ~ mu + b12.time + b13.timedrug + b14.sex + b15.prevoi +
  b16.stratum + b22.drug + b23.sex + b24.prevoi + b25.stratum - 1 +
  f(U11 , model="iid") + f(U12,copy="U11",fixed=FALSE)

mod5 = inla(formula5, family = c("gaussian","exponential"),
  data = joint.data, verbose=TRUE, control.compute=list(dic=TRUE))

##MODEL VI
formula6 = Y ~ mu + b12.time + b13.timedrug + b14.sex + b15.prevoi +
  b16.stratum + b22.drug + b23.sex + b24.prevoi + b25.stratum - 1 +
  f(U11 , model="iid") + f(U12,copy="U11",fixed=FALSE)+ f(U3,model="iid")

mod6 = inla(formula6, family = c("gaussian","exponential"),
  data = joint.data, verbose=TRUE, control.compute=list(dic=TRUE))

##MODEL VII
formula7 = Y ~ mu + b12.time + b13.timedrug + b14.sex + b15.prevoi +
  b16.stratum + b22.drug + b23.sex + b24.prevoi + b25.stratum - 1 +
  f(U11 , model="iid2d", param = c(23,100,100,0), initial = c(-2.7,0.9,-0.22), n=2*N) +
  f(U21, b12.time, copy="U11")

mod7 = inla(formula7, family = c("gaussian","exponential"),
  data = joint.data, verbose=TRUE, control.compute=list(dic=TRUE))

##MODEL VIII
formula8 = Y ~ mu + b12.time + b13.timedrug + b14.sex + b15.prevoi +
  b16.stratum + b22.drug + b23.sex + b24.prevoi + b25.stratum - 1 +
  f(U11 , model="iid2d", param = c(23,100,100,0),
    initial = c(-2.7,0.9,-0.22), n=2*N) +
  f(U21, b12.time, copy="U11") +
  f(U12, copy="U11", fixed = FALSE, param=c(0,0.01), initial = -0.2)

mod8 = inla(formula8, family = c("gaussian","exponential"),
  data = joint.data, verbose=TRUE, control.compute=list(dic=TRUE))


##MODEL IX
formula9 = Y ~ mu + b12.time + b13.timedrug + b14.sex + b15.prevoi +
  b16.stratum + b22.drug + b23.sex + b24.prevoi + b25.stratum - 1 +
  f(U11 , model="iid2d", param = c(23,100,100,0),
    initial = c(-2.7,0.9,-0.22), n=2*N) +
  f(U21, b12.time, copy="U11") +
  f(U22, copy="U11", fixed = FALSE, param=c(0,0.01), initial = -0.2)

mod9 = inla(formula9, family = c("gaussian","exponential"),
  data = joint.data, verbose=TRUE, control.compute=list(dic=TRUE))


##MODEL XI
formula11 = Y ~ mu + b12.time + b13.timedrug + b14.sex + b15.prevoi +
  b16.stratum + b22.drug + b23.sex + b24.prevoi + b25.stratum - 1 +
  f(U11 , model="iid2d", param = c(23,100,100,0),
    initial = c(-2.7,0.9,-0.22), n=2*N) +
  f(U21, b12.time, copy="U11") +
  f(U12, copy="U11", fixed = FALSE, param=c(0,0.01), initial = -0.2) +
  f(U22, copy="U11", fixed = FALSE, param=c(0,0.01), initial = -1.6)

mod11 = inla(formula11, family = c("gaussian","exponential"),
  data = joint.data, verbose=TRUE, control.compute=list(dic=TRUE))
