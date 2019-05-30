## ----sett,echo=F,results='hide',message=FALSE,warning=FALSE--------------
library(knitr)
opts_chunk$set(
fig.path='figs/coreg',
message=FALSE, warning=FALSE
)
options(width=77, prompt = " ", continue = "   ")
library(INLA)
lcall <- inla.getOption('inla.call')
inla.setOption(inla.call='remote')
inla.setOption(num.threads=8)
source('R/spde-tutorial-functions.R')
set.seed(1)


## ----param---------------------------------------------------------------
alpha <- c(-5, 3, 10) ### intercept on reparametrized model
m.var <- (3:5)/10 ### random field marginal variances
kappa <- c(12, 10, 7) ### GRF scales: inverse range parameters
beta <- c(.7, .5, -.5) ### copy par.: reparam. coregionalization par.
n1 <- 99; n2 <- n1+1; n3 <- n2+1 ### number of spatial locations


## ----sloc----------------------------------------------------------------
loc1 <- cbind(runif(n1), runif(n1)) 
loc2 <- cbind(runif(n2), runif(n2)) 
loc3 <- cbind(runif(n3), runif(n3)) 


## ----rfs,results='hide'--------------------------------------------------
z1 <- rMatern(1, rbind(loc1, loc2, loc3), kappa[1], m.var[1])
z2 <- rMatern(1, rbind(loc2, loc3), kappa[2], m.var[2])
z3 <- rMatern(1, loc3, kappa[3], m.var[3])


## ----yyy-----------------------------------------------------------------
e.sd <- c(0.3, 0.2, 0.15)
y1 <- alpha[1] + z1[1:n1] + rnorm(n1, 0, e.sd[1])
y2 <- alpha[2] + beta[1] * z1[n1+1:n2] + z2[1:n2] + 
    rnorm(n2, 0, e.sd[2])
y3 <- alpha[3] + beta[2] * z1[n1+n2+1:n3] + 
    beta[3] * z2[n2+1:n3] + z3 + rnorm(n3, 0, e.sd[3])


## ----mesh----------------------------------------------------------------
mesh <- inla.mesh.2d(rbind(loc1, loc2, loc3), ###loc.domain=locd, 
                     max.edge=c(0.05, 0.2), 
                     offset=c(0.05, 0.3), cutoff=0.01)


## ----eval=F,echo=F,results='hide'----------------------------------------
## mesh$n
## plot(mesh, asp=1)
## points(loc1, pch=15, col=2)
## points(loc2, pch=16, col=3)
## points(loc3, pch=17, col=4)


## ----spde----------------------------------------------------------------
spde <- inla.spde2.pcmatern(
    mesh=mesh, alpha=2, ### mesh and smoothness parameter
    prior.range=c(0.05, 0.01), ### P(practic.range<0.05)=0.01
    prior.sigma=c(1, 0.01)) ### P(sigma>1)=0.01


## ----pcopy---------------------------------------------------------------
hc3 <- hc2 <- hc1 <- list(theta=list(prior='normal', param=c(0,10)))


## ----form----------------------------------------------------------------
form <- y ~ 0 + intercept1 + intercept2 + intercept3 + 
  f(s1, model=spde) + f(s2, model=spde) + f(s3, model=spde) + 
  f(s12, copy="s1", fixed=FALSE, hyper=hc1) + 
  f(s13, copy="s1", fixed=FALSE, hyper=hc2) + 
  f(s23, copy="s2", fixed=FALSE, hyper=hc3) 


## ----stlokA--------------------------------------------------------------
A1 <- inla.spde.make.A(mesh, loc1) 
A2 <- inla.spde.make.A(mesh, loc2) 
A3 <- inla.spde.make.A(mesh, loc3) 


## ----stack---------------------------------------------------------------
stack1 <- inla.stack(
  data=list(y=cbind(as.vector(y1), NA, NA)), A=list(A1), 
  effects=list(list(intercept1=1, s1=1:spde$n.spde))) 
stack2 <- inla.stack(
  data=list(y=cbind(NA, as.vector(y2), NA)), A=list(A2), 
  effects=list(list(intercept2=1, s2=1:spde$n.spde, 
                    s12=1:spde$n.spde)))
stack3 <- inla.stack(
  data=list(y=cbind(NA, NA, as.vector(y3))), A=list(A3), 
  effects=list(list(intercept3=1, s3=1:spde$n.spde, 
                    s13=1:spde$n.spde, 
                    s23=1:spde$n.spde)))
stack <- inla.stack(stack1, stack2, stack3) 


## ----fixnugget-----------------------------------------------------------
eprec <- list(hyper=list(theta=list(prior='pc.prec', 
                                    param=c(1, 0.01))))


## ----initheta------------------------------------------------------------
theta.ini <- c(log(1/e.sd^2), 
               c(log(sqrt(8)/kappa), log(sqrt(m.var)) 
                 )[c(1,4, 2,5, 3,6)], beta)


## ----result,results='hide'-----------------------------------------------
(result <- inla(form, rep('gaussian', 3), 
                data=inla.stack.data(stack), 
                control.family=list(eprec, eprec, eprec), 
                control.predictor=list(A=inla.stack.A(stack)),
                control.mode=list(theta=theta.ini, restart=TRUE),
                control.inla=list(int.strategy='eb')))$cpu

## ----cpu,echo=FALSE------------------------------------------------------
result$cpu

## ----mode----------------------------------------------------------------
result$logfile[grep('Number of function evaluations', result$logfile)] 
round(result$mode$theta, 2) 


## ----intercepts----------------------------------------------------------
round(cbind(true=alpha, result$summary.fix), 2) 


## ----prec----------------------------------------------------------------
round(cbind(true=c(e=e.sd^-2), result$summary.hy[1:3, ]), 4)


## ----fixed---------------------------------------------------------------
round(cbind(true=beta, result$summary.hy[10:12,]), 4)


## ----range---------------------------------------------------------------
round(cbind(true=sqrt(8)/kappa, result$summary.hy[c(4,6,8),]), 3)


## ----rfvar---------------------------------------------------------------
round(cbind(true=m.var^0.5, result$summary.hy[c(5,7,9),]), 3)


## ----zfit,eval=FALSE-----------------------------------------------------
## par(mfrow=c(2,3), mar=c(2.5,2.5,1.5,0.5), mgp=c(1.5,0.5,0))
## plot(drop(A1%*%result$summary.ran$s1$mean), z1[1:n1],
##      xlab='Posterior mean', ylab='Simulated',
##      asp=1, main='z1 in y1'); abline(0:1)
## plot(drop(A2%*%result$summary.ran$s1$mean), z1[n1+1:n2],
##      xlab='Posterior mean', ylab='Simulated',
##      asp=1, main='z1 in y2'); abline(0:1)
## plot(drop(A3%*%result$summary.ran$s1$mean), z1[n1+n2+1:n3],
##      xlab='Posterior mean', ylab='Simulated',
##      asp=1, main='z1 in y3'); abline(0:1)
## plot(drop(A2%*%result$summary.ran$s2$mean), z2[1:n2],
##      xlab='Posterior mean', ylab='Simulated',
##      asp=1, main='z2 in y2'); abline(0:1)
## plot(drop(A3%*%result$summary.ran$s2$mean), z2[n2+1:n3],
##      xlab='Posterior mean', ylab='Simulated',
##      asp=1, main='z2 in y3'); abline(0:1)
## plot(drop(A3%*%result$summary.ran$s3$mean), z3[1:n3],
##      xlab='Posterior mean', ylab='Simulated',
##      asp=1, main='z3 in y3'); abline(0:1)


## ----zfitplot,echo=FALSE,fig.width=10,heigh=4,out.width='0.9\\textwidth'----
par(mfrow=c(2,3), mar=c(2.5,2.5,1.5,0.5), mgp=c(1.5,0.5,0))
plot(drop(A1%*%result$summary.ran$s1$mean), z1[1:n1],
     xlab='Posterior mean', ylab='Simulated', 
     asp=1, main='z1 in y1'); abline(0:1)
plot(drop(A2%*%result$summary.ran$s1$mean), z1[n1+1:n2],
     xlab='Posterior mean', ylab='Simulated', 
     asp=1, main='z1 in y2'); abline(0:1)
plot(drop(A3%*%result$summary.ran$s1$mean), z1[n1+n2+1:n3],
     xlab='Posterior mean', ylab='Simulated', 
     asp=1, main='z1 in y3'); abline(0:1)
plot(drop(A2%*%result$summary.ran$s2$mean), z2[1:n2],
     xlab='Posterior mean', ylab='Simulated', 
     asp=1, main='z2 in y2'); abline(0:1)
plot(drop(A3%*%result$summary.ran$s2$mean), z2[n2+1:n3],
     xlab='Posterior mean', ylab='Simulated', 
     asp=1, main='z2 in y3'); abline(0:1)
plot(drop(A3%*%result$summary.ran$s3$mean), z3[1:n3],
     xlab='Posterior mean', ylab='Simulated', 
     asp=1, main='z3 in y3'); abline(0:1)

