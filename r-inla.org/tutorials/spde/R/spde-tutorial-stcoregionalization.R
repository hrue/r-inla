## ----sett,echo=F,results='hide'------------------------------------------
library(knitr)
opts_chunk$set(
fig.path='figs/stcoreg'
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
rho <- c(.7, .8, .9) ### temporal correlations
n <- 300;  k <- 9  ### number of spatial locations and time points

## ----sloc----------------------------------------------------------------
loc <- cbind(runif(n), runif(n)) 

## ----rfs,results='hide'--------------------------------------------------
x1 <- rMatern(k, loc, kappa[1], m.var[1])
x2 <- rMatern(k, loc, kappa[2], m.var[2])
x3 <- rMatern(k, loc, kappa[3], m.var[3])

## ----st------------------------------------------------------------------
z1 <- x1; z2 <- x2; z3 <- x3
for (j in 2:k) {
    z1[, j] <- rho[1] * z1[,j-1] + sqrt(1-rho[1]^2) * x1[,j]
    z2[, j] <- rho[2] * z2[,j-1] + sqrt(1-rho[2]^2) * x2[,j]
    z3[, j] <- rho[3] * z3[,j-1] + sqrt(1-rho[3]^2) * x3[,j]
}   

## ----yyy-----------------------------------------------------------------
e.sd <- c(0.3, 0.2, 0.15)
y1 <- alpha[1] + z1 + rnorm(n, 0, e.sd[1])
y2 <- alpha[2] + beta[1] * z1 + z2 + rnorm(n, 0, e.sd[2])
y3 <- alpha[3] + beta[2] * z1 + beta[3] * z2 + z3 + 
    rnorm(n, 0, e.sd[3])

## ----mesh----------------------------------------------------------------
mesh <- inla.mesh.2d(loc, max.edge=0.2, 
                     offset=0.1, cutoff=0.1)

## ----eval=F,echo=F,results='hide'----------------------------------------
## mesh$n
## plot(mesh, asp=1)
## points(loc, pch=19)

## ----spde----------------------------------------------------------------
spde <- inla.spde2.pcmatern(
    mesh=mesh, alpha=2, ### mesh and smoothness parameter
    prior.range=c(0.05, 0.01), ### P(practic.range<0.05)=0.01
    prior.sigma=c(1, 0.01)) ### P(sigma>1)=0.01

## ----idx-----------------------------------------------------------------
s1 = s2 = s3 = s12 = s13 = s23 = rep(1:spde$n.spde, times=k)
g1 = g2 = g3 = g12 = g13 = g23 = rep(1:k, each=spde$n.spde)

## ----pbeta---------------------------------------------------------------
rho1p <- list(theta=list(prior='pccor1', param=c(0, 0.9))) 
ctr.g <- list(model='ar1', hyper=rho1p)

## ----pcopy---------------------------------------------------------------
hc3 <- hc2 <- hc1 <- list(theta=list(prior='normal', param=c(0,10)))

## ----form----------------------------------------------------------------
form <- y ~ 0 + intercept1 + intercept2 + intercept3 + 
  f(s1, model=spde, ngroup=k, group=g1, control.group=ctr.g) + 
  f(s2, model=spde, ngroup=k, group=g2, control.group=ctr.g) + 
  f(s3, model=spde, ngroup=k, group=g3, control.group=ctr.g) + 
  f(s12, copy="s1", group=g12, fixed=FALSE, hyper=hc1) + 
  f(s13, copy="s1", group=g13, fixed=FALSE, hyper=hc2) + 
  f(s23, copy="s2", group=g23, fixed=FALSE, hyper=hc3) 

## ----stlokA--------------------------------------------------------------
stloc <- kronecker(matrix(1,k,1), loc) ### rep. coordinates each time
A <- inla.spde.make.A(mesh, stloc, n.group=k, group=rep(1:k, each=n))

## ----stack---------------------------------------------------------------
stack1 <- inla.stack(
  data=list(y=cbind(as.vector(y1), NA, NA)), A=list(A), 
  effects=list(list(intercept1=1, s1=s1, g1=g1))) 
stack2 <- inla.stack(
  data=list(y=cbind(NA, as.vector(y2), NA)), A=list(A), 
  effects=list(list(intercept2=1, s2=s2, g2=g2, 
                    s12=s12, g12=g12))) 
stack3 <- inla.stack(
  data=list(y=cbind(NA, NA, as.vector(y3))), A=list(A), 
  effects=list(list(intercept3=1, s3=s3, g3=g3, 
                    s13=s13, g13=g13, s23=s23, g23=g23))) 
stack <- inla.stack(stack1, stack2, stack3) 

## ----fixnugget-----------------------------------------------------------
eprec <- list(hyper=list(theta=list(prior='pc.prec', 
                                    param=c(1, 0.01))))

## ----initheta------------------------------------------------------------
theta.ini <- c(log(1/e.sd^2), 
               c(log(sqrt(8)/kappa), log(sqrt(m.var)), 
                 qlogis(rho))[c(1,4,7, 2,5,8, 3,6,9)], beta)

## ----result,results='hide'-----------------------------------------------
(result <- inla(form, rep('gaussian', 3), data=inla.stack.data(stack), 
                control.family=list(eprec, eprec, eprec), 
                control.mode=list(theta=theta.ini, restart=TRUE),
                control.inla=list(int.strategy='eb'), 
                control.predictor=list(A=inla.stack.A(stack))))$cpu

## ----cpu,echo=FALSE------------------------------------------------------
result$cpu

## ----mode----------------------------------------------------------------
result$logfile[grep('Number of function evaluations', result$logfile)] 
round(result$mode$theta, 2) 

## ----intercepts----------------------------------------------------------
round(cbind(true=alpha, result$summary.fix), 2) 

## ----prec----------------------------------------------------------------
round(cbind(true=c(e=e.sd^-2), result$summary.hy[1:3, ]), 4)

## ----rho-----------------------------------------------------------------
round(cbind(true=rho, result$summary.hy[c(6,9,12),]), 4) 

## ----fixed---------------------------------------------------------------
round(cbind(true=beta, result$summary.hy[13:15,]), 4)

## ----range---------------------------------------------------------------
round(cbind(true=sqrt(8)/kappa, result$summary.hy[c(4, 7, 10),]), 3)

## ----rfvar---------------------------------------------------------------
round(cbind(true=m.var^0.5, result$summary.hy[c(5, 8, 11),]), 3)

## ----stzfit,eval=FALSE---------------------------------------------------
## par(mfrow=c(1,3), mar=c(2,2,0.5,0.5), mgp=c(1.5,0.5,0))
## plot(drop(A%*%result$summary.ran$s1$mean), as.vector(z1),
##      xlab='', ylab='', asp=1); abline(0:1)
## plot(drop(A%*%result$summary.ran$s2$mean), as.vector(z2),
##      xlab='', ylab='', asp=1); abline(0:1)
## plot(drop(A%*%result$summary.ran$s3$mean), as.vector(z3),
##      xlab='', ylab='', asp=1); abline(0:1)

## ----stzfitplot,echo=FALSE,fig.width=15,heigh=3,out.width='0.97\\textwidth'----
par(mfrow=c(1,3), mar=c(2,2,0.5,0.5), mgp=c(1.5,0.5,0))
plot(drop(A%*%result$summary.ran$s1$mean), as.vector(z1),
     xlab='', ylab='', asp=1); abline(0:1)
plot(drop(A%*%result$summary.ran$s2$mean), as.vector(z2),
     xlab='', ylab='', asp=1); abline(0:1)
plot(drop(A%*%result$summary.ran$s3$mean), as.vector(z3),
     xlab='', ylab='', asp=1); abline(0:1)

