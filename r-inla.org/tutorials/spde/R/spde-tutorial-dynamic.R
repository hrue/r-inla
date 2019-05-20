## ----opts,echo=F,results='hide',message=FALSE,warning=FALSE--------------
library(knitr)
opts_chunk$set(
fig.path='figs/dynamic',
message=FALSE, warning=FALSE
)
options(width=75, prompt = " ", continue = "   ")
library(INLA) 
lcall <- inla.getOption('inla.call') 
inla.setOption(inla.call='remote')
inla.setOption(num.threads=8)
source('R/spde-tutorial-functions.R')


## ----coom----------------------------------------------------------------
n <- 150;   set.seed(1);  coo <- matrix(runif(2*n), n)


## ----sample--------------------------------------------------------------
kappa <- c(10, 12);    sigma2 <- c(1/2, 1/4)
k <- 15;  rho <- c(0.7, 0.5) 
set.seed(2); beta0 <- rMatern(k, coo, kappa[1], sigma2[1]) 
set.seed(3); beta1 <- rMatern(k, coo, kappa[2], sigma2[2]) 
beta0[,1] <- beta0[,1] / (1-rho[1]^2)
beta1[,1] <- beta1[,1] / (1-rho[2]^2)
for (j in 2:k) {
    beta0[, j] <- beta0[,j-1]*rho[1] + beta0[,j] * (1-rho[1]^2)
    beta1[, j] <- beta1[,j-1]*rho[2] + beta1[,j] * (1-rho[2]^2)
}


## ----response------------------------------------------------------------
set.seed(4); hh <- runif(n*k) ### simlate the covariate values
mu.beta <- c(-5, 1);   taue <- 20 
set.seed(5); error <- rnorm(n*k, 0, sqrt(1/taue)) ### error in the observation
length(y <- (mu.beta[1] + beta0) + (mu.beta[2]+beta1)*hh + ### dynamic regression part
           error)


## ----mesh----------------------------------------------------------------
(mesh <- inla.mesh.2d(coo, max.edge=c(0.25), ### coarse mesh
                      offset=c(0.15), cutoff=0.05))$n

## ----vmesh,eval=FALSE,echo=FALSE,results='hide'--------------------------
## mesh$n
## plot(mesh, asp=1)
## points(coo, pch=4, col=2)


## ----spde----------------------------------------------------------------
spde <- inla.spde2.pcmatern(
    mesh=mesh, alpha=2, ### mesh and smoothness parameter
    prior.range=c(0.05, 0.01), ### P(practic.range<0.05)=0.01
    prior.sigma=c(1, 0.01)) ### P(sigma>1)=0.01


## ----idx-----------------------------------------------------------------
i0 <- inla.spde.make.index('i0', spde$n.spde, n.group=k)
i1 <- inla.spde.make.index('i1', spde$n.spde, n.group=k)


## ----spdebuild-----------------------------------------------------------
A0 <- inla.spde.make.A(mesh, cbind(rep(coo[,1], k), rep(coo[,2], k)),
                       group=rep(1:k, each=n))
A1 <- inla.spde.make.A(mesh, cbind(rep(coo[,1], k), rep(coo[,2], k)),
                       group=rep(1:k, each=n), weights=hh)


## ----stky----------------------------------------------------------------
stk.y <- inla.stack(data=list(y=as.vector(y)), tag='y', 
                    A=list(A0, A1, 1), 
                    effects=list(i0, i1, 
                        data.frame(mu1=1, h=hh)))


## ----formula-------------------------------------------------------------
form <- y ~ 0 + mu1 + h + ### to fit mu_beta
    f(i0, model=spde, group=i0.group, control.group=list(model='ar1')) + 
        f(i1, model=spde, group=i1.group, control.group=list(model='ar1'))


## ----theta---------------------------------------------------------------
(theta.ini <- c(log(taue), ## likelihood log precision
                log(sqrt(8)/kappa[1]), ## log range 1
                log(sqrt(sigma2[1])), ## log stdev 1
                log((1+rho[1])/(1-rho[1])), ## inv.logit rho 1
                log(sqrt(8)/kappa[2]), ## log range 1
                log(sqrt(sigma2[2])), ## log stdev 1
                log((1+rho[2])/(1-rho[2]))))## inv.logit rho 2


## ----fittingdyn3---------------------------------------------------------
(res <- inla(form, family='gaussian', data=inla.stack.data(stk.y), 
            control.predictor=list(A=inla.stack.A(stk.y)),
            control.inla=list(int.strategy='eb'), ### no integration wr theta
            control.mode=list(theta=theta.ini, ### initial theta value
                              restart=TRUE)))$cpu 


## ----summarymux----------------------------------------------------------
round(cbind(true=mu.beta, res$summary.fix), 4)


## ----likprec-------------------------------------------------------------
round(c(true=taue, unlist(res$summary.hy[1,])), 3)


## ----hd3pmds, echo=TRUE, eval=FALSE--------------------------------------
## par(mfrow=c(2, 3), mar=c(2.5,2.5,0.3,0.3), mgp=c(1.5,0.5,0))
## for (j in 2:7) {
##     plot(res$marginals.hy[[j]], type='l',
##          xlab=names(res$marginals.hyperpar)[j], ylab='Density')
##     abline(v=c(sqrt(8)/kappa[1], sigma2[1]^0.5, rho[1],
##                sqrt(8)/kappa[2], sigma2[2]^0.5, rho[2])[j-1])
## }


## ----hd3pmdsf, eval=TRUE, echo=FALSE, fig.width=7.5, fig.height=4, out.width='0.97\\textwidth'----
par(mfrow=c(2, 3), mar=c(2.5,2.5,0.3,0.3), mgp=c(1.5,0.5,0)) 
for (j in 2:7) {
    plot(res$marginals.hy[[j]], type='l', 
         xlab=names(res$marginals.hyperpar)[j], ylab='Density')
    abline(v=c(sqrt(8)/kappa[1], sigma2[1]^0.5, rho[1], 
               sqrt(8)/kappa[2], sigma2[2]^0.5, rho[2])[j-1])
}


## ----betas---------------------------------------------------------------
c(beta0=cor(as.vector(beta0), drop(A0%*%res$summary.ran$i0$mean)), 
  beta1=cor(as.vector(beta1), 
      drop(A0%*%res$summary.ran$i1$mean))) ## using A0 to account only for the coeff.

