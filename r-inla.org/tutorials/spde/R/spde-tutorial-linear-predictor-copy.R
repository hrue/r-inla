## ----sett,echo=F,results='hide',message=FALSE,warning=FALSE--------------
library(knitr)
opts_chunk$set(
fig.path='figs/copy',
message=FALSE, warning=FALSE
)
options(width=77, prompt = " ", continue = "   ")
library(INLA)
lcall <- inla.getOption('inla.call')
inla.setOption(inla.call='remote')
inla.setOption(num.threads=8)
source('R/spde-tutorial-functions.R')
set.seed(1)


## ----truebeta------------------------------------------------------------
beta0 = 5
beta1 = 1
beta2 = 0.5
beta3 = 2


## ----s123true------------------------------------------------------------
s123 <- c(0.1, 0.05, 0.15)


## ----kappas2-------------------------------------------------------------
kappab <- 10
sigma2b <- 1


## ----loc-----------------------------------------------------------------
n <- 50
loc <- cbind(runif(n), runif(n))


## ------------------------------------------------------------------------
b <- rMatern(n=1, coords=loc, kappa=10, variance=1)


## ----vizloc,eval=FALSE---------------------------------------------------
## par(mar=c(0,0,0,0))
## plot(loc, asp=1, cex=0.3+2*(b-min(b))/diff(range(b)),
##      pch=19, axes=FALSE); box()


## ----vvizloc,echo=F,results='hide',fig.width=5.5,fig.height=5.5----------
par(mar=c(0,0,0,0))
plot(loc, asp=1, cex=0.3+2*(b-min(b))/diff(range(b)), 
     pch=19, axes=FALSE); box()


## ----covariate-----------------------------------------------------------
x <- runif(n, -1, 1)*sqrt(3)


## ----predictor-----------------------------------------------------------
eta1 <- beta0 + beta1*x + b
eta2 <- beta2*(beta0 + beta1*x)
eta3 <- beta3*eta1


## ----obs-----------------------------------------------------------------
y1 <- rnorm(n, eta1, s123[1])
y2 <- rnorm(n, eta2, s123[2])
y3 <- rnorm(n, eta3, s123[3])


## ----mesh----------------------------------------------------------------
mesh <- inla.mesh.2d(loc.domain=cbind(c(0,1,1,0), c(0,0,1,1)),
                     max.edge=c(0.1, 0.3),
                     offset=c(0.05, 0.35), cutoff=0.05)


## ----As------------------------------------------------------------------
As <- inla.spde.make.A(mesh, loc)


## ----spde----------------------------------------------------------------
spde <- inla.spde2.matern(mesh, alpha=2)


## ----stack1--------------------------------------------------------------
stack1 <- inla.stack(tag='y1',
                     data=list(y=y1), 
                     effects=list(
                         data.frame(beta0=1, beta1=x),
                         s=1:spde$n.spde,
                         e1=1:n),
                     A=list(1, As, 1))


## ----stack01-------------------------------------------------------------
stack01 <- inla.stack(tag='eta1',
                      data=list(y=rep(0,n), offset=-y1),
                      effects=list(
                          s=1:spde$n.spde,
                          list(e1=1:n, 
                               eta1=1:n)),
                      A=list(As, 1))


## ----stack02-------------------------------------------------------------
stack02 <- inla.stack(tag='eta2',
                      data=list(y=rep(0,n), offset=-y1),
                      effects=list(list(e1=1:n, eta2=1:n)),
                      A=list(1))


## ----stack2--------------------------------------------------------------
stack2 <- inla.stack(tag='y2',
                     data=list(y=y2),
                     effects=list(list(eta1c=1:n, e2=1:n)), 
                     A=list(1))


## ----stack3--------------------------------------------------------------
stack3 <- inla.stack(tag='y3',
                     data=list(y=y3), 
                     effects=list(list(eta2c=1:n, e3=1:n)),
                     A=list(1))


## ----stacks--------------------------------------------------------------
stack <- inla.stack(stack1, stack01, stack02, stack2, stack3)


## ----pcprec--------------------------------------------------------------
pcprec <- list(theta=list(prior='pcprec', param=c(0.5, 0.1)))


## ----formula123----------------------------------------------------------
formula123 <- y ~ 0 + beta0 + beta1 + 
    f(s, model=spde) + f(e1, model='iid', hyper=pcprec) +
    f(eta1, model='iid',
      hyper=list(theta=list(initial=-10, fixed=TRUE))) + 
    f(eta2, model='iid',
      hyper=list(theta=list(initial=-10, fixed=TRUE))) + 
    f(eta1c, copy='eta1', fixed=FALSE) +
    f(e2, model='iid', hyper=pcprec) +
    f(eta2c, copy='eta2', fixed=FALSE) +
    f(e3, model='iid', hyper=pcprec)


## ----res123--------------------------------------------------------------
res123 <- inla(formula123,
               data=inla.stack.data(stack),
               offset=offset,
               control.family=list(list(
                   hyper=list(theta=list(initial=10, fixed=TRUE)))),
               control.predictor=list(A=inla.stack.A(stack)))


## ----fixed---------------------------------------------------------------
round(cbind(true=c(beta0,beta1),res123$summary.fixed), 4)


## ----copy12--------------------------------------------------------------
i.b <- match(paste0('Beta for eta', 1:2, 'c'),
             rownames(res123$summary.hyper))
round(cbind(true=c(beta2, beta3), res123$summary.hy[i.b,]), 4)


## ----s123,eval=FALSE-----------------------------------------------------
## i.e123 <- match(paste0('Precision for e', 1:3),
##                 names(res123$marginals.hyper))
## par(mfrow=c(1,3), mar=c(3,3,1,1), mgp=c(1.5, 0.5, 0))
## for (j in 1:3) {
##     plot(inla.tmarginal(function(x) 1/sqrt(x),
##                         res123$marginals.hyperpar[[i.e123[j]]]),
##          type='l', lwd=2, xlab=bquote(sigma[.(j)]^2),
##          ylab='Posterior marginal density')
##     abline(v=s123[j], lwd=2)
## }


## ----vs123,echo=F,results='hide',fig.width=7.5,fig.height=3--------------
i.e123 <- match(paste0('Precision for e', 1:3),
                names(res123$marginals.hyper))
par(mfrow=c(1,3), mar=c(3,3,1,1), mgp=c(1.5, 0.5, 0))
for (j in 1:3) {
    plot(inla.tmarginal(function(x) 1/sqrt(x),
                        res123$marginals.hyperpar[[i.e123[j]]]),
         type='l', lwd=2, xlab=bquote(sigma[.(j)]^2), 
         ylab='Posterior marginal density')
    abline(v=s123[j], lwd=2)
}


## ----rfparams,eval=FALSE-------------------------------------------------
## rfparams <- inla.spde2.result(res123, 's', spde)
## par(mfrow=c(1,3), mar=c(3,3,1,1), mgp=c(1.5, 0.5, 0))
## for (j in 15:17) {
##     plot(rfparams[[j]][[1]], xlab=names(rfparams)[j],
##          type='l', lwd=2, ylab='Posterior marginal density')
##     abline(v=c(10, 1, sqrt(8)/10)[j-14], lwd=2)
## }


## ----vrfparams,echo=F,results='hide',fig.width=7.5,fig.height=3----------
rfparams <- inla.spde2.result(res123, 's', spde)
par(mfrow=c(1,3), mar=c(3,3,1,1), mgp=c(1.5, 0.5, 0))
for (j in 15:17) {
    plot(rfparams[[j]][[1]], xlab=names(rfparams)[j],
         type='l', lwd=2, ylab='Posterior marginal density')
    abline(v=c(10, 1, sqrt(8)/10)[j-14], lwd=2)
}

