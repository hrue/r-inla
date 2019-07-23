## ----sett,echo=FALSE,results='hide',message=FALSE,warning=FALSE----------
library(knitr)
opts_chunk$set(
fig.path='figs/spacetime',
message=FALSE, warning=FALSE
)
options(width=75, prompt = " ", continue = "   ")
library(lattice) 
library(INLA)
lcall <- inla.getOption('inla.call')
## inla.setOption(inla.call='remote')
## inla.setOption(num.threads=4)
source('R/spde-tutorial-functions.R')
load('prmesh1.RData')


## ----data----------------------------------------------------------------
data(PRborder)


## ----k-------------------------------------------------------------------
k <- 12


## ----prprec--------------------------------------------------------------
data(PRprec)
coords <- as.matrix(PRprec[sample(1:nrow(PRprec)), 1:2])


## ----echo=FALSE,results='hide'-------------------------------------------
inla.setOption(inla.call=lcall)


## ----rk------------------------------------------------------------------
params <- c(variance=1, kappa=1) 
set.seed(1)
x.k <- rspde(coords, kappa=params[2], variance=params[1], n=k, 
             mesh=prmesh1, return.attributes=TRUE)
dim(x.k)


## ----beta----------------------------------------------------------------
rho <- 0.7


## ----x-------------------------------------------------------------------
x <- x.k
for (j in 2:k) 
    x[,j] <- rho*x[,j-1] + sqrt(1-rho^2)*x.k[,j]


## ----timevisual,eval=FALSE-----------------------------------------------
## c100 <- rainbow(101)
## par(mfrow=c(4,3), mar=c(0,0,0,0))
## for (j in 1:k)
##   plot(coords, col=c100[round(100*(x[,j]-min(x[,j]))/diff(range(x[,j])))],
##        axes=FALSE, asp=1, pch=19, cex=0.5)


## ----ftimevisual,echo=FALSE,eval=TRUE,fig.width=7.5,fig.height=7---------
c100 <- rainbow(101)
par(mfrow=c(4,3), mar=c(0,0,0,0))
for (j in 1:k)
  plot(coords, col=c100[round(100*(x[,j]-min(x[,j]))/diff(range(x[,j])))], 
       axes=FALSE, asp=1, pch=19, cex=0.5)


## ----categcov------------------------------------------------------------
n <- nrow(coords)
set.seed(2)
table(ccov <- factor(sample(LETTERS[1:3], n*k, replace=TRUE)) )


## ----betacov-------------------------------------------------------------
beta <- -1:1


## ----respst--------------------------------------------------------------
sd.y <- 0.1
y <- beta[unclass(ccov)] + x + rnorm(n*k, 0, sd.y)
tapply(y, ccov, mean)


## ----seldat--------------------------------------------------------------
isel <- sample(1:(n*k), n*k/2) 


## ----dat-----------------------------------------------------------------
dat <- data.frame(y=as.vector(y), w=ccov, 
                  time=rep(1:k, each=n), 
                  xcoo=rep(coords[,1], k), 
                  ycoo=rep(coords[,2], k))[isel, ] 


## ----spde----------------------------------------------------------------
spde <- inla.spde2.pcmatern(
    mesh=prmesh1, alpha=2, ### mesh and smoothness parameter
    prior.range=c(0.5, 0.01), ### P(practic.range<0.05)=0.01
    prior.sigma=c(1, 0.01)) ### P(sigma>1)=0.01


## ----rfindex-------------------------------------------------------------
iset <- inla.spde.make.index('i', n.spde=spde$n.spde, n.group=k)


## ----apred---------------------------------------------------------------
A <- inla.spde.make.A(mesh=prmesh1, 
                      loc=cbind(dat$xcoo, dat$ycoo), 
                      group=dat$time) 


## ----stack---------------------------------------------------------------
sdat <- inla.stack(tag='stdata', data=list(y=dat$y), 
                   A=list(A,1),  effects=list(iset, w=dat$w)) 


## ----hbeta---------------------------------------------------------------
h.spec <- list(theta=list(prior='pccor1', param=c(0, 0.9)))


## ----remote,echo=FALSE---------------------------------------------------
## inla.setOption(inla.call='remote')

## ----ft------------------------------------------------------------------
formulae <- y ~ 0 + w + 
    f(i, model=spde, group=i.group, 
      control.group=list(model='ar1', hyper=h.spec)) 
prec.prior <- list(prior='pc.prec', param=c(1, 0.01))
res <- inla(formulae,  data=inla.stack.data(sdat), 
            control.predictor=list(compute=TRUE, A=inla.stack.A(sdat)), 
            control.family=list(hyper=list(theta=prec.prior)), 
            control.fixed=list(expand.factor.strategy='inla'))


## ----sbeta---------------------------------------------------------------
round(cbind(observed=tapply(dat$y, dat$w, mean), res$summary.fixed), 4) 


## ----rfst,eval=F---------------------------------------------------------
## par(mfrow=c(2,2), mar=c(3,3,1,0.1), mgp=2:0)
## for (j in 1:4) {
##     plot(res$marginals.hyper[[j]], type='l',
##          xlab=names(res$marginals.hyper)[j], ylab='Density')
##     abline(v=c(1/sd.y^2, sqrt(8)/params[1],
##                params[2]^0.5, rho)[j], col=2)
## }


## ----echo=FALSE,fig.width=5.5,fig.height=5.5-----------------------------
par(mfrow=c(2,2), mar=c(3,3,1,0.1), mgp=2:0)
for (j in 1:4) {
    plot(res$marginals.hyper[[j]], type='l', 
         xlab=names(res$marginals.hyper)[j], ylab='Density')
    abline(v=c(1/sd.y^2, sqrt(8)/params[1], 
               params[2]^0.5, rho)[j], col=2)
}


## ----rfidx---------------------------------------------------------------
str(idat <- inla.stack.index(sdat, 'stdata')$data) 


## ----meanrf--------------------------------------------------------------
cor(dat$y, res$summary.linear.predictor$mean[idat])


## ----projgrid------------------------------------------------------------
stepsize <- 4*1/111
nxy <- round(c(diff(range(coords[,1])), diff(range(coords[,2])))/stepsize)
projgrid <- inla.mesh.projector(prmesh1, xlim=range(coords[,1]), 
                                ylim=range(coords[,2]), dims=nxy)


## ----projpmean-----------------------------------------------------------
xmean <- list()
for (j in 1:k)
  xmean[[j]] <- inla.mesh.project(
    projgrid, res$summary.random$i$mean[iset$i.group==j])


## ----inout---------------------------------------------------------------
library(splancs)
xy.in <- inout(projgrid$lattice$loc, cbind(PRborder[,1], PRborder[,2]))


## ----setNA---------------------------------------------------------------
for (j in 1:k)   xmean[[j]][!xy.in] <- NA


## ----strf,eval=FALSE-----------------------------------------------------
## library(gridExtra)
## do.call(function(...) grid.arrange(..., nrow=4),
##         lapply(xmean, levelplot,  xlab='', ylab='',
##                col.regions=topo.colors(16), scale=list(draw=FALSE)))


## ----vstrf,eval=TRUE,echo=FALSE,fig.width=10,fig.height=8.7--------------
library(gridExtra)
do.call(function(...) grid.arrange(..., nrow=4), 
        lapply(xmean, levelplot,  xlab='', ylab='',
               col.regions=topo.colors(16), scale=list(draw=FALSE)))


## ----stkv----------------------------------------------------------------
vdat <- data.frame(r=as.vector(y), w=ccov, t=rep(1:k, each=n), 
                   x=rep(coords[,1], k), y=rep(coords[,2], k))[-isel, ] 
Aval <- inla.spde.make.A(prmesh1, loc=cbind(vdat$x, vdat$y), group=vdat$t) 
stval <- inla.stack(tag='stval', data=list(y=NA), ### set NA in order to predict
                    A=list(Aval,1),  effects=list(iset, w=vdat$w)) 


## ----val-----------------------------------------------------------------
stfull <- inla.stack(sdat, stval) 
vres <- inla(formulae,  data=inla.stack.data(stfull), 
             control.predictor=list(compute=TRUE, A=inla.stack.A(stfull)), 
             control.family=list(hyper=list(theta=prec.prior)), 
             control.fixed=list(expand.factor.strategy='inla'), 
             control.mode=list(theta=res$mode$theta, restart=FALSE))


## ----indval--------------------------------------------------------------
ival <- inla.stack.index(stfull, 'stval')$data 


## ----stvalplot,eval=FALSE------------------------------------------------
## par(mfrow=c(1,1), mar=c(3,3,0.5,0.5), mgp=c(1.75,0.5,0))
## plot(vres$summary.fitted.values$mean[ival], vdat$r,
##      asp=1, xlab='Posterior mean', ylab='Observed')
## abline(0:1, col=gray(.7))


## ----vstval,echo=FALSE,fig.width=5.5,fig.height=5------------------------
par(mfrow=c(1,1), mar=c(3,3,0.5,0.5), mgp=c(1.75,0.5,0))
plot(vres$summary.fitted.values$mean[ival], vdat$r, 
     asp=1, xlab='Posterior mean', ylab='Observed') 
abline(0:1, col=gray(.7)) 


## ----loctime-------------------------------------------------------------
n <- nrow(loc <- unique(as.matrix(PRprec[,1:2])))
time <- sort(runif(n, 0, 1))


## ----stcov---------------------------------------------------------------
stcov <- function(coords, time, kappa.s, kappa.t, variance=1, nu=1) {
    s <- as.matrix(dist(coords))
    t <- as.matrix(dist(time))
    scorr <- exp((1-nu)*log(2) + nu*log(s*kappa.s) - lgamma(nu)) * 
        besselK(s*kappa.s, nu) 
    diag(scorr) <- 1
    return(variance * scorr * exp(-t*kappa.t))
}


## ----sample--------------------------------------------------------------
kappa.s <- 1; kappa.t <- 5; s2 <- 1/2
xx <- crossprod(chol(stcov(loc, time, kappa.s, kappa.t, s2)), 
                rnorm(n))  
beta0 <- -3; tau.error <- 3
y <- beta0 + xx + rnorm(n, 0, sqrt(1/tau.error))


## ----tknots--------------------------------------------------------------
k <- 10
(mesh.t <- inla.mesh.1d(seq(0+0.5/k, 1-0.5/k, length=k)))$loc


## ----rfind---------------------------------------------------------------
iset <- inla.spde.make.index('i', n.spde=spde$n.spde, n.group=k)


## ----apredt--------------------------------------------------------------
A <- inla.spde.make.A(mesh=prmesh1, loc=loc, 
                      group=time, group.mesh=mesh.t) 


## ----stackst-------------------------------------------------------------
sdat <- inla.stack(tag='stdata', data=list(y=y), 
                   A=list(A,1),  effects=list(iset, list(b0=rep(1,n))))


## ----tcor1---------------------------------------------------------------
exp(-kappa.t*diff(mesh.t$loc[1:2]))


## ----cfit----------------------------------------------------------------
formulae <- y ~ 0 + b0 + 
    f(i, model=spde, group=i.group, 
      control.group=list(model='ar1', hyper=h.spec)) 
res <- inla(formulae, data=inla.stack.data(sdat), 
            control.family=list(hyper=list(theta=prec.prior)), 
            control.predictor=list(compute=TRUE, A=inla.stack.A(sdat)))


## ----rfpars--------------------------------------------------------------
round(res$summary.hyper, 4)


## ----stcpost,eval=F------------------------------------------------------
## par(mfrow=c(2,3), mar=c(3,3,1,0.1), mgp=2:0)
## plot(res$marginals.fixed[[1]], type='l',
##      xlab=expression(beta[0]), ylab='Density')
## abline(v=beta0, col=2)
## for (j in 1:4) {
##     plot(res$marginals.hyper[[j]], type='l',
##          xlab=names(res$marginals.hyper)[j], ylab='Density')
##     abline(v=c(tau.error, sqrt(8)/kappa.s, sqrt(s2),
##              exp(-kappa.t*diff(mesh.t$loc[1:2])))[j], col=2)
## }


## ----echo=FALSE,fig.width=7.5,fig.height=5-------------------------------
par(mfrow=c(2,3), mar=c(3,3,1,0.1), mgp=2:0)
plot(res$marginals.fixed[[1]], type='l', 
     xlab=expression(beta[0]), ylab='Density')
abline(v=beta0, col=2)
for (j in 1:4) {
    plot(res$marginals.hyper[[j]], type='l', 
         xlab=names(res$marginals.hyper)[j], ylab='Density')
    abline(v=c(tau.error, sqrt(8)/kappa.s, sqrt(s2), 
             exp(-kappa.t*diff(mesh.t$loc[1:2])))[j], col=2)
}

