## ----sett,echo=F,results='hide',message=FALSE,warning=FALSE--------------
library(knitr)
opts_chunk$set(
fig.path='figs/prefsampl',
message=FALSE, warning=FALSE
)
options(width=75, prompt = " ", continue = "   ")
library(INLA)
lcall <- inla.getOption('inla.call')
inla.setOption(inla.call='remote')
inla.setOption(num.threads=7)
library(splancs)
library(fields)
source('R/spde-tutorial-functions.R')

## ----window--------------------------------------------------------------
library(spatstat)
win <- owin(c(0,3), c(0,3))

## ----gridres-------------------------------------------------------------
spatstat.options(npixel=300)

## ------------------------------------------------------------------------
beta0 <- 3

## ----n-exp---------------------------------------------------------------
exp(beta0) * diff(range(win$x)) * diff(range(win$y))

## ------------------------------------------------------------------------
sigma2x <- 0.2;      kappa <- 2

## ----simulapp,eval=TRUE--------------------------------------------------
library(RandomFields)
set.seed(1)
lg.s <- rLGCP('matern', beta0, 
              var=sigma2x, scale=1/kappa, nu=1, win=win)

## ----xy------------------------------------------------------------------
xy <- cbind(lg.s$x, lg.s$y)[,2:1]

## ----nxy-----------------------------------------------------------------
(n <- nrow(xy))

## ------------------------------------------------------------------------
Lam <- attr(lg.s, 'Lambda')
summary(as.vector(rf.s <- log(Lam$v)))

## ----lgrfpp,eval=F-------------------------------------------------------
## par(mfrow=c(1,1))
## library(fields)
## image.plot(list(x=Lam$yrow, y=Lam$xcol, z=rf.s), main='log-Lambda', asp=1)
## points(xy, pch=19)

## ----echo=F,results='hide',fig.width=5,fig.height=5----------------------
par(mfrow=c(1,1))
library(fields)
image.plot(list(x=Lam$yrow, y=Lam$xcol, z=rf.s), main='log-Lambda', asp=1)
points(xy, pch=19) 

## ----mesh----------------------------------------------------------------
loc.d <- 3*t(matrix(c(0,0,1,0,1,1,0,1,0,0), 2))
(nv <- (mesh <- inla.mesh.2d(
            loc.domain=loc.d, offset=c(.3, 1), 
            max.edge=c(.3, 0.7), cutoff=.05))$n) 

## ----meshplot,eval=FALSE-------------------------------------------------
## par(mar=c(0,0,0,0))
## plot(mesh, asp=1, main='')
## points(xy, col=4, pch=19); lines(loc.d, col=3)

## ----echo=F,results='hide',fig.width=5,fig.height=5----------------------
par(mar=c(0,0,0,0))
plot(mesh, asp=1, main='')
points(xy, col=4, pch=19); lines(loc.d, col=3)

## ----spde----------------------------------------------------------------
spde <- inla.spde2.pcmatern(
    mesh=mesh, alpha=2, ### mesh and smoothness parameter
    prior.range=c(0.05, 0.01), ### P(practic.range<0.05)=0.01
    prior.sigma=c(1, 0.01)) ### P(sigma>1)=0.01

## ------------------------------------------------------------------------
sum(diag(spde$param.inla$M0))

## ----dualmesh------------------------------------------------------------
dmesh <- inla.mesh.dual(mesh)

## ----splocd--------------------------------------------------------------
domainSP <- SpatialPolygons(list(Polygons(
    list(Polygon(loc.d)), '0')))

## ----pols----------------------------------------------------------------
library(rgeos)
sum(w <- sapply(1:length(dmesh), function(i) {
    if (gIntersects(dmesh[i,], domainSP))
        return(gArea(gIntersection(dmesh[i,], domainSP)))
    else return(0)
}))

## ----wsummary------------------------------------------------------------
table(w>0)

## ----dualmeshcode,eval=FALSE---------------------------------------------
## par(mar=c(2,2,1,1), mgp=2:0)
## plot(mesh$loc, asp=1, col=(w==0)+1, pch=19, xlab='', ylab='')
## plot(dmesh, add=TRUE)
## lines(loc.d, col=3)

## ----echo=F,results='hide',fig.width=5,fig.height=5----------------------
par(mar=c(2,2,1,1), mgp=2:0)
plot(mesh$loc, asp=1, col=(w==0)+1, pch=19, xlab='', ylab='') 
plot(dmesh, add=TRUE)
lines(loc.d, col=3)

## ----y01-----------------------------------------------------------------
y.pp <- rep(0:1, c(nv, n))

## ----expected------------------------------------------------------------
e.pp <- c(w, rep(0, n)) 

## ----Aloc----------------------------------------------------------------
lmat <- inla.spde.make.A(mesh, xy)

## ----pp-proj-------------------------------------------------------------
imat <- Diagonal(nv, rep(1, nv))

## ----App-----------------------------------------------------------------
A.pp <- rBind(imat, lmat)

## ----stkpp---------------------------------------------------------------
stk.pp <- inla.stack(data=list(y=y.pp, e=e.pp), 
                     A=list(1,A.pp), tag='pp',
                     effects=list(list(b0=rep(1,nv+n)), list(i=1:nv))) 

## ----ppest---------------------------------------------------------------
pp.res <- inla(y ~ 0 + b0 + f(i, model=spde), 
               family='poisson', data=inla.stack.data(stk.pp), 
               control.predictor=list(A=inla.stack.A(stk.pp)), 
               E=inla.stack.data(stk.pp)$e)

## ----pppars,eval=TRUE----------------------------------------------------
round(pp.res$summary.hyperpar, 4)

## ----pp-viz,eval=F-------------------------------------------------------
## par(mfrow=c(1,3), mar=c(3,3,1,0.3), mgp=c(2,1,0))
## plot(pp.res$marginals.fix[[1]], type='l',
##      xlab=expression(beta[0]), ylab='Density')
## abline(v=beta0, col=2)
## plot(pp.res$marginals.hyperpar[[2]], type='l',
##      xlab=expression(sigma^2), ylab='Density')
## abline(v=sigma2x, col=2)
## plot(pp.res$marginals.hyperpar[[1]], type='l',
##      xlab='Nominal range', ylab='Density')
## abline(v=sqrt(8*1)/kappa, col=2)

## ----echo=F,results='hide',fig.width=10,fig.height=5---------------------
par(mfrow=c(1,3), mar=c(3,3,1,0.3), mgp=c(2,1,0)) 
plot(pp.res$marginals.fix[[1]], type='l', 
     xlab=expression(beta[0]), ylab='Density')
abline(v=beta0, col=2)
plot(pp.res$marginals.hyperpar[[2]], type='l', 
     xlab=expression(sigma^2), ylab='Density')
abline(v=sigma2x, col=2)
plot(pp.res$marginals.hyperpar[[1]], type='l', 
     xlab='Nominal range', ylab='Density')
abline(v=sqrt(8*1)/kappa, col=2)

## ----gridcov-------------------------------------------------------------
y0 <- x0 <- seq(win$xrange[1], win$xrange[2], 
                length=spatstat.options()$npixel)
gridcov <- outer(x0, y0, function(x,y) cos(x) - sin(y-2))

## ----n-exp-cov-----------------------------------------------------------
beta1 <- -0.5
sum(exp(beta0 + beta1*gridcov) * diff(x0[1:2])*diff(y0[1:2]))

## ----simulappc-----------------------------------------------------------
set.seed(1)
lg.s.c <- rLGCP('matern', im(beta0 + beta1*gridcov, xcol=x0, yrow=y0), 
              var=sigma2x, scale=1/kappa, nu=1, win=win)

## ----xyc-----------------------------------------------------------------
(n.c <- nrow(xy.c <- cbind(lg.s.c$x, lg.s.c$y)[,2:1]))

## ----lgrfppc,eval=FALSE--------------------------------------------------
## library(fields)
## par(mfrow=c(1,2), mar=c(2,2,1,1), mgp=c(1,0.5,0))
## image.plot(list(x=x0, y=y0, z=gridcov), main='Covariate', asp=1)
## image.plot(list(x=x0, y=y0, z=log(attr(lg.s.c, 'Lambda')$v)),
##            main='log-Lambda', asp=1)
## points(xy.c, pch=19)

## ----echo=F,results='hide',fig.width=10,fig.height=5---------------------
library(fields)
par(mfrow=c(1,2), mar=c(2,2,1,1), mgp=c(1,0.5,0))
image.plot(list(x=x0, y=y0, z=gridcov), main='Covariate', asp=1)
image.plot(list(x=x0, y=y0, z=log(attr(lg.s.c, 'Lambda')$v)), 
           main='log-Lambda', asp=1)
points(xy.c, pch=19)

## ----collcovar-----------------------------------------------------------
covariate = gridcov[Reduce('cbind', nearest.pixel(
    c(mesh$loc[,1], xy.c[,1]), c(mesh$loc[,2], xy.c[,2]), 
    im(gridcov, x0, y0)))]

## ----datc----------------------------------------------------------------
y.pp.c <- rep(0:1, c(nv, n.c))
e.pp.c <- c(w, rep(0, n.c))

## ----A.c-----------------------------------------------------------------
lmat.c <- inla.spde.make.A(mesh, xy.c)

## ----App.c---------------------------------------------------------------
A.pp.c <- rBind(imat, lmat.c)

## ----stkpp.c-------------------------------------------------------------
stk.pp.c <- inla.stack(data=list(y=y.pp.c, e=e.pp.c), 
                       A=list(1, A.pp.c), tag='pp.c',
                       effects=list(list(b0=1, covariate=covariate), 
                           list(i=1:nv)))

## ----ppest.c-------------------------------------------------------------
pp.c.res <- inla(y ~ 0 + b0 + covariate + f(i, model=spde), 
                 family='poisson', data=inla.stack.data(stk.pp.c), 
                 control.predictor=list(A=inla.stack.A(stk.pp.c)), 
                 E=inla.stack.data(stk.pp.c)$e)

## ----insp-u.c,eval=TRUE--------------------------------------------------
round(pp.c.res$summary.hyperpar, 4)

## ----ppc-viz,eval=F------------------------------------------------------
## par(mfrow=c(2,2), mar=c(3,3,1,0.3), mgp=c(2,1,0))
## plot(pp.c.res$marginals.fix[[1]], type='l', ylab='Density',
##      xlab=expression(beta[0])); abline(v=beta0, col=2)
## plot(pp.c.res$marginals.fix[[2]], type='l', ylab='Density',
##      xlab=expression(beta[1])); abline(v=beta1, col=2)
## plot(pp.c.res$marginals.hyperpar[[2]], type='l', ylab='Density',
##      xlab=expression(sigma)); abline(v=sigma2x^0.5, col=2)
## plot(pp.c.res$marginals.hyperpar[[1]], type='l', ylab='Density',
##      xlab=expression(sqrt(8)/kappa)); abline(v=kappa, col=2)

## ----echo=F,results='hide',fig.width=7.5,fig.height=7.5------------------
par(mfrow=c(2,2), mar=c(3,3,1,0.3), mgp=c(2,1,0)) 
plot(pp.c.res$marginals.fix[[1]], type='l', ylab='Density', 
     xlab=expression(beta[0])); abline(v=beta0, col=2)
plot(pp.c.res$marginals.fix[[2]], type='l', ylab='Density', 
     xlab=expression(beta[1])); abline(v=beta1, col=2)
plot(pp.c.res$marginals.hyperpar[[2]], type='l', ylab='Density', 
     xlab=expression(sigma)); abline(v=sigma2x^0.5, col=2)
plot(pp.c.res$marginals.hyperpar[[1]], type='l', ylab='Density',
     xlab=expression(sqrt(8)/kappa)); abline(v=kappa, col=2)

## ----simulaz-------------------------------------------------------------
summary(z <- log(t(Lam$v)[Reduce(
  'cbind', nearest.pixel(xy[,1], xy[,2], Lam))]))

## ----resp----------------------------------------------------------------
beta0.y <- 10;   beta <- -2;   prec.y <- 16
set.seed(2)
summary(resp <- beta0.y + (z-beta0)/beta + 
        rnorm(length(z), 0, sqrt(1/prec.y)))

## ----rresp---------------------------------------------------------------
stk.u <- inla.stack(data=list(y=resp), A=list(lmat, 1), 
                    effects=list(i=1:nv, b0=rep(1,length(resp))))
u.res <- inla(y ~ 0 + b0 + f(i, model=spde), 
              data=inla.stack.data(stk.u), 
              control.predictor=list(A=inla.stack.A(stk.u)))
round(cbind(True=c(beta0y=beta0.y, prec.y=prec.y), 
            rbind(u.res$summary.fix[, 1:6], u.res$summary.hy[1,])), 4)

## ----insp-u,eval=F-------------------------------------------------------
## par(mfrow=c(1,3), mar=c(3, 3, 0.3, 0.3), mgp=c(2,1,0))
## plot(inla.tmarginal(function(x) sqrt(1/x),
##                     u.res$marginals.hyperpar[[1]]),
##      type='l', ylab='Density', xlab=expression(sigma))
## abline(v=prec.y^-0.5, col=2)
## plot(u.res$marginals.hyperpar[[2]], type='l', ylab='Density',
##      xlab=expression(sqrt(8)/kappa)); abline(v=sqrt(8)/kappa, col=2)
## plot(u.res$marginals.hyperpar[[3]], type='l',
##      xlab=expression(sqrt(sigma^2[x])), ylab='Density')
## abline(v=sigma2x^0.5, col=2)

## ----echo=F,results='hide',fig.width=6,fig.height=2.1--------------------
par(mfrow=c(1,3), mar=c(3, 3, 0.3, 0.3), mgp=c(2,1,0))
plot(inla.tmarginal(function(x) sqrt(1/x), 
                    u.res$marginals.hyperpar[[1]]), 
     type='l', ylab='Density', xlab=expression(sigma)) 
abline(v=prec.y^-0.5, col=2)
plot(u.res$marginals.hyperpar[[2]], type='l', ylab='Density', 
     xlab=expression(sqrt(8)/kappa)); abline(v=sqrt(8)/kappa, col=2)
plot(u.res$marginals.hyperpar[[3]], type='l', 
     xlab=expression(sqrt(sigma^2[x])), ylab='Density')
abline(v=sigma2x^0.5, col=2)

## ----ppstk---------------------------------------------------------------
stk2.y <- inla.stack(data=list(y=cbind(resp,NA), e=rep(0,n)), 
                     A=list(lmat, 1), tag='resp2',
                     effects=list(i=1:nv, b0.y=rep(1,n)))
stk2.pp <- inla.stack(data=list(y=cbind(NA,y.pp), e=e.pp), 
                      A=list(A.pp, 1), tag='pp2',
                      effects=list(j=1:nv, b0.pp=rep(1,nv+n))) 
j.stk <- inla.stack(stk2.y, stk2.pp)

## ----j-res---------------------------------------------------------------
jform <- y ~ 0 + b0.pp + b0.y + 
    f(i, model=spde) + f(j, copy='i', fixed=FALSE)
j.res <- inla(jform, family=c('gaussian', 'poisson'), 
              data=inla.stack.data(j.stk), E=inla.stack.data(j.stk)$e,              
              control.predictor=list(A=inla.stack.A(j.stk)))
round(cbind(True=c(beta0, beta0.y), 
            j.res$summary.fix), 4)

## ----pppost,echo=F,results='hide',fig.width=5.5,fig.height=4-------------
par(mfrow=c(2,3), mar=c(3,3,0.5,0.5), mgp=c(2,0.5,0))
plot(j.res$marginals.fix[[1]], type='l', ylab='Density', 
     xlab=expression(beta[0]), lwd=2) 
lines(pp.res$marginals.fix[[1]], lty=2, lwd=2)
abline(v=beta0, col=2)
plot(j.res$marginals.fix[[2]], type='l', ylab='Density', 
     xlab=expression(beta[0]^y), lwd=2) 
lines(u.res$marginals.fix[[1]], lty=3, lwd=5)
abline(v=beta0.y, col=2)
plot(inla.tmarginal(function(x) 1/x, 
                    j.res$marginals.hy[[1]]), lwd=2, 
     type='l', ylab='Density', xlab=expression(sigma[y]^2)) 
lines(inla.tmarginal(function(x) 1/x, 
                     u.res$marginals.hy[[1]]), lwd=5, lty=3)
abline(v=1/prec.y, col=2)
plot(j.res$marginals.hyperpar[[2]], type='l', xlim=c(0,10),
     xlab=expression(sqrt(8)/kappa), ylab='Density', lwd=2)
lines(pp.res$marginals.hyperpar[[1]], lty=2, lwd=2)
lines(u.res$marginals.hyperpar[[2]], lty=3, lwd=5)
abline(v=sqrt(8)/kappa, col=2)
plot(j.res$marginals.hyperpar[[3]], type='l', lwd=2, xlim=c(0,1),
     xlab=expression(sqrt(sigma[x]^2)), ylab='Density')
lines(pp.res$marginals.hyperpar[[2]], lty=2, lwd=2)
lines(u.res$marginals.hyperpar[[3]], lty=3, lwd=5)
abline(v=sigma2x^0.5, col=2)
plot(j.res$marginals.hy[[4]], type='l', 
     xlab=expression(beta), ylab='Density', lwd=2)
abline(v=beta, col=2)
legend('topright', c('True value', 'Joint', 'Only PP', 'Only Y'), 
       col=c(2,1,1,1), lty=c(1,1,2,3), lwd=c(2,2,3,5), bty='n')

