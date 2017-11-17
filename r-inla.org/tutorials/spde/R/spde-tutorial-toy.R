## ----sett,echo=F,results='hide',message=FALSE,warning=FALSE--------------
library(knitr)
opts_chunk$set(
fig.path='figs/toy',
message=FALSE, warning=FALSE
)
options(width=75, prompt = " ", continue = "   ")
library(lattice) 
library(INLA)
lcall <- inla.getOption('inla.call')
##inla.setOption(inla.call='remote')
##inla.setOption(num.threads=7)

## ----datatoy-------------------------------------------------------------
data(SPDEtoy)

## ----strdata-------------------------------------------------------------
str(SPDEtoy)

## ----buildmesh5----------------------------------------------------------
pl.dom <- cbind(c(0,1,1,0.7,0), c(0,0,0.7,1,1))
mesh5 <- inla.mesh.2d(, pl.dom, max.e=c(0.092, 0.2))

## ----spde2matern---------------------------------------------------------
args(inla.spde2.matern)

## ----spde5def------------------------------------------------------------
spde5 <- inla.spde2.pcmatern(
    mesh=mesh5, alpha=2, ### mesh and smoothness parameter
    prior.range=c(0.3, 0.5), ### P(practic.range<0.3)=0.5
    prior.sigma=c(1, 0.01)) ### P(sigma>1)=0.01

## ----proj2---------------------------------------------------------------
coords <- as.matrix(SPDEtoy[,1:2])
A5 <- inla.spde.make.A(mesh5, loc=coords)

## ----dima1---------------------------------------------------------------
dim(A5)

## ----a5lines-------------------------------------------------------------
table(rowSums(A5>0))

## ----rsum----------------------------------------------------------------
table(rowSums(A5))

## ----colsA---------------------------------------------------------------
table(colSums(A5)>0)

## ----loadmeshes, echo=FALSE----------------------------------------------
for (i in 1:6)
  load(paste('mesh', i, '.RData', sep=''))

## ----eacha1--------------------------------------------------------------
A1 <- inla.spde.make.A(mesh1, loc=coords)
table(rowSums(A1>0))

## ----summarya1-----------------------------------------------------------
table(rowSums(A1))

## ----stackdata1b---------------------------------------------------------
stk5 <- inla.stack(data=list(resp=SPDEtoy$y), A=list(A5,1), 
                   effects=list(i=1:spde5$n.spde, 
                     m=rep(1,nrow(SPDEtoy))), tag='est')

## ----dimA----------------------------------------------------------------
dim(inla.stack.A(stk5))

## ----modelfit------------------------------------------------------------
res5 <- inla(resp ~ 0 + m + f(i, model=spde5), 
             data=inla.stack.data(stk5), 
             control.predictor=list(A=inla.stack.A(stk5)))

## ----beta0summary--------------------------------------------------------
res5$summary.fix

## ----invnuggetsummary----------------------------------------------------
res5$summary.hy[1,]

## ----postnugget----------------------------------------------------------
post.se <- inla.tmarginal(function(x) sqrt(1/x), 
                          res5$marginals.hy[[1]])

## ----summarypostnu-------------------------------------------------------
inla.emarginal(function(x) x, post.se)
inla.qmarginal(c(0.025, 0.5, 0.975), post.se)
inla.hpdmarginal(0.95, post.se)
inla.pmarginal(c(0.5, 0.7), post.se)

## ----variancepost--------------------------------------------------------
res5.field <- inla.spde2.result(res5, 'i', spde5, do.transf=TRUE)

## ----erandf--------------------------------------------------------------
inla.emarginal(function(x) x, res5.field$marginals.kappa[[1]])
inla.emarginal(function(x) x, res5.field$marginals.variance.nominal[[1]])
inla.emarginal(function(x) x, res5.field$marginals.range.nominal[[1]])

## ----pts3----------------------------------------------------------------
pts3 <- rbind(c(.1,.1), c(.5,.55), c(.7,.9))

## ----A5pts3--------------------------------------------------------------
dim(A5pts3 <- inla.spde.make.A(mesh5, loc=pts3))

## ----a5pts3c-------------------------------------------------------------
(jj3 <- which(colSums(A5pts3)>0))
round(A5pts3[, jj3],3)

## ----stk3prd-------------------------------------------------------------
stk5p.rf <- inla.stack(data=list(resp=NA), A=list(A5pts3), 
                       effects=list(i=1:spde5$n.spde), tag='prd5r')

## ----stakfull------------------------------------------------------------
stk5.jp <- inla.stack(stk5, stk5p.rf)

## ----refit---------------------------------------------------------------
res5p <- inla(resp ~ 0 + m + f(i, model=spde5), 
              data=inla.stack.data(stk5.jp), 
              control.predictor=list(A=inla.stack.A(stk5.jp), compute=TRUE))

## ----inddp---------------------------------------------------------------
(indd5p <- inla.stack.index(stk5.jp, tag='prd5r')$data)

## ----postp---------------------------------------------------------------
round(res5p$summary.linear.pred[indd5p,], 4)

## ----marg3p,results='hide'-----------------------------------------------
marg3 <- res5p$marginals.linear[indd5p]

## ----hpdp3---------------------------------------------------------------
inla.hpdmarginal(0.95, marg3[[2]])

## ----meanproj3-----------------------------------------------------------
drop(A5pts3%*%res5$summary.random$i$mean)

## ----projector-----------------------------------------------------------
inla.mesh.project(inla.mesh.projector(mesh5, loc=pts3), 
                  res5$summary.random$i$mean)

## ----sdproj3-------------------------------------------------------------
drop(A5pts3%*%res5$summary.random$i$sd)

## ----sdproj3c------------------------------------------------------------
sqrt(drop((A5pts3^2)%*%(res5$summary.random$i$sd^2)))

## ----grid0---------------------------------------------------------------
pgrid0 <- inla.mesh.projector(mesh5, xlim=0:1, ylim=0:1, dims=c(101,101))

## ----projg---------------------------------------------------------------
prd0.m <- inla.mesh.project(pgrid0,  res5$summary.ran$i$mean)
prd0.s <- inla.mesh.project(pgrid0,  res5$summary.ran$i$s)

## ----stackpresp----------------------------------------------------------
stk5.presp <- inla.stack(data=list(resp=NA), A=list(A5pts3,1), 
                         effects=list(i=1:spde5$n.spde, m=rep(1,3)), 
                         tag='prd5.resp')

## ----rresp---------------------------------------------------------------
stk5.full <- inla.stack(stk5, stk5.presp)
r5presp <- inla(resp ~ 0 + m + f(i, model=spde5), 
                data=inla.stack.data(stk5.full), 
                control.predictor=list(A=inla.stack.A(stk5.full), compute=TRUE))

## ----indd----------------------------------------------------------------
(indd3r <- inla.stack.index(stk5.full, 'prd5.resp')$data)

## ----postd---------------------------------------------------------------
round(r5presp$summary.fitted.values[indd3r,], 3)

## ----margp,results='hide'------------------------------------------------
marg3r <- r5presp$marginals.fitted.values[indd3r]

## ----hpdp----------------------------------------------------------------
inla.hpdmarginal(0.95, marg3r[[2]])

## ----prdmean-------------------------------------------------------------
res5$summary.fix[1,1] + drop(A5pts3%*%res5$summary.random$i$mean)

## ----optlcall,echo=F-----------------------------------------------------
inla.setOption(inla.call=lcall)

## ----cov,results='hide',echo=F-------------------------------------------
q <- inla.spde2.precision(spde5, theta=res5$summary.hyper[2:3,1])
rf.cov <- inla.qinv(q)
dim(rf.cov)
diag(cov3 <- A5pts3%*%rf.cov%*%t(A5pts3))
cov3

## ----interpmean----------------------------------------------------------
summary(rvar <- res5$summary.random$i$sd^2)
sqrt(1^2+res5$summary.fix[1,2]^2 + drop(A5pts3%*%rvar))

## ----prespgrid-----------------------------------------------------------
stkgrid <- inla.stack(data=list(resp=NA), A=list(pgrid0$proj$A,1), 
                      effects=list(i=1:spde5$n.spde,
                        m=rep(1,101*101)), tag='prd.gr')
stk.all <- inla.stack(stk5, stkgrid)
res5g <- inla(resp ~ 0 + m + f(i, model=spde5), 
              data=inla.stack.data(stk.all), 
              control.predictor=list(A=inla.stack.A(stk.all), 
                compute=TRUE), quantiles=NULL, 
              control.results=list(return.marginals.random=FALSE, 
                return.marginals.predictor=FALSE))
res5g$cpu

## ----indgr---------------------------------------------------------------
igr <- inla.stack.index(stk.all, 'prd.gr')$data

## ----pgrid,eval=F--------------------------------------------------------
## library(gridExtra)
## grid.arrange(levelplot(prd0.m, col.regions=topo.colors(99), main='latent field mean',
##                        xlab='', ylab='', scales=list(draw=FALSE)),
##              levelplot(matrix(res5g$summary.fitt[igr,1], 101),
##                        xlab='', ylab='', main='response mean',
##                        col.regions=topo.colors(99), scales=list(draw=FALSE)),
##              levelplot(prd0.s, col.regions=topo.colors(99), main='latent field SD',
##                        xlab='', ylab='', scales=list(draw=FALSE)),
##              levelplot(matrix(res5g$summary.fitt[igr,2], 101),
##                        xlab='', ylab='', main='response SD',
##                        col.regions=topo.colors(99), scales=list(draw=FALSE)),
##              nrow=2)

## ----echo=FALSE,fig.width=7.5,fig.height=7,out.width='0.97\\textwidth'----
library(gridExtra)
grid.arrange(levelplot(prd0.m, col.regions=topo.colors(99), main='latent field mean',
                       xlab='', ylab='', scales=list(draw=FALSE)), 
             levelplot(matrix(res5g$summary.fitt[igr,1], 101), 
                       xlab='', ylab='', main='response mean',
                       col.regions=topo.colors(99), scales=list(draw=FALSE)), 
             levelplot(prd0.s, col.regions=topo.colors(99), main='latent field SD',
                       xlab='', ylab='', scales=list(draw=FALSE)), 
             levelplot(matrix(res5g$summary.fitt[igr,2], 101), 
                       xlab='', ylab='', main='response SD',
                       col.regions=topo.colors(99), scales=list(draw=FALSE)), 
             nrow=2)

## ----meshes--------------------------------------------------------------
lrf <- lres <- l.dat <- l.spde <- l.a <- list()
for (k in 1:6) {
  l.a[[k]] <- inla.spde.make.A(get(paste('mesh', k, sep='')), loc=coords)
  l.spde[[k]] <- inla.spde2.matern(get(paste('mesh', k, sep='')), alpha=2)
  l.dat[[k]] <- list(y=SPDEtoy[,3], i=1:ncol(l.a[[k]]), 
                     m=rep(1, ncol(l.a[[k]])))
  lres[[k]] <- inla(y ~ 0 + m + f(i, model=l.spde[[k]]), 
                    data=l.dat[[k]], control.predictor=list(A=l.a[[k]]))
  lrf[[k]] <- inla.spde2.result(lres[[k]], 'i', l.spde[[k]], do.transf=TRUE)
}

## ----time----------------------------------------------------------------
round(sapply(lres, function(x) x$cpu[2]), 2)

## ----s2marg--------------------------------------------------------------
s2.marg <- lapply(lres, function(m) 
                  inla.tmarginal(function(x) 1/x, m$marginals.hy[[1]]))

## ----truepars------------------------------------------------------------
beta0 <- 10; sigma2e <- 0.3; sigma2u <- 5; kappa <- 7; nu <- 1

## ----lkv-----------------------------------------------------------------
lk.est <- c(beta=9.54, s2e=0.374, s2u=3.32, range=0.336)

## ----compare,eval=F------------------------------------------------------
## rcols <- rainbow(6)##c(rgb(4:1/4,0:3/5,0), c(rgb(0,0:3/5,4:1/4)))
## par(mfrow=c(2,3), mar=c(2.5,2.5,1,.5), mgp=c(1.5,.5,0), las=1)
## 
## xrange <- range(sapply(lres, function(x) range(x$marginals.fix[[1]][,1])))
## yrange <- range(sapply(lres, function(x) range(x$marginals.fix[[1]][,2])))
## plot(lres[[1]]$marginals.fix[[1]], type='l', xlim=xrange, ylim=yrange,
##      xlab=expression(beta[0]), ylab='Density')
## for (k in 1:6)
##   lines(lres[[k]]$marginals.fix[[1]], col=rcols[k], lwd=2)
## abline(v=beta0, lty=2, lwd=2, col=3)
## abline(v=lk.est[1], lty=3, lwd=2, col=3)
## 
## xrange <- range(sapply(s2.marg, function(x) range(x[,1])))
## yrange <- range(sapply(s2.marg, function(x) range(x[,2])))
## plot.default(s2.marg[[1]], type='l', xlim=xrange, ylim=yrange,
##              xlab=expression(sigma[e]^2), ylab='Density')
## for (k in 1:6)
##   lines(s2.marg[[k]], col=rcols[k], lwd=2)
## abline(v=sigma2e, lty=2, lwd=2, col=3)
## abline(v=lk.est[2], lty=3, lwd=2, col=3)
## 
## xrange <- range(sapply(lrf, function(r) range(r$marginals.variance.nominal[[1]][,1])))
## yrange <- range(sapply(lrf, function(r) range(r$marginals.variance.nominal[[1]][,2])))
## plot(lrf[[1]]$marginals.variance.nominal[[1]], type='l',
##      xlim=xrange, ylim=yrange, xlab=expression(sigma[x]^2), ylab='Density')
## for (k in 1:6)
##   lines(lrf[[k]]$marginals.variance.nominal[[1]], col=rcols[k], lwd=2)
## abline(v=sigma2u, lty=2, lwd=2, col=3)
## abline(v=lk.est[3], lty=3, lwd=2, col=3)
## 
## xrange <- range(sapply(lrf, function(r) range(r$marginals.kappa[[1]][,1])))
## yrange <- range(sapply(lrf, function(r) range(r$marginals.kappa[[1]][,2])))
## plot(lrf[[1]]$marginals.kappa[[1]], type='l',
##      xlim=xrange, ylim=yrange, xlab=expression(kappa), ylab='Density')
## for (k in 1:6)
##   lines(lrf[[k]]$marginals.kappa[[1]], col=rcols[k], lwd=2)
## abline(v=kappa, lty=2, lwd=2, col=3)
## abline(v=lk.est[4], lty=3, lwd=2, col=3)
## 
## xrange <- range(sapply(lrf, function(r) range(r$marginals.range.nominal[[1]][,1])))
## yrange <- range(sapply(lrf, function(r) range(r$marginals.range.nominal[[1]][,2])))
## plot(lrf[[1]]$marginals.range.nominal[[1]], type='l',
##      xlim=xrange, ylim=yrange, xlab='nominal range', ylab='Density')
## for (k in 1:6)
##   lines(lrf[[k]]$marginals.range.nominal[[1]], col=rcols[k], lwd=2)
## abline(v=sqrt(8)/kappa, lty=2, lwd=2, col=3)
## abline(v=sqrt(8)/lk.est[4], lty=3, lwd=2, col=3)
## 
## xrange <- range(sapply(lrf, function(r) range(r$marginals.tau[[1]][,1])))
## yrange <- range(sapply(lrf, function(r) range(r$marginals.tau[[1]][,2])))
## plot(lrf[[1]]$marginals.tau[[1]], type='l',
##      xlim=xrange, ylim=yrange, xlab=expression(tau), ylab='Density')
## for (k in 1:6)
##   lines(lrf[[k]]$marginals.tau[[1]], col=rcols[k], lwd=2)
## 
## legend('topright', c(paste('mesh', 1:6, sep=''), 'True', 'Likelihood'),
##        lty=c(rep(1,6), 2, 3), lwd=rep(2, 6), col=c(rcols,3,3), bty='n')

## ----kappai--------------------------------------------------------------
1/kappa

## ----echo=F,fig.width=7.5,fig.height=5-----------------------------------
rcols <- rainbow(6)##c(rgb(4:1/4,0:3/5,0), c(rgb(0,0:3/5,4:1/4)))
par(mfrow=c(2,3), mar=c(2.5,2.5,1,.5), mgp=c(1.5,.5,0), las=1)

xrange <- range(sapply(lres, function(x) range(x$marginals.fix[[1]][,1])))
yrange <- range(sapply(lres, function(x) range(x$marginals.fix[[1]][,2])))
plot(lres[[1]]$marginals.fix[[1]], type='l', xlim=xrange, ylim=yrange, 
     xlab=expression(beta[0]), ylab='Density')
for (k in 1:6)
  lines(lres[[k]]$marginals.fix[[1]], col=rcols[k], lwd=2)
abline(v=beta0, lty=2, lwd=2, col=3) 
abline(v=lk.est[1], lty=3, lwd=2, col=3)

xrange <- range(sapply(s2.marg, function(x) range(x[,1])))
yrange <- range(sapply(s2.marg, function(x) range(x[,2])))
plot.default(s2.marg[[1]], type='l', xlim=xrange, ylim=yrange, 
             xlab=expression(sigma[e]^2), ylab='Density')
for (k in 1:6) 
  lines(s2.marg[[k]], col=rcols[k], lwd=2)
abline(v=sigma2e, lty=2, lwd=2, col=3) 
abline(v=lk.est[2], lty=3, lwd=2, col=3)

xrange <- range(sapply(lrf, function(r) range(r$marginals.variance.nominal[[1]][,1])))
yrange <- range(sapply(lrf, function(r) range(r$marginals.variance.nominal[[1]][,2])))
plot(lrf[[1]]$marginals.variance.nominal[[1]], type='l', 
     xlim=xrange, ylim=yrange, xlab=expression(sigma[x]^2), ylab='Density')
for (k in 1:6)
  lines(lrf[[k]]$marginals.variance.nominal[[1]], col=rcols[k], lwd=2)
abline(v=sigma2u, lty=2, lwd=2, col=3) 
abline(v=lk.est[3], lty=3, lwd=2, col=3)

xrange <- range(sapply(lrf, function(r) range(r$marginals.kappa[[1]][,1])))
yrange <- range(sapply(lrf, function(r) range(r$marginals.kappa[[1]][,2])))
plot(lrf[[1]]$marginals.kappa[[1]], type='l', 
     xlim=xrange, ylim=yrange, xlab=expression(kappa), ylab='Density')
for (k in 1:6)
  lines(lrf[[k]]$marginals.kappa[[1]], col=rcols[k], lwd=2)
abline(v=kappa, lty=2, lwd=2, col=3) 
abline(v=lk.est[4], lty=3, lwd=2, col=3)

xrange <- range(sapply(lrf, function(r) range(r$marginals.range.nominal[[1]][,1])))
yrange <- range(sapply(lrf, function(r) range(r$marginals.range.nominal[[1]][,2])))
plot(lrf[[1]]$marginals.range.nominal[[1]], type='l', 
     xlim=xrange, ylim=yrange, xlab='nominal range', ylab='Density')
for (k in 1:6)
  lines(lrf[[k]]$marginals.range.nominal[[1]], col=rcols[k], lwd=2)
abline(v=sqrt(8)/kappa, lty=2, lwd=2, col=3)
abline(v=sqrt(8)/lk.est[4], lty=3, lwd=2, col=3)

xrange <- range(sapply(lrf, function(r) range(r$marginals.tau[[1]][,1])))
yrange <- range(sapply(lrf, function(r) range(r$marginals.tau[[1]][,2])))
plot(lrf[[1]]$marginals.tau[[1]], type='l', 
     xlim=xrange, ylim=yrange, xlab=expression(tau), ylab='Density')
for (k in 1:6)
  lines(lrf[[k]]$marginals.tau[[1]], col=rcols[k], lwd=2)

legend('topright', c(paste('mesh', 1:6, sep=''), 'True', 'Likelihood'), 
       lty=c(rep(1,6), 2, 3), lwd=rep(2, 6), col=c(rcols,3,3), bty='n')

