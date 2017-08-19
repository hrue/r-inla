## ----settings,include=FALSE,results='hide',message=FALSE,warning=FALSE----
library(knitr)
opts_chunk$set(
fig.path='figs/rain',
message=FALSE, warning=FALSE
)
options(width=75, prompt = " ", continue = "   ")
library(INLA)
lcall <- inla.getOption('inla.call')
inla.setOption(inla.call='remote')
inla.setOption(num.threads=7)
library(lattice) 
library(gridExtra) 

## ----start,results='hide'------------------------------------------------
data(PRprec) 

## ----headdat-------------------------------------------------------------
PRprec[ii <- c(which(is.na(PRprec$A))[which.min(PRprec$La[is.na(PRprec$A)])], 
               which(PRprec$Lo%in%range(PRprec$Lo)), which.max(PRprec$A)), 1:10]

## ----meanjan-------------------------------------------------------------
summary(PRprec$precMean <- rowMeans(PRprec[,3+1:31], na.rm=TRUE) )
table(rowSums(is.na(PRprec[,3+1:31])))

## ----prborder------------------------------------------------------------
data(PRborder) 

## ----precviz,eval=F,results='hide'---------------------------------------
## par(mfrow=c(1,2), mar=c(0,0,2,0))
## plot(PRborder, type='l', asp=1, axes=FALSE, main='Altitude')
## points(PRprec[1:2], col=is.na(PRprec$Alt)+1,
##        cex=ifelse(is.na(PRprec$Alt), 1, .3+PRprec$Alt/1500))
## legend('topright', format(0:4*350), bty='n', pch=1, pt.cex=.3+0:4*35/150)
## lines(PRborder[1034:1078, ], col='cyan')
## 
## plot(PRborder, type='l', asp=1, axes=FALSE,
##      main=paste('Mean of daily accumulated precipitation (mm)'))
## points(PRprec[1:2], cex=0.3+PRprec$precMean/20)
## legend('topright', format(seq(1,21,5)),
##        bty='n', pch=1, pt.cex=0.3+seq(1,21,5)/20)
## points(PRprec[ii, 1:2], pch=3, col=2)
## lines(PRborder[1034:1078, ], col='cyan')

## ----prmap,echo=FALSE,results='hide',fig.width=10,fig.height=4-----------
par(mfrow=c(1,2), mar=c(0,0,2,0)) 
plot(PRborder, type='l', asp=1, axes=FALSE, main='Altitude') 
points(PRprec[1:2], col=is.na(PRprec$Alt)+1, 
       cex=ifelse(is.na(PRprec$Alt), 1, .3+PRprec$Alt/1500)) 
legend('topright', format(0:4*350), bty='n', pch=1, pt.cex=.3+0:4*35/150) 
lines(PRborder[1034:1078, ], col='cyan') 

plot(PRborder, type='l', asp=1, axes=FALSE, 
     main=paste('Mean of daily accumulated precipitation (mm)')) 
points(PRprec[1:2], cex=0.3+PRprec$precMean/20)  
legend('topright', format(seq(1,21,5)), 
       bty='n', pch=1, pt.cex=0.3+seq(1,21,5)/20) 
points(PRprec[ii, 1:2], pch=3, col=2) 
lines(PRborder[1034:1078, ], col='cyan') 

## ----project-------------------------------------------------------------
coords <- as.matrix(PRprec[,1:2]) 
mat.dists <- spDists(coords, PRborder[1034:1078,], longlat=TRUE) 

## ----distseacalc---------------------------------------------------------
PRprec$"seaDist" <- apply(mat.dists, 1, min) 

## ----dispp,eval=F,results='hide'-----------------------------------------
## par(mfrow=c(2,2), mar=c(3,3,0.5,0.5), mgp=c(1.7,.7,0), las=1)
## for (i in c(1:3, ncol(PRprec)))  plot(PRprec[c(i,ncol(PRprec)-1)], cex=.5)

## ----dispfig,echo=F,results='hide',fig.width=7.5,fig.height=5------------
par(mfrow=c(2,2), mar=c(3,3,0.5,0.5), mgp=c(1.7,.7,0), las=1) 
for (i in c(1:3, ncol(PRprec)))  plot(PRprec[c(i,ncol(PRprec)-1)], cex=.5) 

## ----pcprec--------------------------------------------------------------
pcprec <- list(prior='pcprec', param=c(1, 0.01))

## ----prmesh--------------------------------------------------------------
pts.bound <- inla.nonconvex.hull(coords, 0.3, 0.3)
mesh <- inla.mesh.2d(coords, boundary=pts.bound, 
                     max.edge=c(0.3,1), offset=c(1e-5,1.5), cutoff=0.1)

## ----projA---------------------------------------------------------------
A <- inla.spde.make.A(mesh, loc=coords)

## ----spde----------------------------------------------------------------
spde <- inla.spde2.pcmatern(
    mesh=mesh, alpha=2, ### mesh and smoothness parameter
    prior.range=c(0.05, 0.01), ### P(practic.range<0.05)=0.01
    prior.sigma=c(1, 0.01)) ### P(sigma>1)=0.01

## ----stackpr-------------------------------------------------------------
stk.dat <- inla.stack(
    data=list(y=PRprec$precMean), 
    A=list(A,1), tag='dat', 
    effects=list(list(s=1:spde$n.spde), 
                 data.frame(Intercept=1, 
                            gWest=inla.group(coords[,1]), 
                            gSeaDist=inla.group(PRprec$seaDist), 
                            seaDist=PRprec$seaDist))) 

## ----fitwestcoo----------------------------------------------------------
f.west <- y ~ 0 + Intercept + 
    f(gWest, model='rw1', ### first random walk prior 
      scale.model=TRUE, ### scaling this prior
      hyper=list(theta=pcprec)) + ### considering the PC prior
    f(s, model=spde)
r.west <- inla(f.west, family='Gamma', 
               control.compute=list(cpo=TRUE),
               data=inla.stack.data(stk.dat), 
               control.predictor=list(
                   A=inla.stack.A(stk.dat), link=1))

## ----fitdistsea----------------------------------------------------------
f.seaD <- y ~ 0 + Intercept + 
    f(gSeaDist, model='rw1', scale.model=TRUE, 
      hyper=list(theta=pcprec)) + 
    f(s, model=spde)
r.seaD <- inla(f.seaD, family='Gamma', 
               control.compute=list(cpo=TRUE),
               data=inla.stack.data(stk.dat), 
               control.predictor=list(
                   A=inla.stack.A(stk.dat), link=1))

## ----fitdistseal---------------------------------------------------------
f.seaD.l <- y ~ 0 + Intercept + seaDist + 
    f(s, model=spde)
r.seaD.l <- inla(f.seaD.l, family='Gamma', 
                control.compute=list(cpo=TRUE),
                data=inla.stack.data(stk.dat), 
                control.predictor=list(
                    A=inla.stack.A(stk.dat), link=1))

## ----cpos----------------------------------------------------------------
slcpo <- function(m, na.rm=TRUE) 
    -sum(log(m$cpo$cpo), na.rm=na.rm)
c(long=slcpo(r.west), seaD=slcpo(r.seaD), 
  seaD.l=slcpo(r.seaD.l))

## ----betasumary1---------------------------------------------------------
round(r.seaD.l$summary.fixed, 4)

## ----disp----------------------------------------------------------------
round(unlist(r.seaD.l$summary.hy[1,]), 4)

## ----resfield------------------------------------------------------------
round(r.seaD.l$summary.hyper[-1, ], 4)

## ----seacoefs,eval=F-----------------------------------------------------
## par(mfrow=c(2,3), mar=c(3,3.5,0,0), mgp=c(1.5, .5, 0), las=0)
## plot(r.seaD.l$marginals.fix[[1]], type='l',
##      xlab='Intercept', ylab='Density')
## plot(r.seaD.l$marginals.fix[[2]], type='l',
##      xlab='Sea distance coefficient', ylab='Density')
## plot(r.seaD$summary.random[[1]][,1:2], type='l',
##      xlab='Distance to sea (Km)', ylab='Effect')
## abline(h=0, lty=3)
## for (i in c(4,6))
##     lines(r.seaD$summary.random[[1]][,c(1,i)], lty=2)
## abline(h=0)
## for (j in 1:3)
## plot(r.seaD.l$marginals.hy[[j]], type='l',
##      ylab='Density', xlab=names(r.seaD.l$marginals.hy)[j])

## ----seacoefsv,echo=F,results='hide',fig.width=7.5,fig.height=4----------
par(mfrow=c(2,3), mar=c(3,3.5,0,0), mgp=c(1.5, .5, 0), las=0) 
plot(r.seaD.l$marginals.fix[[1]], type='l', 
     xlab='Intercept', ylab='Density') 
plot(r.seaD.l$marginals.fix[[2]], type='l', 
     xlab='Sea distance coefficient', ylab='Density') 
plot(r.seaD$summary.random[[1]][,1:2], type='l', 
     xlab='Distance to sea (Km)', ylab='Effect') 
abline(h=0, lty=3) 
for (i in c(4,6)) 
    lines(r.seaD$summary.random[[1]][,c(1,i)], lty=2) 
abline(h=0)
for (j in 1:3) 
plot(r.seaD.l$marginals.hy[[j]], type='l', 
     ylab='Density', xlab=names(r.seaD.l$marginals.hy)[j])

## ----mod0----------------------------------------------------------------
r0.seaD.l <- inla(y ~ 0 + Intercept + seaDist, 
                  family='Gamma', control.compute=list(cpo=TRUE),
                  data=inla.stack.data(stk.dat), 
                  control.predictor=list(A=inla.stack.A(stk.dat), link=1))

## ----cpo0----------------------------------------------------------------
c(seaD.l=slcpo(r.seaD.l), seaD.l0=slcpo(r0.seaD.l))

## ----stepsize------------------------------------------------------------
(stepsize <- 4*1/111)

## ----ncoords-------------------------------------------------------------
(nxy <- round(c(diff(range(PRborder[,1])), 
                diff(range(PRborder[,2])))/stepsize))

## ----projgrid------------------------------------------------------------
projgrid <- inla.mesh.projector(mesh, xlim=range(PRborder[,1]), 
                                ylim=range(PRborder[,2]), dims=nxy)

## ----projpred------------------------------------------------------------
xmean <- inla.mesh.project(projgrid, r.seaD$summary.random$s$mean)
xsd <- inla.mesh.project(projgrid, r.seaD$summary.random$s$sd)

## ----sp------------------------------------------------------------------
library(splancs)
table(xy.in <- inout(projgrid$lattice$loc, 
                     cbind(PRborder[,1], PRborder[,2])))
xmean[!xy.in] <- xsd[!xy.in] <- NA

## ------------------------------------------------------------------------
Aprd <- projgrid$proj$A[which(xy.in), ]

## ----prdcoo--------------------------------------------------------------
prdcoo <- projgrid$lattice$loc[which(xy.in),]

## ----seaDcov0------------------------------------------------------------
seaDist0 <- apply(spDists(PRborder[1034:1078,], 
                          prdcoo, longlat=TRUE), 2, min)

## ----seaDistk------------------------------------------------------------
seaDist.k <- sort(unique(stk.dat$effects$data$gSeaDist))

## ----gseaD---------------------------------------------------------------
seaDist.b <- (seaDist.k[-1] + seaDist.k[length(seaDist.k)])/2
i0 <- findInterval(seaDist0, seaDist.b)+1
gSeaDist0 <- seaDist.k[i0]

## ----stkprd--------------------------------------------------------------
stk.prd <- inla.stack(
    data=list(y=NA), A=list(Aprd,1), 
    effects=list(s=1:spde$n.spde, 
                 data.frame(Intercept=1, 
                            seaDist=seaDist0)), tag='prd') 
stk.all <- inla.stack(stk.dat, stk.prd)

## ----predfit-------------------------------------------------------------
r2.seaD.l <- inla(f.seaD.l, family='Gamma', 
                  data=inla.stack.data(stk.all), 
                  control.predictor=list(A=inla.stack.A(stk.all), 
                                         compute=TRUE, link=1), 
                  quantiles=NULL, 
                  control.inla=list(strategy='gaussian'), 
                  control.results=list(return.marginals.random=FALSE,
                                       return.marginals.predictor=FALSE), 
                  control.mode=list(theta=r.seaD.l$mode$theta, 
                                    restart=FALSE))

## ----rprevprep-----------------------------------------------------------
id.prd <- inla.stack.index(stk.all, 'prd')$data
sd.prd <- m.prd <- matrix(NA, nxy[1], nxy[2])
m.prd[xy.in] <- r2.seaD.l$summary.fitted.values$mean[id.prd]
sd.prd[xy.in] <- r2.seaD.l$summary.fitted.values$sd[id.prd]

## ----xrain1c,eval=F------------------------------------------------------
## library(gridExtra)
## do.call('grid.arrange',
##         lapply(list(xmean, xsd, m.prd, sd.prd),
##                levelplot, col.regions=terrain.colors(16),
##                xlab='', ylab='', scales=list(draw=FALSE)))

## ----xrain1,echo=F,results='hide',fig.width=10,fig.height=7--------------
library(gridExtra) 
do.call('grid.arrange', 
        lapply(list(xmean, xsd, m.prd, sd.prd), 
               levelplot, col.regions=terrain.colors(16), 
               xlab='', ylab='', scales=list(draw=FALSE)))

## ----dsmesh--------------------------------------------------------------
seaDist.mesh <- apply(
    spDists(PRborder[1034:1078,], 
            mesh$loc[,1:2], longlat=TRUE), 2, min)

## ----stkmesh-------------------------------------------------------------
stk.mesh <- inla.stack(
    tag='mesh', data=list(y=NA), A=list(1,1),  
    effects=list(s=1:spde$n.spde, 
                 data.frame(Intercept=1, seaDist=seaDist.mesh)))
stk.b <- inla.stack(stk.dat, stk.mesh)

## ----fittmesh------------------------------------------------------------
rm.seaD.l <- inla(f.seaD.l, family='Gamma', 
                  data=inla.stack.data(stk.b), 
                  control.predictor=list(A=inla.stack.A(stk.b), 
                                         compute=TRUE, link=1), 
                  quantiles=NULL, 
                  control.results=list(return.marginals.random=FALSE,
                                       return.marginals.predictor=FALSE), 
                  control.compute=list(config=TRUE)) ## need to sample

## ----calll,echo=F,results='hide'-----------------------------------------
inla.setOption(inla.call=lcall)

## ----sampl---------------------------------------------------------------
sampl <- inla.posterior.sample(n=1000, result=rm.seaD.l)

## ----idexs---------------------------------------------------------------
dim(pred.nodes <- exp(sapply(sampl, function(x) 
    x$latent[nrow(PRprec) + 1:nrow(mesh$loc)])))

## ----projsamples---------------------------------------------------------
sd.prd.s <- m.prd.s <- matrix(NA, nxy[1], nxy[2])
m.prd.s[xy.in] <- drop(Aprd %*% rowMeans(pred.nodes))
sd.prd.s[xy.in] <- drop(Aprd %*% apply(pred.nodes, 1, sd))

## ----comparepreds--------------------------------------------------------
cor(as.vector(m.prd.s), as.vector(m.prd), use='p')
cor(log(as.vector(sd.prd.s)), log(as.vector(sd.prd)), use='p')

