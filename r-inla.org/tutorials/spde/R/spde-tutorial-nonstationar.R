## ----settings,echo=FALSE,results='hide',message=FALSE,warning=FALSE------
library(knitr)
opts_chunk$set(
fig.path='figs/nonstationar',
message=FALSE, warning=FALSE
)
options(width=75, prompt = " ", continue = "   ")
library(lattice) 
library(INLA)
library(gridExtra)
lcall <- inla.getOption('inla.call')
###inla.setOption(inla.call='remote')
inla.setOption(num.threads=7)

## ----poly----------------------------------------------------------------
pl01 <- cbind(c(0,1,1,0,0), c(0,0,1,1,0))

## ----mesh----------------------------------------------------------------
(mesh <- inla.mesh.2d(, pl01, cutoff=0.03, 
                      max.edge=c(0.07,.12)))$n

## ----spde----------------------------------------------------------------
spde <- inla.spde2.matern(mesh, 
   B.tau=cbind(0, 1, 0, sin(pi*mesh$loc[,1])),
   B.kappa=cbind(0, 0, 1, 0), 
   theta.prior.mean=rep(0, 3),
   theta.prior.prec=rep(1, 3))

## ----thetas--------------------------------------------------------------
theta1 <- c(-1, 2, -1)
theta2 <- c(-1, 2, 1)

## ----Q-------------------------------------------------------------------
Q1 <- inla.spde2.precision(spde, theta=theta1)
Q2 <- inla.spde2.precision(spde, theta=theta2)

## ----changecalltolocal,echo=FALSE----------------------------------------
inla.setOption(inla.call=lcall)

## ----covs----------------------------------------------------------------
cov1 <- inla.qinv(Q1);         cov2 <- inla.qinv(Q2) 

## ----changetoremodecall,echo=FALSE---------------------------------------
inla.setOption(inla.call='remote')

## ----diagC---------------------------------------------------------------
v1 <- diag(cov1);      v2 <- diag(cov2)
rbind(v1=summary(v1),  v2=summary(v2))

## ----varns,eval=F--------------------------------------------------------
## par(mfrow=c(1,2), mar=c(3,3,.5,.5), mgp=c(2, .7, 0), las=1)
## plot(mesh$loc[,1], v1, ylim=range(v1,v2), las=1,
##      xlab='x-coordinate', ylab='marginal variance')
## points(mesh$loc[,1], v2, col=2)
## i1 <- which((mesh$loc[,1]>0) & (mesh$loc[,1]<1) &
##   (mesh$loc[,2]>0) & (mesh$loc[,2]<1))
## plot(mesh$loc[i1,1], v1[i1], ylim=range(v1[i1],v2[i1]),
##      xlab='x-coordinate', ylab='marginal variance', las=1)
## points(mesh$loc[i1,1], v2[i1], col=2)
## legend('topleft', as.expression(lapply(
##                       c(theta1[3], theta2[3]),
##                       function(x) bquote(theta[3]==.(x)))), col=1:2)

## ----vvarns,echo=FALSE,fig.width=10,heigh=5.5,width='0.97\\textwidth'----
par(mfrow=c(1,2), mar=c(3,3,.5,.5), mgp=c(2, .7, 0), las=1)
plot(mesh$loc[,1], v1, ylim=range(v1,v2), las=1, 
     xlab='x-coordinate', ylab='marginal variance')
points(mesh$loc[,1], v2, col=2)
i1 <- which((mesh$loc[,1]>0) & (mesh$loc[,1]<1) & 
  (mesh$loc[,2]>0) & (mesh$loc[,2]<1))
plot(mesh$loc[i1,1], v1[i1], ylim=range(v1[i1],v2[i1]), 
     xlab='x-coordinate', ylab='marginal variance', las=1)
points(mesh$loc[i1,1], v2[i1], col=2)
legend('topleft', as.expression(lapply(
                      c(theta1[3], theta2[3]), 
                      function(x) bquote(theta[3]==.(x)))), col=1:2)

## ----changecalltolocal2,echo=FALSE---------------------------------------
inla.setOption(inla.call=lcall)

## ----samples-------------------------------------------------------------
sample1 <-  as.vector(inla.qsample(1, Q1, seed=1))
sample2 <-  as.vector(inla.qsample(1, Q2, seed=1))

## ----changetoremodecall2,echo=FALSE--------------------------------------
inla.setOption(inla.call='remote')

## ----ssummary------------------------------------------------------------
tapply(sample1, round(inla.group(mesh$loc[,1], 5),3), var)
tapply(sample2, round(inla.group(mesh$loc[,1], 5),3), var)

## ----plotsamples,eval=F--------------------------------------------------
## proj <- inla.mesh.projector(mesh, xlim=0:1, ylim=0:1)
## grid.arrange(levelplot(inla.mesh.project(proj, field=sample1),
##                        xlab='', ylab='', scale=list(draw=FALSE),
##                        col.regions=topo.colors(100)),
##              levelplot(inla.mesh.project(proj,field=sample2),
##                        xlab='', ylab='', scale=list(draw=FALSE),
##                        col.regions=topo.colors(100)), nrow=1)

## ----vplotsamples,echo=FALSE,fig.width=10,heigh=4.5,width='0.97\\textwidth'----
proj <- inla.mesh.projector(mesh, xlim=0:1, ylim=0:1)
grid.arrange(levelplot(inla.mesh.project(proj, field=sample1), 
                       xlab='', ylab='', scale=list(draw=FALSE),
                       col.regions=topo.colors(100)), 
             levelplot(inla.mesh.project(proj,field=sample2), 
                       xlab='', ylab='', scale=list(draw=FALSE),
                       col.regions=topo.colors(100)), nrow=1)

## ----changecalltolocal3,echo=FALSE---------------------------------------
inla.setOption(inla.call=lcall)

## ----constrs-------------------------------------------------------------
s1r <-  as.vector(inla.qsample(1, Q1, seed=1, constr=spde$f$extraconstr))
s2r <-  as.vector(inla.qsample(1, Q2, seed=1, constr=spde$f$extraconstr))

## ----changetoremodecall3,echo=FALSE--------------------------------------
inla.setOption(inla.call='remote')

## ----comparss------------------------------------------------------------
rbind(s1r=summary(s1r), s2r=summary(s2r))
c(cor1=cor(sample1, s1r), cor2=cor(sample2, s2r))

## ----likehy--------------------------------------------------------------
clik <- list(hyper=list(theta=list(initial=20, fixed=TRUE)))

## ----fit12---------------------------------------------------------------
formula <- y ~ 0 + f(i, model=spde)
fit1 <- inla(formula, control.family=clik, 
             data=data.frame(y=sample1, i=1:mesh$n))
fit2 <- inla(formula, control.family=clik, 
             data=data.frame(y=sample2, i=1:mesh$n))

## ----hy1summaries--------------------------------------------------------
round(cbind(true=theta1, fit1$summary.hy), 4)

## ----hy2summaries--------------------------------------------------------
round(cbind(true=theta2, fit2$summary.hy), 4)

## ------------------------------------------------------------------------
set.seed(2);     n <- 100
loc <- cbind(runif(n), runif(n))

## ----projloc-------------------------------------------------------------
projloc <- inla.mesh.projector(mesh, loc)

## ----projectsamples------------------------------------------------------
x1 <- inla.mesh.project(projloc, sample1)
x2 <- inla.mesh.project(projloc, sample2)

## ----stacks--------------------------------------------------------------
stk1 <- inla.stack(list(y=x1), A=list(projloc$proj$A), tag='d',
                  effects=list(data.frame(i=1:mesh$n)))
stk2 <- inla.stack(list(y=x2), A=list(projloc$proj$A), tag='d',
                  effects=list(data.frame(i=1:mesh$n)))

## ----fitt----------------------------------------------------------------
res1 <- inla(formula, data=inla.stack.data(stk1), control.family=clik, 
             control.predictor=list(compute=TRUE, A=inla.stack.A(stk1)))
res2 <- inla(formula, data=inla.stack.data(stk2), control.family=clik, 
             control.predictor=list(compute=TRUE, A=inla.stack.A(stk2)))

## ----spostth-------------------------------------------------------------
round(cbind(True=theta1, res1$summary.hy), 4)
round(cbind(True=theta2, res2$summary.hy), 4)

## ----projpostrf----------------------------------------------------------
x1.mean <- inla.mesh.project(proj, field=res1$summary.ran$i$mean)
x1.var <- inla.mesh.project(proj, field=res1$summary.ran$i$sd^2)
x2.mean <- inla.mesh.project(proj, field=res2$summary.ran$i$mean)
x2.var <- inla.mesh.project(proj, field=res2$summary.ran$i$sd^2)

## ----visprojp,eval=F-----------------------------------------------------
## do.call(function(...) grid.arrange(..., nrow=2),
##         lapply(list(inla.mesh.project(proj, sample1), x1.mean, x1.var,
##                     inla.mesh.project(proj, sample2), x2.mean, x2.var),
##                levelplot, xlab='', ylab='',
##                col.regions=topo.colors(100), scale=list(draw=FALSE)))

## ----fvisprojp,echo=FALSE,fig.width=10,fig.height=5,width='0.97\\textwidth'----
do.call(function(...) grid.arrange(..., nrow=2), 
        lapply(list(inla.mesh.project(proj, sample1), x1.mean, x1.var, 
                    inla.mesh.project(proj, sample2), x2.mean, x2.var), 
               levelplot, xlab='', ylab='', 
               col.regions=topo.colors(100), scale=list(draw=FALSE)))

