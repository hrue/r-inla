## ----opts,echo=F,results='hide',message=FALSE,warning=FALSE--------------
library(knitr)
opts_chunk$set(
fig.path='figs/survival',
message=FALSE, warning=FALSE
)
options(width=75, prompt = " ", continue = "   ")
library(INLA)
library(fields)
lcall <- inla.getOption('inla.call')
inla.setOption(inla.call='remote')
inla.setOption(num.threads=8)

## ----Leuk----------------------------------------------------------------
data(Leuk);    sapply(Leuk, summary)

## ----mesh----------------------------------------------------------------
loc <- cbind(Leuk$xcoord, Leuk$ycoord)
bnd1 <- inla.nonconvex.hull(loc, convex=0.05)
bnd2 <- inla.nonconvex.hull(loc, convex=0.25)
mesh <- inla.mesh.2d(loc, boundary=list(bnd1, bnd2), 
    max.edge=c(0.05, 0.2), cutoff=0.005)

## ----proj----------------------------------------------------------------
A <- inla.spde.make.A(mesh, loc)

## ----spde----------------------------------------------------------------
spde <- inla.spde2.pcmatern(
    mesh=mesh, alpha=2, ### mesh and smoothness parameter
    prior.range=c(0.05, 0.01), ### P(practic.range<0.05)=0.01
    prior.sigma=c(1, 0.01)) ### P(sigma>1)=0.01

## ----form----------------------------------------------------------------
formula <- inla.surv(time, cens) ~ 0 + a0 +
    sex + age + wbc + tpi +
    f(spatial, model=spde)

## ----stack---------------------------------------------------------------
stk <- inla.stack(data=list(time=Leuk$time, cens=Leuk$cens),
                  A=list(A, 1),
                  effect=list(
                      list(spatial=1:spde$n.spde),
                      data.frame(a0=1, Leuk[,-c(1:4)])))

## ----res,results='hide'--------------------------------------------------
r <- inla(formula, family="weibullsurv", data=inla.stack.data(stk), 
          control.predictor=list(A=inla.stack.A(stk)))

## ----fix-----------------------------------------------------------------
round(r$summary.fix, 4)

## ----hy------------------------------------------------------------------
round(r$summary.hy, 3)

## ----prj-----------------------------------------------------------------
r0 <- diff(range(bbox(nwEngland)[1,]))/diff(range(bbox(nwEngland)[2,]))
prj <- inla.mesh.projector(mesh, xlim=bbox(nwEngland)[1,], 
                           ylim=bbox(nwEngland)[2,], 
                           dims=c(200*r0, 200))

## ----nas-----------------------------------------------------------------
m.spat <- inla.mesh.project(prj, r$summary.ran$spatial$mean)
sd.spat <- inla.mesh.project(prj, r$summary.ran$spatial$sd)
ov <- over(SpatialPoints(prj$lattice$loc), nwEngland)
sd.spat[is.na(ov)] <- m.spat[is.na(ov)] <- NA

## ----wscode,eval=FALSE---------------------------------------------------
## par(mfrow=c(1,2), mar=c(0,0,0,0))
## image.plot(x=prj$x, y=prj$y, z=m.spat, asp=1,
##            xlab='', ylab='', axes=FALSE, horizontal=TRUE)
## plot(nwEngland, add=TRUE)
## image.plot(x=prj$x, y=prj$y, z=sd.spat, asp=1,
##            xlab='', ylab='', axes=FALSE, horizontal=TRUE)
## plot(nwEngland, add=TRUE)

## ----wsmap,echo=FALSE,fig.width=7,fig.height=5.3,out.width='0.97\\textwidth'----
par(mfrow=c(1,2), mar=c(0,0,0,0)) 
image.plot(x=prj$x, y=prj$y, z=m.spat, asp=1, 
           xlab='', ylab='', axes=FALSE, horizontal=TRUE)
plot(nwEngland, add=TRUE)
image.plot(x=prj$x, y=prj$y, z=sd.spat, asp=1, 
           xlab='', ylab='', axes=FALSE, horizontal=TRUE)
plot(nwEngland, add=TRUE)

## ----coxphdat------------------------------------------------------------
formula0 <- inla.surv(time, cens) ~ 0 + a0 + sex + age + wbc + tpi 
cph.leuk <- inla.coxph(formula0, data=data.frame(a0=1, Leuk[, c(1:8)]),
                       control.hazard=list(n.intervals=25))

## ----coxph0,results='hide'-----------------------------------------------
cph.res0 <- inla(formula0, 'coxph', data=data.frame(a0=1,Leuk))

## ------------------------------------------------------------------------
cph.formula <- update(cph.leuk$formula, 
                      '. ~ . + f(spatial, model=spde)')

## ----Acph----------------------------------------------------------------
cph.A <- inla.spde.make.A(mesh, loc=cbind(
   cph.leuk$data$xcoord, cph.leuk$data$ycoord))

## ----stkcph--------------------------------------------------------------
cph.stk <- inla.stack(data=c(list(E=cph.leuk$E), 
                          cph.leuk$data[c('y..coxph')]),
                      A=list(cph.A, 1),
                      effects=list(
                          list(spatial=1:spde$n.spde), 
                          cph.leuk$data[c('baseline.hazard', 'a0',
                                          'age', 'sex', 'wbc', 'tpi')]))
cph.data <- c(inla.stack.data(cph.stk), cph.leuk$data.list)

## ----cphres,results='hide'-----------------------------------------------
cph.res <- inla(cph.formula, family='Poisson', 
                data=cph.data, E=cph.data$E, 
                control.predictor=list(A=inla.stack.A(cph.stk)))

## ----survival------------------------------------------------------------
library(survival)
m0 <- coxph(Surv(time, cens) ~ sex + age + wbc + tpi, Leuk)
cbind(survival=c(NA, coef(m0)),
      r0=cph.res0$summary.fix[,1], r1=cph.res$summary.fix[,1])

## ----corwcph-------------------------------------------------------------
cor(as.vector(m.spat), 
    as.vector(inla.mesh.project(
        prj, cph.res$summary.ran$spatial$mean)), use='p')

cor(log(as.vector(sd.spat)), 
    log(as.vector(inla.mesh.project(
        prj, cph.res$summary.ran$spatial$sd))), use='p')

