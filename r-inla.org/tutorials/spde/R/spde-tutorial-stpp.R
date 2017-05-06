## ----opts,echo=F,results='hide',message=FALSE,warning=FALSE--------------
library(knitr)
opts_chunk$set(
fig.path='figs/stpp',
message=FALSE, warning=FALSE
)
options(width=75, prompt = " ", continue = "   ")
library(INLA)
lcall <- inla.getOption('inla.call')
inla.setOption(inla.call='remote')
inla.setOption(num.threads=7)
source('R/spde-tutorial-functions.R')

## ----sdomain-------------------------------------------------------------
x0 <- seq(0, 4*pi, length=15)
domain <- data.frame(x=c(x0, rev(x0), 0))
domain$y <- c(sin(x0/2)-2, sin(rev(x0/2))+2, sin(0)-2)

## ----sp------------------------------------------------------------------
library(sp)
domainSP <- SpatialPolygons(list(Polygons(list(Polygon(domain)), '0')))

## ----lgcp,results='hide'-------------------------------------------------
library(lgcp) 
ndays <- 15
n <- (xyt <- lgcpSim(
    owin=spatstat:::owin(poly=domain), tlim=c(0,ndays), 
    model.parameters=lgcppars(1,0.5,0.1,0,0.5), cellwidth=0.1,
    spatial.covmodel='matern', covpars=c(nu=1)))$n

## ----tmesh---------------------------------------------------------------
k <- 7; tmesh <- inla.mesh.1d(seq(0, ndays, length=k))

## ----mesh----------------------------------------------------------------
smesh <- inla.mesh.2d(boundary=inla.sp2segment(domainSP),
                      max.edge=1, cutoff=0.3)

## ----stplot--------------------------------------------------------------
par(mfrow=c(2,1), mar=c(1.5,0,0,0), mgp=c(1,0.5,0))
plot(sample(xyt$t,500), rep(1,500), type='h', ylim=0:1,
     xlab='Day', ylab='', axes=FALSE); box(); axis(1)
abline(v=tmesh$loc, col=4, lwd=3)
par(mar=c(0,0,0,0))
plot(smesh, asp=1, main='')
points(xyt$x, xyt$y, cex=0.5, pch=3)

## ----stfit,echo=FALSE,fig.width=10,fig.height=4--------------------------
par(mfrow=c(2,1), mar=c(1.5,0,0,0), mgp=c(1,0.5,0))
plot(sample(xyt$t,500), rep(1,500), type='h', ylim=0:1,
     xlab='Day', ylab='', axes=FALSE); box(); axis(1)
abline(v=tmesh$loc, col=4, lwd=3)
par(mar=c(0,0,0,0))
plot(smesh, asp=1, main='')
points(xyt$x, xyt$y, cex=0.5, pch=3)

## ----voronoi-------------------------------------------------------------
library(deldir)
dd <- deldir(smesh$loc[,1], smesh$loc[,2])
tiles <- tile.list(dd)

## ----sppls---------------------------------------------------------------
polys <- SpatialPolygons(lapply(1:length(tiles), function(i)
    { p <- cbind(tiles[[i]]$x, tiles[[i]]$y)
      n <- nrow(p)
      Polygons(list(Polygon(p[c(1:n, 1),])), i)
  }))

## ----ppinsp--------------------------------------------------------------
area <- factor(over(SpatialPoints(cbind(xyt$x, xyt$y)), 
                    polys), levels=1:length(polys))

## ----tnear---------------------------------------------------------------
t.breaks <- sort(c(tmesh$loc[c(1,k)],
                   tmesh$loc[2:k-1]/2 + tmesh$loc[2:k]/2))
table(time <- factor(findInterval(xyt$t, t.breaks), 
                     levels=1:(length(t.breaks)-1)))

## ----agg-----------------------------------------------------------------
agg.dat <- as.data.frame(table(area, time))
for(j in 1:2) ### set time and area as integer
    agg.dat[[j]] <- as.integer(as.character(agg.dat[[j]])) 
str(agg.dat)

## ----intersect-----------------------------------------------------------
library(rgeos)
sum(w.areas <- sapply(1:length(tiles), function(i)
    { p <- cbind(tiles[[i]]$x, tiles[[i]]$y)
      n <- nrow(p)
      pl <- SpatialPolygons(list(Polygons(list(Polygon(p[c(1:n, 1),])), i)))
      if (gIntersects(pl, domainSP))
          return(gArea(gIntersection(pl, domainSP)))
      else return(0)
  }))

## ----plarea--------------------------------------------------------------
summary(w.areas)

## ----sarea---------------------------------------------------------------
gArea(domainSP)

## ----wtknot--------------------------------------------------------------
(w.t <- diag(inla.mesh.fem(tmesh)$c0))

## ----intensity0----------------------------------------------------------
(i0 <- n / (gArea(domainSP) * diff(range(tmesh$loc))))

## ----stvol---------------------------------------------------------------
summary(e0 <- w.areas[agg.dat$area] * (w.t[agg.dat$time]))

## ----spde----------------------------------------------------------------
A.st <- inla.spde.make.A(smesh, smesh$loc[agg.dat$area,],
                         group=agg.dat$time, mesh.group=tmesh)
spde <- inla.spde2.matern(smesh)
idx <- inla.spde.make.index('s', spde$n.spde, n.group=k)

## ----stack---------------------------------------------------------------
stk <- inla.stack(data=list(y=agg.dat$Freq, exposure=e0), 
                  A=list(A.st, 1), 
                  effects=list(idx, 
                      list(b0=rep(1, nrow(agg.dat)))))

## ----formula-------------------------------------------------------------
formula <- y ~ 0 + b0 + 
    f(s, model=spde, group=s.group, control.group=list(model='ar1'))

## ----fitt----------------------------------------------------------------
res <- inla(formula, family='poisson', 
            data=inla.stack.data(stk), E=exposure, 
            control.predictor=list(A=inla.stack.A(stk)),
            control.inla=list(strategy='gaussian'))

## ----intercept-----------------------------------------------------------
round(cbind(true=log(i0), res$summary.fixed),4)

## ----n-------------------------------------------------------------------
eta.i <- res$summary.fix[1,1] + res$summary.ran$s$mean
c(n=xyt$n, 'E(n)'=sum(rep(w.areas, k)*rep(w.t, each=smesh$n)*exp(eta.i)))

## ----lstsres-------------------------------------------------------------
r0 <- diff(range(domain[,1]))/diff(range(domain[,2]))
prj <- inla.mesh.projector(smesh, xlim=bbox(domainSP)[1,], 
                           ylim=bbox(domainSP)[2,], dims=c(r0*200, 200))
g.no.in <- is.na(over(SpatialPoints(prj$lattice$loc), domainSP))
t.mean <- lapply(1:k, function(j) {
    z <- inla.mesh.project(prj, res$summary.ran$s$mean[idx$s.group==j])
    z[g.no.in] <- NA
    return(z)
})

## ----lstsresp,eval=FALSE-------------------------------------------------
## zlims <- range(unlist(t.mean), na.rm=TRUE)
## library(fields)
## par(mfrow=c(4,2), mar=c(0.1,0.1,0.1,0.1))
## for (j in 1:k) {
##     image(prj$x, prj$y, t.mean[[j]],
##           axes=FALSE, zlim=zlims, col=tim.colors(30))
##     points(xyt$x[time==j], xyt$y[time==j], cex=0.1)
## }
## image.plot(prj$x, prj$y, t.mean[[j]]+1e9, axes=FALSE, zlim=zlims, xlab='',
##            legend.mar=10, legend.width=5, col=tim.colors(30), horizontal=T)

## ----lstsresf,echo=FALSE,fig.width=5,fig.height=5,out.width='0.97\\textwidth'----
zlims <- range(unlist(t.mean), na.rm=TRUE)
library(fields)
par(mfrow=c(4,2), mar=c(0.1,0.1,0.1,0.1))
for (j in 1:k) {
    image(prj$x, prj$y, t.mean[[j]], 
          axes=FALSE, zlim=zlims, col=tim.colors(30))
    points(xyt$x[time==j], xyt$y[time==j], cex=0.1)
}
image.plot(prj$x, prj$y, t.mean[[j]]+1e9, axes=FALSE, zlim=zlims, xlab='', 
           legend.mar=10, legend.width=5, col=tim.colors(30), horizontal=T)

