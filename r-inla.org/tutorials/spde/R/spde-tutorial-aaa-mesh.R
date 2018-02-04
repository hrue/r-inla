## ----sett,include=FALSE,message=FALSE,warning=FALSE----------------------
library(knitr)
opts_chunk$set(
fig.path='figs/mesh',
message=FALSE, warning=FALSE
)
options(width=75, prompt = " ", continue = "   ")
library(INLA)

## ----argsmesh------------------------------------------------------------
args(inla.mesh.2d)

## ----SPDEtoy-------------------------------------------------------------
data(SPDEtoy)
coords <- as.matrix(SPDEtoy[,1:2]) ;   p5 <- coords[1:5,]

## ----domain--------------------------------------------------------------
pl.dom <- cbind(c(0,1,1,0.7,0), c(0,0,0.7,1,1))

## ----mesh5ab-------------------------------------------------------------
m1 <- inla.mesh.2d(p5, max.edge=c(0.5, 0.5)) 
m2 <- inla.mesh.2d(p5, max.edge=c(0.5, 0.5), cutoff=0.1) 
m3 <- inla.mesh.2d(p5, max.edge=c(0.1, 0.5), cutoff=0.1) 
m4 <- inla.mesh.2d(p5, max.edge=c(0.1, 0.5), offset=c(0,-0.65)) 
m5 <- inla.mesh.2d(, pl.dom, max.edge=c(0.3, 0.5), offset=c(0.03, 0.5)) 
m6 <- inla.mesh.2d(, pl.dom, max.edge=c(0.3, 0.5), offset=c(0.03, 0.5), cutoff=0.1)
m7 <- inla.mesh.2d(, pl.dom, max.edge=c(0.3, 0.5), n=5, offset=c(.05,.1)) 
m8 <- inla.mesh.2d(, pl.dom, max.edge=c(.3, 0.5), n=7, offset=c(.01,.3)) 
m9 <- inla.mesh.2d(, pl.dom, max.edge=c(.3, 0.5), n=4, offset=c(.05,.3)) 

## ----vizmesh,eval=FALSE--------------------------------------------------
## par(mfrow=c(3, 3), mar=c(0,0,1,0))
## for (i in 1:9) {
##   plot(pl.dom, type='l', col=3, lwd=2*(i>4), xlim=c(-0.57,1.57),
##        main = paste('m',i,sep=''), asp=1, axes=FALSE)
##   plot(get(paste('m', i, sep='')), add=TRUE)
##   points(p5, pch=19, col=2)
## }

## ----vvizmesh,echo=F,results='hide',fig.width=5.5,fig.height=5.5,out.width='0.97\\textwidth'----
par(mfrow=c(3, 3), mar=c(0,0,1,0))
for (i in 1:9) { 
  plot(pl.dom, type='l', col=3, lwd=2*(i>4), xlim=c(-0.57,1.57), 
       main = paste('m',i,sep=''), asp=1, axes=FALSE)
  plot(get(paste('m', i, sep='')), add=TRUE) 
  points(p5, pch=19, col=2)
}

## ----meshclass-----------------------------------------------------------
class(m1)
names(m1)

## ----n-------------------------------------------------------------------
c(m1$n, m2$n, m3$n, m4$n, m5$n, m6$n, m7$n, m8$n, m9$n)

## ----A1------------------------------------------------------------------
dim(m1$graph$vv)

## ----meshid--------------------------------------------------------------
m1$idx$loc

## ----noncovex------------------------------------------------------------
args(inla.nonconvex.hull)

## ----nonconmesh----------------------------------------------------------
bound1 <- inla.nonconvex.hull(p5)
bound2 <- inla.nonconvex.hull(p5, convex=0.5, concave=-0.15)
bound3 <- inla.nonconvex.hull(p5, concave=0.5)
bound4 <- inla.nonconvex.hull(p5, concave=0.5, resolution=c(20, 20))

m10 <- inla.mesh.2d(boundary=bound1, cutoff=0.05, max.edge=c(.1,.2))
m11 <- inla.mesh.2d(boundary=bound2, cutoff=0.05, max.edge=c(.1,.2))
m12 <- inla.mesh.2d(boundary=bound3, cutoff=0.05, max.edge=c(.1,.2))
m13 <- inla.mesh.2d(boundary=bound4, cutoff=0.05, max.edge=c(.1,.2))

## ----nonconmeshv,eval=FALSE----------------------------------------------
## par(mfrow=c(2,2), mar=c(0,0,1,0))
## for (i in 10:13) {
##    plot(get(paste('m', i, sep='')), asp=1, main='')
##    points(p5, pch=19, col=2); title(main=paste('m', i, sep=''))
## }

## ----vnonconmeshv,echo=F,results='hide',fig.width=5.5,fig.height=5.5,out.width='0.7\\textwidth'----
par(mfrow=c(2,2), mar=c(0,0,1,0))
for (i in 10:13) { 
   plot(get(paste('m', i, sep='')), asp=1, main='') 
   points(p5, pch=19, col=2); title(main=paste('m', i, sep=''))
}

## ----defnc---------------------------------------------------------------
max(diff(range(p5[,1])), diff(range(p5[,2])))*.15

## ----mesh12--------------------------------------------------------------
mesh1 <- inla.mesh.2d(coords, max.edge=c(0.035, 0.1)) 
mesh2 <- inla.mesh.2d(coords, max.edge=c(0.15, 0.2)) 

## ----mesh3---------------------------------------------------------------
mesh3 <- inla.mesh.2d(coords, max.edge=c(0.15, 0.2), cutoff=0.02)

## ----mesh456-------------------------------------------------------------
mesh4 <- inla.mesh.2d(, pl.dom, max.e=c(0.0355, 0.1))
mesh5 <- inla.mesh.2d(, pl.dom, max.e=c(0.092, 0.2))
mesh6 <- inla.mesh.2d(, pl.dom, max.e=c(0.11, 0.2))

## ----nmesh---------------------------------------------------------------
c(mesh1$n, mesh2$n, mesh3$n, mesh4$n, mesh5$n, mesh6$n)

## ----plotmesh1,eval=FALSE------------------------------------------------
## par(mfrow=c(2,3), mar=c(0,0,0,0))
## for (i in 1:6)
##   plot(get(paste('mesh',i,sep='')), asp=1, main='')

## ----vplotmesh1,fig.width=7.5,fig.height=5,echo=F,results='hide',out.width='0.97\\textwidth'----
par(mfrow=c(2,3), mar=c(0,0,0,0)) 
for (i in 1:6) 
  plot(get(paste('mesh',i,sep='')), asp=1, main='')

## ----prrain--------------------------------------------------------------
data(PRprec); dim(PRprec)
PRprec[1:2, 1:10]

## ----prpl----------------------------------------------------------------
data(PRborder); dim(PRborder)

## ----prmeshnch-----------------------------------------------------------
prdomain <- inla.nonconvex.hull(as.matrix(PRprec[,1:2]), 
                                -0.03, -0.05, resolution=c(100,100))

## ----mesh2pr-------------------------------------------------------------
(prmesh1 <- inla.mesh.2d(boundary=prdomain, max.edge=c(.7,.7), 
                         cutoff=0.35, offset=c(-0.05, -0.05)))$n
(prmesh2 <- inla.mesh.2d(boundary=prdomain, max.edge=c(.45,1), cutoff=0.2))$n

## ----prmesh,eval=FALSE---------------------------------------------------
## par(mfrow=c(1,2), mar=c(0,0,0,0))
## plot(prmesh1, asp=1, main='');   lines(PRborder, col=3)
## plot(prmesh2, asp=1, main='');   lines(PRborder, col=3)

## ----vprmesh,echo=F,results='hide',fig.width=7.5,fig.height=3.5,out.width='0.97\\textwidth'----
par(mfrow=c(1,2), mar=c(0,0,0,0))
plot(prmesh1, asp=1, main='');   lines(PRborder, col=3)
plot(prmesh2, asp=1, main='');   lines(PRborder, col=3)

## ----ncmap---------------------------------------------------------------
library(maptools)
nc.fl <- system.file("etc/shapes/sids.shp", package="spdep")[1]
nc.sids <- readShapePoly(nc.fl, ID="FIPSNO", 
                         proj4string=CRS("+proj=longlat +ellps=clrk66"))

## ----uionsp,results='hide',eval=FALSE------------------------------------
## gpclibPermit()

## ----unionSpatialPolygons------------------------------------------------
nc.border <- unionSpatialPolygons(nc.sids, rep(1, nrow(nc.sids)))

## ----ncsegment-----------------------------------------------------------
nc.bdry <- inla.sp2segment(nc.border)

## ----ncmesh--------------------------------------------------------------
(nc.mesh <- inla.mesh.2d(boundary=nc.bdry, cutoff=0.15, 
                         max.edge=c(0.3, 1)))$n

## ----ncmeshplot,eval=FALSE-----------------------------------------------
## par(mar=c(0,0,0,0))
## plot(nc.mesh, asp=1, main='')

## ----ncmeshv,echo=FALSE,fig.width=8,fig.height=3-------------------------
par(mar=c(0,0,0,0)) 
plot(nc.mesh, asp=1, main='')

## ----savemesh,echo=F,results='hide'--------------------------------------
save('mesh1', file='mesh1.RData')
save('mesh2', file='mesh2.RData')
save('mesh3', file='mesh3.RData')
save('mesh4', file='mesh4.RData')
save('mesh5', file='mesh5.RData')
save('mesh6', file='mesh6.RData')
save('prmesh1', file='prmesh1.RData')
save('prmesh2', file='prmesh2.RData')

## ----hexample,eval=FALSE-------------------------------------------------
## pl1 <- Polygon(cbind(c(0,15,15,0,0), c(5,0,20,20,5)), hole=FALSE)
## h1 <- Polygon(cbind(c(5,7,7,5,5), c(7,7,15,15,7)), hole=TRUE)
## pl2 <- Polygon(cbind(c(15,20,20,30,30,15,15), c(10,10,0,0,20,20,10)), hole=FALSE)
## sp <- SpatialPolygons(list(Polygons(list(pl1, h1), '0'), Polygons(list(pl2), '1')))
## par(mar=c(0,0,0,0)); plot(sp) ### to visualize it
## text(c(13, 17, 23), c(3, 12, 3), LETTERS[1:3], cex=3)

## ----plotlake,echo=F,results='hide',fig.width=7,fig.height=5-------------
pl1 <- Polygon(cbind(c(0,15,15,0,0), c(5,0,20,20,5)), hole=FALSE)
h1 <- Polygon(cbind(c(5,7,7,5,5), c(7,7,15,15,7)), hole=TRUE)
pl2 <- Polygon(cbind(c(15,20,20,30,30,15,15), c(10,10,0,0,20,20,10)), hole=FALSE)
sp <- SpatialPolygons(list(Polygons(list(pl1, h1), '0'), Polygons(list(pl2), '1')))
par(mar=c(0,0,0,0)); plot(sp) ### to visualize it
text(c(13, 17, 23), c(3, 12, 3), LETTERS[1:3], cex=3)

## ----hbond---------------------------------------------------------------
bound <- inla.sp2segment(sp)
mesh <- inla.mesh.2d(boundary=bound, max.edge=2)

## ----plothex,echo=F,results='hide',fig.width=7,fig.height=5,out.width='0.7\\textwidth'----
par(mar=c(0,0,0,0), bg=rgb(.7,.5,.5))
plot(sp, col=rgb(.3,.7,.9))
plot(mesh, add=TRUE, lwd=2)
text(c(13, 17, 23), c(3, 12, 3), LETTERS[1:3], cex=3)

