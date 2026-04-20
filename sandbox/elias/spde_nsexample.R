
library(INLA)

inla.setOption(
    inla.mode='experimental')

data(PRprec, package='INLA')

dim(PRprec)
n <- nrow(PRprec)
PRprec[1:3, 1:7]

crs.ll <- CRS('+proj=longlat +ellps=WGS84')
crs.mollkm <- CRS('+proj=moll +units=km')
prec.avg <- SpatialPointsDataFrame(
    spTransform(
        SpatialPoints(
            as.matrix(PRprec[,1:2]),
            crs.ll), crs.mollkm),
    data=data.frame(
        elevation=PRprec[,3],
        prec=rowMeans(PRprec[, 4:368])))

head(prec.avg)
summary(prec.avg)

boundary <- inla.nonconvex.hull(
    coordinates(prec.avg), 20, 20)

mesh <- inla.mesh.2d(
    boundary=boundary,
    max.edge=c(20, 70),
    cutoff=10,
    offset=100,
    n=25)
mesh$n

## projector
A <- inla.spde.make.A(mesh, coordinates(prec.avg))

d.stack <- inla.stack(
    tag='est',
    data=list(y=prec.avg$prec), 
    effects=list(
        b0=rep(1, n), 
        s=1:ncol(A)), 
    A=list(1, A)
)
#################################################
#### consider covariates in the covariance
#################################################

### define the 'covariates' at the mesh nodes
llloc <- coordinates(spTransform(
    SpatialPoints(mesh$loc[,1:2], crs.mollkm), crs.ll))
B0 <- cbind(west=llloc[,1]<(-51.5),
            north=llloc[,2]>(-24.5))

lk0 <- log(8)/2
lt0 <- -log(4*pi)/2

### stationary model
nsspde0 <- inla.spde2.matern(
    mesh=mesh, alpha=2,
    B.kappa=cbind(lk0, 0, 1),
    B.tau  =cbind(lt0, 1, 0))

### one covariate for each
nsspde1 <- inla.spde2.matern(
    mesh=mesh, alpha=2,
    B.kappa=cbind(lk0, 0, 1,       0, -B0[,2]),
    B.tau  =cbind(lt0, 1, 0, -B0[,1],  B0[,2]))

### both covariates for each
nsspde2 <- inla.spde2.matern(
    mesh=mesh, alpha=2,
    B.kappa=cbind(lk0, 0, 1, 0,0, -B0),
    B.tau  =cbind(lt0, 1, 0, -B0,  B0))

fns0 <- y ~ 0 + b0 + f(s, model=nsspde0)
fns1 <- y ~ 0 + b0 + f(s, model=nsspde1)
fns2 <- y ~ 0 + b0 + f(s, model=nsspde2)

ctrc <- list(dic=TRUE, waic=TRUE, cpo=TRUE)

pprec <- list(prior='pc.prec', param=c(1, 0.5))

fitns0 <- inla(
    formula=fns0,
    data=inla.stack.data(d.stack), 
    control.predictor=list(
        A=inla.stack.A(d.stack)),
    control.family=list(hyper=list(theta=pprec)),
    control.compute=ctrc)

fitns1 <- inla(
    formula=fns1,
    data=inla.stack.data(d.stack), 
    control.predictor=list(
        A=inla.stack.A(d.stack)),
    control.family=list(hyper=list(theta=pprec)),
    control.compute=ctrc)

fitns2 <- inla(
    formula=fns2,
    data=inla.stack.data(d.stack), 
    control.predictor=list(
        A=inla.stack.A(d.stack)),
    control.family=list(hyper=list(theta=pprec)),
    control.compute=ctrc)

rbind(fitns0$cpu, fitns1$cpu, fitns2$cpu)

rbind(ns0=fitns0$summary.hy[1,],
      ns1=fitns1$summary.hy[1,],
      ns2=fitns2$summary.hy[1,])

fitns0$summary.hy[-1,]
fitns1$summary.hy[-1,]
fitns2$summary.hy[-1,]

i.ok <- which(!is.na(prec.avg$prec))

rbind(mns0=c(dic=fitns0$dic$dic, waic=fitns0$waic$waic,
            cpo=-sum(log(fitns0$cpo$cpo[i.ok]))),
      mns1=c(dic=fitns1$dic$dic, waic=fitns1$waic$waic,
            cpo=-sum(log(fitns1$cpo$cpo[i.ok]))),
      mns2=c(dic=fitns2$dic$dic, waic=fitns2$waic$waic,
             cpo=-sum(log(fitns2$cpo$cpo[i.ok]))))

