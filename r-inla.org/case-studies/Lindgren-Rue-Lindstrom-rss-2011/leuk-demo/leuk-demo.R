library(rgl)
library(INLA)
source("utils.R")

## Read data:
map = read.table("Leuk.map")
data(Leuk)
Leuk$id = 1:dim(Leuk)[1]

## Build triangular mesh:
mesh = (inla.fmesher.mesh(
                          ## Data locations:
                          cbind(Leuk$xcoord, Leuk$ycoord),
                          ## Where to put the mesh/graph files:
                          prefix="./leuk-graph/mesh.",
                          ## Set to >=0 for visual (not on Windows):
                          x11.delay=NULL,
                          ## Encapsulate data region with a relative margin:
                          cet=c(8,-0.15),
                          ## Refined triangulation,
                          ## minimal angles >=26 degrees,
                          ## exterior maximal edge lengths 0.2,
                          ## interior maximal edge lengths 0.08:
                          rcdt=c(26,0.2,0.08)))
## Store the data-->vertex mapping:
## ( locations.idx[k] is the mesh vertex for data point #k )
Leuk$spatial = mesh$locations.idx
## Create the SPDE/GMRF precision building blocks, 2nd order FEM:
spde = inla.create.spde(prefix=mesh$prefix, fem=2)

## Build the GLM model:
formula = inla.surv(Leuk$time, Leuk$cens) ~ 1 + sex + age + wbc + tpi +
    ## Add the spatial effect model:
    f(spatial, model="spde",
      spde.prefix = spde$prefix, n = spde$n,
      ## Wide prior specification for the log-tau and log-kappa^2 parameters:
      param = c(-2,0.02,4,0.02,NA,NA,NA,NA))

## Run INLA:
r  = (inla(formula, family="weibull",
           data = Leuk,
           ## Prior specification:
           control.data = list(param=c(0.05,0.1)),
           ## Reasonable starting point for the optimisation:
           control.mode = list(theta=c(0,-2,5),restart=TRUE),
           ## We don't need the marginals:
           control.compute = list(return.marginals=FALSE),
           ## We don't need to overoptimise:
           control.inla=list(tolerance=1e-5),
           ## Verbose output, binary I/O, keep the resulting data files:
           verbose=TRUE, keep=TRUE, inla.arg = "-v -b"
           ))

## Extract the SPDE parameters:
tau = exp(r$summary.hyperpar[2,"mean"])
kappa = exp(r$summary.hyperpar[3,"mean"]/2)

## Extract the precision building blocks:
c0 = inla.read.fmesher.file(paste(spde$prefix, "c0", sep=""))
g1 = inla.read.fmesher.file(paste(spde$prefix, "g1", sep=""))
g2 = inla.read.fmesher.file(paste(spde$prefix, "g2", sep=""))

## Build the precision matrix:
Q = tau^2*(kappa^4*c0+2*g1*kappa^2+g2)


## Get a random sample (not used here),
## and the index reordering,
## using an undocumented ad hoc function:
finn=inla.finn(Q,reordering="metis")
## Need to invert the indexing:
neworder = finn$reordering
neworder[finn$reordering] = 1:length(finn$reordering)
## Reorder the matrix:
Q.reordered = Q[ neworder,neworder ]


## Reference point for covariance/correlation comparisons:
ref.s = (which.min((mesh$mesh$s[,1]-mean(range(mesh$mesh$s[,1])))^2 +
                   (mesh$mesh$s[,2]-mean(range(mesh$mesh$s[,2])))^2))

## Calculate covariances (S) and correlations (SS):
S = solve(Q)
SS = diag(1/sqrt(diag(S))) %*% S %*% diag(1/sqrt(diag(S)))
D = as.matrix(dist(mesh$mesh$s))

## Theoretical Matern correlations and covariances:
dd = (0:1000)/1000
SS.theory = (dd*kappa)*besselK(dd*kappa,1)
SS.theory[1] = 1
S.theory = SS.theory/(4*pi*kappa^2)/tau^2




###########################
## Prepare for plotting:

## Resolution for gridded output: (n=100 in the paper)
mapgrid.n = 200
mapgrid.x = (((0:mapgrid.n)/mapgrid.n)*diff(range(mesh$mesh$s[,1]))+
             min(mesh$mesh$s[,1]))
mapgrid.y = (((0:mapgrid.n)/mapgrid.n)*diff(range(mesh$mesh$s[,2]))+
             min(mesh$mesh$s[,2]))
mapgrid.xy = (cbind(rep(mapgrid.x,length(mapgrid.y)),
                    rep(mapgrid.y,1,each=length(mapgrid.x)),
                    0))
## Calculate mapping between triangulation vertices and grid points:
mapgrid = points2mesh.calc(prefix=mesh$prefix, mapgrid.xy)


## Construct greyscale palette function:
my.grey.palette = function (n,...) { return (grey.colors(n,0.05,0.95,...))}
## Use it:
my.palette = my.grey.palette

## Construct map data appropriate for easy plotting:
mm = calc.map(map)


#####################
## Plot results:

## Compare correlations:
dev.new()
plot(D[ref.s,],SS[ref.s,],type="p",pch=20,
     xlab="Distance",ylab="Correlation",
     xlim=c(0.0,1.0),ylim=c(-0.005,1.005))
lines(dd,SS.theory,type="l",
      col=rgb(0.5,0.5,0.5),lwd=2)


## Don't plot the full precision pattern; reduce the dimension first:
Q.coarse = sparse.pattern.coarsen(Q.reordered,4)
## Plot the reordered precision pattern:
dev.new()
image(as.matrix(Q.coarse)>0,
      col=grey.colors(2,start=1,end=0),
      ylim=c(1,0), axes=FALSE, xlim=c(0,1))
box()



## Map resulting posterior mean field to a grid:
plotdata = as.vector(mapgrid$A %*% r$summary.random$spatial[,"mean"])
plotdata[!mapgrid$ok] = NA
plotdata = (matrix(plotdata,length(mapgrid.x),length(mapgrid.y)))
## Plot PM contours:
dev.new()
bbb = (levelplot(row.values=mapgrid.x, column.values=mapgrid.y, x=plotdata,
                 mm=mm, panel=levelplotmap,
                 col.regions=my.palette,
                 xlim=range(mapgrid.x), ylim=range(mapgrid.y), aspect="iso",
                 contour=TRUE, cuts=11, labels=FALSE, pretty=TRUE,
                 xlab="Easting",ylab="Northing"))
print(bbb)


## Map resulting posterior standard deviation field to a grid:
plotdata = as.vector(mapgrid$A %*% r$summary.random$spatial[,"sd"])
plotdata[!mapgrid$ok] = NA
plotdata = (matrix(plotdata,length(mapgrid.x),length(mapgrid.y)))
## Plot std.dev. contours:
dev.new()
bbb = (levelplot(row.values=mapgrid.x, column.values=mapgrid.y, x=plotdata,
                 mm=mm, panel=levelplotmap,
                 col.regions=my.palette,
                 xlim=range(mapgrid.x), ylim=range(mapgrid.y), aspect="iso",
                 contour=TRUE, cuts=11, labels=FALSE, pretty=TRUE,
                 xlab="Easting",ylab="Northing"))
print(bbb)


## Plot data points+map:
rgl.win1 = open3d()
plot.rgl.map(map, color = "black", lwd=2, zoffset=0.1)
plot.fake.points(cbind(Leuk$xcoord, Leuk$ycoord,0.1), radius=0.005,
                 color="black",specular="black")
plot.inla.trimesh(mesh$mesh$tv, mesh$mesh$s,
                  rgb(1,1,1), lwd=1,
                  draw.vertices=FALSE, draw.edges=FALSE,
                  color="white");
view3d(0,0,fov=0,zoom=0.8)

## Plot triangulation+map
rgl.win2 = open3d()
plot.inla.trimesh(mesh$mesh$tv, mesh$mesh$s,
                  rgb(1,1,1), draw.vertices=FALSE, draw.edges=TRUE,
                  edge.color=rgb(0.6,0.6,0.6), lwd=1, color="white")
plot.rgl.map(map, color = "black", lwd=1, zoffset=0.1)
view3d(0,0,fov=0,zoom=0.8)
