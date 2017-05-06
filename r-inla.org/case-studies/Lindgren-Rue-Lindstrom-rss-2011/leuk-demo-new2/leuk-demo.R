## Code for Lindgren, Rue, Lindstrom (2011), partially updated to newer R-INLA
## Note: Will not reproduce exactly the same results as in the paper.
require(rgl)
require(INLA)
require(lattice)
source("utils.R")

## Read data:
map = read.table("Leuk.map")
data(Leuk)
Leuk$id = 1:dim(Leuk)[1]

## Build triangular mesh:
loc = cbind(Leuk$xcoord, Leuk$ycoord)
bnd1 = inla.nonconvex.hull(loc, convex=0.05)
bnd2 = inla.nonconvex.hull(loc, convex=0.25)
mesh = (inla.mesh.2d(
                  ## Data locations to use as location seeds:
                  loc=cbind(Leuk$xcoord, Leuk$ycoord),
                  ## Encapsulate data region:
                  boundary=list(bnd1, bnd2),
                  ## Refined triangulation,
                  ## minimal angles >=26 degrees,
                  ## interior maximal edge lengths 0.05,
                  ## exterior maximal edge lengths 0.2,
                  ## don't add input points closer than 0.05:
                  min.angle=24,
                  max.edge=c(0.05, 0.2),
                  cutoff=0.005,
                  ## Set to >=0 for visual (no effect Windows):
                  plot.delay=0.5
                  ))
## Store the data-->vertex mapping:
## ( mesh$idx$loc[k] is the mesh vertex for data location nr. k )
Leuk$spatial = mesh$idx$loc

## Create the SPDE/GMRF model, (kappa^2-Delta)(tau x) = W:
spde = inla.spde2.matern(mesh, alpha=2)

## Build the GLM model:
formula = inla.surv(Leuk$time, Leuk$cens) ~ 1 + sex + age + wbc + tpi +
    ## Add the spatial effect model:
    f(spatial, model=spde)

## Run INLA:
r  = (inla(formula, family="weibullsurv",
           data = Leuk,
           ## Prior specification:
           control.family = list(param=c(0.05,0.1)),
           ## Reasonable starting point for the optimisation:
           control.mode = list(theta=c(-0.5,-2,2),restart=TRUE),
           ## We don't need the marginals:
           control.compute = list(return.marginals=FALSE),
           ## We don't need to overoptimise:
           control.inla=list(tolerance=1e-5),
           ## Verbose output:
           verbose=TRUE
           ))

## Extract the SPDE parameters:
tau = exp(r$summary.hyperpar[2,"mean"])
kappa = exp(r$summary.hyperpar[3,"mean"])

## Get the precision matrix:
Q = inla.spde2.precision(spde, theta=log(c(tau, kappa)))


## Get a random sample (not used here),
## and the index reordering,
## using an undocumented ad hoc function:
reo=inla.qreordering(Q, reordering="metis")
## Need to invert the indexing:
neworder = reo$reordering
neworder[neworder] = 1:length(neworder)
## Reorder the matrix:
Q.reordered = Q[ neworder,neworder ]


## Reference point for covariance/correlation comparisons:
ref.s = (which.min((mesh$loc[,1]-mean(range(mesh$loc[,1])))^2 +
                   (mesh$loc[,2]-mean(range(mesh$loc[,2])))^2))

## Calculate covariances (S) and correlations (SS):
S = solve(Q)
SS = diag(1/sqrt(diag(S))) %*% S %*% diag(1/sqrt(diag(S)))
D = as.matrix(dist(mesh$loc))

## Theoretical Matern correlations and covariances:
dd = (0:1000)/1000
SS.theory = (dd*kappa)*besselK(dd*kappa,1)
SS.theory[1] = 1
S.theory = SS.theory/(4*pi*kappa^2)/tau^2




###########################
## Prepare for plotting:

## Calculate mapping between triangulation vertices and grid points:
## Resolution for gridded output was dims=c(100,100) in the paper.
proj = inla.mesh.projector(mesh, dims=c(200,200))

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
plotdata = inla.mesh.project(proj, r$summary.random$spatial[,"mean"])
## Plot PM contours:
dev.new()
bbb = (levelplot(row.values=proj$x, column.values=proj$y, x=plotdata,
                 mm=mm, panel=levelplotmap,
                 col.regions=my.palette,
                 xlim=range(proj$x), ylim=range(proj$y), aspect="iso",
                 contour=TRUE, cuts=11, labels=FALSE, pretty=TRUE,
                 xlab="Easting",ylab="Northing"))
print(bbb)


## Map resulting posterior standard deviation field to a grid:
plotdata = inla.mesh.project(proj, r$summary.random$spatial[,"sd"])
## Plot std.dev. contours:
dev.new()
bbb = (levelplot(row.values=proj$x, column.values=proj$y, x=plotdata,
                 mm=mm, panel=levelplotmap,
                 col.regions=my.palette,
                 xlim=range(proj$x), ylim=range(proj$y), aspect="iso",
                 contour=TRUE, cuts=11, labels=FALSE, pretty=TRUE,
                 xlab="Easting",ylab="Northing"))
print(bbb)


## Plot data points+map:
rgl.win1 = open3d()
plot.rgl.map(map, color = "black", lwd=2, zoffset=0.1)
plot.fake.points(cbind(Leuk$xcoord, Leuk$ycoord,0.1), radius=0.005,
                 color="black",specular="black")
plot(mesh, rgl=TRUE, add=TRUE,
     lwd=1, draw.vertices=FALSE, draw.edges=FALSE,
     col="white");
view3d(0,0,fov=0,zoom=0.8)

## Plot triangulation+map
rgl.win2 = open3d()
plot.rgl.map(map, color = "black", lwd=1, zoffset=0.1)
plot(mesh, rgl=TRUE, add=TRUE,
     draw.vertices=FALSE, draw.edges=TRUE,
     edge.color=rgb(0.6,0.6,0.6), lwd=1, col="white")
view3d(0,0,fov=0,zoom=0.8)
