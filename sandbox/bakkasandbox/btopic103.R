## ----setup, include=FALSE------------------------------------------------
rm(list=ls())
knitr::opts_chunk$set(echo = TRUE)

## ---- warning=FALSE, message=FALSE---------------------------------------
library(INLA)
#source('https://haakonbakka.bitbucket.io/functions-barrier-DT.R')
# - if you have an internet connection
source('functions-barrier-DT.R')
# - if you have saved the file locally
set.seed(2016)
set.inla.seed = 2016

## ------------------------------------------------------------------------
smalldist = 0.2
# - the width of the opening in the barrier
width = 0.5
# - The width/thickness of the barrier
p2 = bakka.square.polygon(xlim=c(-1, 5-smalldist/2), 
                          ylim=5+width*c(-.5, .5), ret.SP = T)
p3 = bakka.square.polygon(xlim=c(5+smalldist/2, 11), 
                          ylim=5+width*c(-.5, .5), ret.SP = T)
p = SpatialPolygons(c(p2@polygons, p3@polygons))
# - The total barrier area polygon
plot(p, main="Barrier area polygon")

## ------------------------------------------------------------------------
max.edge.length = 0.4
# - The coarseness of the finite element approximation
# - Corresponds to grid-square width in discretisations
# - - Except that finite element approximations are better
# - Should be compared to size of study area
# - Should be less than a fourth of the estimated (posterior) 
#   spatial range
# - Computational time up to *8 when you halve this value


## ------------------------------------------------------------------------
loc1 = matrix(c(0,0, 10,0, 0,10, 10,10), 4, 2, byrow = T)
# - This defines the extent of the mesh
# - In an application, if you want the mesh to depend on your 
#   data locations, you may use those locations here
seg = inla.sp2segment(p)
# - Transforms a SpatialPolygon to an "inla polygon"
mesh = inla.mesh.2d(loc=loc1, interior = seg, 
                    max.e = max.edge.length, offset=1)
# - The INLA mesh constructor, used for any INLA-SPDE model
mesh = inla.mesh.add.posTri(mesh)
# - Add on mesh$posTri
# - - contains the positions of the triangles


## ------------------------------------------------------------------------
barrier = over(p, SpatialPoints(mesh$posTri), returnList=T)
# - checking which mesh triangles are inside the barrier area
barrier = unlist(barrier)
Omega = dt.Omega(list(barrier, 1:mesh$t), mesh)
Omega.SP = dt.polygon.omega(mesh, Omega)
# - creates polygons for the different areas: Barrier area 
#   and Normal area

## ------------------------------------------------------------------------
plot(mesh, main="Mesh and Omega")
plot(Omega.SP[[1]], add=T, col='grey')
plot(Omega.SP[[2]], add=T, col='lightblue')
plot(mesh, add=T)
points(loc1)

## ------------------------------------------------------------------------
local.plot.field = function(field, ...){
  xlim = c(2, 8); ylim = xlim;
  proj = inla.mesh.projector(mesh, xlim = xlim, 
                             ylim = ylim, dims=c(300, 300))
  # - Can project from the mesh onto a 300x300 grid 
  #   for plots
  field.proj = inla.mesh.project(proj, field)
  # - Do the projection
  image.plot(list(x = proj$x, y=proj$y, z = field.proj), 
             xlim = xlim, ylim = ylim, ...)  
  # - Use image.plot to get nice colors and legend
}
print(mesh$n)
# - This is the appropriate length of the field variable

## ------------------------------------------------------------------------
spde.dt = dt.FEMmatrices(mesh, Omega)
# - Set up the matrices for solving the SPDE
# - May be time consuming

## ------------------------------------------------------------------------
ranges = c(0.3, 3)
# - the first range is for the barrier area
# - - it is not sensitive to the exact value here, 
#     just make it "small"
# - the second range is for the normal area
Q = dt.precision(spde.dt, ranges)
# - the precision matrix for fixed ranges

## ------------------------------------------------------------------------
u = inla.qsample(n=1, Q=Q, seed = set.inla.seed)
u = u[ ,1]
local.plot.field(u, main="The true (simulated) spatial field")
plot(Omega.SP[[1]], add=T, col='grey')
# - Overlay the barrier

## ------------------------------------------------------------------------
num.try = 500 
# - try to sample this number of data locations
loc.try = matrix(runif(num.try*2, min=2, max=8), 
                         num.try, 2)
# - locations sampled inside the barrier will be removed 
#   in a few lines
temp = SpatialPoints(loc.try)
loc.ok = !is.na(over(temp, Omega.SP[[2]]))
# - only allow locations that are not inside the Barrier area
loc.data = loc.try[loc.ok, ]
A.data = inla.spde.make.A(mesh, loc.data)
# - the projector matrix required for any spatial model
# - this matrix can transform the field-defined-on-the-mesh 
#   to the field-defined-on-the-data-locations
c(dim(A.data), mesh$n, nrow(loc.data))
# - shows that the dimensions are correct
u.data = A.data %*% u
# - project the field from the finite element  
#   representation to the data locations
df = data.frame(loc.data)
# - df is the dataframe used for modeling
names(df) = c('locx', 'locy')
sigma.u = 1
# - size of the random effect
# - feel free to change this value
sigma.epsilon = 0.2
# - size of the iid noise in the Gaussian likelihood
# - feel free to change this value
df$y = drop(sigma.u*u.data + sigma.epsilon*rnorm(nrow(df)))
# - sample observations with gaussian noise

## ------------------------------------------------------------------------
summary(df)

## ------------------------------------------------------------------------
stk <- inla.stack(data=list(y=df$y), A=list(A.data, 1),
                  effects=list(s=1:mesh$n, 
                               intercept=rep(1, nrow(df))), 
                  remove.unused = FALSE, tag='est')
# - this is the common stack used in INLA SPDE models
# - see the SPDE-tutorial
# - - http://www.r-inla.org/examples/tutorials/spde-tutorial

## ------------------------------------------------------------------------
model.stat = inla.spde2.matern(mesh)
# - Set up the model component for the spatial SPDE model: 
#   Stationary Matern model
# - I assume you are somewhat familiar with this model

formula <- y ~ 0 + intercept + f(s, model=model.stat)
# - Remove the default intercept
# - - Having it in the stack instead improves the numerical 
#     accuracy of the INLA algorithm
# - Fixed effects + random effects

res.stationary <- inla(formula, data=inla.stack.data(stk),
            control.predictor=list(A = inla.stack.A(stk)),
            family = 'gaussian',
            control.family = list(hyper = list(prec = list(
              prior = "pc.prec", fixed = FALSE, 
              param = c(0.2,0.5)))),
            control.mode=list(restart=T, theta=c(4,-1.7,0.25)))

## ------------------------------------------------------------------------
summary(res.stationary)

## ------------------------------------------------------------------------
local.plot.field(res.stationary$summary.random$s$mean,
          main="Spatial estimate with the stationary model")
# - plot the posterior spatial marginal means
# - we call this the spatial estimate, or the smoothed data
plot(Omega.SP[[1]], add=T, col='grey')
# - Posterior spatial estimate using the stationary model

## ------------------------------------------------------------------------
DT = dt.rgeneric.environment(spde.dt, prior.parameters=c(1, 1), 
                             fixed.ranges = c(0.05, NA))
# - This function creates the environment/variables needed 
#   to run the internal 
#   rgeneric engine inside the C-code in R-INLA
# - We fix the barrier range to a different value than we 
#   used for simulations
# - - Why? It does not matter, as long as it is 'small' 
#     the models are very
#     similar
# - - This shows that you do not need to know the 
#     true 'barrier range'!
# - The prior parameters are the lambdas in the exponential 
#   priors for standard 
#   deviation and inverse-range

filename.rgeneric.environment = 'temp.rgeneric.init.RData'
save('DT', file = filename.rgeneric.environment)
# - save the new environment file
# - this is not an actual R environment, just a 
#   list-of-lists that is loaded 
#   at rgeneric startup

## ------------------------------------------------------------------------
barrier.model = inla.rgeneric.define(dt.rgeneric.dt.model, 
      n=dim(spde.dt$I)[1], ntheta = DT$ntheta, 
      R.init=filename.rgeneric.environment, debug=F) 
# - This contains all the needed functions for the spatial 
#   model component
# - - Most important is the function computing the 
#     precision matrix Q from the hyperparameters theta
# - These functions are called by the C-code in R-INLA 
#   when R-INLA is performing inference
# - - The inla(...) call spawns a C algorithm running 
#     inference, which again
#     spawns an R algorithm running rgeneric

## ------------------------------------------------------------------------
formula2 <- y ~ 0 + intercept + f(s, model=barrier.model)
# - The spatial model component is different from before
# - The rest of the model setup is the same! 
#   (as in the stationary case)
# - - e.g. the inla(...) call below is the same, 
#     only this formula is different

## ------------------------------------------------------------------------
res.barrier = inla(formula2, data=inla.stack.data(stk),
       control.predictor=list(A = inla.stack.A(stk)),
       family = 'gaussian',
       control.family = list(hyper = list(prec = list(
             prior = "pc.prec", fixed = FALSE, 
             param = c(0.2,0.5)))),
       control.mode=list(restart=T, theta=c(3.2, 0.4, 1.6)))

## ------------------------------------------------------------------------
summary(res.barrier)

## ------------------------------------------------------------------------
local.plot.field(res.barrier$summary.random$s$mean, 
                 main="Spatial posterior for Barrier model")
# - plot the posterior spatial marginal means
# - we call this the spatial (smoothing) estimate
plot(Omega.SP[[1]], add=T, col='grey')
# - Posterior spatial estimate using the Barrier model

