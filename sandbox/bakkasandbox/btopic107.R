## ----setup, include=FALSE------------------------------------------------
rm(list=ls())
knitr::opts_chunk$set(echo = TRUE)

## ---- warning=FALSE, message=FALSE---------------------------------------
library(INLA); library(sp); library(fields)

## Load data
#load(file = "https://haakonbakka.bitbucket.io/data/WebSiteData-Archipelago.RData")
# - if you have an internet connection
load(file = "data/WebSiteData-Archipelago.RData")
# - if you have saved the file locally

## What is loaded
# - poly.water is our study area
# - df is our dataframe to be analysed
# - dat is the orginial dataframe

## Load functions
#source('https://haakonbakka.bitbucket.io/functions-barrier-DT.R')
# - if you have an internet connection
source('functions-barrier-DT.R')
# - if you have saved the file locally

set.seed(2016)
set.inla.seed = 2016

## ------------------------------------------------------------------------
max.edge = 0.6
bound.outer = 4.6
mesh = inla.mesh.2d(boundary = poly.water,
                     loc=cbind(df$locx, df$locy),
                    max.edge = c(1,5)*max.edge,
                    cutoff = 0.06,
                    offset = c(max.edge, bound.outer))
plot(mesh, main="Our mesh", lwd=0.5); points(df$locx, df$locy, col="red")
mesh$n

## ------------------------------------------------------------------------
local.plot.field = function(field, mesh, xlim, ylim, ...){
  stopifnot(length(field) == mesh$n)
  # - error when using the wrong mesh
  if (missing(xlim)) xlim = poly.water@bbox[1, ] 
  if (missing(ylim)) ylim = poly.water@bbox[2, ]
  # - choose plotting region to be the same as the study area polygon
  proj = inla.mesh.projector(mesh, xlim = xlim, 
                             ylim = ylim, dims=c(300, 300))
  # - Can project from the mesh onto a 300x300 grid 
  #   for plots
  field.proj = inla.mesh.project(proj, field)
  # - Do the projection
  image.plot(list(x = proj$x, y=proj$y, z = field.proj), 
             xlim = xlim, ylim = ylim, ...)  
}

## ------------------------------------------------------------------------
df$y = df$y.smelt

## ------------------------------------------------------------------------
M = list()
# - the list of all models
M[[1]] = list()
# - the first model
M[[1]]$shortname = "stationary-no-cov"
# - a short description, the name of the model
M[[2]] = list()
M[[2]]$shortname = "stationary-all-cov"
M[[3]] = list()
M[[3]]$shortname = "barrier-no-cov"
M[[4]] = list()
M[[4]]$shortname = "barrier-all-cov"

## ------------------------------------------------------------------------
A.i.s = inla.spde.make.A(mesh, loc=cbind(df$locx, df$locy))

## ---- warning=FALSE, message=FALSE---------------------------------------
 stk <- inla.stack(data=list(y=df$y, e=df$exposure), 
                    effects=list(s=1:mesh$n,
                                 data.frame(m=1, df[ ,5:11]), 
                                 # - m is the intercept
                                 iidx=1:nrow(df)),
                   A=list(A.i.s, 1, 1),
                    remove.unused = FALSE, tag='est')  

## ---- warning=FALSE, message=FALSE---------------------------------------
spde = inla.spde2.pcmatern(mesh, prior.range = c(6, .5), prior.sigma = c(3, 0.01))
# - We put the prior median at approximately diff(range(df$locy))
# - - this is roughly the extent of our study area
# - The prior probability of marginal standard deviation 3 or more is 0.01.

## ------------------------------------------------------------------------
hyper.iid = list(prec = list(prior = 'pc.prec', param = c(3, 0.01))) 
# - the param the same as prior.sigma above, with the same interpretation
M[[1]]$formula = y~ -1+m + f(s, model=spde) + f(iidx, model="iid", hyper=hyper.iid)
# - no covariates (except intercept m)
M[[2]]$formula = as.formula(paste( "y ~ -1 + ",paste(colnames(df)[5:11], collapse = " + ")))
M[[2]]$formula = update(M[[2]]$formula, .~. +m + f(s, model=spde) + f(iidx, model="iid", hyper=hyper.iid))

## ------------------------------------------------------------------------
print(M[[1]])
print(M[[2]])

## ------------------------------------------------------------------------
mesh = inla.mesh.add.posTri(mesh)
# - compute the triangle positions
normal = over(poly.water, SpatialPoints(mesh$posTri), returnList=T)
# - checking which mesh triangles are inside the normal area
normal = unlist(normal)
Omega = dt.Omega(list(normal, 1:mesh$t), mesh)
Omega.SP = dt.polygon.omega(mesh, Omega)

## ------------------------------------------------------------------------
spde.dt = dt.FEMmatrices(mesh, Omega)
lambdas = c(-log(0.01)/3, -log(0.5)*6)
# - the first is log(prob)/exceed, the second log(prob)*exceed
# - the second is exponential for inverse range, therefore multiplication!
DT = dt.rgeneric.environment(spde.dt, 
                             prior.parameters=lambdas, 
                             fixed.ranges = c(NA, max.edge))
filename.rgeneric.environment = 'temp.rgeneric.init.RData'
save('DT', file = filename.rgeneric.environment)
barrier.model = inla.rgeneric.define(dt.rgeneric.dt.model, 
      n=dim(spde.dt$I)[1], ntheta = DT$ntheta, 
      R.init=filename.rgeneric.environment, debug=F) 

## ------------------------------------------------------------------------
M[[3]]$formula = y~ -1+m + f(s, model=barrier.model) + f(iidx, model="iid", hyper=hyper.iid)
# - no covariates (except intercept m)
M[[4]]$formula = as.formula(paste( "y ~ -1 + ",paste(colnames(df)[5:11], collapse = " + ")))
M[[4]]$formula = update(M[[4]]$formula, .~. +m + f(s, model=barrier.model) + f(iidx, model="iid", hyper=hyper.iid))

## ------------------------------------------------------------------------
print(M[[3]])
print(M[[4]])

## ------------------------------------------------------------------------
M[[1]]$init = c(2.499,1.196,-0.668)
# - set these to NULL the first time you run the model
M[[2]]$init = c(1.18,0.321,-0.737)
M[[3]]$init = c(0.606,2.092,-0.589)
M[[4]]$init = c(0.008,1.235,-0.733)

## ---- warning=FALSE, message=FALSE---------------------------------------
for (i in 1:length(M)){
    print(paste("Running:  ", M[[i]]$shortname))
    M[[i]]$res = inla(M[[i]]$formula,
                      data=inla.stack.data(stk),
                      control.predictor=list(A=inla.stack.A(stk)),
                      family="poisson", E = e,
                      control.mode=list(restart=T, theta=M[[i]]$init))  
}

## ------------------------------------------------------------------------
for (i in 1:length(M)){
  print(paste(round(M[[i]]$res$internal.summary.hyperpar$mode, 3), collapse = ','))
}

## ------------------------------------------------------------------------
summary(M[[1]]$res)

## ------------------------------------------------------------------------
summary(M[[2]]$res)

## ------------------------------------------------------------------------
summary(M[[3]]$res)

## ------------------------------------------------------------------------
summary(M[[4]]$res)

## ------------------------------------------------------------------------
#M[[i]]$res$logfile

## ------------------------------------------------------------------------
for (i in c(1,3)) {
  field = M[[i]]$res$summary.random$s$mean + M[[i]]$res$summary.fixed['m', 'mean']
  local.plot.field(field, mesh, main=paste(M[[i]]$shortname), zlim=c(-10.5, 1))
  plot(Omega.SP[[2]], add=T, col="grey")
  points(df$locx, df$locy)
}

## ------------------------------------------------------------------------
for (i in c(2,4)) {
  field = M[[i]]$res$summary.random$s$mean + M[[i]]$res$summary.fixed['m', 'mean']
  local.plot.field(field, mesh, main=paste(M[[i]]$shortname), zlim=c(-9, -4.5))
  plot(Omega.SP[[2]], add=T, col="grey")
  points(df$locx, df$locy)
}

## ------------------------------------------------------------------------
## Set up dataframe with relevant results
res = M[[2]]$res$summary.fixed[ ,c(4,3,5)]
# - all covar, some quantiles
for(i in c(4)) res = rbind(res, M[[i]]$res$summary.fixed[ ,c(4,3,5)])

colnames(res) = c("E", "L", "U")
rownames(res)=NULL
n.covar = nrow(M[[2]]$res$summary.fixed)
res$model = factor(rep(c("MS", "MB"), each=n.covar, 
                   levels = c("MS", "MB")))
res$covar = factor(rep(rownames(M[[2]]$res$summary.fixed), 2))

if(require(ggplot2)) {
  ggplot(res, aes(x = model, y = E)) +
    facet_wrap(~covar, scales = "free_y") +
    geom_point(size = 3) +
    geom_errorbar(aes(ymax = U, ymin = L)) +
    xlab(NULL) + ylab(NULL) 
}

