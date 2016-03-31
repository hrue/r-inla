## ----settings, include=FALSE---------------------------------------------
library(knitr);  library(INLA);  library(fields) 
opts_chunk$set(message=FALSE, warning=FALSE, tidy=FALSE, fig.align='center')
knit_hooks$set(small.mar = function(before, options, envir) {
    if (before) par(mar = c(0.1, 0.1, .1, .1))  # smaller margin on top and right
})
knit_hooks$set(mar3311 = function(before, options, envir) {
    if (before) par(mar = c(3, 3, 1, 1), mgp=c(2, 1, 0))  # smaller margin on top and right
})
set.seed(1) 

## ----locations-----------------------------------------------------------
n = 200  ## number of location points
coo = matrix(runif(2*n), n) ## location points
k <- 10;   s2s <- 0.7 ## RF parameters
R <- s2s*exp(-k*as.matrix(dist(coo))) ## covariance matrix

## ----field simulate------------------------------------------------------
s <- drop(rnorm(n)%*%chol(R)) ## one RF realization

## ----covariate and outcome-----------------------------------------------
x <- runif(n)  ## covariate
beta <- 1:2  ## regression coefficients
lin.pred <- beta[1] + beta[2]*x + s ## linear predictor
s2e <- 0.3 ## error variance (nugget)
y <- lin.pred + rnorm(n, 0, sqrt(s2e))  ## the outcome

## ----mesh, fig.width=4.5, fig.height=4.6, small.mar=TRUE-----------------
r0.1 <- sqrt(0.5 * 8)/k ## distance with correlation around 0.139
mesh <- inla.mesh.2d( ## 2D mesh creator
  loc=coo, ## provided locations 
  max.edge=c(r0.1/3, r0.1), ## maximum edge length (inner, outer): mandatory 
  offset=c(r0.1/3, r0.1*2), ## outer extension
  cutoff=r0.1/10) ## good to have >0
par(mar=c(0,0,1,0))
plot(mesh, asp=1) ## plot the mesh
points(coo, col='red') ## add the points 

## ---- projector, fig.width=8, fig.height=4.5, small.mar=TRUE-------------
image(A <- inla.spde.make.A( ## projector creator
    mesh=mesh, ## provide the mesh
    loc=coo) ### locations where to project the field
    ) ## an 'n' by 'm' projector matrix

## ----SPDE model----------------------------------------------------------
spde <- inla.spde2.matern( ## precision components creator
    mesh=mesh, ## mesh supplied
    alpha=1.5) ## smoothness parameter

## ----data stack----------------------------------------------------------
stk.e <- inla.stack( ## stack creator
  data=list(y=y),  ## response
  effects=list(## two elements:
    data.frame(b0=1, x=x), ## regressor part
    s=1:spde$n.spde),  ## RF index
  A=list(## projector list of each effect
    1, ## for the covariates
    A), ## for the RF
  tag='est') ## tag

## ----fitting-------------------------------------------------------------
formula <- y ~ 0 + b0 + x + ## fixed part
  f(s, model=spde) ## RF term
res <- inla( ## main function in INLA package
  formula, ## model formula
  data=inla.stack.data(stk.e), ## dataset
  control.predictor=list( ## inform projector needed in SPDE models
    A = inla.stack.A(stk.e))) ## projector from the stack data

## ----fixed---------------------------------------------------------------
round(res$summary.fixed, 4) 

## ----nugget, fig.width=3.5, fig.height=3, mar3311=TRUE-------------------
m.prec <- res$marginals.hyperpar$'Precision for the Gaussian observations' ## the marginal
post.s2e <- inla.tmarginal(## function to compute a tranformation 
  function(x) 1/x, ## inverse transformation
  m.prec) ## marginal to be applied
### visualize it
plot(post.s2e, type='l', ylab='Density', 
     xlab=expression(sigma[e]^2))
abline(v=s2e, col=2) ## add value used to generate the data

## ----rf------------------------------------------------------------------
rf <- inla.spde.result( ## function to compute the 'interpretable' parameters
    inla=res, ## the inla() output
    name='s', ## name of RF index set
    spde=spde, ## SPDE model object
    do.transf=TRUE) ## to user scale

## ----parameters, fig.width=10, fig.height=3------------------------------
par(mfrow=c(1,3), mar=c(3,3,0.3,0.3), mgp=c(2,0.5,0))
plot(rf$marginals.var[[1]], ty='l', 
     xlab=expression(sigma[s]^2), yla='Density')
abline(v=s2s, col=2) ## add the true value
plot(rf$marginals.kap[[1]], type='l',
     xlab=expression(kappa), ylab='Density')
abline(v=k, col=2) ## add the true value
plot(rf$marginals.range[[1]], type='l', 
     xlab='range nominal', ylab='Density')
abline(v=sqrt(0.5 * 8)/k, col=2) ## add the 'true' value

## ----project-------------------------------------------------------------
gproj <- inla.mesh.projector( ## projector builder
  mesh, ## mesh used to define the model
  xlim=0:1, ylim=0:1, ## limits where to create the grid
  dims=c(300,300)) ## grid dimension
## project the mean and the SD
g.mean <- inla.mesh.project(gproj, res$summary.random$s$mean)
g.sd <- inla.mesh.project(gproj, res$summary.random$s$sd)

## ----visualize, fig.width=9, fig.height=5--------------------------------
par(mfrow=c(1,2), mar=c(0,0,1,0))
require(fields)
image.plot(g.mean, asp=1, main='RF posterior mean', axes=FALSE, horizontal=TRUE)
image.plot(g.sd, asp=1, main='RF posterior SD', axes=FALSE, horizontal=TRUE)

## ----prediction scenario-------------------------------------------------
tcoo <- rbind(c(0.3,0.3), c(0.5,0.5), c(0.7,0.7))
dim(Ap <- inla.spde.make.A(mesh=mesh, loc=tcoo)) 
x0 <- c(0.5, 0.5, 0.5)

## ----prediction stack----------------------------------------------------
stk.pred <- inla.stack(
  tag='pred', ## will be used to collect the posterior marginals
  data=list(y=NA), ## response set as NA
  effects=list(
    data.frame(x=x0, b0=1), ## covariate scenario
    s=1:spde$n.spde), ## same as before
  A=list(1, Ap)) ## covariate and target locations field projectors

## ----refitting-----------------------------------------------------------
stk.full <- inla.stack(stk.e, stk.pred) ## join both data
p.res <- inla(
  formula, data=inla.stack.data(stk.full), ## using the full data 
  control.predictor=list(compute=TRUE, ## compute the predictor
                         A=inla.stack.A(stk.full)), ## from full data
  control.mode=list(theta=res$mode$theta)) ## use the mode previously found

## ----prediction index----------------------------------------------------
pred.ind <- inla.stack.index( ## stack index extractor function
  stk.full, ## the data stack to be considered
  tag='pred' ## which part of the data to look at
  )$data ## which elements to collect
ypost <- p.res$marginals.fitted.values[pred.ind]

## ----predicted, fig.width=10, fig.height=3-------------------------------
xyl <- apply(Reduce('rbind', ypost), 2, range)
par(mfrow=c(1,3), mar=c(3,3,2,1), mgp=c(2,1,0))
for (j in 1:3) 
  plot(ypost[[j]], type='l', xlim=xyl[,1], ylim=xyl[,2], 
       xlab='y', ylab='density', main=paste0('y', j))

## ----marginals-----------------------------------------------------------
apropos('marginal')

## ----marginals playing---------------------------------------------------
inla.qmarginal(c(0.15, 0.7), ypost[[1]]) ## quantiles
inla.emarginal(function(x) x^2, ypost[[1]]) - ## E(y^2) -
  inla.emarginal(function(x) x, ypost[[1]])^2 ## E(y)^2 to compute the variance
inla.pmarginal(inla.qmarginal(0.3, ypost[[1]]), ypost[[1]]) 
inla.zmarginal(ypost[[1]]) ## posterior summary

