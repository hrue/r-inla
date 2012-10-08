remove(list=ls())
require(geoR)
require(INLA)
require(gstat)
require(fields)

include.elevation = T #include/exclude elevation 


###########################################################################
#DATA PREPARATION
###########################################################################
#-- Load the data --#
data(SIC) #use sic.all

#-- Center the coordinates (with respect to the mean of the border coords) --#                                  
sic.all$coords[,1] = (sic.all$coords[,1]-apply(sic.borders,2,mean)[1])
sic.all$coords[,2] = (sic.all$coords[,2]-apply(sic.borders,2,mean)[2])
sic.borders = apply(sic.borders,2,scale,scale=F)

#-- 100 Stations selected for estimation (as in sic.100) --#
index.for.estimation = which(as.numeric(rownames(sic.all$coords)) %in% as.numeric(rownames(sic.100$coords)))

#-- Data for estimation (100 locations) --#
#-- Take the square root for rain measurements --#
est.coord = sic.all$coords[index.for.estimation,]
est.data =  sqrt(sic.all$data[index.for.estimation])
est.elevation = sic.all$altitude[index.for.estimation]/1000

#-- Data for validation (367 locations) --#
#-- Take the square root for rain measurements --#
val.coord = sic.all$coords[-index.for.estimation,]
val.data = sqrt(sic.all$data[-index.for.estimation])
val.elevation = sic.all$altitude[-index.for.estimation]/1000

#-- Rainfall thematic Map --#
color.map = gray(c(1, .8, .5, 0))
sic.cutoff = quantile(sqrt(sic.all$data))

sic.factor.est = cut(est.data,breaks=sic.cutoff,include.lowest=TRUE,right=TRUE)
q4col.est = color.map[as.numeric(sic.factor.est)]
sic.factor.val = cut(val.data,breaks=sic.cutoff,include.lowest=TRUE,right=TRUE)
q4col.val = color.map[as.numeric(sic.factor.val)]

par(mar=c(0,0,0,0)+0.2)
plot(sic.borders, type="l",axes=F, xlab="",ylab="")
points(est.coord[,"V2"], est.coord[,"V3"], bg=q4col.est, pch=21, cex=1.2, lwd=2)
points(val.coord[,"V2"], val.coord[,"V3"], bg=q4col.val, pch=24, cex=1.2)
legend(x=-190, y=74, pch=rep(21,4), pt.cex=rep(1.2,4),
       pt.bg=rev(color.map), legend=rev(levels(sic.factor.val)), bty="n",cex=1.2)
legend(x=-200, y=74, pch=rep(24,4), pt.cex=rep(1.2,4),
       pt.bg=rev(color.map), legend=c("","","",""), bty="n",cex=1.2)

#-- Elevation surface for prediction on the grid --#
data(sic97) #we need demstd

seq.x.grid = seq(from=demstd@coords[1,1],to=demstd@coords[2,1],
                 length=demstd@grid@cells.dim[1])/1000 #376 (m-->km)
seq.y.grid = seq(from=demstd@coords[1,2],to=demstd@coords[2,2],
                 length=demstd@grid@cells.dim[2])/1000 #253 (m-->km)
pred.grid = expand.grid(seq.x.grid,seq.y.grid)

elevation.grid = matrix(as.matrix(demstd@data), nrow=demstd@grid@cells.dim[1],
                        ncol=demstd@grid@cells.dim[2], byrow=F)
elevation.grid = elevation.grid[,253:1]/1000 #elevation

#-- Plot the elevation surface and monitoring stations--#
par(mar=c(4,4,0.5,1.5))
image.plot(seq.x.grid,seq.y.grid,elevation.grid, col=grey(50:0/50),
          xlab="W-E (km)", ylab="N-S (km)")
lines(sic.borders,lwd=2)


###########################################################################
#SPDE/INLA 
###########################################################################
#-- Create the mesh --#
mesh =
  inla.mesh.create.helper(points=est.coord,
                          points.domain=sic.borders,
                          offset=c(5, 20),
                          max.edge=c(40,100), min.angle=c(21,21))
mesh$n

#-- Plot the triangulation --#
par(mar=c(0,2,2,0))
plot(mesh)
lines(sic.borders,lwd=3)
points(est.coord,pch=19, col=1, cex=1.2)

#-- Construct the SPDE object --#
spde = inla.spde2.matern(mesh=mesh)

#-- Create A matrices (mapping between the mesh nodes and station locations) --#
A.est =
  inla.spde.make.A(mesh, loc=est.coord)
A.val =
  inla.spde.make.A(mesh, loc=val.coord)
A.pred =
  inla.spde.make.A(mesh)

#-- Create index matrix --#
field.indices =
  inla.spde.make.index("field", n.mesh=mesh$n)

#-- Create data stacks --#
if(include.elevation){
  stack.est =
    inla.stack(data=list(rain=est.data),
               A=list(A.est, 1),
               effects=
                 list(c(field.indices,
                        list(Intercept=1)),
                      list(Elevation=est.elevation)),
               tag="est")

  stack.val =
    inla.stack(data=list(rain=NA),
               A=list(A.val,1),
               effects=
                 list(c(field.indices,
                        list(Intercept=1)),
                      list(Elevation=val.elevation)),
               tag="val")

  stack.pred =
    inla.stack(data=list(rain=NA),
               A=list(A.pred),
               effects=
                 list(c(field.indices,
                        list(Intercept=1))),
               tag="pred")

  stack = inla.stack(stack.est, stack.val, stack.pred)
  

  #-- Define the formula --#
  formula <- rain ~ -1 + Intercept + Elevation + f(field, model=spde)

}else{
  #no elevation
  stack.est =
    inla.stack(data=list(rain=est.data),
               A=list(A.est),
               effects=
                 list(c(field.indices,
                        list(Intercept=1))),
               tag="est")
  stack.val =
    inla.stack(data=list(rain=NA),
               A=list(A.val),
               effects=
                 list(c(field.indices,
                        list(Intercept=1))),
               tag="val")

  stack.pred =
    inla.stack(data=list(rain=NA),
               A=list(A.pred),
               effects=
                 list(c(field.indices,
                        list(Intercept=1))),
               tag="pred")

  stack = inla.stack(stack.est, stack.val, stack.pred)

  #-- Define the formula --#
  formula <- rain ~ -1 + Intercept +  f(field, model=spde)

}

#-- Call INLA and get results --#
mod =   inla(formula,
             data=inla.stack.data(stack, spde=spde),
             family="gaussian",
             control.predictor=list(A=inla.stack.A(stack), compute=TRUE),
             control.compute=list(cpo=TRUE, dic=TRUE),
             keep=FALSE, verbose=TRUE)

mod$dic$dic
print(summary(mod))

###########################################################################
#EXTRACT POSTERIOR SUMMARY STATISTICS
###########################################################################

#-- Extract results for fixed effects - covariate coeffs --#
beta = mod$summary.fixed[,"mean"]
beta_sd = mod$summary.fixed[,"sd"]

#-- Extract results for sigma2eps (1/precision)
sigma2e_marg =
  inla.tmarginal(function(x) 1/x,
                 mod$marginals.hyperpar$"Precision for the Gaussian observations")
sigma2e_m1 = inla.emarginal(function(x) x, sigma2e_marg)
sigma2e_m2 = inla.emarginal(function(x) x^2, sigma2e_marg)
sigma2e_stdev = sqrt(sigma2e_m2 - sigma2e_m1^2)
sigma2e_quantiles = inla.qmarginal(c(0.025, 0.5, 0.975), sigma2e_marg)
cat("-----Results for sigma2eps-----\n",
    c(sigma2e_m1, sigma2e_stdev, sigma2e_quantiles),
    "\n-----------------------------")

#-- Extract results for the covariance function parameters --#
mod.field = inla.spde2.result(mod, name="field", spde)

var.nom.marg = mod.field$marginals.variance.nominal[[1]]
var.nom.m1 = inla.emarginal(function(x) x, var.nom.marg)
var.nom.m2 = inla.emarginal(function(x) x^2, var.nom.marg)
var.nom.stdev = sqrt(var.nom.m2 - var.nom.m1^2)
var.nom.quantiles = inla.qmarginal(c(0.025, 0.5, 0.975), var.nom.marg)
cat("-----Results for sigma2omega-----\n",
    c(var.nom.m1, var.nom.stdev, var.nom.quantiles),
    "\n-----------------------------")

range.nom.marg = mod.field$marginals.range.nominal[[1]]
range.nom.m1 = inla.emarginal(function(x) x, range.nom.marg)
range.nom.m2 = inla.emarginal(function(x) x^2, range.nom.marg)
range.nom.stdev = sqrt(range.nom.m2 - range.nom.m1^2)
range.nom.quantiles = inla.qmarginal(c(0.025, 0.5, 0.975), range.nom.marg)
cat("-----Results for the range-----\n",
    c(range.nom.m1, range.nom.stdev, range.nom.quantiles),
    "\n-----------------------------")

###########################################################################
#PREDICTION IN THE VALIDATION STATIONS
###########################################################################
#-- Extract results for the linear predictor (lp) in the validation stations --#
index.val = inla.stack.index(stack,"val")$data
lp.mean.val = mod$summary.linear.predictor[index.val,"mean"]
lp.sd.val = mod$summary.linear.predictor[index.val,"sd"]

#-- Compare observed vs predicted values and compute residuals --#
par(mar=c(4,4,2,2))
plot(val.data, lp.mean.val, xlab="Observed", ylab="Predicted") #sqrt scale

res = val.data - lp.mean.val
cor(val.data, lp.mean.val)
sqrt(mean(res^2)) #RMSE

###########################################################################
#PREDICTION IN THE GRID STATIONS
###########################################################################
#-- Create mapping between mesh and grid --#
proj_grid =
  inla.mesh.projector(mesh,
                      xlim=range(pred.grid[,1]),
                      ylim=range(pred.grid[,2]),
                      dims=c(376,253))


index.pred = inla.stack.index(stack,"pred")$data
lp.mean.pred = mod$summary.linear.predictor[index.pred, "mean"]
lp.sd.pred =  mod$summary.linear.predictor[index.pred, "sd"]
  
#-- Project the linear predictor on the grid --#
lp.mean.grid = inla.mesh.project(proj_grid, lp.mean.pred)
lp.sd.grid = inla.mesh.project(proj_grid, lp.sd.pred)
  
#-- Reconstruct the rain field (on the sqrt scale) --#
if(include.elevation==T){
    rain.mean.grid = lp.mean.grid + elevation.grid * beta[2]
    rain.var.grid = lp.sd.grid^2 + elevation.grid^2 * beta[2]^2
    rain.sd.grid = rain.var.grid^0.5
}else{
    rain.mean.grid = lp.mean.grid
    rain.sd.grid = lp.sd.grid
}
  
#-- Plot the rain field mean (sqrt scale) --#
par(mar=c(4,4,2,2))
image.plot(seq.x.grid,seq.y.grid,rain.mean.grid,
             xlab="W-E (km)", ylab="N-S (km)",col=grey(20:0/20), axes=F)
  contour(seq.x.grid,seq.y.grid,rain.mean.grid, add=T, lwd=2, labcex=1)
  lines(sic.borders,lwd=3)
#-- Plot the rain field sd --#
par(mar=c(4,4,2,2))
image.plot(seq.x.grid,seq.y.grid,rain.sd.grid,
             xlab="W-E (km)", ylab="N-S (km)",col=grey(20:0/20), axes=F)
contour(seq.x.grid,seq.y.grid,rain.sd.grid, add=T,lwd=2, labcex=1)
lines(sic.borders,lwd=3)



