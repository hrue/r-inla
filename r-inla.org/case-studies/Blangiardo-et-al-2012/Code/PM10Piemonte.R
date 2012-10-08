############################################################################################                                  
#                             PM10 CONCENTRATIONS IN PIEMONTE (ITALY)                      #
############################################################################################
remove(list=ls())
require(INLA)
require(fields)
require(abind)

##--- For which day do you want the prediction?
i_day = 122

##--- Path for the Covariates directory
path = "http://www.math.ntnu.no/inla/r-inla.org/case-studies/Cameletti2012/"


##--- Source the function for selecting the covariates for a given day
source(paste(path,"Covariates/Covariates_selector.R",sep=""))


## ################################
## Load the data
## ################################
##--- for the 24 stations and 182 days
Piemonte_data <-read.table(paste(path,"Piemonte_data_byday.csv",sep=""), header=TRUE,sep=",")
coordinates <-read.table(paste(path,"coordinates.csv",sep=""),header=TRUE,sep=",")
rownames(coordinates) = coordinates[,"Station.ID"]

which_date = unique(Piemonte_data$Date)[i_day]
print(paste("**---- You will get a prediction for ", which_date, "---**"))

##--- Borders of Piemonte (in km)
borders <-read.table(paste(path,"Piemonte_borders.csv",sep=""),header=TRUE,sep=",")

## ################################
## Work out how many stations and days there are
## ################################
n_stations <- length(coordinates$Station.ID)
n_data <- length(Piemonte_data$Station.ID)
n_days <- as.integer(n_data/n_stations)

##--- Check that the data is OK
if (n_data %% n_stations != 0) {
  print("The number of data points needs to be an integer multiple of the number of stations!")
  return
}

## ################################
##Standardize covariates and calculate LOG of PM10
## ################################
##--- The covariates are standardised using the mean and std.dev. for
##--- the 24 data sites.
mean_covariates = apply(Piemonte_data[,3:10],2,mean)
sd_covariates = apply(Piemonte_data[,3:10],2,sd)

Piemonte_data[,3:10] =
  scale(Piemonte_data[,3:10],
        mean_covariates, sd_covariates)

Piemonte_data$logPM10 = log(Piemonte_data$PM10)
Piemonte_data$time = rep(1:n_days,each = n_stations)

## ################################
## Load the covariate arrays (each array except for A is 56x72x182)
## ################################
load(url(paste(path,"Covariates/Altitude_GRID.RData",sep=""))) #A; AltitudeGRID
load(url(paste(path,"Covariates/WindSpeed_GRID.RData",sep=""))) #WS; WindSpeedGRID
load(url(paste(path,"Covariates/HMix_GRID.RData",sep=""))) #HMIX; HMixMaxGRID
load(url(paste(path,"Covariates/Emi_GRID.RData",sep=""))) #EMI; EmiGRID
load(url(paste(path,"Covariates/Temp_GRID.RData",sep=""))) #TEMP; Mean_Temp
load(url(paste(path,"Covariates/Prec_GRID.RData",sep=""))) #PREC; Prec

##--- Load the Piemonte grid c(309,529),c(4875,5159),dims=c(56,72)
load(url(paste(path,"Covariates/Piemonte_grid.RData",sep="")))

##--- Extract the standardized covariate for day i_day (you get a 56X72 matrix)
covariate_array_std =
  covariates_selector_funct(i_day, mean_covariates, sd_covariates)

#--- Set to NA the (standardized) altitude values >7
tmp = covariate_array_std[,,1]
index.mountains = which(tmp>7)
tmp[tmp>7] = NA
covariate_array_std[,,1] = tmp

#--- Reshape the 3D array (56x72x8) into a data frame (4032x8) with the 8 covariates on the columns
covariate_matrix_std = data.frame(apply(covariate_array_std,3,function(X) c(t(X))))
colnames(covariate_matrix_std) = colnames(Piemonte_data[,3:10])

## ################################
## Triangulation using borders
## ################################
mesh =
  inla.mesh.create.helper(points=cbind(coordinates$UTMX,
                                       coordinates$UTMY),
                          points.domain=borders,
                          offset=c(10, 140),
                          max.edge=c(50, 1000),
                          min.angle=c(26, 21),
                          cutoff=0,
                          plot.delay=NULL
)

##--- Plot the triangulation
# pdf("triangulation_piemonte.pdf")
# plot(mesh)
# lines(borders, lwd=3)
# points(coordinates$UTMX, coordinates$UTMY, pch=20, cex=2, col=2)
# dev.off()

## ################################
## Make the SPDE object and the formula
## ################################
##--- Construct the SPDE object
spde = inla.spde2.matern(mesh=mesh)

#-- Create A matrices for the estimation part --#
A.est =
  inla.spde.make.A(mesh,
                   loc=as.matrix(coordinates[Piemonte_data$Station.ID, c("UTMX","UTMY")]),
                   group=Piemonte_data$time,
                   n.group=n_days)

field.indices =
  inla.spde.make.index("field",
                       n.mesh=mesh$n,
                       n.group=n_days)

stack.est =
  inla.stack(data=list(logPM10=Piemonte_data$logPM10),
             A=list(A.est, 1),
             effects=
               list(c(field.indices,
                      list(Intercept=1)),
                    list(Piemonte_data[,3:10])),
             tag="est")

#-- Create A matrices for the grid prediction part --#
A.pred =
  inla.spde.make.A(mesh, loc=as.matrix(Piemonte_grid),
                   group=i_day, n.group=n_days)

stack.pred =
  inla.stack(data=list(logPM10=NA),
             A=list(A.pred,1),
             effects=
               list(c(field.indices,
                      list(Intercept=1)),
                    list(covariate_matrix_std)),
                    tag="pred")

#-- Create the "full" stack object (estimation + prediction)
stack = inla.stack(stack.est, stack.pred)


#-- Define the formula --#
formula <- (logPM10 ~ -1 + Intercept + A + UTMX + UTMY + WS + TEMP + HMIX + PREC + EMI + 
  f(field, model=spde, group=field.group, control.group=list(model="ar1")))


## ################################
## Call INLA and get results
## ################################
#Two step procedure:
#Step 1 (compute the mode using stack.est)
mod.mode =
  inla(formula,
       data=inla.stack.data(stack.est, spde=spde),
       family="gaussian",
       control.predictor=list(A=inla.stack.A(stack.est), compute=FALSE),
       control.compute=list(cpo=FALSE),
       keep=FALSE, verbose=TRUE,
       control.inla=list(reordering="metis"))
       
#Step 2 (using result.mode compute the linear predictor in all the grid points in the full stack object)
mod =
  inla(formula,
       data=inla.stack.data(stack, spde=spde),
       family="gaussian",
       control.predictor=list(A=inla.stack.A(stack), compute=TRUE),
       control.compute=list(cpo=TRUE, dic=TRUE),
       control.mode=list(theta=mod.mode$mode$theta, restart=FALSE),
       keep=FALSE, verbose=TRUE,
       control.inla=list(reordering="metis"))
       
print(summary(mod))


## ############################
## Mapping
## ############################
##--- Posterior mean of the linear predictor
index.pred = inla.stack.index(stack,"pred")$data

lp_grid_mean = matrix(
  mod$summary.linear.predictor[index.pred, "mean"], 56, 72, byrow=T)
lp_grid_mean[index.mountains]=NA

lp_grid_sd = matrix(
  mod$summary.linear.predictor[index.pred, "sd"], 56, 72)
lp_grid_sd[index.mountains]=NA

seq.x.grid = seq(range(Piemonte_grid[,1])[1],range(Piemonte_grid[,1])[2],l=56)
seq.y.grid = seq(range(Piemonte_grid[,2])[1],range(Piemonte_grid[,2])[2],l=72)

trellis.par.set(axis.line=list(col=NA))
print(levelplot(x=lp_grid_mean,
                row.values=seq.x.grid,
                column.values=seq.y.grid,
                ylim=c(4875,5159), xlim=c(309,529),
                col.regions=gray(seq(.9,.2,l=100)),
                aspect="iso",
                contour=TRUE, labels=FALSE, pretty=TRUE, 
                xlab="",ylab="", scales = list(draw = FALSE)
))
trellis.focus("panel", 1, 1, highlight=FALSE)
lpoints(borders,col=1,cex=.25)
lpoints(coordinates$UTMX, coordinates$UTMY,col=1,lwd=2,pch=21)
trellis.unfocus()


##-- Extract the posterior marginals for the linear predictor
##-- for computing the (pointwise) exceedance probability 
lp.marginal.val = mod$marginals.linear.predictor[index.pred]
threshold = log(50)
inlaprob  = lapply(X=lp.marginal.val, FUN=function(x) inla.pmarginal(marginal=x,threshold))
tailprob_grid = matrix(1-unlist(inlaprob),56,72, byrow=T)
tailprob_grid[index.mountains] = NA

trellis.par.set(axis.line=list(col=NA))
print(levelplot(x=tailprob_grid,
                row.values=seq.x.grid,
                column.values=seq.y.grid,
                ylim=c(4875,5159), xlim=c(309,529),
                at=seq(0,1,by=.1),
                col.regions=gray(seq(.9,.2,l=100)),
                aspect="iso",
                contour=TRUE, labels=FALSE, pretty=TRUE, 
                xlab="",ylab="", scales = list(draw = FALSE)
                ))
trellis.focus("panel", 1, 1, highlight=FALSE)
lpoints(borders,col=1,cex=.25)
lpoints(coordinates$UTMX, coordinates$UTMY,col=1,lwd=2,pch=21)
trellis.unfocus()


## ############################
## Extract posterior summaries 
## ############################

##--- Results for fixed effects - covariate coeffs -
beta = mod$summary.fixed[,"mean"]
beta_sd = mod$summary.fixed[,"sd"]

##--- ar1 parameter
mod$summary.hyperpar["GroupRho for field",]

##--- sigma2eps (1/precision)
sigma2eps_marg =
  inla.tmarginal(function(x) 1/x,
                 mod$marginals.hyperpar$"Precision for the Gaussian observations")
sigma2eps_m1 = inla.emarginal(function(x) x, sigma2eps_marg)
sigma2eps_m2 = inla.emarginal(function(x) x^2, sigma2eps_marg)
sigma2eps_stdev = sqrt(sigma2eps_m2 - sigma2eps_m1^2)
sigma2eps_quantiles = inla.qmarginal(c(0.025, 0.5, 0.975), sigma2eps_marg)

##--- sigma2omega
result.field = inla.spde.result(mod, "field", spde, do.transform=TRUE)

var.nom.marg = result.field$marginals.variance.nominal[[1]]
var.nom.m1 = inla.emarginal(function(x) x, var.nom.marg)
var.nom.m2 = inla.emarginal(function(x) x^2, var.nom.marg)
var.nom.stdev = sqrt(var.nom.m2 - var.nom.m1^2)
var.nom.quantiles = inla.qmarginal(c(0.025, 0.5, 0.975), var.nom.marg)

##-- range r
range.nom.marg = result.field$marginals.range.nominal[[1]]
range.nom.m1 = inla.emarginal(function(x) x, range.nom.marg)
range.nom.m2 = inla.emarginal(function(x) x^2, range.nom.marg)
range.nom.stdev = sqrt(range.nom.m2 - range.nom.m1^2)
range.nom.quantiles = inla.qmarginal(c(0.025, 0.5, 0.975), range.nom.marg)

print(c(sigma2eps_m1, sigma2eps_stdev, sigma2eps_quantiles))
print(c(var.nom.m1, var.nom.stdev, var.nom.quantiles))
print(c(range.nom.m1, range.nom.stdev, range.nom.quantiles))
print(mod$summary.hyperpar["GroupRho for field",])



