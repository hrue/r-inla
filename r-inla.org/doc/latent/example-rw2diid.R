### Example: Use of the rw2diid-model ###
### BCI-data: The point pattern of Beilschmiedia pendula (bei) with covariates (bei.extra) ##

library(fields)
library(INLA)
library(plotrix)
library(spatstat)

# The observed spatial point pattern
data.ppp = bei
plot(data.ppp, main="")

# The observation window
x.win = data.ppp$window

# Covariates, here elevation and slope
cov1 = bei.extra$elev
cov2 = bei.extra$grad

# Choosing grid resolution, here using a VERY coarse resolution to speed up calculations
nrow = 20
ncol = 40
n = nrow*ncol

# Vector giving the area for each grid cell
Area = rep(max(x.win$xrange)*max(x.win$yrange)/n, n)

# Discretize the observation window
x.grid = quadrats(x.win, ncol, nrow)

# Count the number of points in each grid cell, starting at upper left corner and downwards
y = as.vector(quadratcount(data.ppp, tess = x.grid))

# Extract values for the covariates at the center of each grid cell
func.centers <- function(x.win, nrow, ncol)
{
  eps = c((max(x.win$x)-min(x.win$x))/ncol/2, (max(x.win$y)-min(x.win$y))/nrow/2)
  x.x = matrix(rep(seq(min(x.win$x)+eps[1], max(x.win$x)-eps[1],length=ncol),nrow),nrow,ncol,byrow=T)
  x.y = matrix(rep(seq(max(x.win$y)-eps[2], min(x.win$y)+eps[2],length=nrow),ncol),nrow,ncol)
  x.ppp = ppp(x.x, x.y, window=x.win)
  return(x.ppp)
}
centers = func.centers(x.win, nrow, ncol)
cov1 = lookup.im(cov1, centers$x, centers$y)
cov2 = lookup.im(cov2, centers$x, centers$y)

# Scale covariates, possibly also using log-transform if distribution is skewed
cov1 = scale(c(cov1))
cov2 = scale(log(c(cov2)))

# Prepare dataset
data = list(y=y, cov1=cov1, cov2=cov2, index=seq(1:n))

# Choose parameters for the PC priors for hyperparameters
U.sigma = c(0.10, 0.5, 1.0, 2.0)
alpha.sigma = 0.01
U.phi = 0.5
alpha.phi = 2/3

# Run inla using different values of U.sigma
formula = result = list()
for (j in 1:length(U.sigma)){
  ptm <- proc.time()
  print(j)
  formula[[j]] = y ~ cov1 + cov2 +
    f(index, nrow=nrow, ncol=ncol, model="rw2diid", scale.model=T,
      hyper=list(prec=list(prior="pc.prec", param=c(U.sigma[j], alpha.sigma)),
                 phi =list(prior="pc", param=c(U.phi, alpha.phi))))
  result[[j]] = inla(formula[[j]], family="poisson", data=data, E=Area, 
                     verbose=F, control.compute=list(dic=TRUE))
  print(proc.time()-ptm)
}

# Plot estimated credible intervals for different values of U.sigma
func.plotci <- function(result, U, cov.names){
  n.covariates = length(cov.names)
  y.min = y.mean = y.max = matrix(nrow=length(U), ncol=n.covariates)
  for (j in 1:length(U)){
    y.min[j,] = result[[j]]$summary.fixed[-1,3]
    y.mean[j,]= result[[j]]$summary.fixed[-1,1]
    y.max[j,] = result[[j]]$summary.fixed[-1,5]
  }
  pos1 = seq(0.8, 1.2, length = length(U))
  pos = pos1 
  for (i in 1:(n.covariates-1)) pos = c(pos, pos1+i*1)
  plotCI(pos, y=y.mean, ui=y.max, li=y.min,
         err="y", xlab="", ylab="", xaxt="n", col=1:length(U), lwd=1.5)
  axis(1, at=seq(1:n.covariates), labels=cov.names, cex=2)
  abline(h=0, lwd=2, col=c(1))
}
func.plotci(result, U.sigma, cov.names = c("Elevation","Slope"))

# The posterior means and credible intervals for the hyperparameters. In addition DIC. 
func.summary <- function(result){
  n.res = length(result)
  tab = matrix(nrow=n.res, ncol=8)
  for (j in 1:length(result)){
    tab[j,] = c(U.sigma[j],
                inla.emarginal(function(x) 1/x^0.5, result[[j]]$marginals.hyper[[1]]),
                inla.qmarginal(p=c(0.025,0.975),
                               inla.tmarginal(function(x) 1/x^0.5, result[[j]]$marginals.hyper[[1]])),
                as.numeric(result[[j]]$summary.hyper[2,c(1,3,5)]), 
                result[[j]]$dic$dic)
    colnames(tab) = c("U","sigma","0.025quant","0.975quant","phi","0.025quant","0.975quant","DIC")
  }
  return(tab)
}
func.summary(result)

