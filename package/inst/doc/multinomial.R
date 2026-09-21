## ----setup, include=FALSE-------------------------------------------
set.seed(123)
library(INLA)
inla.setOption(num.threads="1:1")
inla.setOption(smtp="taucs")
if (file.exists("myinit.R")) source("myinit.R")
knitr::opts_chunk$set(echo = TRUE)
knitr::opts_chunk$set(fig.path="figures/multinomial/")

## ----message = FALSE, warning = FALSE-------------------------------
# number of observations 
n = 10

x = rnorm(n, sd = 1)
eta = 1 + x
prob = exp(eta)/(1+exp(eta))

# Sizes
Ntrials = sample(100:1000, n, replace=TRUE)

# Observations
y = rbinom(n, size=Ntrials, prob = prob)

## -------------------------------------------------------------------
lprec = 25

## -------------------------------------------------------------------
sqrt(1/exp(-lprec))

## -------------------------------------------------------------------
mult.model = inla(y ~ -1 + intercept + x,
         data = data.frame(y, x, Ntrials, intercept = rep(1, n)),
         family = "binomial", Ntrials = Ntrials,
         control.predictor = list(hyper = list(prec = list(initial = lprec))),
         control.inla = list(int.strategy="eb"))


## -------------------------------------------------------------------
# initialize
Y = intercept = xx = c()
for(i in 1:n) {
    # corresponding vector of observations
    Y = c(Y, c(y[i], Ntrials[i] - y[i]))
    # corresponding covariates
    xx = c(xx, c(x[i], 0))
    intercept = c(intercept, c(1, 0))
}


## -------------------------------------------------------------------
phi = rep(1:n, each=2)

## -------------------------------------------------------------------
pois.model = inla(Y ~ -1
          + intercept
          + f(phi, model="iid",
              hyper = list(prec = list(initial = -lprec, fixed=TRUE)))
          + xx,
          family="poisson",
          data = data.frame(Y, xx, phi , intercept), 
          control.predictor = list(hyper = list(prec = list(initial = lprec))), 
          control.inla = list(int.strategy="eb"))


## -------------------------------------------------------------------
corr.factor = - sum(dnorm(pois.model$summary.random$phi$mode,
                          sd = sqrt(1/exp(-lprec)), log=TRUE)) + sum(log(Ntrials))
print(mult.model$mlik)
print(pois.model$mlik + corr.factor)

## ----message=F------------------------------------------------------
library(Ecdat)
library(deldir)
library(sp)
library(sf)
library(mvtnorm)
library(gridExtra)
library(mlogit)

## ----message=FALSE--------------------------------------------------
beta = -0.3
deltas = c(1, 4, 3)
gammas = c(0.3, 0.2, 0.4)
param = c(beta, deltas, gammas)
n = 500

set.seed(123)
# alternative specific with generic coefficient beta
X.A = rnorm(n, mean = 30, sd = 2.5)
X.B = rnorm(n, mean = 35, sd = 3.5)
X.C = rnorm(n, mean = 40, sd = 1)

# alternative specific with alternative specific coefficient delta 
W.A = abs(rnorm(n, mean = 1))
W.B = abs(rnorm(n, mean = 1))
W.C = abs(rnorm(n, mean = 1))

# individual specific with alternative specific coefficient gamma 
Z = rnorm(n, 20, 3)


## -------------------------------------------------------------------
Multinom.sample = function(N){
  Y = matrix(NA, ncol = 3, nrow = n)
  for(i in 1:n){
    V.A = beta*X.A[i] + deltas[1]*W.A[i] + gammas[1]*Z[i]  
    V.B = beta*X.B[i] + deltas[2]*W.B[i] + gammas[2]*Z[i] 
    V.C = beta*X.C[i] + deltas[3]*W.C[i] + gammas[3]*Z[i] 
    
    probs = c(V.A, V.B, V.C)
    probs = exp(probs)/sum(exp(probs))
    samp = rmultinom(1, N, prob = probs)
    
    Y[i,] = as.vector(samp)
  }
  colnames(Y ) = c("Y.A", "Y.B", "Y.C")
  return(Y)
}

## -------------------------------------------------------------------
head(Multinom.sample(1), 5)
head(Multinom.sample(100), 5)

## -------------------------------------------------------------------
Y = Multinom.sample(1)
df = data.frame(cbind(Y, X.A, X.B, X.C, W.A, W.B, W.C, Z))

Data.structure = function(df){
  Data = matrix(NA, ncol = 8, nrow = n*3)
  for(i in 1:n){
    # simulated variable
    Data[((i-1)*3+1):(i*3), 1] = c(df$Y.A[i], df$Y.B[i], df$Y.C[i])
    # alternative specific with generic coeff
    Data[((i-1)*3+1):(i*3), 2] = c(df$X.A[i], df$X.B[i], df$X.C[i])
    # alternative specific with alternative coeff
    Data[((i-1)*3+1):(i*3), 3:5] = diag(c(df$W.A[i], df$W.B[i], df$W.C[i]))
    # individual specific with alternative coeff
    Data[((i-1)*3+1):(i*3), 6] = rep(df$Z[i],3)
    # choice situation index
    Data[((i-1)*3+1):(i*3), 7] = rep(i,3)
    # choice alternative index
    Data[((i-1)*3+1):(i*3), 8] = c(1, 2, 3)
  }
  Data = data.frame(Data)
  names(Data) = c('Y', "X","W.A","W.B","W.C","Z",'phi','alt.idx')
  return(Data)
}

round(head(Data.structure(df)),3)

## -------------------------------------------------------------------
formula = Y ~ -1 + X + W.A + W.B + W.C + 
  f(phi, initial = -10, fixed = T) +
  f(alt.idx, Z, fixed = T, constr = T)

## -------------------------------------------------------------------
Data = Data.structure(df)
model = inla(formula, data = Data, family = 'Poisson',
	control.inla=list(int.strategy="eb"))

## -------------------------------------------------------------------
result = rbind(model$summary.fixed[1:5], model$summary.random$alt.idx[2:6])
result = cbind(result, true = param)
row.names(result) = c("beta","delta.A","delta.B","delta.C","gamma.A","gamma.B","gamma.C" )
round(result,3)

## -------------------------------------------------------------------
diff.result = 
  cbind("0.025quant"= diff(model$summary.random$alt.idx$`0.025quant`),
      "0.5quant" = diff(model$summary.random$alt.idx$`0.5quant`),
      "0.975quant" = diff(model$summary.random$alt.idx$`0.975quant`),
      "true" = diff(gammas))
row.names(diff.result) = c("gamma.B - gamma.A", "gamma.C - gamma.B")
round(diff.result,3)

## -------------------------------------------------------------------
random.effect = rnorm(n)
Multinom.sample.rand = function(N, random.effect){
  Y = matrix(NA, ncol = 3, nrow = n)
  for(i in 1:n){
    V.A = beta*X.A[i] + deltas[1]*W.A[i] + gammas[1]*Z[i] + 
      random.effect[i]
    V.B = beta*X.B[i] + deltas[2]*W.B[i] + gammas[2]*Z[i] 
    V.C = beta*X.C[i] + deltas[3]*W.C[i] + gammas[3]*Z[i] 
    
    probs = c(V.A, V.B, V.C)
    probs = exp(probs)/sum(exp(probs))
    samp = rmultinom(1, N, prob = probs)
    
    Y[i,] = as.vector(samp)
  }
  colnames(Y) = c("Y.A", "Y.B", "Y.C")
  return(Y)
}

Y.rand1 = Multinom.sample.rand(100, random.effect)
df.rand1 = data.frame(cbind(Y.rand1, X.A, X.B, X.C, W.A, W.B, W.C, Z))
Data.rand1 = Data.structure(df.rand1)

## -------------------------------------------------------------------
rand.idx = rep(NA, n*3)
rand.idx[seq(1,n*3, by = 3)] = seq(1,n)
Data.rand1$rand.idx = rand.idx
round(head(Data.rand1),3)

## -------------------------------------------------------------------
formula.rand1 = Y ~ -1 + X + W.A + W.B + W.C + 
  f(phi, initial = -10, fixed = T) +
  f(alt.idx, Z, fixed = T, constr = T) + 
  f(rand.idx, model = "iid")  #random effect

## -------------------------------------------------------------------
model.rand1 = inla(formula.rand1, data = Data.rand1, family ='Poisson',
	control.inla=list(int.strategy="eb"))
result.rand1 = rbind(model.rand1$summary.fixed[1:5])
result.rand1 = cbind(result.rand1, true = param[1:4])
row.names(result.rand1) = c("beta","delta.A","delta.B","delta.C")
round(result.rand1,3)

## -------------------------------------------------------------------
diff.result.rand1 = 
  cbind("0.025quant"= diff(model.rand1$summary.random$alt.idx$`0.025quant`),
      "0.5quant" = diff(model.rand1$summary.random$alt.idx$`0.5quant`),
      "0.975quant" = diff(model.rand1$summary.random$alt.idx$`0.975quant`),
      "true" = diff(gammas))
row.names(diff.result.rand1) = c("gamma.B - gamma.A", "gamma.C - gamma.B")
round(diff.result.rand1,3)

## -------------------------------------------------------------------
mean((random.effect - model.rand1$summary.random$rand.idx$`0.5quant`)^2)

## ----fig.height=2.5-------------------------------------------------
#  Generate the random effect
random.walk1 = rep(NA, n)
random.walk1[1] = 0
for(i in 2:n){
  random.walk1[i] = rnorm(1, mean = random.walk1[i-1], sd = 0.1)
}

# Generate data
Y.rw1 = Multinom.sample.rand(100, random.walk1)
df.rw1 = data.frame(cbind(Y.rw1, X.A, X.B, X.C, W.A, W.B, W.C, Z))
Data.rw1 = Data.structure(df.rw1)


## -------------------------------------------------------------------
rw1.idx = rep(NA, n*3)
rw1.idx[seq(1,n*3, by = 3)] = seq(1,n)
Data.rw1$rw1.idx = rw1.idx

formula.rw1 = Y ~ -1 + X + W.A + W.B + W.C + 
  f(phi, initial = -10, fixed = T) +
  f(alt.idx, Z, fixed = T, constr = T) + 
  f(rw1.idx, model = "rw1")  #random walk of order 1

model.rw1 = inla(formula.rw1, family = 'poisson', data = Data.rw1,
	control.inla=list(int.strategy="eb"))

## ----fig.cap=c("Estimate of a random walk, the black line represents the true one and the red line the estimate"), echo=FALSE, fig.height=4----
plot(random.walk1 - mean(random.walk1), type = 'l', 
     main = 'Random walk estimate', ylab = 'values', xlab = 'time')
points(model.rw1$summary.random$rw1.idx$`0.5quant`, type = 'l', col = 'red')

## -------------------------------------------------------------------
x.loc = runif(n)
y.loc = runif(n)
loc = cbind("x.loc" = x.loc, "y.loc" = y.loc)

## ----message=FALSE--------------------------------------------------

voronoi.polygons <- function(x) {
  require(deldir)
  if (.hasSlot(x, 'coords')) {
    crds <- x@coords  
  } else crds <- x
  z <- deldir(crds[,1], crds[,2])
  w <- tile.list(z)
  polys <- vector(mode='list', length=length(w))
  require(sp)
  for (i in seq(along=polys)) {
    pcrds <- cbind(w[[i]]$x, w[[i]]$y)
    pcrds <- rbind(pcrds, pcrds[1,])
    polys[[i]] <- Polygons(list(Polygon(pcrds)), ID=as.character(i))
  }
  SP <- SpatialPolygons(polys)
  voronoi <- SpatialPolygonsDataFrame(SP, 
                                      data=data.frame(x=crds[,1],
                                                      y=crds[,2], 
                                                      row.names=sapply(slot(SP, 'polygons'), 
                                                                       function(x) slot(x, 'ID'))))
}

## ----message=F------------------------------------------------------
# number of regions
n.reg = 50
# corners
boundaries = rbind(c(0,0), c(1,1), c(0,1), c(1,0))
# points
points = rbind(boundaries, cbind(runif(n.reg - 4), runif(n.reg - 4)))
# generate Voronoi polygons
vor = voronoi.polygons(points)

## -------------------------------------------------------------------
# Create an ID column for the regions, it will be our index
vor@data$id = seq(1,n.reg)
# Transform the locations in SpatialPoint
pp = SpatialPoints(loc)
# This vector contains the region in which each point is located
id.samples = over(pp,vor)$id
# Check the number of non epty regions 
check = as.integer(length(unique(id.samples)))
check

## ----fig.cap=c("Partition of the region [0,1]x[0,1], each red point correspond to a different location"), fig.height=5----
# plot the partition and the locations
plot(vor)
points(loc, pch = 21, bg = 2, cex = 0.7)

## ----message = F, warnings = F--------------------------------------
# adjacency matrix
ADJ <- as(st_touches(st_as_sf(vor)), "Matrix")
# Unscaled precision matrix
Q_unscaled <- Matrix::Diagonal(n.reg, colSums(ADJ)) - ADJ * 0.99 # CAR
# marginal inverse standard deviation of each random effect component
SD_inv <- Matrix::Diagonal(n.reg, 1.0^-1)
# Std.dev. from the unscaled precision
SD_scale <- Matrix::Diagonal(n.reg, diag(inla.qinv(Q_unscaled))^0.5)
# covariance matrix
prec.matrix <- (SD_inv %*% SD_scale) %*% Q_unscaled %*% (SD_scale %*% SD_inv)
prec.matrix <- (prec.matrix + t(prec.matrix)) / 2
# generate the random effect and add it to the SpatialPolygonDataFrame
vor@data$rand.eff <- as.numeric(inla.qsample(1, prec.matrix, seed = 1234))


## ----fig.cap=c("Random effect on the partition"), fig.height=4------
spplot(vor, 'rand.eff', col.regions = terrain.colors(32))

## -------------------------------------------------------------------
# generate random effect vector
random.effect = vor@data$rand.eff[id.samples]
# generate sample
Y.spatial = Multinom.sample.rand(N = 100, random.effect)
# construct data set
df.spatial = data.frame(cbind(Y.spatial, X.A, X.B, X.C, W.A, W.B, W.C, Z))
Data.spatial = Data.structure(df.spatial)

## -------------------------------------------------------------------
# create and add index column
rand.eff.idx = rep(NA, n*3)
rand.eff.idx[seq(1,n*3, by = 3)] = id.samples
Data.spatial$rand.eff.idx = rand.eff.idx

# Formula
formula.spatial = Y ~ -1 + X + W.A + W.B + W.C + 
  f(phi, initial = -10, fixed = T) + 
  f(alt.idx, Z, fixed = T, constr = T) + 
  f(rand.eff.idx, model = 'besag', graph = ADJ) #spatial random effect

# Fit the model
model.spatial = inla(formula.spatial, family = "poisson", data = Data.spatial,
	control.inla=list(int.strategy="eb"))

# Check the results
result.spatial = cbind(model.spatial$summary.fixed[3:5], true = param[1:4])
row.names(result.spatial) = c("beta","delta.A","delta.B","delta.C")

diff.result.spatial = 
  cbind("0.025quant"= diff(model.spatial$summary.random$alt.idx$`0.025quant`),
        "0.5quant" = diff(model.spatial$summary.random$alt.idx$`0.5quant`),
        "0.975quant" = diff(model.spatial$summary.random$alt.idx$`0.975quant`),
        "true" = diff(gammas))
row.names(diff.result.spatial) = c("gamma.B - gamma.A", "gamma.C - gamma.B")

round(rbind(result.spatial, diff.result.spatial),3)

## ----fig.cap = c("Comparison between the true random effect and the estimated one"), fig.height=5----
# compute mean square error
mean((vor@data$rand.eff -
        model.spatial$summary.random$rand.eff.idx$`0.5quant`)^2)

# add columns to the SpatialPolygonDataFralme
vor@data$rand.eff.est = model.spatial$summary.random$rand.eff.idx$`0.5quant`
vor@data$residuals = (vor@data$rand.eff-
                        model.spatial$summary.random$rand.eff.idx$`0.5quant`)^2
# plot them
grid.arrange(spplot(vor, 'rand.eff', at = seq(-2.5, 2.5, length.out = 32),
                    col.regions = terrain.colors(32), main = "True"),
             spplot(vor, 'rand.eff.est', at = seq(-2.5,2.5,length.out = 32),
                    col.regions = terrain.colors(32), main = "Estimates"),
             spplot(vor, 'residuals', at = seq(-2.5,2.5,length.out = 32),
                    col.regions = terrain.colors(32), main = "Residuals"),
             nrow = 1)

## -------------------------------------------------------------------
mesh = inla.mesh.2d(loc, max.edge = c(0.02, 0.2), cutoff = 0.02)

## ----fig.cap = c("Triangulation of the study area [0,1]x[0,1], the red dots represent the data points locations"), fig.height=5----
mesh$n
par(mfrow = c(1,1))
plot(mesh)
points(loc, pch = 21, bg = 2)

## -------------------------------------------------------------------
# starting parametrization
size <- min(c(diff(range(mesh$loc[, 1])), diff(range(mesh$loc[, 2]))))
range <- size/5
sigma0 = 1

# new parametrization
kappa0 <- sqrt(8)/range
tau0 <- 1/(sqrt(4 * pi) * kappa0 * sigma0)

# create an spde object
spde <- inla.spde2.matern(mesh, B.tau = cbind(log(tau0), -1, +1),
                          B.kappa = cbind(log(kappa0), 0, -1))

## -------------------------------------------------------------------
Q = inla.spde.precision(spde, theta = c(0,0))
# the number of elements of Q is the square of the number of triangles 
dim(Q)

## -------------------------------------------------------------------
sample = inla.qsample(n = 2, Q, seed = 123)

## -------------------------------------------------------------------
# create A
A = inla.spde.make.A(mesh = mesh,
                      loc = cbind(rep(loc[,1], each = 3), 
                                  rep(loc[,2], each = 3)))

# select rows relative to the first alternatives
A = A*rep(c(1,2,3) == 1, n)

## -------------------------------------------------------------------
random.effect = matrix(drop(A%*%sample[,2]), ncol = 3, byrow = T)[,1]

## ----fig.cap=c("The right figure represents the continuos spatial effect, the left figure represents the value of the random effect on the sample locations" ), eval.after='fig.cap'----
par(mfrow = c(1,2))

# plot the continuos random effect
proj <- inla.mesh.projector(mesh, dims = c(100, 100))
sample_proj = inla.mesh.project(proj, field = sample[,2])
image(proj$x, proj$y, sample_proj , xlim = c(0,1), ylim = c(0,1),
      xlab = '',ylab = '')
contour(proj$x, proj$y, sample_proj, add = T)

# plot the projection on the sample locations
rbPal <- colorRampPalette(heat.colors(100))
Col <- rbPal(100)[as.numeric(cut(random.effect, breaks = 100))]
plot(loc, pch = 20, col = Col, xlab = '', ylab = '')

## -------------------------------------------------------------------
# generate sample
Y.spde = Multinom.sample.rand(N = 100, random.effect)
# construct data set
df.spde = data.frame(cbind(Y.spde, X.A, X.B, X.C, W.A, W.B, W.C, Z))
Data.spde = Data.structure(df.spde)
round(head(Data.spde),3)

## -------------------------------------------------------------------
# index for the random effect
s_index <- inla.spde.make.index(name="spatial.field",
                                n.spde= spde$n.spde)

# stack object
stack = inla.stack(data = list(Y = Data.spde$Y),
                   effects = list(s_index, list(Data.spde[,2:8])),
                   A = list(A, 1))


## -------------------------------------------------------------------
formula.spde = Y ~ -1 + X + W.A + W.B + W.C + 
  f(phi, initial = -10, fixed = T) +
  f(alt.idx, Z, fixed = T, constr = T) +
  f(spatial.field, model = spde, group = spatial.field.group)

init = c(-0.008, -0.093)
model.spde <- inla(formula.spde,
                   data=inla.stack.data(stack, spde = spde),
                   family="poisson",
                   control.predictor=list(A = inla.stack.A(stack), compute=TRUE),
                   control.inla= list(int.strategy = "eb"),
                   control.mode=list(restart=T, theta=init))

round(model.spde$internal.summary.hyperpar$mode, 3)

## -------------------------------------------------------------------
# estimates on the mesh locations
spde.mesh.est = model.spde$summary.random$spatial.field$`0.5quant`
# estimates on the data locations
spde.loc.est = matrix(drop(A%*%spde.mesh.est),
                 ncol=3, byrow=TRUE)[,1]
# Mean Square Error
mean((spde.loc.est - random.effect)^2)

## ----fig.cap=c("Recostruction of the spatial random field")---------
par(mfrow = c(1,2))

# plot the estimates
output_proj = inla.mesh.project(proj, field = spde.mesh.est)
image(proj$x, proj$y, output_proj , xlim = c(0,1), ylim = c(0,1),
      xlab = '',ylab = '', main = "Estimate")

# plot the true one
sample_proj = inla.mesh.project(proj, field = sample[,2])
image(proj$x, proj$y, sample_proj , xlim = c(0,1), ylim = c(0,1),
      xlab = '',ylab = '', main = "True")


## ----message=F------------------------------------------------------
data("Yogurt")
head(Yogurt)

## ----echo = F-------------------------------------------------------
n = nrow(Yogurt)
Data = matrix(NA, ncol = 5, nrow = n*4)
for(i in 1:n){
  choice = rep(0, 4)
  choice[Yogurt$choice[i]] = 1
  Data[((i-1)*4+1):(i*4), 1] = choice
  Data[((i-1)*4+1):(i*4), 2] = c(Yogurt$price.yoplait[i],
                                 Yogurt$price.dannon[i],
                                 Yogurt$price.hiland[i],
                                 Yogurt$price.weight[i])
  Data[((i-1)*4+1):(i*4), 3] = c(Yogurt$feat.yoplait[i],
                                 Yogurt$feat.dannon[i],
                                 Yogurt$feat.hiland[i],
                                 Yogurt$feat.weight[i])
  Data[((i-1)*4+1):(i*4), 4] = rep(Yogurt$id[i], 4)
  Data[((i-1)*4+1):(i*4), 5] = rep(i, 4)
}

Data = data.frame(Data)
names(Data) = c('Y', 'price', 'feat', 'cust.id', 'phi')
Data$price = Data$price/100 
Data$phi = as.factor(Data$phi)
Data$alpha.idx = rep(c("yoplait","dannon", "hiland","weight"), n)

## -------------------------------------------------------------------
head(Data)

## -------------------------------------------------------------------
formula = Y ~ -1 + price + feat + alpha.idx +
  f(phi, initial = -10, fixed = T) 

model = inla(formula, data = Data, family = 'Poisson',
	control.inla=list(int.strategy="eb"))

## -------------------------------------------------------------------
results = model$summary.fixed[,1:5]
rownames(results) = c("price", "feat", "dannon", "hiland", "weight", "yoplait")
round(results,3)

## -------------------------------------------------------------------
Data$price.class = inla.group(Data$price, n = 12, method = "cut")

formula.price = Y ~ -1 + feat + alpha.idx +
  f(price.class, model = "rw1") + 
  f(phi, initial = -10, fixed = T) 

model.price = inla(formula.price, data = Data, family = 'Poisson',
	control.inla=list(int.strategy="eb"))

res = model.price$summary.fixed[,1:5]
rownames(res) = c('feat', 'dannon','hiland','weight','yoplait')
round(model.price$summary.random$price.class[,1:5],3)

## ----fig.cap=c("Relation between the coefficient and the price level")----
plot(model.price$summary.random$price.class$ID,
     model.price$summary.random$price.class$mean, type = 'l', 
     ylab = expression(beta[1]),
     xlab = 'price')

## ----echo = F-------------------------------------------------------
Data = matrix(NA, ncol = 8, nrow = n*4)
for(i in 1:n){
  choice = rep(0, 4)
  choice[Yogurt$choice[i]] = 1
  Data[((i-1)*4+1):(i*4), 1] = choice
  Data[((i-1)*4+1):(i*4), 2] = c(Yogurt$price.yoplait[i],
                                 Yogurt$price.dannon[i],
                                 Yogurt$price.hiland[i],
                                 Yogurt$price.weight[i])
  Data[((i-1)*4+1):(i*4), 3:6] = rbind(c(Yogurt$feat.yoplait[i], NA, NA, NA),
                                       c(NA, Yogurt$feat.dannon[i], NA, NA),
                                       c(NA, NA, Yogurt$feat.hiland[i], NA),
                                       c(NA, NA, NA, Yogurt$feat.weight[i]))
                                
  Data[((i-1)*4+1):(i*4), 7] = rep(Yogurt$id[i], 4)
  Data[((i-1)*4+1):(i*4), 8] = rep(i, 4)
}


Data = data.frame(Data)
names(Data) = c('Y', 'price', 'feat.yoplait','feat.dannon', 'feat.hiland', 'feat.weight', 'cust.id', 'phi')
Data$price = Data$price/100 #(Data$price - mean(Data$price))/sd(Data$price)
Data$phi = as.factor(Data$phi)
Data$alpha.idx = rep(c("yoplait","dannon", "hiland","weight"), n)
head(Data)

## -------------------------------------------------------------------
formula.feat = Y ~ -1 + price + feat.yoplait + feat.dannon + 
  feat.hiland + feat.weight + alpha.idx + 
  f(phi, initial = -10, fixed = T)  
  

model.feat = inla(formula.feat, data = Data, family = 'Poisson',
	control.inla=list(int.strategy="eb"))

## -------------------------------------------------------------------
model.feat$summary.fixed[2:5, 1:5]

