### Simulation study similar to that in section 3 of 
### "Geostatistical inference under preferential sampling" Diggle et al. 2009

####Load the necessary libraries
library(INLA)
library(RandomFields)

### Define the model parameters the grid size and the
###  number of data points to be simulated 
model <- "exponential"   
mean <- 4
variance <- 1.5
nugget <- 0
scale <- 0.15
beta=2 
sd.data=0

nrow=50
ncol=50
x <- seq(0, 1, length.out=nrow) 
y <- seq(0, 1, length.out=ncol)     
n.data = 100
n = nrow*ncol

###simulate the latent field

S0 <- GaussRF(x=x, y=y, model=model, grid=TRUE,
              param=c(mean, variance, nugget, scale))
S0.vec =  inla.matrix2vector(S0)

### map to probability 
p = exp(beta*S0.vec) / max(beta*exp(S0.vec))

###  SAMPLE THE DATA USING DIFFERENT DESIGNS
###Completely random sampling design

### sample points
points.rand = rep(0,n)
points.rand[ sample(1:n, size=n.data, prob=rep(1/n,n),replace = FALSE)] = 1

### observe the marks
y.rand = rep(NA,n)
pp.rand = which(points.rand == 1)
y.rand[pp.rand] = rnorm(n.data, mean = S0.vec[pp.rand], sd=sd.data)

###Preferential sample design
### sample points
points.pref = rep(0,n)
points.pref[ sample(1:n, size=n.data, prob=p, replace = FALSE) ] = 1

### observe the marks
y.pref = rep(NA,n)
pp.pref = which(points.pref == 1)
y.pref[pp.pref] = rnorm(n.data, mean = S0.vec[pp.pref], sd=sd.data)

### Cluster design 

###simulate an independent realization of S
S1 = GaussRF(x=x, y=y, model=model, grid=TRUE,
              param=c(mean, variance, nugget, scale))
S1.vec =  inla.matrix2vector(S1)
p1 = exp(beta*S1.vec) / max(beta*exp(S1.vec))

### sample points
points.clust = rep(0,n)
points.clust[ sample(1:n, size=n.data, prob=p1, replace = FALSE) ] = 1

### observe the marks
y.clust = rep(NA,n)
pp.clust = which(points.clust == 1)
y.clust[pp.clust] = rnorm(n.data, mean = S0.vec[pp.clust], sd=sd.data)


###Plot the data locations
loc.pref=rep(NA,n)
loc.pref[pp.pref]=1
loc.rand=rep(NA,n)
loc.rand[pp.rand]=1
loc.clust=rep(NA,n)
loc.clust[pp.clust]=1
dev.new()
par(mfrow=c(2,2))
image(S0,col=grey(seq(0,1,len=256)))
image(inla.vector2matrix(loc.pref,nrow,ncol),add=T)
title("Preferential sampling")

image(S0,col=grey(seq(0,1,len=256)))
image(inla.vector2matrix(loc.rand,nrow,ncol),add=T)
title("Random sampling")

image(S0,col=grey(seq(0,1,len=256)))
image(inla.vector2matrix(loc.clust,nrow,ncol),add=T)
title("Cluster sampling")

### Fit a simple geostatistical model as in equation (1) in Diggle et
### al as a model for the latent random field we use first a RW2 and
### then a matern field

### Prepare the data sets
data.simple.pref=list(y=y.pref,ii=seq(1,n))
data.simple.rand=list(y=y.rand,ii=seq(1,n))
data.simple.clust=list(y=y.clust,ii=seq(1,n))

###define the model formula (here we use a RW2 prior)
formula.simple = y~f(ii,model="rw2d",nrow=nrow,ncol=ncol)

###fit the three models
simple.model.pref = inla(formula.simple,family="gaussian",data=data.simple.pref,
                          control.predictor = list(compute = TRUE), verbose=T)

simple.model.rand = inla(formula.simple,family="gaussian",data=data.simple.rand,
                         control.predictor = list(compute = TRUE),verbose=T)

simple.model.clust = inla(formula.simple,family="gaussian",data=data.simple.clust,
                          control.predictor = list(compute = TRUE),verbose=T)


###define the model formula (here we use a matern prior and fix the range to be 10)
formula.simple1 = y~f(ii,model="matern2d",nrow=nrow,ncol=ncol,initial=c(1,log(10)),
                      fixed=c(FALSE,TRUE))

###fit the three models
simple.model.pref1 = inla(formula.simple1,family="gaussian",data=data.simple.pref,
                          control.predictor = list(compute = TRUE), verbose=T)

simple.model.rand1 = inla(formula.simple1,family="gaussian",data=data.simple.rand,
                         control.predictor = list(compute = TRUE),verbose=T)

simple.model.clust1 = inla(formula.simple1,family="gaussian",data=data.simple.clust,
                          control.predictor = list(compute = TRUE),verbose=T)

### Model with preferential sampling 
### as in eq (5) of Diggle et al. 
### Prepare the data sets
ii=c(1:n,rep(NA,n))     
jj = c(rep(NA,n), 1:n)
alpha = c(rep(0,n), rep(1,n))
mu = c(rep(1,n), rep(0,n))

### preferential sampling
yy.pref = matrix(NA,2*n,2)
yy.pref[1:n,1] = y.pref
yy.pref[n+1:n,2] = points.pref

data.pref.pref = list(yy=yy.pref,mu=mu,ii=ii,jj=jj,alpha=alpha)

### random sampling
yy.rand = matrix(NA,2*n,2)
yy.rand[1:n,1] = y.rand
yy.rand[n+1:n,2] = points.rand

data.pref.rand = list(yy=yy.rand,mu=mu,ii=ii,jj=jj,alpha=alpha)

### cluster sampling
yy.clust = matrix(NA,2*n,2)
yy.clust[1:n,1] = y.clust
yy.clust[n+1:n,2] = points.clust

data.pref.clust = list(yy=yy.clust,mu=mu,ii=ii,jj=jj,alpha=alpha)

### define the model formula (here with RW prior)
formula.pref = yy ~ alpha + mu + f(ii, model = "rw2d", nrow=nrow, ncol=ncol,
        initial = 5, fixed=c(FALSE), param = c(1,0.001), constr=TRUE, bvalue=1) +
    f(jj, copy="ii", fixed=FALSE, param=c(0,0.1)) -1

pref.model.pref = inla(formula.pref, family = c("gaussian", "poisson"),
        control.family = list(list(initial = 4, fixed=FALSE), list()),
        control.inla= list(strategy = "gaussian"),
        data = data.pref.pref, verbose = TRUE)

pref.model.rand = inla(formula.pref, family = c("gaussian", "poisson"),
        control.family = list(list(initial = 4, fixed=FALSE), list()),
        control.inla= list(strategy = "gaussian"),
        data = data.pref.rand, verbose = TRUE)

pref.model.clust = inla(formula.pref, family = c("gaussian", "poisson"),
        control.family = list(list(initial = 4, fixed=FALSE), list()),
        control.inla= list(strategy = "gaussian"),
        data = data.pref.clust, verbose = TRUE)


### define the model formula (here with Matern prior, range is fixed to be 10)
formula.pref1 = yy ~ alpha + mu + f(ii, model = "matern2d", nrow=nrow, ncol=ncol,
        initial = c(3, log(10)), fixed=c(FALSE,TRUE),
        param = c(1,0.001, NA,NA), constr=TRUE) +
    f(jj, copy="ii", fixed=FALSE, param=c(0,0.1)) -1

pref.model.pref1 = inla(formula.pref1, family = c("gaussian", "poisson"),
        control.family = list(list(initial = 4, fixed=FALSE), list()),
        control.inla= list(strategy = "gaussian"),
        data = data.pref.pref, verbose = TRUE)

pref.model.rand1 = inla(formula.pref1, family = c("gaussian", "poisson"),
        control.family = list(list(initial = 4, fixed=FALSE), list()),
        control.inla= list(strategy = "gaussian"),
        data = data.pref.rand, verbose = TRUE)

pref.model.clust1 = inla(formula.pref1, family = c("gaussian", "poisson"),
        control.family = list(list(initial = 4, fixed=FALSE), list()),
        control.inla= list(strategy = "gaussian"),
        data = data.pref.clust, verbose = TRUE)

