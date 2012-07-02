### In this code we simulate data from a model like that in section 5 of
### Diggle et al. 2009

### Load necessary libraries
library(RandomFields)
library(INLA)


### Define the model parameters the grid size and the 
###  number of data points to be simulated 
model <- "exponential"   
mean <- 4
variance <- 1.5
nugget <- 0
scale <- 0.15
beta=2 
mean0=4
mean1=1
sd.data=0.1

nrow=50
ncol=50
x <- seq(0, 1, length.out=nrow) 
y <- seq(0, 1, length.out=ncol)     
n.data = 100
n = nrow*ncol

### simulate 2 realizations of the latent field
S0 <- GaussRF(x=x, y=y, model=model, grid=TRUE,
              param=c(0, variance, nugget, scale))
S0 = S0 - mean(S0) + mean0
S0.vec =  inla.matrix2vector(S0)

S1 <- GaussRF(x=x, y=y, model=model, grid=TRUE,
              param=c(0, variance, nugget, scale))
S1 = S1 - mean(S1) + mean1
S1.vec =  inla.matrix2vector(S1)

###map to probability
p0=exp(beta*S0.vec) / max(exp(beta*S0.vec))
p1=exp(beta*S1.vec) / max(exp(beta*S1.vec))

### Simulate the data on S0 using a preferential sampling design

### sample points
points0 = rep(0,n)
points0[ sample(1:n, size=n.data, prob=p0, replace = FALSE) ] = 1

### observe the marks
y0 = rep(NA,n)
pp0 = which(points0 == 1)
y0[pp0] = rnorm(n.data, mean = S0.vec[pp0], sd=sd.data)

loc0=rep(NA,n)
loc0[pp0]=1

###     Simulate the data on S1 using a grid design      
###    build the grid containing anout n.data points
nn = floor(sqrt(n.data))
x1 = floor(seq(1,nrow,length.out=nn))
x2 = floor(seq(1,ncol,length.out=nn))
xx = expand.grid(x1,x2)
points1 = rep(0,n)
for(i in 1:dim(xx)[1])
points1[inla.lattice2node(xx[i,1],xx[i,2],nrow,ncol)] = 1

### observe the marks
y1 = rep(NA,n)
pp1 = which(points1 == 1)
y1[pp1] = rnorm(sum(points1), mean = S1.vec[pp1], sd=sd.data)

loc1=rep(NA,n)
loc1[pp1]=1


### Prepare data set for the joint model as in sec 5 of Diggle et al.
yy=matrix(NA,4*n,2)
yy[1:n,1] = y0
yy[n+1:n,1] = y1
yy[2*n+1:n,2] = points0

mu0 = c(rep(1,n),rep(0,3*n))
mu1 = c(rep(0,n),rep(1,n),rep(0,2*n))
alpha = c(rep(0,2*n),rep(1,n),rep(0,n))
ii = c(1:n,1:n,rep(NA,2*n))
jj = c(rep(NA,2*n),1:n,1:n)
replicates = c(rep(1,n),rep(2,n),rep(1,n),rep(2,n))

data = list(yy=yy,mu0=mu0,mu1=mu1,alpha=alpha,ii=ii,jj=jj,replicates=replicates)

### write the formula , model with RW2D latent field
formula = yy ~ alpha + mu0 + mu1 -1 + 
    f(ii, model = "rw2d", nrow=nrow, ncol=ncol,
      replicate=replicates ,constr=TRUE, bvalue=1) +
    f(jj, copy = "ii", replicate=replicates,
      fixed=FALSE, param=c(0,0.1))

res = inla(formula, family = c("gaussian", "poisson"), data = data, 
        control.family = list(list(initial=1), list()),
        control.inla= list(strategy = "gaussian", int.strategy = "eb"),
        control.predictor=list(compute=TRUE),verbose = TRUE, keep=TRUE)

###model with matern2d latent field, the range is fixed to be 10
formula.mat = yy ~ alpha + mu0 + mu1 -1 + 
    f(ii, model = "matern2d", nrow=nrow, ncol=ncol,
      replicate=replicates, initial = c(3, log(10)),
      fixed=c(FALSE,TRUE), constr=TRUE) +
    f(jj, copy="ii", replicate=replicates, fixed=FALSE, param=c(0,0.1))

res.mat = inla(formula.mat, family = c("gaussian", "poisson"), data=data,
        control.family = list(list(initial=1), list()),
        control.inla= list(strategy = "gaussian", int.strategy = "eb"),
        control.predictor=list(compute=TRUE),verbose = TRUE)



### This are to get a better estimate of the hyperparameters of the model
### It takes some time to run this functions

##hyper.res <- inla.hyperpar(res,dz=1,diff.logdens=10,verbose=T)
##hyper.res.mat <- inla.hyperpar(res.mat,dz=1,diff.logdens=10,verbose=T)
##plot(hyper.res)
##plot(hyper.res.mat)


### Plot the estimated beta
par(mfrow=c(2,1))
plot(res$marginals.hyperpar$"Beta1 for jj",type="l",
     main="Estimated Beta for the RW2 model")
plot(res.mat$marginals.hyperpar$"Beta1 for jj",type="l",
      main="Estimated Beta for the Matern model")

### Plot the data and the real S fields together with
### the estimated linear predictor from the different models
dev.new()
par(mfrow=c(3,2))
image(S0,col=grey(seq(0,1,len=256)))
image(inla.vector2matrix(loc0,nrow,ncol),add=T)
title("Latent field S0 with data")

image(S1,col=grey(seq(0,1,len=256)))
image(inla.vector2matrix(loc1,nrow,ncol),add=T)
title("Latent field S1 with data")

image(inla.vector2matrix(res$summary.linear.predictor$mean[1:n],nrow,ncol),
      col=grey(seq(0,1,len=256)))
title("linear predictor for S0. rw2 model")

image(inla.vector2matrix(res$summary.linear.predictor$mean[n+1:n],nrow,ncol),
      col=grey(seq(0,1,len=256)))
title("linear predictor for S1. rw2 model")

image(inla.vector2matrix(res.mat$summary.linear.predictor$mean[1:n],nrow,ncol),
      col=grey(seq(0,1,len=256)))
title("linear predictor for S0. Matern model")

image(inla.vector2matrix(res.mat$summary.linear.predictor$mean[n+1:n],nrow,ncol),
      col=grey(seq(0,1,len=256)))
title("linear predictor for S1. Matern model")
