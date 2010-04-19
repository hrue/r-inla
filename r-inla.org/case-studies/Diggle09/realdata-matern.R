library(INLA)

###read data
lead97=read.table("lead97-new.txt")
names(lead97) <- c("X","Y","data")
lead00=read.table("lead00-new.txt")
names(lead00) <- c("X","Y","data")

###log-scale and scale coordinates
loglead97 <- data.frame(X=lead97$X/100000,Y=lead97$Y/100000,data=log(lead97$data))
loglead00 <- data.frame(X=lead00$X/100000,Y=lead00$Y/100000,data=log(lead00$data))

###define the dimension of the grid
nrow=100
ncol=100
n=nrow*ncol
xrange=range(c(loglead97$X,loglead00$X))
yrange=range(c(loglead97$Y,loglead00$Y))
cutpointsx=seq((xrange[1]-0.1),(xrange[2]+0.1),length.out=nrow)
cutpointsy=seq((yrange[1]-0.1),(yrange[2]+0.1),length.out=ncol)

###work on lead97
grid1.97 <- as.numeric(cut(loglead97$X,cutpointsx))
grid2.97 <- as.numeric(cut(loglead97$Y,cutpointsy))
grid97=cbind(grid1.97,grid2.97)
node97=numeric(dim(loglead97)[1])
for(i in 1: dim(loglead97)[1])
  node97[i]=inla.lattice2node(grid97[i,1],grid97[i,2],nrow=nrow,ncol=ncol)

###work on lead00
grid1.00 <- as.numeric(cut(loglead00$X,cutpointsx))
grid2.00 <- as.numeric(cut(loglead00$Y,cutpointsy))
grid00=cbind(grid1.00,grid2.00)
node00=numeric(dim(loglead00)[1])
for(i in 1: dim(loglead00)[1])
  node00[i]=inla.lattice2node(grid00[i,1],grid00[i,2],nrow=nrow,ncol=ncol)

###Build the data set
y97 = rep(NA,n)
y97[node97] =  loglead97$data
y00 = rep(NA,n)
y00[node00] =  loglead00$data

pois97 = rep(0,n)
pois97[node97] = 1

yy=matrix(NA,4*n,2)
yy[1:n,1] = y97
yy[n+1:n,1 ] = y00
yy[2*n+1:n,2]=pois97

mu0 = c(rep(1,n),rep(0,3*n))
mu1 = c(rep(0,n),rep(1,n),rep(0,2*n))
alpha = c(rep(0,2*n),rep(1,n),rep(0,n))
ii = c(1:n,1:n,rep(NA,2*n))
jj = c(rep(NA,2*n),1:n,1:n)
replicates = c(rep(1,n),rep(2,n),rep(1,n),rep(2,n))

data = list(yy=yy,mu0=mu0,mu1=mu1,alpha=alpha,ii=ii,jj=jj,replicates=replicates)

### define the model and run it
### we choose a matern randon field with range fixed to be 20

formula = yy ~ alpha + mu0 + mu1 -1 + 
    f(ii, model = "matern2d", nrow=nrow, ncol=ncol, replicate=replicates,
      initial = c(2, log(20)), fixed=c(FALSE,TRUE),
      param = c(1,0.001, NA,NA), constr=TRUE) +
    f(jj, copy="ii",replicate=replicates, fixed=FALSE, param=c(0,0.1), initial=0)

res = inla(formula, family = c("gaussian", "poisson"), data =data, verbose=TRUE,
        control.inla= list(strategy = "gaussian"))

### This is to obtain a better estimate of the model hyperparameters
##h.res=inla.hyperpar(res,dz=1,diff.logdens=15,verbose=T)
