nrow=20
ncol=30
n = nrow*ncol

## two covariates
zi.mat = matrix(NA,nrow=nrow,ncol=ncol)
i=1:nrow
for(j in 1:ncol)
    zi.mat[i,j] = rnorm(nrow, mean = i, sd=1)

zj.mat = matrix(NA,nrow=nrow,ncol=ncol)
j=1:ncol
for(i in 1:nrow)
    zj.mat[i,j] = rnorm(ncol, mean = j, sd=1)

## iid noise
noise.mat=matrix(rnorm(nrow*ncol, sd=1),nrow,ncol)

## make simulated data with no spatial component
y.mat = zi.mat + zj.mat + noise.mat

## convert matrices to the internal representation in INLA
y = inla.matrix2vector(y.mat)
zi = inla.matrix2vector(zi.mat)
zj = inla.matrix2vector(zj.mat)
node = 1:n
formula= y ~ 1+zi+zj + f(node, model="matern2d", nu=1, nrow=nrow, ncol=ncol,
                         param=c(NA,NA,1,1))
data=data.frame(y=y,node=node,zi=zi,zj=zj)

## fit the model
result=inla(formula, family="gaussian", data=data, verbose=TRUE,
            control.predictor = list(compute = TRUE))

#plot the posterior mean for `predictor' and compare with the truth
par(mfrow=c(2,1))
image(s.mat)
image(inla.vector2matrix(result$summary.linear.predictor$mean,nrow,ncol))
