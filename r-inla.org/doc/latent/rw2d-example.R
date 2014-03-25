nrow=50
ncol=25
n = nrow*ncol
s.mat=matrix(NA,nrow=nrow,ncol=ncol)
j=1:ncol
for(i in 1:nrow)
    s.mat[i,j] = 0.1*(i+2*j)

## a covariate
z.mat=matrix(runif(nrow*ncol),nrow,ncol)

## noise
noise.mat=matrix(rnorm(nrow*ncol, sd=0.3),nrow,ncol)

## make simulated data
y.mat = s.mat  + 0.5*z.mat + noise.mat

## convert matrices to the internal representation in INLA
y = inla.matrix2vector(y.mat)
z = inla.matrix2vector(z.mat)
node = 1:n
formula= y ~ z + f(node, model="rw2d", nrow=nrow, ncol=ncol)
data=data.frame(y=y,z=z,node=node)
## fit the model
result=inla(formula, family="gaussian", data=data)

#plot the posterior mean for `node' with the truth
dev.new()
INLA:::inla.display.matrix(s.mat)
dev.new()
INLA:::inla.display.matrix(INLA:::inla.vector2matrix(result$summary.random$node$mean,nrow,ncol))
