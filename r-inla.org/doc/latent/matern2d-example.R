nrow=20
ncol=30
n = nrow*ncol
s.noise = 1

zi.mat = matrix(NA,nrow=nrow,ncol=ncol)
i=1:nrow
for(j in 1:ncol)
    zi.mat[i,j] = 3*exp(-(i-j)^2/4)

## iid noise
noise.mat=matrix(rnorm(nrow*ncol, sd=s.noise),nrow,ncol)

## make simulated data with no spatial component
y.mat = zi.mat + noise.mat

## convert matrices to the internal representation in INLA
y = inla.matrix2vector(y.mat)
node = 1:n
formula= y ~ 1+ f(node, model="matern2d", nu=1, nrow=nrow, ncol=ncol,
        hyper = list(range = list(param =c(1, 1),
                                  prior = "loggamma",
                                  initial=1),
                     prec = list(param=c(1, 1))))
data=data.frame(y=y,node=node)

## fit the model
result=inla(formula, family="gaussian", data=data, verbose=TRUE,
        control.predictor = list(compute = TRUE),
        control.family = list(hyper = list(theta = list(initial = log(1/s.noise^2),
                                                 fixed = FALSE))),
        keep=T)

## plot the posterior mean for `predictor' and compare with the truth
dev.new()
INLA:::inla.display.matrix(zi.mat)
dev.new()
INLA:::inla.display.matrix(INLA:::inla.vector2matrix(result$summary.linear.predictor$mean,nrow,ncol))
