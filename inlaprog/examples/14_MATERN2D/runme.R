nrow = 50
ncol = 50

### normal data with no errors (no nugget); latent field was simulated
### with prec = 1, range = 10

##  idx scale y  (ind: zero-based)
yy = read.table("normal-data-no-error.dat")

ind = yy$V1 + 1  
formula = y ~ f(ind, model="matern2d", nrow=50, ncol=50,nu=1)
d = list(y=yy$V3, ind=ind)
res1 = inla(formula, family = "gaussian", scale=yy$V2, data = d, verbose =TRUE,
            control.data = list(initial = 12, fixed= TRUE), keep=TRUE)

image( matrix(res1$summary.random$ind$mean,nrow,ncol))
title("Posterior mean of the latent Matern-field")


### binary data; latent field was simulated with prec = 1, range = 10

## idx Ntrials y (ind: zero-based)
yy = read.table("binary-data.dat")
ind = yy$V1 + 1
formula = y ~ f(ind, model="matern2d", nrow=50, ncol=50, nu=1)
d = list(y=yy$V3, ind=ind)
res2 = inla(formula, family = "binomial", Ntrials=yy$V2, data = d, verbose =TRUE, keep = TRUE)

dev.new()
image( matrix(res1$summary.random$ind$mean,nrow,ncol))
title("Posterior mean of the latent Matern-field")
