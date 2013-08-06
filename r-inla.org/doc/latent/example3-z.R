## Example for implementing the standard mixed model
##    y = X%*%beta + Z%*%b+e
## with b ~N(0,tau * Q) With Q a precision matrix

library(INLA)
set.seed(123)
n = 100
nz = 5
nb = 4
tmp = matrix(rnorm(nz^2), nz, nz)
Qz = tmp %*% t(tmp) 
Z = matrix(runif(n*nz), n, nz)
b = inla.qsample(1, Q=Qz)
beta = 5 + rnorm(nb)
X = matrix(runif(n*nb), n, nb)
X[,1] = 1   # want an intercept

eta = X %*% beta + Z %*% b
y = eta + rnorm(n, sd=0.01)
Qx = diag(nb)

formula = y ~ -1 + X +  f(id.z, model="z", Cmatrix=Qz,  Z=Z)
result = inla(formula, data = list(y=y, id.z = 1:n, X=X))

plot(b, result$summary.random$id.z$mean[-(1:n)], pch=19)
abline(0, 1)
