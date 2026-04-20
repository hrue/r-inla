## the precision parameter in the beta distribution
phi = 5

## generate simulated data
n = 1000
z = rnorm(n, sd=.2)
eta = 1 + z
mu = exp(eta)/(1+exp(eta))
a = mu * phi
b = -mu * phi + phi
y = rbeta(n, a, b)

## this is the censoring
cens <- 0.05
y[y <= cens] <- 0
y[y >= 1-cens] <- 1

## estimate the model
formula = y ~ 1 + z
r = inla(formula, data = data.frame(y, z), family = "beta",
         control.family = list(beta.censor.value = cens))
summary(r)
