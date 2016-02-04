n = 1000
w = runif(n, min = 0.25,  max = 0.75)
phi = 5 * w
z = rnorm(n, sd=0.2)
eta = 1 + z
mu = exp(eta)/(1+exp(eta))
a = mu * phi
b = -mu * phi + phi
y = rbeta(n, a, b)


formula = y ~ 1 + z
r = inla(formula, data = data.frame(y, z, w),
         family = "beta", scale = w)
summary(r)
