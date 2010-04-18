n=100
a = 1
b = 1
z = rnorm(n)
eta = a + b*z
E = sample(1:10, n, replace=TRUE)
lambda = E*exp(eta)
y = rpois(n, lambda = lambda)

data = list(y=y,z=z)
formula = y ~ 1+z
result = inla(formula, family = "poisson", data = data, E=E)
summary(result)
