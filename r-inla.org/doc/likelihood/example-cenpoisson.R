n=100
a = 1
b = 1
z = rnorm(n)
eta = a + b*z
C = 5
E = sample(1:10, n, replace=TRUE)
lambda = E*exp(eta)
y = rpois(n, lambda = lambda)
y[y <= C] = 0

data = list(y=y,z=z)
formula = y ~ 1+z
result = inla(formula, family = "cenpoisson", data = data, E=E,
    control.family = list(cenpoisson.C = C))
summary(result)
