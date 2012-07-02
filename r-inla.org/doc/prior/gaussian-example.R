n=100
a = 1
b = 1
E = rep(1,n)
z = rnorm(n)
eta = a + b*z
mu = E*exp(eta)
size = 15
prob = size/(size + mu)
y = rnbinom(n, size=size, prob = prob)

data = list(y=y,z=z)
formula = y ~ 1+z
result = inla(formula, family = "nbinomial", data = data, E=E,
              control.family = list(prior="gaussian", param = c(0,0.01)))
summary(result)
