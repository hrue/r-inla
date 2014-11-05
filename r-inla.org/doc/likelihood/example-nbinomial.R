n = 1000
x = rnorm(n,  sd = 0.2)
eta = 1 + x
E = runif(n, min = 0, max=10)

mu = E * exp(eta)
size = 3
y = rnbinom(n, size=size, mu=mu)
r = inla(y ~ 1 + x, data = data.frame(y, x, E),
    family = "nbinomial", E=E)

mu = E * exp(eta)
size = E*3
y = rnbinom(n, size=size, mu=mu)
rr = inla(y ~ 1 + x, data = data.frame(y, x, E),
    family = "nbinomial",
    control.family = list(variant = 1),
    E=E)


