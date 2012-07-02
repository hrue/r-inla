library(sn)
n = 1000
z = rnorm(n)
y = z + rsn(n, shape = 2)
formula = y ~ z
r = inla(formula, family = "sn", data = data.frame(z,y),
         control.family = list(sn.shape.max = 5.0))
summary(r)
