library(sn)
n = 1000
z = rnorm(n)
y = z + rsn(n, shape = 3)
formula = y ~ z
r = inla(formula, family = "sn", data = data.frame(z,y))
summary(r)

