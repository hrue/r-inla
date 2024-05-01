n = 1000
x = rnorm(n)
eta = 1 + x
mu = exp(eta)
prec.scale = runif(n, min = 0.5, max = 2)
prec.par = 1.2
a = prec.par * prec.scale
b = mu / (prec.par * prec.scale)
y = rgamma(n, shape = a, scale = b)
r = inla(y ~ 1 + x,  data = data.frame(y, x),
        scale = prec.scale, family = "gamma")
