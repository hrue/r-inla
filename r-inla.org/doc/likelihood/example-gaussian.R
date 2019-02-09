n=100
a = 1
b = 1
z = rnorm(n)
eta = a + b*z
tau = 100
scale = exp(rnorm(n))
prec = scale*tau
y = rnorm(n, mean = eta, sd = 1/sqrt(prec))

data = list(y=y, z=z)
formula = y ~ 1+z
result = inla(formula, family = "gaussian", data = data,
        control.family = list(hyper = list(
                                    prec = list(
                                            prior = "loggamma",
                                            param = c(1.0,0.01),
                                            initial = 2))), 
              scale=scale, keep=TRUE)
summary(result)

## with an offset in the variance
var0 = 1.0 ## fixed offset
var1 = 2.0
v = var0 + var1
s = sqrt(v)
x = rnorm(n)
y = 1 + x + rnorm(n, sd = s)
rr = inla(y ~ x,
         data = data.frame(y, x),
         control.family = list(
             hyper = list(precoffset = list(initial = log(1/var0)))), 
         verbose = TRUE)
summary(rr)
plot(rr$internal.marginals.hyperpar[[1]], type = "l", lwd=3)
abline(v = log(1.0/var1), lwd=3, col = "blue")
