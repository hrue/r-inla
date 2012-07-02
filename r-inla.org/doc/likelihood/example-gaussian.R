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
