n = 300
intercept = 2
x = rnorm(n, sd = 0.2)
beta = 1
eta = intercept + beta * x
alpha = 0.9
y = numeric(n)
E = runif(n, min=1, max=10)
for(i in 1:n) {
    lambda = E[i] * INLA:::inla.qcontpois(exp(eta[i]), alpha = alpha)
    y[i] = rpois(1, lambda)
}

r = inla(y ~ 1 + x,
         data = data.frame(y, x, E),
         family = "qpoisson",
         control.family = list(control.link=list(quantile = alpha)),
         E =E)
summary(r)

