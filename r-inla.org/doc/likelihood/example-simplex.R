## this library is found at
## http://www.commanster.eu/rcode.html
library(rmutil)

n = 1000
x = rnorm(n, sd = 0.2)
eta = 1 + x
mu = exp(eta)/(1+exp(eta))

s = 0.3
y = rsimplex(n, m = mu, s = s)

r = inla(y ~ 1 + x,  data = data.frame(y, x),
        family = "simplex")
## prec = 1/s
summary(r)


