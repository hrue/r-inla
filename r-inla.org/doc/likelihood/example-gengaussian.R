library(gnorm)
n <- 10^4
x <- rnorm(n)
sigma <- 2.0
beta <- 1.5
alpha <- sqrt(gamma(1/beta)/gamma(3/beta)) * sigma
y <- 1 + x + rgnorm(n, alpha = alpha, beta = beta)
r <- inla(y ~ 1 + x,
          data = data.frame(y, x),
          family = "gengaussian",
          control.fixed = list(prec.intercept = 1), 
          control.inla = list(cmin = 0))
summary(r)
