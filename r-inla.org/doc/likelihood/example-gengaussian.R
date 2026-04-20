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
          control.compute = list(cpo = TRUE), 
          control.fixed = list(prec.intercept = 1), 
          control.inla = list(cmin = 0))
summary(r)

n <- 10^5
x <- rnorm(n)
sigma <- 2.0
beta <- 1.5
alpha <- sqrt(gamma(1/beta)/gamma(3/beta)) * sigma
## this is the lin.pred for the quantile
eta.q <- 1 + x  
quantile <- 0.9
## this is the mu/mean/median-parameter in the qgnorm
mu <- eta.q - qgnorm(quantile, alpha = alpha, beta = beta)
y <- rgnorm(n, mu = mu, alpha = alpha, beta = beta)
rr <- inla(y ~ 1 + x,
           data = data.frame(y, x),
           family = "gengaussian",
           control.compute = list(cpo = TRUE), 
           control.family = list(control.link = list(model = "quantile",
                                                     quantile = quantile)), 
           control.fixed = list(prec.intercept = 1), 
           control.inla = list(cmin = 0, int.strategy = "eb"))
summary(r)

