library(tweedie)
library(INLA)

n <- 300
x <- rnorm(n, sd = 0.2)
eta <- 1 + x
mu <- exp(eta)

p <- 1.32
phi <- 2.0
y <- numeric(n)
for(i in 1:n) {
    y[i] <- rtweedie(1, xi = p, mu = mu[i], phi = phi)
}

r <- inla(y ~ 1 + x,
          data = data.frame(y, x),
          ## offset = rep(log(mean(y)), n), 
          family = "tweedie",
          control.family = list(hyper =  list(
                                    theta1 = list(initial = 0), 
                                    theta2 = list(initial = -4,
                                                  prior = "loggamma",
                                                  param = c(100, 100)))), 
          ## control.fixed = list(prec = 0, prec.intercept = 1), 
          ## control.inla = list(cmin = 0, b.strategy = "skip"), 
          ## inla.mode = "experimental", 
          num.threads = "4:1", 
          verbose = T)
summary(r)
