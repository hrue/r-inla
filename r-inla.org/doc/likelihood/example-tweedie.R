library(tweedie)
library(INLA)

n <- 300
x <- rnorm(n, sd = 0.3)
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
          family = "tweedie")
summary(r)
