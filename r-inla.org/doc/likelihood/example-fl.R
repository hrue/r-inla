## example with Poisson likelihood in two ways
n <- 100
x <- rnorm(n)
eta <- 1 + 0.3 * x
y <- rpois(n, exp(eta))

## log-likelihood is in the fl-parameterisation:
## log(y!) + y * eta - 0.5 * 0 * (0 - eta)^2 - 1 * exp(0 + 1 *eta)
C <- cbind(-lfactorial(y), y, 0, 0, 1, 0, 1, 0, 0, 0)

r <- inla(inla.mdata(C) ~ 1 + x, family = "fl", data = list(C = C, x = x))
rr <- inla(y ~ 1 + x, family = "poisson", data = data.frame(y, x))

summary(rr)
summary(r)
r$mlik - rr$mlik
