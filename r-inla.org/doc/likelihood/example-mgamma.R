n <- 10^4
phi <- 2.2
x <- rnorm(n)
m <- exp(1.1 + 1.2 * x)
delta <- (sqrt(phi * (phi + 4.0)) + phi) / 2.0

y <- rgamma(n, shape = 1 + delta, rate = delta / m)
r <- inla(y ~ 1 + x,
          data = data.frame(y, x),
          family = "mgamma",
          control.compute = list(cpo = T), 
          verbose = TRUE)

Y <- inla.surv(y, 1)
rr <- inla(Y ~ 1 + x,
          data = list(Y = Y, x = x),
          family = "mgammasurv",
          control.compute = list(cpo = T), 
          verbose = TRUE)
