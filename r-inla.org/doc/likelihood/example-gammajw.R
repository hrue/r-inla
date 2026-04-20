n <- 300
x <- rnorm(n, sd = 0.3)
eta <- 1 + x
mu <- exp(eta)
y <- rgamma(n, shape = mu, scale = 1)
r <- inla(y ~ 1 + x,
          data = data.frame(y, x),
          family = "gammajw",
          control.compute = list(cpo = TRUE), 
          control.fixed = list(prec.intercept = 0.01), 
          verbose = TRUE)
summary(r)

yy <- inla.surv(y, event = 1)
rr <- inla(yy ~ 1 + x,
          data = list(yy = yy, x = x),
          family = "gammajwsurv",
          control.compute = list(cpo = TRUE), 
          control.fixed = list(prec.intercept = 0.01), 
          verbose = TRUE)
summary(rr)

print(r$summary.fixed - rr$summary.fixed)
