n <- 10^4
phi <- 3
a <- phi
## mu <-  phi / rate
alpha <- 0.85
log.mu <- seq(-5, 5, by = 0.01)

## if exp(eta) is the alpha-quantile,  what is then the mu?
log.q <- log(qgamma(alpha, shape = phi, rate =  phi / exp(log.mu)))
fun <- splinefun(log.q, log.mu)

x <- rnorm(n, sd = 0.3)
eta <- 2 + x
mu <- exp(fun(eta))

## just a check
head(cbind(eta, log(qgamma(alpha, shape = phi, rate = phi / mu))))

y <- rgamma(n, shape = phi, rate = phi / mu)
r <- inla(y ~ 1 + x,
          family = "gamma",
          control.family = list(control.link = list(model = "quantile",
                                                    quantile = alpha)),
          data = data.frame(y, x),
          control.fixed = list(prec.intercept = 1, prec = 1), 
          safe = FALSE,
          verbose = TRUE)
summary(r)

cbind(estimate = c(r$summary.fixed$mean, r$summary.hyperpar$mean),
      true = c(2, 1, phi))
