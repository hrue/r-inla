if (FALSE) {
    inla.setOption(smtp = 'taucs', safe = FALSE, num.threads = "1:1")
}

N <- 200
s <- 0.1
x <- 1:N
eta <- 1 + x / N + sin(x * 0.1) * exp(-2*x/N) 
y <- eta + rnorm(N, sd = s)
m <- 15

Y <- matrix(NA, N + 1, 2)
Y[1:N, 1] <- y
Y[N+1, 2] <- mean(y)
r <- inla(Y ~ -1 +
              ## this model will just define the 'overall level',
              ## but with one value for each i.
              ## We need this as as can then scale this one with
              ## the spline. We add a point with
              ## the second likelihood to lock-it in place
              f(idx,
                model = "rw1",
                scale.model = TRUE,
                constr = FALSE,
                values = 1:N,
                hyper = list(prec = list(initial = 15,
                                         fixed = TRUE))) +
              ## the 'overall level' is scaled by a spline
              f(idx.scopy, scopy = "idx",
                hyper = list(mean = list(param = c(1, 0.1)),
                             slope = list(param = c(0, 0.1)),
                             spline = list(param = c(0, 20))), 
                control.scopy = list(covariate = x, n = m)), 
          ##
          data = list(idx = c(rep(NA, N), 1), 
                      idx.scopy = c(1:N, NA), 
                      x = c(x,1), 
                      m = m),
          family = c("gaussian", "gaussian"), 
          control.family = list(
              list(hyper =
                       list(prec = list(
                                initial = log(1/s^2),
                                fixed = TRUE))),
              list(hyper =
                       list(prec = list(
                                initial = 15, 
                                fixed = TRUE)))), 
          control.compute = list(config = TRUE, residuals = TRUE))

plot(x, y, pch = 19)
## note that the locations are not stored in the results, hence
## we can set them here. This is just for the ease of the output,
## the results are unchanged.
beta <- inla.scopy.summary(r, "idx.scopy", range = c(1, N))
s.mean <- mean(r$summary.random$idx$mean)
lines(beta$x, s.mean * (beta$mean), lwd = 3, col = "blue")
lines(beta$x, s.mean * (beta$mean + 2 * beta$sd), lwd = 2,
      lty = 2, col = "black")
lines(beta$x, s.mean * (beta$mean - 2 * beta$sd), lwd = 2,
      lty = 2, col = "black")
