## Use the spline part of scopy to estimate a spline.
## This example is rather artifical, but illustrate
## well the idea.
N <- 100
s <- 0.1
x <- 1:N
y <- 1 + sin(x * 0.25) * exp(-2*x/N) + rnorm(N, sd = s)
m <- 15
r <- inla(y ~ -1 +
              ##
              ## this model will just define the overall mean level, but
              ## with one value for each i. We need this as as can then
              ## scale this one with the spline
              ##
              f(idx,
                model = "rw1",
                scale.model = TRUE,
                constr = FALSE,
                values = 1:N,
                hyper = list(prec = list(initial = 20, fixed = TRUE))) +
              ##
              ## the 'overall level' is scaled by a spline
              ##
              f(idx.scopy,
                scopy = "idx",
                control.scopy = list(covariate = x,
                                     n = m,
                                     mean = 1,
                                     prec.mean = 1,
                                     prec.betas = 10,
                                     model = "rw2")), 
          ##
          data = list(idx = rep(NA, N), 
                      idx.scopy = 1:N,
                      x = x,
                      m = m),
          ##
          control.family = list(hyper = list(
                                    prec =  list(
                                        initial = log(1/s^2),
                                        fixed = TRUE))),
          ##
          control.compute = list(config = TRUE),
          verbose = TRUE)
plot(x, y, pch = 19)
## note that the locations are not stored in the results, hence we can set them here. This is
## just for the ease of the output, the results are unchanged. 1. we can just rescale
beta <- inla.summary.scopy(r, "idx.scopy", range = c(1, N))
s <- mean(r$summary.random$idx$mean)
lines(beta$x, s * beta$mean, lwd = 3, col = "blue")
lines(beta$x, s * beta$mean + 2 * s * beta$sd, lwd = 1, col = "black")
lines(beta$x, s * beta$mean - 2 * s * beta$sd, lwd = 1, col = "black")
