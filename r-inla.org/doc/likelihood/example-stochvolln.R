n <- 1000
s <- 0.5
rho <- 0.95
x <- s * as.numeric(scale(arima.sim(n, model = list(ar = rho)))) 
v <- exp(x) ## variance
cc <- 0.5 * mean(v)

y <- rnorm(n, mean = cc - 0.5 * v, sd = sqrt(v))

r <- inla(y ~ -1 + f(time, model = "ar1",
                    hyper = list(rho = list(prior = "pc.cor1",
                                            param = c(0.9, 0.5)),
                                 prec = list(prior = "pc.prec",
                                             param = c(1, 0.01),
                                             initial = 0))), 
          data = list(y = y, n = n, time = 1:n), 
          family = "stochvolln",
          verbose = TRUE, safe = F, keep = T)

par(mfrow = c(3, 2))
plot(y, type = "l", main = "y")
plot(x, type = "l", main = "x and x.est")
lines(r$summary.linear.predictor$mean, col = "blue")

m <- inla.tmarginal(function(x) sqrt(1/x), r$marginals.hyperpar[["Precision for time"]])
plot(m, type = "l", lwd = 3, main = "sd")
abline(v=s)

m <- inla.smarginal(r$marginals.hyperpar[["Rho for time"]])
plot(m, type = "l", lwd = 3, main = "rho")
abline(v=rho)

m <- inla.smarginal(r$marginals.hyperpar[["Mean offset for stochvolln"]])
plot(m, type = "l", lwd = 3, main = "c")
abline(v=cc)
