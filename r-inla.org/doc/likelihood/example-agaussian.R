n <- 10
y <- rnorm(n)
Y <- inla.agaussian(y)

r <- inla(Y ~ 1,
          data = list(Y = Y),
          family = "agaussian")
rr <- inla(y ~ 1,
           data = data.frame(y), 
           family = "gaussian")
print(r$mlik - rr$mlik)

inla.dev.new()
par(mfrow = c(1, 2))
plot(r$internal.marginals.hyperpar[[1]], pch = 19, main = "prec")
lines(rr$internal.marginals.hyperpar[[1]], lwd = 3)
plot(r$marginals.fixed$'(Intercept)', pch = 19, main = "intercept")
lines(rr$marginals.fixed$'(Intercept)', lwd = 3)

####################################
####################################

n <- 5
s <- 1:n ## scale the precision
y <- rnorm(n)
yy <- rnorm(n, sd = sqrt(1/s))
Y <- inla.agaussian(rbind(y, yy),
                    rbind(rep(1, n), s))

r <- inla(Y ~ 1,
          data = list(Y = Y),
          family = "agaussian",
          control.compute = list(cpo = TRUE, dic = TRUE))
rr <- inla(yyy ~ 1,
           data = data.frame(yyy = c(y, yy)),
           scale = c(rep(1, n), s), 
           control.compute = list(cpo = TRUE, dic = TRUE), 
           family = "gaussian")
print(r$mlik - rr$mlik)

inla.dev.new()
par(mfrow = c(1, 2))
plot(r$internal.marginals.hyperpar[[1]], pch = 19, main = "prec")
lines(rr$internal.marginals.hyperpar[[1]], lwd = 3)
plot(r$marginals.fixed$'(Intercept)', pch = 19, main = "intercept")
lines(rr$marginals.fixed$'(Intercept)', lwd = 3)
