## ----setup, include=FALSE-------------------------------------------
knitr::opts_chunk$set(echo = TRUE, cache = FALSE)
knitr::opts_chunk$set(fig.path="figures/int-design/")
Design = matrix()
set.seed(123)
library(INLA)
inla.setOption(num.threads="1:1")
inla.setOption(smtp="taucs")
if (file.exists("myinit.R")) source("myinit.R")

## -------------------------------------------------------------------
	opts = control.inla(int.strategy = "user", int.design = Design)

## ----eval=FALSE-----------------------------------------------------
# 	opts = control.inla(int.strategy = "user.std", int.design = Design)

## ----eval=TRUE, echo=TRUE-------------------------------------------
n = 100
rho = 0.9
x = scale(arima.sim(n, model = list(ar = rho)))
y = x + rnorm(n, sd = 0.1)

## -------------------------------------------------------------------
plot(y, xlab = "time", ylab = "value")
lines(x, lwd=2)

## -------------------------------------------------------------------
rho.0 = rho
to.theta = inla.models()$latent$ar1$hyper$theta2$to.theta
rho.0.internal = to.theta(rho.0)

r = inla(y ~ -1 + f(time, model="ar1",
	hyper = list(
		theta1 = list(prior = "loggamma",
		              param = c(1,1)),
	    theta2 = list(initial = rho.0.internal,
		              fixed=TRUE))),
    control.inla = list(int.strategy = "grid"),
	data = data.frame(y, time = 1:n))

sd.0 = r$summary.random$time[1,"sd"]
print(sd.0)

## -------------------------------------------------------------------
summary(r)
nm <- nrow(r$summary.hyperpar)
Design = as.matrix(cbind(r$joint.hyper[, seq_len(nm)], 1))
head(Design)

## -------------------------------------------------------------------
h.rho = 0.01
rho.1.internal = to.theta(rho.0 + h.rho)
rr = inla(y ~ -1 + f(time, model="ar1",
	hyper = list(
		theta1 = list(prior = "loggamma",
		              param = c(1,1)),
	    theta2 = list(initial = rho.1.internal,
		              fixed=TRUE))),
	control.mode = list(result = r, restart=FALSE),	
    data = data.frame(y, time = 1:n),
	control.inla = list(
		int.strategy = "user",
		int.design = Design))
sd.1 = rr$summary.random$time[1,"sd"]
print(sd.1)

## -------------------------------------------------------------------
deriv.1 = (sd.1 - sd.0) / h.rho
print(deriv.1)

## ----eval=FALSE-----------------------------------------------------
# control.inla = list(int.stategy = "user.expert")

## -------------------------------------------------------------------
n = 50
x = rnorm(n)
y = 1 + x + rnorm(n, sd = 0.2)
param = c(1, 0.04)
dz = 0.1
r.std = inla(y ~ 1 + x, data = data.frame(y, x),
             control.inla = list(int.strategy = "grid",
                                 dz = dz,
                                 diff.logdens = 8), 
             control.family = list(
                 hyper = list(
                     prec = list(
                         prior = "loggamma",
                         param = param))))

s = r.std$internal.summary.hyperpar[1,"sd"]
m = r.std$internal.summary.hyperpar[1,"mean"]
theta = m + s*seq(-4, 4, by = dz)
weight = dnorm(theta,  mean = m, sd = s)

r = rep(list(list()), length(theta))
for(k in seq_along(r)) {
    r[[k]] = inla(y ~ 1 + x,
                  control.family = list(
                      hyper = list(
                          prec = list(
                              initial = theta[k],
                              fixed=TRUE))),
                  data = data.frame(y, x))
}
r.merge = inla.merge(r, prob = weight)

r.design = inla(y ~ 1 + x,
                data = data.frame(y, x),
                control.family = list(
                    hyper = list(
                        prec = list(
                            ## the prior here does not really matter, as we will override
                            ## it with the user.expert in any case.
                            prior = "pc.prec",
                            param = c(1, 0.01)))), 
                control.inla = list(int.strategy = "user.expert",
                                    int.design = cbind(theta, weight)))
for(k in 1:2) {
    plot(inla.smarginal(r.std$marginals.fixed[[k]]),
         lwd = 2, lty = 1, type = "l", 
         xlim = inla.qmarginal(c(0.0001, 0.9999), r.std$marginals.fixed[[k]]))
    lines(inla.smarginal(r.design$marginals.fixed[[k]]),
          lwd = 2, col = "blue", lty = 1)
    lines(inla.smarginal(r.merge$marginals.fixed[[k]]),
          lwd = 2, col = "yellow", lty = 1)
}

