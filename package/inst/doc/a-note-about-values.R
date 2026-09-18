## ----setup, include=FALSE, cache=FALSE------------------------------
library(knitr)
knitr::opts_chunk$set(echo = TRUE, warning = FALSE, error = FALSE, message = FALSE, cache =FALSE)
knitr::opts_chunk$set(fig.path="figures/a-note-on-values/")
library(INLA)
inla.setOption(num.threads="1:1")
inla.setOption(smtp="taucs")
if (file.exists("myinit.R")) source("myinit.R")

## ----eval=FALSE-----------------------------------------------------
# y ~ ... + f(year, model = "ar1", ...)

## -------------------------------------------------------------------
year <- 2000:2025

## -------------------------------------------------------------------
values <- sort(unique(year[!is.na(year)]))
values

map.year.to.t <- function(year, values) {
	a <- numeric(length(year))
	for(i in seq_along(year)) 
		a[i] = which(year[i] == values)
	return (a) 
}
yy <- c(2000, 2020, 2025)
print(cbind(year=yy, t=map.year.to.t(yy, values)))

## -------------------------------------------------------------------
year <- c(2000:2020, 2023:2025)

## -------------------------------------------------------------------
values <- sort(unique(year[!is.na(year)]))
values

## -------------------------------------------------------------------
yy <- c(2000, 2020, 2025)
print(cbind(year = yy, t = map.year.to.t(yy, values)))

## ----eval=FALSE-----------------------------------------------------
# time <- c(1.1, 2.2, 3.3, 4.4, 10.0)

## ----eval=FALSE-----------------------------------------------------
# y ~ ... + f(year, model = "ar1", values = 2000:2025)

## ----eval=FALSE-----------------------------------------------------
# y ~ ... + f(year, model = "ar1", values = 1990:2035)

## -------------------------------------------------------------------
values <- sort(unique(values[!is.na(values)]))

## ----plot=TRUE, out.width = "70%", fig.align = "center"-------------
year <- 2000:2025
n <- length(year)
phi <- 0.9
phi.intern <- inla.models()$latent$ar1$hyper$theta2$to.theta(phi)
x <- scale(arima.sim(n, model=list(ar=phi)))
s <- 0.3
y <- x + rnorm(n, sd=s)
plot(year, y, pch=19)

## ----plot=TRUE, out.width = "70%", fig.align = "center"-------------
r <- inla(y ~ -1 + f(year, model="ar1",
	                 hyper = list(prec = list(initial=0, 
						                      fixed=TRUE),
								  rho = list(initial=phi.intern,
									         fixed=TRUE))),
	    family="gaussian",
		control.family = list(hyper = list(
                                  prec = list(initial=log(1/s^2),
								  fixed = TRUE))),
		data = data.frame(y, year))
plot(year, y, pch=19)
lines(year, r$summary.random$year$mean)

## -------------------------------------------------------------------
## this is correct
y[10] <- NA
r <- inla(y ~ -1 + f(year, model="ar1",
	                 hyper = list(prec = list(initial=0, 
						                      fixed=TRUE),
								  rho = list(initial=phi.intern,
									         fixed=TRUE))),
	    family="gaussian",
		control.family = list(hyper = list(
                                  prec = list(initial=log(1/s^2),
								  fixed = TRUE))),
		data = data.frame(y, year))  ## missing data included

## this is wrong
rr <- inla(y ~ -1 + f(year, model="ar1",
	                 hyper = list(prec = list(initial=0, 
						                      fixed=TRUE),
								  rho = list(initial=phi.intern,
									         fixed=TRUE))),
	    family="gaussian",
		control.family = list(hyper = list(
                                  prec = list(initial=log(1/s^2),
								  fixed = TRUE))),
		data = data.frame(y, year)[-10,]) ## missing data removed
## we compare the results
r$mlik - rr$mlik

## -------------------------------------------------------------------
rrr <- inla(y ~ -1 + f(year, model="ar1", 
	                  values = 2000:2025,  ## values added
	                  hyper = list(prec = list(initial=0, 
						                      fixed=TRUE),
								  rho = list(initial=phi.intern,
									         fixed=TRUE))),
	    family="gaussian",
		control.family = list(hyper = list(
                                  prec = list(initial=log(1/s^2),
								  fixed = TRUE))),
		data = data.frame(y, year)[-10,]) ## missing data removed
## we compare the results
r$mlik - rrr$mlik

## ----plot=TRUE, out.width = "70%", fig.align = "center"-------------
values <- 1990:2035
r <- inla(y ~ -1 + f(year, model="ar1",
	                 values = values,  ## values added
	                 hyper = list(prec = list(initial=0, 
						                      fixed=TRUE),
								  rho = list(initial=phi.intern,
									         fixed=TRUE))),
	    family="gaussian",
		control.family = list(hyper = list(
                                  prec = list(initial=log(1/s^2),
								  fixed = TRUE))),
		data = list(y=y, year=year, values=values))
plot(year, y, pch=19, xlim=range(values))
lines(values, r$summary.random$year$mean)

