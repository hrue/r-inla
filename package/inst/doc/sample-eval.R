## ----setup, include=FALSE-------------------------------------------
library(INLA)
inla.setOption(num.threads="1:1")
inla.setOption(smtp="taucs")
inla.setOption(inla.mode="experimental")
if (file.exists("myinit.R")) source("myinit.R")
library(knitr)
set.seed(123)
knitr::opts_chunk$set(echo = TRUE)
opts_chunk$set(size='small'
               , cache=FALSE
               , cache.path='cache/'
               , comment=NA
               , warning=FALSE
               , message=FALSE
               , fig.align='center'
               , fig.path='figures/eval/'
               , fig.pos='H'
               , background='#ffffff'
               , results='hold'
               , eval=TRUE)

## -------------------------------------------------------------------
n <- 100
x <- rnorm(n)
eta <- 1 + x
y <- rnorm(n, mean = eta, sd = 0.1)
r <- inla(y ~ 1 + x,
	data = data.frame(y,x),
    control.compute = list(config=TRUE))

## -------------------------------------------------------------------
samples <- inla.posterior.sample(10000, r)

## -------------------------------------------------------------------
summary(r)

## -------------------------------------------------------------------
fun1 <- function() return (x)

## -------------------------------------------------------------------
eval.fun1 <- inla.posterior.sample.eval(fun1, samples)
str(eval.fun1)

## -------------------------------------------------------------------
hist(eval.fun1[1,], prob=TRUE, n=300)
lines(inla.smarginal(r$marginals.fixed$x), lwd=3)

## -------------------------------------------------------------------
fun2 <- function(x.cov) return (Intercept + x * x.cov)

## -------------------------------------------------------------------
eval.fun2 <- inla.posterior.sample.eval(fun2, samples, x.cov = x)

## -------------------------------------------------------------------
plot(x, apply(eval.fun2, 1, mean))
abline(a=1, b=1) # this is the true curve

## -------------------------------------------------------------------
fun3 <- function(x.cov) return (Predictor - (Intercept + x * x.cov))
eval.fun3 <- inla.posterior.sample.eval(fun3, samples, x.cov = x)
summary(eval.fun3[1,])

## -------------------------------------------------------------------
fun4 <- function() return (theta)
eval.fun4 <- inla.posterior.sample.eval(fun4, samples)
table(eval.fun4[1, ])

## -------------------------------------------------------------------
samples <- inla.posterior.sample(1, r)
fun4 <- function() {
	n <- length(Predictor)
	return (Predictor + rnorm(n, sd = sqrt(1/theta)))
}
eval.fun4 <- inla.posterior.sample.eval(fun4, samples)
plot(x, eval.fun4[,1])

## -------------------------------------------------------------------
r <- inla(y ~ 1 + x,
	family = "sn",
	data = data.frame(y,x),
    control.compute = list(config=TRUE))

## -------------------------------------------------------------------
rownames(r$summary.hyperpar)

## -------------------------------------------------------------------
m <- 100
n <- m^2
## fixed effects
x <- rnorm(n)
xx <- rnorm(n)
## random effects
v <- rnorm(m, sd=0.2)
v.idx <- rep(1:m, each = m)
eta <- 1 + 0.2 * (x + xx) + v
y <- rpois(n, exp(eta))

## -------------------------------------------------------------------
r <- inla(y ~ 1 + x + xx + f(v.idx, model = "iid"),
          data = data.frame(y, x, xx, v.idx),
		  family = "poisson",
          control.compute = list(config = TRUE))
samples <- inla.posterior.sample(10000, r)

## -------------------------------------------------------------------
fun5 <- function(v.index) {
    return (c(Predictor, Predictor - v.idx[v.index]))
}
eval.fun5 <- inla.posterior.sample.eval(fun5, samples, v.index=v.idx)

## -------------------------------------------------------------------
i <- 2
hist(eval.fun5[i,], prob=TRUE, n=300)
lines(density(eval.fun5[n + i,]), col="blue", lwd=3)

## ----eval=FALSE-----------------------------------------------------
# r <- inla(...., control.predictor = list(A=inla.stack.A(...)))

## -------------------------------------------------------------------
n <- 100
m <- 25
## fixed effects
x <- rnorm(n)
xx <- rnorm(n)
## random effects
v <- rnorm(m, sd=0.2)
v.idx <- rep(1:m, each = n %/% m)
eta <- 1 + 0.2 * (x + xx) + v
A <- matrix(rnorm(n^2, sd=sqrt(1/n)), n, n)
eta.star <- A %*% eta
y <- rpois(n, exp(eta.star))

r <- inla(y ~ 1 + x + xx + f(v.idx, model = "iid"),
          data = data.frame(y, x, xx, v.idx),
		  family = "poisson",
	      inla.mode = "experimental",
		  control.predictor = list(A=A),
          control.compute = list(config = TRUE))
samples <- inla.posterior.sample(10000, r)

## -------------------------------------------------------------------
fun6 <- function(v.index) {
    return (c(APredictor, as.numeric(pA %*% (Predictor - v.idx[v.index]))))
}
eval.fun6 <- inla.posterior.sample.eval(fun6, samples, v.index=v.idx)
i <- 2
hist(eval.fun6[i,], prob=TRUE, n=300)
lines(density(eval.fun6[n + i,]), col="blue", lwd=3)

