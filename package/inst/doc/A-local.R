## ----setup, include=FALSE, cache=FALSE----------------------------------------
library(knitr)
knitr::opts_chunk$set(echo = TRUE, warning = FALSE, error = FALSE, message = FALSE, cache =FALSE)
knitr::opts_chunk$set(fig.path="figures/A-local/")
library(INLA)
inla.setOption(num.threads="1:1")
inla.setOption(smtp="taucs")
if (file.exists("myinit.R")) source("myinit.R")

## -----------------------------------------------------------------------------
n <- 100
x <- rnorm(n)
i <- 1:n
im <- c(NA,1:(n-1))
y <- 1 + x + c(0, x[-n]) + rnorm(n)

## -----------------------------------------------------------------------------
r1 <- inla(y ~ 1 + f(i,model="iid") + f(im, copy="i"),
	family="stdnormal",	
	data = data.frame(y,i,im))

## -----------------------------------------------------------------------------
A <- matrix(0,n,n)
A[row(A)-col(A) == 1] <- 1
A <- inla.as.sparse(A)
r2 <- inla(y ~ 1 + f(i, model="iid", A.local=A),
	family="stdnormal",	
	data = list(y=y,i=i,A=A))

## -----------------------------------------------------------------------------
diag(A) <- 1 ## we add also the diagonal
r3 <- inla(y ~ 1 + f(i.NA, model="iid", A.local=A, values=1:n),
	family="stdnormal",	
	data = list(y=y, i.NA=rep(NA,n), A=A, n=n))

## -----------------------------------------------------------------------------
r4 <- inla(y ~ 1 + f(i, w.zero, model="iid", A.local=A),
	family="stdnormal",	
	data = list(y=y, i=i, w.zero=rep(0,n), A=A, n=n))

## -----------------------------------------------------------------------------
as.vector(c(r1$mlik[1,], r2$mlik[1,], r3$mlik[1,], r4$mlik[1,]))

## -----------------------------------------------------------------------------
AA <- rbind(A,A)
r5 <- inla(y ~ 1 + f(i, w.zero, model="iid", A.local=AA),
	family="stdnormal",	
	data = list(y=y, i=i, w.zero=rep(0,n), AA=AA, n=n))
r5$mlik[1,]

## -----------------------------------------------------------------------------
n <- 1000
phi <- 0.9
x <- scale(arima.sim(n=n,  model = list(ar = phi)))
xx <- scale(arima.sim(n=n,  model = list(ar = phi)))

## -----------------------------------------------------------------------------
h <- 10L
nh <- h + 1L
kern <- dnorm(-h:0, sd=h)
kern <- kern / sum(kern)
A <- inla.as.sparse(INLA:::inla.toeplitz(c(kern[nh], rep(0, n-nh), kern[1:h])))

h <- 30L
nh <- h + 1L
kern <- dnorm(-h:0, sd= h)
kern <- kern / sum(kern)
B <- inla.as.sparse(INLA:::inla.toeplitz(c(kern[nh], rep(0, n-nh), kern[1:h])))

## -----------------------------------------------------------------------------
eta <- 5 + as.vector(A %*% x) + as.vector(B %*% xx)
y <- rpois(n, lambda = exp(eta))

## -----------------------------------------------------------------------------
AB <- cbind(A, B)
r <- inla(y ~ 1 + f(idx, w, model="ar1", constr = T, A.local = AB, nrep = 2), 
          data = list(y = y, idx = 1:n, AB = AB, w = rep(0, n)),
          family = "poisson")

## ----plot=T-------------------------------------------------------------------
idx <- 1:n
plot(idx, y, pch=19, ylim = range(pretty(y)), main="Data")

## ----plot=T-------------------------------------------------------------------
plot(idx, eta, col="blue", lwd=2, type="l", ylim = range(pretty(eta)))
lines(idx, r$summary.linear.predictor$mean, col="red", lwd=2)
title("eta.true and E(eta|...)")

## ----plot=T-------------------------------------------------------------------
plot(idx, x, col="blue", lwd=2, type="l", ylim = range(pretty(x)))
lines(idx, r$summary.random$idx$mean[1:n], col="red", lwd=2)
title("AR1 and E(AR1|...)")

## ----plot=T-------------------------------------------------------------------
plot(idx, xx, col="blue", lwd=2, type="l", ylim = range(pretty(x)))
lines(idx, r$summary.random$idx$mean[n + 1:n], col="red", lwd=2)
title("AR1 and E(AR1|...)")

## -----------------------------------------------------------------------------
summary(r)

## -----------------------------------------------------------------------------
n <- 300
m <- 30
Z <- matrix(rnorm(n*m), n, m)
rho <- 0.8
Qz <- toeplitz(rho^(0:(m-1)))
z <- inla.qsample(1, Q=Qz)
eta <- Z %*% z
s <- 0.1
y <- eta + rnorm(n, sd = s)

## -----------------------------------------------------------------------------
r <- inla(y ~ -1 + f(idx, model="z", Z=Z, Cmatrix=Qz, precision=exp(20),
					hyper = list(prec = list(initial = 0, fixed = TRUE))),
        data = list(y=y, idx=1:n, n=n),
        control.family = list(hyper = list(prec = list(initial = log(1/s^2), fixed= TRUE))))

## -----------------------------------------------------------------------------
rr <- inla(y ~ -1 + f(idx, model="generic", A.local=Z, Cmatrix=Qz, values=1:m,
					hyper = list(prec = list(initial = 0, fixed = TRUE))),
        data = list(y=y, idx=rep(NA,n), n=n, m=m),
        control.family = list(hyper = list(prec = list(initial = log(1/s^2), fixed= TRUE))))

## -----------------------------------------------------------------------------
r$mlik - (rr$mlik + 0.5 * determinant(Qz, log=TRUE)$modulus)

## -----------------------------------------------------------------------------
mean(abs(r$summary.random$idx$mean[n + 1:m] - rr$summary.random$idx$mean))
mean(abs(r$summary.random$idx$sd[n + 1:m] / rr$summary.random$idx$sd - 1))

## -----------------------------------------------------------------------------
g <- inla.read.graph(system.file("demodata/germany.graph", package="INLA"))
source(system.file("demodata/Bym-map.R", package="INLA"))
n <- g$n
m <- 4
s <- 0.1
X <- matrix(rnorm(n*m), n, m)
y <- rowSums(X) + rnorm(n, sd = s)
A <- inla.as.sparse(diag(X[, 1]))
for(i in 2:m) {
    A <- cbind(A, inla.as.sparse(diag(X[, i])))
}

r <- inla(y ~ -1 +
              f(idx.na,
                model = "besag",
                scale.model = TRUE,
                constr = FALSE,
                nrep = m,
                graph = g, 
                A.local = A,
                values = 1:n,
                hyper = list(prec = list(prior = "pc.prec",
                                         param = c(0.5, 0.01)))), 
          family = "normal",
          control.family = list(hyper = list(prec = list(initial = log(1/s^2),
                                                         fixed = TRUE))), 
          data = list(n = n, m = m, y = y, X = X, 
			  idx.na = rep(NA, n)))

## -----------------------------------------------------------------------------
beta <- matrix(r$summary.random$idx.na$mean, n,m)
summary(beta)

## ----plot=T-------------------------------------------------------------------
i=1; Bym.map(r$summary.random$idx.na$mean[ (i-1) * n + 1:n])
title(paste("coeff ",  i))

## ----plot=T-------------------------------------------------------------------
i=2; Bym.map(r$summary.random$idx.na$mean[ (i-1) * n + 1:n])
title(paste("coeff ",  i))

## ----plot=T-------------------------------------------------------------------
i=3; Bym.map(r$summary.random$idx.na$mean[ (i-1) * n + 1:n])
title(paste("coeff ",  i))

## ----plot=T-------------------------------------------------------------------
i=4; Bym.map(r$summary.random$idx.na$mean[ (i-1) * n + 1:n])
title(paste("coeff ",  i))

## -----------------------------------------------------------------------------
n <- 10
y <- 1:n
A <- matrix(0, n, n)
A[lower.tri(A, diag = TRUE)] <- 1

## -----------------------------------------------------------------------------
cbind(y, A %*% rep(1,n))

## -----------------------------------------------------------------------------
idx <- 1:n
iidx <- rep(NA, n)
w0 <- rep(0, n)
r <- inla(y ~ -1 +
              f(idx,
                w0, 
                model = "iid",
                values = 1:n, 
                hyper = list(prec = list(initial = -20,
                                         fixed = TRUE))) + 
              f(iidx,
                copy = "idx",
                values = 1:n, 
                A.local = A),
          data = list(idx = idx,
                      iidx = iidx,
                      w0 = w0,
                      A = A,
                      n = n),
          family = "stdnormal")

## -----------------------------------------------------------------------------
round(dig=3, cbind(r$summary.random$idx$mean,
	               r$summary.random$iidx$mean))

## -----------------------------------------------------------------------------
round(dig=3, cbind(y, A %*% r$summary.random$iidx$mean))

## -----------------------------------------------------------------------------
r <- inla(y ~ -1 +
              f(idx,
                w0, 
                model = "iid",
                values = 1:n, 
                hyper = list(prec = list(initial = -20,
                                         fixed = TRUE))) + 
              f(iidx,
                copy = "idx",
                values = 1:n, 
                A.local = A,
                hyper = list(beta = list(initial = 0.5))),
          data = list(idx = idx,
                      iidx = iidx,
                      w0 = w0,
                      A = A,
                      n = n),
          family = "stdnormal")

## -----------------------------------------------------------------------------
round(dig=3, cbind(r$summary.random$idx$mean,
	               r$summary.random$iidx$mean))

