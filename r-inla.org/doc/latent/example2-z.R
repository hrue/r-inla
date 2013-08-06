## This example demonstrate how to use the z-model with intrinsic
## models. The z-model must be proper, which we have to mimic if we
## are using an intrinsic model

## Simulate some data
n = 100
idx = 1:n
x = sin(idx / n * 4 * pi)
s = 0.1
y = x + rnorm(n, sd=s)

## Parameters for the loggamma prior
prior = c(1, 0.001)

## A small constant we add to the diagonal to prevent the model to be
## intrinsic. 
d = 1e-8

## RW1
r = inla(y ~ -1 + f(idx, model="rw1", param=prior,
                    constr=TRUE, diagonal=d),
        data = data.frame(y, idx),
        control.family = list(
                hyper = list(
                        prec = list(
                                initial = log(1/s^2),
                                fixed = TRUE))))
C = toeplitz(c(2, -1, rep(0, n-2)))
C[1, 1] = C[n, n] = 1
## We must add the extra diagonal contribution here, as otherwise it
## applies to the hole model (v, z)
diag(C) = diag(C) + d
Z = diag(n)
rr = inla(y ~ -1 + f(idx, model="z", Z=Z,
                     Cmatrix = C, constr=TRUE, param=prior),
        data = data.frame(y, idx), 
        control.family = list(
                hyper = list(
                        prec = list(
                                initial = log(1/s^2),
                                fixed = TRUE))))
##
par(mfrow=c(2, 2))
plot(idx, r$summary.random$idx$mean)
lines(idx, rr$summary.random$idx$mean[1:n])
title("rw1: idx")
plot(r$internal.marginals.hyperpar[[1]])
lines(rr$internal.marginals.hyperpar[[1]])
title("rw1: log.prec")      

## RW2
r = inla(y ~ -1 + f(idx, model="rw2",
        ## we cannot define the rankdef for the z-model, but it will
        ## be set to 1 as constr=TRUE. so we're using the same here,
        ## even though the rankdef is 0, since we added 'd' on the
        ## diagonal.
        rankdef = 1, 
        param=prior, constr=TRUE, diagonal = d),
        data = data.frame(y, idx),
        control.family = list(
                hyper = list(
                        prec = list(
                                initial = log(1/s^2),
                                fixed = TRUE))))

C = toeplitz(c(2, -1, rep(0, n-3), -1))
C = C[-c(1, n), ]
C = t(C) %*% C
## We must add the extra diagonal contribution here, as otherwise it
## applies to the hole model (v, z)
diag(C) = diag(C) + d
Z = diag(n)
rr = inla(y ~ -1 + f(idx, model="z", Z=Z, Cmatrix = C,
        constr=TRUE, param=prior),
        data = data.frame(y, idx), 
        control.family = list(
                hyper = list(
                        prec = list(
                                initial = log(1/s^2),
                                fixed = TRUE))))
##
plot(idx, r$summary.random$idx$mean)
lines(idx, rr$summary.random$idx$mean[1:n])
title("rw2: idx")
plot(r$internal.marginals.hyperpar[[1]])
lines(rr$internal.marginals.hyperpar[[1]])
title("rw2: log.prec")
      
