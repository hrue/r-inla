## An example demonstrating two ways to implement the model eta = Z*z,
## where z ~ N(0, tau*Q).

## Simulate data
n = 100
m = 10
Z = matrix(rnorm(n*m), n, m)
rho = 0.8
Qz = toeplitz(rho^(0:(m-1)))
prec.fixed = FALSE  ## the precision parameter for z
z = inla.qsample(1, Q=Qz)
eta = Z %*% z
s = 0.1 ## noise stdev
s.fixed = TRUE
y = eta + rnorm(n, sd = s)

## This is normally not needed at all, but it demonstrate how to set
## the 'high precisions': in the z-model, in the A-part of the linear
## predictor, and in the linear predictor iteself.
precision = exp(15)

## The first approach use the z-model.
r = inla(y ~ -1 + f(idx, model="z", Z=Z,
                    precision=precision, 
                    Cmatrix=Qz,
                    hyper = list(
                            prec = list(
                                    initial = 0,
                                    fixed = prec.fixed,
                                    param = c(1, 1)))), 
        data = list(y=y, idx=1:n),
        control.family = list(
                hyper = list(
                        prec = list(
                                initial = log(1/s^2),
                                fixed=s.fixed))), 
        control.predictor = list(
                compute=TRUE,
                precision=precision,
                initial = log(precision)))

## The second one uses the A-matrix
rr = inla(y ~ -1 + f(idx, model="generic",
                     precision = precision, 
                     Cmatrix=Qz,
                     hyper = list(
                             prec = list(
                                     initial = 0,
                                     fixed = prec.fixed,
                                     param = c(1, 1)))), 
        data = list(y=y, idx=1:m),
        control.family = list(
                hyper = list(
                        prec = list(
                                initial = log(1/s^2),
                                fixed=s.fixed))), 
        control.predictor = list(
                compute=TRUE,
                A=Z,
                precision=precision,
                initial = log(precision)))

## Plot some results
par(mfrow=c(2, 2))
plot(r$summary.linear.predictor$mean[1:n], eta,
     main="z-model: (eta.estimated, eta)")
plot(r$summary.linear.predictor$mean[1:n], eta,
     main="generic-model: (eta.estimated, eta)")
plot(r$internal.marginals.hyperpar[[1]],
     main="Prec.param (both)")
lines(rr$internal.marginals.hyperpar[[1]])

## compare (log) marginal likelihood. recall to add the missing part,
## see inla.doc("generic")
print(r$mlik - (rr$mlik + 0.5*log(det(Qz))))


r = inla.hyperpar(r)
rr = inla.hyperpar(rr)
plot(r$internal.marginals.hyperpar[[1]],
     main="Prec.param (improved, both)")
lines(rr$internal.marginals.hyperpar[[1]])

## compare (log) marginal likelihood. recall to add the missing part,
## see inla.doc("generic")
print(r$mlik - (rr$mlik + 0.5*log(det(Qz))))

