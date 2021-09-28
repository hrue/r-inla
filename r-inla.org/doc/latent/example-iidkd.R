n <- 300
m <- 4
N <- m*n
rho <- 0.8

Sigma <- matrix(NA, m, m)
diag(Sigma) <- (1/(1:m))^2
for(i in 1:m) {
    if (i+1 <= m) {
        for (j in (i+1):m) {
            Sigma[i, j] <- Sigma[j, i] <- rho * sqrt(Sigma[i, i]*Sigma[j, j])
        }
    }
}

library(mvtnorm)
yy <- rmvnorm(n, sigma = Sigma)
y <- c()
for(i in 1:m) {
    y <- c(y, yy[, i])
}

r <- inla(y ~ f(i, model = "iidkd", order = m, n=N,
                ## set parameters using 'theta1'.
                ## these are the default parameters. 
                hyper =  list(theta1 = list(
                                  param = c(100, rep(1, m), rep(0, m*(m-1)/2))))), 
         data = data.frame(i = 1:N, y),
         ## fix precision as we have exact observations
         control.family = list(hyper = list(
                                   prec = list(initial = 15, fixed = TRUE))), 
         verbose = TRUE)

## this is how the internal parameters are defined
L <- t(chol(solve(Sigma)))
diag(L) <- log(diag(L))
LL <- t(chol(solve(cov(yy))))
diag(LL) <- log(diag(LL))

## compare the estimated (internal) parameters with MLE and the truth
round(dig = 3, cbind(true = c(diag(L), L[lower.tri(L)]),
                     mle = c(diag(LL), LL[lower.tri(LL)]),
                     inla = r$mode$theta))

## this is how to compute stdev and correlations from the internal parametes.
## these ones are more interpretable. one have to know where the parameters
## are in the list though, but here it is easy...

convert.internal <- function(theta, dim, offset = 0)
{
    ntheta.off.diag <- dim*(dim-1)/2L
    L <- matrix(0, dim, dim)
    diag(L) <- exp(theta[offset + 1:dim])
    L[lower.tri(L)] <- theta[offset + dim + 1:ntheta.off.diag]
    S <- solve(L %*% t(L))
    iSigma <- 1/sqrt(diag(S))
    Cor <- diag(iSigma) %*% S %*% diag(iSigma)
    return(c(sqrt(diag(S)), Cor[lower.tri(Cor)]))
}

xx <- inla.hyperpar.sample(10000, r)
qq <- rowMeans(apply(xx, 1, convert.internal, dim = m))
iSigma <- 1/sqrt(diag(Sigma))
Cor <- diag(iSigma) %*% Sigma %*% diag(iSigma)
round(dig = 3, cbind(inla = c(qq),
                     true = c(sqrt(diag(Sigma)), Cor[lower.tri(Cor)])))
