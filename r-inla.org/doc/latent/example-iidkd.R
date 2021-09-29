library(mvtnorm)

n <- 300
m <- 4
N <- m*n
rho <- 0.8

Sigma <- matrix(NA, m, m)
diag(Sigma) <- (1/(1:m))^2
for(i in 1:m) {
    for (j in 1:m) {
        if (i != j) {
            Sigma[i, j] <- rho^abs(i-j) * sqrt(Sigma[i, i] * Sigma[j, j])
        }
    }
}

y <- c()
yy <- rmvnorm(n, sigma = Sigma)
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
         verbose = FALSE)

## this is how the internal parameters are defined
L <- t(chol(solve(Sigma)))
diag(L) <- log(diag(L))
LL <- t(chol(solve(cov(yy))))
diag(LL) <- log(diag(LL))

## compare the estimated (internal) parameters with MLE and the truth
round(dig = 3, cbind(true = c(diag(L), L[lower.tri(L)]),
                     mle = c(diag(LL), LL[lower.tri(LL)]),
                     inla = r$mode$theta))

## this gives a list of sampled matrices (stdev's and correlations)
xx <- inla.iidkd.sample(10^4, r, "i")
## compute the mean
qq <- matrix(rowMeans(matrix(unlist(xx), nrow = m^2)), m, m)

iSigma <- 1/sqrt(diag(Sigma))
Cor <- diag(iSigma) %*% Sigma %*% diag(iSigma)
round(dig = 3, cbind(inla = c(diag(qq), qq[lower.tri(qq)]), 
                     true = c(sqrt(diag(Sigma)), Cor[lower.tri(Cor)])))
