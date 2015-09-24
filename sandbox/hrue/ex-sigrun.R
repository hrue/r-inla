library(INLA)

cormat = function(n, pacf=0)
{
    ## return the correlation matrix of dimension 'n' for given partial correlation coefficients
    ## pacf
    acf = inla.ar.pacf2acf(pacf, lag.max = n)
    return (toeplitz(acf))
}

n = 20
p = 5
pacf = runif(p, min = -1, max = 1)

sigma0 = cormat(n, pacf)
sigma1 = cormat(n, c(pacf, runif(1, min = -1, max = 1)))

diag(solve(sigma0) %*% sigma1)
