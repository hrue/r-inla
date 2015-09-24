kld = function(sigma0, sigma1)
{
    trace = function(A)
    {
        return (sum(diag(A)))
    }

    logdet = function(A)
    {
        return (as.numeric(determinant(A, log=TRUE)$modulus))
    }

    ## print(trace(solve(sigma0) %*% sigma1))
    
    kl = (0.5 * (trace(solve(sigma0) %*% sigma1)
                 - dim(sigma0)[1]
                 - (logdet(sigma1) - logdet(sigma0))))
    return(pmax(0, kl))
}

dist = function(sigma0, sigma1, normalized = TRUE)
{
    d = sqrt(2*kld(sigma0, sigma1))
    if (normalized) {
        d = d/sqrt(dim(sigma0)[1])
    }
    return (d)
}

cormat = function(n, pacf=0)
{
    ## return the correlation matrix of dimension 'n' for given partial correlation coefficients
    ## pacf
    acf = inla.ar.pacf2acf(pacf, lag.max = n)
    return (toeplitz(acf))
}

dists = function(n, pacf, normalized = TRUE)
{
    g = function(x) 2*exp(x)/(1+exp(x))-1
    
    pacfs = g(seq(-10, 10, len=301))
    if (missing(pacf) || is.null(pacf)) {
        sigma0 = diag(n)
    } else {
        sigma0 = cormat(n, pacf)
    }

    d = c()
    for(pac in pacfs) {
        sigma1 = cormat(n, pacf=c(pacf, pac))
        d = c(d, dist(sigma0, sigma1, normalized))
    }
    return(cbind(pacf=pacfs, d = d))
}

n=20
for(i in 1:10) {

    pacf = runif(3, min = -1, max = 1)
    a = dists(n, pacf)
    if (i == 1) {
        plot(a[, 1], a[, 2],  type="l")
    } else {
        lines(a[, 1],  a[, 2])
    }
}

    
