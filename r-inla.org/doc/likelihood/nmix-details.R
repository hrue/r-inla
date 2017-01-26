if (FALSE) {
    ## single observation
    p = runif(1)
    lambda = sample(1:20, 1) + runif(1)
    N = rpois(1,  lambda)
    y = rbinom(1, N, p)
    nmax = floor(y + lambda*(1-p)/0.01)  # estimate of a large 'n'

    L = 0
    for(n in y:nmax) 
        L = L + dpois(n, lambda) * dbinom(y, n, p)
    L.direct = L

    L = dpois(y, lambda) * dbinom(y, y, p)
    ff = lambda * (1-p)
    fac = 1
    for(i in nmax:1) 
        fac = 1 + fac * ff / i
    L.recursive = L * fac

    print(c(N=N, nmax = nmax, p=p, lambda=lambda, y=y,
            L.recursive = L.recursive,  abs.err = abs(L.direct - L.recursive)))
}

if (TRUE) {
    ## replicated observations
    nrep = 5
    p = runif(1)
    lambda = sample(1:20, 1) + runif(1)
    N = rpois(1,  lambda)
    y = rbinom(nrep, N, p)
    ymax = max(y)
    nmax = floor(ymax + lambda*(1-p)/0.01)  # estimate of a large 'n'

    L = 0
    for(nn in ymax:nmax) 
        L = L + dpois(nn, lambda) * prod(dbinom(y, nn, p))
    L.direct = L

    L = dpois(ymax, lambda) * prod(dbinom(y, ymax, p))
    fac = 1
    for(nn in nmax:(ymax+1))
        fac = 1 + fac * lambda * (1-p)^nrep / nn * prod(nn/(nn-y))
    L.recursive = L * fac

    print(c(N=N,
            nmax = nmax,
            p=p,
            lambda=lambda,
            y=y,
            L.direct = L.direct,
            L.recursive = L.recursive,
            abs.err = abs(L.direct - L.recursive)))
}
