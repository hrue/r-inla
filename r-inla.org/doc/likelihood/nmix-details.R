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

print("Single Poisson")
print(c(N=N, nmax = nmax, p=p, lambda=lambda, y=y,
        L.recursive = L.recursive,  abs.err = abs(L.direct - L.recursive)))

#####################

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

print("Replicated Poisson")
print(c(N=N,
        nmax = nmax,
        p=p,
        lambda=lambda,
        y=y,
        L.direct = L.direct,
        L.recursive = L.recursive,
        abs.err = abs(L.direct - L.recursive)))


#####################################################
## single observation
p = runif(1)
lambda = sample(1:20, 1) + runif(1)
delta = exp(1 + rnorm(1))
N = rnbinom(1,  mu = lambda, size = delta)
y = rbinom(1, N, p)
nmax = floor(y + lambda*(1-p)/0.01)  # estimate of a large 'n'

L = 0
for(n in y:nmax) 
    L = L + dnbinom(n, mu=lambda, size=delta) * dbinom(y, n, p)
L.direct = L

L = dnbinom(y, mu=lambda, size=delta) * dbinom(y, y, p)
q = delta/(delta + lambda)
fac = 1
for(n in nmax:(y+1))
    fac = 1 + fac * (n-1+delta)/(n-y) *(1-q)*(1-p)
L.recursive = L * fac

print("Single nBinom")
print(c(N=N, nmax = nmax, p=p, lambda=lambda, y=y, delta=delta, 
        L.recursive = L.recursive,  abs.err = abs(L.direct - L.recursive)))

#####################################################
## replicated observation
nrep = 5
p = runif(1)
lambda = sample(1:20, 1) + runif(1)
delta = exp(1 + rnorm(1))
N = rnbinom(1,  mu = lambda, size = delta)
y = rbinom(nrep, N, p)
ymax = max(y)
nmax = floor(ymax + lambda*(1-p)/0.01)  # estimate of a large 'n'

L = 0
for(n in ymax:nmax) 
    L = L + dnbinom(n, mu=lambda, size=delta) * prod(dbinom(y, n, p))
L.direct = L

L = dnbinom(ymax, mu=lambda, size=delta) * prod(dbinom(y, ymax, p))
q = delta/(delta + lambda)
fac = 1
for(n in nmax:(ymax+1))
    fac = 1 + fac * (n-1+delta)/n * (1-q) * (1-p)^nrep * prod(n/(n-y))
L.recursive = L * fac

print("Replicated nBinom")
print(c(N=N, nmax = nmax, p=p, lambda=lambda, ymax=ymax, delta=delta, 
        L.recursive = L.recursive,  abs.err = abs(L.direct - L.recursive)))
