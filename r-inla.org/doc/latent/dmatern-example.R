library(INLA)
library(mvtnorm)

# 1D example.  locations are 1,  2,  ...,  'n',  with 'nr' replications
range = 10
n = 50
nr = 20
loc = 1:n
var = 1.0
nu = 0.5

S = matrix(0, n, n)
for(i in 1:n) {
    for(j in i:n) {
        d = sqrt((loc[i] - loc[j])^2)
        S[i, j] = var * INLA:::inla.matern.cf(d, range = range,  nu = nu);
        S[j, i] = S[i, j]
    }
}
y = c(t(rmvnorm(nr, sigma = S)))
r1 = inla(y ~ -1 + f(idx, model="dmatern", locations = loc, replicate = re,
                     ## placing the prior at the correct value,  just for
                     ## demonstration
                     hyper = list(range = list(initial = log(range),
                                               param = c(range, 0.5)))), 
          data = data.frame(y, idx = rep(1:n, nr), re = rep(1:nr, each = n)), 
          family = "gaussian",
          ## just this this at some high value
          control.family = list(hyper = list(
                                    prec = list(initial = 12,  fixed=TRUE))))

if (FALSE)
    plot(r1,  plot.random.effect = FALSE)

# 2D example. Simulate 'n' data in a [0, 1]^2 box, with 'nr' replications
range = 0.2
n = 50
nr = 20
loc = matrix(runif(2*n), ncol = 2, nrow = n)
var = 1.0
nu = 0.5

S = matrix(0, n, n)
for(i in 1:n) {
    for(j in i:n) {
        dif = loc[i, ] - loc[j, ]
        d = sqrt(sum(dif^2))
        S[i, j] = var * INLA:::inla.matern.cf(d, range = range,  nu = nu);
        S[j, i] = S[i, j]
    }
}
y = c(t(rmvnorm(nr, sigma = S)))
r2 = inla(y ~ -1 + f(idx, model="dmatern", locations = loc, replicate = re,
                     ## placing the prior at the correct value,  just for
                     ## demonstration
                     hyper = list(range = list(initial = log(range),
                                               param = c(range, 0.5)))), 
          data = data.frame(y, idx = rep(1:n, nr), re = rep(1:nr, each = n)), 
          family = "gaussian",
          ## just this this at some high value
          control.family = list(hyper = list(
                                    prec = list(initial = 12,  fixed=TRUE))))

if (FALSE) 
    plot(r2,  plot.random.effect = FALSE)
