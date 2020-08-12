p = 4
nsim = 1000
lambda = 5 + (1:(1+p))^2
X = matrix(NA, nsim, p)

for(i in 1:nsim) {
    for(j in 1:p) {
        X[i, j] = inla.pc.rcor0(1, lambda = lambda[j])
    }
}

pairs(t(apply(X, 1, inla.ar.pacf2phi)))

Y = apply(X, 1, inla.ar.pacf2acf, lag.max = 20)
matplot(Y, ylim = c(-1, 1))
