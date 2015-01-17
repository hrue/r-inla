G = function(Alpha, Beta) {
    return (pgamma(Beta, shape=Alpha, rate=1))
}

n = 1000
x = rnorm(n)
eta = 1 + x
alpha = 1.5
T = 1
m = 100
y = numeric(n)
prob = numeric(m+1)

for(i in 1:n) {

    ## compute the discrete probability distribution and
    ## then sample from it
    for(j in 1:m) {
        yy = j-1
        beta = alpha * exp(eta[i])
        prob[j] = (G(yy*alpha, beta*T) -
                   G((yy+1)*alpha, beta*T))
    }
    y[i] = sample(0:m, size=1, prob = prob)
}

r = (inla(y ~ 1 + x,
          data = data.frame(y, x),
          family = "gammacount"))
summary(r)
