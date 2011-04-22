rlogistic = function(n, mean = 0, sd = 1)
{
    p = runif(n)
    A = pi/sqrt(3)
    tauA = A/sd^2
    return ((tauA * mean - log((1-p)/p))/tauA)
}

n = 1000
z = rnorm(n, sd=0.1)
eta = 1 + z
y = rlogistic(n, mean = eta,  sd = 1)

r = inla(y ~ 1 + z,  data = data.frame(y, z), family = "logistic",
        control.compute = list(cpo=TRUE))
