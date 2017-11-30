rkumar = function(n, eta, phi, q=0.5)
{
    kappa = eta
    beta = log(1-q)/log(1-exp(-phi))
    alpha = log(1- (1-q)^(1/beta)) / log(kappa)
    u = runif(n)
    y = (1-u^(1/beta))^(1/alpha)
    return (y)
}

n = 100
q = 0.5
phi = 1
x = rnorm(n, sd = 1)
eta = inla.link.invlogit(1 + x)
y = rkumar(n, eta, phi, q)
r = inla(y ~ 1 + x,
    data = data.frame(y, x),
    family = "qkumar",
    control.family = list(control.link=list(quantile = q)))
summary(r)


