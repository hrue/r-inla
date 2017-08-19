rqloglogistic = function(eta, s, q=0.5)
{
    qs = q/(1-q)
    u = runif(length(eta))
    x = (u/((1-u)*qs))^s * exp(eta)
    return (x)
}

n = 30
q = .10
s = .1
x = rnorm(n, s=0.2)
eta = 1 + 2*x
y = rqloglogistic(eta=eta, s=s, q=q)
r = inla(y ~ 1 + x,
    data = data.frame(y, x),
    family = "qloglogistic",
    control.family = list(quantile = q))
    
summary(r)
