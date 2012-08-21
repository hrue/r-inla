dgpoisson = function(y, mu, phi, p)
{
    a = mu + phi * mu^(p-1.0) * y;
    b = 1. + phi * mu^(p-1.0);
    d = exp(log(mu) + (y-1.0)*log(a) -
            y*log(b) - lfactorial(y) -a/b)
    return (d)
}
rgpoisson = function(n, mu, phi, p)
{
    stopifnot(length(mu) == 1)
    s = sqrt(mu*(1+phi*mu^(p-1)))
    f = 10
    low = as.integer(max(0, mu - f*s))
    high = as.integer(mu + f*s)

    prob = dgpoisson(low:high, mu, phi, p)
    y = sample(low:high, n, replace=TRUE,
            prob = prob)

    return (y)
}

n = 1000
phi = 0
mu = 5
p = 1
y = rgpoisson(n, mu, phi, p)

r = inla(y ~ 1,  data = data.frame(y),
        family = "gpoisson",
        control.family = list(gpoisson.p = p))
