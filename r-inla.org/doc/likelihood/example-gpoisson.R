dgpoisson = function(y, mu, phi, p)
{
    a = mu + phi * mu^(p-1.0) * y;
    b = 1. + phi * mu^(p-1.0);
    d = exp(log(mu) + (y-1.0)*log(a) -
            y*log(b) - lfactorial(y) - a/b)
    return (d)
}
rgpoisson = function(n, mu, phi, p)
{
    stopifnot(length(mu) == 1)
    s = sqrt(mu*(1+phi*mu^(p-1))^2)
    f = 20
    low = as.integer(max(0, mu - f*s))
    high = as.integer(mu + f*s)
    prob = dgpoisson(low:high, mu, phi, p)
    y = sample(low:high, n, replace=TRUE,
            prob = prob)

    return (y)
}

n = 1000
phi = 1
p = 1
mu = exp(1 + 5*(1:n)/n)

y = numeric(n)
for(i in 1:n) {
    y[i] = rgpoisson(1, mu[i], phi, p)
}

idx = (1:n)/n
r = inla(y ~ 1 + idx,  data = data.frame(y, idx),
        family = "gpoisson")
