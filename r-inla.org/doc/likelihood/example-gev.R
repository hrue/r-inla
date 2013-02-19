rgev = function(n=1, xi = 0, mu = 0.0, sd = 1.0) {
    u = runif(n)
    if (xi == 0) {
        x = -log(-log(u))
    } else {
        x = ((-log(u))^(-xi) - 1.0)/xi
    }
    return (x*sd + mu)
}

n = 300
z = rnorm(n)
sd.y = 0.5
xi = 0.2
y = 1+z + rgev(n, xi=xi, sd = sd.y)

r = inla(y ~ 1 + z, data = data.frame(y, z), family = "gev",
        control.family = list(gev.scale.xi = 0.01))
summary(r)

