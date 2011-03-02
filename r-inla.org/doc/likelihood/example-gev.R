rgev = function(n=1, xi = 0, mu = 0.0, sd = 1.0) {
    u = runif(n)
    if (xi == 0) {
        x = -log(-log(u))
    } else {
        x = ((-log(u))^(-xi) - 1.0)/xi
    }
    return (x*sd + mu)
}

n = 100
z = rnorm(n)
sd.y = 0.5
xi = 0
y = 1+z + rgev(n, xi=xi, sd = sd.y)

formula = y ~ 1 + f(inla.group(z), model="rw1")
data = data.frame(y,z)

r = inla(formula, data = data, family = "gev",
        control.data = list(gev.scale.xi = 0.01,
                ## just to show how to set an initial value
                hyper = list(prec=list(initial=2))))

