## simulate an OU-process and estimate its parameters back.
phi = -log(0.95)
sigma = 1
marg.prec = 2*phi/sigma^2
n = 1000
locations = cumsum(sample(c(1, 2, 5, 20),n, replace=TRUE))

## do it sequentially and slow (for clarity)
x = numeric(n)
x[1] = rnorm(1, mean=0, sd = sqrt(1/marg.prec))
for(i in 2:n) {
    delta = locations[i] - locations[i-1]
    x[i] = x[i-1] * exp(-phi * delta) +
        rnorm(1, mean=0,  sd = sqrt(1/marg.prec * (1-exp(-2*phi*delta))))
}

## observe it with a little noise
y = 1 + x + rnorm(n, sd= 0.01)
plot(locations, x, type="l")

formula = y ~  1 + f(locations, model="ou", values=locations)
r = inla(formula, data = data.frame(y, locations))
summary(r)
