## these should give the same results
n = 300
x = runif(n)
eta = 1+x
y = exp(rnorm(n, mean = eta,  sd = 1))
data = list(y=y, event=rep(1, n), x=x)
formula = inla.surv(y, event) ~ 1 + x
r=inla(formula, family ="lognormalsurv", data=data)
summary(r)

data = data.frame(y, x)
formula = y  ~ 1 + x
r=inla(formula, family ="lognormal", data=data)
summary(r)
