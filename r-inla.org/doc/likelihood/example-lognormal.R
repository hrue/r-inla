n = 1000
x = runif(n)
eta = 1+x
y = exp(rnorm(n, mean = eta,  sd = 1))
event = rep(1,n)
data = list(y=y, event=event, x=x)
formula = inla.surv(y, event) ~ 1 + x
r=inla(formula, family ="lognormal", data=data, verbose=T)
