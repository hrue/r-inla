n = 1000
alpha = 2
beta = 2
x = runif(n)
eta = 1+beta*x
lambda = exp(eta)
y = rweibull(n, shape= alpha, scale= lambda^(1/-alpha))
event = rep(1,n)
data = list(y=y, event=event, x=x)
formula=inla.surv(y,event)~ x
model=inla(formula, family ="weibull", data=data, verbose=T)
