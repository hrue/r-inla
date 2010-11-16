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

model=inla(formula, family ="coxph", data=data, verbose=T,
        control.hazard=list(model="rw1", n.intervals=20))
