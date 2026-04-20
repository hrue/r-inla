n = 1000
x = rnorm(n, sd = 0.5)
lambda = exp(1+x)
y = rexp(n, rate=lambda)
event = rep(1,n)
data = list(y=y, event=event, x=x)
formula = inla.surv(y,event) ~ x 
model = inla(formula, family ="exponential.surv", data=data)
summary(model)

formula = y ~ x 
model = inla(formula, family ="exponential", data=data)
summary(model)


n = 1000
x = rnorm(n, sd = 0.5)
lambda = 1/exp(1+x)
yy = rexp(n, rate=lambda)
ys <-  rexp(n, rate = exp(1))
y <- pmin(yy, ys)
event <- as.numeric(ys > yy)

data = list(y=y, event=event, x=x)
summary(inla(inla.surv(y,event) ~ x,
             family ="exponential.surv",
             control.family = list(link = list(model = "neglog")), 
             data=data))

library(survival)
summary(survreg(Surv(y, event) ~ x, data=data, dist="exponential"))
