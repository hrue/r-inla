n = 1000
alpha = 2
beta = 2
rho = 0.5

x = runif(n)
censorTime = runif(n,0,2)
eta = 1+beta*x
lambda = exp(eta)
y = rweibull(n, shape= alpha, scale= lambda^(1/-alpha))
z = rbinom(n,size=1, prob=rho)

censoredEvent = (y > censorTime) | z
yObs = y
yObs[censoredEvent] = censorTime[censoredEvent]
event = as.numeric(!censoredEvent)
data = list(y=inla.surv(yObs, event), x=x)

model=inla(
  y ~ x, 
  family ="weibullcure", 
  data=data,
  control.family = list(hyper=list(
      'log alpha' = list(
          prior='loggamma', param=c(1,1)),
      'logit probability' = list(
          prior='logitbeta', param=c(1,1)))))



summary(model)
