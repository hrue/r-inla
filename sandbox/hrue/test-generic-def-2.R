inla.my.update(b=T)
source("~/p/inla/google-code/inla/rinla/R/xinla.R")
inla.setOption(inla.arg='-b -v -c -t1')
inla.setOption(inla.call='inla.work')

n = 100
rho1=0.9
rho2=-0.9
y = arima.sim(n, model = list(ar = rho1)) * sqrt(1-rho1^2) + arima.sim(n, model = list(ar = rho2)) * sqrt(1-rho2^2)

idx = 1:n
iidx = 1:n
model1 = inla.R.generic.define(inla.R.generic.ar1.model, n=n, ntheta = 2)
model2 = inla.R.generic.define(inla.R.generic.ar1.model, n=n, ntheta = 2)

formula = y ~ -1 + f(idx,  model="ar1",
        hyper = list(prec = list(
                             initial = -1,
                             param = c(1, 1),
                             prior = "loggamma"), 
                rho = list(
                        initial = 2,
                        param=c(2, 1),
                        prior = "normal"))) + 
    f(iidx,  model="ar1",
        hyper = list(prec = list(
                             initial = -1,
                             param = c(1, 1),
                             prior = "loggamma"), 
                rho = list(
                        initial = -2,
                        param=c(-2, 1),
                        prior = "normal")))

c.family =  list(hyper = list(prec = list(initial = 6, fixed=TRUE)))

rr = inla(formula,
        data = data.frame(y, idx, iidx),
        debug=FALSE,
        verbose=TRUE, 
        control.family = c.family
        )

formula = y ~ -1 + f(idx,  model="rgeneric", R.generic = model1)
formula = y ~ -1 + f(idx,  model="rgeneric", R.generic = model1) + f(iidx,  model="rgeneric", R.generic = model2)
r = xinla(formula,
        data = data.frame(y, idx, iidx),
        debug=FALSE,
        verbose=TRUE, 
        control.family = c.family)
