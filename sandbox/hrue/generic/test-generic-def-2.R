inla.my.update(b=T)
source("~/p/inla/google-code/inla/rinla/R/rgeneric.R")
inla.setOption(inla.arg='-b -v')
inla.setOption(inla.call='inla.work')

n = 100
rho1=0.99
rho2=-0.99
y = arima.sim(n, model = list(ar = rho1)) * sqrt(1-rho1^2) + 
    arima.sim(n, model = list(ar = rho2)) * sqrt(1-rho1^2)

idx = 1:n
iidx = 1:n
model = inla.rgeneric.define(inla.rgeneric.ar1.model, n=n, ntheta = 2)
mmodel = inla.rgeneric.define(inla.rgeneric.ar1.model, n=n, ntheta = 2)

formula = y ~ -1 + f(idx,  model="ar1",
        hyper = list(prec = list(
                             initial = 1,
                             param = c(1, 1),
                             prior = "loggamma"), 
                rho = list(
                        initial = 1,
                        param=c(0, 1),
                        prior = "normal"))) + 
    f(iidx,  model="ar1",
        hyper = list(prec = list(
                             initial = 1,
                             param = c(1, 1),
                             prior = "loggamma"), 
                rho = list(
                        initial = -1,
                        param=c(0, 1),
                        prior = "normal")))

c.family =  list(hyper = list(prec = list(initial = 10, fixed=TRUE)))

rr = inla(formula,
        data = data.frame(y, idx, iidx),
        debug=FALSE,
        verbose=TRUE, 
        control.family = c.family
        )

formula = y ~ -1 + f(idx,  model="rgeneric", rgeneric = model) + f(iidx,  model="rgeneric", rgeneric = mmodel)

r = inla(formula,
        data = data.frame(y, idx, iidx),
        debug=FALSE,
        verbose=TRUE, 
        control.family = c.family)

