n = 100
rho=0.9
y = arima.sim(n, model = list(ar = rho)) * sqrt(1-rho^2)
idx = 1:n
model = inla.rgeneric.define(inla.rgeneric.ar1.model, n=n, ntheta = 2)
formula = y ~ -1 + f(idx,  model="rgeneric", rgeneric = model)
r = inla(formula, data = data.frame(y, idx), family = "gaussian")
