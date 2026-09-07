n = 100
s = 0.1
x = rnorm(n)
y = 1 + x + rnorm(n, sd = s)

model = (inla.rgeneric.define(inla.rgeneric.iid.model, n=n, debug=FALSE))
r2 = (inla(y ~ -1 + f(idx, model=model), 
           data = data.frame(y = y, idx = 1:n),
           control.family = list(
               hyper = list(prec = list(initial = 12, fixed=TRUE)))))

r1 = (inla(y ~ -1 + f(idx, model="iid", param = c(1, 1)), 
           data = data.frame(y = y, idx = 1:n),
           control.family = list(
               hyper = list(prec = list(initial=12, fixed=TRUE)))))

print(r1$mlik - r2$mlik)


