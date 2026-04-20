## examples of varying regression

## example 1: time-series
n = 300
z = arima.sim(n,  model = list(ar = 0.9))
z = scale(z)
beta = 1
eta = beta * z
time = 1:n
y = eta + rnorm(n, sd = 0.1)

r = (inla(y ~ 1 + z, data = data.frame(y, z, time)))
r2 = (inla(y ~ 1 + z +
               f(time, z,
                 model = "rw1",
                 constr = TRUE, 
                 scale.model = TRUE, 
                 hyper = list(
                     prec = list(
                         prior = "pc.prec",
                         param = c(1, 0.01)))),
           control.family = list(
               hyper = list(
                   prec = list(
                       prior = "pc.prec",
                       param = c(1, 0.01)))), 
           data = data.frame(y, z, time)))
r3 = (inla(y ~ 1 + z +
                f(time, z,
                  model = "rw2",
                  constr = TRUE, 
                  extraconstr = list(
                      A = matrix(scale(1:n), 1, n),
                      e = 0), 
                  scale.model=TRUE, 
                  hyper = list(
                      prec = list(
                          prior = "pc.prec",
                          param = c(1, 0.001)))),
            control.family = list(
                hyper = list(
                    prec = list(
                        prior = "pc.prec",
                        param = c(1, 0.01)))), 
            data = data.frame(y, z, time)))
r4 = (inla(y ~ 1 + z +
               f(time, z,
                 model = "ar1",
                 constr = TRUE, 
                 hyper = list(
                     prec = list(
                         prior = "pc.prec",
                         param = c(1, 0.001)), 
                     rho = list(
                         prior = "pc.cor1",
                         param = c(0.5, 0.9)))), 
           control.family = list(
               hyper = list(
                   prec = list(
                       prior = "pc.prec",
                       param = c(1, 0.01)))), 
           data = data.frame(y, z, time)))
