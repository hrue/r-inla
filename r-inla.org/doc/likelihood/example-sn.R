library(sn)
set.seed(246)
n = 300
x = rnorm(n, sd = 1)
eta = 1+x
skewness = 0.25
y = numeric(n)
prec <- 100
for(i in 1:n) {
    ## map moments to sn-parameters c(xi, omega, alpha)
    param = INLA:::inla.sn.reparam(moments = c(eta[i], 1/prec, skewness))
    y[i] = rsn(1, xi=param$xi, omega = param$omega, alpha = param$alpha)
}

r = inla(y ~ 1+x,
         family = "sn",
         data = data.frame(y, x),
         control.family = list(
             hyper = list(prec = list(
                              prior = "pc.prec",
                              param = c(3, 0.01)))))
summary(r)
