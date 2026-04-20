library(sn)
n <- 1000
x <- scale(arima.sim(n, model= list(ar = 0.95)))
skew <- 0.2
y <- numeric(n)
for(i in 1:n) {
    variance <- exp(x[i])
    par <- unlist(INLA:::inla.sn.reparam(moments = c(0, variance, skew)))
    y[i] <- rsn(1, dp = par)
}

r = inla(y ~ 1 + f(idx, model="ar1",
                   hyper = list(
                       prec = list(prior = "pc.prec",
                                   param = c(0.5, 0.01)), 
                       rho = list(prior = "pc.cor1",
                                  param = c(0.8, 0.5)))), 
         control.fixed = list(prec.intercept = 1), 
         data = data.frame(y, idx=1:n),
         family = "stochvolsn",
         control.inla = list(cmin = 0, b.strategy="skip"), 
         num.threads = "3:1", 
         verbose = TRUE)
summary(r)
