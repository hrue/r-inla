library(sn)
n = 500
x = rnorm(n)
eta = 1/2 + 2*x
w = runif(n, min = 0.5, max = 2)
prec = 1 * w
skewness = 0.25
y = numeric(n)
for(i in 1:n) {
    param = INLA:::inla.sn.reparam(moments = c(eta[i], 1/prec[i], skewness))
    y[i] = rsn(1, xi=param$xi, omega = param$omega, alpha = param$alpha)
}
r = inla(y ~ 1 + x, family = "sn2", scale = w, data = data.frame(y, x, w),  verbose=T)
summary(r)
