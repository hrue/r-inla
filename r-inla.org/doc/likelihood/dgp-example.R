F = function(y, sigma, xi) 1.0 - (1.0 + xi * (y+1)/sigma)^(-1/xi)
f = function(y, sigma, xi) F(y, sigma, xi) - F(y-1, sigma, xi)

rdgp = function(n, sigma, eta, alpha, xi = 0.001)
{
    if (missing(sigma)) {
        stopifnot(!missing(eta) && !missing(alpha))
        stopifnot(length(eta) == 1)
        sigma = exp(eta) * xi / ((1.0 - alpha)^(-xi) -1.0)
    }
    stopifnot(length(sigma) == 1)
    eps = 1e-7
    y.max = ceiling((eps^(-xi) -1) * sigma/xi)
    return (sample(0:y.max, prob = f(0:y.max, sigma, xi),
                   size=n, replace=TRUE))
}

n = 300
x = runif(n)-0.5
eta = 5+x
alpha = 0.95
xi = 0.3
y = numeric(n)
for(i in 1:n) {
    y[i] = rdgp(1, eta = eta[i], alpha = alpha, xi=xi)
}

r = inla(y ~ 1+x,
         data = data.frame(y, x), 
         family = "dgp",
         control.family = list(
             control.link = list(quantile = alpha),
             hyper = list(tail = list(
                              prior = "pc.gevtail",
                              param = c(7, 0.0, 0.5)))), 
         control.predictor = list(compute=TRUE),
         verbose=TRUE)

summary(r)
plot(r, plot.prior=TRUE)

dev.new()
plot(cbind(r$summary.fitted.values$mean, exp(eta)))
abline(a=0, b=1)
