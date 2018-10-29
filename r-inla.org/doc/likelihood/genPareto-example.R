rgp = function(n, sigma, eta, alpha, xi = 0.001)
{
    if (missing(sigma)) {
        stopifnot(!missing(eta) && !missing(alpha))
        sigma = exp(eta) * xi / ((1.0 - alpha)^(-xi) -1.0)
    }
    return (sigma / xi * (runif(n)^(-xi) -1.0))
}

n = 300
x = runif(n)-0.5
eta = 1+x
alpha = 0.95
xi = 0.3
y = rgp(n, eta = eta, alpha = alpha, xi=xi)

r = inla(y ~ 1+x,
         data = data.frame(y, x), 
         family = "gp",
         control.family = list(control.link = list(quantile = alpha)), 
         control.predictor = list(compute=TRUE),
         verbose=TRUE)

    
rx = range(c(r$summary.fitted.values$mean, exp(eta)))
plot(r$summary.fitted.values$mean, exp(eta),
     xlim = rx, ylim = rx)
abline(a=0,b=1)
