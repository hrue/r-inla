bty = "l"
lwd = 2L

prior.kappa.q = function(q,  u,  alpha,  theta)
{
    if (missing(theta)) {
        theta = -log(alpha)/u
    }
    return ((-theta/log(q))^2)
}


## I want to display the Gamma(1, a) prior, with the same median as
## the new one with (u=1, alpha=0.01).

rx = function(n, a)
{
    prec = rgamma(n, shape = 1, rate = a)
    return (rnorm(n, mean = 0, sd = 1/sqrt(prec)))
}


if (FALSE) {
    ## match the median
    a = optimize(
            function(a) {
                (qgamma(0.5, shape=1, rate = a) - inla.pc.qprec(0.5, 1, 0.01))^2
            }, interval = c(0.001,  1))$minimum
    print(paste("Gamma(1,", a, ") gives the same median as (1,0.01)"))
} else {
    ## match the marginal stdev
    a = optimize(
            function(a) {
                set.seed(1234)
                return ((sd(rx(1000000, a))-0.31)^2)
            }, interval = c(0.01, 0.00001))$minimum
    print(paste("Gamma(1,", a, ") gives the same sd = 0.31"))
}

kappa = exp(seq(-5, 15, len=10000))
prior.new = inla.pc.dprec(kappa, 1, 0.01)
prior.gamma = dgamma(kappa, shape=1, rate = a)

inla.dev.new()
plot(kappa, prior.new,
     xlim = c(0,  500),
     xlab = "Precision",
     ylab = "Density", 
     type = "l", lwd = lwd, bty = bty, lty = 2)
lines(kappa, prior.gamma, lwd = lwd, lty = 1)
dev.print(postscript, file = "precision-priors.ps")

d = 1/sqrt(kappa)
prior.d = prior.new * 2*kappa^(3/2)
prior.gamma.d = prior.gamma * 2 * kappa^(3/2)
inla.dev.new()
plot(d, prior.d,
     xlim = c(0,  2),
     ylim = c(0, 9), 
     xlab = "Distance",
     ylab = "Density", 
     type = "l", lwd = lwd, bty = bty, lty = 2)
lines(d, prior.gamma.d, lwd = lwd, lty = 1)
dev.print(postscript, file = "precision-priors-d.ps")


system("which psfix && psfix prec*.ps")


