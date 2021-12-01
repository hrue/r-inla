rpom = function(alpha, eta) 
{
    ## alpha: the cutpoints. eta: the linear predictor
    F = function(a, eta) pnorm(a-eta)

    ns = length(eta)
    y = numeric(ns)
    nc = length(alpha) + 1

    for(k in 1:ns) {
        p = diff(c(0.0, F(alpha, eta[k]), 1.0))
        y[k] = sample(1:nc, 1, prob = p)
    }
    return (y)
}

library(INLA)
##inla.setOption(inla.mode = "experimental")

n = 3000
nsim = 1E5
x = rnorm(n, sd = 0.3)

eta = x + rnorm(n, sd = 0.2)
alpha = c(-1, 0, 0.5, 1.25)
y = rpom(alpha, eta)
xx = inla.group(x)

r = inla(y ~ -1 +
             f(xx, model="rw2", scale.model=TRUE, 
               hyper = list(prec = list(prior = "pc.prec",
                                        param = c(0.5, 0.01)))) +
             f(idx, model="iid", 
               hyper = list(prec = list(prior = "pc.prec",
                                        param = c(0.5, 0.01)))), 
         data = data.frame(y, x, idx = 1:n, xx),
         family = "pom",
         control.family = list(control.pom = list(cdf = "probit")), 
         ##control.inla = list(cmin=0), 
         control.fixed = list(prec = 1, prec.intercept = 1))

theta = inla.hyperpar.sample(nsim, r, intern=TRUE)
nms = paste(paste0("theta", 1:length(alpha)), "for POM")
sim.alpha = matrix(NA, dim(theta)[1], length(alpha))
for(k in 1:length(alpha)) {
    if (k == 1) {
        sim.alpha[, k] = theta[, nms[1]]
    } else {
        sim.alpha[, k] = sim.alpha[, k-1] + exp(theta[, nms[k]])
    }
}
colnames(sim.alpha) = paste0("alpha", 1:length(alpha))
cbind(truth = alpha,
      estimated = colMeans(sim.alpha),
      "stdev(estimate)" = sqrt(colMeans(sim.alpha^2) - colMeans(sim.alpha)^2))

for(k in 1:length(alpha)) {
    d = density(sim.alpha[, k])
    if (k == 1) {
        plot(d, xlim = 1.2*range(c(sim.alpha)), ylim = c(0, 1.5 * max(d$y)), type="l", lty=k, lwd=2)
    } else {
        lines(d, xlim = range(c(sim.alpha)), lty = k, lwd=2)
    }
    abline(v = alpha[k], lty=k, lwd=2)
}
