rpom = function(alpha, eta) 
{
    ## alpha: the cutpoints. eta: the linear predictor
    F = function(x) 1.0/(1+exp(-x))

    ns = length(eta)
    y = numeric(ns)
    nc = length(alpha) + 1

    for(k in 1:ns) {
        p = diff(c(0.0, F(alpha - eta[k]), 1.0))
        y[k] = sample(1:nc, 1, prob = p)
    }
    return (y)
}

n = 300
nsim = 1E5
x = rnorm(n, sd = 0.3)

eta = x
alpha = c(-1, 0, 0.5)
y = rpom(alpha, eta)
prior.alpha = 3 ## parameter in the Dirichlet prior
r = inla(y ~ -1 + x, data = data.frame(y, x, idx = 1:n), family = "pom", 
         control.family = list(hyper = list(theta1 = list(param = prior.alpha))))
summary(r)

## compute the posterior for the cutpoints 
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
m1 = colMeans(sim.alpha)
m2 = colMeans(sim.alpha^2)
print(cbind(truth = alpha, estimate = m1, stdev = sqrt(m2 - m1^2)))

for(k in 1:length(alpha)) {
    d = density(sim.alpha[, k])
    if (k == 1) {
        plot(d, xlim = 1.2*range(c(sim.alpha)),
             ylim = c(0, 1.5 * max(d$y)), type="l", lty=k, lwd=2)
    } else {
        lines(d, xlim = range(c(sim.alpha)), lty = k, lwd=2)
    }
    abline(v = alpha[k], lty=k, lwd=2)
}
