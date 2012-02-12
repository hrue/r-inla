## overdispersion parameter in the betabinomial
rho = 0.7

n = 1000
z = rnorm(n, sd=0.2)
Ntrials = sample(1:10, n, replace=TRUE)
eta = 1 + z
p.eta = exp(eta)/(1+exp(eta))
a = p.eta * (1-rho)/rho
b = (p.eta * rho - p.eta - rho + 1)/rho
p = rbeta(n, a, b)
y = rbinom(n, Ntrials, p)

formula = y ~ 1 + z
data = data.frame(y, z)
r = inla(formula, data = data,
        family = "betabinomial", Ntrials=Ntrials)
summary(r)

