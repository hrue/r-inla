n = 300
rho = 0.2
z = rnorm(n, sd=0.2)
Ntrials = sample(10:20, n, replace=TRUE)
eta = 1 + z
p.eta = exp(eta)/(1+exp(eta))
a = p.eta * (1-rho)/rho
b = (p.eta * rho - p.eta - rho + 1)/rho
p = rbeta(n, a, b)
y = rbinom(n, Ntrials, p)

formula = y ~ 1 + z
data = data.frame(y, z)
## exact
r = inla(formula, data = data,
        family = "betabinomial", Ntrials=Ntrials)
summary(r)
## approximate
ra = inla(formula, data = data,
        family = "betabinomialna", Ntrials=Ntrials)
summary(ra)



## exact simulation from the approximate model
n = 1000
rho = 0.1
z = rnorm(n, sd=0.4)
Ntrials = sample(1:20, n, replace=TRUE)
eta = 1 + z
p = exp(eta)/(1+exp(eta))
s = runif(n)
m = Ntrials * p
v = Ntrials * p * (1.0 - p) * (1.0 + s * (Ntrials - 1) * rho)
y = rnorm(n, mean = m, sd = sqrt(v))

formula = y ~ 1 + z
data = data.frame(y, z, s)
r = inla(formula, data = data, scale = s, 
         family = "betabinomialna", Ntrials=Ntrials, verbose = TRUE,
         control.inla = list(strategy = "adaptive"))
summary(r)

