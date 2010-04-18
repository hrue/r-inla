n=100
a = 1
b = 1
z = rnorm(n)
eta = a + b*z
Ntrials = sample(c(1,5,10,15), size=n, replace=TRUE)
prob = exp(eta)/(1 + exp(eta))
y = rbinom(n, size=Ntrials, prob = prob)

data = list(y=y,z=z)
formula = y ~ 1+z
result = inla(formula, family = "binomial", data = data, Ntrials=Ntrials)
summary(result)
