N=1000
a = -1
b = 1
z = rnorm(N)
eta = a + b*z
n = sample(c(1,5,10,15), size=N, replace=TRUE)
p = exp(eta)/(1 + exp(eta))
prob = 1.0 - (1-p)^n
k = sample(c(1,5,10,15), size=N, replace=TRUE)
y = rbinom(N, size=k, prob = prob)
    
data = list(y=y,z=z)
formula = y ~ 1+z
result = inla(formula, family = "cbinomial", data = data,
        Ntrials=cbind(k, n), verbose=TRUE)
summary(result)
