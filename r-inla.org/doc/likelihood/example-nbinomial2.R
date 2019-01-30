n = 300
x = rnorm(n, sd = 0.2)
eta = 1 + 1.1*x
p = exp(eta)/(1 + exp(eta))
size = sample(1:5, n, replace=TRUE)
y = rnbinom(n, size = size, prob = p)
r = inla(y ~ 1 + x, family = "nbinomial2", Ntrials = size,
         data = data.frame(y, x, size))
summary(r)
