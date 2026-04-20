n = 300
x = rnorm(n,  sd = 0.2)
eta = 1 + x
E = runif(n, min = 0.5, max=2)
S = runif(n, min = 0.5, max=2)

## variant 0
mu = E * exp(eta)
size = 1
prob = size/(size + mu)
y = rnbinom(n, size, mu=mu)

y.low <- sample(c(1:4, Inf), n, replace = TRUE)
y.high <- y.low + sample(c(1:5, Inf), n, replace = TRUE)

Y <- inla.mdata(cbind(y, y.low, y.high))
r = inla(Y ~ 1 + x,
         data = list(Y = Y, x = x, E = E),
         family = "cennbinomial2",
         control.family = list(variant = 0), 
         E=E, scale = S)
summary(r)
censored = which(y >= y.low & y <= y.high)
print(length(censored)/ n)

## variant 1
mu = E * exp(eta)
size = 1 * E
prob = size/(size + mu)
y = rnbinom(n, size, mu=mu)

y.low <- sample(c(1:4, Inf), n, replace = TRUE)
y.high <- y.low + sample(c(1:5, Inf), n, replace = TRUE)

Y <- inla.mdata(cbind(y, y.low, y.high))
r = inla(Y ~ 1 + x,
         data = list(Y = Y, x = x, E = E),
         family = "cennbinomial2",
         control.family = list(variant = 1), 
         E=E)

summary(r)
censored = which(y >= y.low & y <= y.high)
print(length(censored)/ n)

## variant 2
mu = E * exp(eta)
size = 1 * S 
prob = size/(size + mu)
y = rnbinom(n, size, mu=mu)

y.low <- sample(c(1:4, Inf), n, replace = TRUE)
y.high <- y.low + sample(c(1:5, Inf), n, replace = TRUE)

Y <- inla.mdata(cbind(y, y.low, y.high))
r = inla(Y ~ 1 + x,
         data = list(Y = Y, x = x, E = E, S = S),
         family = "cennbinomial2",
         control.family = list(variant = 2), 
         E=E, scale = S)

summary(r)
censored = which(y >= y.low & y <= y.high)
print(length(censored)/ n)
