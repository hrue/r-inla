n=100
a = 0
b = 1
x = rnorm(n, sd = 0.5)
eta = a + b*x
interval = c(1, 4)
E = sample(1:10, n, replace=TRUE)
lambda = E*exp(eta)
y = rpois(n, lambda = lambda)

censored = (y >= interval[1] & y <= interval[2])
y[censored] = interval[1]

r = (inla(y ~ 1 + x, 
          family = "cenpoisson",
          control.family = list(cenpoisson.I = interval), 
          data = data.frame(y, x),
          E=E))
summary(r)



