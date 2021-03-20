library(flexsurv)
library(INLA)

n <- 1000
alpha <- 1.0
intercept <- 1.1
beta <- 1.2
x <- rnorm(n, sd = 0.2)
eta <- intercept + beta*x
mu <- exp(eta)
event <- rep(1,n)
y <- rgompertz(n, rate = mu, shape = alpha)

r <- inla(y ~ 1 + x,
          family ="gompertz", data=data.frame(y, x))
r.surv <- inla(inla.surv(y, event) ~ 1 + x,
               family ="gompertzsurv", data=data.frame(y, event, x))

## should be 'small'
print(r$mlik - r.surv$mlik)
