## simple example to illustrate the use of 'copy'
set.seed(1234)
N <- 100
n <- 50
u <- scale(rnorm(n))
s <- 0.1

i <- sample(1:n, N, replace = TRUE)
j <- sample(1:n, N, replace = TRUE)
k <- sample(1:n, N, replace = TRUE)
y <- u[i] + u[j] + u[k] + rnorm(N, sd = s)

r <- inla(y ~ -1 +
              f(i, values = 1:n) +
              f(j, copy = "i") +
              f(k, copy = "i"),
          data = data.frame(y, i, j, k),
          control.family = list(hyper = list(
                                prec = list(initial = log(1/s^2),
                                            fixed = TRUE))))
plot(u, r$summary.random$i$mean, pch = 19)
abline(a = 0, b = 1, lwd = 3, col = "blue")

## estimate scaling parameters, assuming
## y <- u[i] + beta.j * u[j] + beta.k * u[k] + rnorm(N, sd = s)
## where the true values are beta.j=1 and beta.k=1

rr <- inla(y ~ -1 +
              f(i, values = 1:n) +
              f(j, copy = "i", hyper = list(
                                   beta = list(fixed = FALSE))) +
              f(k, copy = "i", hyper = list(
                                   beta = list(fixed = FALSE))), 
          data = data.frame(y, i, j, k),
          control.family = list(hyper = list(
                                prec = list(initial = log(1/s^2),
                                            fixed = TRUE))))
rr$summary.hyperpar[,c("mean","sd")]
inla.dev.new()
plot(u, rr$summary.random$i$mean, pch = 19)
abline(a = 0, b = 1, lwd = 3, col = "blue")

## now we assume that we know that beta.k = beta.j, and we estimate
## just beta.j

rrr <- inla(y ~ -1 +
              f(i, values = 1:n) +
              f(j, copy = "i", hyper = list(
                                   beta = list(fixed = FALSE))) +
              f(k, copy = "i", same.as = "j"), 
          data = data.frame(y, i, j, k),
          control.family = list(hyper = list(
                                prec = list(initial = log(1/s^2),
                                            fixed = TRUE))))
rrr$summary.hyperpar[,c("mean","sd")]
inla.dev.new()
plot(u, rrr$summary.random$i$mean, pch = 19)
abline(a = 0, b = 1, lwd = 3, col = "blue")
