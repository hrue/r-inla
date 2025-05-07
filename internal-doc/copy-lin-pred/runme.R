library(INLA)
library(dplyr)

## y1 = intercept.1 + beta.x * x + beta.xx * xx
## y2 = intercept.2 + beta * (beta.x * x + beta.xx * xx)


intercept.1 <- 0
intercept.2 <- 1
beta.x <- 1
beta.xx <- 2
beta <- -1

n <- 1000
x <- rnorm(n)
xx <- rnorm(x)

y1 <- intercept.1 + beta.x * x + beta.xx * xx + rnorm(n)
y2 <- intercept.2 + beta * (beta.x * x + beta.xx * xx) + rnorm(n)

## we need to first set v = beta.x * x + beta.xx * xx

Y <- matrix(NA, 3*n, 3)
Y[1:n, 1] <- 0
data0 <- data.frame(x, xx, v = 1:n, w = rep(-1, n))

Y[n + 1:n, 2] <- y1
data1 <- data.frame(intercept.1 = rep(intercept.1, n), x, xx)

Y[2*n + 1:n, 3] <- y2
data2 <- data.frame(intercept.2 = rep(intercept.2, n), v.copy = 1:n)

model <- Y ~ -1 +
    intercept.1 + intercept.2 +
    x + xx +
    f(v, w, model = "iid", hyper = list(prec = list(initial = -5, fixed = TRUE))) +
    f(v.copy, copy = "v", hyper = list(beta = list(fixed = FALSE)))

data <- as.list(bind_rows(data0, data1, data2))
data$Y <- Y

r <- inla(model,
          data = data,
          family = c("normal", "normal", "normal"),
          control.family = list(
              list(hyper = list(prec = list(initial = 15, fixed = TRUE))),
              list(),
              list()))
summary(r)
