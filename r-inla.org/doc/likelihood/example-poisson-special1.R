n <- 300
a <- 1
b <- 1
p <- 0.2
x <- rnorm(n, sd = 0.2)
mu <- exp(a+b*x)
y.max <- ceiling(max(mu + 10*sqrt(mu)))
y <- numeric(n)

for(i in 1:n) {
    yy <- 1:y.max
    dy <- dpois(yy, lambda = mu[i])
    dy <- dy/sum(dy)
    dy <- dy * (1-p)
    dy[1] <- dy[1] + p
    y[i] <- sample(x = yy, size = 1, prob = dy)
}

r <- inla(y ~ 1 + x,
          data = data.frame(y, x),
          family = "poisson.special1")
summary(r)
