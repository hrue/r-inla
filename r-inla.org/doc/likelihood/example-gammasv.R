n <- 300
phi <- 1.1
s <- runif(n, min = 0.8, max = 1.25)
x <- rnorm(n)
eta <- 1 + 0.2 * x
mu <- exp(eta)
a <- phi * s * mu
b <- phi * s
y <- rgamma(n, shape = a, rate = b)
r <- inla(y ~ 1 + x,
          family = "gammasv",
          scale = s,
          data = data.frame(y, x, s))
summary(r)
