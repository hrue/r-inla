n <- 10^5
## high value of 'size' makes it easy to check the likelihood implementation
size <- sample(100:200, n, replace = TRUE)

beta.p1 <- rnorm(4, sd = 0.5)
beta.p2 <- rnorm(4, sd = 0.5)
beta9 <- rnorm(1, sd = 0.5)
beta <- c(beta.p1, beta.p2, beta9)

Z <- matrix(NA, n, 11)
W <- matrix(NA, n, 2)
Y <- matrix(NA, n, 2)

x <- rnorm(n, sd = 0.5)
xx <- rnorm(n, sd = 0.5)
eta <- numeric(n)

for (i in 1:n) {
    Z[i, ] <- rnorm(11)
    w <- c(rbeta(2, 1, 10), rbeta(1, 10, 1))
    w <- w/sum(w)
    W[i, ] <- w[1:2]

    p1 <- inla.link.invlogit(sum(beta.p1 * Z[i, 1:4]) + beta9 * Z[i, 9])
    p2 <- inla.link.invlogit(sum(beta.p2 * Z[i, 4 + 1:4]) + beta9 * Z[i, 10])

    eta[i] <-  1 + x[i] + xx[i] + beta9 * Z[i, 11]
    p3 <- inla.link.invlogit(eta[i])

    p <- w[1] * p1 + w[2] * p2 + w[3] * p3
    Y[i, ] <- c(rbinom(1, size = size[i], prob = p), size[i])
}

r <- inla(inla.mdata(Y, Z, W) ~ 1 + x + xx, 
          family = "binomialmix",
          data = list(Y = Y, Z = Z, W = W, x = x, xx = xx),
          verbose = TRUE,
          control.inla = list(int.strategy = "eb"))

print(round(dig = 4, cbind(estimate = r$summary.fixed[,"mean"], true = 1)))
print(round(dig = 4, cbind(estimate =  r$summary.hyperpar[,"mean"], true = beta)))

