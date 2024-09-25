n <- 30000

## this makes it to easy to its just to check the likelihood
## implementation
size <- sample(10:20, n, replace = TRUE)

beta1 <- rnorm(3, sd = 0.2)
beta2 <- rnorm(3, sd = 0.2)
beta <- c(beta1, beta2)
Z <- matrix(NA, n, 6)
W <- matrix(NA, n, 2)
Y <- matrix(NA, n, 2)

x <- rnorm(n, sd = 0.1)
xx <- rnorm(n, sd = 0.4)
eta <-  1 + x + xx

for (i in 1:n) {
    Z[i, ] <- rnorm(6)
    w <- c(rbeta(2, 1, 10), rbeta(1, 10, 1))
    w <- w/sum(w)
    W[i, ] <- w[1:2]

    p1 <- inla.link.invlogit(sum(beta1 * Z[i, 1:3]))
    p2 <- inla.link.invlogit(sum(beta2 * Z[i, 3 + 1:3]))
    p3 <- inla.link.invlogit(eta[i])
    p <- w[1] * p1 + w[2] * p2 + w[3] * p3
    Y[i, ] <- c(rbinom(1, size = size[i], prob = p), size[i])
}

r <- inla(inla.mdata(Y, Z, W) ~ 1 + x + xx, 
          family = "binomialmix",
          data = list(Y = Y, Z = Z, W = W, x = x, xx = xx),
          verbose = TRUE)

cbind(estimate = r$summary.fixed[,"mean"], true = 1)
cbind(estimate =  r$summary.hyperpar[,"mean"], true = beta)
