n <- 10^5
size <- sample(5:10, n, replace = TRUE)

m <- 5
m2 <- 2*m
beta.p1 <- rnorm(m, sd = 1/sqrt(m))
beta.p2 <- rnorm(m, sd = 1/sqrt(m))
beta.common <- rnorm(1, sd = 0.5)
beta <- c(beta.p1, beta.p2, beta.common)

Z <- matrix(NA, n, m2+3)
W <- matrix(NA, n, 2)
Y <- matrix(NA, n, 2)

x <- rnorm(n, sd = 0.5)
xx <- rnorm(n, sd = 0.5)
eta <- numeric(n)

for (i in 1:n) {
    Z[i, ] <- rnorm(m2+3)
    w <- c(rbeta(2, 1, 10), rbeta(1, 10, 1))
    w <- w/sum(w)
    W[i, ] <- w[1:2]

    p1 <- inla.link.invlogit(sum(beta.p1 * Z[i, seq_len(m)]) + beta.common * Z[i, m2+1])
    p2 <- inla.link.invlogit(sum(beta.p2 * Z[i, m + seq_len(m)]) + beta.common * Z[i, m2+2])

    eta[i] <-  1 + x[i] + xx[i] + beta.common * Z[i, m2+3]
    p3 <- inla.link.invlogit(eta[i])

    p <- w[1] * p1 + w[2] * p2 + w[3] * p3
    Y[i, ] <- c(rbinom(1, size = size[i], prob = p), size[i])
}

r <- inla(inla.mdata(Y, Z, W) ~ 1 + x + xx, 
          family = "binomialmix",
          data = list(Y = Y, Z = Z, W = W, x = x, xx = xx),
          verbose = TRUE,
          control.inla = list(int.strategy = "eb"))
## r <- inla.rerun(r)

print(round(dig = 4, cbind(estimate = r$summary.fixed[,"mean"], true = 1)))
print(round(dig = 4, cbind(estimate =  r$summary.hyperpar[,"mean"], true = beta)))
