n <- 100
a <- 1
b <- 0.2
z <- rnorm(n)
eta <- a + b*z
E <- runif(n)
lambda <- E * exp(eta)
y <- numeric(n)
for(i in 1:n) {
    while((y[i] <- rpois(1, lambda[i])) == 0) TRUE
}

result <- inla(y ~ 1 + z, family = "nzpoisson",
               data = data.frame(y, z), E=E)
summary(result)
