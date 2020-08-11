n <- 4
A <- matrix(rnorm(n^2), n, n)
b <- rnorm(n)
z <- rnorm(n)
z2x <- function(z) A %*% z
x2z <- function(x) solve(A) %*% x
x <- z2x(z)

f <- function(x) (sum(b*x))^2/2
g <- function(z) f(z2x(z))

library(numDeriv)

grad.f <- grad(f, x)
grad.g <- grad(g, z)

cbind(grad.f,  t(solve(A)) %*% grad.g)

print(round(dig = 4, hessian(f, x)))
print(round(dig = 4, t(solve(A)) %*% hessian(g, z) %*% solve(A)))

