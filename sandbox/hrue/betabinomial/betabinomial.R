## from Gamma(1+z) = z Gamma(z), we can simplify Beta(y+a, n-y+b) / Beta(a, b), and avoid
## special functions

lnormc <- function(y, n, a, b) {
    return (lbeta(y+a, n-y+b) - lbeta(a, b))
}

lfac <- function(n, a) {
    if (n == 0) {
        return (0.0)
    } else {
        return (sum(log((1:n) - 1 + a)))
    }
}

lgamma2 <- function(n, a) {
    if (n == 0) {
        return (lgamma(a))
    } else {
        aa <- a-n
        return (sum(log((1:n) - 1 + aa)) + lgamma(aa))
    }
}

lnormc2 <- function(y, n, a, b) {
    return (lfac(y, a) + lfac(n-y, b) - lfac(n, a+b))
}

a <- runif(1)
b <- runif(1)
y <- sample(1:5, 1)
n <- y + sample(0:3, 1)

c(y, n, a, b)
lnormc(y, n, a, b) - lnormc2(y, n, a, b)
lnormc(y, n, a, b)
lnormc2(y, n, a, b)
