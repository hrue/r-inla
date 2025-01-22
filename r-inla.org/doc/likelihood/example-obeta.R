library(INLA)

robeta <- function(eta, k1, k2, precision)
{
    stopifnot(k1 < k2)
    g <- function(...) inla.link.invlogit(...)
    p <- cbind(1 - g(eta - k1), g(eta - k1) - g(eta - k2), g(eta - k2))
    eta <- as.vector(eta)
    n <- length(eta)
    group <- numeric(n)
    x <- numeric(n)
    for(i in 1:n) {
        group[i] <- sample(1:3, 1, prob = p[i, ])
        if (group[i] == 1) {
            x[i] <- 0
        } else if (group[i] == 3) {
            x[i] <- 1
        } else {
            mu <- g(eta[i])
            a <- mu * precision
            b <- -mu * precision + precision
            x[i] <- rbeta(1, a, b)
        }
    }
    return (x)
}

n <- 10^4
x <- rnorm(n, sd = 0.3)
eta <- 1 + x
k1 <- -2
k2 <- 2
prec <- 2
y <- robeta(eta, k1, k2, prec)

r <- inla(y ~ 1 + x,
          data = data.frame(y, x),
          family = "obeta",
          control.inla = list(cmin = 0), 
          verbose = TRUE)
summary(r)

## produce posteriors for k1 and k2
xx <- inla.hyperpar.sample(10^5, r, intern = FALSE)
loc <- xx[, "offset location-parameter for the obeta observations"]
width <- xx[, "offset width-parameter for the obeta observations"]
k1.mc <-  loc - lwidth
k2.mc <- loc + lwidth
par(mfrow = c(1, 2))
hist(k1.mc, prob=TRUE, n = 100); abline(v=k1, lwd = 5, col = "blue")
hist(k2.mc, prob=TRUE, n = 100); abline(v=k2, lwd = 5, col = "blue")
