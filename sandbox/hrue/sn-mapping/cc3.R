library(sn)
library(numDeriv)

g <- function(x, par) {
    return (dsn(x, xi = par$xi, omega = par$omega, alpha = par$alpha, log = TRUE))
}
gd1 <- function(x, par) {
    return (grad(g, x, par = par, method = "simple"))
}
gd2 <- function(x, par) {
    return (grad(gd1, x, par = par, method = "simple"))
}
gd3 <- function(x, par) {
    return (grad(gd2, x, par = par, method = "simple"))
}


ss <- c()
mm <- c()
d3 <- c()
aa <- c()
m <- 0
for (s in seq(-0.99, 0.99, by = 0.01)) {

    par <- INLA:::inla.sn.reparam(c(0, 1, s))
    m <- optim(m, fn = g, gr = NULL, method = "BFGS",
               control = list(fnscale = -1),
               par)$par

    ss <- c(ss, s)
    mm <- c(mm, m)
    d3 <- c(d3, gd3(m, par))
    aa <- c(aa, par$alpha)
}

plot(d3, ss, ylim=c(-1,1), xlim=c(-5,5),  type = "l", lwd = 3)
inla.dev.new()
plot(d3, aa, ylim = c(-5, 5), xlim=c(-5,5),  type = "l", lwd = 3)
