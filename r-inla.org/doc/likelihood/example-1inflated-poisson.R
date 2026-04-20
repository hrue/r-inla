sim.1poisson <- function(prob, lambda) {
    stopifnot(length(prob) == length(lambda) && length(prob) > 0)
    n <- length(lambda)
    y <- numeric(n)
    event <- (runif(n) < prob)
    idx.one <- which(event)
    idx.non.one <- which(!event)
    y[idx.one] <- 1
    for(i in idx.non.one) {
        y[i] <- 0
        while(y[i] == 0) {
            yy <- rpois(20, lambda = lambda[i])
            if (any(yy > 0)) {
                y[i] <- yy[min(which(yy > 0))]
            }
        }
    }
    return (y)
}

## chose link-function to use for the zero-inflation probability
link.simple <- "logit"
inv.link <- inla.link.invlogit
## link.simple <- "probit"
## inv.link <- inla.link.invprobit
## link.simple <- "cloglog"
## inv.link <- inla.link.invcloglog

n <- 10^4
z <- rnorm(n, sd = 0.3)
x <- rnorm(n, sd = 0.2)
xx <- rnorm(n, sd = 0.3)
zz <- rnorm(n, sd = 0.2)
E <- runif(n, min = 0.1, max = 100)

beta <- c(1, 1.1, 2.1, 0, -2, 1.2, 2.2, 0)
eta2 <- beta[1] + beta[2] * xx + beta[3] * zz + beta[4] * xx * zz
eta1 <- beta[5] + beta[6] * x + beta[7] * z + beta[8] * x * z
prob <- inv.link(eta1)
lambda <- E*exp(eta2)

y <- sim.1poisson(prob, lambda)

## head(data.frame(y, E, x, z, xx, zz))

r <- inla(
    inla.mdata(cbind(y, E), cbind(1, x, z, x*z)) ~ 1 + xx + zz + xx*zz,
    family = "1poisson",
    data = data.frame(y, E, x, z, xx, zz), 
    control.fixed = list(prec = 1, prec.intercept = 1), 
    control.compute = list(cpo = TRUE), 
    control.family = list(link.simple = link.simple,
                          hyper = list(beta1 = list(param = c(0, 1)),
                                       beta2 = list(param = c(0, 1)),
                                       beta3 = list(param = c(0, 1)),
                                       beta4 = list(param = c(0, 1)),
                                       beta5 = list(param = c(0, 1)))))

rr <- inla(
    inla.mdata(cbind(y, E), cbind(1, xx, zz, xx*zz)) ~ 1 + x + z + x*z,
    family = "1poissonS",
    data = data.frame(y, E, x, z, xx, zz), 
    control.fixed = list(prec = 1, prec.intercept = 1), 
    control.compute = list(cpo = TRUE), 
    control.inla = list(cmin = 0), 
    ## in this case we need to define link.simple as the main link
    control.family = list(control.link = list(model = link.simple),
                          hyper = list(beta1 = list(param = c(0, 1)),
                                       beta2 = list(param = c(0, 1)),
                                       beta3 = list(param = c(0, 1)),
                                       beta4 = list(param = c(0, 1)),
                                       beta5 = list(param = c(0, 1))))) 
summary(r)
summary(rr)

res <- cbind("beta" = beta, 
             "1poisson" = c(r$summary.fixed$mean, r$summary.hyperpar$mean),
             "1poissonS" = c(rr$summary.hyperpar$mean, rr$summary.fixed$mean))
res <- cbind(res,
             diff = (res[, 2]-beta), 
             diffS = (res[, 3]-beta),
             "diff/sd" = (res[, 2]-beta) / c(r$summary.fixed$sd, r$summary.hyperpar$sd), 
             "diffS/sd" = (res[, 3]-beta) / c(rr$summary.hyperpar$sd, rr$summary.fixed$sd))
mm <- nrow(res) %/% 2
rownames(res) <- c(paste0("beta", 1:mm, ".poisson"), paste0("beta", 1:mm, ".prob"))
print(round(dig = 2, res))
