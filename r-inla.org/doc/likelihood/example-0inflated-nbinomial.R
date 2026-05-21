sim.nbinom <- function(prob, m, size) 
{
    stopifnot(length(prob) == length(m) &&
              length(size) == length(m) &&
              length(prob) > 0)
    n <- length(m)
    y <- numeric(n)
    event <- (runif(n) < prob)
    idx.zero <- which(event)
    idx.non.zero <- which(!event)
    y[idx.zero] <- 0
    y[idx.non.zero] <- rnbinom(length(idx.non.zero),
                               mu = m[idx.non.zero],
                               size = size[idx.non.zero])
    return (y)
}

## chose link-function to use for the zero-inflation probability
link.simple <- "logit"
inv.link <- inla.link.invlogit

variant <- 0
n <- 10^4
z <- rnorm(n, sd = 0.2)
x <- rnorm(n, sd = 0.2)
xx <- rnorm(n, sd = 0.2)
zz <- rnorm(n, sd = 0.2)
E <- runif(n, min = 0.5, max = 1/0.5)
if (variant == 1) {
    size = 3 * E
} else {
    size <- rep(3, n)
}

xxzz <- 0.2 * scale(xx * zz)
xz <- 0.2 * scale(x*z)
beta <- c(1, 1.1, 2.1, 0, -2, 1.2, 2.2, 0)
eta2 <- beta[1] + beta[2] * xx + beta[3] * zz + beta[4] * xxzz
eta1 <- beta[5] + beta[6] * x + beta[7] * z + beta[8] * xz
prob <- inv.link(eta1)
m <- E*exp(eta2)

ok <- FALSE
while(!ok) {
    y <- sim.nbinom(prob, m, size)
    ok <- !all(y == 0)
}

r <- inla(
    inla.mdata(cbind(y, E), cbind(1, x, z, xz)) ~ 1 + xx + zz + xxzz,
    family = "0nbinomial",
    data = data.frame(y, E, x, z, xz, xx, zz, xxzz), 
    control.fixed = list(prec = 1, prec.intercept = 1), 
    control.compute = list(cpo = TRUE), 
    control.family = list(link.simple = link.simple, variant = variant, 
                          hyper = list(beta1 = list(param = c(0, 1)),
                                       beta2 = list(param = c(0, 1)),
                                       beta3 = list(param = c(0, 1)),
                                       beta4 = list(param = c(0, 1)))))
##r <- inla.rerun(r)
summary(r)

rr <- inla(
    inla.mdata(cbind(y, E), cbind(1, xx, zz, xxzz)) ~ 1 + x + z + xz,
    family = "0nbinomialS",
    data = data.frame(y, E, x, z, xz, xx, zz, xxzz), 
    control.fixed = list(prec = 1, prec.intercept = 1), 
    control.compute = list(cpo = TRUE), 
    control.inla = list(cmin = 0), 
    ## in this case we need to define link.simple as the main link
    control.family = list(control.link = list(model = link.simple),
                          hyper = list(beta1 = list(param = c(0, 1)),
                                       beta2 = list(param = c(0, 1)),
                                       beta3 = list(param = c(0, 1)),
                                       beta4 = list(param = c(0, 1)))))
##rr <- inla.rerun(rr)

summary(r)
summary(rr)

ress <- cbind("beta" = beta, 
             "0nbin" = c(r$summary.fixed$mean, r$summary.hyperpar$mean[-1]),
             "0nbinS" = c(rr$summary.hyperpar$mean[-1], rr$summary.fixed$mean))
res <- cbind(ress,
             diff = (ress[, 2]-beta), 
             diffS = (ress[, 3]-beta),
             "sd" = c(r$summary.fixed$sd, r$summary.hyperpar$sd[-1]), 
             "sdS" = c(rr$summary.hyperpar$sd[-1], rr$summary.fixed$sd), 
             "diff/sd" = (ress[, 2]-beta) / c(r$summary.fixed$sd, r$summary.hyperpar$sd[-1]), 
             "diffS/sd" = (ress[, 3]-beta) / c(rr$summary.hyperpar$sd[-1], rr$summary.fixed$sd))
mm <- nrow(res) %/% 2
rownames(res) <- c(paste0("beta", 1:mm, ".nbin"), paste0("beta", 1:mm, ".prob"))
print(round(dig = 3, res))

round(dig = 3, rbind(r$summary.hyperpar[1,c(1, 2, 6)],
                     rr$summary.hyperpar[1,c(1, 2, 6)]))
plot(r$cpo$pit, rr$cpo$pit)
abline(a=0,b=1,lwd=3,col="blue")
coef(lm(r$cpo$pit ~ rr$cpo$pit))
