sim.binomial <- function(prob, p, size) 
{
    ## - prob=zero-inflation-prob
    ## - binomial(size, p)
    stopifnot(length(prob) == length(p) && length(prob) == length(size)
              && length(prob) > 0)
    n <- length(prob)
    y <- numeric(n)
    event <- (runif(n) < prob)
    idx.zero <- which(event)
    idx.non.zero <- which(!event)
    y[idx.zero] <- 0
    y[idx.non.zero] <- rbinom(length(idx.non.zero),
                              size = size[idx.non.zero],
                              prob = p[idx.non.zero])
    return (y)
}

n <- 1000
z <- rnorm(n, sd = 0.3)
x <- rnorm(n, sd = 0.2)
xx <- rnorm(n, sd = 0.3)
zz <- rnorm(n, sd = 0.2)
Ntrials <- sample(1:10, n, replace = TRUE)

## chose link-function to use for the zero-inflation probability
link.simple <- "logit"
inv.link <- inla.link.invlogit
## link.simple <- "probit"
## inv.link <- inla.link.invprobit
## link.simple <- "cloglog"
## inv.link <- inla.link.invcloglog

beta <- c(1, 1.1, 2.1, 0, -2, 1.2, 2.2, 0)
eta2 <- beta[1] + beta[2] * xx + beta[3] * zz + beta[4] * xx * zz
eta1 <- beta[5] + beta[6] * x + beta[7] * z + beta[8] * x * z
prob <- inv.link(eta1)
p <- 1/(1 + exp(-eta2))

ok <- FALSE
while(!ok) {
    y <- sim.binomial(prob, p, Ntrials)
    ok <- !all(y == 0)
}

## head(data.frame(y, Ntrials, x, z, xx, zz))

r <- inla(
    inla.mdata(cbind(y, Ntrials), cbind(1, x, z, x*z)) ~ 1 + xx + zz + xx*zz,
    family = "0binomial",
    data = data.frame(y, Ntrials, x, z, xx, zz), 
    control.fixed = list(prec = 1, prec.intercept = 1), 
    control.compute = list(cpo = TRUE), 
    control.family = list(link.simple = link.simple, 
                          hyper = list(beta1 = list(param = c(0, 1)),
                                       beta2 = list(param = c(0, 1)),
                                       beta3 = list(param = c(0, 1)),
                                       beta4 = list(param = c(0, 1)),
                                       beta5 = list(param = c(0, 1)))))

rr <- inla(
    inla.mdata(cbind(y, Ntrials), cbind(1, xx, zz, xx*zz)) ~ 1 + x + z + x*z,
    family = "0binomialS",
    data = data.frame(y, Ntrials, x, z, xx, zz), 
    control.fixed = list(prec = 1, prec.intercept = 1), 
    control.compute = list(cpo = TRUE), 
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
             "0binomial" = c(r$summary.fixed$mean, r$summary.hyperpar$mean),
             "0binomialS" = c(rr$summary.hyperpar$mean, rr$summary.fixed$mean))
res <- cbind(res,
             diff = (res[, 2]-beta), 
             diffS = (res[, 3]-beta),
             "diff/sd" = (res[, 2]-beta) / c(r$summary.fixed$sd, r$summary.hyperpar$sd), 
             "diffS/sd" = (res[, 3]-beta) / c(rr$summary.hyperpar$sd, rr$summary.fixed$sd))
mm <- nrow(res) %/% 2
rownames(res) <- c(paste0("beta", 1:mm, ".binomial"), paste0("beta", 1:mm, ".prob"))
print(round(dig = 2, res))
