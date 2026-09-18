## ----setup, include=FALSE-------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  dev = "png",
  dev.args = list(type = "cairo-png"),
  fig.width = 10,
  fig.height = 7, 
  out.width = "69%",
  fig.align = "center"
)
knitr::opts_chunk$set(echo = TRUE)
set.seed(1)
library(INLA)
inla.setOption(smtp="taucs")
inla.setOption(num.threads="1:1")
if (file.exists("myinit.R")) source("myinit.R")

## ----domain---------------------------------------------------------
bb <- rbind(c(0, 10), c(0, 7))

## ----rxy------------------------------------------------------------
rxy <- apply(bb, 1, diff)
rr <- mean(rxy)

## ----loc------------------------------------------------------------
n <- 3000
loc <- cbind(
    runif(n, bb[1, 1], bb[1, 2]),
    runif(n, bb[2, 1], bb[2, 2]))

## ----sfn------------------------------------------------------------
sfn <- function(a, b)
    sin(a - b) + cos(a + b)

## ----sloc-----------------------------------------------------------
u <- sfn(2 * pi * loc[, 1] / rr,
         2 * pi * loc[, 2] / rr) 
summary(u)
par(mar = c(2, 2, 0, 0), mgp = c(1.5, 0.5, 0), las = 1)
plot(loc, asp = 1, pch = 19, cex = 0.5 + (2 + u)/4, xlab = "", ylab = "")

## ----centers--------------------------------------------------------
h <- 1 ## size of the side of each sub-region
group.ce <- as.matrix(expand.grid(
    seq(bb[1, 1] + h/2, bb[1, 2], h),
    seq(bb[2, 1] + h/2, bb[2, 2], h)))
(ng <- nrow(group.ce))

## ----groupid--------------------------------------------------------
group.id <- ceiling(loc[,1]/h) + 
  ceiling(rxy[1]/h) * (ceiling(loc[,2]/h)-1)

## ----xxx------------------------------------------------------------
xxx <- cbind(x1 = runif(n), x2 = runif(n), x3 = runif(n))

## ----x.g------------------------------------------------------------
xxx.g <- aggregate(xxx, list(g = group.id), mean)
str(xxx.g)

## ----y--------------------------------------------------------------
beta.y <- c(5, 0, -3, 3) 
u.g <- tapply(u, group.id, mean)
eta.y.g <- drop(cbind(1, as.matrix(xxx.g[, 2:ncol(xxx.g)])) %*% beta.y) + u.g

## ----s--------------------------------------------------------------
s <- rgamma(n, 7, 3)
summary(s)

## ----ysim-----------------------------------------------------------
sigma.y <- 1
y <- rnorm(n, eta.y.g[group.id], sqrt(sigma.y/s))
summary(y)

## ----ag, warning----------------------------------------------------
library(INLA)
agg <- lapply(1:ng, function(g) {
    ig <- which(group.id==g)
    if(length(ig)>0) 
        return(inla.agaussian(y[ig], s[ig]))
    return(inla.mdata(list(NA, 0, 1, 1, NA))) ### a deal with missing
})
str(YY <- Reduce('rbind', lapply(agg, unlist))) ## five columns matrix

## ----w0-------------------------------------------------------------
beta.w <- c(1, -1, 0, 1)
alpha <- 2
gamma.w <- 0.5
lambdaw <- exp(cbind(1, xxx) %*% beta.w + gamma.w * u)
we0 <- rweibull(n, shape = alpha, scale = 1/lambdaw)
summary(we0)

## ----event----------------------------------------------------------
summary(u.ev <- runif(n, 0.3, 1)) ## censoring factor
table(event <- rbinom(n, 1, 0.5)) ## censored (=0) or event (=1)
summary(we <- ifelse(event == 1, we0, we0 * u.ev)) ## censored outcome

## ----ee-------------------------------------------------------------
summary(ee <- rgamma(n, 10, 2))

## ----po-------------------------------------------------------------
beta.p <- c(2, 1, -1, 0)
gamma.p <- -0.3
delta.p <- exp(cbind(1, xxx) %*% beta.p + gamma.p * u)
po <- rpois(n, delta.p * ee)
summary(po)

## ----mesh-----------------------------------------------------------
mesh <- inla.mesh.2d(
    loc.domain = matrix(bb[c(1,3,3,1,1, 2,2,4,4,2)], ncol = 2), 
    offset = c(.0001, .5) * rr, ## extensions size
    max.edge = c(.05, .3) * rr, 
    cutoff = 0.025 * rr)
mesh$n

## ----vmesh----------------------------------------------------------
par(mar = c(0, 0, 1, 0), mgp = c(1.5, 0.5, 0))
plot(mesh, asp = 1, edge.color = gray(0.5))
points(loc, pch = 19, col = group.id, cex = 0.3 + (2 + u)/4)

## ----spde-----------------------------------------------------------
spde <- inla.spde2.pcmatern(
    mesh,
    prior.range=c(1, 0.1),
    prior.sigma=c(1, 0.1))

## ----formula--------------------------------------------------------
fff <- list(
  inla.mdata(Y),
  inla.surv(w, ev),
  o) ~ 0 +
  a1 + x1y + x2y + x3y + f(s, model = spde) +
  a2 + x1w + x2w + x3w + f(sc1, copy = "s", fixed = FALSE) +
  a3 + x1p + x2p + x3p + f(sc2, copy = "s", fixed = FALSE)

## ----stack1---------------------------------------------------------
stackY <- inla.stack(
    tag = 'y',
    data = list(Y = YY),
    effects= list(
        data.frame(a1 = 1,
                   x1y = xxx.g$x1, 
                   x2y = xxx.g$x2, 
                   x3y = xxx.g$x3),
        s = 1:mesh$n),
    A=list(1,
           inla.spde.make.A(mesh, group.ce))
)

## ----stack2---------------------------------------------------------
stackW <- inla.stack(
    tag = 'w', 
    data = list(w = we,
                ev = event), 
    effects = list(
        data.frame(a2 = 1,
                   x1w = xxx[, 1], 
                   x2w = xxx[, 2], 
                   x3w = xxx[, 3]),
        sc1 = 1:mesh$n),
    A = list(1,
             inla.spde.make.A(mesh, loc))
)

## ----stack3---------------------------------------------------------
stackP <- inla.stack(
    tag = 'p', 
    data = list(o = po,
                E = ee),
    effects = list(
        data.frame(a3 = 1,
                   x1p = xxx[, 1], 
                   x2p = xxx[, 2], 
                   x3p = xxx[, 3]),
        sc2 = 1:mesh$n),
    A = list(1,
             inla.spde.make.A(mesh, loc))
)

## ----jstack---------------------------------------------------------
stacked <- inla.stack(stackY, stackW, stackP)

## ----fit------------------------------------------------------------
result <- inla(
   formula = fff,
    family = c("agaussian", "weibullsurv", "poisson"),
    data = inla.stack.data(stacked),
    E = E, 
    control.predictor = list(
        A=inla.stack.A(stacked)),
    control.family = list(
        list(),
        list(variant = 1),
        list()), 
   control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE)
)
result$cpu

## ----fixef----------------------------------------------------------
round(cbind(
    true = c(beta.y, beta.w, beta.p),
    result$summary.fix), 2)

## ----hpars----------------------------------------------------------
round(cbind(true = c(1/sigma.y^2, alpha, NA, NA, gamma.w, gamma.p), 
            result$summary.hy), 2)


## ----dindex---------------------------------------------------------
idx1 <- inla.stack.index(stacked, tag = "p")$data

## ----fitted---------------------------------------------------------
head(cbind(true = delta.p, 
           result$summary.fitted.values[idx1, ]), 3)

## ----dic------------------------------------------------------------
str(result$dic)
str(result$waic)
str(result$cpo)
result$dic$family.dic ### DIC for each outcome
tapply(result$waic$local.waic, 
       result$dic$family, ## taken from dic family id
       sum) ## WAIC for each outcome
tapply(result$cpo$cpo, 
       result$dic$family, ## taken fron dic family id
        function(x) 
         c(nlCPO = -sum(log(x)))) ## -sum log CPO for each outcome

