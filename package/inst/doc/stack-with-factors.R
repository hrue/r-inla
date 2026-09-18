## ----setup, include=FALSE-------------------------------------------
knitr::opts_chunk$set(echo=TRUE, message=FALSE, warning=FALSE)
knitr::opts_chunk$set(fig.path="figures/stack-with-factors/")
set.seed(123)
library(INLA)
inla.setOption(num.threads="1:1")
inla.setOption(smtp="taucs")
if (file.exists("myinit.R")) source("myinit.R")

## ----dat------------------------------------------------------------
  dat <- expand.grid(a=0:1, b=0:1, d=0:2) ### three discrete
  dat$x <- runif(12) ### one continuous covariate
  dat$y <- 3*dat$a - dat$b + (dat$d-1)*2 + dat$x + rnorm(12,0,0.1)

## ----factor---------------------------------------------------------
  dat$d <- factor(dat$d) 

## ----fit------------------------------------------------------------
	coef(lm(y ~ 1 + a + b + d + x, dat))

## ----fit2-----------------------------------------------------------
	coef(lm(y ~ a + b + d + x, dat))

## ----fit3-----------------------------------------------------------
	coef(lm(y ~ -1 + a + b + d + x, dat))

## ----fit4-----------------------------------------------------------
	coef(lm(y ~ 0 + a + b + d + x, dat))

## ----orderi---------------------------------------------------------
	coef(lm(y ~ b + a + d + x, dat))
	coef(lm(y ~ d + a + b + x, dat))

## ----factor2--------------------------------------------------------
dat$a <- factor(dat$a, levels=0:1, labels=c('1st', '2nd'))
dat$b <- factor(dat$b, levels=1:0, labels=c('2nd', '1st')) ## OBS: reference level changed

## ----reg3f----------------------------------------------------------
coef(lm(y ~ a + b + d + x, dat))

## ----orderf---------------------------------------------------------
	coef(lm(y ~ 0 + a + b + d + x, dat))
	coef(lm(y ~ 0 + b + a + d + x, dat))
	coef(lm(y ~ 0 + d + a + b + x, dat))

## ----mm, results='hide'---------------------------------------------
model.matrix(~a+b+d+x, dat) 

## ----tokyo----------------------------------------------------------
data(Tokyo)
str(Tokyo)

## ----tokyomesh------------------------------------------------------
knots <- seq(1, 367, length = 25)
mesh <- inla.mesh.1d(knots, interval = c(1, 367), degree = 2, boundary = "cyclic")
spde <- inla.spde2.pcmatern(mesh, 
    prior.sigma=c(1, 0.01), ## P(sigma > 1) = 0.01
    prior.range=c(1, 0.01)) ## P(range < 1) = 0.01
A <- inla.spde.make.A(mesh, loc = Tokyo$time)
time.index <- inla.spde.make.index("time", n.spde = spde$n.spde)

## ----addf-----------------------------------------------------------
Tokyo$a <- factor(rbinom(366, 1, 0.5))
Tokyo$b <- factor(rbinom(366, 2, 0.5))
Tokyo$x <- runif(366)

## ----ab-------------------------------------------------------------
abx <- model.matrix(~a+b+x, Tokyo)[, -1]

## ----tokyostack-----------------------------------------------------
stack <- inla.stack(
  data = list(y = Tokyo$y, link = 1, Ntrials = Tokyo$n),
  A = list(A, 1),
  effects = list(time.index, data.frame(mu0=1, abx)),
  tag = "est")
formula <- y ~ 0 + mu0 + a1 + b1 + b2 + x + f(time, model = spde)
data <- inla.stack.data(stack)
result <- inla(formula, family = "binomial", 
              data = data, 
              Ntrials = data$Ntrials,
              control.predictor = list(
                A = inla.stack.A(stack), 
                link = data$link, 
                compute = TRUE))
result$summary.fixed[, 1:5]

## ----preds----------------------------------------------------------
pred.sc <- expand.grid(a1=0:1, b1=0:1, b2=0, x=c(0.5))
pred.sc

## ----predt----------------------------------------------------------
A.pred <- inla.spde.make.A(mesh, loc=rep(180, nrow(pred.sc)))

## ----stackpred------------------------------------------------------
stack.pred <- inla.stack(
  data = list(y = NA, link = 1, Ntrials = 2),
  A = list(A.pred, 1),
  effects = list(time.index, data.frame(mu0=1, pred.sc)),
  tag = "pred")
stack.full <- inla.stack(stack, stack.pred)
data <- inla.stack.data(stack.full)
result <- inla(formula, family = "binomial", 
               data = data, 
               Ntrials = data$Ntrials,
               control.predictor = list(
                 A = inla.stack.A(stack.full), 
                 link = data$link, 
                 compute = TRUE),
               control.mode=list(theta=result$mode$theta, 
                                 restart=FALSE))

## ----predicts-------------------------------------------------------
idx.pred <- inla.stack.index(stack.full, tag='pred')$data
result$summary.fitted.val[idx.pred, 1:5]

