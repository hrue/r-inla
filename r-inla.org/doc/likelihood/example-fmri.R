n <- 300
x <- rnorm(n, sd = 0.3)
df <- 1
prec <- 3
eta <- 1 + x
lambda <- exp(eta)
y <- sqrt(rchisq(n, df = df, ncp = prec * lambda^2) /prec)

r <- inla(y ~ 1 + x, 
          data = data.frame(y, x),
          family = "fmri",
          control.family = list(hyper = list(df = list(initial = df))), 
          control.inla = list(cmin = 0,
                              int.strategy = "eb",
                              strategy = "adaptive"), 
          verbose = TRUE)
summary(r)

## 'cmin=0' seems to be required only for initial values that can give
## 'crazy' values. We can rerun without this re-starting at the prev fit,
## to validate
r$.args$control.inla$cmin <- -Inf
r$.args$control.inla$int.strategy <- "auto"
rr <- inla.rerun(r)
summary(rr)
