n <- 30000
x <- rnorm(n)
eta <- 1 + 0.2 * x
event <- rep(1, n)
E <- runif(n)

offset <- rnorm(n, sd = 0.5)
xx <- rnorm(n)
xxx <- rnorm(n)
eta.c <- offset + 0.3 * xx + 0.5 * xxx

## need two for the censoring
y <- rpois(n, E*exp(eta))
yy <- rpois(n, E*exp(eta))

for(i in 1:n) {
    event[i] <- sample(c(1, 0, 99), 1, prob = c(0.6, 0.1, 0.3))
    if (event[i] == 1) {
        ## y[i] <- y[i]
    } else if (event[i] == 0) {
        y[i] <- min(y[i], yy[i])
    } else {
        prob <- 1/(1+exp(-eta.c[i]))
        if (runif(1) < prob) {
            ## local.event = 1
            ## y[i] <- y[i]
        } else {
            ## local.event = 0
            y[i] <- min(y[i], yy[i])
        }
    }
}

Y <- inla.mdata(y, E, event, offset, xx, xxx)
r <- inla(Y ~ 1 + x,
          data = list(Y = Y, x = x),
          family = "rcpoisson",
          control.family = list(hyper = list(beta1 = list(param = c(0, 1)), 
                                             beta2 = list(param = c(0, 2)))), 
          verbose = !TRUE)
summary(r)
