n <- 3000
x <- rnorm(n)
eta <- 1 + 0.2 * x
event <- rep(1, n)
E <- runif(n)

offset <- rnorm(n, sd = 0.3)
xx <- rnorm(n)
xxx <- rnorm(n)
eta.c <- offset + 0.3 * xx + 0.4 * xxx

y <- numeric(n)
prob <- 1/(1+exp(-eta.c))
event <- sample(c(1, 0), n, prob = c(0.6, 0.4), replace = TRUE)
prob[which(event == 1)] <- 1
y <- rpois(n, prob * E * exp(eta))

Y <- inla.mdata(y, E, event, offset, xx, xxx)
r <- inla(Y ~ 1 + x,
          data = list(Y = Y, x = x),
          family = "tpoisson",
          control.family = list(hyper = list(beta1 = list(param = c(0, 1)), 
                                             beta2 = list(param = c(0, 2)))), 
          verbose = TRUE,
          debug = TRUE)
summary(r)
