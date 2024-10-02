n <- 3000
for (xi in c(-0.3, 0.3)) {
    xi.intern <- inla.models()$likelihood$egp$hyper$theta1$to.theta(xi, interval=c(-0.5, 0.5))
    alpha <- 0.9
    kappa <- 1.1
    kappa.intern <- inla.models()$likelihood$egp$hyper$theta2$to.theta(kappa)
    a <- ((1-alpha^(1/kappa))^(-xi) -1)
    x <- rnorm(n, sd = 0.3)
    eta <- 0.9 + 1.1 * x
    q <- exp(eta)
    sigma <- xi * q / a
    y <- - (sigma / xi) * (1- (1-runif(n)^(1/kappa))^(-xi))

    r <- inla(y ~ 1 + x,
              family = "egp",
              control.family = list(
                  control.link = list(quantile = alpha),
                  hyper = list(tail = list(
                                   ##initial = xi.intern,
                                   fixed = !TRUE, 
                                   prior = "pc.egptail",
                                   param = c(5, -0.5, 0.5)),
                               shape = list(
                                   ##initial = kappa.intern, 
                                   fixed = !TRUE, 
                                   prior = "loggamma",
                                   param = c(100, 100)))), 
              data = data.frame(y, x),
              verbose = !TRUE)
    print(summary(r))
}
