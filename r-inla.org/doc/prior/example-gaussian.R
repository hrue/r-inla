n = 10
y = rnorm(n)
r = inla(y ~ 1, data = data.frame(y),
        control.family = list(
                hyper = list(
                        prec = list(
                                prior = "normal",
                                param = c(0, 1)
                                )
                        )
                )
        )

