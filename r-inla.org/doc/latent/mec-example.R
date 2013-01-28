n = 100
prec.y = 100
prec.obs = 10
prec.x = 1
## true unobserved covariate
x = rnorm(n, sd = 1/sqrt(prec.x)) 
## the observed covariate
xobs = x + rnorm(n, sd = 1/sqrt(prec.obs))
## regression model using the unobserved 'x'
y = 1 + 4*x + rnorm(n, sd = 1/sqrt(prec.y))

## prior parameters
prior.prec = c(1, 0.01)
prior.beta = c(0, 0.1)

formula = y ~ 1 + 
    f(xobs, model="me",
      hyper = list(
              beta = list(
                      param = prior.beta,
                      fixed = FALSE
                      ),
              prec.obs = list(
                      param = prior.prec,
                      initial = log(prec.obs),
                      fixed = TRUE
                      ),
              prec.x = list(
                      param = prior.prec,
                      initial = log(prec.x),
                      fixed = FALSE
                      ),
              mean.x = list(
                      initial = 0,
                      fixed=TRUE
                      )
              )
      )

r = inla(formula,
        data = data.frame(y, xobs), 
        family = "gaussian",
        control.family = list(
                hyper = list(
                        prec = list(param = prior.prec, 
                                initial = log(prec.y),
                                fixed=FALSE
                                )
                        )
                )
        )

summary(r)
