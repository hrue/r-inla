n = 100
beta = 4
prec.y = 1
prec.u = 1
prec.x = 1
## true unobserved covariate
x = rnorm(n, sd = 1/sqrt(prec.x))
## the observed covariate with heteroscedastic scaling
s = runif(n,min=0.5,max=2)
w = x + rnorm(n, sd = 1/sqrt(prec.u*s))
## regression model using the unobserved 'x'
y = 1 + beta*x + rnorm(n, sd = 1/sqrt(prec.y))

## prior parameters
prior.beta = c(0, 0.0001)
prior.prec.u = c(10, 9)
prior.prec.x = c(10, 9)
prior.prec.y = c(10, 9)


formula = y ~ 1 +
f(w, model="mec", scale=s, values=w,
    hyper = list(
        beta = list(
            prior = "gaussian",
            param = prior.beta,
            fixed = FALSE
        ),
        prec.u = list(  
            prior = "loggamma",
            param = prior.prec.u,
            initial = log(prec.u),
            fixed = FALSE
        ),
        prec.x = list(
            prior = "loggamma",
            param = prior.prec.x,
            initial = log(prec.x),
            fixed = FALSE
        ),
        mean.x = list(
            prior = "gaussian",
            initial = 0,
            fixed=TRUE
        )
    )
)

r = inla(formula,
    data = data.frame(y, w, s),
    family = "gaussian",
    control.family = list(
        hyper = list(
            prec = list(param = prior.prec.y,
                initial = log(prec.y),
                fixed=FALSE
            )
        )
    )
)

summary(r)


