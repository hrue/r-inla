n = 100
beta = 2
w = rnorm(n)
prec.u = 1
prec.y = 1
## heteroscedastic scaling
s = runif(n,min=0,max=1)
## true but unobserved covariate
x = w + rnorm(n, sd = 1/sqrt(s*prec.u))
y = 1 + beta*x + rnorm(n, sd = 1/sqrt(prec.y))

## prior parameters
prior.beta = c(0, 0.0001)
prior.prec.u = c(10, 9/prec.u)
prior.prec.y = c(10, 9/prec.y)


formula = y ~ f(w, model="meb", scale=s, values=w,
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
        )
    )
)

r = inla(formula, data = data.frame(y, w, s),
    family = "gaussian",
    control.family = list(
        hyper = list(
            prec = list(param = prior.prec.y,
                fixed = FALSE
            )
        )
    )
)
