n = 1000
beta.x = 2
prec.x = 100
prec.y = 1000 ## fixed

x.true = rnorm(n)
x = x.true + rnorm(n, sd = sqrt(1/prec.x))
y = 1 + beta.x * x.true + rnorm(n, sd = sqrt(1/prec.y))

formula = Y ~ -1 + intercept +
    f(me.fixed.effect, model="iid",
      hyper = list(prec = list(initial = -5, fixed=TRUE)))

intercept = c(rep(1, n), rep(0, n))
me.fixed.effect = rep(1:n, 2)

Y = matrix(NA, 2*n, 2)
Y[1:n, 1] = y
Y[(1:n) + n, 2] = x

r = inla(formula,
        data = list(Y=Y, intercept=intercept,
                me.fixed.effect = me.fixed.effect),
        family = c("gaussian", "mefixedeffect"),
        control.data = list(list(
                hyper = list(
                        prec = list(
                                initial = log(prec.y),
                                fixed=TRUE)
                        )
                ), 
                list()), 
        verbose=TRUE)

