y = rnorm(1)

## using buildt-in prior
a = 1
b = 0.1
hyper = list(prec = list(prior = "loggamma", param = c(a, b)))
r = inla(y ~ 1,  data = data.frame(y),
        control.family = list(hyper = hyper))

## implementing the loggamma-prior using "expression:"
loggamma = "expression:
            a = 1;
            b = 0.1;
            precision = exp(log_precision);
            logdens = log(b^a) - lgamma(a)
                      + (a-1)*log_precision - b*precision;
            log_jacobian = log_precision;
            return(logdens + log_jacobian);"

hyper.new = list(prec = list(prior = loggamma))
r.new = inla(y ~ 1,  data = data.frame(y),
        control.family = list(hyper = hyper.new))

## and we verify that we get the same result...
print(r$summary.hyperpar[1,"mean"])
print(r.new$summary.hyperpar[1,"mean"])
