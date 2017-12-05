n = 1000
alpha = 1.1
beta = 2.2
x = c(scale(runif(n)))
eta = 1+beta*x
lambda = exp(eta)

for(variant in 0:1) {
    y = rweibull(n,
                 shape= alpha,
                 scale= if (variant == 0)
                            lambda^(-1/alpha)
                        else
                            1/lambda)

    print(paste("VARIANT=", variant))
    event = rep(1,n)
    data = list(y=y, event=event, x=x)

    formula=inla.surv(y,event)~ x
    r=inla(formula,
           family ="weibullsurv",
           data=data,
           control.family = list(list(variant = variant)))
    print("SURV")
    print(summary(r))

    formula= y ~ x
    r=inla(formula,
           family ="weibull",
           data=data, 
           control.family = list(list(variant = variant)))
    print("REGRESSION")
    print(summary(r))
}
