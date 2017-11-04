rloglogistic = function(n,  lambda,  alpha, variant=0)
{
    u = runif(n)
    if (variant == 0) {
        y = (lambda/(1.0/u - 1.0))^(1.0/alpha)
    } else if (variant == 1) {
        y = (1.0/(1.0/u -1.0))^(1.0/alpha) / lambda
    } else {
        stop("ERROR")
    }
}
    
n = 1000
alpha = 2.1
x = c(scale(runif(n)))
eta = 1.1+2.2*x
lambda = exp(eta)

for(variant in 0:1) {

    print(paste("variant=", variant))
    y = rloglogistic(n, lambda = lambda,
                     alpha = alpha,
                     variant = variant)

    formula = y ~ 1 + x
    r=inla(formula,
           family ="loglogistic",
           data=data.frame(y, x), 
           control.family = list(variant = variant))
    print("REGRESSION")
    print(summary(r))

    event = rep(1,n)
    formula=inla.surv(y,event) ~ 1 + x
    r=inla(formula,
           family ="loglogisticsurv",
           data = list(y=y, event=event, x=x), 
           control.family = list(variant = variant))
    print("SURVIVAL")
    print(summary(r))
}
