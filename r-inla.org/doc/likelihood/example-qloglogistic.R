lam_loglogistic = function(yq, alpha, q, variant = 0) 
{
    if (variant == 0) {
        lambda = yq^alpha * (1/q-1)
    } else if (variant == 1) {
        lambda = 1/yq * (1/(1/q-1))^(1/alpha)
    } else
        stop("ERR")
    return (lambda)
}
    

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
    
n = 500
alpha = 2.1
x = c(scale(runif(n)))
eta = 1.1+2.2*x
yq = exp(eta)

for(variant in 0:1) {
    for(q in c(0.2, 0.8)) {

        print(paste("variant=", variant, "quantile=", q))
        lambda = lam_loglogistic(yq, alpha, q, variant=variant)
        y = rloglogistic(n,
                         lambda = lambda, 
                         alpha = alpha,
                         variant = variant)
        
        formula = y ~ 1 + x
        rr=inla(formula,
               family ="qloglogistic",
               data=data.frame(y, x), 
               control.family = list(list(variant = variant, control.link = list(quantile = q))))
        print("REGRESSION")
        print(summary(rr))
               
        event = rep(1,n)
        formula=inla.surv(y,event) ~ 1 + x
        r=inla(formula,
               family ="qloglogisticsurv",
               data = list(y=y, event=event, x=x), 
               control.family = list(list(variant = variant, control.link = list(quantile = q))))
        print("SURVIVAL")
        print(summary(r))
    }
}
