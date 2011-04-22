rloglogistic = function(n,  beta,  alpha = 1)
{
    p = runif(n)
    return (beta* (((1-p)/p)^(-1/alpha)))
}
    
n = 1000
alpha = 2
x = runif(n)
eta = 1+x
beta = exp(eta)
y = rloglogistic(n, beta = beta, alpha = alpha)
event = rep(1,n)
data = list(y=y, event=event, x=x)
formula=inla.surv(y,event) ~ x
r=inla(formula, family ="loglogistic", data=data, verbose=T)
