nrep.max = 5
n = 50
y = matrix(NA, n, nrep.max)
x = c()
xx = c()
size = 3
overdispersion = 1/size
intercept = 1
    
for(i in 1:n) {
    local.x = runif(1) - 0.5
    lambda = exp(2 + local.x)
    N = rnbinom(1, mu=lambda, size=size)
    local.xx = runif(1) - 0.5
    eta = intercept + local.xx
    p = exp(eta)/(exp(eta) + 1)
    ## sample the number of replications
    nr = sample(1:nrep.max, 1)
    ## and sample these. note that 'y' is initialized with NA's,
    ## so if nr < nrep.max, then
    ## y[i,(nr+1):nrep.max]  would be NA.
    y[i, 1:nr]= rbinom(nr, size = N, prob = p)
    x = c(x, local.x)
    xx = c(xx, local.xx)
}

Y = inla.mdata(y, 1, x)
r = inla(Y ~ 1 + xx,
         data = list(Y=Y, xx=xx), 
         family = "nmixnb",
         control.fixed = list(prec.intercept=1,  prec=1))

