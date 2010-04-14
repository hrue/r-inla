n = 100
N = 10*n
z = rnorm(N)
rate = exp(1+z)
y = rexp(N, rate=rate)


trunc = 0.2
j = (y > trunc)
y = y[j]
z = z[j]
y = y[1:n]
z= z[1:n]

k = length(y)
low = sort(y)[k*0.1]
upp = sort(y)[k*0.5]

low = low + runif(n)*0.1
upp = upp + runif(n)*0.1

if (any(upp < low))
    stop("upp < low")
    
event = numeric(n)
lower = rep(0,n)
upper = rep(0,n)
trunc = rep(trunc,n)

for(i in 1:n)
{
    if (runif(1) < 0.5)
        e = 1
    else
    {
        e = sample(c(0,2,3), size=1)
        f = 0.8
        if (e == 3)
        {
            lower[i] = max(trunc[i], f*y[i])
            upper[i] = y[i]/f
        }
        else if (e==2)
            upper[i] = y[i]/f
        else if (e==0)
            lower[i] = max(trunc[i],y[i]*f)
        else
            stop("oops")
    }
    event[i] = e
}

A = list(event= event, truncation = trunc, lower = lower, upper = upper, time = y)


inla.update()
d = list(z=z, A=A)
##res = inla(A ~ 1 + z, family = "exponential", data = d, inla.call="inla", verbose = TRUE, keep=TRUE)
res = inla(A ~ 1 + z, family = "weibull", data = d, inla.call="inla", verbose = TRUE, keep=TRUE)
##res = inla(A ~ 1 + z, family = "ps", data = d, inla.call="inla", verbose = TRUE, keep=TRUE)
##h = inla.hyperpar(res)

write(t(cbind(0:(n-1), as.data.frame(A))), ncolumns = 5, file="surv-data.dat")
write(t(cbind(0:(n-1), z)), ncolumns = 2, file="z.dat")
