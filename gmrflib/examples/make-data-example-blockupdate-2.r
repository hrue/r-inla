n = 10
N = n^2
g0 = function(x) return ((x-5.0)/10.0)
g1 = function(x) return (g0(x)^2)
E = 2
mu = 1.0

y = numeric(N)
g0idx = numeric(N)
g1idx = numeric(N)
k=1
for(i in 1:n)
{
    for(j in 1:n)
    {
        y[k] = rpois(1,lambda = E*exp(mu + g0(i) + g1(j) + rnorm(1, sd=0.1)))
        g0idx[k] = i
        g1idx[k] = j
        k = k + 1
    }
}
## convert the indices to base-0.
res = cbind(E,y,g0idx-1,g1idx-1)
write(t(res),file="example-blockupdate-2.dat", ncol=4)



