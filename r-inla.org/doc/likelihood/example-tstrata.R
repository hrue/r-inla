df = 10
n = 100L
nstrata = 5L
ntot = n * nstrata

z = rnorm(ntot)
y = numeric(ntot)
k = 0L
for(i in 1L:nstrata) {
    j = 1L:n
    stdev = i
    y[k + j] = 1 + z[k+j] + rt(n, df=df) / sqrt(df/(df-2)) * stdev
    k = k + n
}

strata = rep(1L:nstrata, each = n)
i = 1L:ntot
formula = y ~ 1 + z

r = inla(formula,  data = data.frame(y, z, strata), family = "tstrata", 
        strata = strata)
