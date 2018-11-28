library(mvtnorm)
n = 300
idx = 1:n
nstrata = 3
strata = sample(1:nstrata, n, replace=TRUE)
nsubject = n %/% nstrata
subject = sample(1:nsubject, n, replace=TRUE)
z = rnorm(n)
gam = c(1, 1 + rnorm(nstrata-1, sd = 0.2))

rho = sqrt(3)/2
Sigma = matrix(c(1/1, NA, NA, 1/2), 2, 2)
Sigma[1,2] = Sigma[2,1] = rho*sqrt(Sigma[1,1]*Sigma[2,2])

ab = rmvnorm(nsubject, sigma=Sigma)
a = ab[,1]
b = ab[,2]
s = 0.01
y = gam[strata] * (a[subject] + z * b[subject]) + rnorm(n, s = 0.01)
r = inla(y ~ -1 + f(idx, model = "intslope",
                    args.intslope = list(subject = subject,
                                         strata = strata,
                                         covariates = z)), 
         data = list(y = y,
                     idx = idx,
                     subject = subject,
                     strata = strata,
                     z = z), 
         control.family = list(hyper = list(
                                   prec = list(initial = log(1/s^2),
                                               fixed=TRUE))))
summary(r)
