## ----setup, include=FALSE-------------------------------------------
library(INLA)
inla.setOption(smtp='taucs')
inla.setOption(num.threads="1:1")
inla.setOption(inla.mode="experimental")
if (file.exists("myinit.R")) source("myinit.R")
knitr::opts_chunk$set(fig.path="figures/old-faq/")
set.seed(123)
source("render_toc.R")

## ----toc, echo=FALSE------------------------------------------------
render_toc("old-faq.Rmd")

## -------------------------------------------------------------------
inla.setOption(inla.mode="classic")

## ----eval=FALSE-----------------------------------------------------
#     expression: statement; statement; ...; return(value)

## ----eval=FALSE-----------------------------------------------------
#     prior.expression = "expression:
#         a = 1;
#         b = 0.1;
#         precision = exp(log_precision);
#         logdens = log(b^a) - lgamma(a)
#            + (a-1)*log_precision - b*precision;
#         log_jacobian = log_precision;
#         return(logdens + log_jacobian);"

## ----eval=FALSE-----------------------------------------------------
#     table: x_1 ... x_n y_1 ... y_n

## -------------------------------------------------------------------
## the loggamma-prior
prior.function = function(log_precision) {
    a = 1;
    b = 0.1;
    precision = exp(log_precision);
    logdens = log(b^a) - lgamma(a) + (a-1)*log_precision - b*precision;
    log_jacobian = log_precision;
    return(logdens + log_jacobian)
}

## implementing the loggamma-prior using "expression:"
prior.expression = "expression:
a = 1;
b = 0.1;
precision = exp(log_precision);
logdens = log(b^a) - lgamma(a)
    + (a-1)*log_precision - b*precision;
log_jacobian = log_precision;
return(logdens + log_jacobian);"

## use suitable support points x
lprec = seq(-10, 10, len=100)
## link the x and corresponding y values into a 
## string which begins with "table:""
prior.table = paste(c("table:", cbind(lprec,
                    prior.function(lprec))), collapse=" ", sep="")

# simulate some data
n = 50
y = rnorm(n)

## use the built-in loggamma prior
r1 = inla(y~1,data = data.frame(y),
control.family = list(hyper = list(prec = list(
                      prior = "loggamma", param = c(1, 0.1)))))

## use the definition using expression
r2 = inla(y~1, data = data.frame(y), 
          control.family = list(hyper = list(
		                        prec = list(prior = prior.expression))))

## use a table of x and y values representing the loggamma prior
r3 = inla(y~1, data = data.frame(y), 
          control.family = list(hyper = list(
		                        prec = list(prior = prior.table))))

print(round(c(r1$mlik[1], r2$mlik[1], r3$mlik[1]), digits=3))

## -------------------------------------------------------------------
n = 100
z = rnorm(n)
eta = 1 + 0.1*z
N = 2

p = inla.link.invlogit(eta)
y = rbinom(n,  size = N, prob = p)
r = inla(y ~ 1 + z,  data = data.frame(y, z), family = "binomial", Ntrials = rep(N, n),
        control.family = list(control.link = list(model="logit")), 
		control.predictor = list(compute=TRUE))

p = inla.link.invprobit(eta)
y = rbinom(n,  size = N, prob = p)
rr = inla(y ~ 1 + z,  data = data.frame(y, z), family = "binomial", Ntrials = rep(N, n),
        control.family = list(control.link = list(model="probit")), 
		control.predictor = list(compute=TRUE))

p = inla.link.invcloglog(eta)
y = rbinom(n,  size = N, prob = p)
rrr = inla(y ~ 1 + z,  data = data.frame(y, z), family = "binomial", Ntrials = rep(N, n),
        control.family = list(control.link = list(model="cloglog")), 
		control.predictor = list(compute=TRUE))        

## -------------------------------------------------------------------
n = 100
n.pred = 10
y = arima.sim(n=n, model=list(ar=0.9))
N = n + n.pred
yy = c(y, rep(NA, n.pred))
i = 1:N
formula = yy ~ f(i, model="ar1")
r = inla(formula, data = data.frame(i,yy),
        control.family = list(initial = 10, fixed=TRUE)) ## no observational noise

## -------------------------------------------------------------------
r$summary.random$i[(n+1):N, c("mean", "sd") ]

## ----plot=TRUE------------------------------------------------------
## simple poisson regression
n = 100
x = sort(runif(n))
eta = 1 + x
lambda = exp(eta)
y = rpois(n, lambda = lambda)

## missing values:
y[1:3] = NA
y[(n-2):n] = NA

## link = 1 is a shortcut for rep(1, n) where n is the appropriate
## length. here '1' is a reference to the first 'family', ie
## 'family[1]'
r = inla(y ~ 1 + x,  family = "poisson",
        data = data.frame(y, x),
        control.predictor = list(link = 1))
plot(exp(eta),type ="l")
points(r$summary.fitted.values$mean, pch=19)

## ----plot=T---------------------------------------------------------
n2 = n %/% 2L
Y = matrix(NA, n, 2)
Y[1:n2, 1] = y[1:n2]
Y[1:n2 + n2, 2] = y[1:n2 + n2]
link = rep(NA, n)
link[which(is.na(y[1:n2]))] = 1
link[n2  + which(is.na(y[1:n2 + n2]))] = 2

r = inla(Y ~ 1 + x,  family = c("poisson", "nbinomial"),
        data = list(Y=Y, x=x),
		control.family= list(list(), list()),
        control.predictor = list(link = link))
plot(exp(eta),type ="l")
points(r$summary.fitted.values$mean, pch=19)

## -------------------------------------------------------------------
## Load the data
data(Tokyo)
summary(Tokyo)

Tokyo$y[300:366] <- NA

## Define the model
formula = y ~ f(time, model="rw2", scale.model=TRUE, 
	            constr=FALSE, cyclic=TRUE,
				hyper = list(prec=list(prior="pc.prec",
	                                   param=c(2,0.01)))) -1 

## We'll get a warning since we have not defined the link argument
result = inla(formula, family="binomial", Ntrials=n, data=Tokyo,
	control.compute = list(return.marginals.predictor = TRUE), 
	control.predictor=list(compute=T))

## need to recompute the fitted values for those with data[i] = NA,
## as the identity link is used.
n = 366
fitted.values.mean = numeric(n)
for(i in 1:366) {
    if (is.na(Tokyo$y[i])) {
        if (FALSE) {
            ## either like this, which is slower
            marg = inla.marginal.transform(
                            function(x) exp(x)/(1+exp(x)),
                            result$marginals.fitted.values[[i]] )
            fitted.values.mean[i] = inla.emarginal(function(x) x, marg)
        } else {
            ## or like this,  which is faster
            fitted.values.mean[i] = inla.emarginal(
                            function(x) exp(x)/(1 +exp(x)),
                            result$marginals.fitted.values[[i]])
        }
    } else {
        fitted.values.mean[i] = result$summary.fitted.values[i,"mean"]
    }
}
plot(fitted.values.mean)

## ----eval=FALSE-----------------------------------------------------
# 4
# 1 2 3 4
# 2 0
# 3 1 1
# 4 1 1

## ----plot=TRUE,fig.width=2,fig.align='center'-----------------------
g = inla.read.graph("4 1 2 3 4 2 0 3 1 1 4 1 1")
if (FALSE && requireNamespace("Rgraphviz") && requireNamespace("graph")) {
  tryCatch(
    plot(g),
    error = function(e) {
      message("Rgraphviz plotting doesn't seem to work on this installation.")
    }
  )
}

## ----eval=FALSE-----------------------------------------------------
# "4 1 2 3 4  2 0 3 1 1 4 1 1"

## ----eval=FALSE-----------------------------------------------------
# formula = y ~ f(idx, model = "besag", graph = "graph.dat")

## ----eval=FALSE-----------------------------------------------------
# formula = y ~ f(idx, model = "besag", graph = "4 1 2 3 4  2 0 3 1 1 4 1 1")

## -------------------------------------------------------------------
C = matrix(c(1, 0, 1, 1, 0, 1, 0, 0, 1, 0, 1, 0, 1, 0, 0, 1),4,4)
C

## -------------------------------------------------------------------
C.sparse= inla.as.sparse(C)
C.sparse

## -------------------------------------------------------------------
str(g)

## ----eval=FALSE-----------------------------------------------------
# formula = y ~ x + f(k, model= <some model>)

## ----eval=F---------------------------------------------------------
# control.fixed = list(expand.factor.strategy = "inla")

## -------------------------------------------------------------------
r = inla(y ~ 1 + x, data = data.frame(y=1:3, x=factor(c("a","b","c"))))
as.matrix(r$model.matrix)

## -------------------------------------------------------------------
r = inla(y ~ 1 + x, data = data.frame(y=1:3,x=factor(c("a","b","c"))), 
         control.fixed = list(expand.factor.strategy="inla"))
as.matrix(r$model.matrix)

## -------------------------------------------------------------------
r = inla(y ~ 1 + x, data = data.frame(y=1:3,x=factor(c("a","b",NA))), 
         control.fixed = list(expand.factor.strategy="inla"))
as.matrix(r$model.matrix)

## ----eval=FALSE-----------------------------------------------------
# control.compute=list(cpo=TRUE)

## ----eval=FALSE-----------------------------------------------------
# improved.result = inla.cpo(result)

## ----eval=FALSE-----------------------------------------------------
# r = inla(..., inla.call = "submit")

## ----eval=FALSE-----------------------------------------------------
# inla.qstat()

## ----eval=FALSE-----------------------------------------------------
# r = inla.qget(r)

## -------------------------------------------------------------------
## A simple model
n = 100
a = rnorm(n)
b = rnorm(n)
idx = 1:n

y = rnorm(n) + a + b
formula = y ~ 1 + a + b + f(idx, model="iid")

## assume we want to compute the posterior for
##
##  2 * beta_a + 3 * beta_b + idx[1] - idx[2]
##
## which we specify as follows (and giving it a unique name)

lc1 = inla.make.lincomb(a=2, b=3, idx = c(1,-1,rep(NA,n-2)))
names(lc1) = "lc1"

## strictly speaking, it is sufficient to use `idx = c(1,-1)', as the
## remaining idx's are not used in any case.

r = inla(formula, data = data.frame(a,b,y),
        ## add the linear combinations here
        lincomb = lc1,
        ## force noise variance to be essiatially zero
        control.family = list(initial=10, fixed=TRUE))
		
## to verify the result, we can compare the mean but the variance and
## marginal cannot be computed from the simpler marginals alone.
lc1.1 = 2 * r$summary.fixed["a", "mean"] + 3 * r$summary.fixed["b",
    "mean"] + r$summary.random$idx$mean[1] -
    r$summary.random$idx$mean[2]
lc1.2= r$summary.lincomb.derived$mean
print(round(c(lc1.1 = lc1.1, lc1.2 = lc1.2), dig=3))

## -------------------------------------------------------------------
## let wa and wb be vectors, and we want to compute the marginals for
## beta_a * wa[i] + beta_b * wb[i], for i=1..m. this is done
## conveniently as follows

m = 10
wa = runif(m)
wb = runif(m)
lc.many = inla.make.lincombs(a = wa, b=wb)

## we can give them names as well, but there are also default names, like
print(names(lc.many))
r = inla(formula, data = data.frame(a,b,y),
        lincomb = lc.many,
        control.family = list(initial=10, fixed=TRUE))
print(round(r$summary.lincomb.derived, dig=3))

## ----eval=FALSE-----------------------------------------------------
# control.inla = list(lincomb.derived.correlation.matrix = TRUE)

## ----eval=FALSE-----------------------------------------------------
# result$misc$lincomb.derived.correlation.matrix

## -------------------------------------------------------------------
n = 100
nPred = 10
phi = 0.9
x = arima.sim(n,  model = list(ar=phi)) * sqrt(1-phi^2)
y = 1 + x + rnorm(n,  sd=0.1)

time = 1:(n + nPred)
Y = c(y, rep(NA, nPred))
formula = Y ~ 1 + f(time, model="ar1")

## make linear combinations which are the nPred linear predictors
B = matrix(NA, nPred, n+nPred)
for(i in 1:nPred) {
    B[i, n+i] = 1
}
lcs = inla.make.lincombs(Predictor = B)

r = inla(formula,  data = data.frame(Y, time),
        control.predictor = list(compute=TRUE),
        lincomb = lcs,
		inla.mode="classic",
        control.inla = list(lincomb.derived.correlation.matrix=TRUE))

print(round(r$misc$lincomb.derived.correlation.matrix,dig=3))

## -------------------------------------------------------------------
n = 100
u = rnorm(n)
eta = 1 + u
p = exp(eta)/(1+exp(eta))
y = rbinom(n, size=1,  prob = p)

idx = 1:n
result = inla(y ~ 1 + f(idx, model="iid",
                        hyper = list(prec = list(prior="pc.prec",
						                         prior = c(1,0.01)))), 
              data =data.frame(y,idx), family = "binomial", 
	          Ntrials = 1)
summary(result)

## -------------------------------------------------------------------
## simple example
set.seed(12345L)
sd.data = 1
n = 2000

## fixed high and low precisions are specified with parameters
fixed.high <- 10
fixed.low <- - 6

## otherwise is the model identical to the previous FAQ-model
x1 = rnorm(n)
eta1 = 1 + x1
x2 = rnorm(n)
eta2 = 2 + 2*eta1 + 2*x2
y1 = rnorm(n, mean=eta1, sd = sd.data)
y2 = rnorm(n, mean=eta2, sd = sd.data)
## the trick is to create a vector 'u' (iid) which is equal to eta1, and then we can copy 'u' to
## create beta*u or beta*eta1. we do this by using 0 = eta1 -u + tiny.noise

formula = Y ~ -1 + intercept1 + X1 + intercept2 +
    f(u, w, model="iid", hyper = list(prec = list(initial = fixed.low, fixed=TRUE))) + 
    f(b.eta2, copy="u", hyper = list(beta = list(fixed = FALSE))) + X2
Y = matrix(NA, 3*n, 3)

## part 1: y1
intercept1 = rep(1, n)
X1 = x1
intercept2 = rep(NA, n)
u = rep(NA, n)
w = rep(NA, n)
b.eta2 = rep(NA, n)
X2 = rep(NA, n)
Y[1:n, 1] = y1
## part 2: 0 = eta1 - u + tiny.noise
intercept1 = c(intercept1,  intercept1)
X1 = c(X1,  x1)
intercept2 = c(intercept2,  rep(NA, n))
u = c(u, 1:n)
w = c(w, rep(-1, n))
b.eta2 = c(b.eta2,  rep(NA, n))
X2 = c(X2,  rep(NA, n))
Y[n + 1:n, 2] = 0
## part 3: y2
intercept1 = c(intercept1,  rep(NA, n))
X1 = c(X1,  rep(NA, n))
intercept2 = c(intercept2, rep(1, n))
u = c(u, rep(NA, n))
w = c(w, rep(NA, n))
b.eta2 = c(b.eta2, 1:n)
X2 = c(X2, x2)
Y[2*n + 1:n, 3] = y2
r = inla(formula,
    data =  list(Y=Y, intercept1=intercept1, X1=X1,
        intercept2=intercept2, u=u, w=w, b.eta2=b.eta2, X2=X2),
    family = c("gaussian", "gaussian", "gaussian"),
    control.family = list(
        list(),
        list(hyper = list(prec = list(initial = fixed.high, fixed=TRUE))),
        list()))
summary(r)

## An alternative specification of the same model applies the inla.stack() function. 
## It is called four times, one time for each submodel, one time for the "copying", 
## and one time for joining the three substacks to one.

## Part 1:
## The first model stack comes in two versions

mod1a.stack <- inla.stack(
	data = list(Y = as.matrix(cbind(y1, NA, NA))), 
	A = list(1, 1),
	effects = list(intercept1 = rep(1,n),
			x1 = x1)
)
		
mod1b.stack <- inla.stack(
	data = list(Y = as.matrix(cbind(y1, NA, NA))), 
	A = list(1),
	effects = list(u = 1:n)
)

## Part 2:		
copy.stack <- inla.stack(
	data = list(Y = as.matrix(cbind(NA, rep(0,n), NA))),
	A = list(1, 1, - 1),
	effects = list(intercept1 = rep(1,n),
			x1 = x1,
			u = 1:n)
)

## Part 3:
mod2.stack <- inla.stack(
	data = list(Y = as.matrix(cbind(NA, NA, y2))), 
	A = list(1, 1, 1),
	effects = list(intercept2 = rep(1,n),
			x2 = x2,
			b.eta2 = 1:n)
)

all.a.stack <- inla.stack(mod1a.stack, copy.stack, mod2.stack)
all.b.stack <- inla.stack(mod1b.stack, copy.stack, mod2.stack)

## The formula is slightly modified because the long vectors X1 and X2 has not been created
## and there is no weight applied as "-1" appears as the effect of u

formula = Y ~ -1 + intercept1 + x1 + intercept2 + x2 +
    f(u, model = "iid", hyper = list(prec = list(initial = fixed.low, fixed = TRUE))) + 
    f(b.eta2, copy = "u", hyper = list(beta = list(fixed = FALSE)))

## The inla() function makes use of the stack
r.new.a = inla(formula,
    data =  inla.stack.data(all.a.stack),
    control.predictor = list(A = inla.stack.A(all.a.stack)),
    family = c("gaussian", "gaussian", "gaussian"),
    control.family = list(
        list(),
        list(hyper = list(prec = list(initial = fixed.high, fixed=TRUE))),
        list()))
summary(r.new.a)

r.new.b = inla(formula,
    data =  inla.stack.data(all.b.stack),
    control.predictor = list(A = inla.stack.A(all.b.stack)),
    family = c("gaussian", "gaussian", "gaussian"),
    control.family = list(
        list(),
        list(hyper = list(prec = list(initial = fixed.high, fixed=TRUE))),
        list()))
summary(r.new.b)

## All specifications give similar and acceptable results both for fixed effects and precisions,
## The version r.new.b seems faster than r.new.a for more complex models with one or more spde-terms


## -------------------------------------------------------------------
inla.list.models()

## ----eval=FALSE-----------------------------------------------------
# y = a + b*w + ...

## ----eval=FALSE-----------------------------------------------------
# f(idx, model="iid2d", n=2*m, ...)

## ----eval=FALSE-----------------------------------------------------
# idx = 1:m
# idx.copy = m + 1:m
# formula = y ~  f(idx, model="iid2d", n=2*m) + f(idx.copy, w, copy="idx") + ....

## -------------------------------------------------------------------
n=1000
i=1:n
j = i
z = rnorm(n)
w = runif(n)
y = z  + 2*z*w + rnorm(n)
formula = y ~ f(i, model="iid",initial=0, fixed=T) +
              f(j, w, copy="i", fixed=FALSE)
r = inla(formula, data = data.frame(i,j,w,y))
summary(r)

## ----eval=FALSE-----------------------------------------------------
# beta = low + (high-low)*exp(beta.local)/(1+exp(beta.local))

## ----eval=FALSE-----------------------------------------------------
# f(idx, model = ..,  replicate = r)

## ----eval=FALSE-----------------------------------------------------
# i = 1:n
# formula = y ~  f(i, model = "iid") + ...

## ----eval=FALSE-----------------------------------------------------
# i = rep(1,n)
# r = 1:n
# formula = y ~  f(i, model="iid", replicate = r) + ...

## -------------------------------------------------------------------
n = 100
y1 = arima.sim(n=n, model=list(ar=c(0.9)))+10
y2 = arima.sim(n=n, model=list(ar=c(0.9)))+20
y3 = arima.sim(n=n, model=list(ar=c(0.9)))+30

formula = y ~ mean -1 + f(i, model="ar1", replicate=r)
y = c(y1,y2,y3)
i = rep(1:n, 3)
r = rep(1:3, each=n)
mean = as.factor(r)
result = inla(formula, family = "gaussian",
              data = data.frame(y, i, mean),
              control.family = list(initial = 12, fixed=TRUE))
summary(result)

## ----eval=FALSE-----------------------------------------------------
# y ~ a + 1

## -------------------------------------------------------------------
##  Simple linear regression with observations with two different
##  variances.
n = 100
N = 2*n
y = numeric(N)
x = rnorm(N)

y[1:n] = 1 + x[1:n] + rnorm(n, sd = 1/sqrt(1))
y[1:n + n] = 1 + x[1:n + n] + rnorm(n, sd = 1/sqrt(2))

Y = matrix(NA, N, 2)
Y[1:n, 1] = y[1:n]
Y[1:n + n, 2] = y[1:n + n]

formula = Y ~ x + 1
result = inla(
        formula,
        data = list(Y=Y, x=x),
        family = c("gaussian", "gaussian"),
        control.family = list(list(prior = "flat", param = numeric()),
                            list()))
summary(result)

## -------------------------------------------------------------------
## Simple example with two types of likelihoods
n = 10
N = 2*n

## common covariates
x = rnorm(n)

## Poisson, depends on x
E1 = runif(n)
y1 = rpois(n, lambda = E1*exp(x))

## Binomial, depends on x
size = sample(1:10, size=n, replace=TRUE)
prob = exp(x)/(1+exp(x))
y2 = rbinom(n, size= size, prob = prob)

## Join them together
Y = matrix(NA, N, 2)
Y[1:n, 1] = y1
Y[1:n + n, 2] = y2

## The E for the Poisson
E = numeric(N)
E[1:n] = E1
E[1:n + n] = NA

## Ntrials for the Binomial
Ntrials = numeric(N)
Ntrials[1:n] = NA
Ntrials[1:n + n] = size

## Duplicate the covariate which is shared
X = numeric(N)
X[1:n] = x
X[1:n + n] = x

## Formula involving Y as a matrix
formula = Y ~ X - 1

## `family' is now
result = inla(formula,
        family = c("poisson", "binomial"),
		control.family=list(list(), list()),
        data = list(Y=Y, X=X),
        E = E, Ntrials = Ntrials)
summary(result)

## ----eval=FALSE-----------------------------------------------------
# X = numeric(N)
# X[1:n] = x
# X[1:n + n] = NA
# 
# XX = numeric(N)
# XX[1:n] = NA
# XX[1:n + n] = xx
# 
# formula = Y  ~ X + XX -1

## -------------------------------------------------------------------
## An example with three independent AR(1)'s with separate means, but
## with the same hyperparameters. These are observed with three
## different likelihoods.

n = 100
x1 = arima.sim(n=n, model=list(ar=c(0.9))) + 0
x2 = arima.sim(n=n, model=list(ar=c(0.9))) + 1
x3 = arima.sim(n=n, model=list(ar=c(0.9))) + 2

## Binomial observations
Nt = 10 + rpois(n,lambda=1)
y1 = rbinom(n, size=Nt, prob = exp(x1)/(1+exp(x1)))

## Poisson observations
Ep = runif(n, min=1, max=10)
y2 = rpois(n, lambda = Ep*exp(x2))

## Gaussian observations
y3 = rnorm(n, mean=x3, sd=0.1)

## stack these in a 3-column matrix with NA's where not observed
y = matrix(NA, 3*n, 3)
y[1:n, 1] = y1
y[n + 1:n, 2] = y2
y[2*n + 1:n, 3] = y3

## define the model
r = c(rep(1,n), rep(2,n), rep(3,n))
rf = as.factor(r)
i = rep(1:n, 3)
formula = y ~ f(i, model="ar1", replicate=r, constr=TRUE) + rf -1
data = list(y=y, i=i, r=r, rf=rf)

## parameters for the binomial and the poisson
Ntrials = rep(NA, 3*n)
Ntrials[1:n] = Nt
E = rep(NA, 3*n)
E[1:n + n] = Ep

result = inla(formula, family = c("binomial", "poisson", "normal"),
              data = data, Ntrials = Ntrials, E = E,
              control.family = list(
                      list(),
                      list(),
                      list()))
summary(result)

## ----eval=FALSE-----------------------------------------------------
# y ~ 1 + z

## ----eval=FALSE-----------------------------------------------------
# eta~ = A %*% eta

## ----eval=FALSE-----------------------------------------------------
# y ~ 1 + z, with addition matrix A

## ----eval=FALSE-----------------------------------------------------
# y    ~ eta~   ## no intercept...
# eta~ = A %*% eta
# eta  = intercept + beta*z

## ----eval=FALSE-----------------------------------------------------
# control.predictor=list(A=A)

## ----eval=FALSE-----------------------------------------------------
# eta~ = A %*% eta + offset.arg
# eta = intercept + beta*z + offset.formula

## -------------------------------------------------------------------
## 'm' is the number of observations of eta*, where eta* = A eta +
## offset.arg, and A is a fixed m x n matrix, and eta has length n. An
## offset in the formula goes into 'eta' whereas an offset in the
## argument of the inla-call, goes into eta*
n = 10
m = 100
offset.formula = 10+ 1:n
offset.arg = 1 + 1:m

## a covariate
z = runif(n)

## the linear predictor eta
eta = 1 + z + offset.formula

## the linear predictor eta* = A eta + offset.arg.
A = matrix(runif(n*m), m, n);
##A = inla.as.sparse(A)  ## sparse is ok
## need 'as.vector', as 'Eta' will be a sparseMatrix if 'A' is sparse
## even if ncol(Eta) = 1
Eta = as.vector(A %*% eta) + offset.arg

s = 1e-6
Y = Eta + rnorm(m, sd=s)

## for a check, we can use several offsets. here, m1=-1 and p1=1, so
## they m1+p1 = 0.
r = inla(Y ~ 1+z + offset(offset.formula) + offset(m1) + offset(p1),
        ## The A-matrix defined here
        control.predictor = list(A = A, compute=TRUE, precision = 1e6),
        ## we need to use a list() as the different lengths of Y
        ## and z
        data = list(Y=Y, z=z,
                m1 = rep(-1, n),
                p1 = rep(1, n),
                offset.formula = offset.formula,
                offset.arg = offset.arg),
        ## this is the offset defined in the argument of inla
        offset = offset.arg,
        ##
        control.family = list(initial = log(1/s^2), fixed=TRUE))
summary(r)

## this should be a small number
print(max(abs(r$summary.linear.predictor$mean - c(Eta, eta))))

## ----eval=FALSE-----------------------------------------------------
# y = intercept + s[j] + 0.5*s[k] + noise

## ----plot=TRUE------------------------------------------------------
n = 100
s = c(-1, 0, 1)
nS = length(s)
j = sample(1L:nS, n, replace=TRUE)
k = j
k[j == 1L] = 2
k[j == 2L] = 3
k[k == 3L] = 1

noise = rnorm(n, sd=0.0001)
y = 1 + s[j] + 0.5*s[k] + noise

## build the formula such that the linear predictor is the intercept
## (index 1) and the 's' term (index 2:(n+1)). then kind of
## 'construct' the model using the A-matrix.
formula = y ~ -1 + intercept + f(idx)
A = matrix(0, n, nS+1L)
for(i in 1L:n) {
  A[i, 1L]        = 1
  A[i, 1L + j[i]] = 1
  A[i, 1L + k[i]] = 0.5
}

data = list(intercept = c(1, rep(NA, nS)), idx = c(NA, 1L:nS))
result = inla(formula, data=data, control.predictor=list(A=A))
## should be a straight line
plot(result$summary.random$idx$mean, s, pch=19)
abline(a=0,b=1)

