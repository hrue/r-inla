## ----setup, include=FALSE-------------------------------------------
set.seed(123)
library(INLA)
inla.setOption(num.threads="1:1")
inla.setOption(smtp="taucs")
inla.setOption(inla.mode="classic")
if (file.exists("myinit.R")) source("myinit.R")
knitr::opts_chunk$set(echo = TRUE, cache = FALSE)
knitr::opts_chunk$set(fig.path="figures/conditional-sampling/")

## -------------------------------------------------------------------
inla.setOption(inla.mode="classic")

## -------------------------------------------------------------------
n = 50
grp.len = 5
x = rnorm(n)
g = rnorm(grp.len)
grp = rep(1:grp.len, each = n %/% grp.len)
y = 1 + x + g[grp] + rnorm(n, sd = 0.01)

## -------------------------------------------------------------------
r = inla(y ~ x + f(grp, model = "iid"),
         data = data.frame(y, x, grp), 
         control.compute = list(config=TRUE)) 

## -------------------------------------------------------------------
r$misc$configs$contents

## -------------------------------------------------------------------
m = sum(r$misc$configs$contents$length)
grp.idx = which(r$misc$configs$contents$tag == "grp")
grp.len = r$misc$configs$contents$length[grp.idx]
A = matrix(0, grp.len, m)
e = matrix(0, grp.len, 1)
for(i in 1:grp.len) {
    A[i, r$misc$configs$contents$start[grp.idx] + i - 1] = 1
    e[i] = 0
}

## -------------------------------------------------------------------
constr = r$misc$configs$constr
if (is.null(constr)) {
    ## nothing there from before
    r$misc$configs$constr = list(
        nc = dim(A)[2],
        A = A,
        e = e)
} else {
    ## create a new one
    r$misc$configs$constr = list(
        nc = constr$nc + dim(A)[2],
        A = rbind(constr$A, A),
        e = rbind(constr$e, e))
}

## -------------------------------------------------------------------
xx = inla.posterior.sample(1000, r)

