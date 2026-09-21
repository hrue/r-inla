## ----setup, include=FALSE-------------------------------------------
library(INLA)
inla.setOption(num.threads="1:1")
inla.setOption(smtp="taucs")
if (file.exists("myinit.R")) source("myinit.R")

library(knitr)
library(formatR)
set.seed(20202)
opts_chunk$set(echo = TRUE)
opts_chunk$set(size='small'
               , cache=FALSE
               , cache.path='cache/scale-model'
               , comment=NA
               , warning=FALSE
               , message=FALSE
               , fig.align='center'
               , fig.path='figures/scale-model/'
               , fig.pos='H')
library(spatstat)
library(mvtnorm)
library(lattice)
library(mgcv)
library(pixmap)
library(numDeriv)
library(fields)
set.seed(123)

## -------------------------------------------------------------------
formula = y ~  f(idx, model = "rw2", scale.model = TRUE, hyper=..)

## ----eval=F---------------------------------------------------------
# f(idx, model = "..", group = g, control.group = list(model = "rw2", scale.model = TRUE))

## -------------------------------------------------------------------
# Functions to calculate the marginal standard deviations of rw1 and rw2
ss.rw1.sdref = function(values,tau) { 
    n = length(values)
    stopifnot(!missing(n)) 
    stopifnot(n > 3) 
    y = NA
    
    formula = y ~ -1 + f(values, model = "rw1", constr = TRUE, 
                         values = values, diagonal = 1e-10, 
                         hyper = list(prec = list(initial = log(tau), fixed = TRUE))) 
    r = inla(formula, data = data.frame(y, values), family = "gaussian",
             control.family = list(fixed = TRUE), 
             control.compute = list(return.marginals = F), verbose = F)
    sd = r$summary.random$values$sd 
    return (sd) 
    }

ss.rw2.sdref = function(values,tau) { 
    n = length(values)
    stopifnot(!missing(n)) 
    stopifnot(n > 3)
    y = NA 
    
    A1 = matrix(1, n, 1) 
    A2 = matrix(0, n, 1) 
    for(i in 1:n) A2[i, ]= i 
    
    A = rbind(c(A1), c(A2)) 
    e = rep(0, 2) 
    extraconstr = list(A = A, e= e)
    
    formula = y ~ -1 + f( values, model="rw2", constr = FALSE,
                          values = values, diagonal = 1e-10, extraconstr = extraconstr,
                          hyper = list(prec = list(initial = log(tau), fixed = TRUE)))
    
    r = inla(formula, data = data.frame(y, values), family = "gaussian", 
             control.family = list(fixed = TRUE),
             control.compute = list(return.marginals = F), verbose=F)
    
    sd = r$summary.random$values$sd 
    return (sd) 
}

# The marginal standard deviations for n nodes
n = 100 
sd.rw1 = ss.rw1.sdref(1:n, 1) 
sd.rw2 = ss.rw2.sdref(1:n, 1)

## ----plot = TRUE, fig.show = 'hold', out.width = '40%', echo = FALSE----
plot(sd.rw1, xlab = "Node i", ylab = "Marginal standard deviation") 
plot(sd.rw2, xlab = "Node i", ylab = "Marginal standard deviation") 

## ----echo=TRUE, eval=F----------------------------------------------
# U = sqrt(b/qgamma(alpha, a, 1))

## ----echo=F---------------------------------------------------------
func.u = function(a, b, alpha, sigma.ref) { 
    upper.limit = sqrt(b)*sigma.ref/sqrt(qgamma(alpha, a, 1)) 
    return(upper.limit) 
}
a = 1 
b = 5*10^{-5} 
alpha = 0.001 
upper.limit = func.u(a, b, alpha, 1) 

## -------------------------------------------------------------------
formula = y ~  f(x, model = "rw2", scale.model = TRUE, hyper=...)


## ----echo = F-------------------------------------------------------
func <-function(x, t)
{
    f.x = x/t + sin(2*pi* x/t) - 0.5
    return(f.x)
}

f.est <- function(y, x, a, b, option)
{
    formula = y ~ 1 + f(x, model="rw2",
                    hyper = list(prec = list(prior = "loggamma", param = c(a,b))), 
                    scale.model = option)
    result = inla(formula,family = "gaussian", data = data.frame(y,x))
    f.est = result$summary.random$x$mean
    return(f.est)
}

## ----echo = F-------------------------------------------------------
### Same values of the function, same observations, but  defined on different intervals ###
set.seed(89236)
x2 = 0:n
x1 = x2/100
x3 = x2*10

sigma = 0.5
f.x = func(x2,length(x2))
y = f.x + rnorm(length(f.x), 0, sigma)

a = 1
b = 0.00005

mat1 = cbind(f.x, f.est(y, x1, a, b, T), f.est(y, x1, a, b, F))
mat2 = cbind(f.x, f.est(y, x2, a, b, T), f.est(y, x2, a, b, F))
mat3 = cbind(f.x, f.est(y, x3, a, b, T), f.est(y, x3, a, b, F))

# Generalized marginal variances  
v1 = exp(mean(log(ss.rw2.sdref(x1, 1)^2)))
v2 = exp(mean(log(ss.rw2.sdref(x2, 1)^2)))
v3 = exp(mean(log(ss.rw2.sdref(x3, 1)^2)))

u1 = func.u(a, b, alpha,sqrt(v1))
u2 = func.u(a, b, alpha,sqrt(v2))
u3 = func.u(a, b, alpha,sqrt(v3))

## ----eval=F---------------------------------------------------------
# result = inla(formula, family = "gaussian", data = data.frame(y, x))

## ----eval=F---------------------------------------------------------
# result$summary.random$x$mean

## ----plot=TRUE, echo = FALSE, fig.show = 'hold', out.width = '30%'----
r=range(c(mat1,mat2,mat3))
matplot(x1, mat1, type="lll", lty=c(1,2,4),col=c(1,4,2),xlab="x",ylab="f(x)",ylim=r,lwd=3) 
matplot(x2, mat2, type="lll", lty=c(1,2,4),col=c(1,4,2),xlab="x",ylab="f(x)",ylim=r,lwd=3)
matplot(x3, mat3, type="lll", lty=c(1,2,4),col=c(1,4,2),xlab="x",ylab="f(x)",ylim=r,lwd=3)

## ----echo=TRUE, eval=F----------------------------------------------
# hyper.prec = list(prec = list(param = c(1, 0.01)))

## ----echo=TRUE, eval=F----------------------------------------------
# data(Munich)
# g = system.file("demodata/munich.graph", package="INLA")
# 
# ## Note that here we want to have an estimator of the effect of year
# ## also the for years where we have no observation, therefore we give a
# ## vector with all possible values assumed by the covariate year, that
# ## is seq(1918,2001)
# 
# formula = rent ~ -1 +
#    f(location, model = "besag", graph = g, initial = 1, hyper = hyper.prec,
#      scale.model = option) +
#    f(year, model = "rw2", values = seq(1918, 2001), hyper = hyper.prec,
#      scale.model = option) +
#    f(floor.size, model = "rw2", hyper = hyper.prec, scale.model = option) +
#    Gute.Wohnlage + Beste.Wohnlage + Keine.Wwv + Keine.Zh +
#    Kein.Badkach  + Besond.Bad + Gehobene.Kueche +
#    zim1 + zim2 + zim3 + zim4 + zim5 + zim6
# 

## ----echo=TRUE, eval=F----------------------------------------------
# mod = inla(formula, data = Munich, control.fixed=list(prec=1))

## ----echo=FALSE-----------------------------------------------------
data(Munich)
g = system.file("demodata/munich.graph", package="INLA")

func.munich <- function(a, b, option)
{
    hyper.prec = list(prec = list(param = c(a, b)))
    formula = rent ~ f(location, model = "besag", graph = g,initial=1, 
                       hyper = hyper.prec, scale.model = option) +
      f(year, model = "rw2", values = seq(1918, 2001), 
        hyper = hyper.prec, scale.model = option) +
      f(floor.size, model = "rw2", hyper = hyper.prec, scale.model = option) +
      Gute.Wohnlage + Beste.Wohnlage + Keine.Wwv + Keine.Zh + 
      Kein.Badkach  + Besond.Bad + Gehobene.Kueche + 
      zim1 + zim2 + zim3 + zim4 + zim5 + zim6 -1
    
    return (inla(formula, data = Munich, control.fixed = list(prec=1, prec.intercept = 1)))
}


## ----echo=FALSE-----------------------------------------------------
a = 1
b = 0.01
u.munich = func.u(a, b, alpha, 1)
mod = func.munich(a, b, FALSE)
mod.scaled = func.munich(a, b, TRUE)

## ----echo=F---------------------------------------------------------
alpha = 0.001
u.vec = c(0.001,0.1,10,30)
a.vec = c(1,1,1,1)
b.vec = u.vec^2*qgamma(alpha,a.vec,1)

option=TRUE
c1 = func.munich(a.vec[1], b.vec[1], option)
c2 = func.munich(a.vec[2], b.vec[2], option)
c3 = func.munich(a.vec[3], b.vec[3], option)
c4 = func.munich(a.vec[4], b.vec[4], option)

## ----plot=TRUE, echo = F, fig.show = 'hold', out.width = '40%'------
mat = cbind(mod$summary.random$year$mean, mod.scaled$summary.random$year$mean)
matplot(mod$summary.random$year$ID, mat, type="ll", lty=c(2,1), col=c(2,1), 
        xlab="Year of construction", ylab="", lwd=2.5)

mat = cbind(mod$summary.random$floor.size$mean, mod.scaled$summary.random$floor.size$mean)
matplot(mod$summary.random$floor.size$ID, mat, type="ll", lty=c(2,1), col=c(2,1),
        xlab="Floor size",ylab="",lwd=2.5)

## ----plot=TRUE, fig.show = 'hold', out.width = '40%', echo = FALSE----
mat = cbind(c1$summary.random$year$mean,
            c2$summary.random$year$mean,
            c3$summary.random$year$mean,
            c4$summary.random$year$mean)
matplot(mod$summary.random$year$ID,mat, type="llll", lty=c(1,2,3,4), col=c(1,2,3,4),
        xlab="Year of construction", ylab="", lwd=2.5)
legend("topleft", paste("U= ",u.vec), lty=1:4, col=1:4, lwd=2.5)

mat = cbind(c1$summary.random$floor.size$mean,
            c2$summary.random$floor.size$mean,
            c3$summary.random$floor.size$mean,
            c4$summary.random$floor.size$mean)
matplot(mod$summary.random$floor.size$ID,mat, type="llll", lty=c(1,2,3,4), col=c(1,2,3,4),
        xlab="Floor size", ylab="", lwd=2.5)
legend("topright", paste("U= ",u.vec), lty=1:4, col=1:4, lwd=2.5)

## -------------------------------------------------------------------
inla.setOption(scale.model.default = TRUE)

