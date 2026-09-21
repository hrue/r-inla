## ----setup, include=FALSE-------------------------------------------
set.seed(123)
library(knitr)
knitr::opts_chunk$set(echo = TRUE, results = 'show', 
	fig.keep = 'high', fig.align = 'center', fig.pos = 'h')
knitr::opts_chunk$set(fig.path="figures/bgev/")

library(INLA)
inla.setOption(smtp="taucs")
inla.setOption(num.threads="1:1")
inla.setOption(inla.mode="experimental")
if (file.exists("myinit.R")) source("myinit.R")

## ----echo=F, message=F, warning=FALSE-------------------------------
colFmt = function(x,color){
  outputFormat = opts_knit$get("rmarkdown.pandoc.to")
  if(outputFormat == 'latex')
    paste("\\textcolor[HTML]{",color,"}{",x,"}",sep="")
  else if(outputFormat == 'html')
    paste("<font color='",color,"'>",x,"</font>",sep="")
  else
    x
}

library(INLA)
library(evd)
giveme.gev.par = function(q, sbeta, alpha, beta, xi) 
{
  .mu = function(q, sbeta, alpha, beta, xi) {
    a = -log(1-beta/2)
    b = -log(beta/2)
    c = -log(alpha)
    if (all(xi > 0.0)) {
      tmp0 = (c^(-xi) - 1)/xi
      tmp1 = a^(-xi)
      tmp2 = b^(-xi)
      dbeta = (tmp1 - tmp2)/xi
      return(q - (sbeta/dbeta) * tmp0)
    } else if (all(xi == 0.0)) {
      dbeta = log(b) - log(a)
      tmp0 = log(c)
      return(q + (sbeta/dbeta) * tmp0)
    } else {
      stop("mixed case not implemented")
    }
  }
  
  .sigma = function(q, sbeta, alpha, beta, xi) {
    a = -log(1-beta/2)
    b = -log(beta/2)
    if (all(xi > 0.0)) {
      tmp1 = a^(-xi)
      tmp2 = b^(-xi)
      dbeta = (tmp1 - tmp2)/xi
      return(sbeta/dbeta)
    } else if (all(xi == 0.0)) {
      dbeta = log(b) - log(a)
      return(sbeta/dbeta)
    } else {
      stop("mixed case not implemented")
    }
  }
  
  return(list(mu = .mu(q, sbeta, alpha, beta, xi),
              sigma = .sigma(q, sbeta, alpha, beta, xi),
              xi = xi))
}

map.tail = function(x, interval, inverse = FALSE) {
  if (!inverse) {
    return (interval[1] + (interval[2] - interval[1]) * exp(x)/(1.0 + exp(x)))
  } else {
    return (log((x-interval[1])/(interval[2]-x)))
  }
}

## ----pressure, echo=FALSE, fig.cap="\\label{fig:distributionH.pdf} bGEV distribution (H, black) constructed from distributions F (Frechét, red), G (Gumbel, green) and Beta weight function $p$ (purple). The shaded area is the mixing area, where F and G are merged.", out.width = '50%'----
knitr::include_graphics("fig_distribution1.png")

## -------------------------------------------------------------------
n = 1000
x = rnorm(n, sd=0.5) # we generate values for x from a N(0,0.5^2) dist.
eta.x = 1 + 0.4*x

## -------------------------------------------------------------------
spread = 0.3
tail = 0.1

## -------------------------------------------------------------------
p.alpha = 0.5
p.beta = 0.25

## -------------------------------------------------------------------
par = giveme.gev.par(q = eta.x, sbeta = spread, alpha = p.alpha, beta = p.beta, 
                     xi = tail)
y = numeric(n)
for(i in 1:n) 
  y[i] = rgev(1, loc = par$mu[i], scale = par$sigma, shape = par$xi)

## -------------------------------------------------------------------
hyper.spread = list(initial = 1,
                    fixed=FALSE,
                    prior = "loggamma",
                    param = c(3, 3))

## -------------------------------------------------------------------
tail.interval = c(0, 0.5)
tail.intern = map.tail(tail, tail.interval, inverse=TRUE)

## -------------------------------------------------------------------
hyper.tail = list(initial = tail.intern, 
                  prior = "pc.gevtail",
                  param = c(7, tail.interval), 
                  fixed= FALSE)


## -------------------------------------------------------------------
hyper.tail = list(initial = if (tail == 0.0) -Inf else tail.intern, 
                  prior = "pc.gevtail",
                  param = c(7, tail.interval), 
                  fixed= if (tail == 0.0) TRUE else FALSE)


## -------------------------------------------------------------------

hyper.bgev = list(spread = hyper.spread, 
                  tail = hyper.tail)

## -------------------------------------------------------------------
control.bgev = list(q.location = p.alpha,
                    q.spread = p.beta,
                    # quantile levels for the mixing part
                    q.mix= c(0.05, 0.20), 
                    # the Beta(s1, s2) mixing distribution parameters. 
                    # Hard-coded for the moment: s1=s2=5
                    beta.ab = 5)

## -------------------------------------------------------------------
null.matrix = matrix(nrow = n, ncol= 0)
spread.x = null.matrix
tail.x = null.matrix

## -------------------------------------------------------------------
data.bgev = data.frame(y = y, intercept = 1, x = x, spread.x = spread.x, tail.x = tail.x)
formula = inla.mdata(y, spread.x, tail.x) ~ -1 + intercept + x

## ----message=F, eval = T, warning=FALSE-----------------------------
r1 = inla(formula,
         family = "bgev",
         data = data.bgev,
         control.family = list(hyper = hyper.bgev,
                               control.bgev = control.bgev),
         control.predictor = list(compute = TRUE),
	     control.fixed = list(prec=1000),
         control.compute = list(cpo = TRUE),
         control.inla = list(int.strategy = "eb"),
         verbose=FALSE, safe=TRUE)

## -------------------------------------------------------------------
round(r1$summary.fixed,4)
round(r1$summary.hyperpar,4)

## -------------------------------------------------------------------
n = 1000
x1 = rnorm(n)
eta.x = 1 + 0.4*x1

## -------------------------------------------------------------------
x2 = rnorm(n, sd= 0.2)
s.x = exp(0.1 + 0.3*x2)
x3 = runif(n,-0.25,1)
t.x = 0.1 + 0.2*x3
tail.intern = map.tail(t.x, tail.interval, inverse=TRUE) # internal xi

## ----eval = T-------------------------------------------------------
par = giveme.gev.par(q = eta.x, sbeta = s.x, alpha = p.alpha, beta = p.beta, 
                     xi = t.x)
y = numeric(n)
for(i in 1:n) 
  y[i] = rgev(1, loc = par$mu[i], scale = par$sigma[i], shape = par$xi[i])

## -------------------------------------------------------------------
hyper.beta1 = hyper.beta2 = list(prior = "normal",
                                 param = c(0, 300),
                                 initial = 0)

## -------------------------------------------------------------------

hyper.bgev = list(spread = hyper.spread, 
                  tail = hyper.tail,
                  beta1 = hyper.beta1,
                  beta2 = hyper.beta2)

## -------------------------------------------------------------------
spread.x = x2
tail.x = x3
formula = inla.mdata(y, spread.x, tail.x) ~ -1 + intercept + x

## -------------------------------------------------------------------
data.bgev = data.frame(y = y, intercept = 1, x = x, spread.x = spread.x, tail.x = tail.x)

## ----message=F, eval = T, warning=FALSE-----------------------------
r2 = inla(formula,
         family = "bgev",
         data = data.bgev,
         control.family = list(hyper = hyper.bgev,
                               control.bgev = control.bgev),
         control.predictor = list(compute = TRUE),
	     control.fixed = list(prec=1000),
         control.compute = list(cpo = TRUE),
         control.inla = list(int.strategy = "eb"),
         verbose=FALSE, safe=TRUE)

## ----eval = T-------------------------------------------------------
round(r2$summary.fixed,4)
round(r2$summary.hyperpar,4)

## ----eval = T-------------------------------------------------------
n = 1000
x = rnorm(n)
z1 = seq(0, 6, length.out = n)
z2 = 1:n
p = 2 # AR order
pacf = runif(p)
phi = inla.ar.pacf2phi(pacf)
eta.x = 1 + 0.4*x + sin(z1) + c(scale(arima.sim(n, model = list(ar = phi))))

## -------------------------------------------------------------------
x2 = rnorm(n, sd = 0.2)
x4 = rnorm(n, sd = 0.2)
s.x = exp(0.1 + 0.3*x2 + x4)
x3 = runif(n,-0.2, 0.2)
t.x = 0.1 + 0.2*x3
tail.intern = map.tail(t.x, tail.interval, inverse=TRUE) # internal xi

## ----eval = T-------------------------------------------------------
par = giveme.gev.par(q = eta.x, sbeta = s.x, alpha = p.alpha, beta = p.beta, 
                     xi = t.x)
y = numeric(n)
for(i in 1:n) 
  y[i] = rgev(1, loc = par$mu[i], scale = par$sigma[i], shape = par$xi[i])

## -------------------------------------------------------------------
hyper.beta1 = hyper.beta2 = hyper.beta3 = list(prior = "normal",
                                               param = c(0, 300),
                                               initial = 0)
hyper.bgev = list(spread = hyper.spread, 
                  tail = hyper.tail,
                  beta1 = hyper.beta1,
                  beta2 = hyper.beta2,
                  beta3 = hyper.beta3)

## ----message=F, eval = T, warning=FALSE-----------------------------
spread.x = x2
spread.xx = x4
tail.x = x3
# With this change of variable it is easier to keep track of the effect 
# of the covariates in each parameter, but it is not needed.
formula = inla.mdata(y, cbind(spread.x, spread.xx), tail.x) ~ -1 + intercept + x + 
	f(z1, model = "rw1", scale.model=TRUE, constr=TRUE,
		hyper = list(prec = list(prior = "pc.prec", 
		param = c(0.1, 0.01)))) +
    f(z2, model = 'ar', order = 1, 
		   hyper=list(prec=list(prior="pc.prec", constr=TRUE,
		                        param=c(0.1,0.01)),
                      pacf1=list(param=c(0.5,0.8),
                      pacf2=list(param=c(0.5,0.8)))))
data.bgev = data.frame(y = y, intercept = 1, x = x, z1 = z1, z2 = z2,
                       spread.x = spread.x, spread.xx = spread.xx, tail.x = tail.x)
r3 = inla(formula,
         family = "bgev",
         data = data.bgev,
         control.family = list(hyper = hyper.bgev,
                               control.bgev = control.bgev),
         control.predictor = list(compute = TRUE),
	     control.fixed = list(prec=1000,prec.intercept=1000),
         control.compute = list(cpo = TRUE),
         control.inla = list(int.strategy = "eb",
	                         cmin=0,
	                         b.strategy="keep"),
         verbose=FALSE, safe=TRUE)

## ----eval = T-------------------------------------------------------
round(r3$summary.fixed,4)
round(r3$summary.hyperpar,4)

## ----eval=F---------------------------------------------------------
# library(INLA)
# inla.list.models('latent')

## -------------------------------------------------------------------
library(evd)
giveme.gev.par = function(q, sbeta, alpha, beta, xi) 
{
  .mu = function(q, sbeta, alpha, beta, xi) {
    a = -log(1-beta/2)
    b = -log(beta/2)
    c = -log(alpha)
    if (all(xi > 0.0)) {
      tmp0 = (c^(-xi) - 1)/xi
      tmp1 = a^(-xi)
      tmp2 = b^(-xi)
      dbeta = (tmp1 - tmp2)/xi
      return(q - (sbeta/dbeta) * tmp0)
    } else if (all(xi == 0.0)) {
      dbeta = log(b) - log(a)
      tmp0 = log(c)
      return(q + (sbeta/dbeta) * tmp0)
    } else {
      stop("mixed case not implemented")
    }
  }
  
  .sigma = function(q, sbeta, alpha, beta, xi) {
    a = -log(1-beta/2)
    b = -log(beta/2)
    if (all(xi > 0.0)) {
      tmp1 = a^(-xi)
      tmp2 = b^(-xi)
      dbeta = (tmp1 - tmp2)/xi
      return(sbeta/dbeta)
    } else if (all(xi == 0.0)) {
      dbeta = log(b) - log(a)
      return(sbeta/dbeta)
    } else {
      stop("mixed case not implemented")
    }
  }
  
  return(list(mu = .mu(q, sbeta, alpha, beta, xi),
              sigma = .sigma(q, sbeta, alpha, beta, xi),
              xi = xi))
}


## -------------------------------------------------------------------
map.tail = function(x, interval, inverse = FALSE) {
  if (!inverse) {
    return (interval[1] + (interval[2] - interval[1]) * exp(x)/(1.0 + exp(x)))
  } else {
    return (log((x-interval[1])/(interval[2]-x)))
  }
}

