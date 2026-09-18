## ----setup, include=FALSE-------------------------------------------
set.seed(1234)
library(INLA)
inla.setOption(num.threads="1:1")
inla.setOption(smtp="taucs")
inla.setOption(inla.mode="classic")
if (file.exists("myinit.R")) source("myinit.R")
library(ggplot2)
library(sn)
library(knitr)
set.seed(123)
knitr::opts_chunk$set(echo = TRUE)
opts_chunk$set(size='small'
               , cache=FALSE
               , cache.path='cache/'
               , comment=NA
               , warning=FALSE
               , message=FALSE
               , fig.align='center'
               , fig.path='figures/jmarginal/'
               , fig.pos='H'
               , background='#ffffff'
               , results='hold'
               , eval=TRUE)

## ----eval=F---------------------------------------------------------
# inla.setOption(inla.mode="classic")

## ----eval=F---------------------------------------------------------
# inla.setOption(inla.mode="compact")

## -------------------------------------------------------------------
n = 100
x = rnorm(n, mean = 1, sd = 0.3)
xx = rnorm(n, mean = -2, sd = 0.3)
y = rpois(n, lambda = exp(x + xx))
r = inla(y ~ 1 + x + xx,
         data = data.frame(y, x, xx), 
	     family = "poisson")

## -------------------------------------------------------------------
selection = list(Predictor = 3:4, x = 1, xx = 1)

## -------------------------------------------------------------------
rs = inla(y ~ 1 + x + xx,
          data = data.frame(y, x, xx), 
	      family = "poisson",
	      inla.mode = "classic",
          control.compute = list(return.marginals.predictor = TRUE), 
	      control.predictor = list(compute = TRUE),
		  selection = selection)

## -------------------------------------------------------------------
#summary(rs$selection)
print(rs$selection)

## -------------------------------------------------------------------
ns = 10000
xs = inla.rjmarginal(ns, rs) ## or rs$selection
str(xs)

## -------------------------------------------------------------------
pairs(t(xs$samples[, 1:3000]), cex = 0.1)

## -------------------------------------------------------------------
hist(xs$samples["Predictor:3",], n = 300, prob = TRUE, 
     main = 'Histogram of Predictor:3', xlab = 'Samples 
     from the joint and linear predictor marginal (black straight line)')
lines(inla.smarginal(rs$marginals.linear.predictor[[3]]), lwd = 3)

## -------------------------------------------------------------------
dxs = inla.1djmarginal(jmarginal = rs$selection)
str(dxs)

## -------------------------------------------------------------------
ggplot(data = data.frame(y = xs$samples["Predictor:3",]), aes(y, after_stat(density), colour = "True")) +
  stat_density(alpha = .1) +
  geom_line(data = as.data.frame(dxs$`Predictor:3`), aes(x = x, y = y, colour = "Deterministic"))+
  labs(title= '"Marginal Results for Predictor:3"', x='x', y='Density') +
  scale_colour_manual("", 
                      breaks = c("True", "Deterministic"),
                      values = c("black", "red")) 

ggplot(data = data.frame(y = xs$samples["Predictor:4",]), aes(y, after_stat(density), colour = "True")) +
  stat_density(alpha = .1) +
  geom_line(data = as.data.frame(dxs$`Predictor:4`), aes(x = x, y = y, colour = "Deterministic"))+
  labs(title= '"Marginal Results for Predictor:4"', x='x', y='Density') +
  scale_colour_manual("", 
                      breaks = c("True", "Deterministic"),
                      values = c("black", "red")) 

ggplot(data = data.frame(y = xs$samples["x:1",]), aes(y, after_stat(density), colour = "True")) +
  stat_density(alpha = .1) +
  geom_line(data = as.data.frame(dxs$`x:1`), aes(x = x, y = y, colour = "Deterministic"))+
  labs(title= '"Marginal Results for x:1"', x='x', y='Density') +
  scale_colour_manual("", 
                      breaks = c("True", "Deterministic"),
                      values = c("black", "red")) 

ggplot(data = data.frame(y = xs$samples["xx:1",]), aes(y, after_stat(density), colour = "True")) +
  stat_density(alpha = .1) +
  geom_line(data = as.data.frame(dxs$`xx:1`), aes(x = x, y = y, colour = "Deterministic"))+
  labs(title= '"Marginal Results for xx:1"', x='x', y='Density') +
  scale_colour_manual("", 
                      breaks = c("True", "Deterministic"),
                      values = c("black", "red"))  

## -------------------------------------------------------------------
summary(rs$selection)

## -------------------------------------------------------------------
A = matrix(c(1,1,0,0,0,0,1,1), nrow = 2, ncol = 4, byrow = T)
rownames(A) <- c("Test1", "Test2")
A

## -------------------------------------------------------------------
m = inla.tjmarginal(jmarginal = rs, A = A)
m
class(m)


## -------------------------------------------------------------------
dxs.lin = inla.1djmarginal(jmarginal = m)
str(dxs.lin)

fun1 <- function(...) {Predictor[1]+Predictor[2]}
fun2 <- function(...) {x+xx}

xs.lin1 = inla.rjmarginal.eval(fun1, xs)
xs.lin2 = inla.rjmarginal.eval(fun2, xs)

ggplot(data = data.frame(y = xs.lin1[1, ]), aes(y, after_stat(density), colour = "True")) +
  stat_density(alpha = .1) +
  geom_line(data = as.data.frame(dxs.lin$Test1), aes(x = x, y = y, colour = "Deterministic"))+
  labs(title= '"Marginal Results for Lin1"', x='x', y='Density') +
  scale_colour_manual("", 
                      breaks = c("True", "Deterministic"),
                      values = c("black", "red")) 

ggplot(data = data.frame(y = xs.lin2[1, ]), aes(y, after_stat(density), colour = "True")) +
  stat_density(alpha = .1) +
  geom_line(data = as.data.frame(dxs.lin$Test2), aes(x = x, y = y, colour = "Deterministic"))+
  labs(title= '"Marginal Results for Lin2"', x='x', y='Density') +
  scale_colour_manual("", 
                      breaks = c("True", "Deterministic"),
                      values = c("black", "red")) 

## -------------------------------------------------------------------
summary(m)

## -------------------------------------------------------------------
fun.exp <- function(x) exp(x)

fun5 <- function(...) {exp(x)}
fun6 <- function(...) {exp(xx)}
fun7 <- function(...) {exp(x+xx)}

tdx <- inla.tmarginal(fun = fun.exp, marginal = dxs$`x:1`)
tdxx <- inla.tmarginal(fun = fun.exp, marginal = dxs$`xx:1`)
tdx.lin <- inla.tmarginal(fun = fun.exp, marginal = dxs.lin$Test2)

tx = inla.rjmarginal.eval(fun5, xs)
txx = inla.rjmarginal.eval(fun6, xs)
tx.lin = inla.rjmarginal.eval(fun7, xs)

ggplot(data = data.frame(y = tx[1, ]), aes(y, after_stat(density), colour = "True")) +
  stat_density(alpha = .1) +
  geom_line(data = as.data.frame(tdx), aes(x = x, y = y, colour = "Deterministic"))+
  labs(title= '"Marginal Results for exp(x:1)"', x='x', y='Density') +
  scale_colour_manual("", 
                      breaks = c("True", "Deterministic"),
                      values = c("black", "red")) 

ggplot(data = data.frame(y = txx[1, ]), aes(y, after_stat(density), colour = "True")) +
  stat_density(alpha = .1) +
  geom_line(data = as.data.frame(tdxx), aes(x = x, y = y, colour = "Deterministic"))+
  labs(title= '"Marginal Results for exp(xx:1)"', x='x', y='Density') +
  scale_colour_manual("", 
                      breaks = c("True", "Deterministic"),
                      values = c("black", "red")) 

ggplot(data = data.frame(y = tx.lin[1, ]), aes(y, after_stat(density), colour = "True")) +
  stat_density(alpha = .1) +
  geom_line(data = as.data.frame(tdx.lin), aes(x = x, y = y, colour = "Deterministic"))+
  labs(title= '"Marginal Results for exp(x:1+xx:1)"', x='x', y='Density') +
  scale_colour_manual("", 
                      breaks = c("True", "Deterministic"),
                      values = c("black", "red")) 

## -------------------------------------------------------------------
expx = inla.zmarginal(marginal = tdx, silent = TRUE)
expxx = inla.zmarginal(marginal = tdxx, silent = TRUE)
expx.lin = inla.zmarginal(marginal = tdx.lin, silent = TRUE)

exp.summaries = rbind(expx, expxx, expx.lin)
exp.summaries

