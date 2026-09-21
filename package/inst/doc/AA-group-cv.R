## ----setup, include=FALSE-----------------------------------------------------
library(INLA)
library(knitr)
library(formatR)
library(fastGHQuad)
library(pracma)
inla.setOption(smtp="taucs")
inla.setOption(num.threads="1:1")
if (file.exists("myinit.R")) source("myinit.R")
set.seed(20202)
knitr::opts_chunk$set(echo = TRUE)
opts_chunk$set(size='small'
               , cache=FALSE
               , cache.path='cache/group-cv'
               , comment=NA
               , warning=FALSE
               , message=FALSE
               , fig.align='center'
               , fig.path='figures/group-cv/'
               , fig.pos='H'
               , background='#ffffff'
               , results='hold'
               , eval=TRUE)

# gh.x = c(-6.995680e+00,-6.275079e+00,-5.673961e+00,-5.133596e+00,-4.631560e+00,-4.156272e+00,-3.700743e+00,-3.260321e+00,-2.831680e+00,-2.412318e+00,-2.000259e+00,-1.593886e+00,-1.191827e+00,-7.928770e-01,-3.959427e-01,-2.722499e-16,3.959427e-01,7.928770e-01,1.191827e+00,1.593886e+00,2.000259e+00,2.412318e+00,2.831680e+00,3.260321e+00,3.700743e+00,4.156272e+00,4.631560e+00,5.133596e+00,5.673961e+00,6.275079e+00,6.995680e+00)
# gh.w = c(4.618968e-22,5.110609e-18,5.899556e-15,1.860374e-12,2.352492e-10,1.461199e-08,5.043713e-07,1.049860e-05,1.395209e-04,1.233683e-03,7.482800e-03,3.184723e-02,9.671795e-02,2.121328e-01,3.387727e-01,3.957786e-01,3.387727e-01,2.121328e-01,9.671795e-02,3.184723e-02,7.482800e-03,1.233683e-03,1.395209e-04,1.049860e-05,5.043713e-07,1.461199e-08,2.352492e-10,1.860374e-12,5.899556e-15,5.110609e-18,4.618968e-22)

gh.x = fastGHQuad::gaussHermiteData(51)$x
gh.w = fastGHQuad::gaussHermiteData(51)$w

inla.lik.pred = function(y,args = NULL,args.string,eta.mean,eta.sd,eta.skewness){
    res = dnorm(y,mean = eta.mean,sd = sqrt(eta.sd^2 + 1/exp(args$theta)))
    return(res)
}

## ----message=FALSE, warning=FALSE,echo = FALSE--------------------------------
n_rat = 10
n_each = 10
n = n_rat*n_each
id = rep(1:n_rat, each = n_each)
mu_rat = rep(x = rnorm(n_rat,mean = 30,sd = 1),each = n_each)
sd_measure = .1
y = rnorm(n = n, mean = mu_rat,sd = sd_measure)
plot(id,y,main = "Data against its group id")

## ----message = FALSE, warning=FALSE-------------------------------------------
data = data.frame(y = c(y,rep(NA,n_rat+1)),id = c(id,1:n_rat,n_rat+1))
formula = y ~ 1 + f(id,model = "iid")
res = inla(formula = formula,
           data = data,
           family = "normal"
           )

## ----message = FALSE, warning=FALSE, echo = TRUE, eval=FALSE------------------
# ##
# ## this example needs to be fixed so it does not depend on the removed
# ## inla.likelihood* feature
# ##
# res = inla(formula = formula,
#            data = data.frame(y = c(y,rep(NA,n_rat+1)), id = c(id,1:n_rat,n_rat+1)),
#            family = "normal",
#            control.compute = list(config = TRUE,likelihood.info = T)
#            )
# nconfig = res$misc$configs$nconfig
# weights = unlist(lapply(1:nconfig,FUN = function(i){res$misc$configs$config[[i]]$log.posterior}))
# weights = exp(weights)/sum(exp(weights))
# eps = 0.001
# q = seq(eps,1-eps,eps)
# preddensity = vector(mode = "list",length = n_rat+1)
# for(rat_idx in c(1,11)){
#     predictor_idx = n + rat_idx
#     mode_idx = which.max(weights)
#     args = INLA:::inla.likelihood.parser(res$misc$configs$config[[1]]$arg.str[1])
#     for(mode_idx in mode_idx){
#         args$theta = res$misc$configs$config[[mode_idx]]$theta[1]
#         mode_eta_mu = res$misc$configs$config[[mode_idx]]$Predictor[predictor_idx,"mean"]
#         mode_eta_sd = sqrt(res$misc$configs$config[[mode_idx]]$Predictor[predictor_idx,"variance"])
#         eta_eval = mode_eta_mu + seq(-5,5)*mode_eta_sd
#         args$linear.predictor = eta_eval
#         qfunc = inla.likelihood(type = "q",args)
#         y_eval = c()
#         for(i in 1:length(q)){
#             y_eval = c(y_eval,qfunc(q[i]))
#         }
#     }
#     y_eval = sort(y_eval)
# 
#     preddensity_conditional = list(mode = "list",length = nconfig)
#     preddensity[[rat_idx]]$y_eval = y_eval
#     preddensity[[rat_idx]]$density = numeric(length(y_eval))
#     for(config_idx in 1:nconfig){
#         args$theta = res$misc$configs$config[[config_idx]]$theta[1]
#         preddensity_conditional[[config_idx]] = inla.lik.pred(y_eval,
#                                                   args = args,
#                                                   eta.mean = res$misc$configs$config[[config_idx]]$Predictor[predictor_idx,"mean"],
#                                                   eta.sd = sqrt(res$misc$configs$config[[config_idx]]$Predictor[predictor_idx,"variance"])
#                                                   )
#         preddensity[[rat_idx]]$density  = preddensity[[rat_idx]]$density +  weights[config_idx]*preddensity_conditional[[config_idx]]
#     }
# }
# plot(preddensity[[1]]$y_eval,preddensity[[1]]$density,
#      main = "Predictive Density of A Measure from Rat 1",
#      xlab = "weight",
#      ylab = "density",
#      type="l")
# plot(preddensity[[11]]$y_eval,preddensity[[11]]$density,
#      main = "Predictive Density of A Measure from An Unobserved Rat",
#      xlab = "weight",
#      ylab = "density",
#      type="l")

## ----message = FALSE, warning=FALSE-------------------------------------------
loocv_res = inla.group.cv(result = res)
ULOOCV = mean(log(loocv_res$cv[1:n]))

groups = lapply(1:n,FUN = function(i){which(id == id[i])})
lgocv_res = inla.group.cv(result = res,groups = groups)
ULGOCV = mean(log(lgocv_res$cv[1:n]))

print(paste0("ULOOCV is ", round(ULOOCV,5),"."))
print(paste0("ULGOCV is ", round(ULGOCV,5),"."))

## ----message = FALSE, warning=FALSE-------------------------------------------
lgocv_auto_res = inla.group.cv(result = res,num.level.sets = 1)
ULGOCV_auto = mean(log(lgocv_auto_res$cv[1:n]))
print(paste0("ULGOCV_auto is ", round(ULGOCV_auto,5),"."))
print(paste0("The group for data point 1 is:"))
print(paste(lgocv_auto_res$groups[[1]]$idx,collapse=" "))

## -----------------------------------------------------------------------------
num.level.sets = 3
# We intentionally create some replicated correlations since it is not rare to have 
#|corr(eta_i,eta_j)| == |corr(eta_i,eta_k)|, where j != k.
C_i = c(1,1,0.9,0.9,0.8,0.8,-0.1,-0.1,0,0)
# L_i contains the correlation levels.
L_i = sort(unique(abs(C_i)),decreasing = T)
# Now we construct I_i
I_i = c()
for (l in 1:num.level.sets){
  I_i = c(I_i,which(abs(C_i) == L_i[l]))
}
print(paste(I_i,collapse=" "))

## ----message=FALSE, warning=FALSE,echo = FALSE--------------------------------
n_rat = 10
n_each = 10
n = n_rat*n_each
id = rep(1:n_rat, each = n_each)
x = rnorm(n)
mu_rat = rep(x = rnorm(n_rat,mean = 30,sd = 2),each = n_each) + x
sd_measure = .1
y = rnorm(n = n, mean = mu_rat,sd = sd_measure)
plot(id,y,main = "Data against its group id")
plot(x,y,main = "Data against x")

## ----message = FALSE, warning=FALSE-------------------------------------------
data = data.frame(y = y,id = id,id2 = id,x = x)
formula = y ~ 1 + x + f(id, model = "iid") + f(id2, W = x,model = "iid")
res = inla(formula = formula,
           data = data,
           family = "normal",
           control.compute = list(config = TRUE,likelihood.info = T)
           )

## ----message = FALSE, warning=FALSE-------------------------------------------
formula_alt =  y ~ 1 + x + f(id, model = "iid")
res_alt = inla(formula = formula_alt,
           data = data,
           family = "normal",
           control.compute = list(config = TRUE,likelihood.info = T)
           )

## ----message = FALSE, warning=FALSE-------------------------------------------
lgocv_res = inla.group.cv(result = res,num.level.sets = 1)
lgocv_res_alt = inla.group.cv(result = res_alt,group.cv = lgocv_res)

## ----warning = FALSE----------------------------------------------------------
lgocv_res = inla.group.cv(result = res,num.level.sets = 1)
print(paste(lgocv_res$groups[[1]]$idx,collapse=" "))
print(paste(lgocv_res$groups[[2]]$idx,collapse=" "))
print(paste(lgocv_res$groups[[3]]$idx,collapse=" "))

## ----warning = FALSE----------------------------------------------------------
lgocv_res = inla.group.cv(result = res,num.level.sets = 1,strategy = "prior",keep = "id")
print(paste(lgocv_res$groups[[1]]$idx,collapse=" "))
print(paste(lgocv_res$groups[[2]]$idx,collapse=" "))
print(paste(lgocv_res$groups[[3]]$idx,collapse=" "))

## ----warning = FALSE----------------------------------------------------------
lgocv_res = inla.group.cv(result = res,num.level.sets = 1,select = c(12,53))

## ----warning = FALSE----------------------------------------------------------
groups = vector(mode = "list",length = n)
groups[[12]] = which(id==id[12])
groups[[53]] = which(id==id[53])
lgocv_res = inla.group.cv(result = res,groups = groups)

## ----warning = FALSE----------------------------------------------------------
control.gcpo = list(enable = TRUE,num.level.sets = 1)
res = inla(formula = formula,
           data = data,
           family = "normal",
           control.compute = list(control.gcpo = control.gcpo,config = TRUE)
           )

## ----warning = FALSE----------------------------------------------------------
control.gcpo = list(enable = TRUE,num.level.sets = 1,strategy = "prior",keep = "id")
res = inla(formula = formula,
           data = data,
           family = "normal",
           control.compute = list(control.gcpo = control.gcpo,config = TRUE)
           )
# 1. choose the testing point
testing.point = 1
y_observed = y[testing.point]
# 2. get the hyperparameter density
nconfigs = res$misc$configs$nconfig
weights = unlist(lapply(X = 1:nconfigs,
	FUN = function(config.idx){
    res$misc$configs$config[[config.idx]]$log.posterior 
    + res$misc$configs$config[[config.idx]]$gcpodens.moments[testing.point,3]
    }))
weights = exp(weights)/sum(exp(weights))
#3. get samples of y_i|y_{-I_i},theta
n = 1e6
n_theta = as.numeric(rmultinom(n = 1,size = n,prob = weights))
y_sample = c()
for(config.idx in 1:nconfigs){
    if(n_theta[config.idx]>0){
		eta = rnorm(n_theta[config.idx],
			mean = res$misc$configs$config[[config.idx]]$gcpodens.moments[testing.point,1],
			sd = sqrt(res$misc$configs$config[[config.idx]]$gcpodens.moments[testing.point,2]))
         sd = sqrt(1/exp(res$misc$configs$config[[config.idx]]$theta[1]))
         y_sample_config = rnorm(n = n_theta[config.idx],mean = eta,sd = sd)
         y_sample = c(y_sample,y_sample_config)
    }
}
#4 We can use y_sample to compute the empirical cumulative function then evaluate CRPS
F_y = function(y){sum(y_sample<=y)/n} #inefficient but straightforward...
yy = seq(min(y_sample)-10,max(y_sample)+10,0.01)
FF = unlist(lapply(yy,F_y))
CRPS = pracma::trapz(yy,(FF-(yy>=y_observed))^2)

