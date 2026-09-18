## ----setup, include=FALSE-------------------------------------------
set.seed(123)
library(INLA)
inla.setOption(num.threads="1:1")
inla.setOption(smtp="taucs")
if (file.exists("myinit.R")) source("myinit.R")
library(knitr)
library(rmarkdown)
knitr::opts_chunk$set(echo=TRUE, cache=FALSE, message=FALSE,warning=FALSE)
knitr::opts_chunk$set(fig.path="figures/conditional-logit/")

## ----formula, include=TRUE,eval=F-----------------------------------
# fisher.dat <- readRDS(system.file("demodata/data_fisher2.rds", package
# = "INLA"))
# fisher.dat$id1 <- fisher.dat$id
# fisher.dat$dist_cent <- scale(fisher.dat$dist_cent)
# 
# formula.inla <- y ~ sex + landuse + dist_cent +
#    f(stratum,model="iid",hyper=list(theta = list(initial=log(1e-6),fixed=T))) +
#    f(id1,dist_cent, model="iid")

## ----inla_call, include=TRUE,eval=F---------------------------------
# r.inla <- inla(formula.inla, family ="Poisson", data=fisher.dat)

