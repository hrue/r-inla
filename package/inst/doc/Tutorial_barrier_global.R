## ----setup, include=FALSE-------------------------------------------
set.seed(123)
library(INLA)
##inla.setOption(num.threads="4:1")
##inla.setOption(smtp="taucs")
##if (file.exists("myinit.R")) source("myinit.R")
library(knitr)
library(rmarkdown)
knitr::opts_chunk$set(echo=TRUE, cache=FALSE, message=FALSE, warning=FALSE)
knitr::opts_chunk$set(fig.path="figures/barrier-global/")

