library(INLA)
library(JointModel)
data1 <- prostate
data2 <- dropout
ng <- nrow(data1)
ns <- nrow(data2)
data1 = data1[order(data1$VisitTime),] #Just for "nice" graphs

## prepare the response variable
library(INLA)
y.longPSA <- c(data1$logPSA.postRT , rep(NA,ng),rep(NA, ns))
y.etaPSA <- c(rep(NA,ng),rep(0,ng),rep(NA,ns))
y.survPSA <- inla.surv(time = c(rep(NA, ng),rep(NA,ng), data2$DropTime),
                       event = c(rep(NA, ng),rep(NA,ng),data2$Status))
Yjoint <- list(y.longPSA, y.etaPSA,y.survPSA)
N <- length(unique(prostate$ID))
linear.covariatePSA <- data.frame(mu = as.factor(c(rep(1, ng),rep(1, ng),rep(2, ns))),
                                  b13.PSAbase = c(data1$logPSA.base,data1$logPSA.base, rep(0, ns)))
random.covariatePSA <- list(V1 = c(data1$VisitTime,rep(NA,ng),rep(NA,ns)),
                            u = c(rep(NA,ng),data1$ID,rep(NA,ns)),
                            w = c(rep(NA,ng),rep(-1,ng),rep(NA,ns)),
                            b.eta = c(rep(NA,ng),rep(NA,ng),data2$ID2))
## Entire linear predictor
joint.dataPSA <- c(linear.covariatePSA, random.covariatePSA)
joint.dataPSA$Y <- Yjoint
saveRDS(joint.dataPSA,  file = "jointdataPSA.rds")
