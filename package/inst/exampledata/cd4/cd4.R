library(INLA)
library(JMbayes)

data("aids")
data1 <- aids
data1$Time=data1$Time/max(max(data1$Time),max(data1$obstime))
data1$obstime=data1$obstime/max(max(data1$Time),max(data1$obstime))
datas <- data1[data1$obstime==0,]
datal <- data1[,c(1,4:12)]
ns <- nrow(datas)
nl <- nrow(datal)

y.long <- c(datal$CD4,rep(NA, ns))
y.surv <- inla.surv(time = c(rep(NA, nl),datas$Time ), event = c(rep(NA, nl),datas$death))
Yjoint <- list(y.long, y.surv)
N <- length(unique(data1$patient))

fixed.covariate <- data.frame(mu = as.factor(c(rep(1, nl), rep(2, ns))),
                              l.drug = as.factor(c(as.factor(datal$drug), rep(NA, ns))),
                              l.gender = as.factor(c(as.factor(datal$gender), rep(NA, ns))),
                              l.prevOI = as.factor(c(as.factor(datal$prevOI) , rep(NA, ns))),
                              l.AZT = as.factor(c(as.factor(datal$AZT) , rep(NA, ns))),
                              s.drug = as.factor(c( rep(NA, nl),as.factor(datas$drug))),
                              s.gender = as.factor(c(rep(NA, nl),as.factor(datas$gender))),
                              s.prevOI = as.factor(c(rep(NA, nl),as.factor(datas$prevOI))),
                              s.AZT = as.factor(c(rep(NA, nl),as.factor(datas$AZT))),
                              l.time=c(datal$obstime, rep(0, ns)),
                              s.time=c(rep(0, nl),datas$Time ))
random.covariate <- list(U11 = c(datal$patient,rep(NA, ns)),
                         U21 = c(datal$patient,rep(NA, ns)),
                         U12=c(rep(NA,nl),datas$patient),
                         U22=c(rep(NA,nl),datas$patient))
joint.dataCD4  <-  c(fixed.covariate,random.covariate)
joint.dataCD4$Y  <-  Yjoint #What I need to "call" is joint.dataCD4

saveRDS(joint.dataCD4, file = "jointdataCD4.rds")
