data(BivMetaAnalysis)

## formulation 1
formula <- Y~f(diid,model="2diidwishart",param=c(4,1,2,0.1))+lag.tp+ lag.tn+ ct.tp+ ct.tn+ mr.tp+ mr.tn -1
model=inla(formula,family="binomial", data=BivMetaAnalysis, Ntrials=N,keep=TRUE)
h = inla.hyperpar(model)


## we know that diid = 1:n. assigne the odd numbers to part0 and the even ones to part1
n = dim(BivMetaAnalysis)[1]
k = rep(NA,n)
k[ seq(1,n,by = 2) ] = 1:(n/2)
BivMetaAnalysis2 = cbind(BivMetaAnalysis, "diid.part0" = k)

k = rep(NA,n)
k[ seq(2,n,by = 2) ] = 1:(n/2)
BivMetaAnalysis2 = cbind(BivMetaAnalysis2, "diid.part1" = k)

## now we can do formulation 2 of the problem.
formula2 <- Y~f(diid.part0,model="2diidwishartpart0",param=c(4,1,2,0.1))+f(diid.part1,model="2diidwishartpart1")+lag.tp+ lag.tn+ ct.tp+ ct.tn+ mr.tp+ mr.tn -1
model2=inla(formula2,family="binomial", data=BivMetaAnalysis2, Ntrials=N)

h2 = inla.hyperpar(model2)
