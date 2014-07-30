library (INLA)

## Crossed Random Effects - Salamander

load("salam.RData")
## organize data into a form suitable for logistic regression
dat0=data.frame("y"=c(salam$y), "fW"=as.integer(salam$x[,"W/R"]==1 | salam$x[,"W/W"]==1), 
    "mW"=as.integer(salam$x[,"R/W"]==1 | salam$x[,"W/W"]==1), 
    "WW"=as.integer(salam$x[,"W/W"]==1 ) )
## add salamander id
id = t( apply(salam$z, 1, function(x) {
        tmp = which (x==1)
        tmp[2] = tmp[2] - 20
        tmp
    }) ) 
## ids are suitable for model A and C, but not B
id.modA = rbind(id, id+40, id+20)
colnames (id.modA) = c("f.modA","m.modA")
dat0=cbind (dat0, id.modA, group=1)
dat0$experiment=as.factor(rep(1:3, each=120))
dat0$group=as.factor(dat0$group)

salamander = dat0
salamander.e1 = subset (dat0, dat0$experiment==1)
salamander.e2 = subset (dat0, dat0$experiment==2)
salamander.e3 = subset (dat0, dat0$experiment==3)

# salamander.e1
salamander.e1.inla.fit = inla(y~fW+mW+WW + f(f.modA, model="iid", param=c(1,.622)) + f(m.modA, model="iid", param=c(1,.622)), 
        family="binomial", data=salamander.e1, Ntrials=rep(1,nrow(salamander.e1)))
salamander.e1.hyperpar = inla.hyperpar (salamander.e1.inla.fit)
summary(salamander.e1.inla.fit)
summary(salamander.e1.hyperpar)

inla.emarginal(function(x) 1/x^.5, salamander.e1.hyperpar$marginals[[1]])
inla.emarginal(function(x) 1/x^.5, salamander.e1.hyperpar$marginals[[2]])

## same for salamander.e2 and salamander.e3
