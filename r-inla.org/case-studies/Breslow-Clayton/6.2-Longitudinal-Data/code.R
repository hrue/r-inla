library (INLA)

## visit should be 1 to 4

## Visit is created by glmmAK, it corresponds to Breslow and Clayton's
## Visit/10, because the codes are -3,-1,1,3.
require (glmmAK)
data(epilepticBC)
epil = epilepticBC
epil$id2=epil$id
epil$rand=1:nrow(epil)
epil$V4=epil$visit==4
epil$newid=rep(1:(nrow(epil)/4), each=4)

epil.inla.fit.1 = inla(Seizure ~ Base + Trt + I(Base*Trt) + Age + V4 +
        f(id,model="iid",param=c(2, 1.140), diagonal=0), data=epil,
        family="poisson" )

epil.hyperpar.1 = inla.hyperpar(epil.inla.fit.1)
summary(epil.inla.fit.1)
inla.emarginal(function(x) 1/x^.5, epil.hyperpar.1$marginals[[1]])

epil.inla.fit.2 = inla(Seizure ~ Base + Trt + I(Base*Trt) + Age + V4 +
        f(id,model="iid",param=c(2, 1.240), diagonal=0) +
        f(rand,model="iid",param=c(2, 1.140), diagonal=0), data=epil,
        family="poisson" )

epil.hyperpar.2 = inla.hyperpar(epil.inla.fit.2)
summary(epil.inla.fit.2)
inla.emarginal(function(x) 1/x^.5, epil.hyperpar.2$marginals[[1]])

epil.inla.fit.3 = inla(Seizure ~ Base + Trt + I(Base*Trt) + Age +
        Visit +f(id, model="2diidwishartpart0", param=c(5, 2.277904,
        1.692047, 0), diagonal=0) + f(id2, Visit,
        model="2diidwishartpart1", diagonal = 0), data=epil,
        family="poisson" )

epil.hyperpar.3 = inla.hyperpar(epil.inla.fit.3)
summary(epil.inla.fit.3)
inla.emarginal(function(x) 1/x^.5, epil.hyperpar.3$marginals[[1]])


