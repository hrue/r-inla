library(INLA)
data(Seeds)

seeds.inla.fit.1 = inla(r ~ x1 + x2 + f(plate, model="iid",
        param=c(.5, .0164)), data=Seeds, family="binomial", Ntrials=n )

seeds.hyperpar.1 = inla.hyperpar(seeds.inla.fit.1)
summary(seeds.inla.fit.1)
inla.emarginal(function(x) 1/x^.5, seeds.hyperpar.1$marginals[[1]])

seeds.inla.fit.2 = inla(r ~ x1 + x2 + I(x1*x2) + f(plate, model="iid",
        param=c(.5, .0164)), data=Seeds, family="binomial", Ntrials=n )

seeds.hyperpar.2 = inla.hyperpar(seeds.inla.fit.2)
summary(seeds.inla.fit.2)
inla.emarginal(function(x) 1/x^.5, seeds.hyperpar.2$marginals[[1]])
