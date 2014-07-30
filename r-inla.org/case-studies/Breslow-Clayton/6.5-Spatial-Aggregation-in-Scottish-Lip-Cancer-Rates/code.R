library (INLA)

# Spatial Aggregation - Scotland

data(Scotland)
Scotland$Region2=Scotland$Region

scotland.inla.fit.1 = inla(Counts ~ 1 + f(Region, model="iid", param=c(1, .0014)),
        data=Scotland, family="poisson", E=E )
scotland.hyperpar.1 = inla.hyperpar (scotland.inla.fit.1)
summary(scotland.inla.fit.1)
inla.emarginal(function(x) 1/x^.5, scotland.hyperpar.1$marginals[[1]])

scotland.inla.fit.2 = inla(Counts ~ 1 + I(X/10) +
        f(Region, model="iid", param=c(1, .0014)),
        data=Scotland, family="poisson", E=E )
scotland.hyperpar.2 = inla.hyperpar (scotland.inla.fit.2)
summary(scotland.inla.fit.2)
inla.emarginal(function(x) 1/x^.5, scotland.hyperpar.2$marginals[[1]])

scotland.inla.fit.3 = inla(Counts ~ 1 +
        f(Region, model="iid", param=c(1, .0014)) +
        f(Region2, model="besag", graph="scotland.graph", param=c(1, .2/.59)),
        data=Scotland, family="poisson", E=E )
scotland.hyperpar.3 = inla.hyperpar (scotland.inla.fit.3)
summary(scotland.inla.fit.3)
inla.emarginal(function(x) 1/x^.5, scotland.hyperpar.3$marginals[[1]])

scotland.inla.fit.4 = inla(Counts ~ 1 + I(X/10) +
        f(Region, model="iid", param=c(1, .0014)) +
        f(Region2, model="besag", graph="scotland.graph", param=c(1, .4/.59)),
        data=Scotland, family="poisson", E=E )
scotland.hyperpar.4 = inla.hyperpar (scotland.inla.fit.4)
summary(scotland.inla.fit.4)
inla.emarginal(function(x) 1/x^.5, scotland.hyperpar.4$marginals[[1]])
