library (INLA)

# Log Odds Ratio - Oxford

load("oxford.Rdata")
oxford$birth.year = oxford$birth.year - 54 ## center

formula = cases ~ as.factor(birth.year) + exposed +
    I(exposed * birth.year) + I(exposed * (birth.year**2-22)) +
    f(birth.year, exposed, model="iid", param=c(.5, .0164))

oxford.inla.fit = inla(formula, data=oxford, family="binomial", Ntrials=total)

oxford.hyperpar = inla.hyperpar (oxford.inla.fit)
summary(oxford.inla.fit)
inla.emarginal(function(x) 1/x^.5, oxford.hyperpar$marginals[[1]])
