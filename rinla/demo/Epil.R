data(Epil)
my.center = function(x) (x - mean(x))

## make centered covariates
Epil$CTrt    = my.center(Epil$Trt)
Epil$ClBase4 = my.center(log(Epil$Base/4))
Epil$CV4     = my.center(Epil$V4)
Epil$ClAge   = my.center(log(Epil$Age))
Epil$CBT     = my.center(Epil$Trt*Epil$ClBase4)

formula = y ~ ClBase4 + CTrt + CBT+ ClAge + CV4 +
          f(Ind, model="iid") + f(rand,model="iid")

result = inla(formula,family="poisson", data = Epil)
summary(result)
plot(result)

hyper = inla.hyperpar(result)
summary(hyper)
plot(hyper)

