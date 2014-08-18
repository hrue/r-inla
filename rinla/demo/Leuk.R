data(Leuk)
## get the function to plot the map
source(system.file("demodata/Leuk-map.R", package="INLA"))
g = system.file("demodata/Leuk.graph", package="INLA")

formula = inla.surv(Leuk$time, Leuk$cens) ~ sex + age +
    f(inla.group(wbc), model="rw2")+
    f(inla.group(tpi), model="rw2")+
    f(district,model="besag",graph = g)

result = inla(formula, family="coxph", data=Leuk)
summary(result)
Leuk.map(result$summary.random$district$mean)
plot(result)
