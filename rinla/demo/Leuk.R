data(Leuk)

## graph file
g <- system.file("demodata/Leuk.graph", package="INLA")

## model formula
formula <- inla.surv(Leuk$time, Leuk$cens) ~ sex + age +
    f(inla.group(wbc), model="rw2")+
    f(inla.group(tpi), model="rw2")+
    f(district, model = "besag", graph = g)

## fit the model
result <- inla(formula, family="coxph", data=Leuk)
summary(result)

if(FALSE) {## old code for the map
  
  ## get the function to plot the map
  source(system.file("demodata/Leuk-map.R", package="INLA"))
  Leuk.map(result$summary.random$district$mean)
  
} 

## use the new sf map
library(ggplot2)
ggplot() + theme_minimal() + 
  geom_sf(aes(fill = result$summary.random$district$mean), 
          data = NEmap) +
    labs(fill = "log frailty")

#plot(result)

## Note: check the demo/LeukSPDE.R for continuous spatial domain frailty
