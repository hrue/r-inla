## Redo the Leuk-example using more senible PC-priors
data(Leuk)

## graph file
g <- system.file("demodata/Leuk.graph", package="INLA")

## pc-prior parameters
U <- 0.3/0.31
alpha <- 0.01
pc.param <- list(prec = list(prior = "pc.prec", param = c(U, alpha)))

## model formula
formula <- inla.surv(Leuk$time, Leuk$cens) ~ sex + age +
    f(inla.group(wbc), model="rw2",
      scale.model = TRUE, hyper = pc.param)+
    f(inla.group(tpi), model="rw2",
      scale.model = TRUE, hyper = pc.param)+
    f(district,model="besag",graph = g,
      scale.model = TRUE, hyper = pc.param)

## fit the model 
result <- inla(formula, family="coxph", data=Leuk,
        control.hazard = list(scale.model = TRUE, hyper = pc.param))
summary(result)

if(FALSE) {## old code for the map
  
  ## get the function to plot the map
  source(system.file("demodata/Leuk-map.R", package="INLA"))
  Leuk.map(result$summary.random$district$mean)
  
}

## use the new sf map
library(ggplot2)

ggplot() + theme_minimal() + 
  geom_sf(aes(fill = exp(result$summary.random$district$mean)), 
          data = NEmap) +
    labs(fill = "frailty")

#plot(result)

## Note: check the demo/LeukSPDE.R for continuous spatial domain frailty
