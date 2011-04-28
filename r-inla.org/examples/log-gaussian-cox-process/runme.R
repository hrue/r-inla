library(INLA)

## the 'i' and 'j' in these datafiles are 0-based, so we need to add
## 1L to them later
data = read.table("data-full.dat",  col.names = c("i", "j",  "E", "y"))
altitude = read.table("altitude-full.dat",  col.names = c("i", "j",  "altitude"))
gradient = read.table("gradient-full.dat",  col.names = c("i", "j",  "gradient"))

## all indices in these files have the same ordering, which makes it
## easier.
nrow = 101
ncol = 201

data$idx.data = inla.lattice2node(data$i + 1L, data$j + 1L,  nrow,  ncol)
data$altitude = altitude$altitude
data$gradient = gradient$gradient

## add a better initial value...
formula = y ~ 1 + f(idx.data,  model="rw2d", nrow=nrow, ncol=ncol,
        hyper = list(prec = list(initial = 0))) + altitude + gradient

r = inla(formula,  family = "poisson",  E=E, data = data, verbose=TRUE,
        control.inla = list(strategy = "gaussian"))

inla.display.matrix( matrix(r$summary.random$idx.data$mean, nrow, ncol) )


