nrow=101
ncol=201

## all these files has the same ordering, index = 0,1,2,....
gradient = read.table("gradient-full.dat")$V2
altitude = read.table("altitude-full.dat")$V2
E = read.table("data-full.dat")$V2
y = read.table("data-full.dat")$V3

ind = 1:(nrow*ncol)
formula = y ~ 1 + gradient + altitude + f(ind, model="rw2d", ncol=ncol, nrow=nrow)

d = list(y=y, gradient = gradient, altitude = altitude)
res = inla(formula, family = "poisson", data=d, E = E, verbose = TRUE,
           control.inla = list(strategy = "gaussian"),
           control.predictor = list(fixed = FALSE))

