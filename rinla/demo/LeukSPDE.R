INLA::inla.setOption(fmesher.evolution.warn = TRUE)
INLA::inla.setOption(fmesher.evolution.verbosity = "stop")

## Redo the Leuk-example using continuous domain frailty
library(INLA)

data(Leuk)
ls()

#############################################################
## model without spatial frailty
model1 <- inla.surv(Leuk$time, Leuk$cens) ~ 0 +
    age + sex + wbc + tpi 

## defining the baseline model, and its hyperparameter prior
pprior <- list(
    theta = list(
        prior = "pc.prec",
        param = c(0.5, 0.01)))
hmodel <- list(
    model = 'rw1',
    n.intervals = 50,
    constr = FALSE, 
    scale.model = TRUE,
    hyper = pprior
)

## data expansion to use it later
## expand the data (internally done in inla)
data.expanded <- inla.coxph(
    formula = model1,
    data = Leuk, 
    control.hazard = hmodel)
str(data.expanded, 1)
data.expanded$formula
names(data.expanded$data)

## fit the model
fit1 <- inla(
    formula = data.expanded$formula, 
    data = c(data.expanded$data,
             data.expanded$data.list), 
    family = "poisson",
    E = E..coxph
)

fit1$summary.hyper[, 1:2]

fit1$summary.fixed[, 1:2]

fit1$mlik

if(FALSE)
    plot(fit1)

rplot <- function(r, ...) {
    n <- nrow(r)
    y <- c(r[, 4], r[n:1, 6], r[1, 4])
    mc <- match.call()
    if(any(names(mc) == 'ylim')) {
        plot(r$ID, r$mean, type = 'l',
             ...)        
    } else { 
        plot(r$ID, r$mean, type = 'l', 
             ylim = range(y),
             ...)
    }
    polygon(r$ID[c(1:n, n:1, 1)], y, 
            col = gray(0.5, 0.5),
            border = gray(0.5, 0.5))
}

par(mfrow = c(1, 1), mar = c(3, 3, 1, 0.1), mgp = c(2, 1, 0), las = 1, bty = 'n')
rplot(fit1$summary.random$baseline, 
      main = 'Baseline with no spatial frailty',
      xlab = 'time', ylab = 'baseline')

### Model with spatial frailty

## get info from the map
library(sf)
if(!any(ls()=="NEmap")) ## older INLA versions does not contains this
    NEmap <- st_sfc(lapply(nwEngland@polygons[[1]]@Polygons, function(p)
        st_polygon(list(p@coords))))

bb <- matrix(st_bbox(NEmap), 2)
bb
rxy <- apply(bb, 1, diff)
rxy
rr <- mean(rxy)

## Create a grid to be used later
g0 <- st_make_grid(
    NEmap,
    cellsize = rr/100
)
g1 <- st_intersection(g0, st_buffer(NEmap, rr/100))

## take only the full pixels
ipoly <- which(sapply(g1, class)[2, ] == "POLYGON")
grid <- st_sf(
    data.frame(ID1 = ipoly),
    geom = st_sfc(g1[ipoly])
)

## non-convex boundaries around the map
library(fmesher)
bnd1 <- fm_nonconvex_hull_inla(
    grid,
    convex = rr * 0.05,
    concave = rr * 0.1)

## build a mesh accounting for boundaries
mesh <- fm_mesh_2d(
    boundary = bnd1, 
    max.edge = rr * c(0.05, 0.15),
    cutoff = rr * 0.02)

mesh$n

plot(mesh)
plot(st_geometry(NEmap), add = TRUE)

## SPDE random effect design matrix
Alocs <- inla.spde.make.A(
    mesh = mesh,
    loc = cbind(data.expanded$data$xcoord,
                data.expanded$data$ycoord)
)
stopifnot(sum(Alocs) == nrow(data.expanded$data))

## SPDE model and hyper-prior definition
rr
spde <- inla.spde2.pcmatern(
    mesh = mesh,
    prior.range = c(rr*0.05, 0.01),
    prior.sigma = c(0.5, 0.01)
)

## model formula
model2 <- update(
    data.expanded$formula, 
    .~. + f(spatial, model = spde, A.local = Alocs))

## fit the model 
fit2 <- inla(
    formula = model2,
    data = c(data.frame(data.expanded$data, spatial = NA),
             data.expanded$data.list), 
    family = "poisson", E = E..coxph
)

fit2$summary.hyperpar[, 1:2]
fit2$summary.fixed[, 1:2]

data.frame(m1 = fit1$mlik,
           m2 = fit2$mlik)

## plot(fit2)

par(mfrow = c(2, 1), mar = c(3, 3, 1, 0.1), mgp = c(2, 1, 0), las = 1, bty = 'n')
rplot(fit1$summary.random$baseline, ylim = c(-12, -7),
      main = 'Baseline with no spatial frailty',
      xlab = 'time', ylab = 'baseline')
rplot(fit2$summary.random$baseline, ylim = c(-12, -7),
      main = 'Baseline under spatial frailty',
      xlab = 'time', ylab = 'baseline')

## evaluate the spatial frailty at the grid locations
gproj <- fm_evaluator(mesh = mesh, loc = st_centroid(grid))

grid$frailty.mean <- fm_evaluate(gproj, field = fit2$summary.random$spatial$mean)
grid$frailty.sd <- fm_evaluate(gproj, field = fit2$summary.random$spatial$sd)

if(FALSE)
    plot(grid[, 3:4], border = 'transparent') 

library(ggplot2)
library(ggpubr)

gg0 <- ggplot(grid) + theme_minimal() 
gg1 <- geom_sf(data = NEmap, fill = 'transparent')

ggarrange(
    gg0 + geom_sf(aes(fill = frailty.mean), color = 'transparent') +
    scale_fill_distiller(
        direction = 1,
        palette = "RdBu"
    ) + gg1,
    gg0 + geom_sf(aes(fill = frailty.sd), color = 'transparent') +
    scale_fill_distiller(
        direction = 1,
        palette = "PuBu"
    ) + gg1 + geom_point(
                  aes(x = xcoord, y = ycoord),
                  data = Leuk, cex = 0.1) + xlab("") + ylab("")
)
