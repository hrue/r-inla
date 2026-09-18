## ----setup, include=FALSE, cache=FALSE------------------------------
knitr::opts_chunk$set(echo = TRUE, warning = FALSE, error = FALSE, message = FALSE)

## ----libs, warning=FALSE, message=FALSE, errors=FALSE---------------
# libraries
library(maps)
library(ggplot2)
library(sf)
library(terra)
library(tidyterra) # raster plotting
library(tidyr)
library(scales)
library(inlabru)
library(INLA)
library(fmesher)
library(dplyr)
# Note: the 'splancs' package also needs to be installed,
# but doesn't need to be loaded

# set option
select <- dplyr::select
options(scipen = 99999)
options(max.print = 99999)
options(stringsAsFactors = FALSE)

## ----echo=F---------------------------------------------------------
inla.setOption(num.threads = "1:1")
inla.setOption(smtp = "taucs")
inla.setOption(inla.mode = "experimental")
if (file.exists("myinit.R")) source("myinit.R")

## ----proj, warning=F, message=FALSE---------------------------------
# define a crs
epsg6703km <- paste(
  "+proj=aea +lat_0=23 +lon_0=-96 +lat_1=29.5",
  "+lat_2=45.5 +x_0=0 +y_0=0 +datum=NAD83",
  "+units=km +no_defs"
)

# make a base map
states <- maps::map("state", plot = FALSE, fill = TRUE) %>%
  sf::st_as_sf() %>%
  filter(ID %in% c(
    "texas", "oklahoma", "kansas", "missouri",
    "arkansas", "louisiana"
  )) %>%
  sf::st_make_valid() %>%
  sf::st_transform(epsg6703km) %>%
  sf::st_make_valid()

## ----data-init, eval=T, cache=FALSE---------------------------------
# get data 
## download_dat <- read.csv(paste0(
##      "https://raw.github.com/tmeeha/inlaSVCBC",
##      "/master/code/modeling_data.csv"))
download_dat <- read.csv("data/modeling_data.csv.gz")

# select subset
count_dat <- download_dat %>%
  select(
    circle, bcr, state, year, std_yr, count, log_hrs,
    lon, lat, obs
  ) %>%
  mutate(year = year + 1899) %>%
  filter(
    state %in% c(
      "TEXAS", "OKLAHOMA", "KANSAS", "MISSOURI",
      "ARKANSAS", "LOUISIANA"
    ),
    year >= 1987
  )

## ----data-subset, eval=TRUE-----------------------------------------
count_dat <- count_dat %>%
  mutate(site_idx = as.numeric(factor(paste(circle, lon, lat)))) %>%
  group_by(site_idx) %>%
  mutate(n_years = n()) %>%
  filter(n_years >= 20) %>%
  ungroup() %>%
  mutate(
    std_yr = year - max(year),
    obs = seq_len(nrow(.)),
    site_idx = as.numeric(factor(paste(circle, lon, lat))),
    year_idx = as.numeric(factor(year)),
    site_year_idx = as.numeric(factor(paste(circle, lon, lat, year)))
  ) %>%
  st_as_sf(coords = c("lon", "lat"), crs = 4326, remove = FALSE) %>%
  st_transform(epsg6703km) %>%
  mutate(
    easting = st_coordinates(.)[, 1],
    northing = st_coordinates(.)[, 2]
  ) %>%
  arrange(circle, year)

## ----data_map, fig.pos='b', fig.align='center'----------------------
# map it
ggplot() +
  geom_sf(
    data = count_dat %>% filter(year_idx %in% seq(1, 30, 3)),
    aes(col = log(count + 1))
  ) +
  geom_sf(data = states, fill = NA) +
  coord_sf(datum = NA) +
  facet_wrap(~year) +
  scale_color_distiller(palette = "Spectral") +
  theme_bw()

## ----spat_dat-------------------------------------------------------
# make a set of distinct study sites for mapping
site_map <- count_dat %>%
  select(circle, easting, northing) %>%
  distinct() %>%
  select(circle, easting, northing)

# save coordinates as matrices for later
unique_coords <- count_dat %>%
  select(circle, easting, northing) %>%
  distinct() %>%
  select(easting, northing) %>%
  st_drop_geometry() %>%
  as.matrix()
all_coords <- as.matrix(st_coordinates(count_dat))

## ----mesh, cache=FALSE----------------------------------------------
# make a hull and mesh for spatial model
hull <- fm_nonconvex_hull(
  count_dat,
  convex = 200,
  concave = 350
)
mesh <- fm_mesh_2d(
  boundary = hull, max.edge = c(100, 600), # km inside and outside
  cutoff = 50, offset = c(100, 300),
  crs = fm_crs(count_dat)
) # cutoff is min edge

## ----mesh-plot------------------------------------------------------
# plot it
ggplot() +
  gg(data = mesh) +
  geom_sf(data = site_map, col = "darkgreen", size = 1) +
  geom_sf(data = states, fill = NA) +
  theme_bw() +
  labs(x = "", y = "")

## ----spde-----------------------------------------------------------
# make spde
spde <- inla.spde2.pcmatern(
  mesh = mesh,
  prior.range = c(500, 0.5),
  prior.sigma = c(1, 0.5)
)

## ----indices--------------------------------------------------------
# make index sets
alpha_idx <- inla.spde.make.index(name = "alpha", n.spde = mesh$n)
eps_idx <- inla.spde.make.index(name = "eps", n.spde = mesh$n)
tau_idx <- inla.spde.make.index(name = "tau", n.spde = mesh$n)

## ----projector------------------------------------------------------
# make projector matrices
A_alpha <- inla.spde.make.A(mesh = mesh, loc = count_dat)
A_eps <- inla.spde.make.A(
  mesh = mesh, loc = count_dat,
  weights = count_dat$log_hrs
) # note weights argument
A_tau <- inla.spde.make.A(
  mesh = mesh, loc = count_dat,
  weights = count_dat$std_yr
) # note weights argument

## ----stack1---------------------------------------------------------
# stack observed data
stack_fit <- inla.stack(
  tag = "obs",
  data = list(count = as.vector(count_dat$count)), # response from data frame
  effects = list(data.frame(
    intercept = 1,
    kappa = count_dat$site_idx
  ), # predictors from data frame
  alpha = alpha_idx, # or index sets if spatial
  eps = eps_idx,
  tau = tau_idx
  ),
  A = list(
    1, # a value of 1 is given for non-spatial terms
    A_alpha,
    A_eps,
    A_tau
  )
)

## ----prior----------------------------------------------------------
# iid prior
pc_prec <- list(prior = "pcprec", param = c(1, 0.1))

## ----form-----------------------------------------------------------
# formula
svc_form <- count ~ -1 +
  f(kappa, model = "iid", constr = TRUE, hyper = list(prec = pc_prec)) +
  f(alpha, model = spde) +
  f(eps, model = spde) +
  f(tau, model = spde)

## ----run, eval=T, cache=FALSE---------------------------------------
res <- inla(svc_form,
  family = "nbinomial",
  data = inla.stack.data(stack_fit), # stack the stack
  control.predictor = list(
    A = inla.stack.A(stack_fit),
    compute = FALSE
  ), # must define A
  control.compute = list(waic = TRUE, cpo = TRUE, config = TRUE),
  control.inla = list(strategy = "adaptive", int.strategy = "eb"),
  verbose = FALSE
)

## ----load_res, echo=F, eval=F---------------------------------------
# load("spatial_cbc/code/spde_vignette_v3.RData")

## ----vc-------------------------------------------------------------
# view results
res$summary.hyperpar[-1, c(1, 2)]

## ----alpha_res------------------------------------------------------
summary(exp(res$summary.random$alp$"0.5quant")) # exp(alpha) posterior median

## ----eps_res--------------------------------------------------------
summary(res$summary.random$eps$mean) # epsilon

## ----tau_res--------------------------------------------------------
summary((exp(res$summary.random$tau$"0.5quant") - 1) * 100) # (exp(tau)-1)*100

## ----map_grid-------------------------------------------------------
# get easting and northing limits
bbox <- fm_bbox(hull[[1]])
grd_dims <- round(c(x = diff(bbox[[1]]), y = diff(bbox[[2]])) / 25)

# make mesh projector to get model summaries from the mesh to the mapping grid
mesh_proj <- fm_evaluator(
  mesh,
  xlim = bbox[[1]], ylim = bbox[[2]], dims = grd_dims
)

## ----map_vals-------------------------------------------------------
# pull data
kappa <- data.frame(
  median = exp(res$summary.random$kappa$"0.5quant"),
  range95 = exp(res$summary.random$kappa$"0.975quant") -
    exp(res$summary.random$kappa$"0.025quant")
)
alph <- data.frame(
  median = exp(res$summary.random$alpha$"0.5quant"),
  range95 = exp(res$summary.random$alpha$"0.975quant") -
    exp(res$summary.random$alpha$"0.025quant")
)
epsi <- data.frame(
  median = res$summary.random$eps$"0.5quant",
  range95 = (res$summary.random$eps$"0.975quant" -
    res$summary.random$eps$"0.025quant")
)
taus <- data.frame(
  median = (exp(res$summary.random$tau$"0.5quant") - 1) * 100,
  range95 = (exp(res$summary.random$tau$"0.975quant") -
    exp(res$summary.random$tau$"0.025quant")) * 100
)

# loop to get estimates on a mapping grid
pred_grids <- lapply(
  list(alpha = alph, epsilon = epsi, tau = taus),
  function(x) as.matrix(fm_evaluate(mesh_proj, x))
)

## ----map_rast-------------------------------------------------------
# make a terra raster stack with the posterior median and range95
out_stk <- rast()
for (j in 1:3) {
  mean_j <- cbind(expand.grid(x = mesh_proj$x, y = mesh_proj$y),
    Z = c(matrix(pred_grids[[j]][, 1], grd_dims[1]))
  )
  mean_j <- rast(mean_j, crs = epsg6703km)
  range95_j <- cbind(expand.grid(X = mesh_proj$x, Y = mesh_proj$y),
    Z = c(matrix(pred_grids[[j]][, 2], grd_dims[1]))
  )
  range95_j <- rast(range95_j, crs = epsg6703km)
  out_j <- c(mean_j, range95_j)
  terra::add(out_stk) <- out_j
}
names(out_stk) <- c(
  "alpha_median", "alpha_range95", "epsilon_median",
  "epsilon_range95", "tau_median", "tau_range95"
)
out_stk <- terra::mask(
  out_stk,
  terra::vect(sf::st_union(states)),
  updatevalue = NA,
  touches = FALSE
)

## ----svc_plot_construction, warning=F, message=F--------------------
make_plot_field <- function(data_stk, scale_label) {
  ggplot(states) +
    geom_sf(fill = NA) +
    coord_sf(datum = NA) +
    geom_spatraster(data = data_stk) +
    labs(x = "", y = "") +
    scale_fill_distiller(
      scale_label,
      palette = "Spectral",
      na.value = "transparent"
    ) +
    theme_bw() +
    geom_sf(fill = NA)
}
make_plot_site <- function(data, scale_label) {
  ggplot(states) +
    geom_sf() +
    coord_sf(datum = NA) +
    geom_sf(data = data, size = 1, mapping = aes(colour = value)) +
    scale_colour_distiller(scale_label, palette = "Spectral") +
    labs(x = "", y = "") +
    theme_bw() +
    geom_sf(fill = NA)
}

# medians
# fields alpha_s, epsilon_s, tau_s
pa <- make_plot_field(
  data_stk = out_stk[["alpha_median"]],
  scale_label = "posterior\nmedian\nexp(alpha_s)"
)
pe <- make_plot_field(
  data_stk = out_stk[["epsilon_median"]],
  scale_label = "posterior\nmedian\nepsilon_s"
)
pt <- make_plot_field(
  data_stk = out_stk[["tau_median"]],
  scale_label = "posterior\nmedian\n100(exp(tau_s)-1)"
)
# sites kappa_s
ps <- make_plot_site(
  data = cbind(site_map, data.frame(value = kappa$median)),
  scale_label = "posterior\nmedian\nexp(kappa_s)"
)
# range95
# fields alpha_s, epsilon_s, tau_s
pa_range95 <- make_plot_field(
  data_stk = out_stk[["alpha_range95"]],
  scale_label = "posterior\nrange95\nexp(alpha_s)"
)
pe_range95 <- make_plot_field(
  data_stk = out_stk[["epsilon_range95"]],
  scale_label = "posterior\nrange95\nepsilon_s"
)
pt_range95 <- make_plot_field(
  data_stk = out_stk[["tau_range95"]],
  scale_label = "posterior\nrange95\n100(exp(tau_s)-1)"
)
# sites kappa_s
ps_range95 <- make_plot_site(
  data = cbind(site_map, data.frame(value = kappa$range95)),
  scale_label = "posterior\nrange95\nexp(kappa_s)"
)

## ----svc_plots, warning=F, message=F, fig.width=8, fig.align='center'----
# plot together
multiplot(ps, pa, pe, pt, cols = 2)

## ----svc_plots_range95, warning=F, message=F, fig.width=8, fig.align='center'----
# plot together
multiplot(ps_range95, pa_range95, pe_range95, pt_range95, cols = 2)

