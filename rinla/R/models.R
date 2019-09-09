## Export: inla.models

## this is how to define a must-be-enabled...
## status = "changed:Oct.25.2017", 


`inla.models.section.latent` = function()
{
    return
    list(latent =
             list(
                 linear = list(
                     doc = "Alternative interface to an fixed effect", 
                     hyper = list(
                     ),
                     constr = FALSE,
                     nrow.ncol = FALSE,
                     augmented = FALSE,
                     aug.factor = 1L,
                     aug.constr = NULL,
                     n.div.by = NULL,
                     n.required = FALSE,
                     set.default.values = FALSE,
                     pdf = "linear"
                 ),

                 iid = list(
                     doc = "Gaussian random effects in dim=1", 
                     hyper = list(
                         theta = list(
                             hyperid =  1001, 
                             name = "log precision",
                             short.name = "prec",
                             prior = "loggamma",
                             param = c(1, 0.00005),
                             initial = 4,
                             fixed = FALSE,
                             to.theta = function(x) log(x),
                             from.theta = function(x) exp(x)
                         )
                     ),
                     constr = FALSE,
                     nrow.ncol = FALSE,
                     augmented = FALSE,
                     aug.factor = 1L,
                     aug.constr = NULL,
                     n.div.by = NULL,
                     n.required = FALSE,
                     set.default.values = FALSE,
                     pdf = "indep"
                 ),

                 mec = list(
                     doc = "Classical measurement error model", 
                     hyper = list(
                         theta1 = list(
                             hyperid =  2001,
                             name = "beta",
                             short.name = "b",
                             prior = "gaussian",
                             param = c(1, 0.001),
                             initial = 1,
                             fixed = FALSE,
                             to.theta = function(x) x,
                             from.theta = function(x) x
                         ), 
                         theta2 = list(
                             hyperid =  2002,
                             name = "prec.u",
                             short.name = "prec",
                             prior = "loggamma",
                             param = c(1, 0.0001),
                             initial = log(1/0.0001),
                             fixed = TRUE,
                             to.theta = function(x) log(x),
                             from.theta = function(x) exp(x)
                         ), 
                         theta3 = list(
                             hyperid =  2003,
                             name = "mean.x",
                             short.name = "mu.x",
                             prior = "gaussian",
                             param = c(0, 0.0001),
                             initial = 0,
                             fixed = TRUE,
                             to.theta = function(x) x,
                             from.theta = function(x) x
                         ), 
                         theta4 = list(
                             hyperid =  2004,
                             name = "prec.x",
                             short.name = "prec.x",
                             prior = "loggamma",
                             param = c(1, 10000),
                             initial = log(1/10000),
                             fixed = TRUE,
                             to.theta = function(x) log(x),
                             from.theta = function(x) exp(x)
                         )
                     ),
                     constr = FALSE,
                     nrow.ncol = FALSE,
                     augmented = FALSE,
                     aug.factor = 1L,
                     aug.constr = NULL,
                     n.div.by = NULL,
                     n.required = FALSE,
                     set.default.values = FALSE,
                     pdf = "mec"
                 ),

                 meb = list(
                     doc = "Berkson measurement error model", 
                     hyper = list(
                         theta1 = list(
                             hyperid =  3001,
                             name = "beta",
                             short.name = "b",
                             prior = "gaussian",
                             param = c(1, 0.001),
                             initial = 1,
                             fixed = FALSE,
                             to.theta = function(x) x,
                             from.theta = function(x) x
                         ), 
                         theta2 = list(
                             hyperid =  3002,
                             name = "prec.u",
                             short.name = "prec",
                             prior = "loggamma",
                             param = c(1, 0.0001),
                             initial = log(1000),
                             fixed = FALSE,
                             to.theta = function(x) log(x),
                             from.theta = function(x) exp(x)
                         )
                     ),
                     constr = FALSE,
                     nrow.ncol = FALSE,
                     augmented = FALSE,
                     aug.factor = 1L,
                     aug.constr = NULL,
                     n.div.by = NULL,
                     n.required = FALSE,
                     set.default.values = FALSE,
                     pdf = "meb"
                 ),

                 rgeneric = list(
                     doc = "Generic latent model spesified using R", 
                     hyper = list(), 
                     constr = FALSE,
                     nrow.ncol = FALSE,
                     augmented = FALSE,
                     aug.factor = 1L,
                     aug.constr = NULL,
                     n.div.by = NULL,
                     n.required = TRUE,
                     set.default.values = TRUE,
                     status = "experimental", 
                     pdf = "rgeneric"
                 ),

                 rw1 = list(
                     doc = "Random walk of order 1", 
                     hyper = list(
                         theta = list(
                             hyperid =  4001,
                             name = "log precision",
                             short.name = "prec",
                             prior = "loggamma",
                             param = c(1, 0.00005),
                             initial = 4,
                             fixed = FALSE,
                             to.theta = function(x) log(x),
                             from.theta = function(x) exp(x)
                         )
                     ),
                     constr = TRUE,
                     nrow.ncol = FALSE,
                     augmented = FALSE,
                     aug.factor = 1L,
                     aug.constr = NULL,
                     n.div.by = NULL,
                     n.required = FALSE,
                     set.default.values = FALSE,
                     min.diff = 1E-5, 
                     pdf = "rw1"
                 ),

                 rw2 = list(
                     doc = "Random walk of order 2", 
                     hyper = list(
                         theta = list(
                             hyperid =  5001,
                             name = "log precision",
                             short.name = "prec",
                             prior = "loggamma",
                             param = c(1, 0.00005),
                             initial = 4,
                             fixed = FALSE,
                             to.theta = function(x) log(x),
                             from.theta = function(x) exp(x)
                         )
                     ),
                     constr = TRUE,
                     nrow.ncol = FALSE,
                     augmented = FALSE,
                     aug.factor = 1L,
                     aug.constr = NULL,
                     n.div.by = NULL,
                     n.required = FALSE,
                     set.default.values = FALSE,
                     min.diff = 1E-3, 
                     pdf = "rw2"
                 ),

                 crw2 = list(
                     doc = "Exact solution to the random walk of order 2", 
                     hyper = list(
                         theta = list(
                             hyperid =  6001,
                             name = "log precision",
                             short.name = "prec",
                             prior = "loggamma",
                             param = c(1, 0.00005),
                             initial = 4,
                             fixed = FALSE,
                             to.theta = function(x) log(x),
                             from.theta = function(x) exp(x)
                         )
                     ),
                     constr = TRUE,
                     nrow.ncol = FALSE,
                     augmented = FALSE,
                     aug.factor = 2L,
                     aug.constr = 1L,
                     n.div.by = NULL,
                     n.required = FALSE,
                     set.default.values = FALSE,
                     min.diff = 1E-3, 
                     pdf = "crw2"
                 ),

                 seasonal = list(
                     doc = "Seasonal model for time series", 
                     hyper = list(
                         theta = list(
                             hyperid =  7001,
                             name = "log precision",
                             short.name = "prec",
                             prior = "loggamma",
                             param = c(1, 0.00005),
                             initial = 4,
                             fixed = FALSE,
                             to.theta = function(x) log(x),
                             from.theta = function(x) exp(x)
                         )
                     ),
                     constr = FALSE,
                     nrow.ncol = FALSE,
                     augmented = FALSE,
                     aug.factor = 1L,
                     aug.constr = NULL,
                     n.div.by = NULL,
                     n.required = FALSE,
                     set.default.values = FALSE,
                     pdf = "seasonal"
                 ),

                 besag = list(
                     doc = "The Besag area model (CAR-model)", 
                     hyper = list(
                         theta = list(
                             hyperid =  8001,
                             name = "log precision",
                             short.name = "prec",
                             prior = "loggamma",
                             param = c(1, 0.00005),
                             initial = 4,
                             fixed = FALSE,
                             to.theta = function(x) log(x),
                             from.theta = function(x) exp(x)
                         )
                     ),
                     constr = TRUE,
                     nrow.ncol = FALSE,
                     augmented = FALSE,
                     aug.factor = 1L,
                     aug.constr = NULL,
                     n.div.by = NULL,
                     n.required = TRUE,
                     set.default.values = TRUE,
                     pdf = "besag"
                 ),

                 besag2 = list(
                     doc = "The shared Besag model", 
                     hyper = list(
                         theta1 = list(
                             hyperid =  9001,
                             name = "log precision",
                             short.name = "prec",
                             prior = "loggamma",
                             param = c(1, 0.00005),
                             initial = 4,
                             fixed = FALSE,
                             to.theta = function(x) log(x),
                             from.theta = function(x) exp(x)
                         ),
                         theta2 = list(
                             hyperid =  9002,
                             name = "scaling parameter",
                             short.name = "a",
                             prior = "loggamma",
                             param = c(10, 10),
                             initial = 0,
                             fixed = FALSE,
                             to.theta = function(x) log(x),
                             from.theta = function(x) exp(x)
                         )
                     ),
                     constr = FALSE,
                     nrow.ncol = FALSE,
                     augmented = FALSE,
                     aug.factor = 1L,
                     aug.constr = c(1L, 2L),
                     n.div.by = 2L,
                     n.required = TRUE,
                     set.default.values = TRUE,
                     pdf = "besag2"
                 ),

                 bym = list(
                     doc = "The BYM-model (Besag-York-Mollier model)", 
                     hyper = list(
                         theta1 = list(
                             hyperid =  10001,
                             name = "log unstructured precision",
                             short.name = "prec.unstruct",
                             prior = "loggamma",
                             param = c(1, 0.0005),
                             initial = 4,
                             fixed = FALSE,
                             to.theta = function(x) log(x),
                             from.theta = function(x) exp(x)
                         ),
                         theta2 = list(
                             hyperid =  10002,
                             name = "log spatial precision",
                             short.name = "prec.spatial",
                             prior = "loggamma",
                             param = c(1, 0.0005),
                             initial = 4,
                             fixed = FALSE,
                             to.theta = function(x) log(x),
                             from.theta = function(x) exp(x)
                         )
                     ),
                     constr = TRUE,
                     nrow.ncol = FALSE,
                     augmented = TRUE,
                     aug.factor = 2L,
                     aug.constr = 2L,
                     n.div.by = NULL,
                     n.required = TRUE,
                     set.default.values = TRUE,
                     pdf = "bym"
                 ),

                 bym2 = list(
                     doc = "The BYM-model with the PC priors", 
                     hyper = list(
                         theta1 = list(
                             hyperid =  11001,
                             name = "log precision",
                             short.name = "prec",
                             prior = "pc.prec",
                             param = c(1, .01),
                             initial = 4,
                             fixed = FALSE,
                             to.theta = function(x) log(x),
                             from.theta = function(x) exp(x)
                         ),
                         theta2 = list(
                             hyperid =  11002,
                             name = "logit phi",
                             short.name = "phi",
                             prior = "pc",
                             param = c(0.5, 0.5),
                             initial = -3,
                             fixed = FALSE,
                             to.theta = function(x) log(x/(1-x)),
                             from.theta = function(x) exp(x)/(1+exp(x))
                         )
                     ),
                     constr = TRUE,
                     nrow.ncol = FALSE,
                     augmented = TRUE,
                     aug.factor = 2L,
                     aug.constr = 2L,
                     n.div.by = NULL,
                     n.required = TRUE,
                     set.default.values = TRUE,
                     status = "experimental", 
                     pdf = "bym2"
                 ),

                 besagproper = list(
                     doc = "A proper version of the Besag model", 
                     hyper = list(
                         theta1 = list(
                             hyperid =  12001,
                             name = "log precision",
                             short.name = "prec",
                             prior = "loggamma",
                             param = c(1, 0.0005),
                             initial = 2,
                             fixed = FALSE,
                             to.theta = function(x) log(x),
                             from.theta = function(x) exp(x)
                         ),
                         theta2 = list(
                             hyperid =  12002,
                             name = "log diagonal",
                             short.name = "diag",
                             prior = "loggamma",
                             param = c(1, 1),
                             initial = 1,
                             fixed = FALSE,
                             to.theta = function(x) log(x),
                             from.theta = function(x) exp(x)
                         )
                     ),
                     constr = FALSE,
                     nrow.ncol = FALSE,
                     augmented = FALSE,
                     aug.factor = 1L,
                     aug.constr = NULL,
                     n.div.by = NULL,
                     n.required = TRUE,
                     set.default.values = TRUE,
                     status = "experimental", 
                     pdf = "besagproper"
                 ),

                 besagproper2 = list(
                     doc = "An alternative proper version of the Besag model", 
                     hyper = list(
                         theta1 = list(
                             hyperid =  13001,
                             name = "log precision",
                             short.name = "prec",
                             prior = "loggamma",
                             param = c(1, 0.0005),
                             initial = 2,
                             fixed = FALSE,
                             to.theta = function(x) log(x),
                             from.theta = function(x) exp(x)
                         ),
                         theta2 = list(
                             hyperid =  13002,
                             name = "logit lambda",
                             short.name = "lambda",
                             prior = "gaussian",
                             param = c(0, 0.45),
                             initial = 3,
                             fixed = FALSE,
                             to.theta = function(x) log(x/(1-x)),
                             from.theta = function(x) exp(x)/(1+exp(x))
                         )
                     ),
                     constr = FALSE,
                     nrow.ncol = FALSE,
                     augmented = FALSE,
                     aug.factor = 1L,
                     aug.constr = NULL,
                     n.div.by = NULL,
                     n.required = TRUE,
                     set.default.values = TRUE,
                     status = "experimental", 
                     pdf = "besagproper2"
                 ),

                 fgn = list(
                     doc = "Fractional Gaussian noise model", 
                     hyper = list(
                         theta1 = list(
                             hyperid =  13101,
                             name = "log precision",
                             short.name = "prec",
                             prior = "pc.prec",
                             param = c(3, 0.01),
                             initial = 1,
                             fixed = FALSE,
                             to.theta = function(x) log(x),
                             from.theta = function(x) exp(x)
                         ),
                         theta2 = list(
                             hyperid =  13102,
                             name = "logit H",
                             short.name = "H",
                             prior = "pcfgnh",
                             param = c(0.9, 0.1),
                             initial = 2,
                             fixed = FALSE,
                             to.theta = function(x) log((2*x-1)/(2*(1-x))), 
                             from.theta = function(x) 0.5 + 0.5*exp(x)/(1+exp(x))
                         )
                     ),
                     constr = FALSE,
                     nrow.ncol = FALSE,
                     augmented = TRUE,
                     aug.factor = 5L,
                     aug.constr = 1L,
                     n.div.by = NULL,
                     n.required = FALSE,
                     set.default.values = TRUE,
                     order.default = 4L,     ## default order for approximation
                     order.defined = 3L:4L,  ## the list of orders which are implemented
                     pdf = "fgn"
                 ),

                 fgn2 = list(
                     doc = "Fractional Gaussian noise model (alt 2)", 
                     hyper = list(
                         theta1 = list(
                             hyperid =  13111,
                             name = "log precision",
                             short.name = "prec",
                             prior = "pc.prec",
                             param = c(3, 0.01),
                             initial = 1,
                             fixed = FALSE,
                             to.theta = function(x) log(x),
                             from.theta = function(x) exp(x)
                         ),
                         theta2 = list(
                             hyperid =  13112,
                             name = "logit H",
                             short.name = "H",
                             prior = "pcfgnh",
                             param = c(0.9, 0.1),
                             initial = 2,
                             fixed = FALSE,
                             to.theta = function(x) log((2*x-1)/(2*(1-x))), 
                             from.theta = function(x) 0.5 + 0.5*exp(x)/(1+exp(x))
                         )
                     ),
                     constr = FALSE,
                     nrow.ncol = FALSE,
                     augmented = TRUE,
                     aug.factor = 4L,
                     aug.constr = 1L,
                     n.div.by = NULL,
                     n.required = FALSE,
                     set.default.values = TRUE,
                     order.default = 4L,     ## default order for approximation
                     order.defined = 3L:4L,  ## the list of orders which are implemented
                     pdf = "fgn"
                 ),

                 ar1 = list(
                     doc = "Auto-regressive model of order 1 (AR(1))", 
                     hyper = list(
                         theta1 = list(
                             hyperid =  14001,
                             name = "log precision",
                             short.name = "prec",
                             prior = "loggamma",
                             param = c(1, 0.00005),
                             initial = 4,
                             fixed = FALSE,
                             to.theta = function(x) log(x),
                             from.theta = function(x) exp(x)
                         ),
                         theta2 = list(
                             hyperid =  14002,
                             name = "logit lag one correlation",
                             short.name = "rho",
                             prior = "normal",
                             param = c(0, 0.15),
                             initial = 2,
                             fixed = FALSE,
                             to.theta = function(x) log((1+x)/(1-x)),
                             from.theta = function(x) 2*exp(x)/(1+exp(x))-1
                         ), 
                         theta3 = list(
                             hyperid =  14003,
                             name = "mean",
                             short.name = "mean",
                             prior = "normal",
                             param = c(0, 1),
                             initial = 0,
                             fixed = TRUE,
                             to.theta = function(x) x, 
                             from.theta = function(x) x
                         )
                     ),
                     constr = FALSE,
                     nrow.ncol = FALSE,
                     augmented = FALSE,
                     aug.factor = 1L,
                     aug.constr = NULL,
                     n.div.by = NULL,
                     n.required = FALSE,
                     set.default.values = FALSE,
                     pdf = "ar1"
                 ),

                 ar1c = list(
                     doc = "Auto-regressive model of order 1 w/covariates", 
                     hyper = list(
                         theta1 = list(
                             hyperid =  14101,
                             name = "log precision",
                             short.name = "prec",
                             prior = "pc.prec",
                             param = c(1, 0.01),
                             initial = 4,
                             fixed = FALSE,
                             to.theta = function(x) log(x),
                             from.theta = function(x) exp(x)
                         ),
                         theta2 = list(
                             hyperid =  14102,
                             name = "logit lag one correlation",
                             short.name = "rho",
                             prior = "pc.cor0",
                             param = c(0.5, 0.5),
                             initial = 2,
                             fixed = FALSE,
                             to.theta = function(x) log((1+x)/(1-x)),
                             from.theta = function(x) 2*exp(x)/(1+exp(x))-1
                         )
                     ),
                     constr = FALSE,
                     nrow.ncol = FALSE,
                     augmented = FALSE,
                     aug.factor = 1L,
                     aug.constr = NULL,
                     n.div.by = NULL,
                     n.required = FALSE,
                     set.default.values = TRUE,
                     status = "experimental", 
                     pdf = "ar1c"
                 ),

                 ar = list(
                     doc = "Auto-regressive model of order p (AR(p))", 
                     ## to many parameters here, but ...
                     hyper = list(
                         theta1 = list(
                             hyperid =  15001,
                             name = "log precision",
                             short.name = "prec",
                             initial = 4,
                             fixed = FALSE,
                             prior = "pc.prec",
                             param = c(3, 0.01),
                             to.theta = function(x) log(x),
                             from.theta = function(x) exp(x)
                         ),
                         theta2 = list(
                             hyperid =  15002,
                             name = "pacf1",
                             short.name = "pacf1",
                             initial = 1,
                             fixed = FALSE,
                             prior = "pc.cor0",
                             param = c(0.5, 0.5), 
                             to.theta = function(x) log((1+x)/(1-x)),
                             from.theta = function(x) 2*exp(x)/(1+exp(x))-1
                         ),
                         theta3 = list(
                             hyperid =  15003,
                             name = "pacf2",
                             short.name = "pacf2",
                             initial = 0,
                             fixed = FALSE,
                             prior = "pc.cor0",
                             param = c(0.5, 0.4), 
                             to.theta = function(x) log((1+x)/(1-x)),
                             from.theta = function(x) 2*exp(x)/(1+exp(x))-1
                         ),
                         theta4 = list(
                             hyperid =  15004,
                             name = "pacf3",
                             short.name = "pacf3",
                             initial = 0,
                             fixed = FALSE,
                             prior = "pc.cor0",
                             param = c(0.5, 0.3), 
                             to.theta = function(x) log((1+x)/(1-x)),
                             from.theta = function(x) 2*exp(x)/(1+exp(x))-1
                         ),
                         theta5 = list(
                             hyperid =  15005,
                             name = "pacf4",
                             short.name = "pacf4",
                             initial = 0,
                             fixed = FALSE,
                             prior = "pc.cor0",
                             param = c(0.5, 0.2), 
                             to.theta = function(x) log((1+x)/(1-x)),
                             from.theta = function(x) 2*exp(x)/(1+exp(x))-1
                         ),
                         theta6 = list(
                             hyperid =  15006,
                             name = "pacf5",
                             short.name = "pacf5",
                             initial = 0,
                             fixed = FALSE,
                             prior = "pc.cor0",
                             param = c(0.5, 0.1), 
                             to.theta = function(x) log((1+x)/(1-x)),
                             from.theta = function(x) 2*exp(x)/(1+exp(x))-1
                         ),
                         theta7 = list(
                             hyperid =  15007,
                             name = "pacf6",
                             short.name = "pacf6",
                             initial = 0,
                             fixed = FALSE,
                             prior = "pc.cor0",
                             param = c(0.5, 0.1), 
                             to.theta = function(x) log((1+x)/(1-x)),
                             from.theta = function(x) 2*exp(x)/(1+exp(x))-1
                         ),
                         theta8 = list(
                             hyperid =  15008,
                             name = "pacf7",
                             short.name = "pacf7",
                             initial = 0,
                             fixed = FALSE,
                             prior = "pc.cor0",
                             param = c(0.5, 0.1), 
                             to.theta = function(x) log((1+x)/(1-x)),
                             from.theta = function(x) 2*exp(x)/(1+exp(x))-1
                         ),
                         theta9 = list(
                             hyperid =  15009,
                             name = "pacf8",
                             short.name = "pacf8",
                             initial = 0,
                             fixed = FALSE,
                             prior = "pc.cor0",
                             param = c(0.5, 0.1), 
                             to.theta = function(x) log((1+x)/(1-x)),
                             from.theta = function(x) 2*exp(x)/(1+exp(x))-1
                         ),
                         theta10 = list(
                             hyperid =  15010,
                             name = "pacf9",
                             short.name = "pacf9",
                             initial = 0,
                             fixed = FALSE,
                             prior = "pc.cor0",
                             param = c(0.5, 0.1), 
                             to.theta = function(x) log((1+x)/(1-x)),
                             from.theta = function(x) 2*exp(x)/(1+exp(x))-1
                         ),
                         theta11 = list(
                             hyperid =  15011,
                             name = "pacf10",
                             short.name = "pacf10",
                             initial = 0,
                             fixed = FALSE,
                             prior = "pc.cor0",
                             param = c(0.5, 0.1), 
                             to.theta = function(x) log((1+x)/(1-x)),
                             from.theta = function(x) 2*exp(x)/(1+exp(x))-1
                         )
                     ),
                     constr = FALSE,
                     nrow.ncol = FALSE,
                     augmented = FALSE,
                     aug.factor = 1L,
                     aug.constr = NULL,
                     n.div.by = NULL,
                     n.required = FALSE,
                     set.default.values = FALSE,
                     pdf = "ar"
                 ),

                 ou = list(
                     doc = "The Ornstein-Uhlenbeck process", 
                     hyper = list(
                         theta1 = list(
                             hyperid =  16001,
                             name = "log precision",
                             short.name = "prec",
                             prior = "loggamma",
                             param = c(1, 0.00005),
                             initial = 4,
                             fixed = FALSE,
                             to.theta = function(x) log(x),
                             from.theta = function(x) exp(x)
                         ),
                         theta2 = list(
                             hyperid =  16002,
                             name = "log phi",
                             short.name = "phi",
                             prior = "normal",
                             param = c(0, .2),
                             initial = -1,
                             fixed = FALSE,
                             to.theta = function(x) log(x), 
                             from.theta = function(x) exp(x)
                         )
                     ),
                     constr = FALSE,
                     nrow.ncol = FALSE,
                     augmented = FALSE,
                     aug.factor = 1L,
                     aug.constr = NULL,
                     n.div.by = NULL,
                     n.required = FALSE,
                     set.default.values = FALSE,
                     pdf = "ou"
                 ),

                 intslope = list(
                     doc = "Intecept-slope model with Wishart-prior", 
                     hyper = list(
                         theta1 = list(
                             hyperid =  16101,
                             name = "log precision1",
                             short.name = "prec1",
                             initial = 4,
                             fixed = FALSE,
                             prior = "wishart2d",
                             param = c(4, 1, 1, 0),
                             to.theta = function(x) log(x),
                             from.theta = function(x) exp(x)
                         ),
                         theta2 = list(
                             hyperid =  16102,
                             name = "log precision2",
                             short.name = "prec2",
                             initial = 4,
                             fixed = FALSE,
                             prior = "none",
                             param = numeric(),
                             to.theta = function(x) log(x),
                             from.theta = function(x) exp(x)
                         ),
                         theta3 = list(
                             hyperid =  16103,
                             name = "logit correlation",
                             short.name = "cor",
                             initial = 4,
                             fixed = FALSE,
                             prior = "none",
                             param = numeric(),
                             to.theta = function(x) log((1+x)/(1-x)),
                             from.theta = function(x) 2*exp(x)/(1+exp(x))-1
                         ), 
                         theta4 = list(
                             hyperid =  16104,
                             name = "gamma1",
                             short.name = "g1",
                             initial = 1,
                             fixed = TRUE,   ## YES.
                             prior = "normal",
                             param = c(1, 36),
                             to.theta = function(x) x,
                             from.theta = function(x) x
                         ),
                         theta5 = list(
                             hyperid =  16105,
                             name = "gamma2",
                             short.name = "g2",
                             initial = 1,
                             fixed = FALSE,
                             prior = "normal",
                             param = c(1, 36),
                             to.theta = function(x) x,
                             from.theta = function(x) x
                         ),
                         theta6 = list(
                             hyperid =  16106,
                             name = "gamma3",
                             short.name = "g3",
                             initial = 1,
                             fixed = FALSE,
                             prior = "normal",
                             param = c(1, 36),
                             to.theta = function(x) x,
                             from.theta = function(x) x
                         ),
                         theta7 = list(
                             hyperid =  16107,
                             name = "gamma4",
                             short.name = "g4",
                             initial = 1,
                             fixed = FALSE,
                             prior = "normal",
                             param = c(1, 36),
                             to.theta = function(x) x,
                             from.theta = function(x) x
                         ),
                         theta8 = list(
                             hyperid =  16108,
                             name = "gamma5",
                             short.name = "g5",
                             initial = 1,
                             fixed = FALSE,
                             prior = "normal",
                             param = c(1, 36),
                             to.theta = function(x) x,
                             from.theta = function(x) x
                         ),
                         theta9 = list(
                             hyperid =  16109,
                             name = "gamma6",
                             short.name = "g6",
                             initial = 1,
                             fixed = FALSE,
                             prior = "normal",
                             param = c(1, 36),
                             to.theta = function(x) x,
                             from.theta = function(x) x
                         ),
                         theta10 = list(
                             hyperid =  16110,
                             name = "gamma7",
                             short.name = "g7",
                             initial = 1,
                             fixed = FALSE,
                             prior = "normal",
                             param = c(1, 36),
                             to.theta = function(x) x,
                             from.theta = function(x) x
                         ),
                         theta11 = list(
                             hyperid =  16111,
                             name = "gamma8",
                             short.name = "g8",
                             initial = 1,
                             fixed = FALSE,
                             prior = "normal",
                             param = c(1, 36),
                             to.theta = function(x) x,
                             from.theta = function(x) x
                         ),
                         theta12 = list(
                             hyperid =  16112,
                             name = "gamma9",
                             short.name = "g9",
                             initial = 1,
                             fixed = FALSE,
                             prior = "normal",
                             param = c(1, 36),
                             to.theta = function(x) x,
                             from.theta = function(x) x
                         ),
                         theta13 = list(
                             hyperid =  16113,
                             name = "gamma10",
                             short.name = "g10",
                             initial = 1,
                             fixed = FALSE,
                             prior = "normal",
                             param = c(1, 36),
                             to.theta = function(x) x,
                             from.theta = function(x) x
                         )
                     ),
                     constr = FALSE,
                     nrow.ncol = FALSE,
                     augmented = FALSE,
                     aug.factor = 1L,
                     aug.constr = NULL,
                     n.div.by = NULL,
                     n.required = FALSE,
                     set.default.values = TRUE,
                     status = "experimental", 
                     pdf = "intslope"
                 ),

                 generic = list(
                     doc = "A generic model", 
                     hyper = list(
                         theta = list(
                             hyperid =  17001,
                             name = "log precision",
                             short.name = "prec",
                             prior = "loggamma",
                             param = c(1, 0.00005),
                             initial = 4,
                             fixed = FALSE,
                             to.theta = function(x) log(x),
                             from.theta = function(x) exp(x)
                         )
                     ),
                     constr = FALSE,
                     nrow.ncol = FALSE,
                     augmented = FALSE,
                     aug.factor = 1L,
                     aug.constr = NULL,
                     n.div.by = NULL,
                     n.required = TRUE,
                     set.default.values = TRUE,
                     pdf = "generic0"
                 ),

                 generic0 = list(
                     doc = "A generic model (type 0)", 
                     hyper = list(
                         theta = list(
                             hyperid =  18001,
                             name = "log precision",
                             short.name = "prec",
                             prior = "loggamma",
                             param = c(1, 0.00005),
                             initial = 4,
                             fixed = FALSE,
                             to.theta = function(x) log(x),
                             from.theta = function(x) exp(x)
                         )
                     ),
                     constr = FALSE,
                     nrow.ncol = FALSE,
                     augmented = FALSE,
                     aug.factor = 1L,
                     aug.constr = NULL,
                     n.div.by = NULL,
                     n.required = TRUE,
                     set.default.values = TRUE,
                     pdf = "generic0"
                 ),

                 generic1 = list(
                     doc = "A generic model (type 1)", 
                     hyper = list(
                         theta1 = list(
                             hyperid =  19001,
                             name = "log precision",
                             short.name = "prec",
                             prior = "loggamma",
                             param = c(1, 0.00005),
                             initial = 4,
                             fixed = FALSE,
                             to.theta = function(x) log(x),
                             from.theta = function(x) exp(x)
                         ),
                         theta2 = list(
                             hyperid =  19002,
                             name = "beta",
                             short.name = "beta",
                             initial = 2,
                             fixed = FALSE,
                             prior = "gaussian",
                             param = c(0, 0.1),
                             to.theta = function(x) log(x/(1-x)),
                             from.theta = function(x) exp(x)/(1+exp(x))
                         )
                     ),
                     constr = FALSE,
                     nrow.ncol = FALSE,
                     augmented = FALSE,
                     aug.factor = 1L,
                     aug.constr = NULL,
                     n.div.by = NULL,
                     n.required = TRUE,
                     set.default.values = TRUE,
                     pdf = "generic1"
                 ),

                 generic2 = list(
                     doc = "A generic model (type 2)", 
                     hyper = list(
                         theta1 = list(
                             hyperid =  20001,
                             name = "log precision cmatrix",
                             short.name = "prec",
                             initial = 4,
                             fixed = FALSE,
                             prior = "loggamma",
                             param = c(1, 0.00005),
                             to.theta = function(x) log(x),
                             from.theta = function(x) exp(x)
                         ),
                         theta2 = list(
                             hyperid =  20002,
                             name = "log precision random",
                             short.name = "prec.random",
                             initial = 4,
                             fixed = FALSE,
                             prior = "loggamma",
                             param = c(1, 0.001),
                             to.theta = function(x) log(x),
                             from.theta = function(x) exp(x)
                         )
                     ),
                     constr = FALSE,
                     nrow.ncol = FALSE,
                     augmented = FALSE,
                     aug.factor = 2L,
                     aug.constr = 2L,
                     n.div.by = NULL,
                     n.required = TRUE,
                     set.default.values = TRUE,
                     pdf = "generic2"
                 ),

                 generic3 = list(
                     doc = "A generic model (type 3)", 
                     hyper = list(
                         theta1 = list(
                             hyperid =  21001,
                             name = "log precision1",
                             short.name = "prec1",
                             initial = 4,
                             fixed = FALSE,
                             prior = "loggamma",
                             param = c(1, 0.00005),
                             to.theta = function(x) log(x),
                             from.theta = function(x) exp(x)
                         ),
                         theta2 = list(
                             hyperid =  21002,
                             name = "log precision2",
                             short.name = "prec2",
                             initial = 4,
                             fixed = FALSE,
                             prior = "loggamma",
                             param = c(1, 0.00005),
                             to.theta = function(x) log(x),
                             from.theta = function(x) exp(x)
                         ),
                         theta3 = list(
                             hyperid =  21003,
                             name = "log precision3",
                             short.name = "prec3",
                             initial = 4,
                             fixed = FALSE,
                             prior = "loggamma",
                             param = c(1, 0.00005),
                             to.theta = function(x) log(x),
                             from.theta = function(x) exp(x)
                         ),
                         theta4 = list(
                             hyperid =  21004,
                             name = "log precision4",
                             short.name = "prec4",
                             initial = 4,
                             fixed = FALSE,
                             prior = "loggamma",
                             param = c(1, 0.00005),
                             to.theta = function(x) log(x),
                             from.theta = function(x) exp(x)
                         ),
                         theta5 = list(
                             hyperid =  21005,
                             name = "log precision5",
                             short.name = "prec5",
                             initial = 4,
                             fixed = FALSE,
                             prior = "loggamma",
                             param = c(1, 0.00005),
                             to.theta = function(x) log(x),
                             from.theta = function(x) exp(x)
                         ),
                         theta6 = list(
                             hyperid =  21006,
                             name = "log precision6",
                             short.name = "prec6",
                             initial = 4,
                             fixed = FALSE,
                             prior = "loggamma",
                             param = c(1, 0.00005),
                             to.theta = function(x) log(x),
                             from.theta = function(x) exp(x)
                         ),
                         theta7 = list(
                             hyperid =  21007,
                             name = "log precision7",
                             short.name = "prec7",
                             initial = 4,
                             fixed = FALSE,
                             prior = "loggamma",
                             param = c(1, 0.00005),
                             to.theta = function(x) log(x),
                             from.theta = function(x) exp(x)
                         ),
                         theta8 = list(
                             hyperid =  21008,
                             name = "log precision8",
                             short.name = "prec8",
                             initial = 4,
                             fixed = FALSE,
                             prior = "loggamma",
                             param = c(1, 0.00005),
                             to.theta = function(x) log(x),
                             from.theta = function(x) exp(x)
                         ),
                         theta9 = list(
                             hyperid =  21009,
                             name = "log precision9",
                             short.name = "prec9",
                             initial = 4,
                             fixed = FALSE,
                             prior = "loggamma",
                             param = c(1, 0.00005),
                             to.theta = function(x) log(x),
                             from.theta = function(x) exp(x)
                         ),
                         theta10 = list(
                             hyperid =  21010,
                             name = "log precision10",
                             short.name = "prec10",
                             initial = 4,
                             fixed = FALSE,
                             prior = "loggamma",
                             param = c(1, 0.00005),
                             to.theta = function(x) log(x),
                             from.theta = function(x) exp(x)
                         ),
                         theta11 = list(
                             hyperid =  21011,
                             name = "log precision common",
                             short.name = "prec.common",
                             initial = 0, ## yes!
                             fixed = TRUE, ## yes!
                             prior = "loggamma",
                             param = c(1, 0.00005),
                             to.theta = function(x) log(x),
                             from.theta = function(x) exp(x)
                         )
                     ),
                     constr = FALSE,
                     nrow.ncol = FALSE,
                     augmented = FALSE,
                     aug.factor = 1L,
                     aug.constr = NULL,
                     n.div.by = NULL,
                     n.required = TRUE,
                     set.default.values = TRUE,
                     status = "experimental", 
                     pdf = "generic3"
                 ),

                 spde = list(
                     doc = "A SPDE model", 
                     ## this will be redone anyway soon....
                     hyper = list(
                         theta1 = list(
                             hyperid =  22001,
                             name = "theta.T",
                             short.name = "T",
                             initial = 2,
                             fixed = FALSE,
                             prior = "normal",
                             param = c(0, 1),
                             to.theta = function(x) x,
                             from.theta = function(x) x
                         ),
                         theta2 = list(
                             hyperid =  22002,
                             name = "theta.K",
                             short.name = "K",
                             initial = -2,
                             fixed = FALSE,
                             prior = "normal",
                             param = c(0, 1),
                             to.theta = function(x) x,
                             from.theta = function(x) x
                         ),
                         theta3 = list(
                             hyperid =  22003,
                             name = "theta.KT",
                             short.name = "KT",
                             initial = 0,
                             fixed = FALSE,
                             prior = "normal",
                             param = c(0, 1),
                             to.theta = function(x) x,
                             from.theta = function(x) x
                         ),
                         theta4 = list(
                             hyperid =  22004,
                             name = "theta.OC",
                             short.name = "OC",
                             initial = -20,
                             fixed = TRUE,
                             prior = "normal",
                             param = c(0, 0.2),
                             to.theta = function(x) log(x/(1-x)),
                             from.theta = function(x) exp(x)/(1+exp(x))
                         )
                     ),
                     constr = FALSE,
                     nrow.ncol = FALSE,
                     augmented = FALSE,
                     aug.factor = 1L,
                     aug.constr = NULL,
                     n.div.by = NULL,
                     n.required = TRUE,
                     set.default.values = TRUE,
                     pdf = "spde"
                 ),

                 spde2 = list(
                     doc = "A SPDE2 model", 
                     ## to many parameters here, but ...
                     hyper = list(
                         theta1 = list(
                             hyperid =  23001,
                             name = "theta1",
                             short.name = "t1",
                             initial = 0,
                             fixed = FALSE,
                             prior = "mvnorm",
                             param = c(1, 1),
                             to.theta = function(x) x,
                             from.theta = function(x) x
                         ),
                         theta2 = list(
                             hyperid =  23002,
                             name = "theta2",
                             short.name = "t2",
                             initial = 0,
                             fixed = FALSE,
                             prior = "none",
                             param = numeric(0),
                             to.theta = function(x) x,
                             from.theta = function(x) x
                         ),
                         theta3 = list(
                             hyperid =  23003,
                             name = "theta3",
                             short.name = "t3",
                             initial = 0,
                             fixed = FALSE,
                             prior = "none",
                             param = numeric(0),
                             to.theta = function(x) x,
                             from.theta = function(x) x
                         ),
                         theta4 = list(
                             hyperid =  23004,
                             name = "theta4",
                             short.name = "t4",
                             initial = 0,
                             fixed = FALSE,
                             prior = "none",
                             param = numeric(0),
                             to.theta = function(x) x,
                             from.theta = function(x) x
                         ),
                         theta5 = list(
                             hyperid =  23005,
                             name = "theta5",
                             short.name = "t5",
                             initial = 0,
                             fixed = FALSE,
                             prior = "none",
                             param = numeric(0),
                             to.theta = function(x) x,
                             from.theta = function(x) x
                         ),
                         theta6 = list(
                             hyperid =  23006,
                             name = "theta6",
                             short.name = "t6",
                             initial = 0,
                             fixed = FALSE,
                             prior = "none",
                             param = numeric(0),
                             to.theta = function(x) x,
                             from.theta = function(x) x
                         ),
                         theta7 = list(
                             hyperid =  23007,
                             name = "theta7",
                             short.name = "t7",
                             initial = 0,
                             fixed = FALSE,
                             prior = "none",
                             param = numeric(0),
                             to.theta = function(x) x,
                             from.theta = function(x) x
                         ),
                         theta8 = list(
                             hyperid =  23008,
                             name = "theta8",
                             short.name = "t8",
                             initial = 0,
                             fixed = FALSE,
                             prior = "none",
                             param = numeric(0),
                             to.theta = function(x) x,
                             from.theta = function(x) x
                         ),
                         theta9 = list(
                             hyperid =  23009,
                             name = "theta9",
                             short.name = "t9",
                             initial = 0,
                             fixed = FALSE,
                             prior = "none",
                             param = numeric(0),
                             to.theta = function(x) x,
                             from.theta = function(x) x
                         ),
                         theta10 = list(
                             hyperid =  23010,
                             name = "theta10",
                             short.name = "t10",
                             initial = 0,
                             fixed = FALSE,
                             prior = "none",
                             param = numeric(0),
                             to.theta = function(x) x,
                             from.theta = function(x) x
                         ),
                         theta11 = list(
                             hyperid =  23011,
                             name = "theta11",
                             short.name = "t11",
                             initial = 0,
                             fixed = FALSE,
                             prior = "none",
                             param = numeric(0),
                             to.theta = function(x) x,
                             from.theta = function(x) x
                         ),
                         theta12 = list(
                             hyperid =  23012,
                             name = "theta12",
                             short.name = "t12",
                             initial = 0,
                             fixed = FALSE,
                             prior = "none",
                             param = numeric(0),
                             to.theta = function(x) x,
                             from.theta = function(x) x
                         ),
                         theta13 = list(
                             hyperid =  23013,
                             name = "theta13",
                             short.name = "t13",
                             initial = 0,
                             fixed = FALSE,
                             prior = "none",
                             param = numeric(0),
                             to.theta = function(x) x,
                             from.theta = function(x) x
                         ),
                         theta14 = list(
                             hyperid =  23014,
                             name = "theta14",
                             short.name = "t14",
                             initial = 0,
                             fixed = FALSE,
                             prior = "none",
                             param = numeric(0),
                             to.theta = function(x) x,
                             from.theta = function(x) x
                         ),
                         theta15 = list(
                             hyperid =  23015,
                             name = "theta15",
                             short.name = "t15",
                             initial = 0,
                             fixed = FALSE,
                             prior = "none",
                             param = numeric(0),
                             to.theta = function(x) x,
                             from.theta = function(x) x
                         ),
                         theta16 = list(
                             hyperid =  23016,
                             name = "theta16",
                             short.name = "t16",
                             initial = 0,
                             fixed = FALSE,
                             prior = "none",
                             param = numeric(0),
                             to.theta = function(x) x,
                             from.theta = function(x) x
                         ),
                         theta17 = list(
                             hyperid =  23017,
                             name = "theta17",
                             short.name = "t17",
                             initial = 0,
                             fixed = FALSE,
                             prior = "none",
                             param = numeric(0),
                             to.theta = function(x) x,
                             from.theta = function(x) x
                         ),
                         theta18 = list(
                             hyperid =  23018,
                             name = "theta18",
                             short.name = "t18",
                             initial = 0,
                             fixed = FALSE,
                             prior = "none",
                             param = numeric(0),
                             to.theta = function(x) x,
                             from.theta = function(x) x
                         ),
                         theta19 = list(
                             hyperid =  23019,
                             name = "theta19",
                             short.name = "t19",
                             initial = 0,
                             fixed = FALSE,
                             prior = "none",
                             param = numeric(0),
                             to.theta = function(x) x,
                             from.theta = function(x) x
                         ),
                         theta20 = list(
                             hyperid =  23020,
                             name = "theta20",
                             short.name = "t20",
                             initial = 0,
                             fixed = FALSE,
                             prior = "none",
                             param = numeric(0),
                             to.theta = function(x) x,
                             from.theta = function(x) x
                         ),
                         theta21 = list(
                             hyperid =  23021,
                             name = "theta21",
                             short.name = "t21",
                             initial = 0,
                             fixed = FALSE,
                             prior = "none",
                             param = numeric(0),
                             to.theta = function(x) x,
                             from.theta = function(x) x
                         ),
                         theta22 = list(
                             hyperid =  23022,
                             name = "theta22",
                             short.name = "t22",
                             initial = 0,
                             fixed = FALSE,
                             prior = "none",
                             param = numeric(0),
                             to.theta = function(x) x,
                             from.theta = function(x) x
                         ),
                         theta23 = list(
                             hyperid =  23023,
                             name = "theta23",
                             short.name = "t23",
                             initial = 0,
                             fixed = FALSE,
                             prior = "none",
                             param = numeric(0),
                             to.theta = function(x) x,
                             from.theta = function(x) x
                         ),
                         theta24 = list(
                             hyperid =  23024,
                             name = "theta24",
                             short.name = "t24",
                             initial = 0,
                             fixed = FALSE,
                             prior = "none",
                             param = numeric(0),
                             to.theta = function(x) x,
                             from.theta = function(x) x
                         ),
                         theta25 = list(
                             hyperid =  23025,
                             name = "theta25",
                             short.name = "t25",
                             initial = 0,
                             fixed = FALSE,
                             prior = "none",
                             param = numeric(0),
                             to.theta = function(x) x,
                             from.theta = function(x) x
                         ),
                         theta26 = list(
                             hyperid =  23026,
                             name = "theta26",
                             short.name = "t26",
                             initial = 0,
                             fixed = FALSE,
                             prior = "none",
                             param = numeric(0),
                             to.theta = function(x) x,
                             from.theta = function(x) x
                         ),
                         theta27 = list(
                             hyperid =  23027,
                             name = "theta27",
                             short.name = "t27",
                             initial = 0,
                             fixed = FALSE,
                             prior = "none",
                             param = numeric(0),
                             to.theta = function(x) x,
                             from.theta = function(x) x
                         ),
                         theta28 = list(
                             hyperid =  23028,
                             name = "theta28",
                             short.name = "t28",
                             initial = 0,
                             fixed = FALSE,
                             prior = "none",
                             param = numeric(0),
                             to.theta = function(x) x,
                             from.theta = function(x) x
                         ),
                         theta29 = list(
                             hyperid =  23029,
                             name = "theta29",
                             short.name = "t29",
                             initial = 0,
                             fixed = FALSE,
                             prior = "none",
                             param = numeric(0),
                             to.theta = function(x) x,
                             from.theta = function(x) x
                         ),
                         theta30 = list(
                             hyperid =  23030,
                             name = "theta30",
                             short.name = "t30",
                             initial = 0,
                             fixed = FALSE,
                             prior = "none",
                             param = numeric(0),
                             to.theta = function(x) x,
                             from.theta = function(x) x
                         ),
                         theta31 = list(
                             hyperid =  23031,
                             name = "theta31",
                             short.name = "t31",
                             initial = 0,
                             fixed = FALSE,
                             prior = "none",
                             param = numeric(0),
                             to.theta = function(x) x,
                             from.theta = function(x) x
                         ),
                         theta32 = list(
                             hyperid =  23032,
                             name = "theta32",
                             short.name = "t32",
                             initial = 0,
                             fixed = FALSE,
                             prior = "none",
                             param = numeric(0),
                             to.theta = function(x) x,
                             from.theta = function(x) x
                         ),
                         theta33 = list(
                             hyperid =  23033,
                             name = "theta33",
                             short.name = "t33",
                             initial = 0,
                             fixed = FALSE,
                             prior = "none",
                             param = numeric(0),
                             to.theta = function(x) x,
                             from.theta = function(x) x
                         ),
                         theta34 = list(
                             hyperid =  23034,
                             name = "theta34",
                             short.name = "t34",
                             initial = 0,
                             fixed = FALSE,
                             prior = "none",
                             param = numeric(0),
                             to.theta = function(x) x,
                             from.theta = function(x) x
                         ),
                         theta35 = list(
                             hyperid =  23035,
                             name = "theta35",
                             short.name = "t35",
                             initial = 0,
                             fixed = FALSE,
                             prior = "none",
                             param = numeric(0),
                             to.theta = function(x) x,
                             from.theta = function(x) x
                         ),
                         theta36 = list(
                             hyperid =  23036,
                             name = "theta36",
                             short.name = "t36",
                             initial = 0,
                             fixed = FALSE,
                             prior = "none",
                             param = numeric(0),
                             to.theta = function(x) x,
                             from.theta = function(x) x
                         ),
                         theta37 = list(
                             hyperid =  23037,
                             name = "theta37",
                             short.name = "t37",
                             initial = 0,
                             fixed = FALSE,
                             prior = "none",
                             param = numeric(0),
                             to.theta = function(x) x,
                             from.theta = function(x) x
                         ),
                         theta38 = list(
                             hyperid =  23038,
                             name = "theta38",
                             short.name = "t38",
                             initial = 0,
                             fixed = FALSE,
                             prior = "none",
                             param = numeric(0),
                             to.theta = function(x) x,
                             from.theta = function(x) x
                         ),
                         theta39 = list(
                             hyperid =  23039,
                             name = "theta39",
                             short.name = "t39",
                             initial = 0,
                             fixed = FALSE,
                             prior = "none",
                             param = numeric(0),
                             to.theta = function(x) x,
                             from.theta = function(x) x
                         ),
                         theta40 = list(
                             hyperid =  23040,
                             name = "theta40",
                             short.name = "t40",
                             initial = 0,
                             fixed = FALSE,
                             prior = "none",
                             param = numeric(0),
                             to.theta = function(x) x,
                             from.theta = function(x) x
                         ),
                         theta41 = list(
                             hyperid =  23041,
                             name = "theta41",
                             short.name = "t41",
                             initial = 0,
                             fixed = FALSE,
                             prior = "none",
                             param = numeric(0),
                             to.theta = function(x) x,
                             from.theta = function(x) x
                         ),
                         theta42 = list(
                             hyperid =  23042,
                             name = "theta42",
                             short.name = "t42",
                             initial = 0,
                             fixed = FALSE,
                             prior = "none",
                             param = numeric(0),
                             to.theta = function(x) x,
                             from.theta = function(x) x
                         ),
                         theta43 = list(
                             hyperid =  23043,
                             name = "theta43",
                             short.name = "t43",
                             initial = 0,
                             fixed = FALSE,
                             prior = "none",
                             param = numeric(0),
                             to.theta = function(x) x,
                             from.theta = function(x) x
                         ),
                         theta44 = list(
                             hyperid =  23044,
                             name = "theta44",
                             short.name = "t44",
                             initial = 0,
                             fixed = FALSE,
                             prior = "none",
                             param = numeric(0),
                             to.theta = function(x) x,
                             from.theta = function(x) x
                         ),
                         theta45 = list(
                             hyperid =  23045,
                             name = "theta45",
                             short.name = "t45",
                             initial = 0,
                             fixed = FALSE,
                             prior = "none",
                             param = numeric(0),
                             to.theta = function(x) x,
                             from.theta = function(x) x
                         ),
                         theta46 = list(
                             hyperid =  23046,
                             name = "theta46",
                             short.name = "t46",
                             initial = 0,
                             fixed = FALSE,
                             prior = "none",
                             param = numeric(0),
                             to.theta = function(x) x,
                             from.theta = function(x) x
                         ),
                         theta47 = list(
                             hyperid =  23047,
                             name = "theta47",
                             short.name = "t47",
                             initial = 0,
                             fixed = FALSE,
                             prior = "none",
                             param = numeric(0),
                             to.theta = function(x) x,
                             from.theta = function(x) x
                         ),
                         theta48 = list(
                             hyperid =  23048,
                             name = "theta48",
                             short.name = "t48",
                             initial = 0,
                             fixed = FALSE,
                             prior = "none",
                             param = numeric(0),
                             to.theta = function(x) x,
                             from.theta = function(x) x
                         ),
                         theta49 = list(
                             hyperid =  23049,
                             name = "theta49",
                             short.name = "t49",
                             initial = 0,
                             fixed = FALSE,
                             prior = "none",
                             param = numeric(0),
                             to.theta = function(x) x,
                             from.theta = function(x) x
                         ),
                         theta50 = list(
                             hyperid =  23050,
                             name = "theta50",
                             short.name = "t50",
                             initial = 0,
                             fixed = FALSE,
                             prior = "none",
                             param = numeric(0),
                             to.theta = function(x) x,
                             from.theta = function(x) x
                         ),
                         theta51 = list(
                             hyperid =  23051,
                             name = "theta51",
                             short.name = "t51",
                             initial = 0,
                             fixed = FALSE,
                             prior = "none",
                             param = numeric(0),
                             to.theta = function(x) x,
                             from.theta = function(x) x
                         ),
                         theta52 = list(
                             hyperid =  23052,
                             name = "theta52",
                             short.name = "t52",
                             initial = 0,
                             fixed = FALSE,
                             prior = "none",
                             param = numeric(0),
                             to.theta = function(x) x,
                             from.theta = function(x) x
                         ),
                         theta53 = list(
                             hyperid =  23053,
                             name = "theta53",
                             short.name = "t53",
                             initial = 0,
                             fixed = FALSE,
                             prior = "none",
                             param = numeric(0),
                             to.theta = function(x) x,
                             from.theta = function(x) x
                         ),
                         theta54 = list(
                             hyperid =  23054,
                             name = "theta54",
                             short.name = "t54",
                             initial = 0,
                             fixed = FALSE,
                             prior = "none",
                             param = numeric(0),
                             to.theta = function(x) x,
                             from.theta = function(x) x
                         ),
                         theta55 = list(
                             hyperid =  23055,
                             name = "theta55",
                             short.name = "t55",
                             initial = 0,
                             fixed = FALSE,
                             prior = "none",
                             param = numeric(0),
                             to.theta = function(x) x,
                             from.theta = function(x) x
                         ),
                         theta56 = list(
                             hyperid =  23056,
                             name = "theta56",
                             short.name = "t56",
                             initial = 0,
                             fixed = FALSE,
                             prior = "none",
                             param = numeric(0),
                             to.theta = function(x) x,
                             from.theta = function(x) x
                         ),
                         theta57 = list(
                             hyperid =  23057,
                             name = "theta57",
                             short.name = "t57",
                             initial = 0,
                             fixed = FALSE,
                             prior = "none",
                             param = numeric(0),
                             to.theta = function(x) x,
                             from.theta = function(x) x
                         ),
                         theta58 = list(
                             hyperid =  23058,
                             name = "theta58",
                             short.name = "t58",
                             initial = 0,
                             fixed = FALSE,
                             prior = "none",
                             param = numeric(0),
                             to.theta = function(x) x,
                             from.theta = function(x) x
                         ),
                         theta59 = list(
                             hyperid =  23059,
                             name = "theta59",
                             short.name = "t59",
                             initial = 0,
                             fixed = FALSE,
                             prior = "none",
                             param = numeric(0),
                             to.theta = function(x) x,
                             from.theta = function(x) x
                         ),
                         theta60 = list(
                             hyperid =  23060,
                             name = "theta60",
                             short.name = "t60",
                             initial = 0,
                             fixed = FALSE,
                             prior = "none",
                             param = numeric(0),
                             to.theta = function(x) x,
                             from.theta = function(x) x
                         ),
                         theta61 = list(
                             hyperid =  23061,
                             name = "theta61",
                             short.name = "t61",
                             initial = 0,
                             fixed = FALSE,
                             prior = "none",
                             param = numeric(0),
                             to.theta = function(x) x,
                             from.theta = function(x) x
                         ),
                         theta62 = list(
                             hyperid =  23062,
                             name = "theta62",
                             short.name = "t62",
                             initial = 0,
                             fixed = FALSE,
                             prior = "none",
                             param = numeric(0),
                             to.theta = function(x) x,
                             from.theta = function(x) x
                         ),
                         theta63 = list(
                             hyperid =  23063,
                             name = "theta63",
                             short.name = "t63",
                             initial = 0,
                             fixed = FALSE,
                             prior = "none",
                             param = numeric(0),
                             to.theta = function(x) x,
                             from.theta = function(x) x
                         ),
                         theta64 = list(
                             hyperid =  23064,
                             name = "theta64",
                             short.name = "t64",
                             initial = 0,
                             fixed = FALSE,
                             prior = "none",
                             param = numeric(0),
                             to.theta = function(x) x,
                             from.theta = function(x) x
                         ),
                         theta65 = list(
                             hyperid =  23065,
                             name = "theta65",
                             short.name = "t65",
                             initial = 0,
                             fixed = FALSE,
                             prior = "none",
                             param = numeric(0),
                             to.theta = function(x) x,
                             from.theta = function(x) x
                         ),
                         theta66 = list(
                             hyperid =  23066,
                             name = "theta66",
                             short.name = "t66",
                             initial = 0,
                             fixed = FALSE,
                             prior = "none",
                             param = numeric(0),
                             to.theta = function(x) x,
                             from.theta = function(x) x
                         ),
                         theta67 = list(
                             hyperid =  23067,
                             name = "theta67",
                             short.name = "t67",
                             initial = 0,
                             fixed = FALSE,
                             prior = "none",
                             param = numeric(0),
                             to.theta = function(x) x,
                             from.theta = function(x) x
                         ),
                         theta68 = list(
                             hyperid =  23068,
                             name = "theta68",
                             short.name = "t68",
                             initial = 0,
                             fixed = FALSE,
                             prior = "none",
                             param = numeric(0),
                             to.theta = function(x) x,
                             from.theta = function(x) x
                         ),
                         theta69 = list(
                             hyperid =  23069,
                             name = "theta69",
                             short.name = "t69",
                             initial = 0,
                             fixed = FALSE,
                             prior = "none",
                             param = numeric(0),
                             to.theta = function(x) x,
                             from.theta = function(x) x
                         ),
                         theta70 = list(
                             hyperid =  23070,
                             name = "theta70",
                             short.name = "t70",
                             initial = 0,
                             fixed = FALSE,
                             prior = "none",
                             param = numeric(0),
                             to.theta = function(x) x,
                             from.theta = function(x) x
                         ),
                         theta71 = list(
                             hyperid =  23071,
                             name = "theta71",
                             short.name = "t71",
                             initial = 0,
                             fixed = FALSE,
                             prior = "none",
                             param = numeric(0),
                             to.theta = function(x) x,
                             from.theta = function(x) x
                         ),
                         theta72 = list(
                             hyperid =  23072,
                             name = "theta72",
                             short.name = "t72",
                             initial = 0,
                             fixed = FALSE,
                             prior = "none",
                             param = numeric(0),
                             to.theta = function(x) x,
                             from.theta = function(x) x
                         ),
                         theta73 = list(
                             hyperid =  23073,
                             name = "theta73",
                             short.name = "t73",
                             initial = 0,
                             fixed = FALSE,
                             prior = "none",
                             param = numeric(0),
                             to.theta = function(x) x,
                             from.theta = function(x) x
                         ),
                         theta74 = list(
                             hyperid =  23074,
                             name = "theta74",
                             short.name = "t74",
                             initial = 0,
                             fixed = FALSE,
                             prior = "none",
                             param = numeric(0),
                             to.theta = function(x) x,
                             from.theta = function(x) x
                         ),
                         theta75 = list(
                             hyperid =  23075,
                             name = "theta75",
                             short.name = "t75",
                             initial = 0,
                             fixed = FALSE,
                             prior = "none",
                             param = numeric(0),
                             to.theta = function(x) x,
                             from.theta = function(x) x
                         ),
                         theta76 = list(
                             hyperid =  23076,
                             name = "theta76",
                             short.name = "t76",
                             initial = 0,
                             fixed = FALSE,
                             prior = "none",
                             param = numeric(0),
                             to.theta = function(x) x,
                             from.theta = function(x) x
                         ),
                         theta77 = list(
                             hyperid =  23077,
                             name = "theta77",
                             short.name = "t77",
                             initial = 0,
                             fixed = FALSE,
                             prior = "none",
                             param = numeric(0),
                             to.theta = function(x) x,
                             from.theta = function(x) x
                         ),
                         theta78 = list(
                             hyperid =  23078,
                             name = "theta78",
                             short.name = "t78",
                             initial = 0,
                             fixed = FALSE,
                             prior = "none",
                             param = numeric(0),
                             to.theta = function(x) x,
                             from.theta = function(x) x
                         ),
                         theta79 = list(
                             hyperid =  23079,
                             name = "theta79",
                             short.name = "t79",
                             initial = 0,
                             fixed = FALSE,
                             prior = "none",
                             param = numeric(0),
                             to.theta = function(x) x,
                             from.theta = function(x) x
                         ),
                         theta80 = list(
                             hyperid =  23080,
                             name = "theta80",
                             short.name = "t80",
                             initial = 0,
                             fixed = FALSE,
                             prior = "none",
                             param = numeric(0),
                             to.theta = function(x) x,
                             from.theta = function(x) x
                         ),
                         theta81 = list(
                             hyperid =  23081,
                             name = "theta81",
                             short.name = "t81",
                             initial = 0,
                             fixed = FALSE,
                             prior = "none",
                             param = numeric(0),
                             to.theta = function(x) x,
                             from.theta = function(x) x
                         ),
                         theta82 = list(
                             hyperid =  23082,
                             name = "theta82",
                             short.name = "t82",
                             initial = 0,
                             fixed = FALSE,
                             prior = "none",
                             param = numeric(0),
                             to.theta = function(x) x,
                             from.theta = function(x) x
                         ),
                         theta83 = list(
                             hyperid =  23083,
                             name = "theta83",
                             short.name = "t83",
                             initial = 0,
                             fixed = FALSE,
                             prior = "none",
                             param = numeric(0),
                             to.theta = function(x) x,
                             from.theta = function(x) x
                         ),
                         theta84 = list(
                             hyperid =  23084,
                             name = "theta84",
                             short.name = "t84",
                             initial = 0,
                             fixed = FALSE,
                             prior = "none",
                             param = numeric(0),
                             to.theta = function(x) x,
                             from.theta = function(x) x
                         ),
                         theta85 = list(
                             hyperid =  23085,
                             name = "theta85",
                             short.name = "t85",
                             initial = 0,
                             fixed = FALSE,
                             prior = "none",
                             param = numeric(0),
                             to.theta = function(x) x,
                             from.theta = function(x) x
                         ),
                         theta86 = list(
                             hyperid =  23086,
                             name = "theta86",
                             short.name = "t86",
                             initial = 0,
                             fixed = FALSE,
                             prior = "none",
                             param = numeric(0),
                             to.theta = function(x) x,
                             from.theta = function(x) x
                         ),
                         theta87 = list(
                             hyperid =  23087,
                             name = "theta87",
                             short.name = "t87",
                             initial = 0,
                             fixed = FALSE,
                             prior = "none",
                             param = numeric(0),
                             to.theta = function(x) x,
                             from.theta = function(x) x
                         ),
                         theta88 = list(
                             hyperid =  23088,
                             name = "theta88",
                             short.name = "t88",
                             initial = 0,
                             fixed = FALSE,
                             prior = "none",
                             param = numeric(0),
                             to.theta = function(x) x,
                             from.theta = function(x) x
                         ),
                         theta89 = list(
                             hyperid =  23089,
                             name = "theta89",
                             short.name = "t89",
                             initial = 0,
                             fixed = FALSE,
                             prior = "none",
                             param = numeric(0),
                             to.theta = function(x) x,
                             from.theta = function(x) x
                         ),
                         theta90 = list(
                             hyperid =  23090,
                             name = "theta90",
                             short.name = "t90",
                             initial = 0,
                             fixed = FALSE,
                             prior = "none",
                             param = numeric(0),
                             to.theta = function(x) x,
                             from.theta = function(x) x
                         ),
                         theta91 = list(
                             hyperid =  23091,
                             name = "theta91",
                             short.name = "t91",
                             initial = 0,
                             fixed = FALSE,
                             prior = "none",
                             param = numeric(0),
                             to.theta = function(x) x,
                             from.theta = function(x) x
                         ),
                         theta92 = list(
                             hyperid =  23092,
                             name = "theta92",
                             short.name = "t92",
                             initial = 0,
                             fixed = FALSE,
                             prior = "none",
                             param = numeric(0),
                             to.theta = function(x) x,
                             from.theta = function(x) x
                         ),
                         theta93 = list(
                             hyperid =  23093,
                             name = "theta93",
                             short.name = "t93",
                             initial = 0,
                             fixed = FALSE,
                             prior = "none",
                             param = numeric(0),
                             to.theta = function(x) x,
                             from.theta = function(x) x
                         ),
                         theta94 = list(
                             hyperid =  23094,
                             name = "theta94",
                             short.name = "t94",
                             initial = 0,
                             fixed = FALSE,
                             prior = "none",
                             param = numeric(0),
                             to.theta = function(x) x,
                             from.theta = function(x) x
                         ),
                         theta95 = list(
                             hyperid =  23095,
                             name = "theta95",
                             short.name = "t95",
                             initial = 0,
                             fixed = FALSE,
                             prior = "none",
                             param = numeric(0),
                             to.theta = function(x) x,
                             from.theta = function(x) x
                         ),
                         theta96 = list(
                             hyperid =  23096,
                             name = "theta96",
                             short.name = "t96",
                             initial = 0,
                             fixed = FALSE,
                             prior = "none",
                             param = numeric(0),
                             to.theta = function(x) x,
                             from.theta = function(x) x
                         ),
                         theta97 = list(
                             hyperid =  23097,
                             name = "theta97",
                             short.name = "t97",
                             initial = 0,
                             fixed = FALSE,
                             prior = "none",
                             param = numeric(0),
                             to.theta = function(x) x,
                             from.theta = function(x) x
                         ),
                         theta98 = list(
                             hyperid =  23098,
                             name = "theta98",
                             short.name = "t98",
                             initial = 0,
                             fixed = FALSE,
                             prior = "none",
                             param = numeric(0),
                             to.theta = function(x) x,
                             from.theta = function(x) x
                         ),
                         theta99 = list(
                             hyperid =  23099,
                             name = "theta99",
                             short.name = "t99",
                             initial = 0,
                             fixed = FALSE,
                             prior = "none",
                             param = numeric(0),
                             to.theta = function(x) x,
                             from.theta = function(x) x
                         ),
                         theta100 = list(
                             hyperid =  23100,
                             name = "theta100",
                             short.name = "t100",
                             initial = 0,
                             fixed = FALSE,
                             prior = "none",
                             param = numeric(0),
                             to.theta = function(x) x,
                             from.theta = function(x) x
                         )
                     ),
                     constr = FALSE,
                     nrow.ncol = FALSE,
                     augmented = FALSE,
                     aug.factor = 1L,
                     aug.constr = NULL,
                     n.div.by = NULL,
                     n.required = TRUE,
                     set.default.values = TRUE,
                     pdf = "spde2"
                 ),

                 spde3 = list(
                     doc = "A SPDE3 model", 
                     ## to many parameters here, but ...
                     hyper = list(
                         theta1 = list(
                             hyperid =  24001,
                             name = "theta1",
                             short.name = "t1",
                             initial = 0,
                             fixed = FALSE,
                             prior = "mvnorm",
                             param = c(1, 1),
                             to.theta = function(x) x,
                             from.theta = function(x) x
                         ),
                         theta2 = list(
                             hyperid =  24002,
                             name = "theta2",
                             short.name = "t2",
                             initial = 0,
                             fixed = FALSE,
                             prior = "none",
                             param = numeric(0),
                             to.theta = function(x) x,
                             from.theta = function(x) x
                         ),
                         theta3 = list(
                             hyperid =  24003,
                             name = "theta3",
                             short.name = "t3",
                             initial = 0,
                             fixed = FALSE,
                             prior = "none",
                             param = numeric(0),
                             to.theta = function(x) x,
                             from.theta = function(x) x
                         ),
                         theta4 = list(
                             hyperid =  24004,
                             name = "theta4",
                             short.name = "t4",
                             initial = 0,
                             fixed = FALSE,
                             prior = "none",
                             param = numeric(0),
                             to.theta = function(x) x,
                             from.theta = function(x) x
                         ),
                         theta5 = list(
                             hyperid =  24005,
                             name = "theta5",
                             short.name = "t5",
                             initial = 0,
                             fixed = FALSE,
                             prior = "none",
                             param = numeric(0),
                             to.theta = function(x) x,
                             from.theta = function(x) x
                         ),
                         theta6 = list(
                             hyperid =  24006,
                             name = "theta6",
                             short.name = "t6",
                             initial = 0,
                             fixed = FALSE,
                             prior = "none",
                             param = numeric(0),
                             to.theta = function(x) x,
                             from.theta = function(x) x
                         ),
                         theta7 = list(
                             hyperid =  24007,
                             name = "theta7",
                             short.name = "t7",
                             initial = 0,
                             fixed = FALSE,
                             prior = "none",
                             param = numeric(0),
                             to.theta = function(x) x,
                             from.theta = function(x) x
                         ),
                         theta8 = list(
                             hyperid =  24008,
                             name = "theta8",
                             short.name = "t8",
                             initial = 0,
                             fixed = FALSE,
                             prior = "none",
                             param = numeric(0),
                             to.theta = function(x) x,
                             from.theta = function(x) x
                         ),
                         theta9 = list(
                             hyperid =  24009,
                             name = "theta9",
                             short.name = "t9",
                             initial = 0,
                             fixed = FALSE,
                             prior = "none",
                             param = numeric(0),
                             to.theta = function(x) x,
                             from.theta = function(x) x
                         ),
                         theta10 = list(
                             hyperid =  24010,
                             name = "theta10",
                             short.name = "t10",
                             initial = 0,
                             fixed = FALSE,
                             prior = "none",
                             param = numeric(0),
                             to.theta = function(x) x,
                             from.theta = function(x) x
                         ),
                         theta11 = list(
                             hyperid =  24011,
                             name = "theta11",
                             short.name = "t11",
                             initial = 0,
                             fixed = FALSE,
                             prior = "none",
                             param = numeric(0),
                             to.theta = function(x) x,
                             from.theta = function(x) x
                         ),
                         theta12 = list(
                             hyperid =  24012,
                             name = "theta12",
                             short.name = "t12",
                             initial = 0,
                             fixed = FALSE,
                             prior = "none",
                             param = numeric(0),
                             to.theta = function(x) x,
                             from.theta = function(x) x
                         ),
                         theta13 = list(
                             hyperid =  24013,
                             name = "theta13",
                             short.name = "t13",
                             initial = 0,
                             fixed = FALSE,
                             prior = "none",
                             param = numeric(0),
                             to.theta = function(x) x,
                             from.theta = function(x) x
                         ),
                         theta14 = list(
                             hyperid =  24014,
                             name = "theta14",
                             short.name = "t14",
                             initial = 0,
                             fixed = FALSE,
                             prior = "none",
                             param = numeric(0),
                             to.theta = function(x) x,
                             from.theta = function(x) x
                         ),
                         theta15 = list(
                             hyperid =  24015,
                             name = "theta15",
                             short.name = "t15",
                             initial = 0,
                             fixed = FALSE,
                             prior = "none",
                             param = numeric(0),
                             to.theta = function(x) x,
                             from.theta = function(x) x
                         ),
                         theta16 = list(
                             hyperid =  24016,
                             name = "theta16",
                             short.name = "t16",
                             initial = 0,
                             fixed = FALSE,
                             prior = "none",
                             param = numeric(0),
                             to.theta = function(x) x,
                             from.theta = function(x) x
                         ),
                         theta17 = list(
                             hyperid =  24017,
                             name = "theta17",
                             short.name = "t17",
                             initial = 0,
                             fixed = FALSE,
                             prior = "none",
                             param = numeric(0),
                             to.theta = function(x) x,
                             from.theta = function(x) x
                         ),
                         theta18 = list(
                             hyperid =  24018,
                             name = "theta18",
                             short.name = "t18",
                             initial = 0,
                             fixed = FALSE,
                             prior = "none",
                             param = numeric(0),
                             to.theta = function(x) x,
                             from.theta = function(x) x
                         ),
                         theta19 = list(
                             hyperid =  24019,
                             name = "theta19",
                             short.name = "t19",
                             initial = 0,
                             fixed = FALSE,
                             prior = "none",
                             param = numeric(0),
                             to.theta = function(x) x,
                             from.theta = function(x) x
                         ),
                         theta20 = list(
                             hyperid =  24020,
                             name = "theta20",
                             short.name = "t20",
                             initial = 0,
                             fixed = FALSE,
                             prior = "none",
                             param = numeric(0),
                             to.theta = function(x) x,
                             from.theta = function(x) x
                         ),
                         theta21 = list(
                             hyperid =  24021,
                             name = "theta21",
                             short.name = "t21",
                             initial = 0,
                             fixed = FALSE,
                             prior = "none",
                             param = numeric(0),
                             to.theta = function(x) x,
                             from.theta = function(x) x
                         ),
                         theta22 = list(
                             hyperid =  24022,
                             name = "theta22",
                             short.name = "t22",
                             initial = 0,
                             fixed = FALSE,
                             prior = "none",
                             param = numeric(0),
                             to.theta = function(x) x,
                             from.theta = function(x) x
                         ),
                         theta23 = list(
                             hyperid =  24023,
                             name = "theta23",
                             short.name = "t23",
                             initial = 0,
                             fixed = FALSE,
                             prior = "none",
                             param = numeric(0),
                             to.theta = function(x) x,
                             from.theta = function(x) x
                         ),
                         theta24 = list(
                             hyperid =  24024,
                             name = "theta24",
                             short.name = "t24",
                             initial = 0,
                             fixed = FALSE,
                             prior = "none",
                             param = numeric(0),
                             to.theta = function(x) x,
                             from.theta = function(x) x
                         ),
                         theta25 = list(
                             hyperid =  24025,
                             name = "theta25",
                             short.name = "t25",
                             initial = 0,
                             fixed = FALSE,
                             prior = "none",
                             param = numeric(0),
                             to.theta = function(x) x,
                             from.theta = function(x) x
                         ),
                         theta26 = list(
                             hyperid =  24026,
                             name = "theta26",
                             short.name = "t26",
                             initial = 0,
                             fixed = FALSE,
                             prior = "none",
                             param = numeric(0),
                             to.theta = function(x) x,
                             from.theta = function(x) x
                         ),
                         theta27 = list(
                             hyperid =  24027,
                             name = "theta27",
                             short.name = "t27",
                             initial = 0,
                             fixed = FALSE,
                             prior = "none",
                             param = numeric(0),
                             to.theta = function(x) x,
                             from.theta = function(x) x
                         ),
                         theta28 = list(
                             hyperid =  24028,
                             name = "theta28",
                             short.name = "t28",
                             initial = 0,
                             fixed = FALSE,
                             prior = "none",
                             param = numeric(0),
                             to.theta = function(x) x,
                             from.theta = function(x) x
                         ),
                         theta29 = list(
                             hyperid =  24029,
                             name = "theta29",
                             short.name = "t29",
                             initial = 0,
                             fixed = FALSE,
                             prior = "none",
                             param = numeric(0),
                             to.theta = function(x) x,
                             from.theta = function(x) x
                         ),
                         theta30 = list(
                             hyperid =  24030,
                             name = "theta30",
                             short.name = "t30",
                             initial = 0,
                             fixed = FALSE,
                             prior = "none",
                             param = numeric(0),
                             to.theta = function(x) x,
                             from.theta = function(x) x
                         ),
                         theta31 = list(
                             hyperid =  24031,
                             name = "theta31",
                             short.name = "t31",
                             initial = 0,
                             fixed = FALSE,
                             prior = "none",
                             param = numeric(0),
                             to.theta = function(x) x,
                             from.theta = function(x) x
                         ),
                         theta32 = list(
                             hyperid =  24032,
                             name = "theta32",
                             short.name = "t32",
                             initial = 0,
                             fixed = FALSE,
                             prior = "none",
                             param = numeric(0),
                             to.theta = function(x) x,
                             from.theta = function(x) x
                         ),
                         theta33 = list(
                             hyperid =  24033,
                             name = "theta33",
                             short.name = "t33",
                             initial = 0,
                             fixed = FALSE,
                             prior = "none",
                             param = numeric(0),
                             to.theta = function(x) x,
                             from.theta = function(x) x
                         ),
                         theta34 = list(
                             hyperid =  24034,
                             name = "theta34",
                             short.name = "t34",
                             initial = 0,
                             fixed = FALSE,
                             prior = "none",
                             param = numeric(0),
                             to.theta = function(x) x,
                             from.theta = function(x) x
                         ),
                         theta35 = list(
                             hyperid =  24035,
                             name = "theta35",
                             short.name = "t35",
                             initial = 0,
                             fixed = FALSE,
                             prior = "none",
                             param = numeric(0),
                             to.theta = function(x) x,
                             from.theta = function(x) x
                         ),
                         theta36 = list(
                             hyperid =  24036,
                             name = "theta36",
                             short.name = "t36",
                             initial = 0,
                             fixed = FALSE,
                             prior = "none",
                             param = numeric(0),
                             to.theta = function(x) x,
                             from.theta = function(x) x
                         ),
                         theta37 = list(
                             hyperid =  24037,
                             name = "theta37",
                             short.name = "t37",
                             initial = 0,
                             fixed = FALSE,
                             prior = "none",
                             param = numeric(0),
                             to.theta = function(x) x,
                             from.theta = function(x) x
                         ),
                         theta38 = list(
                             hyperid =  24038,
                             name = "theta38",
                             short.name = "t38",
                             initial = 0,
                             fixed = FALSE,
                             prior = "none",
                             param = numeric(0),
                             to.theta = function(x) x,
                             from.theta = function(x) x
                         ),
                         theta39 = list(
                             hyperid =  24039,
                             name = "theta39",
                             short.name = "t39",
                             initial = 0,
                             fixed = FALSE,
                             prior = "none",
                             param = numeric(0),
                             to.theta = function(x) x,
                             from.theta = function(x) x
                         ),
                         theta40 = list(
                             hyperid =  24040,
                             name = "theta40",
                             short.name = "t40",
                             initial = 0,
                             fixed = FALSE,
                             prior = "none",
                             param = numeric(0),
                             to.theta = function(x) x,
                             from.theta = function(x) x
                         ),
                         theta41 = list(
                             hyperid =  24041,
                             name = "theta41",
                             short.name = "t41",
                             initial = 0,
                             fixed = FALSE,
                             prior = "none",
                             param = numeric(0),
                             to.theta = function(x) x,
                             from.theta = function(x) x
                         ),
                         theta42 = list(
                             hyperid =  24042,
                             name = "theta42",
                             short.name = "t42",
                             initial = 0,
                             fixed = FALSE,
                             prior = "none",
                             param = numeric(0),
                             to.theta = function(x) x,
                             from.theta = function(x) x
                         ),
                         theta43 = list(
                             hyperid =  24043,
                             name = "theta43",
                             short.name = "t43",
                             initial = 0,
                             fixed = FALSE,
                             prior = "none",
                             param = numeric(0),
                             to.theta = function(x) x,
                             from.theta = function(x) x
                         ),
                         theta44 = list(
                             hyperid =  24044,
                             name = "theta44",
                             short.name = "t44",
                             initial = 0,
                             fixed = FALSE,
                             prior = "none",
                             param = numeric(0),
                             to.theta = function(x) x,
                             from.theta = function(x) x
                         ),
                         theta45 = list(
                             hyperid =  24045,
                             name = "theta45",
                             short.name = "t45",
                             initial = 0,
                             fixed = FALSE,
                             prior = "none",
                             param = numeric(0),
                             to.theta = function(x) x,
                             from.theta = function(x) x
                         ),
                         theta46 = list(
                             hyperid =  24046,
                             name = "theta46",
                             short.name = "t46",
                             initial = 0,
                             fixed = FALSE,
                             prior = "none",
                             param = numeric(0),
                             to.theta = function(x) x,
                             from.theta = function(x) x
                         ),
                         theta47 = list(
                             hyperid =  24047,
                             name = "theta47",
                             short.name = "t47",
                             initial = 0,
                             fixed = FALSE,
                             prior = "none",
                             param = numeric(0),
                             to.theta = function(x) x,
                             from.theta = function(x) x
                         ),
                         theta48 = list(
                             hyperid =  24048,
                             name = "theta48",
                             short.name = "t48",
                             initial = 0,
                             fixed = FALSE,
                             prior = "none",
                             param = numeric(0),
                             to.theta = function(x) x,
                             from.theta = function(x) x
                         ),
                         theta49 = list(
                             hyperid =  24049,
                             name = "theta49",
                             short.name = "t49",
                             initial = 0,
                             fixed = FALSE,
                             prior = "none",
                             param = numeric(0),
                             to.theta = function(x) x,
                             from.theta = function(x) x
                         ),
                         theta50 = list(
                             hyperid =  24050,
                             name = "theta50",
                             short.name = "t50",
                             initial = 0,
                             fixed = FALSE,
                             prior = "none",
                             param = numeric(0),
                             to.theta = function(x) x,
                             from.theta = function(x) x
                         ),
                         theta51 = list(
                             hyperid =  24051,
                             name = "theta51",
                             short.name = "t51",
                             initial = 0,
                             fixed = FALSE,
                             prior = "none",
                             param = numeric(0),
                             to.theta = function(x) x,
                             from.theta = function(x) x
                         ),
                         theta52 = list(
                             hyperid =  24052,
                             name = "theta52",
                             short.name = "t52",
                             initial = 0,
                             fixed = FALSE,
                             prior = "none",
                             param = numeric(0),
                             to.theta = function(x) x,
                             from.theta = function(x) x
                         ),
                         theta53 = list(
                             hyperid =  24053,
                             name = "theta53",
                             short.name = "t53",
                             initial = 0,
                             fixed = FALSE,
                             prior = "none",
                             param = numeric(0),
                             to.theta = function(x) x,
                             from.theta = function(x) x
                         ),
                         theta54 = list(
                             hyperid =  24054,
                             name = "theta54",
                             short.name = "t54",
                             initial = 0,
                             fixed = FALSE,
                             prior = "none",
                             param = numeric(0),
                             to.theta = function(x) x,
                             from.theta = function(x) x
                         ),
                         theta55 = list(
                             hyperid =  24055,
                             name = "theta55",
                             short.name = "t55",
                             initial = 0,
                             fixed = FALSE,
                             prior = "none",
                             param = numeric(0),
                             to.theta = function(x) x,
                             from.theta = function(x) x
                         ),
                         theta56 = list(
                             hyperid =  24056,
                             name = "theta56",
                             short.name = "t56",
                             initial = 0,
                             fixed = FALSE,
                             prior = "none",
                             param = numeric(0),
                             to.theta = function(x) x,
                             from.theta = function(x) x
                         ),
                         theta57 = list(
                             hyperid =  24057,
                             name = "theta57",
                             short.name = "t57",
                             initial = 0,
                             fixed = FALSE,
                             prior = "none",
                             param = numeric(0),
                             to.theta = function(x) x,
                             from.theta = function(x) x
                         ),
                         theta58 = list(
                             hyperid =  24058,
                             name = "theta58",
                             short.name = "t58",
                             initial = 0,
                             fixed = FALSE,
                             prior = "none",
                             param = numeric(0),
                             to.theta = function(x) x,
                             from.theta = function(x) x
                         ),
                         theta59 = list(
                             hyperid =  24059,
                             name = "theta59",
                             short.name = "t59",
                             initial = 0,
                             fixed = FALSE,
                             prior = "none",
                             param = numeric(0),
                             to.theta = function(x) x,
                             from.theta = function(x) x
                         ),
                         theta60 = list(
                             hyperid =  24060,
                             name = "theta60",
                             short.name = "t60",
                             initial = 0,
                             fixed = FALSE,
                             prior = "none",
                             param = numeric(0),
                             to.theta = function(x) x,
                             from.theta = function(x) x
                         ),
                         theta61 = list(
                             hyperid =  24061,
                             name = "theta61",
                             short.name = "t61",
                             initial = 0,
                             fixed = FALSE,
                             prior = "none",
                             param = numeric(0),
                             to.theta = function(x) x,
                             from.theta = function(x) x
                         ),
                         theta62 = list(
                             hyperid =  24062,
                             name = "theta62",
                             short.name = "t62",
                             initial = 0,
                             fixed = FALSE,
                             prior = "none",
                             param = numeric(0),
                             to.theta = function(x) x,
                             from.theta = function(x) x
                         ),
                         theta63 = list(
                             hyperid =  24063,
                             name = "theta63",
                             short.name = "t63",
                             initial = 0,
                             fixed = FALSE,
                             prior = "none",
                             param = numeric(0),
                             to.theta = function(x) x,
                             from.theta = function(x) x
                         ),
                         theta64 = list(
                             hyperid =  24064,
                             name = "theta64",
                             short.name = "t64",
                             initial = 0,
                             fixed = FALSE,
                             prior = "none",
                             param = numeric(0),
                             to.theta = function(x) x,
                             from.theta = function(x) x
                         ),
                         theta65 = list(
                             hyperid =  24065,
                             name = "theta65",
                             short.name = "t65",
                             initial = 0,
                             fixed = FALSE,
                             prior = "none",
                             param = numeric(0),
                             to.theta = function(x) x,
                             from.theta = function(x) x
                         ),
                         theta66 = list(
                             hyperid =  24066,
                             name = "theta66",
                             short.name = "t66",
                             initial = 0,
                             fixed = FALSE,
                             prior = "none",
                             param = numeric(0),
                             to.theta = function(x) x,
                             from.theta = function(x) x
                         ),
                         theta67 = list(
                             hyperid =  24067,
                             name = "theta67",
                             short.name = "t67",
                             initial = 0,
                             fixed = FALSE,
                             prior = "none",
                             param = numeric(0),
                             to.theta = function(x) x,
                             from.theta = function(x) x
                         ),
                         theta68 = list(
                             hyperid =  24068,
                             name = "theta68",
                             short.name = "t68",
                             initial = 0,
                             fixed = FALSE,
                             prior = "none",
                             param = numeric(0),
                             to.theta = function(x) x,
                             from.theta = function(x) x
                         ),
                         theta69 = list(
                             hyperid =  24069,
                             name = "theta69",
                             short.name = "t69",
                             initial = 0,
                             fixed = FALSE,
                             prior = "none",
                             param = numeric(0),
                             to.theta = function(x) x,
                             from.theta = function(x) x
                         ),
                         theta70 = list(
                             hyperid =  24070,
                             name = "theta70",
                             short.name = "t70",
                             initial = 0,
                             fixed = FALSE,
                             prior = "none",
                             param = numeric(0),
                             to.theta = function(x) x,
                             from.theta = function(x) x
                         ),
                         theta71 = list(
                             hyperid =  24071,
                             name = "theta71",
                             short.name = "t71",
                             initial = 0,
                             fixed = FALSE,
                             prior = "none",
                             param = numeric(0),
                             to.theta = function(x) x,
                             from.theta = function(x) x
                         ),
                         theta72 = list(
                             hyperid =  24072,
                             name = "theta72",
                             short.name = "t72",
                             initial = 0,
                             fixed = FALSE,
                             prior = "none",
                             param = numeric(0),
                             to.theta = function(x) x,
                             from.theta = function(x) x
                         ),
                         theta73 = list(
                             hyperid =  24073,
                             name = "theta73",
                             short.name = "t73",
                             initial = 0,
                             fixed = FALSE,
                             prior = "none",
                             param = numeric(0),
                             to.theta = function(x) x,
                             from.theta = function(x) x
                         ),
                         theta74 = list(
                             hyperid =  24074,
                             name = "theta74",
                             short.name = "t74",
                             initial = 0,
                             fixed = FALSE,
                             prior = "none",
                             param = numeric(0),
                             to.theta = function(x) x,
                             from.theta = function(x) x
                         ),
                         theta75 = list(
                             hyperid =  24075,
                             name = "theta75",
                             short.name = "t75",
                             initial = 0,
                             fixed = FALSE,
                             prior = "none",
                             param = numeric(0),
                             to.theta = function(x) x,
                             from.theta = function(x) x
                         ),
                         theta76 = list(
                             hyperid =  24076,
                             name = "theta76",
                             short.name = "t76",
                             initial = 0,
                             fixed = FALSE,
                             prior = "none",
                             param = numeric(0),
                             to.theta = function(x) x,
                             from.theta = function(x) x
                         ),
                         theta77 = list(
                             hyperid =  24077,
                             name = "theta77",
                             short.name = "t77",
                             initial = 0,
                             fixed = FALSE,
                             prior = "none",
                             param = numeric(0),
                             to.theta = function(x) x,
                             from.theta = function(x) x
                         ),
                         theta78 = list(
                             hyperid =  24078,
                             name = "theta78",
                             short.name = "t78",
                             initial = 0,
                             fixed = FALSE,
                             prior = "none",
                             param = numeric(0),
                             to.theta = function(x) x,
                             from.theta = function(x) x
                         ),
                         theta79 = list(
                             hyperid =  24079,
                             name = "theta79",
                             short.name = "t79",
                             initial = 0,
                             fixed = FALSE,
                             prior = "none",
                             param = numeric(0),
                             to.theta = function(x) x,
                             from.theta = function(x) x
                         ),
                         theta80 = list(
                             hyperid =  24080,
                             name = "theta80",
                             short.name = "t80",
                             initial = 0,
                             fixed = FALSE,
                             prior = "none",
                             param = numeric(0),
                             to.theta = function(x) x,
                             from.theta = function(x) x
                         ),
                         theta81 = list(
                             hyperid =  24081,
                             name = "theta81",
                             short.name = "t81",
                             initial = 0,
                             fixed = FALSE,
                             prior = "none",
                             param = numeric(0),
                             to.theta = function(x) x,
                             from.theta = function(x) x
                         ),
                         theta82 = list(
                             hyperid =  24082,
                             name = "theta82",
                             short.name = "t82",
                             initial = 0,
                             fixed = FALSE,
                             prior = "none",
                             param = numeric(0),
                             to.theta = function(x) x,
                             from.theta = function(x) x
                         ),
                         theta83 = list(
                             hyperid =  24083,
                             name = "theta83",
                             short.name = "t83",
                             initial = 0,
                             fixed = FALSE,
                             prior = "none",
                             param = numeric(0),
                             to.theta = function(x) x,
                             from.theta = function(x) x
                         ),
                         theta84 = list(
                             hyperid =  24084,
                             name = "theta84",
                             short.name = "t84",
                             initial = 0,
                             fixed = FALSE,
                             prior = "none",
                             param = numeric(0),
                             to.theta = function(x) x,
                             from.theta = function(x) x
                         ),
                         theta85 = list(
                             hyperid =  24085,
                             name = "theta85",
                             short.name = "t85",
                             initial = 0,
                             fixed = FALSE,
                             prior = "none",
                             param = numeric(0),
                             to.theta = function(x) x,
                             from.theta = function(x) x
                         ),
                         theta86 = list(
                             hyperid =  24086,
                             name = "theta86",
                             short.name = "t86",
                             initial = 0,
                             fixed = FALSE,
                             prior = "none",
                             param = numeric(0),
                             to.theta = function(x) x,
                             from.theta = function(x) x
                         ),
                         theta87 = list(
                             hyperid =  24087,
                             name = "theta87",
                             short.name = "t87",
                             initial = 0,
                             fixed = FALSE,
                             prior = "none",
                             param = numeric(0),
                             to.theta = function(x) x,
                             from.theta = function(x) x
                         ),
                         theta88 = list(
                             hyperid =  24088,
                             name = "theta88",
                             short.name = "t88",
                             initial = 0,
                             fixed = FALSE,
                             prior = "none",
                             param = numeric(0),
                             to.theta = function(x) x,
                             from.theta = function(x) x
                         ),
                         theta89 = list(
                             hyperid =  24089,
                             name = "theta89",
                             short.name = "t89",
                             initial = 0,
                             fixed = FALSE,
                             prior = "none",
                             param = numeric(0),
                             to.theta = function(x) x,
                             from.theta = function(x) x
                         ),
                         theta90 = list(
                             hyperid =  24090,
                             name = "theta90",
                             short.name = "t90",
                             initial = 0,
                             fixed = FALSE,
                             prior = "none",
                             param = numeric(0),
                             to.theta = function(x) x,
                             from.theta = function(x) x
                         ),
                         theta91 = list(
                             hyperid =  24091,
                             name = "theta91",
                             short.name = "t91",
                             initial = 0,
                             fixed = FALSE,
                             prior = "none",
                             param = numeric(0),
                             to.theta = function(x) x,
                             from.theta = function(x) x
                         ),
                         theta92 = list(
                             hyperid =  24092,
                             name = "theta92",
                             short.name = "t92",
                             initial = 0,
                             fixed = FALSE,
                             prior = "none",
                             param = numeric(0),
                             to.theta = function(x) x,
                             from.theta = function(x) x
                         ),
                         theta93 = list(
                             hyperid =  24093,
                             name = "theta93",
                             short.name = "t93",
                             initial = 0,
                             fixed = FALSE,
                             prior = "none",
                             param = numeric(0),
                             to.theta = function(x) x,
                             from.theta = function(x) x
                         ),
                         theta94 = list(
                             hyperid =  24094,
                             name = "theta94",
                             short.name = "t94",
                             initial = 0,
                             fixed = FALSE,
                             prior = "none",
                             param = numeric(0),
                             to.theta = function(x) x,
                             from.theta = function(x) x
                         ),
                         theta95 = list(
                             hyperid =  24095,
                             name = "theta95",
                             short.name = "t95",
                             initial = 0,
                             fixed = FALSE,
                             prior = "none",
                             param = numeric(0),
                             to.theta = function(x) x,
                             from.theta = function(x) x
                         ),
                         theta96 = list(
                             hyperid =  24096,
                             name = "theta96",
                             short.name = "t96",
                             initial = 0,
                             fixed = FALSE,
                             prior = "none",
                             param = numeric(0),
                             to.theta = function(x) x,
                             from.theta = function(x) x
                         ),
                         theta97 = list(
                             hyperid =  24097,
                             name = "theta97",
                             short.name = "t97",
                             initial = 0,
                             fixed = FALSE,
                             prior = "none",
                             param = numeric(0),
                             to.theta = function(x) x,
                             from.theta = function(x) x
                         ),
                         theta98 = list(
                             hyperid =  24098,
                             name = "theta98",
                             short.name = "t98",
                             initial = 0,
                             fixed = FALSE,
                             prior = "none",
                             param = numeric(0),
                             to.theta = function(x) x,
                             from.theta = function(x) x
                         ),
                         theta99 = list(
                             hyperid =  24099,
                             name = "theta99",
                             short.name = "t99",
                             initial = 0,
                             fixed = FALSE,
                             prior = "none",
                             param = numeric(0),
                             to.theta = function(x) x,
                             from.theta = function(x) x
                         ),
                         theta100 = list(
                             hyperid =  24100,
                             name = "theta100",
                             short.name = "t100",
                             initial = 0,
                             fixed = FALSE,
                             prior = "none",
                             param = numeric(0),
                             to.theta = function(x) x,
                             from.theta = function(x) x
                         )
                     ),
                     constr = FALSE,
                     nrow.ncol = FALSE,
                     augmented = FALSE,
                     aug.factor = 1L,
                     aug.constr = NULL,
                     n.div.by = NULL,
                     n.required = TRUE,
                     set.default.values = TRUE,
                     pdf = "spde3"
                 ),

                 iid1d = list(
                     doc = "Gaussian random effect in dim=1 with Wishart prior", 
                     hyper = list(
                         theta = list(
                             hyperid =  25001,
                             name = "precision",
                             short.name = "prec",
                             initial = 4,
                             fixed = FALSE,
                             prior = "wishart1d",
                             ## this is the
                             ## corresponding
                             ## wishart prior for
                             ## the (a, b) gamma
                             ## parameters, so
                             ## default we get the
                             ## same.
                             param = c(2*1, 2*0.00005),
                             to.theta = function(x) log(x),
                             from.theta = function(x) exp(x)
                         )
                     ),
                     constr = FALSE,
                     nrow.ncol = FALSE,
                     augmented = FALSE,
                     aug.factor = 1L,
                     aug.constr = NULL,
                     n.div.by = NULL,
                     n.required = FALSE, ## TRUE is not
                     ## needed here.
                     set.default.values = TRUE,
                     pdf = "iid123d"
                 ),

                 iid2d = list(
                     doc = "Gaussian random effect in dim=2 with Wishart prior", 
                     hyper = list(
                         theta1 = list(
                             hyperid =  26001,
                             name = "log precision1",
                             short.name = "prec1",
                             initial = 4,
                             fixed = FALSE,
                             prior = "wishart2d",
                             param = c(4, 1, 1, 0),
                             to.theta = function(x) log(x),
                             from.theta = function(x) exp(x)
                         ),
                         theta2 = list(
                             hyperid =  26002,
                             name = "log precision2",
                             short.name = "prec2",
                             initial = 4,
                             fixed = FALSE,
                             prior = "none",
                             param = numeric(),
                             to.theta = function(x) log(x),
                             from.theta = function(x) exp(x)
                         ),
                         theta3 = list(
                             hyperid =  26003,
                             name = "logit correlation",
                             short.name = "cor",
                             initial = 4,
                             fixed = FALSE,
                             prior = "none",
                             param = numeric(),
                             to.theta = function(x) log((1+x)/(1-x)),
                             from.theta = function(x) 2*exp(x)/(1+exp(x))-1
                         )
                     ),
                     constr = FALSE,
                     nrow.ncol = FALSE,
                     augmented = TRUE,
                     aug.factor = 1L,
                     aug.constr = c(1L, 2L),
                     n.div.by = 2L,
                     n.required = TRUE,
                     set.default.values = TRUE,
                     pdf = "iid123d"
                 ),

                 iid3d = list(
                     doc = "Gaussian random effect in dim=3 with Wishart prior", 
                     hyper = list(
                         theta1 = list(
                             hyperid =  27001,
                             name = "log precision1",
                             short.name = "prec1",
                             initial = 4,
                             fixed = FALSE,
                             prior = "wishart3d",
                             param = c(7, 1, 1, 1, 0, 0, 0),
                             to.theta = function(x) log(x),
                             from.theta = function(x) exp(x)
                         ),
                         theta2 = list(
                             hyperid =  27002,
                             name = "log precision2",
                             short.name = "prec2",
                             initial = 4,
                             fixed = FALSE,
                             prior = "none",
                             param = numeric(),
                             to.theta = function(x) log(x),
                             from.theta = function(x) exp(x)
                         ),
                         theta3 = list(
                             hyperid =  27003,
                             name = "log precision3",
                             short.name = "prec3",
                             initial = 4,
                             fixed = FALSE,
                             prior = "none",
                             param = numeric(),
                             to.theta = function(x) log(x),
                             from.theta = function(x) exp(x)
                         ),
                         theta4= list(
                             hyperid =  27004,
                             name = "logit correlation12",
                             short.name = "cor12",
                             initial = 0,
                             fixed = FALSE,
                             prior = "none",
                             param = numeric(),
                             to.theta = function(x) log((1+x)/(1-x)),
                             from.theta = function(x) 2*exp(x)/(1+exp(x))-1
                         ),
                         theta5 = list(
                             hyperid =  27005,
                             name = "logit correlation13",
                             short.name = "cor13",
                             initial = 0,
                             fixed = FALSE,
                             prior = "none",
                             param = numeric(),
                             to.theta = function(x) log((1+x)/(1-x)),
                             from.theta = function(x) 2*exp(x)/(1+exp(x))-1
                         ),
                         theta6 = list(
                             hyperid =  27006,
                             name = "logit correlation23",
                             short.name = "cor23",
                             initial = 0,
                             fixed = FALSE,
                             prior = "none",
                             param = numeric(),
                             to.theta = function(x) log((1+x)/(1-x)),
                             from.theta = function(x) 2*exp(x)/(1+exp(x))-1
                         )
                     ),
                     constr = FALSE,
                     nrow.ncol = FALSE,
                     augmented = TRUE,
                     aug.factor = 1L,
                     aug.constr = c(1L, 2L, 3L),
                     n.div.by = 3L,
                     n.required = TRUE,
                     set.default.values = TRUE,
                     pdf = "iid123d"
                 ),

                 iid4d = list(
                     doc = "Gaussian random effect in dim=4 with Wishart prior", 
                     hyper = list(
                         theta1 = list(
                             hyperid =  28001,
                             name = "log precision1",
                             short.name = "prec1",
                             initial = 4,
                             fixed = FALSE,
                             prior = "wishart4d",
                             param = c(11, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0),
                             to.theta = function(x) log(x),
                             from.theta = function(x) exp(x)
                         ),
                         theta2 = list(
                             hyperid =  28002,
                             name = "log precision2",
                             short.name = "prec2",
                             initial = 4,
                             fixed = FALSE,
                             prior = "none",
                             param = numeric(),
                             to.theta = function(x) log(x),
                             from.theta = function(x) exp(x)
                         ),
                         theta3 = list(
                             hyperid =  28003,
                             name = "log precision3",
                             short.name = "prec3",
                             initial = 4,
                             fixed = FALSE,
                             prior = "none",
                             param = numeric(),
                             to.theta = function(x) log(x),
                             from.theta = function(x) exp(x)
                         ),
                         theta4 = list(
                             hyperid =  28004,
                             name = "log precision4",
                             short.name = "prec4",
                             initial = 4,
                             fixed = FALSE,
                             prior = "none",
                             param = numeric(),
                             to.theta = function(x) log(x),
                             from.theta = function(x) exp(x)
                         ),
                         theta5= list(
                             hyperid =  28005,
                             name = "logit correlation12",
                             short.name = "cor12",
                             initial = 0,
                             fixed = FALSE,
                             prior = "none",
                             param = numeric(),
                             to.theta = function(x) log((1+x)/(1-x)),
                             from.theta = function(x) 2*exp(x)/(1+exp(x))-1
                         ),
                         theta6 = list(
                             hyperid =  28006,
                             name = "logit correlation13",
                             short.name = "cor13",
                             initial = 0,
                             fixed = FALSE,
                             prior = "none",
                             param = numeric(),
                             to.theta = function(x) log((1+x)/(1-x)),
                             from.theta = function(x) 2*exp(x)/(1+exp(x))-1
                         ),
                         theta7 = list(
                             hyperid =  28007,
                             name = "logit correlation14",
                             short.name = "cor14",
                             initial = 0,
                             fixed = FALSE,
                             prior = "none",
                             param = numeric(),
                             to.theta = function(x) log((1+x)/(1-x)),
                             from.theta = function(x) 2*exp(x)/(1+exp(x))-1
                         ),
                         theta8 = list(
                             hyperid =  28008,
                             name = "logit correlation23",
                             short.name = "cor23",
                             initial = 0,
                             fixed = FALSE,
                             prior = "none",
                             param = numeric(),
                             to.theta = function(x) log((1+x)/(1-x)),
                             from.theta = function(x) 2*exp(x)/(1+exp(x))-1
                         ),
                         theta9 = list(
                             hyperid =  28009,
                             name = "logit correlation24",
                             short.name = "cor24",
                             initial = 0,
                             fixed = FALSE,
                             prior = "none",
                             param = numeric(),
                             to.theta = function(x) log((1+x)/(1-x)),
                             from.theta = function(x) 2*exp(x)/(1+exp(x))-1
                         ),
                         theta10 = list(
                             hyperid =  28010,
                             name = "logit correlation34",
                             short.name = "cor34",
                             initial = 0,
                             fixed = FALSE,
                             prior = "none",
                             param = numeric(),
                             to.theta = function(x) log((1+x)/(1-x)),
                             from.theta = function(x) 2*exp(x)/(1+exp(x))-1
                         )
                     ),
                     constr = FALSE,
                     nrow.ncol = FALSE,
                     augmented = TRUE,
                     aug.factor = 1L,
                     aug.constr = c(1L, 2L, 3L, 4L),
                     n.div.by = 4L,
                     n.required = TRUE,
                     set.default.values = TRUE,
                     pdf = "iid123d"
                 ),

                 iid5d = list(
                     doc = "Gaussian random effect in dim=5 with Wishart prior", 
                     hyper = list(
                         theta1 = list(
                             hyperid =  29001,
                             name = "log precision1",
                             short.name = "prec1",
                             initial = 4,
                             fixed = FALSE,
                             prior = "wishart5d",
                             param = c(16, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
                             to.theta = function(x) log(x),
                             from.theta = function(x) exp(x)
                         ),
                         theta2 = list(
                             hyperid =  29002,
                             name = "log precision2",
                             short.name = "prec2",
                             initial = 4,
                             fixed = FALSE,
                             prior = "none",
                             param = numeric(),
                             to.theta = function(x) log(x),
                             from.theta = function(x) exp(x)
                         ),
                         theta3 = list(
                             hyperid =  29003,
                             name = "log precision3",
                             short.name = "prec3",
                             initial = 4,
                             fixed = FALSE,
                             prior = "none",
                             param = numeric(),
                             to.theta = function(x) log(x),
                             from.theta = function(x) exp(x)
                         ),
                         theta4 = list(
                             hyperid =  29004,
                             name = "log precision4",
                             short.name = "prec4",
                             initial = 4,
                             fixed = FALSE,
                             prior = "none",
                             param = numeric(),
                             to.theta = function(x) log(x),
                             from.theta = function(x) exp(x)
                         ),
                         theta5 = list(
                             hyperid =  29005,
                             name = "log precision5",
                             short.name = "prec5",
                             initial = 4,
                             fixed = FALSE,
                             prior = "none",
                             param = numeric(),
                             to.theta = function(x) log(x),
                             from.theta = function(x) exp(x)
                         ),
                         theta6= list(
                             hyperid =  29006,
                             name = "logit correlation12",
                             short.name = "cor12",
                             initial = 0,
                             fixed = FALSE,
                             prior = "none",
                             param = numeric(),
                             to.theta = function(x) log((1+x)/(1-x)),
                             from.theta = function(x) 2*exp(x)/(1+exp(x))-1
                         ),
                         theta7 = list(
                             hyperid =  29007,
                             name = "logit correlation13",
                             short.name = "cor13",
                             initial = 0,
                             fixed = FALSE,
                             prior = "none",
                             param = numeric(),
                             to.theta = function(x) log((1+x)/(1-x)),
                             from.theta = function(x) 2*exp(x)/(1+exp(x))-1
                         ),
                         theta8 = list(
                             hyperid =  29008,
                             name = "logit correlation14",
                             short.name = "cor14",
                             initial = 0,
                             fixed = FALSE,
                             prior = "none",
                             param = numeric(),
                             to.theta = function(x) log((1+x)/(1-x)),
                             from.theta = function(x) 2*exp(x)/(1+exp(x))-1
                         ),
                         theta9 = list(
                             hyperid =  29009,
                             name = "logit correlation15",
                             short.name = "cor15",
                             initial = 0,
                             fixed = FALSE,
                             prior = "none",
                             param = numeric(),
                             to.theta = function(x) log((1+x)/(1-x)),
                             from.theta = function(x) 2*exp(x)/(1+exp(x))-1
                         ),
                         theta10 = list(
                             hyperid =  29010,
                             name = "logit correlation23",
                             short.name = "cor23",
                             initial = 0,
                             fixed = FALSE,
                             prior = "none",
                             param = numeric(),
                             to.theta = function(x) log((1+x)/(1-x)),
                             from.theta = function(x) 2*exp(x)/(1+exp(x))-1
                         ),
                         theta11 = list(
                             hyperid =  29011,
                             name = "logit correlation24",
                             short.name = "cor24",
                             initial = 0,
                             fixed = FALSE,
                             prior = "none",
                             param = numeric(),
                             to.theta = function(x) log((1+x)/(1-x)),
                             from.theta = function(x) 2*exp(x)/(1+exp(x))-1
                         ),
                         theta12 = list(
                             hyperid =  29012,
                             name = "logit correlation25",
                             short.name = "cor25",
                             initial = 0,
                             fixed = FALSE,
                             prior = "none",
                             param = numeric(),
                             to.theta = function(x) log((1+x)/(1-x)),
                             from.theta = function(x) 2*exp(x)/(1+exp(x))-1
                         ),
                         theta13 = list(
                             hyperid =  29013,
                             name = "logit correlation34",
                             short.name = "cor34",
                             initial = 0,
                             fixed = FALSE,
                             prior = "none",
                             param = numeric(),
                             to.theta = function(x) log((1+x)/(1-x)),
                             from.theta = function(x) 2*exp(x)/(1+exp(x))-1
                         ),
                         theta14 = list(
                             hyperid =  29014,
                             name = "logit correlation35",
                             short.name = "cor35",
                             initial = 0,
                             fixed = FALSE,
                             prior = "none",
                             param = numeric(),
                             to.theta = function(x) log((1+x)/(1-x)),
                             from.theta = function(x) 2*exp(x)/(1+exp(x))-1
                         ),
                         theta15 = list(
                             hyperid =  29015,
                             name = "logit correlation45",
                             short.name = "cor45",
                             initial = 0,
                             fixed = FALSE,
                             prior = "none",
                             param = numeric(),
                             to.theta = function(x) log((1+x)/(1-x)),
                             from.theta = function(x) 2*exp(x)/(1+exp(x))-1
                         )
                     ),
                     constr = FALSE,
                     nrow.ncol = FALSE,
                     augmented = TRUE,
                     aug.factor = 1L,
                     aug.constr = c(1L, 2L, 3L, 4L, 5L),
                     n.div.by = 5L,
                     n.required = TRUE,
                     set.default.values = TRUE,
                     pdf = "iid123d"
                 ),

                 "2diid" = list(
                     doc = "(This model is obsolute)", 
                     hyper = list(
                         theta1 = list(
                             hyperid =  30001,
                             name = "log precision1",
                             short.name = "prec1",
                             initial = 4,
                             fixed = FALSE,
                             prior = "loggamma",
                             param = c(1, 0.00005),
                             to.theta = function(x) log(x),
                             from.theta = function(x) exp(x)
                         ),
                         theta2 = list(
                             hyperid =  30002,
                             name = "log precision2",
                             short.name = "prec2",
                             initial = 4,
                             fixed = FALSE,
                             prior = "loggamma",
                             param = c(1, 0.00005),
                             to.theta = function(x) log(x),
                             from.theta = function(x) exp(x)
                         ),
                         theta3 = list(
                             hyperid =  30003,
                             name = "correlation",
                             short.name = "cor",
                             initial = 4,
                             fixed = FALSE,
                             prior = "normal",
                             param = c(0, 0.15),
                             to.theta = function(x) log((1+x)/(1-x)),
                             from.theta = function(x) 2*exp(x)/(1+exp(x))-1
                         )
                     ),
                     constr = FALSE,
                     nrow.ncol = FALSE,
                     augmented = FALSE,
                     aug.factor = 1L,
                     aug.constr = c(1L, 2L),
                     n.div.by = 2L,
                     n.required = TRUE,
                     set.default.values = TRUE,
                     pdf = "iid123d"
                 ),

                 z = list(
                     doc = "The z-model in a classical mixed model formulation", 
                     hyper = list(
                         theta = list(
                             hyperid =  31001,
                             name = "log precision",
                             short.name = "prec",
                             initial = 4,
                             fixed = FALSE,
                             prior = "loggamma",
                             param = c(1, 0.00005),
                             to.theta = function(x) log(x),
                             from.theta = function(x) exp(x)
                         )
                     ),
                     constr = FALSE,
                     nrow.ncol = FALSE,
                     augmented = FALSE,
                     aug.factor = 1L,
                     aug.constr = NULL,
                     n.div.by = NULL,
                     n.required = TRUE,
                     set.default.values = TRUE,
                     pdf = "z", 
                     status = "experimental"
                 ),

                 rw2d = list(
                     doc = "Thin-plate spline model", 
                     hyper = list(
                         theta = list(
                             hyperid =  32001,
                             name = "log precision",
                             short.name = "prec",
                             initial = 4,
                             fixed = FALSE,
                             prior = "loggamma",
                             param = c(1, 0.00005),
                             to.theta = function(x) log(x),
                             from.theta = function(x) exp(x)
                         )
                     ),
                     constr = TRUE,
                     nrow.ncol = TRUE,
                     augmented = FALSE,
                     aug.factor = 1L,
                     aug.constr = NULL,
                     n.div.by = NULL,
                     n.required = FALSE,
                     set.default.values = TRUE,
                     pdf = "rw2d"
                 ),

                 rw2diid = list(
                     doc = "Thin-plate spline with iid noise", 
                     hyper = list(
                         theta1 = list(
                             hyperid =  33001,
                             name = "log precision",
                             short.name = "prec",
                             prior = "pc.prec",
                             param = c(1, .01),
                             initial = 4,
                             fixed = FALSE,
                             to.theta = function(x) log(x),
                             from.theta = function(x) exp(x)
                         ),
                         theta2 = list(
                             hyperid =  33002,
                             name = "logit phi",
                             short.name = "phi",
                             prior = "pc",
                             param = c(0.5, 0.5),
                             initial = 3,
                             fixed = FALSE,
                             to.theta = function(x) log(x/(1-x)),
                             from.theta = function(x) exp(x)/(1+exp(x))
                         )
                     ),
                     constr = TRUE,
                     nrow.ncol = TRUE,
                     augmented = TRUE,
                     aug.factor = 2L,
                     aug.constr = 2L,
                     n.div.by = NULL,
                     n.required = FALSE,
                     set.default.values = TRUE,
                     status = "experimental", 
                     pdf = "rw2diid"
                 ),

                 slm = list(
                     doc = "Spatial lag model", 
                     hyper = list(
                         theta1 = list(
                             hyperid =  34001,
                             name = "log precision",
                             short.name = "prec",
                             initial = 4,
                             fixed = FALSE,
                             prior = "loggamma",
                             param = c(1, 0.00005),
                             to.theta = function(x) log(x),
                             from.theta = function(x) exp(x)
                         ), 
                         theta2 = list(
                             hyperid =  34002,
                             name = "rho",
                             short.name = "rho",
                             initial = 0,
                             fixed = FALSE,
                             prior = "normal",
                             param = c(0, 10),
                             to.theta = function(x) log(x/(1-x)), 
                             from.theta = function(x) 1/(1+exp(-x))
                         )
                     ),
                     constr = FALSE,
                     nrow.ncol = FALSE,
                     augmented = FALSE,
                     aug.factor = 1L,
                     aug.constr = NULL,
                     n.div.by = NULL,
                     n.required = TRUE,
                     set.default.values = TRUE,
                     pdf = "slm", 
                     status = "experimental"
                 ),

                 matern2d = list(
                     doc = "Matern covariance function on a regular grid", 
                     hyper = list(
                         theta1 = list(
                             hyperid =  35001,
                             name = "log precision",
                             short.name = "prec",
                             initial = 4,
                             fixed = FALSE,
                             prior = "loggamma",
                             param = c(1, 0.00005),
                             to.theta = function(x) log(x),
                             from.theta = function(x) exp(x)
                         ),
                         theta2 = list(
                             hyperid =  35002,
                             name = "log range",
                             short.name = "range",
                             initial = 2,
                             fixed = FALSE,
                             prior = "loggamma",
                             param = c(1, 0.01),
                             to.theta = function(x) log(x),
                             from.theta = function(x) exp(x)
                         )
                     ),
                     constr = FALSE,
                     nrow.ncol = TRUE,
                     augmented = FALSE,
                     aug.factor = 1L,
                     aug.constr = NULL,
                     n.div.by = NULL,
                     n.required = FALSE,
                     set.default.values = TRUE,
                     pdf = "matern2d"
                 ),

                 dmatern = list(
                     doc = "Dense Matern field", 
                     hyper = list(
                         theta1 = list(
                             hyperid =  35101,
                             name = "log precision",
                             short.name = "prec",
                             initial = 3,
                             fixed = FALSE,
                             prior = "pc.prec",
                             param = c(1, 0.01),
                             to.theta = function(x) log(x),
                             from.theta = function(x) exp(x)
                         ),
                         theta2 = list(
                             hyperid =  35102,
                             name = "log range",
                             short.name = "range",
                             initial = 0,
                             fixed = FALSE,
                             prior = "pc.range",
                             param = c(1, 0.5),
                             to.theta = function(x) log(x),
                             from.theta = function(x) exp(x)
                         ), 
                         theta3 = list(
                             hyperid =  35103,
                             name = "log nu",
                             short.name = "nu",
                             initial = log(0.5),
                             fixed = TRUE,
                             prior = "loggamma",
                             param = c(0.5, 1),
                             to.theta = function(x) log(x),
                             from.theta = function(x) exp(x)
                         )
                     ),
                     constr = FALSE,
                     nrow.ncol = FALSE,
                     augmented = FALSE,
                     aug.factor = 1L,
                     aug.constr = NULL,
                     n.div.by = NULL,
                     n.required = TRUE,
                     set.default.values = TRUE,
                     status = "experimental", 
                     pdf = "dmatern"
                 ),

                 copy = list(
                     doc = "Create a copy of a model component", 
                     hyper = list(
                         theta = list(
                             hyperid =  36001,
                             name = "beta",
                             short.name = "b",
                             initial = 1,
                             fixed = TRUE,
                             prior = "normal",
                             param = c(1, 10),
                             to.theta = function(x, REPLACE.ME.low, REPLACE.ME.high) {
                             if (all(is.infinite(c(low, high))) || low == high) {
                                 return (x)
                             } else if (all(is.finite(c(low, high)))) {
                                 stopifnot(low < high)
                                 return (log( - (low - x)/(high -x)))
                             } else if (is.finite(low) && is.infinite(high) && high > low) {
                                 return (log(x-low))
                             } else {
                                 stop("Condition not yet implemented")
                             }
                         }, 
                         from.theta = function(x, REPLACE.ME.low, REPLACE.ME.high) {
                             if (all(is.infinite(c(low, high))) || low == high) {
                                 return (x)
                             } else if (all(is.finite(c(low, high)))) {
                                 stopifnot(low < high)
                                 return (low + exp(x)/(1+exp(x)) * (high - low))
                             } else if (is.finite(low) && is.infinite(high) && high > low) {
                                 return (low + exp(x))
                             } else {
                                 stop("Condition not yet implemented")
                             }
                         }
                         )
                     ),
                     constr = FALSE,
                     nrow.ncol = FALSE,
                     augmented = FALSE,
                     aug.factor = 1L,
                     aug.constr = NULL,
                     n.div.by = NULL,
                     n.required = FALSE,
                     set.default.values = FALSE,
                     pdf = NA
                 ), 

                 clinear = list(
                     doc = "Constrained linear effect", 
                     hyper = list(
                         theta = list(
                             hyperid =  37001,
                             name = "beta",
                             short.name = "b",
                             initial = 1,
                             fixed = FALSE,
                             prior = "normal",
                             param = c(1, 10),
                             to.theta = function(x, REPLACE.ME.low, REPLACE.ME.high) {
                             if (all(is.infinite(c(low, high))) || low == high) {
                                 stopifnot(low < high)
                                 return (x)
                             } else if (all(is.finite(c(low, high)))) {
                                 stopifnot(low < high)
                                 return (log( - (low - x)/(high -x)))
                             } else if (is.finite(low) && is.infinite(high) && high > low) {
                                 return (log(x-low))
                             } else {
                                 stop("Condition not yet implemented")
                             }
                         }, 
                         from.theta = function(x, REPLACE.ME.low, REPLACE.ME.high) {
                             if (all(is.infinite(c(low, high))) || low == high) {
                                 stopifnot(low < high)
                                 return (x)
                             } else if (all(is.finite(c(low, high)))) {
                                 stopifnot(low < high)
                                 return (low + exp(x)/(1+exp(x)) * (high - low))
                             } else if (is.finite(low) && is.infinite(high) && high > low) {
                                 return (low + exp(x))
                             } else {
                                 stop("Condition not yet implemented")
                             }
                         }
                         )
                     ),
                     constr = FALSE,
                     nrow.ncol = FALSE,
                     augmented = FALSE,
                     aug.factor = 1L,
                     aug.constr = NULL,
                     n.div.by = NULL,
                     n.required = FALSE,
                     set.default.values = FALSE,
                     pdf = "clinear"
                 ), 

                 sigm = list(
                     doc = "Sigmoidal effect of a covariate", 
                     hyper = list(
                         theta1 = list(
                             hyperid =  38001,
                             name = "beta",
                             short.name = "b",
                             initial = 1,
                             fixed = FALSE,
                             prior = "normal",
                             param = c(1, 10),
                             to.theta = function(x) x, 
                             from.theta = function(x) x
                         ),
                         theta2 = list(
                             hyperid =  38002,
                             name = "loghalflife",
                             short.name = "halflife",
                             initial = 3,
                             fixed = FALSE,
                             prior = "loggamma",
                             param = c(3, 1),
                             to.theta = function(x) log(x),
                             from.theta = function(x) exp(x)
                         ), 
                         theta3 = list(
                             hyperid =  38003,
                             name = "logshape",
                             short.name = "shape",
                             initial = 0,
                             fixed = FALSE,
                             prior = "loggamma",
                             param = c(10, 10),
                             to.theta = function(x) log(x),
                             from.theta = function(x) exp(x)
                         )
                     ),
                     constr = FALSE,
                     nrow.ncol = FALSE,
                     augmented = FALSE,
                     aug.factor = 1L,
                     aug.constr = NULL,
                     n.div.by = NULL,
                     n.required = FALSE,
                     set.default.values = FALSE,
                     status = "experimental", 
                     pdf = "sigm"
                 ), 

                 revsigm = list(
                     doc = "Reverse sigmoidal effect of a covariate", 
                     hyper = list(
                         theta1 = list(
                             hyperid =  39001,
                             name = "beta",
                             short.name = "b",
                             initial = 1,
                             fixed = FALSE,
                             prior = "normal",
                             param = c(1, 10),
                             to.theta = function(x) x, 
                             from.theta = function(x) x
                         ),
                         theta2 = list(
                             hyperid =  39002,
                             name = "loghalflife",
                             short.name = "halflife",
                             initial = 3,
                             fixed = FALSE,
                             prior = "loggamma",
                             param = c(3, 1),
                             to.theta = function(x) log(x),
                             from.theta = function(x) exp(x)
                         ), 
                         theta3 = list(
                             hyperid =  39003,
                             name = "logshape",
                             short.name = "shape",
                             initial = 0,
                             fixed = FALSE,
                             prior = "loggamma",
                             param = c(10, 10),
                             to.theta = function(x) log(x),
                             from.theta = function(x) exp(x)
                         )
                     ),
                     constr = FALSE,
                     nrow.ncol = FALSE,
                     augmented = FALSE,
                     aug.factor = 1L,
                     aug.constr = NULL,
                     n.div.by = NULL,
                     n.required = FALSE,
                     set.default.values = FALSE,
                     status = "experimental", 
                     pdf = "sigm"
                 ),

                 log1exp = list(
                     doc = "A nonlinear model of a covariate", 
                     hyper = list(
                         theta1 = list(
                             hyperid =  39011,
                             name = "beta",
                             short.name = "b",
                             initial = 1,
                             fixed = FALSE,
                             prior = "normal",
                             param = c(0, 1),
                             to.theta = function(x) x, 
                             from.theta = function(x) x
                         ),
                         theta2 = list(
                             hyperid =  39012,
                             name = "alpha",
                             short.name = "a",
                             initial = 0,
                             fixed = FALSE,
                             prior = "normal",
                             param = c(0, 1),
                             to.theta = function(x) x,
                             from.theta = function(x) x
                         ), 
                         theta3 = list(
                             hyperid =  39013,
                             name = "gamma",
                             short.name = "g",
                             initial = 0,
                             fixed = FALSE,
                             prior = "normal",
                             param = c(0, 1),
                             to.theta = function(x) x,
                             from.theta = function(x) x
                         )
                     ),
                     constr = FALSE,
                     nrow.ncol = FALSE,
                     augmented = FALSE,
                     aug.factor = 1L,
                     aug.constr = NULL,
                     n.div.by = NULL,
                     n.required = FALSE,
                     set.default.values = FALSE,
                     status = "experimental", 
                     pdf = "log1exp"
                 ),
                 
                 logdist = list(
                     doc = "A nonlinear model of a covariate", 
                     hyper = list(
                         theta1 = list(
                             hyperid =  39021,
                             name = "beta",
                             short.name = "b",
                             initial = 1,
                             fixed = FALSE,
                             prior = "normal",
                             param = c(0, 1),
                             to.theta = function(x) x, 
                             from.theta = function(x) x
                         ),
                         theta2 = list(
                             hyperid =  39022,
                             name = "alpha1",
                             short.name = "a1",
                             initial = 0,
                             fixed = FALSE,
                             prior = "loggamma",
                             param = c(0.1, 1),
                             to.theta = function(x) log(x),
                             from.theta = function(x) exp(x)
                         ), 
                         theta3 = list(
                             hyperid =  39023,
                             name = "alpha2",
                             short.name = "a2",
                             initial = 0,
                             fixed = FALSE,
                             prior = "loggamma",
                             param = c(0.1, 1),
                             to.theta = function(x) log(x),
                             from.theta = function(x) exp(x)
                         )
                     ),
                     constr = FALSE,
                     nrow.ncol = FALSE,
                     augmented = FALSE,
                     aug.factor = 1L,
                     aug.constr = NULL,
                     n.div.by = NULL,
                     n.required = FALSE,
                     set.default.values = FALSE,
                     status = "experimental", 
                     pdf = "logdist"
                 )
                 ##
             )
         )
}

`inla.models.section.group` = function()
{
    ## prevent a warning with R CMD check
    ngroup = NULL

    return
    list(group =
             list(
                 exchangeable = list(
                     doc = "Exchangeable correlations", 
                     hyper = list(
                         theta = list(
                             hyperid =  40001,
                             name = "logit correlation",
                             short.name = "rho",
                             initial = 1,
                             fixed = FALSE,
                             prior = "normal",
                             param = c(0, 0.2),
                             to.theta = function(x, REPLACE.ME.ngroup) log((1+x*(ngroup-1))/(1-x)), 
                             from.theta = function(x, REPLACE.ME.ngroup) (exp(x)-1)/(exp(x) + ngroup -1)
                             )
                         )
                     ),

                 exchangeablepos = list(
                     doc = "Exchangeable positive correlations", 
                     hyper = list(
                         theta = list(
                             hyperid =  40101,
                             name = "logit correlation",
                             short.name = "rho",
                             initial = 1,
                             fixed = FALSE,
                             prior = "pc.cor0",
                             param = c(0.5, 0.5),
                             to.theta = function(x) log(x/(1-x)), 
                             from.theta = function(x) exp(x)/(1+exp(x))
                             )
                         )
                     ),

                 ar1 = list(
                     doc = "AR(1) correlations", 
                     hyper = list(
                         theta = list(
                             hyperid =  41001,
                             name = "logit correlation",
                             short.name = "rho",
                             initial = 2,
                             fixed = FALSE,
                             prior = "normal",
                             param = c(0, 0.15),
                             to.theta = function(x) log((1+x)/(1-x)),
                             from.theta = function(x) 2*exp(x)/(1+exp(x))-1
                             )
                         )
                     ),
                 
                 ar = list(
                     doc = "AR(p) correlations", 
                     ## to many parameters here, but ...
                     hyper = list(
                         theta1 = list(
                             hyperid =  42001,
                             name = "log precision",
                             short.name = "prec",
                             initial = 0,
                             fixed = TRUE,
                             prior = "pc.prec",
                             param = c(3, 0.01),
                             to.theta = function(x) log(x),
                             from.theta = function(x) exp(x)
                             ),
                         theta2 = list(
                             hyperid =  42002,
                             name = "pacf1",
                             short.name = "pacf1",
                             initial = 2,
                             fixed = FALSE,
                             prior = "pc.cor0",
                             param = c(0.5, 0.5), 
                             to.theta = function(x) log((1+x)/(1-x)),
                             from.theta = function(x) 2*exp(x)/(1+exp(x))-1
                             ),
                         theta3 = list(
                             hyperid =  42003,
                             name = "pacf2",
                             short.name = "pacf2",
                             initial = 0,
                             fixed = FALSE,
                             prior = "pc.cor0",
                             param = c(0.5, 0.4),
                             to.theta = function(x) log((1+x)/(1-x)),
                             from.theta = function(x) 2*exp(x)/(1+exp(x))-1
                             ),
                         theta4 = list(
                             hyperid =  42004,
                             name = "pacf3",
                             short.name = "pacf3",
                             initial = 0,
                             fixed = FALSE,
                             prior = "pc.cor0",
                             param = c(0.5, 0.3),
                             to.theta = function(x) log((1+x)/(1-x)),
                             from.theta = function(x) 2*exp(x)/(1+exp(x))-1
                             ),
                         theta5 = list(
                             hyperid =  42005,
                             name = "pacf4",
                             short.name = "pacf4",
                             initial = 0,
                             fixed = FALSE,
                             prior = "pc.cor0",
                             param = c(0.5, 0.2),
                             to.theta = function(x) log((1+x)/(1-x)),
                             from.theta = function(x) 2*exp(x)/(1+exp(x))-1
                             ),
                         theta6 = list(
                             hyperid =  42006,
                             name = "pacf5",
                             short.name = "pacf5",
                             initial = 0,
                             fixed = FALSE,
                             prior = "pc.cor0",
                             param = c(0.5, 0.1),
                             to.theta = function(x) log((1+x)/(1-x)),
                             from.theta = function(x) 2*exp(x)/(1+exp(x))-1
                             ),
                         theta7 = list(
                             hyperid =  42007,
                             name = "pacf6",
                             short.name = "pacf6",
                             initial = 0,
                             fixed = FALSE,
                             prior = "pc.cor0",
                             param = c(0.5, 0.1),
                             to.theta = function(x) log((1+x)/(1-x)),
                             from.theta = function(x) 2*exp(x)/(1+exp(x))-1
                             ),
                         theta8 = list(
                             hyperid =  42008,
                             name = "pacf7",
                             short.name = "pacf7",
                             initial = 0,
                             fixed = FALSE,
                             prior = "pc.cor0",
                             param = c(0.5, 0.1),
                             to.theta = function(x) log((1+x)/(1-x)),
                             from.theta = function(x) 2*exp(x)/(1+exp(x))-1
                             ),
                         theta9 = list(
                             hyperid =  42009,
                             name = "pacf8",
                             short.name = "pacf8",
                             initial = 0,
                             fixed = FALSE,
                             prior = "pc.cor0",
                             param = c(0.5, 0.1),
                             to.theta = function(x) log((1+x)/(1-x)),
                             from.theta = function(x) 2*exp(x)/(1+exp(x))-1
                             ),
                         theta10 = list(
                             hyperid =  42010,
                             name = "pacf9",
                             short.name = "pacf9",
                             initial = 0,
                             fixed = FALSE,
                             prior = "pc.cor0",
                             param = c(0.5, 0.1),
                             to.theta = function(x) log((1+x)/(1-x)),
                             from.theta = function(x) 2*exp(x)/(1+exp(x))-1
                             ),
                         theta11 = list(
                             hyperid =  42011,
                             name = "pacf10",
                             short.name = "pacf10",
                             initial = 0,
                             fixed = FALSE,
                             prior = "pc.cor0",
                             param = c(0.5, 0.1),
                             to.theta = function(x) log((1+x)/(1-x)),
                             from.theta = function(x) 2*exp(x)/(1+exp(x))-1
                             )
                         )
                     ),

                 rw1 = list(
                     doc = "Random walk of order 1", 
                     hyper = list(
                         theta = list(
                             hyperid =  43001,
                             name = "log precision",
                             short.name = "prec",
                             prior = "loggamma",
                             param = c(1, 0.00005),
                             initial = 0,
                             fixed = TRUE,
                             to.theta = function(x) log(x),
                             from.theta = function(x) exp(x)
                             )
                         )
                     ), 

                 rw2 = list(
                     doc = "Random walk of order 2", 
                     hyper = list(
                         theta = list(
                             hyperid =  44001,
                             name = "log precision",
                             short.name = "prec",
                             prior = "loggamma",
                             param = c(1, 0.00005),
                             initial = 0,
                             fixed = TRUE,
                             to.theta = function(x) log(x),
                             from.theta = function(x) exp(x)
                             )
                         )
                     ), 

                 besag = list(
                     doc = "Besag model", 
                     hyper = list(
                         theta = list(
                             hyperid =  45001,
                             name = "log precision",
                             short.name = "prec",
                             prior = "loggamma",
                             param = c(1, 0.00005),
                             initial = 0,
                             fixed = TRUE,
                             to.theta = function(x) log(x),
                             from.theta = function(x) exp(x)
                             )
                         )
                     ),
                 iid = list(
                     doc = "Independent model", 
                     hyper = list(
                         theta = list(
                             hyperid =  46001,
                             name = "log precision",
                             short.name = "prec",
                             prior = "loggamma",
                             param = c(1, 0.00005),
                             initial = 0,
                             fixed = TRUE,
                             to.theta = function(x) log(x),
                             from.theta = function(x) exp(x)
                             )
                         )
                     ) 
                 )
         )
}

`inla.models.section.mix` = function()
{
    return
    list(mix =
             list(
                 gaussian = list(
                     doc = "Gaussian mixture", 
                     hyper = list(
                         theta = list(
                             hyperid =  47001,
                             name = "log precision",
                             short.name = "prec",
                             prior = "pc.prec",
                             param = c(1, 0.01),
                             initial = 0,
                             fixed = FALSE,
                             to.theta = function(x) log(x),
                             from.theta = function(x) exp(x)
                         )
                     )
                 ), 
                 loggamma = list(
                     ## > a = 1/inla.pc.rgamma(100000, lambda = 4.8)
                     ## > x = log(rgamma(length(a), shape=a, rate = a))
                     ## > sd (x)
                     ## [1] 0.3344429082
                     doc = "LogGamma mixture", 
                     hyper = list(
                         theta = list(
                             hyperid =  47101,
                             name = "log precision",
                             short.name = "prec",
                             prior = "pc.mgamma",
                             param = 4.8, 
                             initial = 4,
                             fixed = FALSE,
                             to.theta = function(x) log(x),
                             from.theta = function(x) exp(x)
                         )
                     )
                 ), 
                 mloggamma = list(
                     ## > a = 1/inla.pc.rgamma(100000, lambda = 4.8)
                     ## > x = log(rgamma(length(a), shape=a, rate = a))
                     ## > sd (x)
                     ## [1] 0.3344429082
                     doc = "Minus-LogGamma mixture", 
                     hyper = list(
                         theta = list(
                             hyperid =  47201,
                             name = "log precision",
                             short.name = "prec",
                             prior = "pc.mgamma",
                             param = 4.8, 
                             initial = 4,
                             fixed = FALSE,
                             to.theta = function(x) log(x),
                             from.theta = function(x) exp(x)
                         )
                     )
                 )
             )
         )
}

`inla.models.section.link` = function()
{
    return
    list(link =
             list(
                 default = list(
                     doc = "The default link", 
                     hyper = list()
                 ), 
                 cloglog = list(
                     doc = "The complementary log-log link", 
                     hyper = list()
                 ), 
                 loglog = list(
                     doc = "The log-log link", 
                     hyper = list()
                 ), 
                 identity = list(
                     doc = "The identity link", 
                     hyper = list()
                 ), 
                 inverse = list(
                     doc = "The inverse link", 
                     hyper = list()
                 ), 
                 log = list(
                     doc = "The log-link", 
                     hyper = list()
                 ), 
                 neglog = list(
                     doc = "The negative log-link", 
                     hyper = list()
                 ), 
                 logit = list(
                     doc = "The logit-link", 
                     hyper = list()
                 ), 
                 probit = list(
                     doc = "The probit-link", 
                     hyper = list()
                 ), 
                 cauchit = list(
                     doc = "The cauchit-link", 
                     hyper = list()
                 ), 
                 tan = list(
                     doc = "The tan-link", 
                     hyper = list()
                 ), 
                 quantile = list(
                     doc = "The quantile-link", 
                     hyper = list()
                 ), 
                 pquantile = list(
                     doc = "The population quantile-link", 
                     hyper = list()
                 ), 
                 sslogit = list(
                     doc = "Logit link with sensitivity and specificity", 
                     hyper = list(
                         theta1 = list(
                             hyperid =  48001,
                             name = "sensitivity",
                             short.name = "sens", 
                             prior = "logitbeta",
                             param = c(10, 5),
                             initial = 1,
                             fixed = FALSE,
                             to.theta = function(x) log(x/(1-x)),
                             from.theta = function(x) exp(x)/(1+exp(x))
                             ),
                         theta2 = list(
                             hyperid =  48002,
                             name = "specificity",
                             short.name = "spec",
                             prior = "logitbeta",
                             param = c(10, 5),
                             initial = 1,
                             fixed = FALSE,
                             to.theta = function(x) log(x/(1-x)),
                             from.theta = function(x) exp(x)/(1+exp(x))
                             )), 
                     pdf = NA
                     ), 

                 logoffset = list(
                     doc = "Log-link with an offset", 
                     ## variant = 0, a+exp(...), a>0
                     ## variant = 1, a-exp(...), a>0
                     hyper = list(
                         theta = list(
                             hyperid =  49001,
                             name = "beta",
                             short.name = "b",
                             prior = "normal",
                             param = c(0, 100),
                             initial = 0,
                             fixed = TRUE,
                             to.theta = function(x) log(x),
                             from.theta = function(x) exp(x)
                             )),
                     pdf = "logoffset"
                     ), 

                 logitoffset = list(
                     doc = "Logit-link with an offset", 
                     hyper = list(
                         theta = list(
                             hyperid =  49011,
                             name = "prob",
                             short.name = "p",
                             prior = "normal",
                             param = c(-1, 100),
                             initial = -1,
                             fixed = FALSE,
                             to.theta = function(x) log(x/(1-x)),
                             from.theta = function(x) exp(x)/(1+exp(x))
                             )),
                     status = "experimental", 
                     pdf = "logitoffset"
                     ), 
                 robit = list(
                     doc = "Robit link", 
                     hyper = list(
                         theta = list(
                             hyperid =  49021,
                             name = "log degrees of freedom",
                             short.name = "dof",
                             initial = log(5),
                             fixed = TRUE,
                             prior = "pc.dof",
                             param = c(50, 0.5),
                             to.theta = function(x) log(x-2),
                             from.theta = function(x) 2+exp(x)
                             )
                         ),
                     status = "experimental", 
                     pdf = "robit"
                     ), 
                 sn = list(
                     doc = "Skew-normal link", 
                     hyper = list(
                         theta = list(
                             hyperid =  49031,
                             name = "alpha",
                             short.name = "alpha",
                             initial = 0,
                             fixed = TRUE,
                             prior = "pc.sn",
                             param = 40,
                             to.theta = function(x) x,
                             from.theta = function(x) x
                             )
                         ),
                     status = "experimental", 
                     pdf = "linksn"
                     ), 
                 test1 = list(
                     doc = "A test1-link function (experimental)", 
                     hyper = list(
                         theta = list(
                             hyperid =  50001,
                             name = "beta",
                             short.name = "b",
                             prior = "normal",
                             param = c(0, 100),
                             initial = 0,
                             fixed = FALSE,
                             to.theta = function(x) x,
                             from.theta = function(x) x
                             )),
                     pdf = NA
                     ),

                 special1 = list(
                     doc = "A special1-link function (experimental)", 
                     hyper = list(
                         theta1 = list(
                             hyperid =  51001,
                             name = "log precision",
                             short.name = "prec",
                             initial = 0,
                             fixed = FALSE,
                             prior = "loggamma",
                             param = c(1, 1),
                             to.theta = function(x) x,
                             from.theta = function(x) x
                             ),
                         theta2 = list(
                             hyperid =  51002,
                             name = "beta1",
                             short.name = "beta1",
                             initial = 0,
                             fixed = FALSE,
                             prior = "mvnorm",
                             param = c(0, 100), 
                             to.theta = function(x) x, 
                             from.theta = function(x) x
                             ),
                         theta3 = list(
                             hyperid =  51003,
                             name = "beta2",
                             short.name = "beta2",
                             initial = 0,
                             fixed = FALSE,
                             prior = "none",
                             param = numeric(0),
                             to.theta = function(x) x, 
                             from.theta = function(x) x
                             ),
                         theta4 = list(
                             hyperid =  51004,
                             name = "beta3",
                             short.name = "beta3",
                             initial = 0,
                             fixed = FALSE,
                             prior = "none",
                             param = numeric(0),
                             to.theta = function(x) x, 
                             from.theta = function(x) x
                             ),
                         theta5 = list(
                             hyperid =  51005,
                             name = "beta4",
                             short.name = "beta4",
                             initial = 0,
                             fixed = FALSE,
                             prior = "none",
                             param = numeric(0),
                             to.theta = function(x) x, 
                             from.theta = function(x) x
                             ),
                         theta6 = list(
                             hyperid =  51006,
                             name = "beta5",
                             short.name = "beta5",
                             initial = 0,
                             fixed = FALSE,
                             prior = "none",
                             param = numeric(0),
                             to.theta = function(x) x, 
                             from.theta = function(x) x
                             ),
                         theta7 = list(
                             hyperid =  51007,
                             name = "beta6",
                             short.name = "beta6",
                             initial = 0,
                             fixed = FALSE,
                             prior = "none",
                             param = numeric(0),
                             to.theta = function(x) x, 
                             from.theta = function(x) x
                             ),
                         theta8 = list(
                             hyperid =  51008,
                             name = "beta7",
                             short.name = "beta7",
                             initial = 0,
                             fixed = FALSE,
                             prior = "none",
                             param = numeric(0),
                             to.theta = function(x) x, 
                             from.theta = function(x) x
                             ),
                         theta9 = list(
                             hyperid =  51009,
                             name = "beta8",
                             short.name = "beta8",
                             initial = 0,
                             fixed = FALSE,
                             prior = "none",
                             param = numeric(0),
                             to.theta = function(x) x, 
                             from.theta = function(x) x
                             ),
                         theta10 = list(
                             hyperid =  51010,
                             name = "beta9",
                             short.name = "beta9",
                             initial = 0,
                             fixed = FALSE,
                             prior = "none",
                             param = numeric(0),
                             to.theta = function(x) x, 
                             from.theta = function(x) x
                             ),
                         theta11 = list(
                             hyperid =  51011,
                             name = "beta10",
                             short.name = "beta10",
                             initial = 0,
                             fixed = FALSE,
                             prior = "none",
                             param = numeric(0),
                             to.theta = function(x) x, 
                             from.theta = function(x) x
                             )
                         ),
                     pdf = NA
                     ),
                 
                 special2 = list(
                     doc = "A special2-link function (experimental)", 
                     hyper = list(
                         theta = list(
                             hyperid =  52001,
                             name = "beta",
                             short.name = "b",
                             prior = "normal",
                             param = c(0, 10),
                             initial = 0,
                             fixed = FALSE,
                             to.theta = function(x) x,
                             from.theta = function(x) x
                             )),
                     pdf = NA
                     )
                 )
         )
}

`inla.models.section.predictor` = function()
{
    return
    list(predictor =
             list(
                 predictor = list(
                     doc = "(not used)", 
                     hyper = list(
                         theta = list(
                             hyperid =  53001,
                             name = "log precision",
                             short.name = "prec",
                             initial = 12,
                             fixed = TRUE,
                             prior = "loggamma",
                             param = c(1, 0.00001),
                             to.theta = function(x) log(x),
                             from.theta = function(x) exp(x)
                             )
                         )
                     )
                 )
         )
}

`inla.models.section.hazard` = function()
{
    return
    list(hazard =
             list(
                 rw1 = list(
                     doc = "A random walk of order 1 for the log-hazard", 
                     hyper = list(
                         theta = list(
                             hyperid =  54001,
                             name = "log precision",
                             short.name = "prec",
                             initial = 4,
                             fixed = FALSE,
                             prior = "loggamma",
                             param = c(1, 0.00005),
                             to.theta = function(x) log(x),
                             from.theta = function(x) exp(x)
                             )
                         )
                     ),
                 rw2 = list(
                     doc = "A random walk of order 2 for the log-hazard", 
                     hyper = list(
                         theta = list(
                             hyperid =  55001,
                             name = "log precision",
                             short.name = "prec",
                             initial = 4,
                             fixed = FALSE,
                             prior = "loggamma",
                             param = c(1, 0.00005),
                             to.theta = function(x) log(x),
                             from.theta = function(x) exp(x)
                             )
                         )
                     )
                 )
         )
}

`inla.models.section.likelihood` = function()
{
    return
    list(likelihood =
             list(
                 ## the first non-default link-function is the default one.
                 poisson = list(
                     doc = "The Poisson likelihood", 
                     hyper = list(
                         ),
                     survival = FALSE,
                     discrete = TRUE,
                     link = c("default", "log", "logoffset", "quantile", "test1", "special1", "special2"),
                     pdf = "poisson"
                     ),
                 
                 contpoisson = list(
                     doc = "The Cont Poisson likelihood", 
                     hyper = list(
                         ),
                     survival = FALSE,
                     discrete = TRUE,
                     link = c("default", "log"), 
                     pdf = "contpoisson"
                     ),
                 
                 qcontpoisson = list(
                     doc = "The quantile Cont Poisson likelihood", 
                     hyper = list(
                         ),
                     survival = FALSE,
                     discrete = TRUE,
                     link = c("default", "log"), 
                     pdf = "qcontpoisson"
                     ),
                 
                 cenpoisson = list(
                     doc = "Then censored Poisson likelihood", 
                     hyper = list(
                         ),
                     survival = FALSE,
                     discrete = TRUE,
                     link = c("default", "log", "logoffset", "test1", "special1", "special2"),
                     status = "experimental", 
                     pdf = "cenpoisson"
                     ),

                 gpoisson = list(
                     doc = "The generalized Poisson likelihood", 
                     hyper = list(
                         theta1 = list(
                             hyperid =  56001,
                             name = "overdispersion",
                             short.name = "phi",
                             initial = 0,
                             fixed = FALSE,
                             prior = "loggamma",
                             param = c(1, 1),
                             to.theta = function(x) log(x), 
                             from.theta = function(x) exp(x)
                             ), 
                         theta2 = list(
                             hyperid =  56002,
                             name = "p",
                             short.name = "p",
                             initial = 1,
                             fixed = TRUE,
                             prior = "normal",
                             ## use a tight prior, as this can go very wrong if set to weak
                             param = c(1, 100),
                             to.theta = function(x) x, 
                             from.theta = function(x) x
                             )
                         ),
                     survival = FALSE,
                     discrete = TRUE,
                     link = c("default", "log", "logoffset"),
                     pdf = "gpoisson",
                     status = "experimental"
                     ),

                 binomial = list(
                     doc = "The Binomial likelihood", 
                     hyper = list(
                         ),
                     survival = FALSE,
                     discrete = TRUE,
                     link = c("default", "logit", "cauchit", "probit", "cloglog", "loglog",
                              "log", "sslogit", "logitoffset", "quantile", "pquantile", "robit", "sn"),
                     pdf = "binomial"
                     ),

                 pom = list(
                     doc = "Likelihood for the proportional odds model", 
                     hyper = list(
                         theta1 = list(
                             hyperid =  57101,
                             name = "theta1",
                             short.name = "theta1",
                             initial = NA,
                             fixed = FALSE,
                             prior = "dirichlet",
                             param = 3.0, 
                             to.theta = function(x) x, 
                             from.theta = function(x) x
                             ), 
                         theta2 = list(
                             hyperid =  57102,
                             name = "theta2",
                             short.name = "theta2",
                             initial = NA,
                             fixed = FALSE,
                             prior = "none",
                             param = numeric(0), 
                             to.theta = function(x) log(x), 
                             from.theta = function(x) exp(x)
                             ), 
                         theta3 = list(
                             hyperid =  57103,
                             name = "theta3",
                             short.name = "theta3",
                             initial = NA,
                             fixed = FALSE,
                             prior = "none",
                             param = numeric(0), 
                             to.theta = function(x) log(x), 
                             from.theta = function(x) exp(x)
                             ), 
                         theta4 = list(
                             hyperid =  57104,
                             name = "theta4",
                             short.name = "theta4",
                             initial = NA,
                             fixed = FALSE,
                             prior = "none",
                             param = numeric(0), 
                             to.theta = function(x) log(x), 
                             from.theta = function(x) exp(x) 
                             ), 
                         theta5 = list(
                             hyperid =  57105,
                             name = "theta5",
                             short.name = "theta5",
                             initial = NA,
                             fixed = FALSE,
                             prior = "none",
                             param = numeric(0), 
                             to.theta = function(x) log(x), 
                             from.theta = function(x) exp(x) 
                             ), 
                         theta6 = list(
                             hyperid =  57106,
                             name = "theta6",
                             short.name = "theta6",
                             initial = NA,
                             fixed = FALSE,
                             prior = "none",
                             param = numeric(0), 
                             to.theta = function(x) log(x), 
                             from.theta = function(x) exp(x) 
                             ), 
                         theta7 = list(
                             hyperid =  57107,
                             name = "theta7",
                             short.name = "theta7",
                             initial = NA,
                             fixed = FALSE,
                             prior = "none",
                             param = numeric(0), 
                             to.theta = function(x) log(x), 
                             from.theta = function(x) exp(x) 
                             ), 
                         theta8 = list(
                             hyperid =  57108,
                             name = "theta8",
                             short.name = "theta8",
                             initial = NA,
                             fixed = FALSE,
                             prior = "none",
                             param = numeric(0), 
                             to.theta = function(x) log(x), 
                             from.theta = function(x) exp(x) 
                             ), 
                         theta9 = list(
                             hyperid =  57109,
                             name = "theta9",
                             short.name = "theta9",
                             initial = NA,
                             fixed = FALSE,
                             prior = "none",
                             param = numeric(0), 
                             to.theta = function(x) log(x), 
                             from.theta = function(x) exp(x) 
                             ), 
                         theta10 = list(
                             hyperid =  57110,
                             name = "theta10",
                             short.name = "theta10",
                             initial = NA,
                             fixed = FALSE,
                             prior = "none",
                             param = numeric(0), 
                             to.theta = function(x) log(x), 
                             from.theta = function(x) exp(x) 
                             )
                         ),
                     status = "experimental", 
                     survival = FALSE,
                     discrete = TRUE,
                     link = c("default", "identity"), 
                     pdf = "pom"
                     ),

                 gev2 = list(
                     doc = "The Generalized Extreme Value likelihood (2nd variant)", 
                     hyper = list(
                         theta1 = list(
                             hyperid =  57201,
                             name = "spread",
                             short.name = "sd",
                             initial = 0,
                             fixed = FALSE,
                             prior = "loggamma",
                             param = c(1, 3),
                             to.theta = function(x) log(x),
                             from.theta = function(x) exp(x)
                         ),
                         theta2 = list(
                             hyperid =  57202,
                             name = "tail",
                             short.name = "tail",
                             initial = -4,
                             fixed = FALSE,
                             prior = "pc.gevtail",
                             param = c(7, 0.0, 0.5), 
                             to.theta = function(x, interval = c(REPLACE.ME.low, REPLACE.ME.high)) log(-(interval[1] - x)/(interval[2] - x)), 
                             from.theta = function(x, interval = c(REPLACE.ME.low, REPLACE.ME.high)) interval[1] + (interval[2]-interval[1]) * exp(x)/(1.0 + exp(x))
                         ), 
                         theta3 = list(
                             hyperid =  57203,
                             name = "beta1",
                             short.name = "beta1",
                             initial = NA,
                             fixed = FALSE,
                             prior = "normal",
                             param = c(0, 300), 
                             to.theta = function(x) x, 
                             from.theta = function(x) x
                         ), 
                         theta4 = list(
                             hyperid =  57204,
                             name = "beta2",
                             short.name = "beta2",
                             initial = NA,
                             fixed = FALSE,
                             prior = "normal",
                             param = c(0, 300), 
                             to.theta = function(x) x, 
                             from.theta = function(x) x
                         ), 
                         theta5 = list(
                             hyperid =  57205,
                             name = "beta3",
                             short.name = "beta3",
                             initial = NA,
                             fixed = FALSE,
                             prior = "normal",
                             param = c(0, 300), 
                             to.theta = function(x) x, 
                             from.theta = function(x) x
                         ), 
                         theta6 = list(
                             hyperid =  57206,
                             name = "beta4",
                             short.name = "beta4",
                             initial = NA,
                             fixed = FALSE,
                             prior = "normal",
                             param = c(0, 300), 
                             to.theta = function(x) x, 
                             from.theta = function(x) x
                         ), 
                         theta7 = list(
                             hyperid =  57207,
                             name = "beta5",
                             short.name = "beta5",
                             initial = NA,
                             fixed = FALSE,
                             prior = "normal",
                             param = c(0, 300), 
                             to.theta = function(x) x, 
                             from.theta = function(x) x
                         ), 
                         theta8 = list(
                             hyperid =  57208,
                             name = "beta6",
                             short.name = "beta6",
                             initial = NA,
                             fixed = FALSE,
                             prior = "normal",
                             param = c(0, 300), 
                             to.theta = function(x) x, 
                             from.theta = function(x) x
                         ), 
                         theta9 = list(
                             hyperid =  57209,
                             name = "beta7",
                             short.name = "beta7",
                             initial = NA,
                             fixed = FALSE,
                             prior = "normal",
                             param = c(0, 300), 
                             to.theta = function(x) x, 
                             from.theta = function(x) x
                         ), 
                         theta10 = list(
                             hyperid =  57210,
                             name = "beta8",
                             short.name = "beta8",
                             initial = NA,
                             fixed = FALSE,
                             prior = "normal",
                             param = c(0, 300), 
                             to.theta = function(x) x, 
                             from.theta = function(x) x
                         ), 
                         theta11 = list(
                             hyperid =  57211,
                             name = "beta9",
                             short.name = "beta9",
                             initial = NA,
                             fixed = FALSE,
                             prior = "normal",
                             param = c(0, 300), 
                             to.theta = function(x) x, 
                             from.theta = function(x) x
                         ), 
                         theta12 = list(
                             hyperid =  57212,
                             name = "beta10",
                             short.name = "beta",
                             initial = NA,
                             fixed = FALSE,
                             prior = "normal",
                             param = c(0, 300), 
                             to.theta = function(x) x, 
                             from.theta = function(x) x
                         )
                     ), 
                     status = "experimental", 
                     survival = FALSE,
                     discrete = FALSE,
                     link = c("default", "identity"), 
                     pdf = "gev2"
                     ),

                 gamma = list(
                     doc = "The Gamma likelihood", 
                     hyper = list(
                         theta = list(
                             hyperid =  58001,
                             name = "precision parameter",
                             short.name = "prec",
                             initial = log(100),
                             fixed = FALSE,
                             prior = "loggamma",
                             param = c(1, 0.01),
                             to.theta = function(x) log(x), 
                             from.theta = function(x) exp(x)
                             )
                         ),
                     survival = FALSE,
                     discrete = FALSE,
                     link = c("default", "log", "quantile"),
                     pdf = "gamma"
                     ),

                 gammasurv = list(
                     doc = "The Gamma likelihood (survival)", 
                     hyper = list(
                         theta = list(
                             hyperid =  58101,
                             name = "precision parameter",
                             short.name = "prec",
                             initial = log(1),
                             fixed = FALSE,
                             prior = "loggamma",
                             param = c(1, 0.01),
                             to.theta = function(x) log(x), 
                             from.theta = function(x) exp(x)
                             )
                         ),
                     survival = TRUE,
                     discrete = FALSE,
                     status = "experimental", 
                     link = c("default", "log", "quantile"),
                     pdf = "gammasurv"
                     ),

                 gammacount = list(
                     doc = "A Gamma generalisation of the Poisson likelihood", 
                     hyper = list(
                         theta = list(
                             hyperid =  59001,
                             name = "log alpha",
                             short.name = "alpha",
                             initial = log(1.0),
                             fixed = FALSE,
                             prior = "pc.gammacount",
                             param = 3,
                             to.theta = function(x) log(x), 
                             from.theta = function(x) exp(x)
                             )
                         ),
                     survival = FALSE,
                     discrete = FALSE,
                     link = c("default", "log"),
                     status = "experimental", 
                     pdf = "gammacount"
                     ),

                 qkumar = list(
                     doc = "A quantile version of the Kumar likelihood", 
                     hyper = list(
                         theta = list(
                             hyperid =  60001,
                             name = "precision parameter",
                             short.name = "prec",
                             initial = 0,
                             fixed = FALSE,
                             prior = "loggamma",
                             param = c(1, 0.001),
                             to.theta = function(x) log(x), 
                             from.theta = function(x) exp(x)
                         )
                     ), 
                     survival = FALSE,
                     discrete = FALSE,
                     link = c("default", "logit", "cauchit"),
                     pdf = "qkumar"
                 ),

                 qloglogistic = list(
                     doc = "A quantile loglogistic likelihood", 
                     hyper = list(
                         theta = list(
                             hyperid =  60011,
                             name = "log alpha",
                             short.name = "alpha",
                             initial = 1,
                             fixed = FALSE,
                             prior = "loggamma",
                             param = c(25, 25),
                             to.theta = function(x) log(x), 
                             from.theta = function(x) exp(x)
                         )
                     ), 
                     survival = FALSE,
                     discrete = FALSE,
                     link = c("default", "log", "neglog"), 
                     pdf = "qloglogistic"
                 ),

                 qloglogisticsurv = list(
                     doc = "A quantile loglogistic likelihood (survival)", 
                     hyper = list(
                         theta = list(
                             hyperid =  60021,
                             name = "log alpha",
                             short.name = "alpha",
                             initial = 1,
                             fixed = FALSE,
                             prior = "loggamma",
                             param = c(25, 25),
                             to.theta = function(x) log(x), 
                             from.theta = function(x) exp(x)
                         )
                     ), 
                     survival = TRUE,
                     discrete = FALSE,
                     link = c("default", "log", "neglog"), 
                     pdf = "qloglogistic"
                 ),

                 beta = list(
                     doc = "The Beta likelihood", 
                     hyper = list(
                         theta = list(
                             hyperid =  61001,
                             name = "precision parameter",
                             short.name = "phi",
                             initial = log(10),
                             fixed = FALSE,
                             prior = "loggamma",
                             param = c(1, 0.1),
                             to.theta = function(x) log(x), 
                             from.theta = function(x) exp(x)
                             )
                         ),
                     survival = FALSE,
                     discrete = FALSE,
                     link = c("default", "logit", "cauchit", "probit", "cloglog", "loglog"),
                     pdf = "beta"
                     ),

                 betabinomial = list(
                     doc = "The Beta-Binomial likelihood", 
                     hyper = list(
                         theta = list(
                             hyperid =  62001,
                             name = "overdispersion",
                             short.name = "rho",
                             initial = 0,
                             fixed = FALSE,
                             prior = "gaussian",
                             param = c(0.0, 0.4),
                             to.theta = function(x) log(x/(1-x)),
                             from.theta = function(x) exp(x)/(1+exp(x))
                             )
                         ),
                     survival = FALSE,
                     discrete = TRUE,
                     link = c("default", "logit", "cauchit", "probit", "cloglog", "loglog", "robit", "sn"),
                     pdf = "betabinomial"
                     ),

                 cbinomial = list(
                     doc = "The clustered Binomial likelihood", 
                     hyper = list(),
                     survival = FALSE,
                     discrete = TRUE,
                     link = c("default", "logit", "cauchit", "probit", "cloglog", "loglog", "robit", "sn"),
                     status = "experimental", 
                     pdf = "cbinomial"
                     ),

                 nbinomial = list(
                     doc = "The negBinomial likelihood", 
                     hyper = list(
                         theta = list(
                             hyperid =  63001,
                             name = "size",
                             short.name = "size",
                             initial = log(10),
                             fixed = FALSE,
                             prior = "pc.mgamma",
                             param = 7,
                             to.theta = function(x) log(x),
                             from.theta = function(x) exp(x)
                             )
                         ),
                     survival = FALSE,
                     discrete = TRUE,
                     link = c("default", "log", "logoffset", "quantile"),
                     pdf = "nbinomial"
                     ),

                 nbinomial2 = list(
                     doc = "The negBinomial2 likelihood", 
                     hyper = list(), 
                     survival = FALSE,
                     discrete = TRUE,
                     link = c("default", "logit", "cauchit", "probit", "cloglog", "loglog"), 
                     pdf = "nbinomial"
                     ),

                 simplex = list(
                     doc = "The simplex likelihood", 
                     hyper = list(
                         theta = list(
                             hyperid =  64001,
                             name = "log precision",
                             short.name = "prec",
                             initial = 4,
                             fixed = FALSE,
                             prior = "loggamma",
                             param = c(1, 0.00005),
                             to.theta = function(x) log(x),
                             from.theta = function(x) exp(x)
                             )
                         ),
                     survival = FALSE,
                     discrete = FALSE,
                     link = c("default", "logit", "cauchit", "probit", "cloglog", "loglog"),
                     pdf = "simplex"
                     ),

                 gaussian = list(
                     doc = "The Gaussian likelihoood", 
                     hyper = list(
                         theta1 = list(
                             hyperid =  65001,
                             name = "log precision",
                             short.name = "prec",
                             initial = 4,
                             fixed = FALSE,
                             prior = "loggamma",
                             param = c(1, 0.00005),
                             to.theta = function(x) log(x),
                             from.theta = function(x) exp(x)
                             ), 
                         theta2 = list(
                             hyperid =  65002,
                             name = "log precision offset",
                             short.name = "precoffset",
                             initial = (-2.0 * log(.Machine$double.eps)), 
                             fixed = TRUE,
                             prior = "none",
                             param = numeric(), 
                             to.theta = function(x) log(x),
                             from.theta = function(x) exp(x)
                             )
                         ),
                     survival = FALSE,
                     discrete = FALSE,
                     link = c("default", "identity", "logit", "cauchit", "log", "logoffset"),
                     pdf = "gaussian"
                     ),

                 circularnormal = list(
                     doc = "The circular Gaussian likelihoood", 
                     hyper = list(
                         theta = list(
                             hyperid =  67001,
                             name = "log precision parameter",
                             short.name = "prec",
                             initial = 2,
                             fixed = FALSE,
                             prior = "loggamma",
                             param = c(1, 0.01),
                             to.theta = function(x) log(x),
                             from.theta = function(x) exp(x)
                             )
                         ),
                     survival = FALSE,
                     discrete = FALSE,
                     link = c("default", "tan"),
                     pdf = "circular-normal",
                     status = "experimental"
                     ),
                 
                 wrappedcauchy = list(
                     doc = "The wrapped Cauchy likelihoood", 
                     hyper = list(
                         theta = list(
                             hyperid =  68001,
                             name = "log precision parameter",
                             short.name = "prec",
                             initial = 2,
                             fixed = FALSE,
                             prior = "loggamma",
                             param = c(1, 0.005),
                             to.theta = function(x) log(x/(1-x)), 
                             from.theta = function(x) exp(x)/(1+exp(x))
                             )
                         ),
                     survival = FALSE,
                     discrete = FALSE,
                     link = c("default", "tan"),
                     pdf = "wrapped-cauchy",
                     status = "disabled"
                     ),
                 
                 iidgamma = list(
                     doc = "(experimental)", 
                     hyper = list(
                         theta1 = list(
                             hyperid =  69001,
                             name = "logshape",
                             short.name = "shape",
                             initial = 0,
                             fixed = FALSE,
                             prior = "loggamma",
                             param = c(100, 100),
                             to.theta = function(x) log(x),
                             from.theta = function(x) exp(x)
                             ),
                         theta2 = list(
                             hyperid =  69002,
                             name = "lograte",
                             short.name = "rate",
                             initial = 0,
                             fixed = FALSE,
                             prior = "loggamma",
                             param = c(100, 100),
                             to.theta = function(x) log(x),
                             from.theta = function(x) exp(x)
                             )
                         ),
                     survival = FALSE,
                     discrete = FALSE,
                     link = c("default", "identity"),
                     pdf = "iidgamma",
                     status = "experimental"
                     ),

                 iidlogitbeta = list(
                     doc = "(experimental)", 
                     hyper = list(
                         theta1 = list(
                             hyperid =  70001,
                             name = "log.a",
                             short.name = "a",
                             initial = 1,
                             fixed = FALSE,
                             prior = "loggamma",
                             param = c(1, 1),
                             to.theta = function(x) log(x),
                             from.theta = function(x) exp(x)
                             ),
                         theta2 = list(
                             hyperid =  70002,
                             name = "log.b",
                             short.name = "b",
                             initial = 1,
                             fixed = FALSE,
                             prior = "loggamma",
                             param = c(1, 1),
                             to.theta = function(x) log(x),
                             from.theta = function(x) exp(x)
                             )
                         ),
                     survival = FALSE,
                     discrete = FALSE,
                     link = c("default", "logit"),
                     pdf = "iidlogitbeta",
                     status = "experimental"
                     ),

                 loggammafrailty = list(
                     doc = "(experimental)", 
                     hyper = list(
                         theta = list(
                             hyperid =  71001,
                             name = "log precision",
                             short.name = "prec",
                             initial = 4,
                             fixed = FALSE,
                             prior = "loggamma",
                             param = c(1, 0.00005),
                             to.theta = function(x) log(x),
                             from.theta = function(x) exp(x)
                             )
                         ),
                     survival = FALSE,
                     discrete = FALSE,
                     link = c("default", "identity"),
                     pdf = "loggammafrailty", 
                     status = "experimental"
                     ),

                 logistic = list(
                     doc = "The Logistic likelihoood", 
                     hyper = list(
                         theta = list(
                             hyperid =  72001,
                             name = "log precision",
                             short.name = "prec",
                             initial = 1,
                             fixed = FALSE,
                             prior = "loggamma",
                             param = c(1, 0.00005),
                             to.theta = function(x) log(x),
                             from.theta = function(x) exp(x)
                             )
                         ),
                     survival = FALSE,
                     discrete = FALSE,
                     link = c("default", "identity"),
                     pdf = "logistic"
                     ),

                 skewnormal = list(
                     doc = "The Skew-Normal likelihoood", 
                     hyper = list(
                         theta1 = list(
                             hyperid =  73001,
                             name = "log inverse scale",
                             short.name = "iscale",
                             initial = 4,
                             fixed = FALSE,
                             prior = "loggamma",
                             param = c(1, 0.00005), 
                             to.theta = function(x) log(x),
                             from.theta = function(x) exp(x)
                         ),
                         theta2 = list(
                             hyperid =  73002,
                             name = "logit skewness",
                             short.name = "skew",
                             initial = 4,
                             fixed = FALSE,
                             prior = "gaussian",
                             param = c(0, 10), 
                             to.theta = function(x, shape.max = 1) log((1+x/shape.max)/(1-x/shape.max)),
                             from.theta = function(x, shape.max = 1) shape.max*(2*exp(x)/(1+exp(x))-1)
                         )
                     ),
                     survival = FALSE,
                     discrete = FALSE,
                     link = c("default", "identity"),
                     pdf = "sn"
                 ),

                 sn = list(
                     doc = "The Skew-Normal likelihoood", 
                     hyper = list(
                         theta1 = list(
                             hyperid =  74001,
                             name = "log inverse scale",
                             short.name = "iscale",
                             initial = 4,
                             fixed = FALSE,
                             prior = "loggamma",
                             param = c(1, 0.00005), 
                             to.theta = function(x) log(x),
                             from.theta = function(x) exp(x)
                         ),
                         theta2 = list(
                             hyperid =  74002,
                             name = "logit skewness",
                             short.name = "skew",
                             initial = 0,
                             fixed = FALSE,
                             prior = "gaussian",
                             param = c(0, 10),
                             to.theta = function(x, shape.max = 1) log((1+x/shape.max)/(1-x/shape.max)),
                             from.theta = function(x, shape.max = 1) shape.max*(2*exp(x)/(1+exp(x))-1)
                             )
                         ),
                     survival = FALSE,
                     discrete = FALSE,
                     link = c("default", "identity"),
                     pdf = "sn"
                     ),

                 sn2 = list(
                     doc = "The Skew-Normal likelihoood (alt param)", 
                     hyper = list(
                         theta1 = list(
                             hyperid =  75001,
                             name = "log precision",
                             short.name = "prec",
                             initial = 1,
                             fixed = FALSE,
                             prior = "loggamma",
                             param = c(1, 0.00005), 
                             to.theta = function(x) log(x),
                             from.theta = function(x) exp(x)
                         ),
                         theta2 = list(
                             hyperid =  75002,
                             name = "logit skewness",
                             short.name = "skew",
                             initial = 0,
                             fixed = FALSE,
                             prior = "gaussian",
                             param = c(0, 10),
                             to.theta = function(x, shape.max = 1) log((1+x/shape.max)/(1-x/shape.max)),
                             from.theta = function(x, shape.max = 1) shape.max*(2*exp(x)/(1+exp(x))-1)
                         )
                     ),
                     survival = FALSE,
                     discrete = FALSE,
                     link = c("default", "identity"),
                     status = "experimental", 
                     pdf = "sn2"
                 ),

                 gev = list(
                     doc = "The Generalized Extreme Value likelihood", 
                     hyper = list(
                         theta1 = list(
                             hyperid =  76001,
                             name = "log precision",
                             short.name = "prec",
                             initial = 4,
                             fixed = FALSE,
                             prior = "loggamma",
                             param = c(1, 0.00005),
                             to.theta = function(x) log(x),
                             from.theta = function(x) exp(x)
                             ),
                         theta2 = list(
                             hyperid =  76002,
                             name = "tail parameter",
                             short.name = "tail",
                             initial = 0,
                             fixed = FALSE,
                             prior = "gaussian",
                             param = c(0, 25),
                             to.theta = function(x) x,
                             from.theta = function(x) x
                             )
                         ),
                     survival = FALSE,
                     discrete = FALSE,
                     link = c("default", "identity"),
                     status = "experimental", 
                     pdf = "gev"
                     ),

                 lognormal = list(
                     doc = "The log-Normal likelihood", 
                     hyper = list(
                         theta = list(
                             hyperid =  77101,
                             name = "log precision",
                             short.name = "prec",
                             initial = 0,
                             fixed = FALSE,
                             prior = "loggamma",
                             param = c(1, 0.00005),
                             to.theta = function(x) log(x),
                             from.theta = function(x) exp(x)
                             )
                         ),
                     survival = FALSE,
                     discrete = FALSE,
                     link = c("default", "identity"),
                     pdf = "lognormal"
                     ),

                 lognormalsurv = list(
                     doc = "The log-Normal likelihood (survival)", 
                     hyper = list(
                         theta = list(
                             hyperid =  78001,
                             name = "log precision",
                             short.name = "prec",
                             initial = 0,
                             fixed = FALSE,
                             prior = "loggamma",
                             param = c(1, 0.00005),
                             to.theta = function(x) log(x),
                             from.theta = function(x) exp(x)
                             )
                         ),
                     survival = TRUE,
                     discrete = FALSE,
                     link = c("default", "identity"),
                     pdf = "lognormal"
                     ),

                 exponential = list(
                     doc = "The Exponential likelihood", 
                     hyper = list(
                         ),
                     survival = FALSE,
                     discrete = FALSE,
                     link = c("default", "log"),
                     pdf = "exponential"
                     ),

                 exponentialsurv = list(
                     doc = "The Exponential likelihood (survival)", 
                     hyper = list(
                         ),
                     survival = TRUE,
                     discrete = FALSE,
                     link = c("default", "log", "neglog"),
                     pdf = "exponential"
                     ),

                 coxph = list(
                     doc = "Cox-proportional hazard likelihood", 
                     hyper = list(
                         ),
                     survival = TRUE,
                     discrete = TRUE,
                     link = c("default", "log", "neglog"),
                     pdf = "coxph"
                     ),

                 weibull = list(
                     doc = "The Weibull likelihood", 
                     ## variant=0: lambda*y^alpha
                     ## variant=1: (lambda*y)^alpha
                     hyper = list(
                         theta = list(
                             hyperid =  79001,
                             name = "log alpha",
                             short.name = "alpha",
                             initial = 0.1,
                             fixed = FALSE,
                             prior = "pc.alphaw",
                             param = c(5),
                             ## the 'sc' constant is defined in inla.h, and must be the same.
                             ## I know, this is hard-coded for the moment. Should be a generic
                             ## way of doing this...
                             to.theta = function(x, sc = 0.1) log(x)/sc,
                             from.theta = function(x, sc = 0.1) exp(sc*x)
                             )
                         ),
                     survival = FALSE,
                     discrete = FALSE,
                     link = c("default", "log", "neglog", "quantile"),
                     pdf = "weibull"
                     ),

                 weibullsurv = list(
                     doc = "The Weibull likelihood (survival)", 
                     ## variant=0: lambda*y^alpha
                     ## variant=1: (lambda*y)^alpha
                     hyper = list(
                         theta = list(
                             hyperid =  79101,
                             name = "log alpha",
                             short.name = "alpha",
                             initial = 0.1,
                             fixed = FALSE,
                             prior = "pc.alphaw",
                             param = c(5),
                             ## the 'sc' constant is defined in inla.h, and must be the same.
                             ## I know, this is hard-coded for the moment. Should be a generic
                             ## way of doing this...
                             to.theta = function(x, sc = 0.1) log(x)/sc,
                             from.theta = function(x, sc = 0.1) exp(sc*x)
                             )
                         ),
                     survival = TRUE,
                     discrete = FALSE,
                     link = c("default", "log", "neglog", "quantile"),
                     pdf = "weibull"
                     ),

                 loglogistic = list(
                     doc = "The loglogistic likelihood", 
                     hyper = list(
                         theta = list(
                             hyperid =  80001,
                             name = "log alpha",
                             short.name = "alpha",
                             initial = 1,
                             fixed = FALSE,
                             prior = "loggamma",
                             param = c(25, 25),
                             to.theta = function(x) log(x),
                             from.theta = function(x) exp(x)
                             )
                         ),
                     survival = FALSE,
                     discrete = FALSE,
                     link = c("default", "log", "neglog"),
                     pdf = "loglogistic"
                     ),

                 loglogisticsurv = list(
                     doc = "The loglogistic likelihood (survival)", 
                     hyper = list(
                         theta = list(
                             hyperid =  80011,
                             name = "log alpha",
                             short.name = "alpha",
                             initial = 1,
                             fixed = FALSE,
                             prior = "loggamma",
                             param = c(25, 25),
                             to.theta = function(x) log(x),
                             from.theta = function(x) exp(x)
                             )
                         ),
                     survival = TRUE,
                     discrete = FALSE,
                     link = c("default", "log", "neglog"),
                     pdf = "loglogistic"
                     ),

                 weibullcure = list(
                     doc = "The Weibull-cure likelihood (survival)", 
                     hyper = list(
                         theta1 = list(
                             hyperid =  81001,
                             name = "log alpha",
                             short.name = "a",
                             initial = 0.1,
                             fixed = FALSE,
                             prior = "pc.alphaw",
                             param = c(5),
                             ## the 'sc' constant is defined in inla.h, and must be the same.
                             ## I know, this is hard-coded for the moment. Should be a generic
                             ## way of doing this...
                             to.theta = function(x, sc = 0.1) log(x)/sc,
                             from.theta = function(x, sc = 0.1) exp(sc*x)
                             ),
                         theta2 = list(
                             hyperid =  81002,
                             name = "logit probability",
                             short.name = "prob",
                             initial = -1,
                             fixed = FALSE,
                             prior = "gaussian",
                             param = c(-1, 0.2),
                             to.theta = function(x) log(x/(1-x)),
                             from.theta = function(x) exp(x)/(1+exp(x))
                             )
                         ),
                     survival = TRUE,
                     discrete = FALSE,
                     link = c("default", "log", "neglog"),
                     pdf = "weibullcure"
                     ),

                 stochvol = list(
                     doc = "The Gaussian stochvol likelihood", 
                     hyper = list(
                         theta = list(
                             hyperid =  82001,
                             name = "log precision",
                             short.name = "prec",
                             initial = 500, ## yes, this is correct
                             fixed = TRUE, ## yes, this is correct
                             prior = "loggamma",
                             param = c(1, 0.005),
                             to.theta = function(x) log(x),
                             from.theta = function(x) exp(x)
                             )
                         ),
                     survival = FALSE,
                     discrete = FALSE,
                     link = c("default", "log"),
                     pdf = "stochvolgaussian"
                     ),

                 stochvolt = list(
                     doc = "The Student-t stochvol likelihood", 
                     hyper = list(
                         theta = list(
                             hyperid =  83001,
                             name = "log degrees of freedom",
                             short.name = "dof",
                             initial = 4,
                             fixed = FALSE,
                             prior = "pc.dof",
                             param = c(15, 0.5),
                             to.theta = function(x) log(x-2),
                             from.theta = function(x) 2+exp(x)
                             )
                         ),
                     survival = FALSE,
                     discrete = FALSE,
                     link = c("default", "log"),
                     pdf = "stochvolt"
                     ),

                 stochvolnig = list(
                     doc = "The Normal inverse Gaussian stochvol likelihood", 
                     hyper = list(
                         theta1 = list(
                             hyperid =  84001,
                             name = "skewness",
                             short.name = "skew",
                             initial = 0,
                             fixed = FALSE,
                             prior = "gaussian",
                             param = c(0, 10),
                             to.theta = function(x) x,
                             from.theta = function(x) x
                             ),
                         theta2 = list(
                             hyperid =  84002,
                             name = "shape",
                             short.name = "shape",
                             initial = 0,
                             fixed = FALSE,
                             prior = "loggamma",
                             param = c(1, 0.5),
                             to.theta = function(x) log(x-1),
                             from.theta = function(x) 1+exp(x)
                             )
                         ),
                     survival = FALSE,
                     discrete = FALSE,
                     link = c("default", "log"),
                     pdf = "stochvolnig"
                     ),

                 zeroinflatedpoisson0 = list(
                     doc = "Zero-inflated Poisson, type 0", 
                     hyper = list(
                         theta = list(
                             hyperid =  85001,
                             name = "logit probability",
                             short.name = "prob",
                             initial = -1,
                             fixed = FALSE,
                             prior = "gaussian",
                             param = c(-1, 0.2),
                             to.theta = function(x) log(x/(1-x)),
                             from.theta = function(x) exp(x)/(1+exp(x))
                             )
                         ),
                     survival = FALSE,
                     discrete = FALSE,
                     link = c("default", "log"),
                     pdf = "zeroinflated"
                     ),

                 zeroinflatedpoisson1 = list(
                     doc = "Zero-inflated Poisson, type 1", 
                     hyper = list(
                         theta = list(
                             hyperid =  86001,
                             name = "logit probability",
                             short.name = "prob",
                             initial = -1,
                             fixed = FALSE,
                             prior = "gaussian",
                             param = c(-1, 0.2),
                             to.theta = function(x) log(x/(1-x)),
                             from.theta = function(x) exp(x)/(1+exp(x))
                             )
                         ),
                     survival = FALSE,
                     discrete = FALSE,
                     link = c("default", "log"),
                     pdf = "zeroinflated"
                     ),

                 zeroinflatedpoisson2 = list(
                     doc = "Zero-inflated Poisson, type 2", 
                     hyper = list(
                         theta = list(
                             hyperid =  87001,
                             name = "log alpha",
                             short.name = "a",
                             initial = log(2),
                             fixed = FALSE,
                             prior = "gaussian",
                             param = c(log(2), 1),
                             to.theta = function(x) log(x),
                             from.theta = function(x) exp(x)
                             )
                         ),
                     survival = FALSE,
                     discrete = FALSE,
                     link = c("default", "log"),
                     pdf = "zeroinflated"
                     ),

                 zeroinflatedbetabinomial0 = list(
                     doc = "Zero-inflated Beta-Binomial, type 0", 
                     hyper = list(
                         theta1 = list(
                             hyperid =  88001,
                             name = "overdispersion",
                             short.name = "rho",
                             initial = 0,
                             fixed = FALSE,
                             prior = "gaussian",
                             param = c(0.0, 0.4),
                             to.theta = function(x) log(x/(1-x)),
                             from.theta = function(x) exp(x)/(1+exp(x))
                             ), 
                         theta2 = list(
                             hyperid =  88002,
                             name = "logit probability",
                             short.name = "prob",
                             initial = -1,
                             fixed = FALSE,
                             prior = "gaussian",
                             param = c(-1, 0.2),
                             to.theta = function(x) log(x/(1-x)),
                             from.theta = function(x) exp(x)/(1+exp(x))
                             )
                         ),
                     survival = FALSE,
                     discrete = TRUE,
                     link = c("default", "logit", "cauchit", "probit", "cloglog", "loglog", "robit", "sn"),
                     pdf = "zeroinflated"
                     ),

                 zeroinflatedbetabinomial1 = list(
                     doc = "Zero-inflated Beta-Binomial, type 1", 
                     hyper = list(
                         theta1 = list(
                             hyperid =  89001,
                             name = "overdispersion",
                             short.name = "rho",
                             initial = 0,
                             fixed = FALSE,
                             prior = "gaussian",
                             param = c(0.0, 0.4),
                             to.theta = function(x) log(x/(1-x)),
                             from.theta = function(x) exp(x)/(1+exp(x))
                             ), 
                         theta2 = list(
                             hyperid =  89002,
                             name = "logit probability",
                             short.name = "prob",
                             initial = -1,
                             fixed = FALSE,
                             prior = "gaussian",
                             param = c(-1, 0.2),
                             to.theta = function(x) log(x/(1-x)),
                             from.theta = function(x) exp(x)/(1+exp(x))
                             )
                         ),
                     survival = FALSE,
                     discrete = TRUE,
                     link = c("default", "logit", "cauchit", "probit", "cloglog", "loglog", "robit", "sn"),
                     pdf = "zeroinflated"
                     ),

                 zeroinflatedbinomial0 = list(
                     doc = "Zero-inflated Binomial, type 0", 
                     hyper = list(
                         theta = list(
                             hyperid =  90001,
                             name = "logit probability",
                             short.name = "prob",
                             initial = -1,
                             fixed = FALSE,
                             prior = "gaussian",
                             param = c(-1, 0.2),
                             to.theta = function(x) log(x/(1-x)),
                             from.theta = function(x) exp(x)/(1+exp(x))
                             )
                         ),
                     survival = FALSE,
                     discrete = FALSE,
                     link = c("default", "logit", "cauchit", "probit", "cloglog", "loglog", "robit", "sn"),
                     pdf = "zeroinflated"
                     ),

                 zeroinflatedbinomial1 = list(
                     doc = "Zero-inflated Binomial, type 1", 
                     hyper = list(
                         theta = list(
                             hyperid =  91001,
                             name = "logit probability",
                             short.name = "prob",
                             initial = -1,
                             fixed = FALSE,
                             prior = "gaussian",
                             param = c(-1, 0.2),
                             to.theta = function(x) log(x/(1-x)),
                             from.theta = function(x) exp(x)/(1+exp(x))
                             )
                         ),
                     survival = FALSE,
                     discrete = FALSE,
                     link = c("default", "logit", "cauchit", "probit", "cloglog", "loglog", "robit", "sn"),
                     pdf = "zeroinflated"
                     ),

                 zeroinflatedbinomial2 = list(
                     doc = "Zero-inflated Binomial, type 2", 
                     hyper = list(
                         theta = list(
                             hyperid =  92001,
                             name = "alpha",
                             short.name = "alpha",
                             initial = -1,
                             fixed = FALSE,
                             prior = "gaussian",
                             param = c(-1, 0.2),
                             to.theta = function(x) log(x),
                             from.theta = function(x) exp(x)
                             )
                         ),
                     survival = FALSE,
                     discrete = FALSE,
                     link = c("default", "logit", "cauchit", "probit", "cloglog", "loglog", "robit", "sn"),
                     pdf = "zeroinflated"
                     ),

                 zeroninflatedbinomial2 = list(
                     doc = "Zero and N inflated binomial, type 2", 
                     hyper = list(
                         theta1 = list(
                             hyperid =  93001,
                             name = "alpha1",
                             short.name = "alpha1",
                             initial = -1,
                             fixed = FALSE,
                             prior = "gaussian",
                             param = c(-1, 0.2),
                             to.theta = function(x) log(x),
                             from.theta = function(x) exp(x)
                             ),
                         theta2 = list(
                             hyperid =  93002,
                             name = "alpha2",
                             short.name = "alpha2",
                             initial = -1,
                             fixed = FALSE,
                             prior = "gaussian",
                             param = c(-1, 0.2),
                             to.theta = function(x) log(x),
                             from.theta = function(x) exp(x)
                             )
                         ),
                     survival = FALSE,
                     discrete = FALSE,
                     link = c("default", "logit", "cauchit", "probit", "cloglog", "loglog", "robit", "sn"),
                     pdf = NA
                     ),

                 zeroninflatedbinomial3 = list(
                     doc = "Zero and N inflated binomial, type 3", 
                     hyper = list(
                         theta1 = list(
                             hyperid =  93101,
                             name = "alpha0",
                             short.name = "alpha0",
                             initial = 1,
                             fixed = FALSE,
                             prior = "loggamma",
                             param = c(1, 1),
                             to.theta = function(x) log(x),
                             from.theta = function(x) exp(x)
                             ),
                         theta2 = list(
                             hyperid =  93102,
                             name = "alphaN",
                             short.name = "alphaN",
                             initial = 1,
                             fixed = FALSE,
                             prior = "loggamma",
                             param = c(1, 1),
                             to.theta = function(x) log(x),
                             from.theta = function(x) exp(x)
                             )
                         ),
                     status = "experimental", 
                     survival = FALSE,
                     discrete = FALSE,
                     link = c("default", "logit", "cauchit", "probit", "cloglog", "loglog", "robit", "sn"),
                     pdf = "zeroinflated"
                     ),

                 zeroinflatedbetabinomial2 = list(
                     doc = "Zero inflated Beta-Binomial, type 2", 
                     hyper = list(
                         theta1 = list(
                             hyperid =  94001,
                             name = "log alpha",
                             short.name = "a",
                             initial = log(2),
                             fixed = FALSE,
                             prior = "gaussian",
                             param = c(log(2), 1),
                             to.theta = function(x) log(x),
                             from.theta = function(x) exp(x)
                             ),
                         theta2 = list(
                             hyperid =  94002,
                             name = "beta",
                             short.name = "b",
                             initial = log(1),
                             fixed = FALSE,
                             prior = "gaussian",
                             param = c(0, 1),
                             to.theta = function(x) log(x),
                             from.theta = function(x) exp(x)
                             )
                         ),
                     survival = FALSE,
                     discrete = FALSE,
                     link = c("default", "logit", "cauchit", "probit", "cloglog", "loglog", "robit", "sn"),
                     pdf = "zeroinflated"
                     ),

                 zeroinflatednbinomial0 = list(
                     doc = "Zero inflated negBinomial, type 0", 
                     hyper = list(
                         theta1 = list(
                             hyperid =  95001,
                             name = "log size",
                             short.name = "size",
                             initial = log(10),
                             fixed = FALSE,
                             prior = "pc.mgamma",
                             param = 7,
                             to.theta = function(x) log(x),
                             from.theta = function(x) exp(x)
                             ),
                         theta2 = list(
                             hyperid =  95002,
                             name = "logit probability",
                             short.name = "prob",
                             initial = -1,
                             fixed = FALSE,
                             prior = "gaussian",
                             param = c(-1, 0.2),
                             to.theta = function(x) log(x/(1-x)),
                             from.theta = function(x) exp(x)/(1+exp(x))
                             )
                         ),
                     survival = FALSE,
                     discrete = FALSE,
                     link = c("default", "log"),
                     pdf = "zeroinflated"
                     ),

                 zeroinflatednbinomial1 = list(
                     doc = "Zero inflated negBinomial, type 1", 
                     hyper = list(
                         theta1 = list(
                             hyperid =  96001,
                             name = "log size",
                             short.name = "size",
                             initial = log(10),
                             fixed = FALSE,
                             prior = "pc.mgamma",
                             param = 7,
                             to.theta = function(x) log(x),
                             from.theta = function(x) exp(x)
                             ),
                         theta2 = list(
                             hyperid =  96002,
                             name = "logit probability",
                             short.name = "prob",
                             initial = -1,
                             fixed = FALSE,
                             prior = "gaussian",
                             param = c(-1, 0.2),
                             to.theta = function(x) log(x/(1-x)),
                             from.theta = function(x) exp(x)/(1+exp(x))
                             )
                         ),
                     survival = FALSE,
                     discrete = FALSE,
                     link = c("default", "log"),
                     pdf = "zeroinflated"
                     ),

                 zeroinflatednbinomial1strata2 = list(
                     doc = "Zero inflated negBinomial, type 1, strata 2", 
                     hyper = list(
                         theta1 = list(
                             hyperid =  97001,
                             name = "log size",
                             short.name = "size",
                             initial = log(10),
                             fixed = FALSE,
                             prior = "pc.mgamma",
                             param = 7,
                             to.theta = function(x) log(x),
                             from.theta = function(x) exp(x)
                         ),
                         theta2 = list(
                             hyperid =  97002,
                             name = "logit probability 1",
                             short.name = "prob1",
                             initial = -1,
                             fixed = FALSE,
                             prior = "gaussian",
                             param = c(-1, 0.2),
                             to.theta = function(x) log(x/(1-x)),
                             from.theta = function(x) exp(x)/(1+exp(x))
                         ), 
                         theta3 = list(
                             hyperid =  97003,
                             name = "logit probability 2",
                             short.name = "prob2",
                             initial = -1,
                             fixed = FALSE,
                             prior = "gaussian",
                             param = c(-1, 0.2),
                             to.theta = function(x) log(x/(1-x)),
                             from.theta = function(x) exp(x)/(1+exp(x))
                         ), 
                         theta4 = list(
                             hyperid =  97004,
                             name = "logit probability 3",
                             short.name = "prob3",
                             initial = -1,
                             fixed = TRUE,
                             prior = "gaussian",
                             param = c(-1, 0.2),
                             to.theta = function(x) log(x/(1-x)),
                             from.theta = function(x) exp(x)/(1+exp(x))
                         ), 
                         theta5 = list(
                             hyperid =  97005,
                             name = "logit probability 4",
                             short.name = "prob4",
                             initial = -1,
                             fixed = TRUE,
                             prior = "gaussian",
                             param = c(-1, 0.2),
                             to.theta = function(x) log(x/(1-x)),
                             from.theta = function(x) exp(x)/(1+exp(x))
                         ), 
                         theta6 = list(
                             hyperid =  97006,
                             name = "logit probability 5",
                             short.name = "prob5",
                             initial = -1,
                             fixed = TRUE,
                             prior = "gaussian",
                             param = c(-1, 0.2),
                             to.theta = function(x) log(x/(1-x)),
                             from.theta = function(x) exp(x)/(1+exp(x))
                         ), 
                         theta7 = list(
                             hyperid =  97007,
                             name = "logit probability 6",
                             short.name = "prob6",
                             initial = -1,
                             fixed = TRUE,
                             prior = "gaussian",
                             param = c(-1, 0.2),
                             to.theta = function(x) log(x/(1-x)),
                             from.theta = function(x) exp(x)/(1+exp(x))
                         ), 
                         theta8 = list(
                             hyperid =  97008,
                             name = "logit probability 7",
                             short.name = "prob7",
                             initial = -1,
                             fixed = TRUE,
                             prior = "gaussian",
                             param = c(-1, 0.2),
                             to.theta = function(x) log(x/(1-x)),
                             from.theta = function(x) exp(x)/(1+exp(x))
                         ), 
                         theta9 = list(
                             hyperid =  97009,
                             name = "logit probability 8",
                             short.name = "prob8",
                             initial = -1,
                             fixed = TRUE,
                             prior = "gaussian",
                             param = c(-1, 0.2),
                             to.theta = function(x) log(x/(1-x)),
                             from.theta = function(x) exp(x)/(1+exp(x))
                         ), 
                         theta10 = list(
                             hyperid =  97010,
                             name = "logit probability 9",
                             short.name = "prob9",
                             initial = -1,
                             fixed = TRUE,
                             prior = "gaussian",
                             param = c(-1, 0.2),
                             to.theta = function(x) log(x/(1-x)),
                             from.theta = function(x) exp(x)/(1+exp(x))
                         ), 
                         theta11 = list(
                             hyperid =  97011,
                             name = "logit probability 10",
                             short.name = "prob10",
                             initial = -1,
                             fixed = TRUE,
                             prior = "gaussian",
                             param = c(-1, 0.2),
                             to.theta = function(x) log(x/(1-x)),
                             from.theta = function(x) exp(x)/(1+exp(x))
                         )
                     ), 
                     status = "experimental", 
                     survival = FALSE,
                     discrete = FALSE,
                     link = c("default", "log"),
                     pdf = "zeroinflated"
                     ),

                 zeroinflatednbinomial1strata3 = list(
                     doc = "Zero inflated negBinomial, type 1, strata 3", 
                     hyper = list(
                         theta1 = list(
                             hyperid =  98001,
                             name = "logit probability",
                             short.name = "prob",
                             initial = -1,
                             fixed = FALSE,
                             prior = "gaussian",
                             param = c(-1, 0.2),
                             to.theta = function(x) log(x/(1-x)),
                             from.theta = function(x) exp(x)/(1+exp(x))
                         ), 
                         theta2 = list(
                             hyperid =  98002,
                             name = "log size 1",
                             short.name = "size1",
                             initial = log(10),
                             fixed = FALSE,
                             prior = "pc.mgamma",
                             param = 7,
                             to.theta = function(x) log(x),
                             from.theta = function(x) exp(x)
                         ),
                         theta3 = list(
                             hyperid =  98003,
                             name = "log size 2",
                             short.name = "size2",
                             initial = log(10),
                             fixed = FALSE,
                             prior = "pc.mgamma",
                             param = 7,
                             to.theta = function(x) log(x),
                             from.theta = function(x) exp(x)
                         ),
                         theta4 = list(
                             hyperid =  98004,
                             name = "log size 3",
                             short.name = "size3",
                             initial = log(10),
                             fixed = TRUE,
                             prior = "pc.mgamma",
                             param = 7,
                             to.theta = function(x) log(x),
                             from.theta = function(x) exp(x)
                         ),
                         theta5 = list(
                             hyperid =  98005, 
                             name = "log size 4",
                             short.name = "size4",
                             initial = log(10),
                             fixed = TRUE,
                             prior = "pc.mgamma",
                             param = 7,
                             to.theta = function(x) log(x),
                             from.theta = function(x) exp(x)
                         ),
                         theta6 = list(
                             hyperid =  98006,
                             name = "log size 5",
                             short.name = "size5",
                             initial = log(10),
                             fixed = TRUE,
                             prior = "pc.mgamma",
                             param = 7,
                             to.theta = function(x) log(x),
                             from.theta = function(x) exp(x)
                         ),
                         theta7 = list(
                             hyperid =  98007,
                             name = "log size 6",
                             short.name = "size6",
                             initial = log(10),
                             fixed = TRUE,
                             prior = "pc.mgamma",
                             param = 7,
                             to.theta = function(x) log(x),
                             from.theta = function(x) exp(x)
                         ),
                         theta8 = list(
                             hyperid =  98008,
                             name = "log size 7",
                             short.name = "size7",
                             initial = log(10),
                             fixed = TRUE,
                             prior = "pc.mgamma",
                             param = 7,
                             to.theta = function(x) log(x),
                             from.theta = function(x) exp(x)
                         ),
                         theta9 = list(
                             hyperid =  98009,
                             name = "log size 8",
                             short.name = "size8",
                             initial = log(10),
                             fixed = TRUE,
                             prior = "pc.mgamma",
                             param = 7,
                             to.theta = function(x) log(x),
                             from.theta = function(x) exp(x)
                         ),
                         theta10 = list(
                             hyperid =  98010,
                             name = "log size 9",
                             short.name = "size9",
                             initial = log(10),
                             fixed = TRUE,
                             prior = "pc.mgamma",
                             param = 7,
                             to.theta = function(x) log(x),
                             from.theta = function(x) exp(x)
                         ),
                         theta11 = list(
                             hyperid =  98011,
                             name = "log size 10",
                             short.name = "size10",
                             initial = log(10),
                             fixed = TRUE,
                             prior = "pc.mgamma",
                             param = 7,
                             to.theta = function(x) log(x),
                             from.theta = function(x) exp(x)
                         )
                     ),
                     status = "experimental", 
                     survival = FALSE,
                     discrete = FALSE,
                     link = c("default", "log"),
                     pdf = "zeroinflated"
                     ),

                 zeroinflatednbinomial2 = list(
                     doc = "Zero inflated negBinomial, type 2", 
                     hyper = list(
                         theta1 = list(
                             hyperid =  99001,
                             name = "log size",
                             short.name = "size",
                             initial = log(10),
                             fixed = FALSE,
                             prior = "pc.mgamma",
                             param = 7,
                             to.theta = function(x) log(x),
                             from.theta = function(x) exp(x)
                             ),
                         theta2 = list(
                             hyperid =  99002,
                             name = "log alpha",
                             short.name = "a",
                             initial = log(2),
                             fixed = FALSE,
                             prior = "gaussian",
                             param = c(2, 1),
                             to.theta = function(x) log(x),
                             from.theta = function(x) exp(x)
                             )
                         ),
                     survival = FALSE,
                     discrete = FALSE,
                     link = c("default", "log"),
                     pdf = "zeroinflated"
                     ),

                 t = list(
                     doc = "Student-t likelihood", 
                     hyper = list(
                         theta1 = list(
                             hyperid =  100001,
                             name = "log precision",
                             short.name = "prec",
                             initial = 0,
                             fixed = FALSE,
                             prior = "loggamma",
                             param = c(1, 0.00005),
                             to.theta = function(x) log(x),
                             from.theta = function(x) exp(x)
                             ),
                         theta2 = list(
                             hyperid =  100002,
                             name = "log degrees of freedom",
                             short.name = "dof",
                             initial = 5,
                             fixed = FALSE,
                             prior = "pc.dof",
                             param = c(15, 0.5),
                             to.theta = function(x) log(x-2),
                             from.theta = function(x) 2+exp(x)
                             )
                         ),
                     survival = FALSE,
                     discrete = FALSE,
                     link = c("default", "identity"),
                     pdf = "student-t"
                     ),

                 tstrata = list(
                     doc = "A stratified version of the Student-t likelihood", 
                     hyper = list(
                         theta1 = list(
                             hyperid =  101001,
                             name = "log degrees of freedom",
                             short.name = "dof",
                             initial = 4,
                             fixed = FALSE,
                             prior = "pc.dof",
                             param = c(15, 0.5),
                             to.theta = function(x) log(x-5),
                             from.theta = function(x) 5+exp(x)
                             ),
                         theta2 = list(
                             hyperid =  101002,
                             name = "log precision1",
                             short.name = "prec1",
                             initial = 2,
                             fixed = FALSE,
                             prior = "loggamma",
                             param = c(1, 0.00005),
                             to.theta = function(x) log(x),
                             from.theta = function(x) exp(x)
                             ),
                         theta3 = list(
                             hyperid =  101003,
                             name = "log precision2",
                             short.name = "prec2",
                             initial = 2,
                             fixed = FALSE,
                             prior = "loggamma",
                             param = c(1, 0.00005),
                             to.theta = function(x) log(x),
                             from.theta = function(x) exp(x)
                             ),
                         theta4 = list(
                             hyperid =  101004,
                             name = "log precision3",
                             short.name = "prec3",
                             initial = 2,
                             fixed = FALSE,
                             prior = "loggamma",
                             param = c(1, 0.00005),
                             to.theta = function(x) log(x),
                             from.theta = function(x) exp(x)
                             ),
                         theta5 = list(
                             hyperid =  101005,
                             name = "log precision4",
                             short.name = "prec4",
                             initial = 2,
                             fixed = FALSE,
                             prior = "loggamma",
                             param = c(1, 0.00005),
                             to.theta = function(x) log(x),
                             from.theta = function(x) exp(x)
                             ),
                         theta6 = list(
                             hyperid =  101006,
                             name = "log precision5",
                             short.name = "prec5",
                             initial = 2,
                             fixed = FALSE,
                             prior = "loggamma",
                             param = c(1, 0.00005),
                             to.theta = function(x) log(x),
                             from.theta = function(x) exp(x)
                             ),
                         theta7 = list(
                             hyperid =  101007,
                             name = "log precision6",
                             short.name = "prec6",
                             initial = 2,
                             fixed = FALSE,
                             prior = "loggamma",
                             param = c(1, 0.00005),
                             to.theta = function(x) log(x),
                             from.theta = function(x) exp(x)
                             ),
                         theta8 = list(
                             hyperid =  101008,
                             name = "log precision7",
                             short.name = "prec7",
                             initial = 2,
                             fixed = FALSE,
                             prior = "loggamma",
                             param = c(1, 0.00005),
                             to.theta = function(x) log(x),
                             from.theta = function(x) exp(x)
                             ),
                         theta9 = list(
                             hyperid =  101009,
                             name = "log precision8",
                             short.name = "prec8",
                             initial = 2,
                             fixed = FALSE,
                             prior = "loggamma",
                             param = c(1, 0.00005),
                             to.theta = function(x) log(x),
                             from.theta = function(x) exp(x)
                             ),
                         theta10 = list(
                             hyperid =  101010,
                             name = "log precision9",
                             short.name = "prec9",
                             initial = 2,
                             fixed = FALSE,
                             prior = "loggamma",
                             param = c(1, 0.00005),
                             to.theta = function(x) log(x),
                             from.theta = function(x) exp(x)
                             ),
                         theta11 = list(
                             hyperid =  101011,
                             name = "log precision10",
                             short.name = "prec10",
                             initial = 2,
                             fixed = FALSE,
                             prior = "loggamma",
                             param = c(1, 0.00005),
                             to.theta = function(x) log(x),
                             from.theta = function(x) exp(x)
                             )
                         ),
                     survival = FALSE,
                     discrete = FALSE,
                     link = c("default", "identity"),
                     pdf = "tstrata"
                     ),

                 nmix = list(
                     doc = "Binomial-Poisson mixture", 
                     hyper = list(
                         theta1 = list(
                             hyperid =  101101,
                             name = "beta1",
                             short.name = "beta1",
                             initial = log(10),
                             fixed = FALSE,
                             prior = "normal",
                             param = c(0, 0.5),
                             to.theta = function(x) x, 
                             from.theta = function(x) x
                         ),
                         theta2 = list(
                             hyperid =  101102,
                             name = "beta2",
                             short.name = "beta2",
                             initial = 0,
                             fixed = FALSE,
                             prior = "normal",
                             param = c(0, 1),
                             to.theta = function(x) x, 
                             from.theta = function(x) x
                         ),
                         theta3 = list(
                             hyperid =  101103,
                             name = "beta3",
                             short.name = "beta3",
                             initial = 0,
                             fixed = FALSE,
                             prior = "normal",
                             param = c(0, 1),
                             to.theta = function(x) x, 
                             from.theta = function(x) x
                         ),
                         theta4 = list(
                             hyperid =  101104,
                             name = "beta4",
                             short.name = "beta4",
                             initial = 0,
                             fixed = FALSE,
                             prior = "normal",
                             param = c(0, 1),
                             to.theta = function(x) x, 
                             from.theta = function(x) x
                         ),
                         theta5 = list(
                             hyperid =  101105,
                             name = "beta5",
                             short.name = "beta5",
                             initial = 0,
                             fixed = FALSE,
                             prior = "normal",
                             param = c(0, 1),
                             to.theta = function(x) x, 
                             from.theta = function(x) x
                         ), 
                         theta6 = list(
                             hyperid =  101106,
                             name = "beta6",
                             short.name = "beta6",
                             initial = 0,
                             fixed = FALSE,
                             prior = "normal",
                             param = c(0, 1),
                             to.theta = function(x) x, 
                             from.theta = function(x) x
                         ), 
                         theta7 = list(
                             hyperid =  101107,
                             name = "beta7",
                             short.name = "beta7",
                             initial = 0,
                             fixed = FALSE,
                             prior = "normal",
                             param = c(0, 1),
                             to.theta = function(x) x, 
                             from.theta = function(x) x
                         ),
                         theta8 = list(
                             hyperid =  101108,
                             name = "beta8",
                             short.name = "beta8",
                             initial = 0,
                             fixed = FALSE,
                             prior = "normal",
                             param = c(0, 1),
                             to.theta = function(x) x, 
                             from.theta = function(x) x
                         ), 
                         theta9 = list(
                             hyperid =  101109,
                             name = "beta9",
                             short.name = "beta9",
                             initial = 0,
                             fixed = FALSE,
                             prior = "normal",
                             param = c(0, 1),
                             to.theta = function(x) x, 
                             from.theta = function(x) x
                         ),
                         theta10 = list(
                             hyperid =  101110,
                             name = "beta10",
                             short.name = "beta10",
                             initial = 0,
                             fixed = FALSE,
                             prior = "normal",
                             param = c(0, 1),
                             to.theta = function(x) x, 
                             from.theta = function(x) x
                         ), 
                         theta11 = list(
                             hyperid =  101111,
                             name = "beta11",
                             short.name = "beta11",
                             initial = 0,
                             fixed = FALSE,
                             prior = "normal",
                             param = c(0, 1),
                             to.theta = function(x) x, 
                             from.theta = function(x) x
                         ), 
                         theta12 = list(
                             hyperid =  101112,
                             name = "beta12",
                             short.name = "beta12",
                             initial = 0,
                             fixed = FALSE,
                             prior = "normal",
                             param = c(0, 1),
                             to.theta = function(x) x, 
                             from.theta = function(x) x
                         ), 
                         theta13 = list(
                             hyperid =  101113,
                             name = "beta13",
                             short.name = "beta13",
                             initial = 0,
                             fixed = FALSE,
                             prior = "normal",
                             param = c(0, 1),
                             to.theta = function(x) x, 
                             from.theta = function(x) x
                         ), 
                         theta14 = list(
                             hyperid =  101114,
                             name = "beta14",
                             short.name = "beta14",
                             initial = 0,
                             fixed = FALSE,
                             prior = "normal",
                             param = c(0, 1),
                             to.theta = function(x) x, 
                             from.theta = function(x) x
                         ), 
                         theta15 = list(
                             hyperid =  101115,
                             name = "beta15",
                             short.name = "beta15",
                             initial = 0,
                             fixed = FALSE,
                             prior = "normal",
                             param = c(0, 1),
                             to.theta = function(x) x, 
                             from.theta = function(x) x
                         )
                     ),
                     status = "experimental", 
                     survival = FALSE,
                     discrete = TRUE,
                     link = c("default", "logit", "probit"),
                     pdf = "nmix"
                 ),

                 nmixnb = list(
                     doc = "NegBinomial-Poisson mixture", 
                     hyper = list(
                         theta1 = list(
                             hyperid =  101121,
                             name = "beta1",
                             short.name = "beta1",
                             initial = log(10),
                             fixed = FALSE,
                             prior = "normal",
                             param = c(0, 0.5),
                             to.theta = function(x) x, 
                             from.theta = function(x) x
                         ),
                         theta2 = list(
                             hyperid =  101122,
                             name = "beta2",
                             short.name = "beta2",
                             initial = 0,
                             fixed = FALSE,
                             prior = "normal",
                             param = c(0, 1),
                             to.theta = function(x) x, 
                             from.theta = function(x) x
                         ),
                         theta3 = list(
                             hyperid =  101123,
                             name = "beta3",
                             short.name = "beta3",
                             initial = 0,
                             fixed = FALSE,
                             prior = "normal",
                             param = c(0, 1),
                             to.theta = function(x) x, 
                             from.theta = function(x) x
                         ),
                         theta4 = list(
                             hyperid =  101124,
                             name = "beta4",
                             short.name = "beta4",
                             initial = 0,
                             fixed = FALSE,
                             prior = "normal",
                             param = c(0, 1),
                             to.theta = function(x) x, 
                             from.theta = function(x) x
                         ),
                         theta5 = list(
                             hyperid =  101125,
                             name = "beta5",
                             short.name = "beta5",
                             initial = 0,
                             fixed = FALSE,
                             prior = "normal",
                             param = c(0, 1),
                             to.theta = function(x) x, 
                             from.theta = function(x) x
                         ), 
                         theta6 = list(
                             hyperid =  101126,
                             name = "beta6",
                             short.name = "beta6",
                             initial = 0,
                             fixed = FALSE,
                             prior = "normal",
                             param = c(0, 1),
                             to.theta = function(x) x, 
                             from.theta = function(x) x
                         ), 
                         theta7 = list(
                             hyperid =  101127,
                             name = "beta7",
                             short.name = "beta7",
                             initial = 0,
                             fixed = FALSE,
                             prior = "normal",
                             param = c(0, 1),
                             to.theta = function(x) x, 
                             from.theta = function(x) x
                         ), 
                         theta8 = list(
                             hyperid =  101128,
                             name = "beta8",
                             short.name = "beta8",
                             initial = 0,
                             fixed = FALSE,
                             prior = "normal",
                             param = c(0, 1),
                             to.theta = function(x) x, 
                             from.theta = function(x) x
                         ), 
                         theta9 = list(
                             hyperid =  101129,
                             name = "beta9",
                             short.name = "beta9",
                             initial = 0,
                             fixed = FALSE,
                             prior = "normal",
                             param = c(0, 1),
                             to.theta = function(x) x, 
                             from.theta = function(x) x
                         ), 
                         theta10 = list(
                             hyperid =  101130,
                             name = "beta10",
                             short.name = "beta10",
                             initial = 0,
                             fixed = FALSE,
                             prior = "normal",
                             param = c(0, 1),
                             to.theta = function(x) x, 
                             from.theta = function(x) x
                         ), 
                         theta11 = list(
                             hyperid =  101131,
                             name = "beta11",
                             short.name = "beta11",
                             initial = 0,
                             fixed = FALSE,
                             prior = "normal",
                             param = c(0, 1),
                             to.theta = function(x) x, 
                             from.theta = function(x) x
                         ), 
                         theta12 = list(
                             hyperid =  101132,
                             name = "beta12",
                             short.name = "beta12",
                             initial = 0,
                             fixed = FALSE,
                             prior = "normal",
                             param = c(0, 1),
                             to.theta = function(x) x, 
                             from.theta = function(x) x
                         ), 
                         theta13 = list(
                             hyperid =  101133,
                             name = "beta13",
                             short.name = "beta13",
                             initial = 0,
                             fixed = FALSE,
                             prior = "normal",
                             param = c(0, 1),
                             to.theta = function(x) x, 
                             from.theta = function(x) x
                         ), 
                         theta14 = list(
                             hyperid =  101134,
                             name = "beta14",
                             short.name = "beta14",
                             initial = 0,
                             fixed = FALSE,
                             prior = "normal",
                             param = c(0, 1),
                             to.theta = function(x) x, 
                             from.theta = function(x) x
                         ), 
                         theta15 = list(
                             hyperid =  101135,
                             name = "beta15",
                             short.name = "beta15",
                             initial = 0,
                             fixed = FALSE,
                             prior = "normal",
                             param = c(0, 1),
                             to.theta = function(x) x, 
                             from.theta = function(x) x
                         ), 
                         theta16 = list(
                             hyperid =  101136,
                             name = "overdispersion",
                             short.name = "overdispersion",
                             initial = 0,
                             fixed = FALSE,
                             prior = "pc.gamma",
                             param = 7,
                             to.theta = function(x) log(x), 
                             from.theta = function(x) exp(x)
                         )
                     ),
                     status = "experimental", 
                     survival = FALSE,
                     discrete = TRUE,
                     link = c("default", "logit", "probit"),
                     pdf = "nmixnb"
                 ),

                 gp = list(
                     doc = "Generalized Pareto likelihood", 
                     hyper = list(
                         theta = list(
                             hyperid =  101201,
                             name = "shape",
                             short.name = "xi",
                             initial = log(0.1),
                             fixed = FALSE,
                             prior = "loggamma",
                             param = c(1, 15), 
                             to.theta = function(x) log(x), 
                             from.theta = function(x) exp(x)
                         )
                     ),
                     status = "experimental", 
                     survival = FALSE,
                     discrete = TRUE,
                     link = c("default", "quantile"), 
                     pdf = "genPareto"
                 ),

                 logperiodogram = list(
                     doc = "Likelihood for the log-periodogram", 
                     hyper = list(
                     ),
                     survival = FALSE,
                     discrete = FALSE,
                     link = c("default", "identity"),
                     pdf = NA
                 )
             )
         )
}

`inla.models.section.prior` = function()
{
    return
    list(prior =
             list(
                 normal = list(
                     doc = "Normal prior", 
                     nparameters = 2L,
                     pdf = "gaussian"
                 ),
                 gaussian = list(
                     doc = "Gaussian prior", 
                     nparameters = 2L,
                     pdf = "gaussian"
                 ),
                 wishart1d = list(
                     doc = "Wishart prior dim=1", 
                     nparameters = 2L,
                     pdf = "iid123d"
                 ),
                 wishart2d = list(
                     doc = "Wishart prior dim=2", 
                     nparameters = 4L,
                     pdf = "iid123d"
                 ),
                 wishart3d = list(
                     doc = "Wishart prior dim=3", 
                     nparameters = 7L,
                     pdf = "iid123d"
                 ),
                 wishart4d = list(
                     doc = "Wishart prior dim=4", 
                     nparameters = 11L,
                     pdf = "iid123d"
                 ),
                 wishart5d = list(
                     doc = "Wishart prior dim=5", 
                     nparameters = 16L,
                     pdf = "iid123d"
                 ),
                 loggamma = list(
                     doc = "Log-Gamma prior", 
                     nparameters = 2L,
                     pdf = "prior-loggamma"
                 ),
                 gamma = list(
                     doc = "Gamma prior", 
                     nparameters = 2L,
                     pdf = "prior-loggamma"
                 ),
                 minuslogsqrtruncnormal = list(
                     doc = "(obsolete)", 
                     nparameters = 2L,
                     pdf = "prior-logtnorm"
                 ),
                 logtnormal = list(
                     doc = "Truncated Normal prior", 
                     nparameters = 2L,
                     pdf = "prior-logtnorm"
                 ),
                 logtgaussian = list(
                     doc = "Truncated Gaussian prior", 
                     nparameters = 2L,
                     pdf = "prior-logtnorm"
                 ),
                 flat=list(
                     doc = "A constant prior", 
                     nparameters = 0L,
                     pdf = "various-flat"
                 ),
                 logflat=list(
                     doc = "A constant prior for log(theta)", 
                     nparameters = 0L,
                     pdf = "various-flat"
                 ),
                 logiflat=list(
                     doc = "A constant prior for log(1/theta)", 
                     nparameters = 0L,
                     pdf = "various-flat"
                 ),
                 mvnorm = list(
                     doc = "A multivariate Normal prior", 
                     nparameters = -1L,
                     pdf = "mvnorm"
                 ),
                 pc.alphaw = list(
                     doc = "PC prior for alpha in Weibull",
                     nparameters = 1L,
                     pdf = "pc.alphaw"
                 ), 
                 pc.ar = list(
                     doc = "PC prior for the AR(p) model", 
                     nparameters = 1L,
                     pdf = "pc.ar"
                 ),
                 dirichlet = list(
                     doc ="Dirichlet prior",
                     nparameters = 1L,
                     pdf = "dirichlet"
                 ),

                 ## this is the 'no prior needed' prior
                 none = list(
                     doc = "No prior", 
                     nparameters = 0L),

                 ## this is the 'flag an error if used' prior
                 invalid = list(
                     doc = "Void prior", 
                     nparameters = 0L),

                 betacorrelation = list(
                     doc = "Beta prior for the correlation", 
                     nparameters = 2L,
                     pdf = "betacorrelation"
                 ),

                 logitbeta = list(
                     doc = "Logit prior for a probability", 
                     nparameters = 2L,
                     pdf = "logitbeta"
                 ),

                 pc.prec = list(
                     doc = "PC prior for log(precision)", 
                     nparameters = 2L,
                     pdf = "pc.prec"
                 ),
                 
                 pc.dof = list(
                     doc = "PC prior for log(dof-2)", 
                     nparameters = 2L,
                     pdf = "pc.dof"
                 ),
                 
                 pc.cor0 = list(
                     doc = "PC prior correlation, basemodel cor=0", 
                     nparameters = 2L,
                     pdf = "pc.cor0"
                 ),
                 
                 pc.cor1 = list(
                     doc = "PC prior correlation, basemodel cor=1", 
                     nparameters = 2L,
                     pdf = "pc.cor1"
                 ),
                 
                 pc.fgnh = list(
                     doc = "PC prior for the Hurst parameter in FGN", 
                     nparameters = 2L,
                     pdf = "pc.fgnh"
                 ), 
                 ## experimental prior from GA
                 pc.spde.GA = list(
                     doc = "(experimental)", 
                     nparameters = 4L,
                     pdf = NA
                 ), 
                 
		 pc.matern = list(
                     doc = "PC prior for the Matern SPDE", 
                     nparameters = 3L,
                     pdf = NA
                 ), 
                 
		 pc.range = list(
                     doc = "PC prior for the range in the Matern SPDE", 
                     nparameters = 2L,
                     pdf = NA
                 ), 

		 pc.sn = list(
                     doc = "PC prior for the skew-normal", 
                     nparameters = 1L,
                     pdf = "pc.sn"
                 ), 

		 pc.gamma = list(
                     doc = "PC prior for a Gamma parameter", 
                     nparameters = 1L,
                     pdf = "pc.gamma"
                 ), 

		 pc.mgamma = list(
                     doc = "PC prior for a Gamma parameter", 
                     nparameters = 1L,
                     pdf = "pc.gamma"
                 ), 

		 pc.gammacount = list(
                     doc = "PC prior for the GammaCount likelihood", 
                     nparameters = 1L,
                     pdf = "pc.gammacount"
                 ), 

		 pc.gevtail = list(
                     doc = "PC prior for the tail in the GEV likelihood", 
                     nparameters = 3L,
                     pdf = "pc.gevtail"
                 ), 

                 ## this is the generic one, which is case-spesific and possibly adaptive
                 pc = list(
                     doc = "Generic PC prior", 
                     nparameters = 2L,
                     pdf = NA
                 ), 

                 ref.ar = list(
                     doc = "Reference prior for the AR(p) model, p<=3", 
                     nparameters = 0L,
                     pdf = NA
                 ),
                 
                 pom = list(
                     doc = "#classes-dependent prior for the POM model", 
                     nparameters = 0L,
                     pdf = "pom"
                 ),
                 
                 jeffreystdf = list(
                     doc = "Jeffreys prior for the doc", 
                     nparameters = 0L,
                     pdf = "jeffreystdf"
                 ),

                 "expression:" = list(
                     doc = "A generic prior defined using expressions", 
                     nparameters = -1L,
                     pdf = "expression"
                 ), 

                 "table:" = list(
                     doc = "A generic tabulated prior", 
                     nparameters = -1L,
                     pdf = "table"
                 )
             )
         )
}

`inla.models.section.wrapper` = function()
{
    return
    list(wrapper =
             list(
                 joint = list(
                     doc = "(experimental)", 
                     hyper = list(
                         theta = list(
                             hyperid = 102001,
                             name = "log precision",
                             short.name = "prec",
                             initial = 0,
                             fixed = TRUE,
                             prior = "loggamma",
                             param = c(1, 0.00005),
                             to.theta = function(x) log(x),
                             from.theta = function(x) exp(x)
                             )
                         ),
                     constr = FALSE,
                     nrow.ncol = FALSE,
                     augmented = FALSE,
                     aug.factor = 1L,
                     aug.constr = NULL,
                     n.div.by = NULL,
                     n.required = FALSE,
                     set.default.values = FALSE,
                     pdf = NA
                     )
                 )
         )
}

`inla.models` = function()
{
    ## this is not very clean solution, but for the moment is ok. the
    ## inla.models() function takes just to much time!!!

    envir = inla.get.inlaEnv()
    
    if (exists("inla.models", envir = envir) &&
        exists("hgid", envir = envir) &&
        get("hgid", envir = envir) == inla.version("hgid")) {

        return (get("inla.models", envir = envir))

    } else {
        ## have to split it, as option keep.source has an upper limit...
        models = c(
            inla.models.section.latent(),
            inla.models.section.group(),
            inla.models.section.mix(),
            inla.models.section.link(),
            inla.models.section.predictor(),
            inla.models.section.hazard(),
            inla.models.section.likelihood(),
            inla.models.section.prior(),
            inla.models.section.wrapper())
        ## set "read.only" attribute for the `hyper' at those elements
        ## that cannot be changed.
        for (section in names(models)) {
            for (model in names(models[[section]])) {
                a = models[[section]][[model]]
                if (!is.null(a$hyper) && length(a$hyper) > 0) {
                    for (theta in names(a$hyper)) {
                        for (elm in names(a$hyper[[theta]])) {
                            val = FALSE
                            if (elm == "prior") {
                                if (a$hyper[[theta]][[elm]] == "none" ||
                                    a$hyper[[theta]][[elm]] == "wishart1d" ||
                                    a$hyper[[theta]][[elm]] == "wishart2d" ||
                                    a$hyper[[theta]][[elm]] == "wishart3d" ||
                                    a$hyper[[theta]][[elm]] == "wishart4d" ||
                                    a$hyper[[theta]][[elm]] == "wishart5d") {
                                    val = TRUE
                                }
                            }
                            if (elm == "to.theta" || elm == "from.theta")
                                val = TRUE
                            ## append the new attribute, so the code
                            ## depends on if the object has an attribute
                            ## from before or not.
                            if (is.null(attributes(models[[section]][[model]]$hyper[[theta]][[elm]]))) {
                                attributes(models[[section]][[model]]$hyper[[theta]][[elm]]) =
                                    list("inla.read.only" = val)
                            } else {
                                attributes(models[[section]][[model]]$hyper[[theta]][[elm]]) =
                                    c("inla.read.only" = val,
                                      attributes(models[[section]][[model]]$hyper[[theta]][[elm]]))
                            }
                        }
                    }
                }
            }
        }

        assign("inla.models", models, envir = envir)
        assign("hgid", inla.version("hgid"), envir = envir)

        return (models)
    }

    stop("This should not happen")
}

`inla.is.model` = function(model, section = NULL, 
        stop.on.error = TRUE, ignore.case = FALSE)
{
    mm = inla.models()
    if (is.null(section)) {
        stop("No section given; please fix...")
    }
    section = match.arg(section, names(mm))
    models = names((mm[ names(mm) == section ])[[1]])

    if (is.character(model) && length(model) > 0) {
        if (ignore.case) {
            m = tolower(model)
            ms = tolower(models)
        } else {
            m = model
            ms = models
        }

        ret = c()
        for(i in 1L:length(m)) {

            if (is.element(m[i], ms)) {
                ret[i] = TRUE
            } else {
                if (stop.on.error) {
                    print("Valid models are:")
                    print(models)
                    stop(paste(c("\n\tUnknown name [", model[i], "]\n", "\tValid choices are: ",
                                 models), sep=" ", collapse=" "))
                } else {
                    ret[i] = FALSE
                }
            }
        }
    } else {
        if (stop.on.error)
            stop(paste("\n\tModel [", model, "] is not of type character()\n", sep=""))
    }
    return (ret)
}

`inla.model.properties` = function(
        model,
        section = NULL, 
        stop.on.error = TRUE,
        ignore.case = FALSE)
{
    if (is.null(section)) {
        stop("No section given; please fix...")
    }
    mm = inla.models()
    section = match.arg(section, names(mm))
    m = inla.model.properties.generic(inla.trim.family(model),
        (mm[names(mm) == section])[[1]],
        stop.on.error, ignore.case,
        ## would like to know the section for a possible warning/error if the status says so.
        section = section)

    if (is.null(m)) {
        return (NULL)
    }

    return (m)
}

`inla.model.properties.generic` = function(model, models, stop.on.error = TRUE, ignore.case = TRUE, section = "UKNOWN")
{
    ## argument 'section' is for 'status' only...
    
    m = ifelse(ignore.case, tolower(model), model)
    if (ignore.case) {
        ms = tolower(names(models))
    } else {
        ms = names(models)
    }
    ms = inla.trim.family(ms)
    k = grep(paste("^", m, "$", sep=""), ms)
    if (length(k) == 0L) {
        if (stop.on.error) {
            stop(paste("Name '", model, "' is unknown. Valid choices are: ", inla.paste(ms), ".", sep=""))
        }
        return (NULL)
    } else {
        ## if 'status' is set, then issue a warning/error depending on
        ## 'status'. do this the first time only if status is
        ## 'experimental'.
        status = models[[k]]$status
        if (is.null(status)) {
            ## do nothing; all ok.
        } else {
            status.core = strsplit(status, ":")[[1]][1]
            stopifnot(any(inla.strcasecmp(status.core, c("experimental", "disabled", "changed"))))
            envir = inla.get.inlaEnv()
            var = paste("processed.status.for.model.", model, ".in.section.", section, sep="")

            if (inla.strcasecmp(status.core, "experimental")) {
                if (!(exists(var, envir = envir) && get(var, envir = envir))) {
                    assign(var, TRUE, envir = envir)
                    msg = paste0("Model '", model, "' in section '", section, "' is marked as '", status, 
                        "'; changes may appear at any time.",
                        "\n  ",
                        "Use this model with extra care!!! Further warnings are disabled.")
                    warning(msg)
                } else {
                    ## the warning is already given; do nothing
                }
            } else if (inla.strcasecmp(status.core, "disabled")) {
                assign(var, TRUE, envir = envir)
                msg = paste("Model '", model, "' in section '", section, "' is marked as '",
                            status, ".\n",
                            "Usage is either not recommended and unsupported.",
                            "\n", sep="")
                var = paste("enable.model.", section, ".", model, sep="")
                if (!(exists(var, envir = envir) && get(var, envir = envir))) {
                    msg = paste0(c(msg, paste("  You can enable this model setting variable '", var,
                        "'\n  to 'TRUE' in environment INLA:::inla.get.inlaEnv().\n",
                        "  If you chose to do so, you are on your own.")))
                    stop(msg)
                }
            } else if (inla.strcasecmp(status.core, "changed")) {
                assign(var, TRUE, envir = envir)
                msg = paste0("Model '", model, "' in section '", section, "' is marked as '",
                            status, ".\n",
                            "  There have been a change in the model definition, which is not backward compatible.\n",
                            "  Please refer to the documentation before proceeeding.")
                var = paste("enable.model.", section, ".", model, sep="")
                if (!(exists(var, envir = envir) && get(var, envir = envir))) {
                    msg = paste0(c(msg, paste0("\n  You can bypass this check setting variable '", var,
                        "'\n  to 'TRUE' in environment INLA:::inla.get.inlaEnv().\n")))
                    stop(msg)
                }
            }
        }
        
        return (models[[k]])
    }
}

`inla.model.validate.link.function` = function(model, link)
{
    valid.links = inla.model.properties(model, "likelihood")$link

    stopifnot(!is.null(valid.links))
    stopifnot(length(valid.links) >= 2L)
    stopifnot(valid.links[1L] == "default")

    link = tolower(link)
    if (is.element(link, valid.links)) {
        ## this is the convention: the default link is the second
        ## entry in the list. the first entry is always "default"
        if (link == "default") {
            link = valid.links[2L]
        }
    } else {
        stop(inla.paste(c("Link function `", link, "' is not valid or yet implemented.",
                          "\n",
                          "Valid ones are: ", inla.paste(valid.links), "\n"), sep=""))
    }

    return (link)
}
