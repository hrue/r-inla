## Export: inla.models

`inla.models.section.latent` = function()
{
    return
    list(latent =
         list(

             linear = list(
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
                 hyper = list(
                     theta = list(
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
                 hyper = list(
                     theta1 = list(
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
                 status = "experimental", 
                 pdf = "mec"
                 ),

             meb = list(
                 hyper = list(
                     theta1 = list(
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
                 status = "experimental", 
                 pdf = "meb"
                 ),

             rgeneric = list(
                 hyper = list(), 
                 constr = FALSE,
                 nrow.ncol = FALSE,
                 augmented = FALSE,
                 aug.factor = 1L,
                 aug.constr = NULL,
                 n.div.by = NULL,
                 n.required = FALSE,
                 set.default.values = FALSE,
                 status = "experimental", 
                 pdf = "rgeneric"
                 ),

             rw1 = list(
                 hyper = list(
                     theta = list(
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
                 pdf = "rw1"
                 ),

             rw2 = list(
                 hyper = list(
                     theta = list(
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
                 pdf = "rw2"
                 ),

             crw2 = list(
                 hyper = list(
                     theta = list(
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
                 pdf = "crw2"
                 ),

             seasonal = list(
                 hyper = list(
                     theta = list(
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
                 hyper = list(
                     theta = list(
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
                 hyper = list(
                     theta1 = list(
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
                 hyper = list(
                     theta1 = list(
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
                 hyper = list(
                     theta1 = list(
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
                         name = "logit phi",
                         short.name = "phi",
                         prior = "pc",
                         param = c(0.5, -1),
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
                 hyper = list(
                     theta1 = list(
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
                 hyper = list(
                     theta1 = list(
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

             ar1 = list(
                 hyper = list(
                     theta1 = list(
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
                         name = "logit lag one correlation",
                         short.name = "rho",
                         prior = "normal",
                         param = c(0, 0.15),
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
                 set.default.values = FALSE,
                 pdf = "ar1"
                 ),

             ar = list(
                 ## to many parameters here, but ...
                 hyper = list(
                     theta1 = list(
                         name = "log precision",
                         short.name = "prec",
                         initial = 4,
                         fixed = FALSE,
                         prior = "pc.prec",
                         param = c(1, 0.01),
                         to.theta = function(x) log(x),
                         from.theta = function(x) exp(x)
                         ),
                     theta2 = list(
                         name = "pacf1",
                         short.name = "pacf1",
                         initial = 1,
                         fixed = FALSE,
                         prior = "pc.rho0",
                         param = c(0.5, 0.5), 
                         to.theta = function(x) log((1+x)/(1-x)),
                         from.theta = function(x) 2*exp(x)/(1+exp(x))-1
                         ),
                     theta3 = list(
                         name = "pacf2",
                         short.name = "pacf2",
                         initial = 0,
                         fixed = FALSE,
                         prior = "pc.rho0",
                         param = c(0.5, 0.4), 
                         to.theta = function(x) log((1+x)/(1-x)),
                         from.theta = function(x) 2*exp(x)/(1+exp(x))-1
                         ),
                     theta4 = list(
                         name = "pacf3",
                         short.name = "pacf3",
                         initial = 0,
                         fixed = FALSE,
                         prior = "pc.rho0",
                         param = c(0.5, 0.3), 
                         to.theta = function(x) log((1+x)/(1-x)),
                         from.theta = function(x) 2*exp(x)/(1+exp(x))-1
                         ),
                     theta5 = list(
                         name = "pacf4",
                         short.name = "pacf4",
                         initial = 0,
                         fixed = FALSE,
                         prior = "pc.rho0",
                         param = c(0.5, 0.2), 
                         to.theta = function(x) log((1+x)/(1-x)),
                         from.theta = function(x) 2*exp(x)/(1+exp(x))-1
                         ),
                     theta6 = list(
                         name = "pacf5",
                         short.name = "pacf5",
                         initial = 0,
                         fixed = FALSE,
                         prior = "pc.rho0",
                         param = c(0.5, 0.1), 
                         to.theta = function(x) log((1+x)/(1-x)),
                         from.theta = function(x) 2*exp(x)/(1+exp(x))-1
                         ),
                     theta7 = list(
                         name = "pacf6",
                         short.name = "pacf6",
                         initial = 0,
                         fixed = FALSE,
                         prior = "pc.rho0",
                         param = c(0.5, 0.1), 
                         to.theta = function(x) log((1+x)/(1-x)),
                         from.theta = function(x) 2*exp(x)/(1+exp(x))-1
                         ),
                     theta8 = list(
                         name = "pacf7",
                         short.name = "pacf7",
                         initial = 0,
                         fixed = FALSE,
                         prior = "pc.rho0",
                         param = c(0.5, 0.1), 
                         to.theta = function(x) log((1+x)/(1-x)),
                         from.theta = function(x) 2*exp(x)/(1+exp(x))-1
                         ),
                     theta9 = list(
                         name = "pacf8",
                         short.name = "pacf8",
                         initial = 0,
                         fixed = FALSE,
                         prior = "pc.rho0",
                         param = c(0.5, 0.1), 
                         to.theta = function(x) log((1+x)/(1-x)),
                         from.theta = function(x) 2*exp(x)/(1+exp(x))-1
                         ),
                     theta10 = list(
                         name = "pacf9",
                         short.name = "pacf9",
                         initial = 0,
                         fixed = FALSE,
                         prior = "pc.rho0",
                         param = c(0.5, 0.1), 
                         to.theta = function(x) log((1+x)/(1-x)),
                         from.theta = function(x) 2*exp(x)/(1+exp(x))-1
                         ),
                     theta11 = list(
                         name = "pacf10",
                         short.name = "pacf10",
                         initial = 0,
                         fixed = FALSE,
                         prior = "pc.rho0",
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
                 status = "experimental", 
                 pdf = "ar"
                 ),

             ou = list(
                 hyper = list(
                     theta1 = list(
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

             generic = list(
                 hyper = list(
                     theta = list(
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
                 hyper = list(
                     theta = list(
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
                 hyper = list(
                     theta1 = list(
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
                 hyper = list(
                     theta1 = list(
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
                 hyper = list(
                     theta1 = list(
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
                         name = "log precision common",
                         short.name = "prec.common",
                         initial = 0,                  ## yes!
                         fixed = TRUE,                 ## yes!
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
                 ## this will be redone anyway soon....
                 hyper = list(
                     theta1 = list(
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
                 ## to many parameters here, but ...
                 hyper = list(
                     theta1 = list(
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
                 ## to many parameters here, but ...
                 hyper = list(
                     theta1 = list(
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
                 hyper = list(
                     theta = list(
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
                 hyper = list(
                     theta1 = list(
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
                 hyper = list(
                     theta1 = list(
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
                 hyper = list(
                     theta1 = list(
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
                 hyper = list(
                     theta1 = list(
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
                 hyper = list(
                     theta1 = list(
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
                 hyper = list(
                     theta = list(
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
                 hyper = list(
                     theta = list(
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
                 hyper = list(
                     theta1 = list(
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
                         name = "logit phi",
                         short.name = "phi",
                         prior = "pc",
                         param = c(0.5, -1),
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
                 hyper = list(
                     theta1 = list(
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
                 hyper = list(
                     theta1 = list(
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

             copy = list(
                 hyper = list(
                     theta = list(
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
                 hyper = list(
                     theta = list(
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
                 hyper = list(
                     theta1 = list(
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
                 hyper = list(
                     theta1 = list(
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
                 hyper = list(
                     theta = list(
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

             ar1 = list(
                 hyper = list(
                     theta = list(
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
                 ## to many parameters here, but ...
                 hyper = list(
                     theta1 = list(
                         name = "log precision",
                         short.name = "prec",
                         initial = 0,
                         fixed = TRUE,
                         prior = "loggamma",
                         param = c(1, 0.00005),
                         to.theta = function(x) log(x),
                         from.theta = function(x) exp(x)
                         ),
                     theta2 = list(
                         name = "pacf1",
                         short.name = "pacf1",
                         initial = 2,
                         fixed = FALSE,
                         prior = "mvnorm",
                         param = c(0, 0.15), ## same as for AR1
                         to.theta = function(x) log((1+x)/(1-x)),
                         from.theta = function(x) 2*exp(x)/(1+exp(x))-1
                         ),
                     theta3 = list(
                         name = "pacf2",
                         short.name = "pacf2",
                         initial = 0,
                         fixed = FALSE,
                         prior = "none",
                         param = numeric(0),
                         to.theta = function(x) log((1+x)/(1-x)),
                         from.theta = function(x) 2*exp(x)/(1+exp(x))-1
                         ),
                     theta4 = list(
                         name = "pacf3",
                         short.name = "pacf3",
                         initial = 0,
                         fixed = FALSE,
                         prior = "none",
                         param = numeric(0),
                         to.theta = function(x) log((1+x)/(1-x)),
                         from.theta = function(x) 2*exp(x)/(1+exp(x))-1
                         ),
                     theta5 = list(
                         name = "pacf4",
                         short.name = "pacf4",
                         initial = 0,
                         fixed = FALSE,
                         prior = "none",
                         param = numeric(0),
                         to.theta = function(x) log((1+x)/(1-x)),
                         from.theta = function(x) 2*exp(x)/(1+exp(x))-1
                         ),
                     theta6 = list(
                         name = "pacf5",
                         short.name = "pacf5",
                         initial = 0,
                         fixed = FALSE,
                         prior = "none",
                         param = numeric(0),
                         to.theta = function(x) log((1+x)/(1-x)),
                         from.theta = function(x) 2*exp(x)/(1+exp(x))-1
                         ),
                     theta7 = list(
                         name = "pacf6",
                         short.name = "pacf6",
                         initial = 0,
                         fixed = FALSE,
                         prior = "none",
                         param = numeric(0),
                         to.theta = function(x) log((1+x)/(1-x)),
                         from.theta = function(x) 2*exp(x)/(1+exp(x))-1
                         ),
                     theta8 = list(
                         name = "pacf7",
                         short.name = "pacf7",
                         initial = 0,
                         fixed = FALSE,
                         prior = "none",
                         param = numeric(0),
                         to.theta = function(x) log((1+x)/(1-x)),
                         from.theta = function(x) 2*exp(x)/(1+exp(x))-1
                         ),
                     theta9 = list(
                         name = "pacf8",
                         short.name = "pacf8",
                         initial = 0,
                         fixed = FALSE,
                         prior = "none",
                         param = numeric(0),
                         to.theta = function(x) log((1+x)/(1-x)),
                         from.theta = function(x) 2*exp(x)/(1+exp(x))-1
                         ),
                     theta10 = list(
                         name = "pacf9",
                         short.name = "pacf9",
                         initial = 0,
                         fixed = FALSE,
                         prior = "none",
                         param = numeric(0),
                         to.theta = function(x) log((1+x)/(1-x)),
                         from.theta = function(x) 2*exp(x)/(1+exp(x))-1
                         ),
                     theta11 = list(
                         name = "pacf10",
                         short.name = "pacf10",
                         initial = 0,
                         fixed = FALSE,
                         prior = "none",
                         param = numeric(0),
                         to.theta = function(x) log((1+x)/(1-x)),
                         from.theta = function(x) 2*exp(x)/(1+exp(x))-1
                         )
                     )
                 ),

             rw1 = list(
                 hyper = list(
                     theta = list(
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
                 hyper = list(
                     theta = list(
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
                 hyper = list(
                     theta = list(
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
             I = list(
                 hyper = list(
                     theta = list(
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
                 hyper = list(
                     theta = list(
                         name = "log precision",
                         short.name = "prec",
                         prior = "loggamma",
                         param = c(1, 0.01),
                         initial = 0,
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
             default = list(hyper = list()),
             cloglog = list(hyper = list()), 
             identity = list(hyper = list()), 
             log = list(hyper = list()), 
             logit = list(hyper = list()), 
             probit = list(hyper = list()), 
             tan = list(hyper = list()), 
             sslogit = list(
                 hyper = list(
                     theta1 = list(
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
                 ## variant = 0, a+exp(...),  a>0
                 ## variant = 1, a-exp(...),  a>0
                 hyper = list(
                     theta = list(
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

             test1 = list(
                 hyper = list(
                     theta = list(
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
                 hyper = list(
                     theta1 = list(
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
                 hyper = list(
                     theta = list(
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
                 hyper = list(
                     theta = list(
                         name = "log precision",
                         short.name = "prec",
                         initial = 11,
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
                 hyper = list(
                     theta = list(
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
                 hyper = list(
                     theta = list(
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
                 hyper = list(
                     ),
                 survival = FALSE,
                 discrete = TRUE,
                 link = c("default", "log", "logoffset", "test1", "special1", "special2"),
                 pdf = "poisson"
                 ),

             gpoisson = list(
                 hyper = list(
                     theta1 = list(
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
                 hyper = list(
                     ),
                 survival = FALSE,
                 discrete = TRUE,
                 link = c("default", "logit", "probit", "cloglog", "log", "sslogit"),
                 pdf = "binomial"
                 ),

             testbinomial1 = list(
                 hyper = list(
                     theta1 = list(
                         name = "sensitivity",
                         short.name = "s",
                         initial = 3,
                         fixed = FALSE,
                         prior = "logitbeta",
                         param = c(2, 1),
                         to.theta = function(x) log(x/(1-x)),
                         from.theta = function(x) exp(x)/(1+exp(x))
                         ), 
                     theta2 = list(
                         name = "specificity",
                         short.name = "e",
                         initial = 3,
                         fixed = FALSE,
                         prior = "logitbeta",
                         param = c(2, 1),
                         to.theta = function(x) log(x/(1-x)),
                         from.theta = function(x) exp(x)/(1+exp(x))
                         )
                     ),
                 status = "experimental", 
                 survival = FALSE,
                 discrete = TRUE,
                 link = c("default", "logit", "probit", "cloglog", "log"),
                 pdf = "testbinomial1"
                 ),

             gamma = list(
                 hyper = list(
                     theta = list(
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
                 link = c("default", "log"),
                 pdf = "gamma"
                 ),

             gammacount = list(
                 hyper = list(
                     theta = list(
                         name = "log alpha",
                         short.name = "alpha",
                         initial = log(1),
                         fixed = FALSE,
                         prior = "loggamma",
                         param = c(10, 10),
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

             kumar = list(
                 hyper = list(
                     theta1 = list(
                         name = "precision parameter",
                         short.name = "prec",
                         initial = 0,
                         fixed = FALSE,
                         prior = "loggamma",
                         param = c(1, 0.001),
                         to.theta = function(x) log(x), 
                         from.theta = function(x) exp(x)
                         ), 
                     theta2 = list(
                         name = "quantile",
                         short.name = "q",
                         initial = 0.5,
                         fixed = TRUE,
                         prior = "invalid",
                         param = numeric(0),
                         to.theta = function(x) x, 
                         from.theta = function(x) x
                         )
                     ),
                 survival = FALSE,
                 discrete = FALSE,
                 link = c("default", "logit"),
                 pdf = "kumar"
                 ),

             beta = list(
                 hyper = list(
                     theta = list(
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
                 link = c("default", "logit", "probit", "cloglog"),
                 pdf = "beta"
                 ),

             betabinomial = list(
                 hyper = list(
                     theta = list(
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
                 link = c("default", "logit", "probit", "cloglog"),
                 pdf = "betabinomial"
                 ),

             cbinomial = list(
                 hyper = list(
                     ),
                 survival = FALSE,
                 discrete = TRUE,
                 link = c("default", "logit", "probit", "cloglog"),
                 status = "experimental", 
                 pdf = "cbinomial"
                 ),

             nbinomial = list(
                 hyper = list(
                     theta = list(
                         name = "size",
                         short.name = "size",
                         initial = log(10),
                         fixed = FALSE,
                         prior = "loggamma",
                         param = c(1, 1),
                         to.theta = function(x) log(x),
                         from.theta = function(x) exp(x)
                         )
                     ),
                 survival = FALSE,
                 discrete = TRUE,
                 link = c("default", "log", "logoffset"),
                 pdf = "nbinomial"
                 ),

             simplex = list(
                 hyper = list(
                     theta = list(
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
                 link = c("default", "logit", "probit", "cloglog"),
                 pdf = "simplex"
                 ),

             gaussian = list(
                 hyper = list(
                     theta = list(
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
                 link = c("default", "identity", "logit", "log", "logoffset"),
                 pdf = "gaussian"
                 ),

             normal = list(
                 hyper = list(
                     theta = list(
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
                 pdf = "gaussian"
                 ),

             circularnormal = list(
                 hyper = list(
                     theta = list(
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
                 hyper = list(
                     theta = list(
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
                 hyper = list(
                     theta1 = list(
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
                 hyper = list(
                     theta1 = list(
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
                 hyper = list(
                     theta = list(
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
                 pdf = "loggammafrailty"
                 ),

             logistic = list(
                 hyper = list(
                     theta = list(
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
                 hyper = list(
                     theta1 = list(
                         name = "inverse.scale",
                         short.name = "iscale",
                         initial = 4,
                         fixed = FALSE,
                         prior = "loggamma",
                         param = c(1, 0.00005)
                         ),
                     theta2 = list(
                         name = "skewness",
                         short.name = "skew",
                         initial = 4,
                         fixed = FALSE,
                         prior = "gaussian",
                         param = c(0, 10)
                         )
                     ),
                 survival = FALSE,
                 discrete = FALSE,
                 link = c("default", "identity"),
                 pdf = "sn"
                 ),

             sn = list(
                 hyper = list(
                     theta1 = list(
                         name = "log inverse scale",
                         short.name = "iscale",
                         initial = 4,
                         fixed = FALSE,
                         prior = "loggamma",
                         param = c(1, 0.00005)
                         ),
                     theta2 = list(
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
                 hyper = list(
                     theta1 = list(
                         name = "log precision",
                         short.name = "prec",
                         initial = 1,
                         fixed = FALSE,
                         prior = "loggamma",
                         param = c(1, 0.00005)
                         ),
                     theta2 = list(
                         name = "logit skewness",
                         short.name = "skew",
                         initial = 0,
                         fixed = FALSE,
                         prior = "gaussian",
                         param = c(0, 10),
                         to.theta = function(x) log((1+x)/(1-x)),
                         from.theta = function(x) (2*exp(x)/(1+exp(x))-1)
                         )
                     ),
                 survival = FALSE,
                 discrete = FALSE,
                 link = c("default", "identity"),
                 status = "experimental", 
                 pdf = "sn2"
                 ),

             gev = list(
                 hyper = list(
                     theta1 = list(
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
                         name = "gev parameter",
                         short.name = "gev",
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

             laplace = list(
                 hyper = list(
                     theta = list(
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
                 status = "disabled", 
                 pdf = "laplace"
                 ),

             lognormal = list(
                 hyper = list(
                     theta = list(
                         name = "log precision",
                         short.name = "prec",
                         initial = 2,
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
                 hyper = list(
                     ),
                 survival = TRUE,
                 discrete = FALSE,
                 link = c("default", "log"),
                 pdf = "exponential"
                 ),

             coxph = list(
                 hyper = list(
                     ),
                 survival = TRUE,
                 discrete = TRUE,
                 link = c("default", "log"),
                 pdf = "coxph"
                 ),

             weibull = list(
                 hyper = list(
                     theta = list(
                         name = "log alpha",
                         short.name = "a",
                         initial = 0,
                         fixed = FALSE,
                         prior = "loggamma",
                         param = c(25, 25),
                         to.theta = function(x) log(x),
                         from.theta = function(x) exp(x)
                         )
                     ),
                 survival = TRUE,
                 discrete = FALSE,
                 link = c("default", "log"),
                 pdf = "weibull"
                 ),

             loglogistic = list(
                 hyper = list(
                     theta = list(
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
                 link = c("default", "log"),
                 pdf = "loglogistic"
                 ),

             weibullcure = list(
                 hyper = list(
                     theta1 = list(
                         name = "log alpha",
                         short.name = "a",
                         initial = 4,
                         fixed = FALSE,
                         prior = "loggamma",
                         param = c(25, 25),
                         to.theta = function(x) log(x),
                         from.theta = function(x) exp(x)
                         ),
                     theta2 = list(
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
                 link = c("default", "log"),
                 pdf = NA
                 ),

             stochvol = list(
                 hyper = list(
                     theta = list(
                         name = "log precision",
                         short.name = "prec",
                         initial = 500,             ## yes, this is correct
                         fixed = TRUE,              ## yes, this is correct
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
                 hyper = list(
                     theta = list(
                         name = "log degrees of freedom",
                         short.name = "dof",
                         initial = 4,
                         fixed = FALSE,
                         prior = "loggamma",
                         param = c(1, 0.5),
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
                 hyper = list(
                     theta1 = list(
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
                 hyper = list(
                     theta = list(
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
                 hyper = list(
                     theta = list(
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
                 hyper = list(
                     theta = list(
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
                 hyper = list(
                     theta1 = list(
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
                 link = c("default", "logit", "probit", "cloglog"),
                 pdf = "zeroinflated"
                 ),

             zeroinflatedbetabinomial1 = list(
                 hyper = list(
                     theta1 = list(
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
                 link = c("default", "logit", "probit", "cloglog"),
                 pdf = "zeroinflated"
                 ),

             zeroinflatedbinomial0 = list(
                 hyper = list(
                     theta = list(
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
                 link = c("default", "logit", "probit", "cloglog"),
                 pdf = "zeroinflated"
                 ),

             zeroinflatedbinomial1 = list(
                 hyper = list(
                     theta = list(
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
                 link = c("default", "logit", "probit", "cloglog"),
                 pdf = "zeroinflated"
                 ),

             zeroinflatedbinomial2 = list(
                 hyper = list(
                     theta = list(
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
                 link = c("default", "logit", "probit", "cloglog"),
                 pdf = "zeroinflated"
                 ),

             zeroninflatedbinomial2 = list(
                 hyper = list(
                     theta1 = list(
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
                 link = c("default", "logit", "probit", "cloglog"),
                 pdf = NA
                 ),

             zeroinflatedbetabinomial2 = list(
                 hyper = list(
                     theta1 = list(
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
                 link = c("default", "logit", "probit", "cloglog"),
                 pdf = "zeroinflated"
                 ),

             zeroinflatednbinomial0 = list(
                 hyper = list(
                     theta1 = list(
                         name = "log size",
                         short.name = "size",
                         initial = log(10),
                         fixed = FALSE,
                         prior = "loggamma",
                         param = c(1, 1),
                         to.theta = function(x) log(x),
                         from.theta = function(x) exp(x)
                         ),
                     theta2 = list(
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
                 hyper = list(
                     theta1 = list(
                         name = "log size",
                         short.name = "size",
                         initial = log(10),
                         fixed = FALSE,
                         prior = "loggamma",
                         param = c(1, 1),
                         to.theta = function(x) log(x),
                         from.theta = function(x) exp(x)
                         ),
                     theta2 = list(
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
                 hyper = list(
                     theta1 = list(
                         name = "log size",
                         short.name = "size",
                         initial = log(10),
                         fixed = FALSE,
                         prior = "loggamma",
                         param = c(1, 1),
                         to.theta = function(x) log(x),
                         from.theta = function(x) exp(x)
                         ),
                     theta2 = list(
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
                         name = "logit probability 2",
                         short.name = "prob2",
                         initial = -1,
                         fixed = FALSE,
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
                 hyper = list(
                     theta1 = list(
                         name = "log size 1",
                         short.name = "size1",
                         initial = log(10),
                         fixed = FALSE,
                         prior = "loggamma",
                         param = c(1, 1),
                         to.theta = function(x) log(x),
                         from.theta = function(x) exp(x)
                         ),
                     theta2 = list(
                         name = "log size 2",
                         short.name = "size2",
                         initial = log(10),
                         fixed = FALSE,
                         prior = "loggamma",
                         param = c(1, 1),
                         to.theta = function(x) log(x),
                         from.theta = function(x) exp(x)
                         ),
                     theta3 = list(
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
                 status = "experimental", 
                 survival = FALSE,
                 discrete = FALSE,
                 link = c("default", "log"),
                 pdf = "zeroinflated"
                 ),

             zeroinflatednbinomial2 = list(
                 hyper = list(
                     theta1 = list(
                         name = "log size",
                         short.name = "size",
                         initial = log(10),
                         fixed = FALSE,
                         prior = "loggamma",
                         param = c(1, 1),
                         to.theta = function(x) log(x),
                         from.theta = function(x) exp(x)
                         ),
                     theta2 = list(
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
                 hyper = list(
                     theta1 = list(
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
                         name = "log degrees of freedom",
                         short.name = "dof",
                         initial = 5,
                         fixed = FALSE,
                         prior = "loggamma",
                         param = c(1, 0.5),
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
                 hyper = list(
                     theta1 = list(
                         name = "log degrees of freedom",
                         short.name = "dof",
                         initial = 4,
                         fixed = FALSE,
                         prior = "loggamma",
                         param = c(1, 0.01),
                         to.theta = function(x) log(x-5),
                         from.theta = function(x) 5+exp(x)
                         ),
                     theta2 = list(
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

             logperiodogram = list(
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
                 nparameters = 2L,
                 pdf = "gaussian"
                 ),
             gaussian = list(
                 nparameters = 2L,
                 pdf = "gaussian"
                 ),
             wishart1d = list(
                 nparameters = 2L,
                 pdf = "iid123d"
                 ),
             wishart2d = list(
                 nparameters = 4L,
                 pdf = "iid123d"
                 ),
             wishart3d = list(
                 nparameters = 7L,
                 pdf = "iid123d"
                 ),
             wishart4d = list(
                 nparameters = 11L,
                 pdf = "iid123d"
                 ),
             wishart5d = list(
                 nparameters = 16L,
                 pdf = "iid123d"
                 ),
             loggamma = list(
                 nparameters = 2L,
                 pdf = "prior-loggamma"
                 ),
             minuslogsqrtruncnormal = list(
                 nparameters = 2L,
                 pdf = "prior-logtnorm"
                 ),
             logtnormal = list(
                 nparameters = 2L,
                 pdf = "prior-logtnorm"
                 ),
             logtgaussian = list(
                 nparameters = 2L,
                 pdf = "prior-logtnorm"
                 ),
             flat=list(
                 nparameters = 0L,
                 pdf = "various-flat"
                 ),
             logflat=list(
                 nparameters = 0L,
                 pdf = "various-flat"
                 ),
             logiflat=list(
                 nparameters = 0L,
                 pdf = "various-flat"
                 ),

             mvnorm = list(
                 nparameters = -1L,
                 pdf = "mvnorm"
                 ),

             pc.ar = list(
                 nparameters = 1L,
                 pdf = "pc.ar"
                 ),

             ## this is the 'no prior needed' prior
             none = list(nparameters = 0L),

             ## this is the 'flag an error if used' prior
             invalid = list(nparameters = 0L),

             betacorrelation = list(
                 nparameters = 2L,
                 pdf = "betacorrelation"
                 ),

             logitbeta = list(
                 nparameters = 2L,
                 pdf = "logitbeta"
                 ),

             pc.prec = list(
                 nparameters = 2L,
                 pdf = "pc.prec"
                 ),
             
             pc.dof = list(
                 nparameters = 2L,
                 pdf = "pc.dof"
                 ),
             
             pc.rho0 = list(
                 nparameters = 2L,
                 pdf = "pc.rho0"
                 ),
             
             pc.rho1 = list(
                 nparameters = 2L,
                 pdf = "pc.rho1"
                 ),
             
             ## experimental prior from GA
             pc.spde.GA = list(
                 nparameters = 4L,
                 pdf = NA
                 ), 

             ## this is the generic one,  which is case-spesific and possibly adaptive
             pc = list(
                 nparameters = 2L,
                 pdf = NA
                 ), 

             ref.ar = list(
                 nparameters = 0L,
                 pdf = NA
                 ),
             
             jeffreystdf = list(
                 nparameters = 0L,
                 pdf = "jeffreystdf"
                 ),

             "expression:" = list(
                 nparameters = -1L,
                 pdf = "expression"
                 ), 

             "table:" = list(
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
                 hyper = list(
                     theta = list(
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

`inla.is.model` = function(model, section = names(inla.models()),
        stop.on.error = TRUE, ignore.case = FALSE)
{
    section = match.arg(section)
    mm = inla.models()
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
        section = c("..invalid.model..", names(inla.models())),
        stop.on.error = TRUE,
        ignore.case = FALSE)
{
    section = match.arg(section)
    if (section == "..invalid.model..")
        stop("No section given; please fix...")

    mm = inla.models()
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
            stop(paste("Name '", model,  "' is unknown. Valid choices are: ", inla.paste(ms), ".", sep=""))
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
            stopifnot(any(inla.strcasecmp(status, c("experimental", "disabled"))))
            envir = inla.get.inlaEnv()
            var = paste("processed.status.for.model.", model, ".in.section.", section, sep="")

            if (inla.strcasecmp(status, "experimental")) {
                if (!(exists(var, envir = envir) && get(var, envir = envir))) {
                    assign(var, TRUE, envir = envir)
                    msg = paste("Model '", model, "' in section '", section, "' is marked as '", status, 
                        "'; changes may appear at any time.",
                        "\n  ",
                        "Use this model with extra care!!! Further warnings are disabled.", sep="")
                    warning(msg)
                } else {
                    ## the warning is already given; do nothing
                }
            } else if (inla.strcasecmp(status, "disabled")) {
                assign(var, TRUE, envir = envir)
                msg = paste("Model '", model, "' in section '", section, "' is marked as '", status,
                    "'.\n  Usage is not recommended and unsupported.\n", sep="")
                var = paste("enable.model.", section, ".", model, sep="")
                if (!(exists(var, envir = envir) && get(var, envir = envir))) {
                    msg = paste(c(msg, paste("  You can enable it setting variable '", var,
                        "'\n  to 'TRUE' in environment INLA:::inla.get.inlaEnv().\n",
                        "  If you chose to do so, you are on your own.", sep="")))
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
