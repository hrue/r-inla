## models and its hyperparameters are defined here

`inla.models` = function()
{
    return (
            list(
                 ## latent models
                 latent = list(

                         linear = list(
                                 hyper = list(
                                         ),
                                 constr = FALSE,
                                 nrow.ncol = FALSE,
                                 augmented = FALSE,
                                 aug.factor = 1,
                                 aug.constr = NULL,
                                 n.div.by = NULL,
                                 n.required = FALSE,
                                 set.default.values = FALSE
                                 ), 
                                 
                         iid = list(
                                 hyper = list(
                                         theta = list(
                                                 name = "precision",
                                                 short.name = "prec",
                                                 initial = 4,
                                                 fixed = FALSE,
                                                 prior = "loggamma",
                                                 param = c(1, 0.0001)
                                                 )
                                         ), 
                                 constr = FALSE,
                                 nrow.ncol = FALSE,
                                 augmented = FALSE,
                                 aug.factor = 1,
                                 aug.constr = NULL,
                                 n.div.by = NULL,
                                 n.required = FALSE,
                                 set.default.values = FALSE
                                 ), 

                         rw1 = list(
                                 hyper = list(
                                         theta = list(
                                                 name = "precision",
                                                 short.name = "prec",
                                                 initial = 4,
                                                 fixed = FALSE,
                                                 prior = "loggamma",
                                                 param = c(1, 0.0001)
                                                 )
                                         ), 
                                 constr = TRUE,
                                 nrow.ncol = FALSE,
                                 augmented = FALSE,
                                 aug.factor = 1,
                                 aug.constr = NULL,
                                 n.div.by = NULL,
                                 n.required = FALSE,
                                 set.default.values = FALSE
                                 ), 

                         rw2 = list(
                                 hyper = list(
                                         theta = list(
                                                 name = "precision",
                                                 short.name = "prec",
                                                 initial = 4,
                                                 fixed = FALSE,
                                                 prior = "loggamma",
                                                 param = c(1, 0.0001)
                                                 )
                                         ), 
                                 constr = TRUE,
                                 nrow.ncol = FALSE,
                                 augmented = FALSE,
                                 aug.factor = 1,
                                 aug.constr = NULL,
                                 n.div.by = NULL,
                                 n.required = FALSE,
                                 set.default.values = FALSE
                                 ), 

                         crw2 = list(
                                 hyper = list(
                                         theta = list(
                                                 name = "precision",
                                                 short.name = "prec",
                                                 initial = 4,
                                                 fixed = FALSE,
                                                 prior = "loggamma",
                                                 param = c(1, 0.0001)
                                                 )
                                         ), 
                                 constr = TRUE,
                                 nrow.ncol = FALSE,
                                 augmented = FALSE,
                                 aug.factor = 2,
                                 aug.constr = 1,
                                 n.div.by = NULL,
                                 n.required = FALSE,
                                 set.default.values = FALSE
                                 ), 

                         seasonal = list(
                                 hyper = list(
                                         theta = list(
                                                 name = "precision",
                                                 short.name = "prec",
                                                 initial = 4,
                                                 fixed = FALSE,
                                                 prior = "loggamma",
                                                 param = c(1, 0.0001)
                                                 )
                                         ), 
                                 constr = FALSE,
                                 nrow.ncol = FALSE,
                                 augmented = FALSE,
                                 aug.factor = 1,
                                 aug.constr = NULL,
                                 n.div.by = NULL,
                                 n.required = FALSE,
                                 set.default.values = FALSE
                                 ), 

                         besag = list(
                                 hyper = list(
                                         theta = list(
                                                 name = "precision",
                                                 short.name = "prec",
                                                 initial = 4,
                                                 fixed = FALSE,
                                                 prior = "loggamma",
                                                 param = c(1, 0.0001)
                                                 )
                                         ), 
                                 constr = TRUE,
                                 nrow.ncol = FALSE,
                                 augmented = FALSE,
                                 aug.factor = 1,
                                 aug.constr = NULL,
                                 n.div.by = NULL,
                                 n.required = TRUE,
                                 set.default.values = TRUE
                                 ), 

                         besag2 = list(
                                 hyper = list(
                                         theta1 = list(
                                                 name = "precision",
                                                 short.name = "prec",
                                                 initial = 4,
                                                 fixed = FALSE,
                                                 prior = "loggamma",
                                                 param = c(1, 0.0001)
                                                 ), 
                                         theta2 = list(
                                                 name = "a",
                                                 short.name = "a",
                                                 initial = 4,
                                                 fixed = FALSE,
                                                 prior = "normal",
                                                 param = c(0, 6.25)
                                                 )
                                         ), 
                                 constr = FALSE,
                                 nrow.ncol = FALSE,
                                 augmented = FALSE,
                                 aug.factor = 1,
                                 aug.constr = c(1, 2),
                                 n.div.by = 2,
                                 n.required = TRUE,
                                 set.default.values = TRUE
                                 ), 

                         bym = list(
                                 hyper = list(
                                         theta1 = list(
                                                 name = "precision iid",
                                                 short.name = "prec.iid",
                                                 initial = 4,
                                                 fixed = FALSE,
                                                 prior = "loggamma",
                                                 param = c(1, 0.0001)
                                                 ), 
                                         theta2 = list(
                                                 name = "precision spatial",
                                                 short.name = "prec.spatial",
                                                 initial = 4,
                                                 fixed = FALSE,
                                                 prior = "normal",
                                                 param = c(0, 0.0001)
                                                 )
                                         ), 
                                 constr = TRUE,
                                 nrow.ncol = FALSE,
                                 augmented = TRUE,
                                 aug.factor = 2,
                                 aug.constr = 2,
                                 n.div.by = NULL,
                                 n.required = TRUE,
                                 set.default.values = TRUE
                                 ), 

                         ar1 = list(
                                 hyper = list(
                                         theta1 = list(
                                                 name = "precision",
                                                 short.name = "prec",
                                                 initial = 4,
                                                 fixed = FALSE,
                                                 prior = "loggamma",
                                                 param = c(1, 0.0001)
                                                 ), 
                                         theta2 = list(
                                                 name = "lag-one correlation",
                                                 short.name = "rho",
                                                 initial = 2,
                                                 fixed = FALSE,
                                                 prior = "normal",
                                                 param = c(0, 0.2)
                                                 )
                                         ), 
                                 constr = FALSE,
                                 nrow.ncol = FALSE,
                                 augmented = FALSE,
                                 aug.factor = 1,
                                 aug.constr = NULL,
                                 n.div.by = NULL,
                                 n.required = FALSE,
                                 set.default.values = FALSE
                                 ), 
                         
                         generic = list(
                                 hyper = list(
                                         theta = list(
                                                 name = "precision",
                                                 short.name = "prec",
                                                 initial = 4,
                                                 fixed = FALSE,
                                                 prior = "loggamma",
                                                 param = c(1, 0.0001)
                                                 )
                                         ), 
                                 constr = FALSE,
                                 nrow.ncol = FALSE,
                                 augmented = FALSE,
                                 aug.factor = 1,
                                 aug.constr = NULL,
                                 n.div.by = NULL,
                                 n.required = TRUE,
                                 set.default.values = TRUE
                                 ), 
                         
                         generic0 = list(
                                 hyper = list(
                                         theta = list(
                                                 name = "precision",
                                                 short.name = "prec",
                                                 initial = 4,
                                                 fixed = FALSE,
                                                 prior = "loggamma",
                                                 param = c(1, 0.0001)
                                                 )
                                         ), 
                                 constr = FALSE,
                                 nrow.ncol = FALSE,
                                 augmented = FALSE,
                                 aug.factor = 1,
                                 aug.constr = NULL,
                                 n.div.by = NULL,
                                 n.required = TRUE,
                                 set.default.values = TRUE
                                 ), 
                         
                         generic1 = list(
                                 hyper = list(
                                         theta1 = list(
                                                 name = "precision",
                                                 short.name = "prec",
                                                 initial = 4,
                                                 fixed = FALSE,
                                                 prior = "loggamma",
                                                 param = c(1, 0.0001)
                                                 ), 
                                         theta2 = list(
                                                 name = "beta",
                                                 short.name = "beta",
                                                 initial = 2,
                                                 fixed = FALSE,
                                                 prior = "gaussian",
                                                 param = c(0, 0.1)
                                                 )
                                         ), 
                                 constr = FALSE,
                                 nrow.ncol = FALSE,
                                 augmented = FALSE,
                                 aug.factor = 1,
                                 aug.constr = NULL,
                                 n.div.by = NULL,
                                 n.required = TRUE,
                                 set.default.values = TRUE
                                 ), 

                         generic2 = list(
                                 hyper = list(
                                         theta1 = list(
                                                 name = "precision-cmatrix",
                                                 short.name = "prec",
                                                 initial = 4,
                                                 fixed = FALSE,
                                                 prior = "loggamma",
                                                 param = c(1, 0.0001)
                                                 ), 
                                         theta2 = list(
                                                 name = "h2",
                                                 short.name = "h2",
                                                 initial = 4,
                                                 fixed = FALSE,
                                                 prior = "loggamma",
                                                 param = c(0, 0.0001)
                                                 )
                                         ), 
                                 constr = FALSE,
                                 nrow.ncol = FALSE,
                                 augmented = FALSE,
                                 aug.factor = 2,
                                 aug.constr = 2,
                                 n.div.by = NULL,
                                 n.required = TRUE,
                                 set.default.values = TRUE
                                 ), 

                         spde = list(
                                 hyper = list(
                                         theta1 = list(
                                                 name = "theta.T",
                                                 short.name = "T",
                                                 initial = -3,
                                                 fixed = FALSE,
                                                 prior = "normal",
                                                 param = c(0, 1)
                                                 ),
                                         theta2 = list(
                                                 name = "theta.K",
                                                 short.name = "K",
                                                 initial = 3,
                                                 fixed = FALSE,
                                                 prior = "normal",
                                                 param = c(0, 1)
                                                 ), 
                                         theta3 = list(
                                                 name = "theta.KT",
                                                 short.name = "KT",
                                                 initial = 0,
                                                 fixed = FALSE,
                                                 prior = "normal",
                                                 param = c(0, 1)
                                                 ), 
                                         theta4 = list(
                                                 name = "theta.OC",
                                                 short.name = "OC",
                                                 initial = -20,
                                                 fixed = TRUE,
                                                 prior = "flat",
                                                 param = numeric()
                                                 )
                                         ),
                                 constr = FALSE,
                                 nrow.ncol = FALSE,
                                 augmented = FALSE,
                                 aug.factor = 1,
                                 aug.constr = NULL,
                                 n.div.by = NULL,
                                 n.required = TRUE,
                                 set.default.values = TRUE
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
                                                 ## parameters
                                                 param = c(2*1, 2*0.0001)
                                                 )
                                         ), 
                                 constr = FALSE,
                                 nrow.ncol = FALSE,
                                 augmented = FALSE,
                                 aug.factor = 1,
                                 aug.constr = NULL,
                                 n.div.by = NULL,
                                 n.required = TRUE,
                                 set.default.values = TRUE
                                 ), 

                         iid2d = list(
                                 hyper = list(
                                         theta1 = list(
                                                 name = "precision1",
                                                 short.name = "prec1",
                                                 initial = 4,
                                                 fixed = FALSE,
                                                 prior = "wishart2d",
                                                 param = c(4, 1, 1, 0)
                                                 ), 
                                         theta2 = list(
                                                 name = "precision2",
                                                 short.name = "prec2",
                                                 initial = 4,
                                                 fixed = FALSE,
                                                 prior = "none",
                                                 param = numeric()
                                                 ), 
                                         theta3 = list(
                                                 name = "correlation",
                                                 short.name = "cor",
                                                 initial = 4,
                                                 fixed = FALSE,
                                                 prior = "none",
                                                 param = numeric()
                                                 )
                                         ), 
                                 constr = FALSE,
                                 nrow.ncol = FALSE,
                                 augmented = TRUE,
                                 aug.factor = 1,
                                 aug.constr = c(1, 2),
                                 n.div.by = 2,
                                 n.required = TRUE,
                                 set.default.values = TRUE
                                 ), 

                         iid3d = list(
                                 hyper = list(
                                         theta1 = list(
                                                 name = "precision1",
                                                 short.name = "prec1",
                                                 initial = 4,
                                                 fixed = FALSE,
                                                 prior = "wishart3d",
                                                 param = c(7, 1, 1, 1, 0, 0, 0)
                                                 ), 
                                         theta2 = list(
                                                 name = "precision2",
                                                 short.name = "prec2",
                                                 initial = 4,
                                                 fixed = FALSE,
                                                 prior = "none",
                                                 param = numeric()
                                                 ), 
                                         theta3 = list(
                                                 name = "precision3",
                                                 short.name = "prec3",
                                                 initial = 4,
                                                 fixed = FALSE,
                                                 prior = "none",
                                                 param = numeric()
                                                 ), 
                                         theta4= list(
                                                 name = "correlation12",
                                                 short.name = "cor12",
                                                 initial = 0,
                                                 fixed = FALSE,
                                                 prior = "none",
                                                 param = numeric()
                                                 ), 
                                         theta5 = list(
                                                 name = "correlation13",
                                                 short.name = "cor13",
                                                 initial = 0,
                                                 fixed = FALSE,
                                                 prior = "none",
                                                 param = numeric()
                                                 ), 
                                         theta6 = list(
                                                 name = "correlation23",
                                                 short.name = "cor23",
                                                 initial = 0,
                                                 fixed = FALSE,
                                                 prior = "none",
                                                 param = numeric()
                                                 )
                                         ), 
                                 constr = FALSE,
                                 nrow.ncol = FALSE,
                                 augmented = TRUE,
                                 aug.factor = 1,
                                 aug.constr = c(1, 2, 3),
                                 n.div.by = 3,
                                 n.required = TRUE,
                                 set.default.values = TRUE
                                 ), 

                         "2diid" = list(
                                 hyper = list(
                                         theta1 = list(
                                                 name = "precision1",
                                                 short.name = "prec1",
                                                 initial = 4,
                                                 fixed = FALSE,
                                                 prior = "loggamma",
                                                 param = c(1, 0.0001)
                                                 ), 
                                         theta2 = list(
                                                 name = "precision2",
                                                 short.name = "prec2",
                                                 initial = 4,
                                                 fixed = FALSE,
                                                 prior = "loggamma",
                                                 param = c(1, 0.0001)
                                                 ), 
                                         theta3 = list(
                                                 name = "correlation",
                                                 short.name = "cor",
                                                 initial = 4,
                                                 fixed = FALSE,
                                                 prior = "normal",
                                                 param = c(0, 0.2)
                                                 )
                                         ), 
                                 constr = FALSE,
                                 nrow.ncol = FALSE,
                                 augmented = FALSE,
                                 aug.factor = 1,
                                 aug.constr = c(1, 2),
                                 n.div.by = 2,
                                 n.required = TRUE,
                                 set.default.values = TRUE
                                 ), 
                         
                         z = list(
                                 hyper = list(
                                         theta = list(
                                                 name = "precision",
                                                 short.name = "prec",
                                                 initial = 4,
                                                 fixed = FALSE,
                                                 prior = "loggamma",
                                                 param = c(1, 0.0001)
                                                 )
                                         ), 
                                 constr = FALSE,
                                 nrow.ncol = FALSE,
                                 augmented = FALSE,
                                 aug.factor = 1,
                                 aug.constr = NULL,
                                 n.div.by = NULL,
                                 n.required = FALSE,
                                 set.default.values = FALSE
                                 ), 
                         
                         rw2d = list(
                                 hyper = list(
                                         theta = list(
                                                 name = "precision",
                                                 short.name = "prec",
                                                 initial = 4,
                                                 fixed = FALSE,
                                                 prior = "loggamma",
                                                 param = c(1, 0.0001)
                                                 )
                                         ), 
                                 constr = TRUE,
                                 nrow.ncol = TRUE,
                                 augmented = FALSE,
                                 aug.factor = 1,
                                 aug.constr = NULL,
                                 n.div.by = NULL,
                                 n.required = FALSE,
                                 set.default.values = TRUE
                                 ), 

                         matern2d = list(
                                 hyper = list(
                                         theta1 = list(
                                                 name = "precision",
                                                 short.name = "prec",
                                                 initial = 4,
                                                 fixed = FALSE,
                                                 prior = "loggamma",
                                                 param = c(1, 0.0001)
                                                 ), 
                                         theta2 = list(
                                                 name = "range",
                                                 short.name = "range",
                                                 initial = 2,
                                                 fixed = FALSE,
                                                 prior = "loggamma",
                                                 param = c(1, 0.01)
                                                 )
                                         ), 
                                 constr = TRUE,
                                 nrow.ncol = TRUE,
                                 augmented = FALSE,
                                 aug.factor = 1,
                                 aug.constr = NULL,
                                 n.div.by = NULL,
                                 n.required = FALSE,
                                 set.default.values = TRUE
                                 ), 

                         copy = list(
                                 hyper = list(
                                         theta = list(
                                                 name = "beta",
                                                 short.name = "beta",
                                                 initial = 1,
                                                 fixed = TRUE,
                                                 prior = "normal",
                                                 param = c(1, 10)
                                                 )
                                         ), 
                                 constr = FALSE,
                                 nrow.ncol = FALSE,
                                 augmented = FALSE,
                                 aug.factor = 1,
                                 aug.constr = NULL,
                                 n.div.by = NULL,
                                 n.required = FALSE,
                                 set.default.values = FALSE
                                 )
                         ##
                         ),


                 ## models for group
                 group = list(

                         exchangeable = list(
                                 hyper = list(
                                         theta = list(
                                                 name = "correlation",
                                                 short.name = "rho",
                                                 initial = 1,
                                                 fixed = FALSE,
                                                 prior = "normal",
                                                 param = c(0, 0.2)
                                                 )
                                         )
                                 ),

                         ar1 = list(
                                 hyper = list(
                                         theta = list(
                                                 name = "correlation",
                                                 short.name = "rho",
                                                 initial = 2,
                                                 fixed = FALSE,
                                                 prior = "normal",
                                                 param = c(0, 0.2)
                                                 )
                                         )
                                 )
                         ), 

                 predictor = list(

                         predictor = list(
                                 hyper = list(
                                         theta = list(
                                                 name = "precision",
                                                 short.name = "prec",
                                                 initial = 11,
                                                 fixed = TRUE,
                                                 prior = "loggamma",
                                                 param = c(1, 0.00001)
                                                 )
                                         )
                                 )
                         ), 

                 hazard = list(

                         rw1 = list(
                                 hyper = list(
                                         theta = list(
                                                 name = "precision",
                                                 short.name = "prec",
                                                 initial = 4,
                                                 fixed = FALSE,
                                                 prior = "loggamma",
                                                 param = c(1, 0.0001)
                                                 )
                                         )
                                 ), 

                         rw2 = list(
                                 hyper = list(
                                         theta = list(
                                                 name = "precision",
                                                 short.name = "prec",
                                                 initial = 4,
                                                 fixed = FALSE,
                                                 prior = "loggamma",
                                                 param = c(1, 0.0001)
                                                 )
                                         )
                                 )
                         ), 
                         
                 ## likelihood models
                 likelihood = list(

                         poisson = list(
                                 hyper = list(
                                         ),
                                 survival = FALSE,
                                 discrete = TRUE
                                 ),

                         binomial = list(
                                 hyper = list(
                                         ),
                                 survival = FALSE,
                                 discrete = TRUE
                                 ), 

                         nbinomial = list(
                                 hyper = list(
                                         theta = list(
                                                 name = "size",
                                                 short.name = "size",
                                                 initial = log(10),
                                                 fixed = FALSE,
                                                 prior = "loggamma",
                                                 param = c(1, 100)
                                                 )
                                         ), 
                                 survival = FALSE,
                                 discrete = TRUE
                                 ),

                         exponential = list(
                                 hyper = list(
                                         ), 
                                 survival = TRUE,
                                 discrete = FALSE
                                 ),

                         coxph = list(
                                 hyper = list(
                                         ),
                                 survival = TRUE,
                                 discrete = TRUE
                                 ),

                         gaussian = list(
                                 hyper = list(
                                         theta = list(
                                                 name = "precision",
                                                 short.name = "prec",
                                                 initial = 4,
                                                 fixed = FALSE,
                                                 prior = "loggamma",
                                                 param = c(1, 0.0001)
                                                 )
                                         ), 
                                 survival = FALSE,
                                 discrete = FALSE
                                 ),

                         normal = list(
                                 hyper = list(
                                         theta = list(
                                                 name = "precision",
                                                 short.name = "prec",
                                                 initial = 4,
                                                 fixed = FALSE,
                                                 prior = "loggamma",
                                                 param = c(1, 0.0001)
                                                 )
                                         ), 
                                 survival = FALSE,
                                 discrete = FALSE
                                 ),

                         skewnormal = list(
                                 hyper = list(
                                         theta1 = list(
                                                 name = "inverse.scale",
                                                 short.name = "iscale",
                                                 initial = 4,
                                                 fixed = FALSE,
                                                 prior = "loggamma",
                                                 param = c(1, 0.0001)
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
                                 discrete = FALSE
                                 ),

                         sn = list(
                                 hyper = list(
                                         theta1 = list(
                                                 name = "inverse.scale",
                                                 short.name = "iscale",
                                                 initial = 4,
                                                 fixed = FALSE,
                                                 prior = "loggamma",
                                                 param = c(1, 0.0001)
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
                                 discrete = FALSE
                                 ),

                         gev = list(
                                 hyper = list(
                                         theta1 = list(
                                                 name = "precision",
                                                 short.name = "prec",
                                                 initial = 4,
                                                 fixed = FALSE,
                                                 prior = "loggamma",
                                                 param = c(1, 0.00001)
                                                 ),
                                         theta2 = list(
                                                 name = "GEVparameter",
                                                 short.name = "gev",
                                                 initial = 0,
                                                 fixed = FALSE, 
                                                 prior = "gaussian",
                                                 param = c(0, 6.25)
                                                 )
                                         ),
                                 survival = FALSE,
                                 discrete = FALSE
                                 ),

                         laplace = list(
                                 hyper = list(
                                         theta = list(
                                                 name = "precision",
                                                 short.name = "prec",
                                                 initial = 4,
                                                 fixed = FALSE,
                                                 prior = "loggamma",
                                                 param = c(1, 0.0001)
                                                 )
                                         ), 
                                 survival = FALSE,
                                 discrete = FALSE
                                 ),

                         weibull = list(
                                 hyper = list(
                                         theta = list(
                                                 name = "alpha",
                                                 short.name = "a",
                                                 initial = 4,
                                                 fixed = FALSE,
                                                 prior = "loggamma",
                                                 param = c(25, 25)
                                                 )
                                         ), 
                                 survival = TRUE,
                                 discrete = FALSE
                                 ),

                         weibullcure = list(
                                 hyper = list(
                                         theta1 = list(
                                                 name = "alpha",
                                                 short.name = "a",
                                                 initial = 4,
                                                 fixed = FALSE,
                                                 prior = "loggamma",
                                                 param = c(25, 25)
                                                 ), 
                                         theta2 = list(
                                                 name = "probability",
                                                 short.name = "p",
                                                 initial = 0,
                                                 fixed = FALSE,
                                                 prior = "gaussian",
                                                 param = c(0, 1)
                                                 )
                                         ), 
                                 survival = TRUE,
                                 discrete = FALSE
                                 ),

                         stochvol = list(
                                 hyper = list(
                                         ),
                                 survival = FALSE,
                                 discrete = FALSE
                                 ), 
                         
                         stochvolt = list(
                                 hyper = list(
                                         theta = list(
                                                 name = "degrees of freedom",
                                                 short.name = "dof",
                                                 initial = 4,
                                                 fixed = FALSE,
                                                 prior = "loggamma",
                                                 param = c(1, 0.5)
                                                 )
                                         ), 
                                 survival = FALSE,
                                 discrete = FALSE
                                 ),

                         stochvolnig = list(
                                 hyper = list(
                                         theta1 = list(
                                                 name = "skewness",
                                                 short.name = "skew",
                                                 initial = 0,
                                                 fixed = FALSE,
                                                 prior = "gaussian",
                                                 param = c(0, 10)
                                                 ),
                                         theta2 = list(
                                                 name = "shape",
                                                 short.name = "shape",
                                                 initial = 0,
                                                 fixed = FALSE,
                                                 prior = "loggamma",
                                                 param = c(1, 0.5)
                                                 )
                                         ), 
                                 survival = FALSE,
                                 discrete = FALSE
                                 ),

                         zeroinflatedpoisson0 = list(
                                 hyper = list(
                                         theta = list(
                                                 name = "probability",
                                                 short.name = "p",
                                                 initial = 0,
                                                 fixed = FALSE,
                                                 prior = "gaussian",
                                                 param = c(0, 1)
                                                 )
                                         ),
                                 survival = FALSE,
                                 discrete = FALSE
                                 ),

                         zeroinflatedpoisson1 = list(
                                 hyper = list(
                                         theta = list(
                                                 name = "probability",
                                                 short.name = "p",
                                                 initial = 0,
                                                 fixed = FALSE,
                                                 prior = "gaussian",
                                                 param = c(0, 1)
                                                 )
                                         ),
                                 survival = FALSE,
                                 discrete = FALSE
                                 ),

                         zeroinflatedpoisson2 = list(
                                 hyper = list(
                                         theta = list(
                                                 name = "probability",
                                                 short.name = "p",
                                                 initial = 0,
                                                 fixed = FALSE,
                                                 prior = "gaussian",
                                                 param = c(0, 1)
                                                 )
                                         ),
                                 survival = FALSE,
                                 discrete = FALSE
                                 ),

                         zeroinflatedbinomial0 = list(
                                 hyper = list(
                                         theta = list(
                                                 name = "probability",
                                                 short.name = "p",
                                                 initial = 0,
                                                 fixed = FALSE,
                                                 prior = "gaussian",
                                                 param = c(0, 1)
                                                 )
                                         ),
                                 survival = FALSE,
                                 discrete = FALSE
                                 ),
                         
                         zeroinflatedbinomial1 = list(
                                 hyper = list(
                                         theta = list(
                                                 name = "probability",
                                                 short.name = "p",
                                                 initial = 0,
                                                 fixed = FALSE,
                                                 prior = "gaussian",
                                                 param = c(0, 1)
                                                 )
                                         ),
                                 survival = FALSE,
                                 discrete = FALSE
                                 ),
                         
                         zeroinflatedbinomial2 = list(
                                 hyper = list(
                                         theta = list(
                                                 name = "probability",
                                                 short.name = "p",
                                                 initial = 0,
                                                 fixed = FALSE,
                                                 prior = "gaussian",
                                                 param = c(0, 1)
                                                 )
                                         ),
                                 survival = FALSE,
                                 discrete = FALSE
                                 ),
                         
                         zeroinflatedbetabinomial2 = list(
                                 hyper = list(
                                         theta1 = list(
                                                 name = "alpha",
                                                 short.name = "a",
                                                 initial = log(2),
                                                 fixed = FALSE,
                                                 prior = "gaussian",
                                                 param = c(0, 1)
                                                 ), 
                                         theta2 = list(
                                                 name = "beta",
                                                 short.name = "b",
                                                 initial = log(1),
                                                 fixed = FALSE,
                                                 prior = "gaussian",
                                                 param = c(0, 1)
                                                 )
                                         ),
                                 survival = FALSE,
                                 discrete = FALSE
                                 ),
                         
                         zeroinflatednbinomial0 = list(
                                 hyper = list(
                                         theta1 = list(
                                                 name = "size",
                                                 short.name = "size",
                                                 initial = log(10),
                                                 fixed = FALSE,
                                                 prior = "loggamma",
                                                 param = c(1, 0.01)
                                                 ), 
                                         theta2 = list(
                                                 name = "probability",
                                                 short.name = "p",
                                                 initial = -1,
                                                 fixed = FALSE,
                                                 prior = "gaussian",
                                                 param = c(0, 1)
                                                 )
                                         ),
                                 survival = FALSE,
                                 discrete = FALSE
                                 ),

                         zeroinflatednbinomial1 = list(
                                 hyper = list(
                                         theta1 = list(
                                                 name = "size",
                                                 short.name = "size",
                                                 initial = log(10),
                                                 fixed = FALSE,
                                                 prior = "loggamma",
                                                 param = c(1, 0.01)
                                                 ), 
                                         theta2 = list(
                                                 name = "probability",
                                                 short.name = "p",
                                                 initial = -1,
                                                 fixed = FALSE,
                                                 prior = "gaussian",
                                                 param = c(0, 1)
                                                 )
                                         ),
                                 survival = FALSE,
                                 discrete = FALSE
                                 ),
                         
                         zeroinflatednbinomial2 = list(
                                 hyper = list(
                                         theta1 = list(
                                                 name = "size",
                                                 short.name = "size",
                                                 initial = log(10),
                                                 fixed = FALSE,
                                                 prior = "loggamma",
                                                 param = c(1, 0.01)
                                                 ), 
                                         theta2 = list(
                                                 name = "alpha",
                                                 short.name = "a",
                                                 initial = log(2),
                                                 fixed = FALSE,
                                                 prior = "gaussian",
                                                 param = c(0, 1)
                                                 )
                                         ),
                                 survival = FALSE,
                                 discrete = FALSE
                                 ),
                         
                         t = list(
                                 hyper = list(
                                         theta1 = list(
                                                 name = "precision",
                                                 short.name = "prec",
                                                 initial = 4,
                                                 fixed = FALSE,
                                                 prior = "loggamma",
                                                 param = c(1, 0.0001)
                                                 ), 
                                         theta2 = list(
                                                 name = "degrees of freedom",
                                                 short.name = "dof",
                                                 initial = 4,
                                                 fixed = FALSE,
                                                 prior = "loggamma",
                                                 param = c(1, 0.5)
                                                 )
                                         ),
                                 survival = FALSE,
                                 discrete = FALSE
                                 ),
                         
                         logperiodogram = list(
                                 hyper = list(
                                         ),
                                 survival = FALSE,
                                 discrete = FALSE
                                 )
                         ),
                 
                 prior = list(
                         ## these are the same
                         normal = list(nparameters = 2),
                         gaussian = list(nparameters = 2),

                         wishart = list(nparameters = NA),
                         wishart1d = list(nparameters = 2),
                         wishart2d = list(nparameters = 4),
                         wishart3d = list(nparameters = 7),
                         loggamma = list(nparameters = 2),

                         ## these are the same
                         minuslogsqrtruncnormal = list(nparameters = 1),
                         logtnormal = list(nparameters = 1),
                         logtgaussian = list(nparameters = 1),

                         flat=list(nparameters = 0),
                         logflat=list(nparameters = 0),
                         logiflat=list(nparameters = 0), 

                         ## this is the 'no prior needed' prior
                         none = list(nparameters = 0)
                         ),

                 ## models that are defined through the wrapper.
                 wrapper = list(
                         joint = list(
                                 hyper = list(
                                         theta = list(
                                                 name = "precision",
                                                 short.name = "prec",
                                                 initial = 0,
                                                 fixed = TRUE,
                                                 prior = "loggamma",
                                                 param = c(1, 0.0001)
                                                 )
                                         ), 
                                 constr = FALSE,
                                 nrow.ncol = FALSE,
                                 augmented = FALSE,
                                 aug.factor = 1,
                                 aug.constr = NULL,
                                 n.div.by = NULL,
                                 n.required = FALSE,
                                 set.default.values = FALSE
                                 )
                         )
                 )
            )
}

`inla.is.model` = function(model, section = names(inla.models()), 
        stop.on.error = TRUE, ignore.case = FALSE)
{
    section = match.arg(section)
    models = names((inla.models()[ names(inla.models()) == section ])[[1]])

    if (is.character(model) && length(model) > 0) {
        if (ignore.case) {
            m = tolower(model)
            ms = tolower(models)
        } else {
            m = model
            ms = models
        }

        ret = c()
        for(i in 1:length(m)) {
        
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

    m = inla.model.properties.generic(inla.trim.family(model),
            (inla.models()[names(inla.models()) == section])[[1]],
            stop.on.error, ignore.case)

    if (is.null(m)) {
        return (NULL)
    }
    
    return (m)
}

`inla.model.properties.generic` = function(model, models, stop.on.error = TRUE, ignore.case = TRUE)
{
    ans = c()
    for(mm in model) {
        m = ifelse(ignore.case, tolower(mm), mm)
        if (ignore.case) {
            ms = tolower(names(models))
        } else {
            ms = names(models)
        }
        if (inla.is.generic(mm, stop.on.error, models, ignore.case)) {
            k = grep(paste("^", m, "$", sep=""), ms)
            ans = c(ans, list(models[[k]]))
        } else {
            ans = c(ans, list(NA))
        }
    }
    ## treat the case length(ans) == 1 specially. do not need a list
    ## of list then.
    if (length(ans) == 1) {
        return (ans[[1]])
    } else {
        return (ans)
    }
}

`inla.models.generate.Rd` = function(file = NULL)
{
    ## this function is used to generate the man-page for the models

    if (!is.null(file))
        sink(file)
    
    tab = "\t"
    tab2 = paste(tab, "\t",  sep="")
    tab3 = paste(tab2, "\t",  sep="")
    tab4 = paste(tab3, "\t",  sep="")
    tab5 = paste(tab4, "\t",  sep="")
    tab6 = paste(tab5, "\t",  sep="")
    tab7 = paste(tab6, "\t",  sep="")
    tab8 = paste(tab7, "\t",  sep="")
    

    cat("
\\name{inla.models}
\\alias{inla.models}
\\title{Valid models in INLA}
\\description{
This page describe the models implemented in \\code{inla}, divided into sections: ")
    cat(inla.paste(names(inla.models()), sep=", "), ".\n}\n")

    cat("
\\usage{
inla.models()
}\n")
    
    cat("\\value{
", tab, "Valid sections are: ")
    cat(inla.paste(names(inla.models()),  sep=", "),  "\n")
    cat("
", tab, "\\describe{
")

    for(section in names(inla.models())) {
        cat(tab2, "\\item{Section `", section, "'.}{\n",sep="")
        cat(tab3, "Valid models in this section are:\n")
        cat(tab3, "\\describe{\n")
        
        for(model in names(inla.models()[[section]])) {
            cat(tab4, "\\item{Model `", model, "'.}{",  sep="")

            if (section != "prior") {
                nhyper = length(inla.models()[[section]][[model]]$hyper)
                cat("Number of hyperparmeters are ", nhyper, ".\n", sep="")
                if (nhyper > 0) {
                    cat(tab5, "\\describe{\n")
                    
                    h = inla.models()[[section]][[model]]$hyper
                    nhyper = length(h)
                    if (nhyper > 0) {
                        for(nh in 1:nhyper) {
                            nm = names(h)[nh]
                            cat(tab6, "\\item{Hyperparameter `", nm, "'}{\n", sep="")
                            cat(tab7, "\\describe{\n")

                            for(m in names(h[[nm]])) {
                                mval = h[[nm]][[m]]
                                if (is.null(mval))
                                    mval = "NULL"
                                cat(tab7, "\\item{", m, " = }{`", inla.paste(mval), "'}\n", sep="")
                            }

                            cat(tab7, "}\n", sep="")
                            cat(tab6, "}\n", sep="")
                        }
                        cat(tab5, "}\n")
                    }

                    h = inla.models()[[section]][[model]]
                    ms = names(h)
                    ms = ms[ ms != "hyper" ]
                    
                    cat(tab5, "\\describe{\n")
                    cat(tab6, "\\item{Properties:}{\n", sep="")
                    cat(tab7, "\\describe{\n", sep="")
                    
                    for(m in ms) {
                        mval = h[[m]]
                        if (is.null(mval))
                            mval = "NULL"
                        cat(tab8, "\\item{", m, " = }{`", inla.paste(mval), "'}\n", sep="")
                    }
                    
                    cat(tab7, "}\n")
                    cat(tab6, "}\n")
                    cat(tab5, "}\n")
                }
            } else if (section == "prior") {
                np = inla.models()[[section]][[model]]$nparameters
                cat(tab5, "Number of parameters in the prior =", np, "\n")
            } else {
                stop("Should not happen.")
            }
                
            cat(tab4, "}\n")
        }
        cat(tab3, "}\n")
        cat(tab2, "}\n")
    }

    cat(tab, "} \n")
    cat("} \n")

    if (!is.null(file))
        sink()

}

`inla.models.generate.tex` = function(a.list = inla.models(), file = NULL)
{
    my.doit.recursively = function(a.list, level = 0L)
    {
        if (level == 0L)
            cat("\\begin{description}\n")

        tab = inla.paste(rep("\t", level + 1L))

        for(nm in names(a.list)) {
            if (is.list(a.list[[nm]])) {
                if (length(a.list[[nm]]) > 0) {
                    cat(tab, "\\item[", nm, "]\\ ", "\n", sep="")
                    cat(tab, "\\begin{description}\n")
                    my.doit.recursively(a.list[[nm]], level + 1L)
                    cat(tab, "\\end{description}\n")
                } else {
                    cat(tab, "\\item[", nm, "]\\ ", "\n", sep="")
                }
            } else {
                val = a.list[[nm]]
                cat(tab, "\\item[", nm, "]", as.character(enquote(a.list[[nm]]))[-1], "\n")
            }
        }

        if (level == 0L)
            cat("\\end{description}\n")
    }

    ##
    ##
    
    if (!is.null(file))
        sink(file)
    
    my.doit.recursively(a.list)
    
    if (!is.null(file))
        sink()
}
