## models and its hyperparameters are defined here

`inla.models` = function()
{
    ## this is not very clean solution, but for the moment is ok. the
    ## inla.models() function takes just to much time!!!
    
    envir = .GlobalEnv
    if (exists("...inla.hash.inla.models",  envir = envir) &&
        exists("...inla.hash.inla.models.hgid",  envir = envir) &&
        get("...inla.hash.inla.models.hgid", envir = envir) == inla.version(hgid = TRUE)) {

        return (get("...inla.hash.inla.models", envir = envir))

    } else {
        models = list(
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
                                aug.factor = 1,
                                aug.constr = NULL,
                                n.div.by = NULL,
                                n.required = FALSE,
                                set.default.values = FALSE
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
                                aug.factor = 1,
                                aug.constr = NULL,
                                n.div.by = NULL,
                                n.required = FALSE,
                                set.default.values = FALSE
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
                                aug.factor = 1,
                                aug.constr = NULL,
                                n.div.by = NULL,
                                n.required = FALSE,
                                set.default.values = FALSE
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
                                aug.factor = 2,
                                aug.constr = 1,
                                n.div.by = NULL,
                                n.required = FALSE,
                                set.default.values = FALSE
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
                                aug.factor = 1,
                                aug.constr = NULL,
                                n.div.by = NULL,
                                n.required = FALSE,
                                set.default.values = FALSE
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
                                aug.factor = 1,
                                aug.constr = NULL,
                                n.div.by = NULL,
                                n.required = TRUE,
                                set.default.values = TRUE
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
                                aug.factor = 1,
                                aug.constr = c(1, 2),
                                n.div.by = 2,
                                n.required = TRUE,
                                set.default.values = TRUE
                                ), 

                        bym = list(
                                hyper = list(
                                        theta1 = list(
                                                name = "log unstructured precision",
                                                short.name = "prec.unstruct",
                                                prior = "loggamma",
                                                param = c(1, 0.001), 
                                                initial = 4,
                                                fixed = FALSE,
                                                to.theta = function(x) log(x),
                                                from.theta = function(x) exp(x)
                                                ), 
                                        theta2 = list(
                                                name = "log spatial precision",
                                                short.name = "prec.spatial",
                                                prior = "normal",
                                                param = c(0, 0.00005), 
                                                initial = 4,
                                                fixed = FALSE,
                                                to.theta = function(x) log(x),
                                                from.theta = function(x) exp(x)
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
                                aug.factor = 1,
                                aug.constr = NULL,
                                n.div.by = NULL,
                                n.required = FALSE,
                                set.default.values = FALSE
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
                                aug.factor = 1,
                                aug.constr = NULL,
                                n.div.by = NULL,
                                n.required = TRUE,
                                set.default.values = TRUE
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
                                aug.factor = 1,
                                aug.constr = NULL,
                                n.div.by = NULL,
                                n.required = TRUE,
                                set.default.values = TRUE
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
                                                param = c(0, 0.001), 
                                                to.theta = function(x) log(x),
                                                from.theta = function(x) exp(x)
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
                                ## this will be redone anyway soon....
                                hyper = list(
                                        theta1 = list(
                                                name = "theta.T",
                                                short.name = "T",
                                                initial = 2,
                                                fixed = FALSE,
                                                prior = "normal",
                                                param = c(0, 1)
                                                ),
                                        theta2 = list(
                                                name = "theta.K",
                                                short.name = "K",
                                                initial = -2,
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
                                                ## corresponding wishart
                                                ## prior for the (a, b)
                                                ## gamma parameters
                                                param = c(2*1, 2*0.00005),
                                                to.theta = function(x) log(x),
                                                from.theta = function(x) exp(x)
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
                                aug.factor = 1,
                                aug.constr = c(1, 2),
                                n.div.by = 2,
                                n.required = TRUE,
                                set.default.values = TRUE
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
                                aug.factor = 1,
                                aug.constr = c(1, 2, 3),
                                n.div.by = 3,
                                n.required = TRUE,
                                set.default.values = TRUE
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
                                aug.factor = 1,
                                aug.constr = c(1, 2),
                                n.div.by = 2,
                                n.required = TRUE,
                                set.default.values = TRUE
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
                                                from.theta = function(x) exp(x)                                                                                                                          )
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
                                aug.factor = 1,
                                aug.constr = NULL,
                                n.div.by = NULL,
                                n.required = FALSE,
                                set.default.values = TRUE
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
                                                short.name = "b",
                                                initial = 1,
                                                fixed = TRUE,
                                                prior = "normal",
                                                param = c(1, 10),
                                                to.theta = function(x, low = -Inf, high = Inf) {
                                                    if (low == -Inf && high == Inf) {
                                                        return (x)
                                                    } else {
                                                        stopifnot((low != -Inf) &&  (high != Inf) && (low < high))
                                                        return (log( - (low - x)/(high -x))) 
                                                    }
                                                }, 
                                                from.theta = function(x, low = -Inf, high = Inf) {
                                                    if (low == -Inf && high == Inf) {
                                                        return (x)
                                                    } else {
                                                        stopifnot((low != -Inf) &&  (high != Inf) && (low < high))
                                                        return (low + exp(x)/(1+exp(x)) * (high - low))
                                                    }
                                                }
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
                                                name = "logit correlation",
                                                short.name = "rho",
                                                initial = 1,
                                                fixed = FALSE,
                                                prior = "normal",
                                                param = c(0, 0.2),
                                                to.theta = function(x, ngroup=2) log((1+x*(ngroup-1))/(1-x)), 
                                                from.theta = function(x, ngroup=2) (exp(x)-1)/(exp(x) + ngroup -1)
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
                                )
                        ), 

                predictor = list(

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
                        ), 

                hazard = list(

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
                                                param = c(1, 100), 
                                                to.theta = function(x) log(x), 
                                                from.theta = function(x) exp(x)
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
                                discrete = FALSE
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
                                discrete = FALSE
                                ),

                        sqrnormal = list(
                                hyper = list(
                                        theta1 = list(
                                                name = "log precision",
                                                short.name = "prec",
                                                initial = 4,
                                                fixed = FALSE,
                                                prior = "loggamma",
                                                param = c(1, 0.00005)
                                                ),
                                        theta2 = list(
                                                name = "scale",
                                                short.name = "scale",
                                                initial = 1,
                                                fixed = FALSE,
                                                prior = "gaussian",
                                                param = c(1, 100)
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
                                discrete = FALSE
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
                                                initial = 4,
                                                fixed = FALSE,
                                                prior = "gaussian",
                                                param = c(0, 10), 
                                                to.theta = function(x, shape.max = 1) log((1+x/shape.max)/(1-x/shape.max)), 
                                                from.theta = function(x, shape.max = 1) shape.max*(2*exp(x)/(1+exp(x))-1)
                                                )
                                        ), 
                                survival = FALSE,
                                discrete = FALSE
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
                                                param = c(0, 6.25),
                                                to.theta = function(x) x, 
                                                from.theta = function(x) x
                                                )
                                        ),
                                survival = FALSE,
                                discrete = FALSE
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
                                discrete = FALSE
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
                                discrete = FALSE
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
                                discrete = FALSE
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
                                discrete = FALSE
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
                                discrete = FALSE
                                ),

                        zeroinflatedpoisson2 = list(
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
                                discrete = FALSE
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
                                discrete = FALSE
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
                                discrete = FALSE
                                ),
                         
                        zeroinflatedbinomial2 = list(
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
                                discrete = FALSE
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
                                discrete = FALSE
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
                                discrete = FALSE
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
                                discrete = FALSE
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
                                discrete = FALSE
                                ),
                         
                        t = list(
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
                        none = list(nparameters = 0),

                        betacorrelation = list(nparameters = 2)
                        ),

                ## models that are defined through the wrapper.
                wrapper = list(
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
                                aug.factor = 1,
                                aug.constr = NULL,
                                n.div.by = NULL,
                                n.required = FALSE,
                                set.default.values = FALSE
                                )
                        )
                )

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

    
        ...inla.hash.inla.models <<- models
        ...inla.hash.inla.models.hgid <<- inla.version(hgid=TRUE)

        return (models)
    }

    stop("Should not happen")
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

    mm = inla.models()
    m = inla.model.properties.generic(inla.trim.family(model),
            (mm[names(mm) == section])[[1]],
            stop.on.error, ignore.case)

    if (is.null(m)) {
        return (NULL)
    }
    
    return (m)
}

`inla.model.properties.generic` = function(model, models, stop.on.error = TRUE, ignore.case = TRUE)
{
    m = ifelse(ignore.case, tolower(model), model)
    if (ignore.case) {
        ms = tolower(names(models))
    } else {
        ms = names(models)
    }
    k = grep(paste("^", m, "$", sep=""), ms)
    if (length(k) == 0L) {
        return (NULL)
    } else {
        return (models[[k]])
    }
}
