## this is how to define a must-be-enabled...
## status = "changed:Oct.25.2017",

`inla.special.number` <- function() {
    return (1048576.0)
}

`inla.models.section.latent` <- function() {
    list(
        latent =
            list(
                linear = list(
                    doc = "Alternative interface to an fixed effect",
                    hyper = list(),
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
                            hyperid = 1001,
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
                            hyperid = 2001,
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
                            hyperid = 2002,
                            name = "prec.u",
                            short.name = "prec",
                            prior = "loggamma",
                            param = c(1, 0.0001),
                            initial = log(1 / 0.0001),
                            fixed = TRUE,
                            to.theta = function(x) log(x),
                            from.theta = function(x) exp(x)
                        ),
                        theta3 = list(
                            hyperid = 2003,
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
                            hyperid = 2004,
                            name = "prec.x",
                            short.name = "prec.x",
                            prior = "loggamma",
                            param = c(1, 10000),
                            initial = log(1 / 10000),
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
                            hyperid = 3001,
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
                            hyperid = 3002,
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
                    doc = "Generic latent model specified using R",
                    hyper = list(),
                    constr = FALSE,
                    nrow.ncol = FALSE,
                    augmented = FALSE,
                    aug.factor = 1L,
                    aug.constr = NULL,
                    n.div.by = NULL,
                    n.required = TRUE,
                    set.default.values = TRUE,
                    pdf = "rgeneric"
                ),

                cgeneric = list(
                    doc = "Generic latent model specified using C",
                    hyper = list(),
                    constr = FALSE,
                    nrow.ncol = FALSE,
                    augmented = FALSE,
                    aug.factor = 1L,
                    aug.constr = NULL,
                    n.div.by = NULL,
                    n.required = TRUE,
                    set.default.values = TRUE,
                    pdf = "rgeneric" ## Yes
                ),

                rw1 = list(
                    doc = "Random walk of order 1",
                    hyper = list(
                        theta = list(
                            hyperid = 4001,
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
                            hyperid = 5001,
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
                            hyperid = 6001,
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
                            hyperid = 7001,
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
                            hyperid = 8001,
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
                            hyperid = 9001,
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
                            hyperid = 9002,
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
                    constr = TRUE,
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
                            hyperid = 10001,
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
                            hyperid = 10002,
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
                            hyperid = 11001,
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
                            hyperid = 11002,
                            name = "logit phi",
                            short.name = "phi",
                            prior = "pc",
                            param = c(0.5, 0.5),
                            initial = -3,
                            fixed = FALSE,
                            to.theta = function(x) log(x / (1 - x)),
                            from.theta = function(x) exp(x) / (1 + exp(x))
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
                    pdf = "bym2"
                ),

                besagproper = list(
                    doc = "A proper version of the Besag model",
                    hyper = list(
                        theta1 = list(
                            hyperid = 12001,
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
                            hyperid = 12002,
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
                    pdf = "besagproper"
                ),

                besagproper2 = list(
                    doc = "An alternative proper version of the Besag model",
                    hyper = list(
                        theta1 = list(
                            hyperid = 13001,
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
                            hyperid = 13002,
                            name = "logit lambda",
                            short.name = "lambda",
                            prior = "gaussian",
                            param = c(0, 0.45),
                            initial = 3,
                            fixed = FALSE,
                            to.theta = function(x) log(x / (1 - x)),
                            from.theta = function(x) exp(x) / (1 + exp(x))
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
                    pdf = "besagproper2"
                ),

                fgn = list(
                    doc = "Fractional Gaussian noise model",
                    hyper = list(
                        theta1 = list(
                            hyperid = 13101,
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
                            hyperid = 13102,
                            name = "logit H",
                            short.name = "H",
                            prior = "pcfgnh",
                            param = c(0.9, 0.1),
                            initial = 2,
                            fixed = FALSE,
                            to.theta = function(x) log((2 * x - 1) / (2 * (1 - x))),
                            from.theta = function(x) 0.5 + 0.5 * exp(x) / (1 + exp(x))
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
                    order.default = 4L, ## default order for approximation
                    order.defined = 3L:4L, ## the list of orders which are implemented
                    pdf = "fgn"
                ),

                fgn2 = list(
                    doc = "Fractional Gaussian noise model (alt 2)",
                    hyper = list(
                        theta1 = list(
                            hyperid = 13111,
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
                            hyperid = 13112,
                            name = "logit H",
                            short.name = "H",
                            prior = "pcfgnh",
                            param = c(0.9, 0.1),
                            initial = 2,
                            fixed = FALSE,
                            to.theta = function(x) log((2 * x - 1) / (2 * (1 - x))),
                            from.theta = function(x) 0.5 + 0.5 * exp(x) / (1 + exp(x))
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
                    order.default = 4L, ## default order for approximation
                    order.defined = 3L:4L, ## the list of orders which are implemented
                    pdf = "fgn"
                ),

                ar1 = list(
                    doc = "Auto-regressive model of order 1 (AR(1))",
                    hyper = list(
                        theta1 = list(
                            hyperid = 14001,
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
                            hyperid = 14002,
                            name = "logit lag one correlation",
                            short.name = "rho",
                            prior = "normal",
                            param = c(0, 0.15),
                            initial = 2,
                            fixed = FALSE,
                            to.theta = function(x) log((1 + x) / (1 - x)),
                            from.theta = function(x) 2 * exp(x) / (1 + exp(x)) - 1
                        ),
                        theta3 = list(
                            hyperid = 14003,
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
                            hyperid = 14101,
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
                            hyperid = 14102,
                            name = "logit lag one correlation",
                            short.name = "rho",
                            prior = "pc.cor0",
                            param = c(0.5, 0.5),
                            initial = 2,
                            fixed = FALSE,
                            to.theta = function(x) log((1 + x) / (1 - x)),
                            from.theta = function(x) 2 * exp(x) / (1 + exp(x)) - 1
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
                    pdf = "ar1c"
                ),

                ar = list(
                    doc = "Auto-regressive model of order p (AR(p))",
                    ## to many parameters here, but ...
                    hyper = list(
                        theta1 = list(
                            hyperid = 15001,
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
                            hyperid = 15002,
                            name = "pacf1",
                            short.name = "pacf1",
                            initial = 1,
                            fixed = FALSE,
                            prior = "pc.cor0",
                            param = c(0.5, 0.5),
                            to.theta = function(x) log((1 + x) / (1 - x)),
                            from.theta = function(x) 2 * exp(x) / (1 + exp(x)) - 1
                        ),
                        theta3 = list(
                            hyperid = 15003,
                            name = "pacf2",
                            short.name = "pacf2",
                            initial = 0,
                            fixed = FALSE,
                            prior = "pc.cor0",
                            param = c(0.5, 0.4),
                            to.theta = function(x) log((1 + x) / (1 - x)),
                            from.theta = function(x) 2 * exp(x) / (1 + exp(x)) - 1
                        ),
                        theta4 = list(
                            hyperid = 15004,
                            name = "pacf3",
                            short.name = "pacf3",
                            initial = 0,
                            fixed = FALSE,
                            prior = "pc.cor0",
                            param = c(0.5, 0.3),
                            to.theta = function(x) log((1 + x) / (1 - x)),
                            from.theta = function(x) 2 * exp(x) / (1 + exp(x)) - 1
                        ),
                        theta5 = list(
                            hyperid = 15005,
                            name = "pacf4",
                            short.name = "pacf4",
                            initial = 0,
                            fixed = FALSE,
                            prior = "pc.cor0",
                            param = c(0.5, 0.2),
                            to.theta = function(x) log((1 + x) / (1 - x)),
                            from.theta = function(x) 2 * exp(x) / (1 + exp(x)) - 1
                        ),
                        theta6 = list(
                            hyperid = 15006,
                            name = "pacf5",
                            short.name = "pacf5",
                            initial = 0,
                            fixed = FALSE,
                            prior = "pc.cor0",
                            param = c(0.5, 0.1),
                            to.theta = function(x) log((1 + x) / (1 - x)),
                            from.theta = function(x) 2 * exp(x) / (1 + exp(x)) - 1
                        ),
                        theta7 = list(
                            hyperid = 15007,
                            name = "pacf6",
                            short.name = "pacf6",
                            initial = 0,
                            fixed = FALSE,
                            prior = "pc.cor0",
                            param = c(0.5, 0.1),
                            to.theta = function(x) log((1 + x) / (1 - x)),
                            from.theta = function(x) 2 * exp(x) / (1 + exp(x)) - 1
                        ),
                        theta8 = list(
                            hyperid = 15008,
                            name = "pacf7",
                            short.name = "pacf7",
                            initial = 0,
                            fixed = FALSE,
                            prior = "pc.cor0",
                            param = c(0.5, 0.1),
                            to.theta = function(x) log((1 + x) / (1 - x)),
                            from.theta = function(x) 2 * exp(x) / (1 + exp(x)) - 1
                        ),
                        theta9 = list(
                            hyperid = 15009,
                            name = "pacf8",
                            short.name = "pacf8",
                            initial = 0,
                            fixed = FALSE,
                            prior = "pc.cor0",
                            param = c(0.5, 0.1),
                            to.theta = function(x) log((1 + x) / (1 - x)),
                            from.theta = function(x) 2 * exp(x) / (1 + exp(x)) - 1
                        ),
                        theta10 = list(
                            hyperid = 15010,
                            name = "pacf9",
                            short.name = "pacf9",
                            initial = 0,
                            fixed = FALSE,
                            prior = "pc.cor0",
                            param = c(0.5, 0.1),
                            to.theta = function(x) log((1 + x) / (1 - x)),
                            from.theta = function(x) 2 * exp(x) / (1 + exp(x)) - 1
                        ),
                        theta11 = list(
                            hyperid = 15011,
                            name = "pacf10",
                            short.name = "pacf10",
                            initial = 0,
                            fixed = FALSE,
                            prior = "pc.cor0",
                            param = c(0.5, 0.1),
                            to.theta = function(x) log((1 + x) / (1 - x)),
                            from.theta = function(x) 2 * exp(x) / (1 + exp(x)) - 1
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
                            hyperid = 16001,
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
                            hyperid = 16002,
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
                            hyperid = 16101,
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
                            hyperid = 16102,
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
                            hyperid = 16103,
                            name = "logit correlation",
                            short.name = "cor",
                            initial = 4,
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(),
                            to.theta = function(x) log((1 + x) / (1 - x)),
                            from.theta = function(x) 2 * exp(x) / (1 + exp(x)) - 1
                        ),
                        theta4 = list(
                            hyperid = 16104,
                            name = "gamma1",
                            short.name = "g1",
                            initial = 1,
                            fixed = TRUE,
                            prior = "normal",
                            param = c(1, 36),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta5 = list(
                            hyperid = 16105,
                            name = "gamma2",
                            short.name = "g2",
                            initial = 1,
                            fixed = TRUE,
                            prior = "normal",
                            param = c(1, 36),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta6 = list(
                            hyperid = 16106,
                            name = "gamma3",
                            short.name = "g3",
                            initial = 1,
                            fixed = TRUE,
                            prior = "normal",
                            param = c(1, 36),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta7 = list(
                            hyperid = 16107,
                            name = "gamma4",
                            short.name = "g4",
                            initial = 1,
                            fixed = TRUE,
                            prior = "normal",
                            param = c(1, 36),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta8 = list(
                            hyperid = 16108,
                            name = "gamma5",
                            short.name = "g5",
                            initial = 1,
                            fixed = TRUE,
                            prior = "normal",
                            param = c(1, 36),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta9 = list(
                            hyperid = 16109,
                            name = "gamma6",
                            short.name = "g6",
                            initial = 1,
                            fixed = TRUE,
                            prior = "normal",
                            param = c(1, 36),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta10 = list(
                            hyperid = 16110,
                            name = "gamma7",
                            short.name = "g7",
                            initial = 1,
                            fixed = TRUE,
                            prior = "normal",
                            param = c(1, 36),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta11 = list(
                            hyperid = 16111,
                            name = "gamma8",
                            short.name = "g8",
                            initial = 1,
                            fixed = TRUE,
                            prior = "normal",
                            param = c(1, 36),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta12 = list(
                            hyperid = 16112,
                            name = "gamma9",
                            short.name = "g9",
                            initial = 1,
                            fixed = TRUE,
                            prior = "normal",
                            param = c(1, 36),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta13 = list(
                            hyperid = 16113,
                            name = "gamma10",
                            short.name = "g10",
                            initial = 1,
                            fixed = TRUE,
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
                    pdf = "intslope"
                ),

                generic = list(
                    doc = "A generic model",
                    hyper = list(
                        theta = list(
                            hyperid = 17001,
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
                            hyperid = 18001,
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
                            hyperid = 19001,
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
                            hyperid = 19002,
                            name = "beta",
                            short.name = "beta",
                            initial = 2,
                            fixed = FALSE,
                            prior = "gaussian",
                            param = c(0, 0.1),
                            to.theta = function(x) log(x / (1 - x)),
                            from.theta = function(x) exp(x) / (1 + exp(x))
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
                            hyperid = 20001,
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
                            hyperid = 20002,
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
                            hyperid = 21001,
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
                            hyperid = 21002,
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
                            hyperid = 21003,
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
                            hyperid = 21004,
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
                            hyperid = 21005,
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
                            hyperid = 21006,
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
                            hyperid = 21007,
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
                            hyperid = 21008,
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
                            hyperid = 21009,
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
                            hyperid = 21010,
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
                            hyperid = 21011,
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
                    pdf = "generic3"
                ),

                spde = list(
                    doc = "A SPDE model",
                    ## this will be redone anyway soon....
                    hyper = list(
                        theta1 = list(
                            hyperid = 22001,
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
                            hyperid = 22002,
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
                            hyperid = 22003,
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
                            hyperid = 22004,
                            name = "theta.OC",
                            short.name = "OC",
                            initial = -20,
                            fixed = TRUE,
                            prior = "normal",
                            param = c(0, 0.2),
                            to.theta = function(x) log(x / (1 - x)),
                            from.theta = function(x) exp(x) / (1 + exp(x))
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
                            hyperid = 23001,
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
                            hyperid = 23002,
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
                            hyperid = 23003,
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
                            hyperid = 23004,
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
                            hyperid = 23005,
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
                            hyperid = 23006,
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
                            hyperid = 23007,
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
                            hyperid = 23008,
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
                            hyperid = 23009,
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
                            hyperid = 23010,
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
                            hyperid = 23011,
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
                            hyperid = 23012,
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
                            hyperid = 23013,
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
                            hyperid = 23014,
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
                            hyperid = 23015,
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
                            hyperid = 23016,
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
                            hyperid = 23017,
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
                            hyperid = 23018,
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
                            hyperid = 23019,
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
                            hyperid = 23020,
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
                            hyperid = 23021,
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
                            hyperid = 23022,
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
                            hyperid = 23023,
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
                            hyperid = 23024,
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
                            hyperid = 23025,
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
                            hyperid = 23026,
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
                            hyperid = 23027,
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
                            hyperid = 23028,
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
                            hyperid = 23029,
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
                            hyperid = 23030,
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
                            hyperid = 23031,
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
                            hyperid = 23032,
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
                            hyperid = 23033,
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
                            hyperid = 23034,
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
                            hyperid = 23035,
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
                            hyperid = 23036,
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
                            hyperid = 23037,
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
                            hyperid = 23038,
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
                            hyperid = 23039,
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
                            hyperid = 23040,
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
                            hyperid = 23041,
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
                            hyperid = 23042,
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
                            hyperid = 23043,
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
                            hyperid = 23044,
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
                            hyperid = 23045,
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
                            hyperid = 23046,
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
                            hyperid = 23047,
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
                            hyperid = 23048,
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
                            hyperid = 23049,
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
                            hyperid = 23050,
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
                            hyperid = 23051,
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
                            hyperid = 23052,
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
                            hyperid = 23053,
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
                            hyperid = 23054,
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
                            hyperid = 23055,
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
                            hyperid = 23056,
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
                            hyperid = 23057,
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
                            hyperid = 23058,
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
                            hyperid = 23059,
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
                            hyperid = 23060,
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
                            hyperid = 23061,
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
                            hyperid = 23062,
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
                            hyperid = 23063,
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
                            hyperid = 23064,
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
                            hyperid = 23065,
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
                            hyperid = 23066,
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
                            hyperid = 23067,
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
                            hyperid = 23068,
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
                            hyperid = 23069,
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
                            hyperid = 23070,
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
                            hyperid = 23071,
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
                            hyperid = 23072,
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
                            hyperid = 23073,
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
                            hyperid = 23074,
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
                            hyperid = 23075,
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
                            hyperid = 23076,
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
                            hyperid = 23077,
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
                            hyperid = 23078,
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
                            hyperid = 23079,
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
                            hyperid = 23080,
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
                            hyperid = 23081,
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
                            hyperid = 23082,
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
                            hyperid = 23083,
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
                            hyperid = 23084,
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
                            hyperid = 23085,
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
                            hyperid = 23086,
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
                            hyperid = 23087,
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
                            hyperid = 23088,
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
                            hyperid = 23089,
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
                            hyperid = 23090,
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
                            hyperid = 23091,
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
                            hyperid = 23092,
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
                            hyperid = 23093,
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
                            hyperid = 23094,
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
                            hyperid = 23095,
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
                            hyperid = 23096,
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
                            hyperid = 23097,
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
                            hyperid = 23098,
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
                            hyperid = 23099,
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
                            hyperid = 23100,
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
                            hyperid = 24001,
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
                            hyperid = 24002,
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
                            hyperid = 24003,
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
                            hyperid = 24004,
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
                            hyperid = 24005,
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
                            hyperid = 24006,
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
                            hyperid = 24007,
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
                            hyperid = 24008,
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
                            hyperid = 24009,
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
                            hyperid = 24010,
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
                            hyperid = 24011,
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
                            hyperid = 24012,
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
                            hyperid = 24013,
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
                            hyperid = 24014,
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
                            hyperid = 24015,
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
                            hyperid = 24016,
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
                            hyperid = 24017,
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
                            hyperid = 24018,
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
                            hyperid = 24019,
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
                            hyperid = 24020,
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
                            hyperid = 24021,
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
                            hyperid = 24022,
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
                            hyperid = 24023,
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
                            hyperid = 24024,
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
                            hyperid = 24025,
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
                            hyperid = 24026,
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
                            hyperid = 24027,
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
                            hyperid = 24028,
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
                            hyperid = 24029,
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
                            hyperid = 24030,
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
                            hyperid = 24031,
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
                            hyperid = 24032,
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
                            hyperid = 24033,
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
                            hyperid = 24034,
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
                            hyperid = 24035,
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
                            hyperid = 24036,
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
                            hyperid = 24037,
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
                            hyperid = 24038,
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
                            hyperid = 24039,
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
                            hyperid = 24040,
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
                            hyperid = 24041,
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
                            hyperid = 24042,
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
                            hyperid = 24043,
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
                            hyperid = 24044,
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
                            hyperid = 24045,
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
                            hyperid = 24046,
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
                            hyperid = 24047,
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
                            hyperid = 24048,
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
                            hyperid = 24049,
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
                            hyperid = 24050,
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
                            hyperid = 24051,
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
                            hyperid = 24052,
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
                            hyperid = 24053,
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
                            hyperid = 24054,
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
                            hyperid = 24055,
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
                            hyperid = 24056,
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
                            hyperid = 24057,
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
                            hyperid = 24058,
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
                            hyperid = 24059,
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
                            hyperid = 24060,
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
                            hyperid = 24061,
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
                            hyperid = 24062,
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
                            hyperid = 24063,
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
                            hyperid = 24064,
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
                            hyperid = 24065,
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
                            hyperid = 24066,
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
                            hyperid = 24067,
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
                            hyperid = 24068,
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
                            hyperid = 24069,
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
                            hyperid = 24070,
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
                            hyperid = 24071,
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
                            hyperid = 24072,
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
                            hyperid = 24073,
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
                            hyperid = 24074,
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
                            hyperid = 24075,
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
                            hyperid = 24076,
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
                            hyperid = 24077,
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
                            hyperid = 24078,
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
                            hyperid = 24079,
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
                            hyperid = 24080,
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
                            hyperid = 24081,
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
                            hyperid = 24082,
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
                            hyperid = 24083,
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
                            hyperid = 24084,
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
                            hyperid = 24085,
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
                            hyperid = 24086,
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
                            hyperid = 24087,
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
                            hyperid = 24088,
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
                            hyperid = 24089,
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
                            hyperid = 24090,
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
                            hyperid = 24091,
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
                            hyperid = 24092,
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
                            hyperid = 24093,
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
                            hyperid = 24094,
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
                            hyperid = 24095,
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
                            hyperid = 24096,
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
                            hyperid = 24097,
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
                            hyperid = 24098,
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
                            hyperid = 24099,
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
                            hyperid = 24100,
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
                            hyperid = 25001,
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
                            param = c(2 * 1, 2 * 0.00005),
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
                            hyperid = 26001,
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
                            hyperid = 26002,
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
                            hyperid = 26003,
                            name = "logit correlation",
                            short.name = "cor",
                            initial = 4,
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(),
                            to.theta = function(x) log((1 + x) / (1 - x)),
                            from.theta = function(x) 2 * exp(x) / (1 + exp(x)) - 1
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
                            hyperid = 27001,
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
                            hyperid = 27002,
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
                            hyperid = 27003,
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
                            hyperid = 27004,
                            name = "logit correlation12",
                            short.name = "cor12",
                            initial = 0,
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(),
                            to.theta = function(x) log((1 + x) / (1 - x)),
                            from.theta = function(x) 2 * exp(x) / (1 + exp(x)) - 1
                        ),
                        theta5 = list(
                            hyperid = 27005,
                            name = "logit correlation13",
                            short.name = "cor13",
                            initial = 0,
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(),
                            to.theta = function(x) log((1 + x) / (1 - x)),
                            from.theta = function(x) 2 * exp(x) / (1 + exp(x)) - 1
                        ),
                        theta6 = list(
                            hyperid = 27006,
                            name = "logit correlation23",
                            short.name = "cor23",
                            initial = 0,
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(),
                            to.theta = function(x) log((1 + x) / (1 - x)),
                            from.theta = function(x) 2 * exp(x) / (1 + exp(x)) - 1
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
                            hyperid = 28001,
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
                            hyperid = 28002,
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
                            hyperid = 28003,
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
                            hyperid = 28004,
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
                            hyperid = 28005,
                            name = "logit correlation12",
                            short.name = "cor12",
                            initial = 0,
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(),
                            to.theta = function(x) log((1 + x) / (1 - x)),
                            from.theta = function(x) 2 * exp(x) / (1 + exp(x)) - 1
                        ),
                        theta6 = list(
                            hyperid = 28006,
                            name = "logit correlation13",
                            short.name = "cor13",
                            initial = 0,
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(),
                            to.theta = function(x) log((1 + x) / (1 - x)),
                            from.theta = function(x) 2 * exp(x) / (1 + exp(x)) - 1
                        ),
                        theta7 = list(
                            hyperid = 28007,
                            name = "logit correlation14",
                            short.name = "cor14",
                            initial = 0,
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(),
                            to.theta = function(x) log((1 + x) / (1 - x)),
                            from.theta = function(x) 2 * exp(x) / (1 + exp(x)) - 1
                        ),
                        theta8 = list(
                            hyperid = 28008,
                            name = "logit correlation23",
                            short.name = "cor23",
                            initial = 0,
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(),
                            to.theta = function(x) log((1 + x) / (1 - x)),
                            from.theta = function(x) 2 * exp(x) / (1 + exp(x)) - 1
                        ),
                        theta9 = list(
                            hyperid = 28009,
                            name = "logit correlation24",
                            short.name = "cor24",
                            initial = 0,
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(),
                            to.theta = function(x) log((1 + x) / (1 - x)),
                            from.theta = function(x) 2 * exp(x) / (1 + exp(x)) - 1
                        ),
                        theta10 = list(
                            hyperid = 28010,
                            name = "logit correlation34",
                            short.name = "cor34",
                            initial = 0,
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(),
                            to.theta = function(x) log((1 + x) / (1 - x)),
                            from.theta = function(x) 2 * exp(x) / (1 + exp(x)) - 1
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
                            hyperid = 29001,
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
                            hyperid = 29002,
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
                            hyperid = 29003,
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
                            hyperid = 29004,
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
                            hyperid = 29005,
                            name = "log precision5",
                            short.name = "prec5",
                            initial = 4,
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(),
                            to.theta = function(x) log(x),
                            from.theta = function(x) exp(x)
                        ),
                        theta6 = list(
                            hyperid = 29006,
                            name = "logit correlation12",
                            short.name = "cor12",
                            initial = 0,
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(),
                            to.theta = function(x) log((1 + x) / (1 - x)),
                            from.theta = function(x) 2 * exp(x) / (1 + exp(x)) - 1
                        ),
                        theta7 = list(
                            hyperid = 29007,
                            name = "logit correlation13",
                            short.name = "cor13",
                            initial = 0,
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(),
                            to.theta = function(x) log((1 + x) / (1 - x)),
                            from.theta = function(x) 2 * exp(x) / (1 + exp(x)) - 1
                        ),
                        theta8 = list(
                            hyperid = 29008,
                            name = "logit correlation14",
                            short.name = "cor14",
                            initial = 0,
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(),
                            to.theta = function(x) log((1 + x) / (1 - x)),
                            from.theta = function(x) 2 * exp(x) / (1 + exp(x)) - 1
                        ),
                        theta9 = list(
                            hyperid = 29009,
                            name = "logit correlation15",
                            short.name = "cor15",
                            initial = 0,
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(),
                            to.theta = function(x) log((1 + x) / (1 - x)),
                            from.theta = function(x) 2 * exp(x) / (1 + exp(x)) - 1
                        ),
                        theta10 = list(
                            hyperid = 29010,
                            name = "logit correlation23",
                            short.name = "cor23",
                            initial = 0,
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(),
                            to.theta = function(x) log((1 + x) / (1 - x)),
                            from.theta = function(x) 2 * exp(x) / (1 + exp(x)) - 1
                        ),
                        theta11 = list(
                            hyperid = 29011,
                            name = "logit correlation24",
                            short.name = "cor24",
                            initial = 0,
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(),
                            to.theta = function(x) log((1 + x) / (1 - x)),
                            from.theta = function(x) 2 * exp(x) / (1 + exp(x)) - 1
                        ),
                        theta12 = list(
                            hyperid = 29012,
                            name = "logit correlation25",
                            short.name = "cor25",
                            initial = 0,
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(),
                            to.theta = function(x) log((1 + x) / (1 - x)),
                            from.theta = function(x) 2 * exp(x) / (1 + exp(x)) - 1
                        ),
                        theta13 = list(
                            hyperid = 29013,
                            name = "logit correlation34",
                            short.name = "cor34",
                            initial = 0,
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(),
                            to.theta = function(x) log((1 + x) / (1 - x)),
                            from.theta = function(x) 2 * exp(x) / (1 + exp(x)) - 1
                        ),
                        theta14 = list(
                            hyperid = 29014,
                            name = "logit correlation35",
                            short.name = "cor35",
                            initial = 0,
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(),
                            to.theta = function(x) log((1 + x) / (1 - x)),
                            from.theta = function(x) 2 * exp(x) / (1 + exp(x)) - 1
                        ),
                        theta15 = list(
                            hyperid = 29015,
                            name = "logit correlation45",
                            short.name = "cor45",
                            initial = 0,
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(),
                            to.theta = function(x) log((1 + x) / (1 - x)),
                            from.theta = function(x) 2 * exp(x) / (1 + exp(x)) - 1
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

                iidkd = list(
                    doc = "Gaussian random effect in dim=k with Wishart prior",
                    hyper = list(
                        theta1 = list(
                            hyperid = 29101,
                            name = "theta1",
                            short.name = "theta1",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "wishartkd",
                            param = c(24, rep(inla.special.number(), (24*25)/2)), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta2 = list(
                            hyperid = 29102,
                            name = "theta2",
                            short.name = "theta2",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta3 = list(
                            hyperid = 29103,
                            name = "theta3",
                            short.name = "theta3",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta4 = list(
                            hyperid = 29104,
                            name = "theta4",
                            short.name = "theta4",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta5 = list(
                            hyperid = 29105,
                            name = "theta5",
                            short.name = "theta5",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta6 = list(
                            hyperid = 29106,
                            name = "theta6",
                            short.name = "theta6",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta7 = list(
                            hyperid = 29107,
                            name = "theta7",
                            short.name = "theta7",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta8 = list(
                            hyperid = 29108,
                            name = "theta8",
                            short.name = "theta8",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta9 = list(
                            hyperid = 29109,
                            name = "theta9",
                            short.name = "theta9",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta10 = list(
                            hyperid = 29110,
                            name = "theta10",
                            short.name = "theta10",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta11 = list(
                            hyperid = 29111,
                            name = "theta11",
                            short.name = "theta11",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta12 = list(
                            hyperid = 29112,
                            name = "theta12",
                            short.name = "theta12",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta13 = list(
                            hyperid = 29113,
                            name = "theta13",
                            short.name = "theta13",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta14 = list(
                            hyperid = 29114,
                            name = "theta14",
                            short.name = "theta14",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta15 = list(
                            hyperid = 29115,
                            name = "theta15",
                            short.name = "theta15",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta16 = list(
                            hyperid = 29116,
                            name = "theta16",
                            short.name = "theta16",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta17 = list(
                            hyperid = 29117,
                            name = "theta17",
                            short.name = "theta17",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta18 = list(
                            hyperid = 29118,
                            name = "theta18",
                            short.name = "theta18",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta19 = list(
                            hyperid = 29119,
                            name = "theta19",
                            short.name = "theta19",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta20 = list(
                            hyperid = 29120,
                            name = "theta20",
                            short.name = "theta20",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta21 = list(
                            hyperid = 29121,
                            name = "theta21",
                            short.name = "theta21",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta22 = list(
                            hyperid = 29122,
                            name = "theta22",
                            short.name = "theta22",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta23 = list(
                            hyperid = 29123,
                            name = "theta23",
                            short.name = "theta23",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta24 = list(
                            hyperid = 29124,
                            name = "theta24",
                            short.name = "theta24",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta25 = list(
                            hyperid = 29125,
                            name = "theta25",
                            short.name = "theta25",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta26 = list(
                            hyperid = 29126,
                            name = "theta26",
                            short.name = "theta26",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta27 = list(
                            hyperid = 29127,
                            name = "theta27",
                            short.name = "theta27",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta28 = list(
                            hyperid = 29128,
                            name = "theta28",
                            short.name = "theta28",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta29 = list(
                            hyperid = 29129,
                            name = "theta29",
                            short.name = "theta29",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta30 = list(
                            hyperid = 29130,
                            name = "theta30",
                            short.name = "theta30",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta31 = list(
                            hyperid = 29131,
                            name = "theta31",
                            short.name = "theta31",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta32 = list(
                            hyperid = 29132,
                            name = "theta32",
                            short.name = "theta32",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta33 = list(
                            hyperid = 29133,
                            name = "theta33",
                            short.name = "theta33",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta34 = list(
                            hyperid = 29134,
                            name = "theta34",
                            short.name = "theta34",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta35 = list(
                            hyperid = 29135,
                            name = "theta35",
                            short.name = "theta35",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta36 = list(
                            hyperid = 29136,
                            name = "theta36",
                            short.name = "theta36",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta37 = list(
                            hyperid = 29137,
                            name = "theta37",
                            short.name = "theta37",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta38 = list(
                            hyperid = 29138,
                            name = "theta38",
                            short.name = "theta38",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta39 = list(
                            hyperid = 29139,
                            name = "theta39",
                            short.name = "theta39",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta40 = list(
                            hyperid = 29140,
                            name = "theta40",
                            short.name = "theta40",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta41 = list(
                            hyperid = 29141,
                            name = "theta41",
                            short.name = "theta41",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta42 = list(
                            hyperid = 29142,
                            name = "theta42",
                            short.name = "theta42",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta43 = list(
                            hyperid = 29143,
                            name = "theta43",
                            short.name = "theta43",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta44 = list(
                            hyperid = 29144,
                            name = "theta44",
                            short.name = "theta44",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta45 = list(
                            hyperid = 29145,
                            name = "theta45",
                            short.name = "theta45",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta46 = list(
                            hyperid = 29146,
                            name = "theta46",
                            short.name = "theta46",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta47 = list(
                            hyperid = 29147,
                            name = "theta47",
                            short.name = "theta47",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta48 = list(
                            hyperid = 29148,
                            name = "theta48",
                            short.name = "theta48",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta49 = list(
                            hyperid = 29149,
                            name = "theta49",
                            short.name = "theta49",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta50 = list(
                            hyperid = 29150,
                            name = "theta50",
                            short.name = "theta50",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta51 = list(
                            hyperid = 29151,
                            name = "theta51",
                            short.name = "theta51",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta52 = list(
                            hyperid = 29152,
                            name = "theta52",
                            short.name = "theta52",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta53 = list(
                            hyperid = 29153,
                            name = "theta53",
                            short.name = "theta53",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta54 = list(
                            hyperid = 29154,
                            name = "theta54",
                            short.name = "theta54",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta55 = list(
                            hyperid = 29155,
                            name = "theta55",
                            short.name = "theta55",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta56 = list(
                            hyperid = 29156,
                            name = "theta56",
                            short.name = "theta56",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta57 = list(
                            hyperid = 29157,
                            name = "theta57",
                            short.name = "theta57",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta58 = list(
                            hyperid = 29158,
                            name = "theta58",
                            short.name = "theta58",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta59 = list(
                            hyperid = 29159,
                            name = "theta59",
                            short.name = "theta59",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta60 = list(
                            hyperid = 29160,
                            name = "theta60",
                            short.name = "theta60",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta61 = list(
                            hyperid = 29161,
                            name = "theta61",
                            short.name = "theta61",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta62 = list(
                            hyperid = 29162,
                            name = "theta62",
                            short.name = "theta62",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta63 = list(
                            hyperid = 29163,
                            name = "theta63",
                            short.name = "theta63",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta64 = list(
                            hyperid = 29164,
                            name = "theta64",
                            short.name = "theta64",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta65 = list(
                            hyperid = 29165,
                            name = "theta65",
                            short.name = "theta65",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta66 = list(
                            hyperid = 29166,
                            name = "theta66",
                            short.name = "theta66",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta67 = list(
                            hyperid = 29167,
                            name = "theta67",
                            short.name = "theta67",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta68 = list(
                            hyperid = 29168,
                            name = "theta68",
                            short.name = "theta68",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta69 = list(
                            hyperid = 29169,
                            name = "theta69",
                            short.name = "theta69",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta70 = list(
                            hyperid = 29170,
                            name = "theta70",
                            short.name = "theta70",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta71 = list(
                            hyperid = 29171,
                            name = "theta71",
                            short.name = "theta71",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta72 = list(
                            hyperid = 29172,
                            name = "theta72",
                            short.name = "theta72",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta73 = list(
                            hyperid = 29173,
                            name = "theta73",
                            short.name = "theta73",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta74 = list(
                            hyperid = 29174,
                            name = "theta74",
                            short.name = "theta74",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta75 = list(
                            hyperid = 29175,
                            name = "theta75",
                            short.name = "theta75",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta76 = list(
                            hyperid = 29176,
                            name = "theta76",
                            short.name = "theta76",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta77 = list(
                            hyperid = 29177,
                            name = "theta77",
                            short.name = "theta77",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta78 = list(
                            hyperid = 29178,
                            name = "theta78",
                            short.name = "theta78",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta79 = list(
                            hyperid = 29179,
                            name = "theta79",
                            short.name = "theta79",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta80 = list(
                            hyperid = 29180,
                            name = "theta80",
                            short.name = "theta80",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta81 = list(
                            hyperid = 29181,
                            name = "theta81",
                            short.name = "theta81",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta82 = list(
                            hyperid = 29182,
                            name = "theta82",
                            short.name = "theta82",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta83 = list(
                            hyperid = 29183,
                            name = "theta83",
                            short.name = "theta83",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta84 = list(
                            hyperid = 29184,
                            name = "theta84",
                            short.name = "theta84",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta85 = list(
                            hyperid = 29185,
                            name = "theta85",
                            short.name = "theta85",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta86 = list(
                            hyperid = 29186,
                            name = "theta86",
                            short.name = "theta86",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta87 = list(
                            hyperid = 29187,
                            name = "theta87",
                            short.name = "theta87",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta88 = list(
                            hyperid = 29188,
                            name = "theta88",
                            short.name = "theta88",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta89 = list(
                            hyperid = 29189,
                            name = "theta89",
                            short.name = "theta89",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta90 = list(
                            hyperid = 29190,
                            name = "theta90",
                            short.name = "theta90",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta91 = list(
                            hyperid = 29191,
                            name = "theta91",
                            short.name = "theta91",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta92 = list(
                            hyperid = 29192,
                            name = "theta92",
                            short.name = "theta92",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta93 = list(
                            hyperid = 29193,
                            name = "theta93",
                            short.name = "theta93",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta94 = list(
                            hyperid = 29194,
                            name = "theta94",
                            short.name = "theta94",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta95 = list(
                            hyperid = 29195,
                            name = "theta95",
                            short.name = "theta95",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta96 = list(
                            hyperid = 29196,
                            name = "theta96",
                            short.name = "theta96",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta97 = list(
                            hyperid = 29197,
                            name = "theta97",
                            short.name = "theta97",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta98 = list(
                            hyperid = 29198,
                            name = "theta98",
                            short.name = "theta98",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta99 = list(
                            hyperid = 29199,
                            name = "theta99",
                            short.name = "theta99",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta100 = list(
                            hyperid = 29200,
                            name = "theta100",
                            short.name = "theta100",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta101 = list(
                            hyperid = 29201,
                            name = "theta101",
                            short.name = "theta101",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta102 = list(
                            hyperid = 29202,
                            name = "theta102",
                            short.name = "theta102",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta103 = list(
                            hyperid = 29203,
                            name = "theta103",
                            short.name = "theta103",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta104 = list(
                            hyperid = 29204,
                            name = "theta104",
                            short.name = "theta104",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta105 = list(
                            hyperid = 29205,
                            name = "theta105",
                            short.name = "theta105",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta106 = list(
                            hyperid = 29206,
                            name = "theta106",
                            short.name = "theta106",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta107 = list(
                            hyperid = 29207,
                            name = "theta107",
                            short.name = "theta107",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta108 = list(
                            hyperid = 29208,
                            name = "theta108",
                            short.name = "theta108",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta109 = list(
                            hyperid = 29209,
                            name = "theta109",
                            short.name = "theta109",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta110 = list(
                            hyperid = 29210,
                            name = "theta110",
                            short.name = "theta110",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta111 = list(
                            hyperid = 29211,
                            name = "theta111",
                            short.name = "theta111",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta112 = list(
                            hyperid = 29212,
                            name = "theta112",
                            short.name = "theta112",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta113 = list(
                            hyperid = 29213,
                            name = "theta113",
                            short.name = "theta113",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta114 = list(
                            hyperid = 29214,
                            name = "theta114",
                            short.name = "theta114",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta115 = list(
                            hyperid = 29215,
                            name = "theta115",
                            short.name = "theta115",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta116 = list(
                            hyperid = 29216,
                            name = "theta116",
                            short.name = "theta116",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta117 = list(
                            hyperid = 29217,
                            name = "theta117",
                            short.name = "theta117",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta118 = list(
                            hyperid = 29218,
                            name = "theta118",
                            short.name = "theta118",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta119 = list(
                            hyperid = 29219,
                            name = "theta119",
                            short.name = "theta119",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta120 = list(
                            hyperid = 29220,
                            name = "theta120",
                            short.name = "theta120",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta121 = list(
                            hyperid = 29221,
                            name = "theta121",
                            short.name = "theta121",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta122 = list(
                            hyperid = 29222,
                            name = "theta122",
                            short.name = "theta122",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta123 = list(
                            hyperid = 29223,
                            name = "theta123",
                            short.name = "theta123",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta124 = list(
                            hyperid = 29224,
                            name = "theta124",
                            short.name = "theta124",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta125 = list(
                            hyperid = 29225,
                            name = "theta125",
                            short.name = "theta125",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta126 = list(
                            hyperid = 29226,
                            name = "theta126",
                            short.name = "theta126",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta127 = list(
                            hyperid = 29227,
                            name = "theta127",
                            short.name = "theta127",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta128 = list(
                            hyperid = 29228,
                            name = "theta128",
                            short.name = "theta128",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta129 = list(
                            hyperid = 29229,
                            name = "theta129",
                            short.name = "theta129",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta130 = list(
                            hyperid = 29230,
                            name = "theta130",
                            short.name = "theta130",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta131 = list(
                            hyperid = 29231,
                            name = "theta131",
                            short.name = "theta131",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta132 = list(
                            hyperid = 29232,
                            name = "theta132",
                            short.name = "theta132",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta133 = list(
                            hyperid = 29233,
                            name = "theta133",
                            short.name = "theta133",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta134 = list(
                            hyperid = 29234,
                            name = "theta134",
                            short.name = "theta134",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta135 = list(
                            hyperid = 29235,
                            name = "theta135",
                            short.name = "theta135",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta136 = list(
                            hyperid = 29236,
                            name = "theta136",
                            short.name = "theta136",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta137 = list(
                            hyperid = 29237,
                            name = "theta137",
                            short.name = "theta137",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta138 = list(
                            hyperid = 29238,
                            name = "theta138",
                            short.name = "theta138",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta139 = list(
                            hyperid = 29239,
                            name = "theta139",
                            short.name = "theta139",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta140 = list(
                            hyperid = 29240,
                            name = "theta140",
                            short.name = "theta140",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta141 = list(
                            hyperid = 29241,
                            name = "theta141",
                            short.name = "theta141",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta142 = list(
                            hyperid = 29242,
                            name = "theta142",
                            short.name = "theta142",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta143 = list(
                            hyperid = 29243,
                            name = "theta143",
                            short.name = "theta143",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta144 = list(
                            hyperid = 29244,
                            name = "theta144",
                            short.name = "theta144",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta145 = list(
                            hyperid = 29245,
                            name = "theta145",
                            short.name = "theta145",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta146 = list(
                            hyperid = 29246,
                            name = "theta146",
                            short.name = "theta146",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta147 = list(
                            hyperid = 29247,
                            name = "theta147",
                            short.name = "theta147",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta148 = list(
                            hyperid = 29248,
                            name = "theta148",
                            short.name = "theta148",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta149 = list(
                            hyperid = 29249,
                            name = "theta149",
                            short.name = "theta149",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta150 = list(
                            hyperid = 29250,
                            name = "theta150",
                            short.name = "theta150",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta151 = list(
                            hyperid = 29251,
                            name = "theta151",
                            short.name = "theta151",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta152 = list(
                            hyperid = 29252,
                            name = "theta152",
                            short.name = "theta152",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta153 = list(
                            hyperid = 29253,
                            name = "theta153",
                            short.name = "theta153",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta154 = list(
                            hyperid = 29254,
                            name = "theta154",
                            short.name = "theta154",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta155 = list(
                            hyperid = 29255,
                            name = "theta155",
                            short.name = "theta155",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta156 = list(
                            hyperid = 29256,
                            name = "theta156",
                            short.name = "theta156",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta157 = list(
                            hyperid = 29257,
                            name = "theta157",
                            short.name = "theta157",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta158 = list(
                            hyperid = 29258,
                            name = "theta158",
                            short.name = "theta158",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta159 = list(
                            hyperid = 29259,
                            name = "theta159",
                            short.name = "theta159",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta160 = list(
                            hyperid = 29260,
                            name = "theta160",
                            short.name = "theta160",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta161 = list(
                            hyperid = 29261,
                            name = "theta161",
                            short.name = "theta161",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta162 = list(
                            hyperid = 29262,
                            name = "theta162",
                            short.name = "theta162",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta163 = list(
                            hyperid = 29263,
                            name = "theta163",
                            short.name = "theta163",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta164 = list(
                            hyperid = 29264,
                            name = "theta164",
                            short.name = "theta164",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta165 = list(
                            hyperid = 29265,
                            name = "theta165",
                            short.name = "theta165",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta166 = list(
                            hyperid = 29266,
                            name = "theta166",
                            short.name = "theta166",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta167 = list(
                            hyperid = 29267,
                            name = "theta167",
                            short.name = "theta167",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta168 = list(
                            hyperid = 29268,
                            name = "theta168",
                            short.name = "theta168",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta169 = list(
                            hyperid = 29269,
                            name = "theta169",
                            short.name = "theta169",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta170 = list(
                            hyperid = 29270,
                            name = "theta170",
                            short.name = "theta170",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta171 = list(
                            hyperid = 29271,
                            name = "theta171",
                            short.name = "theta171",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta172 = list(
                            hyperid = 29272,
                            name = "theta172",
                            short.name = "theta172",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta173 = list(
                            hyperid = 29273,
                            name = "theta173",
                            short.name = "theta173",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta174 = list(
                            hyperid = 29274,
                            name = "theta174",
                            short.name = "theta174",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta175 = list(
                            hyperid = 29275,
                            name = "theta175",
                            short.name = "theta175",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta176 = list(
                            hyperid = 29276,
                            name = "theta176",
                            short.name = "theta176",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta177 = list(
                            hyperid = 29277,
                            name = "theta177",
                            short.name = "theta177",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta178 = list(
                            hyperid = 29278,
                            name = "theta178",
                            short.name = "theta178",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta179 = list(
                            hyperid = 29279,
                            name = "theta179",
                            short.name = "theta179",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta180 = list(
                            hyperid = 29280,
                            name = "theta180",
                            short.name = "theta180",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta181 = list(
                            hyperid = 29281,
                            name = "theta181",
                            short.name = "theta181",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta182 = list(
                            hyperid = 29282,
                            name = "theta182",
                            short.name = "theta182",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta183 = list(
                            hyperid = 29283,
                            name = "theta183",
                            short.name = "theta183",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta184 = list(
                            hyperid = 29284,
                            name = "theta184",
                            short.name = "theta184",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta185 = list(
                            hyperid = 29285,
                            name = "theta185",
                            short.name = "theta185",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta186 = list(
                            hyperid = 29286,
                            name = "theta186",
                            short.name = "theta186",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta187 = list(
                            hyperid = 29287,
                            name = "theta187",
                            short.name = "theta187",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta188 = list(
                            hyperid = 29288,
                            name = "theta188",
                            short.name = "theta188",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta189 = list(
                            hyperid = 29289,
                            name = "theta189",
                            short.name = "theta189",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta190 = list(
                            hyperid = 29290,
                            name = "theta190",
                            short.name = "theta190",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta191 = list(
                            hyperid = 29291,
                            name = "theta191",
                            short.name = "theta191",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta192 = list(
                            hyperid = 29292,
                            name = "theta192",
                            short.name = "theta192",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta193 = list(
                            hyperid = 29293,
                            name = "theta193",
                            short.name = "theta193",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta194 = list(
                            hyperid = 29294,
                            name = "theta194",
                            short.name = "theta194",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta195 = list(
                            hyperid = 29295,
                            name = "theta195",
                            short.name = "theta195",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta196 = list(
                            hyperid = 29296,
                            name = "theta196",
                            short.name = "theta196",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta197 = list(
                            hyperid = 29297,
                            name = "theta197",
                            short.name = "theta197",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta198 = list(
                            hyperid = 29298,
                            name = "theta198",
                            short.name = "theta198",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta199 = list(
                            hyperid = 29299,
                            name = "theta199",
                            short.name = "theta199",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta200 = list(
                            hyperid = 29300,
                            name = "theta200",
                            short.name = "theta200",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta201 = list(
                            hyperid = 29301,
                            name = "theta201",
                            short.name = "theta201",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta202 = list(
                            hyperid = 29302,
                            name = "theta202",
                            short.name = "theta202",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta203 = list(
                            hyperid = 29303,
                            name = "theta203",
                            short.name = "theta203",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta204 = list(
                            hyperid = 29304,
                            name = "theta204",
                            short.name = "theta204",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta205 = list(
                            hyperid = 29305,
                            name = "theta205",
                            short.name = "theta205",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta206 = list(
                            hyperid = 29306,
                            name = "theta206",
                            short.name = "theta206",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta207 = list(
                            hyperid = 29307,
                            name = "theta207",
                            short.name = "theta207",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta208 = list(
                            hyperid = 29308,
                            name = "theta208",
                            short.name = "theta208",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta209 = list(
                            hyperid = 29309,
                            name = "theta209",
                            short.name = "theta209",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta210 = list(
                            hyperid = 29310,
                            name = "theta210",
                            short.name = "theta210",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta211 = list(
                            hyperid = 29311,
                            name = "theta211",
                            short.name = "theta211",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta212 = list(
                            hyperid = 29312,
                            name = "theta212",
                            short.name = "theta212",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta213 = list(
                            hyperid = 29313,
                            name = "theta213",
                            short.name = "theta213",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta214 = list(
                            hyperid = 29314,
                            name = "theta214",
                            short.name = "theta214",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta215 = list(
                            hyperid = 29315,
                            name = "theta215",
                            short.name = "theta215",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta216 = list(
                            hyperid = 29316,
                            name = "theta216",
                            short.name = "theta216",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta217 = list(
                            hyperid = 29317,
                            name = "theta217",
                            short.name = "theta217",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta218 = list(
                            hyperid = 29318,
                            name = "theta218",
                            short.name = "theta218",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta219 = list(
                            hyperid = 29319,
                            name = "theta219",
                            short.name = "theta219",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta220 = list(
                            hyperid = 29320,
                            name = "theta220",
                            short.name = "theta220",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta221 = list(
                            hyperid = 29321,
                            name = "theta221",
                            short.name = "theta221",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta222 = list(
                            hyperid = 29322,
                            name = "theta222",
                            short.name = "theta222",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta223 = list(
                            hyperid = 29323,
                            name = "theta223",
                            short.name = "theta223",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta224 = list(
                            hyperid = 29324,
                            name = "theta224",
                            short.name = "theta224",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta225 = list(
                            hyperid = 29325,
                            name = "theta225",
                            short.name = "theta225",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta226 = list(
                            hyperid = 29326,
                            name = "theta226",
                            short.name = "theta226",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta227 = list(
                            hyperid = 29327,
                            name = "theta227",
                            short.name = "theta227",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta228 = list(
                            hyperid = 29328,
                            name = "theta228",
                            short.name = "theta228",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta229 = list(
                            hyperid = 29329,
                            name = "theta229",
                            short.name = "theta229",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta230 = list(
                            hyperid = 29330,
                            name = "theta230",
                            short.name = "theta230",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta231 = list(
                            hyperid = 29331,
                            name = "theta231",
                            short.name = "theta231",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta232 = list(
                            hyperid = 29332,
                            name = "theta232",
                            short.name = "theta232",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta233 = list(
                            hyperid = 29333,
                            name = "theta233",
                            short.name = "theta233",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta234 = list(
                            hyperid = 29334,
                            name = "theta234",
                            short.name = "theta234",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta235 = list(
                            hyperid = 29335,
                            name = "theta235",
                            short.name = "theta235",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta236 = list(
                            hyperid = 29336,
                            name = "theta236",
                            short.name = "theta236",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta237 = list(
                            hyperid = 29337,
                            name = "theta237",
                            short.name = "theta237",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta238 = list(
                            hyperid = 29338,
                            name = "theta238",
                            short.name = "theta238",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta239 = list(
                            hyperid = 29339,
                            name = "theta239",
                            short.name = "theta239",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta240 = list(
                            hyperid = 29340,
                            name = "theta240",
                            short.name = "theta240",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta241 = list(
                            hyperid = 29341,
                            name = "theta241",
                            short.name = "theta241",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta242 = list(
                            hyperid = 29342,
                            name = "theta242",
                            short.name = "theta242",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta243 = list(
                            hyperid = 29343,
                            name = "theta243",
                            short.name = "theta243",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta244 = list(
                            hyperid = 29344,
                            name = "theta244",
                            short.name = "theta244",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta245 = list(
                            hyperid = 29345,
                            name = "theta245",
                            short.name = "theta245",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta246 = list(
                            hyperid = 29346,
                            name = "theta246",
                            short.name = "theta246",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta247 = list(
                            hyperid = 29347,
                            name = "theta247",
                            short.name = "theta247",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta248 = list(
                            hyperid = 29348,
                            name = "theta248",
                            short.name = "theta248",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta249 = list(
                            hyperid = 29349,
                            name = "theta249",
                            short.name = "theta249",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta250 = list(
                            hyperid = 29350,
                            name = "theta250",
                            short.name = "theta250",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta251 = list(
                            hyperid = 29351,
                            name = "theta251",
                            short.name = "theta251",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta252 = list(
                            hyperid = 29352,
                            name = "theta252",
                            short.name = "theta252",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta253 = list(
                            hyperid = 29353,
                            name = "theta253",
                            short.name = "theta253",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta254 = list(
                            hyperid = 29354,
                            name = "theta254",
                            short.name = "theta254",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta255 = list(
                            hyperid = 29355,
                            name = "theta255",
                            short.name = "theta255",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta256 = list(
                            hyperid = 29356,
                            name = "theta256",
                            short.name = "theta256",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta257 = list(
                            hyperid = 29357,
                            name = "theta257",
                            short.name = "theta257",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta258 = list(
                            hyperid = 29358,
                            name = "theta258",
                            short.name = "theta258",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta259 = list(
                            hyperid = 29359,
                            name = "theta259",
                            short.name = "theta259",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta260 = list(
                            hyperid = 29360,
                            name = "theta260",
                            short.name = "theta260",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta261 = list(
                            hyperid = 29361,
                            name = "theta261",
                            short.name = "theta261",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta262 = list(
                            hyperid = 29362,
                            name = "theta262",
                            short.name = "theta262",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta263 = list(
                            hyperid = 29363,
                            name = "theta263",
                            short.name = "theta263",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta264 = list(
                            hyperid = 29364,
                            name = "theta264",
                            short.name = "theta264",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta265 = list(
                            hyperid = 29365,
                            name = "theta265",
                            short.name = "theta265",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta266 = list(
                            hyperid = 29366,
                            name = "theta266",
                            short.name = "theta266",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta267 = list(
                            hyperid = 29367,
                            name = "theta267",
                            short.name = "theta267",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta268 = list(
                            hyperid = 29368,
                            name = "theta268",
                            short.name = "theta268",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta269 = list(
                            hyperid = 29369,
                            name = "theta269",
                            short.name = "theta269",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta270 = list(
                            hyperid = 29370,
                            name = "theta270",
                            short.name = "theta270",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta271 = list(
                            hyperid = 29371,
                            name = "theta271",
                            short.name = "theta271",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta272 = list(
                            hyperid = 29372,
                            name = "theta272",
                            short.name = "theta272",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta273 = list(
                            hyperid = 29373,
                            name = "theta273",
                            short.name = "theta273",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta274 = list(
                            hyperid = 29374,
                            name = "theta274",
                            short.name = "theta274",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta275 = list(
                            hyperid = 29375,
                            name = "theta275",
                            short.name = "theta275",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta276 = list(
                            hyperid = 29376,
                            name = "theta276",
                            short.name = "theta276",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta277 = list(
                            hyperid = 29377,
                            name = "theta277",
                            short.name = "theta277",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta278 = list(
                            hyperid = 29378,
                            name = "theta278",
                            short.name = "theta278",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta279 = list(
                            hyperid = 29379,
                            name = "theta279",
                            short.name = "theta279",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta280 = list(
                            hyperid = 29380,
                            name = "theta280",
                            short.name = "theta280",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta281 = list(
                            hyperid = 29381,
                            name = "theta281",
                            short.name = "theta281",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta282 = list(
                            hyperid = 29382,
                            name = "theta282",
                            short.name = "theta282",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta283 = list(
                            hyperid = 29383,
                            name = "theta283",
                            short.name = "theta283",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta284 = list(
                            hyperid = 29384,
                            name = "theta284",
                            short.name = "theta284",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta285 = list(
                            hyperid = 29385,
                            name = "theta285",
                            short.name = "theta285",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta286 = list(
                            hyperid = 29386,
                            name = "theta286",
                            short.name = "theta286",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta287 = list(
                            hyperid = 29387,
                            name = "theta287",
                            short.name = "theta287",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta288 = list(
                            hyperid = 29388,
                            name = "theta288",
                            short.name = "theta288",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta289 = list(
                            hyperid = 29389,
                            name = "theta289",
                            short.name = "theta289",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta290 = list(
                            hyperid = 29390,
                            name = "theta290",
                            short.name = "theta290",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta291 = list(
                            hyperid = 29391,
                            name = "theta291",
                            short.name = "theta291",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta292 = list(
                            hyperid = 29392,
                            name = "theta292",
                            short.name = "theta292",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta293 = list(
                            hyperid = 29393,
                            name = "theta293",
                            short.name = "theta293",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta294 = list(
                            hyperid = 29394,
                            name = "theta294",
                            short.name = "theta294",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta295 = list(
                            hyperid = 29395,
                            name = "theta295",
                            short.name = "theta295",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta296 = list(
                            hyperid = 29396,
                            name = "theta296",
                            short.name = "theta296",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta297 = list(
                            hyperid = 29397,
                            name = "theta297",
                            short.name = "theta297",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta298 = list(
                            hyperid = 29398,
                            name = "theta298",
                            short.name = "theta298",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta299 = list(
                            hyperid = 29399,
                            name = "theta299",
                            short.name = "theta299",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta300 = list(
                            hyperid = 29400,
                            name = "theta300",
                            short.name = "theta300",
                            initial = inla.special.number(),
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        )
                    ), 
                    constr = FALSE,
                    nrow.ncol = FALSE,
                    augmented = TRUE,
                    aug.factor = 1L,
                    aug.constr = c(1L, 2L, 3L, 4L, 5L, 6L, 7L, 8L, 9L, 10L, 11L, 12L, 13L, 14L, 15L, 16L, 17L, 18L, 19L, 20L, 21L, 22L, 23L, 24L),
                    n.div.by = -1,
                    n.required = TRUE,
                    set.default.values = TRUE,
                    pdf = "iidkd"
                ),

                "2diid" = list(
                    doc = "(This model is obsolute)",
                    hyper = list(
                        theta1 = list(
                            hyperid = 30001,
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
                            hyperid = 30002,
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
                            hyperid = 30003,
                            name = "correlation",
                            short.name = "cor",
                            initial = 4,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(0, 0.15),
                            to.theta = function(x) log((1 + x) / (1 - x)),
                            from.theta = function(x) 2 * exp(x) / (1 + exp(x)) - 1
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
                            hyperid = 31001,
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
                    pdf = "z"
                ),

                rw2d = list(
                    doc = "Thin-plate spline model",
                    hyper = list(
                        theta = list(
                            hyperid = 32001,
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
                            hyperid = 33001,
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
                            hyperid = 33002,
                            name = "logit phi",
                            short.name = "phi",
                            prior = "pc",
                            param = c(0.5, 0.5),
                            initial = 3,
                            fixed = FALSE,
                            to.theta = function(x) log(x / (1 - x)),
                            from.theta = function(x) exp(x) / (1 + exp(x))
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
                    pdf = "rw2diid"
                ),

                slm = list(
                    doc = "Spatial lag model",
                    hyper = list(
                        theta1 = list(
                            hyperid = 34001,
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
                            hyperid = 34002,
                            name = "rho",
                            short.name = "rho",
                            initial = 0,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(0, 10),
                            to.theta = function(x) log(x / (1 - x)),
                            from.theta = function(x) 1 / (1 + exp(-x))
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
                    pdf = "slm"
                ),

                matern2d = list(
                    doc = "Matern covariance function on a regular grid",
                    hyper = list(
                        theta1 = list(
                            hyperid = 35001,
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
                            hyperid = 35002,
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
                            hyperid = 35101,
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
                            hyperid = 35102,
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
                            hyperid = 35103,
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
                    pdf = "dmatern"
                ),

                copy = list(
                    doc = "Create a copy of a model component",
                    hyper = list(
                        theta = list(
                            hyperid = 36001,
                            name = "beta",
                            short.name = "b",
                            initial = 0.0, ## adaptive: if (fixed) initial=1.0 else initial=0.1
                            fixed = TRUE,
                            prior = "normal",
                            param = c(1, 10),
                            to.theta = function(x, REPLACE.ME.low, REPLACE.ME.high) {
                                if (all(is.infinite(c(low, high))) || low == high) {
                                    return(x)
                                } else if (all(is.finite(c(low, high)))) {
                                    stopifnot(low < high)
                                    return(log(-(low - x) / (high - x)))
                                } else if (is.finite(low) && is.infinite(high) && high > low) {
                                    return(log(x - low))
                                } else {
                                    stop("Condition not yet implemented")
                                }
                            },
                            from.theta = function(x, REPLACE.ME.low, REPLACE.ME.high) {
                                if (all(is.infinite(c(low, high))) || low == high) {
                                    return(x)
                                } else if (all(is.finite(c(low, high)))) {
                                    stopifnot(low < high)
                                    return(low + exp(x) / (1 + exp(x)) * (high - low))
                                } else if (is.finite(low) && is.infinite(high) && high > low) {
                                    return(low + exp(x))
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
                    pdf = "copy"
                ),

                scopy = list(
                    doc = "Create a scopy of a model component",
                    hyper = list(
                        theta1 = list(
                            hyperid = 36101,
                            name = "mean",
                            short.name = "mean",
                            initial = 1.0, 
                            fixed = FALSE,
                            prior = "normal", 
                            param = c(1, 10),
                            to.theta = function(x) x, 
                            from.theta = function(x) x
                        ),
                        theta2 = list(
                            hyperid = 36102,
                            name = "slope",
                            short.name = "slope",
                            initial = 0,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(0, 10), 
                            to.theta = function(x) x, 
                            from.theta = function(x) x
                        ),
                        theta3 = list(
                            hyperid = 36103,
                            name = "spline.theta1",
                            short.name = "spline", ## yes do not use 'spline1' here
                            initial = 0,
                            fixed = FALSE,
                            prior = "laplace",
                            param = c(0, 10), 
                            to.theta = function(x) x, 
                            from.theta = function(x) x
                        ), 
                        theta4 = list(
                            hyperid = 36104,
                            name = "spline.theta2",
                            short.name = "spline2",
                            initial = 0,
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(0), 
                            to.theta = function(x) x, 
                            from.theta = function(x) x
                        ), 
                        theta5 = list(
                            hyperid = 36105,
                            name = "spline.theta3",
                            short.name = "spline3",
                            initial = 0,
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(0), 
                            to.theta = function(x) x, 
                            from.theta = function(x) x
                        ), 
                        theta6 = list(
                            hyperid = 36106,
                            name = "spline.theta4",
                            short.name = "spline4",
                            initial = 0,
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(0), 
                            to.theta = function(x) x, 
                            from.theta = function(x) x
                        ), 
                        theta7 = list(
                            hyperid = 36107,
                            name = "spline.theta5",
                            short.name = "spline5",
                            initial = 0,
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(0), 
                            to.theta = function(x) x, 
                            from.theta = function(x) x
                        ), 
                        theta8 = list(
                            hyperid = 36108,
                            name = "spline.theta6",
                            short.name = "spline6",
                            initial = 0,
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(0), 
                            to.theta = function(x) x, 
                            from.theta = function(x) x
                        ), 
                        theta9 = list(
                            hyperid = 36109,
                            name = "spline.theta7",
                            short.name = "spline7",
                            initial = 0,
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(0), 
                            to.theta = function(x) x, 
                            from.theta = function(x) x
                        ), 
                        theta10 = list(
                            hyperid = 36110,
                            name = "spline.theta8",
                            short.name = "spline8",
                            initial = 0,
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(0), 
                            to.theta = function(x) x, 
                            from.theta = function(x) x
                        ), 
                        theta11 = list(
                            hyperid = 36111,
                            name = "spline.theta9",
                            short.name = "spline9",
                            initial = 0,
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(0), 
                            to.theta = function(x) x, 
                            from.theta = function(x) x
                        ), 
                        theta12 = list(
                            hyperid = 36112,
                            name = "spline.theta10",
                            short.name = "spline10",
                            initial = 0,
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(0), 
                            to.theta = function(x) x, 
                            from.theta = function(x) x
                        ), 
                        theta13 = list(
                            hyperid = 36113,
                            name = "spline.theta11",
                            short.name = "spline11",
                            initial = 0,
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(0), 
                            to.theta = function(x) x, 
                            from.theta = function(x) x
                        ), 
                        theta14 = list(
                            hyperid = 36114,
                            name = "spline.theta12",
                            short.name = "spline12",
                            initial = 0,
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(0), 
                            to.theta = function(x) x, 
                            from.theta = function(x) x
                        ), 
                        theta15 = list(
                            hyperid = 36115,
                            name = "spline.theta13",
                            short.name = "spline13",
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
                    n.required = FALSE,
                    set.default.values = FALSE,
                    pdf = "scopy"
                ),

                clinear = list(
                    doc = "Constrained linear effect",
                    hyper = list(
                        theta = list(
                            hyperid = 37001,
                            name = "beta",
                            short.name = "b",
                            initial = 1,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(1, 10),
                            to.theta = function(x, REPLACE.ME.low, REPLACE.ME.high) {
                                if (all(is.infinite(c(low, high))) || low == high) {
                                    stopifnot(low < high)
                                    return(x)
                                } else if (all(is.finite(c(low, high)))) {
                                    stopifnot(low < high)
                                    return(log(-(low - x) / (high - x)))
                                } else if (is.finite(low) && is.infinite(high) && high > low) {
                                    return(log(x - low))
                                } else {
                                    stop("Condition not yet implemented")
                                }
                            },
                            from.theta = function(x, REPLACE.ME.low, REPLACE.ME.high) {
                                if (all(is.infinite(c(low, high))) || low == high) {
                                    stopifnot(low < high)
                                    return(x)
                                } else if (all(is.finite(c(low, high)))) {
                                    stopifnot(low < high)
                                    return(low + exp(x) / (1 + exp(x)) * (high - low))
                                } else if (is.finite(low) && is.infinite(high) && high > low) {
                                    return(low + exp(x))
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
                            hyperid = 38001,
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
                            hyperid = 38002,
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
                            hyperid = 38003,
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
                    pdf = "sigm"
                ),

                revsigm = list(
                    doc = "Reverse sigmoidal effect of a covariate",
                    hyper = list(
                        theta1 = list(
                            hyperid = 39001,
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
                            hyperid = 39002,
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
                            hyperid = 39003,
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
                    pdf = "sigm"
                ),

                log1exp = list(
                    doc = "A nonlinear model of a covariate",
                    hyper = list(
                        theta1 = list(
                            hyperid = 39011,
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
                            hyperid = 39012,
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
                            hyperid = 39013,
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
                    pdf = "log1exp"
                ),

                logdist = list(
                    doc = "A nonlinear model of a covariate",
                    hyper = list(
                        theta1 = list(
                            hyperid = 39021,
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
                            hyperid = 39022,
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
                            hyperid = 39023,
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
                    pdf = "logdist"
                )
                ##
            )
    )
}

`inla.models.section.group` <- function() {
    ## prevent a warning with R CMD check
    ngroup <- NULL

    list(
        group =
            list(
                exchangeable = list(
                    doc = "Exchangeable correlations",
                    hyper = list(
                        theta = list(
                            hyperid = 40001,
                            name = "logit correlation",
                            short.name = "rho",
                            initial = 1,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(0, 0.2),
                            to.theta = function(x, REPLACE.ME.ngroup) log((1 + x * (ngroup - 1)) / (1 - x)),
                            from.theta = function(x, REPLACE.ME.ngroup) (exp(x) - 1) / (exp(x) + ngroup - 1)
                        )
                    )
                ),

                exchangeablepos = list(
                    doc = "Exchangeable positive correlations",
                    hyper = list(
                        theta = list(
                            hyperid = 40101,
                            name = "logit correlation",
                            short.name = "rho",
                            initial = 1,
                            fixed = FALSE,
                            prior = "pc.cor0",
                            param = c(0.5, 0.5),
                            to.theta = function(x) log(x / (1 - x)),
                            from.theta = function(x) exp(x) / (1 + exp(x))
                        )
                    )
                ),

                ar1 = list(
                    doc = "AR(1) correlations",
                    hyper = list(
                        theta = list(
                            hyperid = 41001,
                            name = "logit correlation",
                            short.name = "rho",
                            initial = 2,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(0, 0.15),
                            to.theta = function(x) log((1 + x) / (1 - x)),
                            from.theta = function(x) 2 * exp(x) / (1 + exp(x)) - 1
                        )
                    )
                ),

                ar = list(
                    doc = "AR(p) correlations",
                    ## to many parameters here, but ...
                    hyper = list(
                        theta1 = list(
                            hyperid = 42001,
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
                            hyperid = 42002,
                            name = "pacf1",
                            short.name = "pacf1",
                            initial = 2,
                            fixed = FALSE,
                            prior = "pc.cor0",
                            param = c(0.5, 0.5),
                            to.theta = function(x) log((1 + x) / (1 - x)),
                            from.theta = function(x) 2 * exp(x) / (1 + exp(x)) - 1
                        ),
                        theta3 = list(
                            hyperid = 42003,
                            name = "pacf2",
                            short.name = "pacf2",
                            initial = 0,
                            fixed = FALSE,
                            prior = "pc.cor0",
                            param = c(0.5, 0.4),
                            to.theta = function(x) log((1 + x) / (1 - x)),
                            from.theta = function(x) 2 * exp(x) / (1 + exp(x)) - 1
                        ),
                        theta4 = list(
                            hyperid = 42004,
                            name = "pacf3",
                            short.name = "pacf3",
                            initial = 0,
                            fixed = FALSE,
                            prior = "pc.cor0",
                            param = c(0.5, 0.3),
                            to.theta = function(x) log((1 + x) / (1 - x)),
                            from.theta = function(x) 2 * exp(x) / (1 + exp(x)) - 1
                        ),
                        theta5 = list(
                            hyperid = 42005,
                            name = "pacf4",
                            short.name = "pacf4",
                            initial = 0,
                            fixed = FALSE,
                            prior = "pc.cor0",
                            param = c(0.5, 0.2),
                            to.theta = function(x) log((1 + x) / (1 - x)),
                            from.theta = function(x) 2 * exp(x) / (1 + exp(x)) - 1
                        ),
                        theta6 = list(
                            hyperid = 42006,
                            name = "pacf5",
                            short.name = "pacf5",
                            initial = 0,
                            fixed = FALSE,
                            prior = "pc.cor0",
                            param = c(0.5, 0.1),
                            to.theta = function(x) log((1 + x) / (1 - x)),
                            from.theta = function(x) 2 * exp(x) / (1 + exp(x)) - 1
                        ),
                        theta7 = list(
                            hyperid = 42007,
                            name = "pacf6",
                            short.name = "pacf6",
                            initial = 0,
                            fixed = FALSE,
                            prior = "pc.cor0",
                            param = c(0.5, 0.1),
                            to.theta = function(x) log((1 + x) / (1 - x)),
                            from.theta = function(x) 2 * exp(x) / (1 + exp(x)) - 1
                        ),
                        theta8 = list(
                            hyperid = 42008,
                            name = "pacf7",
                            short.name = "pacf7",
                            initial = 0,
                            fixed = FALSE,
                            prior = "pc.cor0",
                            param = c(0.5, 0.1),
                            to.theta = function(x) log((1 + x) / (1 - x)),
                            from.theta = function(x) 2 * exp(x) / (1 + exp(x)) - 1
                        ),
                        theta9 = list(
                            hyperid = 42009,
                            name = "pacf8",
                            short.name = "pacf8",
                            initial = 0,
                            fixed = FALSE,
                            prior = "pc.cor0",
                            param = c(0.5, 0.1),
                            to.theta = function(x) log((1 + x) / (1 - x)),
                            from.theta = function(x) 2 * exp(x) / (1 + exp(x)) - 1
                        ),
                        theta10 = list(
                            hyperid = 42010,
                            name = "pacf9",
                            short.name = "pacf9",
                            initial = 0,
                            fixed = FALSE,
                            prior = "pc.cor0",
                            param = c(0.5, 0.1),
                            to.theta = function(x) log((1 + x) / (1 - x)),
                            from.theta = function(x) 2 * exp(x) / (1 + exp(x)) - 1
                        ),
                        theta11 = list(
                            hyperid = 42011,
                            name = "pacf10",
                            short.name = "pacf10",
                            initial = 0,
                            fixed = FALSE,
                            prior = "pc.cor0",
                            param = c(0.5, 0.1),
                            to.theta = function(x) log((1 + x) / (1 - x)),
                            from.theta = function(x) 2 * exp(x) / (1 + exp(x)) - 1
                        )
                    )
                ),

                rw1 = list(
                    doc = "Random walk of order 1",
                    hyper = list(
                        theta = list(
                            hyperid = 43001,
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
                            hyperid = 44001,
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
                            hyperid = 45001,
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
                            hyperid = 46001,
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

`inla.models.section.scopy` <- function() {
    ## the priors for the overall mean and the precision is given elsewhere. this is just to
    ## give the allowed models. 
    list(
        scopy =
            list(
                rw1 = list(
                    doc = "Random walk of order 1",
                    hyper = list()
                ),

                rw2 = list(
                    doc = "Random walk of order 2",
                    hyper = list()
                )
            )
    )
}

`inla.models.section.mix` <- function() {
    list(
        mix =
            list(
                gaussian = list(
                    doc = "Gaussian mixture",
                    hyper = list(
                        theta = list(
                            hyperid = 47001,
                            name = "log precision",
                            short.name = "prec",
                            output.name = "Precision for the Gaussian observations", 
                            output.name.intern = "Log precision for the Gaussian observations", 
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
                            hyperid = 47101,
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
                            hyperid = 47201,
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

`inla.models.section.link` <- function() {
    list(
        link =
            list(
                default = list(
                    doc = "The default link",
                    hyper = list()
                ),

                cloglog = list(
                    doc = "The complementary log-log link",
                    hyper = list()
                ),

                ccloglog = list(
                    doc = "The complement complementary log-log link",
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

                loga = list(
                    doc = "The loga-link",
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
                            hyperid = 48001,
                            name = "sensitivity",
                            short.name = "sens",
                            prior = "logitbeta",
                            param = c(10, 5),
                            initial = 1,
                            fixed = FALSE,
                            to.theta = function(x) log(x / (1 - x)),
                            from.theta = function(x) exp(x) / (1 + exp(x))
                        ),
                        theta2 = list(
                            hyperid = 48002,
                            name = "specificity",
                            short.name = "spec",
                            prior = "logitbeta",
                            param = c(10, 5),
                            initial = 1,
                            fixed = FALSE,
                            to.theta = function(x) log(x / (1 - x)),
                            from.theta = function(x) exp(x) / (1 + exp(x))
                        )
                    ),
                    status = "disabled",
                    pdf = NA
                ),

                logoffset = list(
                    doc = "Log-link with an offset",
                    ## variant = 0, a+exp(...), a>0
                    ## variant = 1, a-exp(...), a>0
                    hyper = list(
                        theta = list(
                            hyperid = 49001,
                            name = "beta",
                            short.name = "b",
                            prior = "normal",
                            param = c(0, 100),
                            initial = 0,
                            fixed = TRUE,
                            to.theta = function(x) log(x),
                            from.theta = function(x) exp(x)
                        )
                    ),
                    pdf = "logoffset"
                ),

                logitoffset = list(
                    doc = "Logit-link with an offset",
                    hyper = list(
                        theta = list(
                            hyperid = 49011,
                            name = "prob",
                            short.name = "p",
                            prior = "normal",
                            param = c(-1, 100),
                            initial = -1,
                            fixed = FALSE,
                            to.theta = function(x) log(x / (1 - x)),
                            from.theta = function(x) exp(x) / (1 + exp(x))
                        )
                    ),
                    pdf = "logitoffset"
                ),

                robit = list(
                    doc = "Robit link",
                    hyper = list(
                        theta = list(
                            hyperid = 49021,
                            name = "log degrees of freedom",
                            short.name = "dof",
                            initial = log(5),
                            fixed = TRUE,
                            prior = "pc.dof",
                            param = c(50, 0.5),
                            to.theta = function(x) log(x - 2),
                            from.theta = function(x) 2 + exp(x)
                        )
                    ),
                    pdf = "robit"
                ),

                sn = list(
                    doc = "Skew-normal link",
                    hyper = list(
                        theta1 = list(
                            hyperid = 49031,
                            name = "skewness",
                            short.name = "skew",
                            initial = 0.00123456789,
                            fixed = FALSE,
                            prior = "pc.sn",
                            param = 10,
                            ## This value defined by GMRFLib_SN_SKEWMAX
                            to.theta = function(x, skew.max = 0.988) log((1 + x / skew.max) / (1 - x / skew.max)),
                            from.theta = function(x, skew.max = 0.988) skew.max * (2 * exp(x) / (1 + exp(x)) - 1)
                        ),
                        theta2 = list(
                            hyperid = 49032,
                            name = "intercept",
                            short.name = "intercept",
                            initial = 0.0,
                            fixed = FALSE,
                            prior = "linksnintercept",
                            ## (mean, prec) in the corresponding N(mean, prec)
                            ## prior in the probit case
                            param = c(0, 0),
                            to.theta = function(x) log(x / (1 - x)),
                            from.theta = function(x) exp(x) / (1 + exp(x))
                        )
                    ),
                    pdf = "linksn"
                ),

                gev = list(
                    doc = "GEV link",
                    hyper = list(
                        theta1 = list(
                            hyperid = 49033,
                            name = "tail",
                            short.name = "xi",
                            initial = -3,
                            fixed = FALSE,
                            prior = "pc.gevtail",
                            param = c(7, 0.0, 0.5),
                            to.theta = function(x, interval = c(REPLACE.ME.low, REPLACE.ME.high)) log(-(interval[1] - x) / (interval[2] - x)),
                            from.theta = function(x, interval = c(REPLACE.ME.low, REPLACE.ME.high)) interval[1] + (interval[2] - interval[1]) * exp(x) / (1.0 + exp(x))
                        ),
                        theta2 = list(
                            hyperid = 49034,
                            name = "intercept",
                            short.name = "intercept",
                            initial = 0.0,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(0, 1), 
                            to.theta = function(x) log(x / (1 - x)),
                            from.theta = function(x) 1 / (1 + exp(-x))
                        )
                    ),
                    pdf = "gev"
                ),

                cgev = list(
                    doc = "Complement GEV link",
                    hyper = list(
                        theta1 = list(
                            hyperid = 49035,
                            name = "tail",
                            short.name = "xi",
                            initial = -3,
                            fixed = FALSE,
                            prior = "pc.gevtail",
                            param = c(7, 0.0, 0.5),
                            to.theta = function(x, interval = c(REPLACE.ME.low, REPLACE.ME.high)) log(-(interval[1] - x) / (interval[2] - x)),
                            from.theta = function(x, interval = c(REPLACE.ME.low, REPLACE.ME.high)) interval[1] + (interval[2] - interval[1]) * exp(x) / (1.0 + exp(x))
                        ),
                        theta2 = list(
                            hyperid = 49036,
                            name = "intercept",
                            short.name = "intercept",
                            initial = 0.0,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(0, 1), 
                            to.theta = function(x) log(x / (1 - x)),
                            from.theta = function(x) 1 / (1 + exp(-x))
                        )
                    ),
                    pdf = "cgev"
                ),

                powerlogit = list(
                    doc = "Power logit link",
                    hyper = list(
                        theta1 = list(
                            hyperid = 49131,
                            name = "power",
                            short.name = "power",
                            initial = 0.00123456789,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(0, 10),
                            to.theta = function(x) log(x),
                            from.theta = function(x) exp(x)
                        ),
                        theta2 = list(
                            hyperid = 49132,
                            name = "intercept",
                            short.name = "intercept",
                            initial = 0.0,
                            fixed = FALSE,
                            prior = "logitbeta",
                            param = c(1, 1),
                            to.theta = function(x) log(x / (1 - x)),
                            from.theta = function(x) exp(x) / (1 + exp(x))
                        )
                    ),
                    pdf = "powerlogit"
                ),

                test1 = list(
                    doc = "A test1-link function (experimental)",
                    hyper = list(
                        theta = list(
                            hyperid = 50001,
                            name = "beta",
                            short.name = "b",
                            prior = "normal",
                            param = c(0, 100),
                            initial = 0,
                            fixed = FALSE,
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        )
                    ),
                    pdf = NA
                ),

                special1 = list(
                    doc = "A special1-link function (experimental)",
                    hyper = list(
                        theta1 = list(
                            hyperid = 51001,
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
                            hyperid = 51002,
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
                            hyperid = 51003,
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
                            hyperid = 51004,
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
                            hyperid = 51005,
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
                            hyperid = 51006,
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
                            hyperid = 51007,
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
                            hyperid = 51008,
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
                            hyperid = 51009,
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
                            hyperid = 51010,
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
                            hyperid = 51011,
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
                            hyperid = 52001,
                            name = "beta",
                            short.name = "b",
                            prior = "normal",
                            param = c(0, 10),
                            initial = 0,
                            fixed = FALSE,
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        )
                    ),
                    pdf = NA
                )
            )
    )
}

`inla.models.section.predictor` <- function() {
    list(
        predictor =
            list(
                predictor = list(
                    doc = "(do not use)",
                    hyper = list(
                        theta = list(
                            hyperid = 53001,
                            name = "log precision",
                            short.name = "prec",
                            initial = log(1/0.001^2),
                            ## do not change
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

`inla.models.section.hazard` <- function() {
    list(
        hazard =
            list(
                rw1 = list(
                    doc = "A random walk of order 1 for the log-hazard",
                    hyper = list(
                        theta = list(
                            hyperid = 54001,
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
                            hyperid = 55001,
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
                iid = list(
                    doc = "An iid model for the log-hazard",
                    hyper = list(
                        theta = list(
                            hyperid = 55501,
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

`inla.models.section.likelihood` <- function() {
    list(
        likelihood =
            list(
                ## the first non-default link-function is the default one.
                ## the first non-default link-function is the default one.
                fl = list(
                    doc = "The fl likelihood",
                    hyper = list(),
                    survival = FALSE,
                    discrete = TRUE,
                    link = c("default", "identity"), 
                    status = "experimental", 
                    pdf = "fl"
                ),

                poisson = list(
                    doc = "The Poisson likelihood",
                    hyper = list(),
                    survival = FALSE,
                    discrete = TRUE,
                    link = c("default", "log", "logoffset", "quantile", "test1", "special1", "special2"),
                    pdf = "poisson"
                ),

                npoisson = list(
                    doc = "The Normal approximation to the Poisson likelihood",
                    hyper = list(),
                    survival = FALSE,
                    discrete = TRUE,
                    link = c("default", "log", "logoffset"), 
                    pdf = "poisson"
                ),

                nzpoisson = list(
                    doc = "The nzPoisson likelihood",
                    hyper = list(),
                    survival = FALSE,
                    discrete = TRUE,
                    link = c("default", "log", "logoffset"), 
                    pdf = "nzpoisson"
                ),

                xpoisson = list(
                    doc = "The Poisson likelihood (expert version)",
                    hyper = list(),
                    survival = FALSE,
                    discrete = TRUE,
                    link = c("default", "log", "logoffset", "quantile", "test1", "special1", "special2"),
                    pdf = "poisson"
                ),

                cenpoisson = list(
                    doc = "Then censored Poisson likelihood",
                    hyper = list(),
                    survival = FALSE,
                    discrete = TRUE,
                    link = c("default", "log", "logoffset", "test1", "special1", "special2"),
                    pdf = "cenpoisson"
                ),

                cenpoisson2 = list(
                    doc = "Then censored Poisson likelihood (version 2)",
                    hyper = list(),
                    survival = FALSE,
                    discrete = TRUE,
                    link = c("default", "log", "logoffset", "test1", "special1", "special2"),
                    pdf = "cenpoisson2"
                ),

                gpoisson = list(
                    doc = "The generalized Poisson likelihood",
                    hyper = list(
                        theta1 = list(
                            hyperid = 56001,
                            name = "overdispersion",
                            short.name = "phi",
                            output.name = "Overdispersion for gpoisson",
                            output.name.intern = "Log overdispersion for gpoisson", 
                            initial = 0,
                            fixed = FALSE,
                            prior = "loggamma",
                            param = c(1, 1),
                            to.theta = function(x) log(x),
                            from.theta = function(x) exp(x)
                        ),
                        theta2 = list(
                            hyperid = 56002,
                            name = "p",
                            short.name = "p",
                            output.name = "Parameter p for gpoisson", 
                            output.name.intern = "Parameter p_intern for gpoisson", 
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
                    pdf = "gpoisson"
                ),

                poisson.special1 = list(
                    doc = "The Poisson.special1 likelihood",
                    hyper = list(
                        theta = list(
                            hyperid = 56100,
                            name = "logit probability",
                            short.name = "prob",
                            output.name = "one-probability parameter for poisson.special1",
                            output.name.intern = "intern one-probability parameter for poisson.special1",
                            initial = -1,
                            fixed = FALSE,
                            prior = "gaussian",
                            param = c(-1, 0.2),
                            to.theta = function(x) log(x / (1 - x)),
                            from.theta = function(x) exp(x) / (1 + exp(x))
                        )
                    ),
                    survival = FALSE,
                    discrete = TRUE,
                    link = c("default", "log"),
                    pdf = "poisson-special"
                ),

                "0poisson" = list(
                    doc = "New 0-inflated Poisson",
                    hyper = list(
                        theta1 = list(
                            hyperid = 56201,
                            name = "beta1",
                            short.name = "beta1",
                            output.name = "beta1 for 0poisson observations",
                            output.name.intern = "beta1 for 0poisson observations",
                            initial = -4,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(-4, 10),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta2 = list(
                            hyperid = 56202,
                            name = "beta2",
                            short.name = "beta2",
                            output.name = "beta2 for 0poisson observations",
                            output.name.intern = "beta2 for 0poisson observations",
                            initial = 0,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(0, 10),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta3 = list(
                            hyperid = 56203,
                            name = "beta3",
                            short.name = "beta3",
                            output.name = "beta3 for 0poisson observations",
                            output.name.intern = "beta3 for 0poisson observations",
                            initial = 0,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(0, 10),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta4 = list(
                            hyperid = 56204,
                            name = "beta4",
                            short.name = "beta4",
                            output.name = "beta4 for 0poisson observations",
                            output.name.intern = "beta4 for 0poisson observations",
                            initial = 0,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(0, 10),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta5 = list(
                            hyperid = 56205,
                            name = "beta5",
                            short.name = "beta5",
                            output.name = "beta5 for 0poisson observations",
                            output.name.intern = "beta5 for 0poisson observations",
                            initial = 0,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(0, 10),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta6 = list(
                            hyperid = 56206,
                            name = "beta6",
                            short.name = "beta6",
                            output.name = "beta6 for 0poisson observations",
                            output.name.intern = "beta6 for 0poisson observations",
                            initial = 0,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(0, 10),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta7 = list(
                            hyperid = 56207,
                            name = "beta7",
                            short.name = "beta7",
                            output.name = "beta7 for 0poisson observations",
                            output.name.intern = "beta7 for 0poisson observations",
                            initial = 0,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(0, 10),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta8 = list(
                            hyperid = 56208,
                            name = "beta8",
                            short.name = "beta8",
                            output.name = "beta8 for 0poisson observations",
                            output.name.intern = "beta8 for 0poisson observations",
                            initial = 0,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(0, 10),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta9 = list(
                            hyperid = 56209,
                            name = "beta9",
                            short.name = "beta9",
                            output.name = "beta9 for 0poisson observations",
                            output.name.intern = "beta9 for 0poisson observations",
                            initial = 0,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(0, 10),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta10 = list(
                            hyperid = 56210,
                            name = "beta10",
                            short.name = "beta10",
                            output.name = "beta10 for 0poisson observations",
                            output.name.intern = "beta10 for 0poisson observations",
                            initial = 0,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(0, 10),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        )
                    ),
                    survival = FALSE,
                    discrete = TRUE,
                    link = c("default", "log", "quantile"), 
                    link.simple = c("default", "logit", "cauchit", "probit", "cloglog", "ccloglog"), 
                    pdf = "0inflated"
                ), 

                "0poissonS" = list(
                    doc = "New 0-inflated Poisson Swap",
                    hyper = list(
                        theta1 = list(
                            hyperid = 56301,
                            name = "beta1",
                            short.name = "beta1",
                            output.name = "beta1 for 0poissonS observations",
                            output.name.intern = "beta1 for 0poissonS observations",
                            initial = -4,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(-4, 10),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta2 = list(
                            hyperid = 56302,
                            name = "beta2",
                            short.name = "beta2",
                            output.name = "beta2 for 0poissonS observations",
                            output.name.intern = "beta2 for 0poissonS observations",
                            initial = 0,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(0, 10),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta3 = list(
                            hyperid = 56303,
                            name = "beta3",
                            short.name = "beta3",
                            output.name = "beta3 for 0poissonS observations",
                            output.name.intern = "beta3 for 0poissonS observations",
                            initial = 0,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(0, 10),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta4 = list(
                            hyperid = 56304,
                            name = "beta4",
                            short.name = "beta4",
                            output.name = "beta4 for 0poissonS observations",
                            output.name.intern = "beta4 for 0poissonS observations",
                            initial = 0,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(0, 10),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta5 = list(
                            hyperid = 56305,
                            name = "beta5",
                            short.name = "beta5",
                            output.name = "beta5 for 0poissonS observations",
                            output.name.intern = "beta5 for 0poissonS observations",
                            initial = 0,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(0, 10),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta6 = list(
                            hyperid = 56306,
                            name = "beta6",
                            short.name = "beta6",
                            output.name = "beta6 for 0poissonS observations",
                            output.name.intern = "beta6 for 0poissonS observations",
                            initial = 0,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(0, 10),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta7 = list(
                            hyperid = 56307,
                            name = "beta7",
                            short.name = "beta7",
                            output.name = "beta7 for 0poissonS observations",
                            output.name.intern = "beta7 for 0poissonS observations",
                            initial = 0,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(0, 10),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta8 = list(
                            hyperid = 56308,
                            name = "beta8",
                            short.name = "beta8",
                            output.name = "beta8 for 0poissonS observations",
                            output.name.intern = "beta8 for 0poissonS observations",
                            initial = 0,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(0, 10),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta9 = list(
                            hyperid = 56309,
                            name = "beta9",
                            short.name = "beta9",
                            output.name = "beta9 for 0poissonS observations",
                            output.name.intern = "beta9 for 0poissonS observations",
                            initial = 0,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(0, 10),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta10 = list(
                            hyperid = 56310,
                            name = "beta10",
                            short.name = "beta10",
                            output.name = "beta10 for 0poissonS observations",
                            output.name.intern = "beta10 for 0poissonS observations",
                            initial = 0,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(0, 10),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        )
                    ),
                    survival = FALSE,
                    discrete = TRUE,
                    link = c("default", "logit", "loga", "cauchit", "probit", "cloglog", "ccloglog", "loglog",
                             "log", "sslogit", "logitoffset", "quantile", "pquantile", "robit", "sn",
                             "powerlogit"),
                    link.simple = c("default", "log"), 
                    pdf = "0inflated"
                ), 

                bell = list(
                    doc = "The Bell likelihood",
                    hyper = list(),
                    survival = FALSE,
                    discrete = TRUE,
                    link = c("default", "log"), 
                    pdf = "bell"
                ),

                "0binomial" = list(
                    doc = "New 0-inflated Binomial",
                    hyper = list(
                        theta1 = list(
                            hyperid = 56401,
                            name = "beta1",
                            short.name = "beta1",
                            output.name = "beta1 for 0binomial observations", 
                            output.name.intern = "beta1 for 0binomial observations", 
                            initial = -4,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(-4, 10),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta2 = list(
                            hyperid = 56402,
                            name = "beta2",
                            short.name = "beta2",
                            output.name = "beta2 for 0binomial observations", 
                            output.name.intern = "beta2 for 0binomial observations", 
                            initial = 0,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(0, 10),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta3 = list(
                            hyperid = 56403,
                            name = "beta3",
                            short.name = "beta3",
                            output.name = "beta3 for 0binomial observations", 
                            output.name.intern = "beta3 for 0binomial observations", 
                            initial = 0,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(0, 10),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta4 = list(
                            hyperid = 56404,
                            name = "beta4",
                            short.name = "beta4",
                            output.name = "beta4 for 0binomial observations", 
                            output.name.intern = "beta4 for 0binomial observations", 
                            initial = 0,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(0, 10),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta5 = list(
                            hyperid = 56405,
                            name = "beta5",
                            short.name = "beta5",
                            output.name = "beta5 for 0binomial observations", 
                            output.name.intern = "beta5 for 0binomial observations", 
                            initial = 0,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(0, 10),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta6 = list(
                            hyperid = 56406,
                            name = "beta6",
                            short.name = "beta6",
                            output.name = "beta6 for 0binomial observations", 
                            output.name.intern = "beta6 for 0binomial observations", 
                            initial = 0,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(0, 10),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta7 = list(
                            hyperid = 56407,
                            name = "beta7",
                            short.name = "beta7",
                            output.name = "beta7 for 0binomial observations", 
                            output.name.intern = "beta7 for 0binomial observations", 
                            initial = 0,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(0, 10),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta8 = list(
                            hyperid = 56408,
                            name = "beta8",
                            short.name = "beta8",
                            output.name = "beta8 for 0binomial observations", 
                            output.name.intern = "beta8 for 0binomial observations", 
                            initial = 0,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(0, 10),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta9 = list(
                            hyperid = 56409,
                            name = "beta9",
                            short.name = "beta9",
                            output.name = "beta9 for 0binomial observations", 
                            output.name.intern = "beta9 for 0binomial observations", 
                            initial = 0,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(0, 10),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta10 = list(
                            hyperid = 56410,
                            name = "beta10",
                            short.name = "beta10",
                            output.name = "beta10 for 0binomial observations", 
                            output.name.intern = "beta10 for 0binomial observations", 
                            initial = 0,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(0, 10),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        )
                    ),
                    survival = FALSE,
                    discrete = TRUE,
                    link = c("default", "logit", "loga", "cauchit", "probit", "cloglog", "ccloglog", "loglog", "log"), 
                    link.simple = c("default", "logit", "cauchit", "probit", "cloglog", "ccloglog"), 
                    pdf = "0inflated"
                ), 

                "0binomialS" = list(
                    doc = "New 0-inflated Binomial Swap",
                    hyper = list(
                        theta1 = list(
                            hyperid = 56501,
                            name = "beta1",
                            short.name = "beta1",
                            output.name = "beta1 for 0binomialS observations", 
                            output.name.intern = "beta1 for 0binomialS observations", 
                            initial = -4,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(-4, 10),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta2 = list(
                            hyperid = 56502,
                            name = "beta2",
                            short.name = "beta2",
                            output.name = "beta2 for 0binomialS observations", 
                            output.name.intern = "beta2 for 0binomialS observations", 
                            initial = 0,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(0, 10),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta3 = list(
                            hyperid = 56503,
                            name = "beta3",
                            short.name = "beta3",
                            output.name = "beta3 for 0binomialS observations", 
                            output.name.intern = "beta3 for 0binomialS observations", 
                            initial = 0,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(0, 10),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta4 = list(
                            hyperid = 56504,
                            name = "beta4",
                            short.name = "beta4",
                            output.name = "beta4 for 0binomialS observations", 
                            output.name.intern = "beta4 for 0binomialS observations", 
                            initial = 0,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(0, 10),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta5 = list(
                            hyperid = 56505,
                            name = "beta5",
                            short.name = "beta5",
                            output.name = "beta5 for 0binomialS observations", 
                            output.name.intern = "beta5 for 0binomialS observations", 
                            initial = 0,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(0, 10),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta6 = list(
                            hyperid = 56506,
                            name = "beta6",
                            short.name = "beta6",
                            output.name = "beta6 for 0binomialS observations", 
                            output.name.intern = "beta6 for 0binomialS observations", 
                            initial = 0,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(0, 10),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta7 = list(
                            hyperid = 56507,
                            name = "beta7",
                            short.name = "beta7",
                            output.name = "beta7 for 0binomialS observations", 
                            output.name.intern = "beta7 for 0binomialS observations", 
                            initial = 0,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(0, 10),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta8 = list(
                            hyperid = 56508,
                            name = "beta8",
                            short.name = "beta8",
                            output.name = "beta8 for 0binomialS observations", 
                            output.name.intern = "beta8 for 0binomialS observations", 
                            initial = 0,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(0, 10),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta9 = list(
                            hyperid = 56509,
                            name = "beta9",
                            short.name = "beta9",
                            output.name = "beta9 for 0binomialS observations", 
                            output.name.intern = "beta9 for 0binomialS observations", 
                            initial = 0,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(0, 10),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta10 = list(
                            hyperid = 56510,
                            name = "beta10",
                            short.name = "beta10",
                            output.name = "beta10 for 0binomialS observations", 
                            output.name.intern = "beta10 for 0binomialS observations", 
                            initial = 0,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(0, 10),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        )
                    ),
                    survival = FALSE,
                    discrete = TRUE,
                    link = c("default", "logit", "loga", "cauchit", "probit", "cloglog",
                             "ccloglog", "loglog", "log"),
                    link.simple = c("default", "logit", "cauchit", "probit", "cloglog", "ccloglog"), 
                    pdf = "0inflated"
                ), 

                binomial = list(
                    doc = "The Binomial likelihood",
                    hyper = list(),
                    survival = FALSE,
                    discrete = TRUE,
                    link = c("default", "logit", "loga", "cauchit", "probit", "cloglog",
                             "ccloglog", "loglog", "log", "sslogit", "logitoffset", "quantile",
                             "pquantile", "robit", "sn", "powerlogit", "gev", "cgev"),
                    pdf = "binomial"
                ),

                xbinomial = list(
                    doc = "The Binomial likelihood (expert version)",
                    hyper = list(),
                    survival = FALSE,
                    discrete = TRUE,
                    link = c("default", "logit", "loga", "cauchit", "probit", "cloglog", "ccloglog", "loglog",
                             "log", "sslogit", "logitoffset", "quantile", "pquantile", "robit", "sn",
                             "powerlogit", "gev", "cgev"),
                    pdf = "binomial"
                ),

                pom = list(
                    doc = "Likelihood for the proportional odds model",
                    hyper = list(
                        theta1 = list(
                            hyperid = 57101,
                            name = "theta1",
                            short.name = "theta1",
                            output.name = "theta1 for POM", 
                            output.name.intern = "theta1 for POM", 
                            initial = NA,
                            fixed = FALSE,
                            prior = "dirichlet",
                            param = 3.0,
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta2 = list(
                            hyperid = 57102,
                            name = "theta2",
                            short.name = "theta2",
                            output.name = "theta2 for POM", 
                            output.name.intern = "theta2 for POM", 
                            initial = NA,
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(0),
                            to.theta = function(x) log(x),
                            from.theta = function(x) exp(x)
                        ),
                        theta3 = list(
                            hyperid = 57103,
                            name = "theta3",
                            short.name = "theta3",
                            output.name = "theta3 for POM", 
                            output.name.intern = "theta3 for POM", 
                            initial = NA,
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(0),
                            to.theta = function(x) log(x),
                            from.theta = function(x) exp(x)
                        ),
                        theta4 = list(
                            hyperid = 57104,
                            name = "theta4",
                            short.name = "theta4",
                            output.name = "theta4 for POM", 
                            output.name.intern = "theta4 for POM", 
                            initial = NA,
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(0),
                            to.theta = function(x) log(x),
                            from.theta = function(x) exp(x)
                        ),
                        theta5 = list(
                            hyperid = 57105,
                            name = "theta5",
                            short.name = "theta5",
                            output.name = "theta5 for POM", 
                            output.name.intern = "theta5 for POM", 
                            initial = NA,
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(0),
                            to.theta = function(x) log(x),
                            from.theta = function(x) exp(x)
                        ),
                        theta6 = list(
                            hyperid = 57106,
                            name = "theta6",
                            short.name = "theta6",
                            output.name = "theta6 for POM", 
                            output.name.intern = "theta6 for POM", 
                            initial = NA,
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(0),
                            to.theta = function(x) log(x),
                            from.theta = function(x) exp(x)
                        ),
                        theta7 = list(
                            hyperid = 57107,
                            name = "theta7",
                            short.name = "theta7",
                            output.name = "theta7 for POM", 
                            output.name.intern = "theta7 for POM", 
                            initial = NA,
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(0),
                            to.theta = function(x) log(x),
                            from.theta = function(x) exp(x)
                        ),
                        theta8 = list(
                            hyperid = 57108,
                            name = "theta8",
                            short.name = "theta8",
                            output.name = "theta8 for POM", 
                            output.name.intern = "theta8 for POM", 
                            initial = NA,
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(0),
                            to.theta = function(x) log(x),
                            from.theta = function(x) exp(x)
                        ),
                        theta9 = list(
                            hyperid = 57109,
                            name = "theta9",
                            short.name = "theta9",
                            output.name = "theta9 for POM", 
                            output.name.intern = "theta9 for POM", 
                            initial = NA,
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(0),
                            to.theta = function(x) log(x),
                            from.theta = function(x) exp(x)
                        ),
                        theta10 = list(
                            hyperid = 57110,
                            name = "theta10",
                            short.name = "theta10",
                            output.name = "theta10 for POM", 
                            output.name.intern = "theta10 for POM", 
                            initial = NA,
                            fixed = FALSE,
                            prior = "none",
                            param = numeric(0),
                            to.theta = function(x) log(x),
                            from.theta = function(x) exp(x)
                        )
                    ),
                    survival = FALSE,
                    discrete = TRUE,
                    link = c("default", "identity"),
                    pdf = "pom"
                ),

                bgev = list(
                    doc = "The blended Generalized Extreme Value likelihood",
                    hyper = list(
                        theta1 = list(
                            hyperid = 57201,
                            name = "spread",
                            short.name = "sd",
                            output.name = "spread for BGEV observations",
                            output.name.intern = "log spread for BGEV observations",
                            initial = 0,
                            fixed = FALSE,
                            prior = "loggamma",
                            param = c(1, 3),
                            to.theta = function(x) log(x),
                            from.theta = function(x) exp(x)
                        ),
                        theta2 = list(
                            hyperid = 57202,
                            name = "tail",
                            short.name = "xi",
                            output.name = "tail for BGEV observations",
                            output.name.intern = "intern tail for BGEV observations",
                            initial = -4,
                            fixed = FALSE,
                            prior = "pc.gevtail",
                            param = c(7, 0.0, 0.5),
                            to.theta = function(x, interval = c(REPLACE.ME.low, REPLACE.ME.high)) log(-(interval[1] - x) / (interval[2] - x)),
                            from.theta = function(x, interval = c(REPLACE.ME.low, REPLACE.ME.high)) interval[1] + (interval[2] - interval[1]) * exp(x) / (1.0 + exp(x))
                        ),
                        theta3 = list(
                            hyperid = 57203,
                            name = "beta1",
                            short.name = "beta1",
                            output.name = "MUST BE FIXED",
                            output.name.intern = "MUST BE FIXED",
                            initial = NA,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(0, 300),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta4 = list(
                            hyperid = 57204,
                            name = "beta2",
                            short.name = "beta2",
                            output.name = "MUST BE FIXED",
                            output.name.intern = "MUST BE FIXED",
                            initial = NA,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(0, 300),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta5 = list(
                            hyperid = 57205,
                            name = "beta3",
                            short.name = "beta3",
                            output.name = "MUST BE FIXED",
                            output.name.intern = "MUST BE FIXED",
                            initial = NA,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(0, 300),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta6 = list(
                            hyperid = 57206,
                            name = "beta4",
                            short.name = "beta4",
                            output.name = "MUST BE FIXED",
                            output.name.intern = "MUST BE FIXED",
                            initial = NA,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(0, 300),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta7 = list(
                            hyperid = 57207,
                            name = "beta5",
                            short.name = "beta5",
                            output.name = "MUST BE FIXED",
                            output.name.intern = "MUST BE FIXED",
                            initial = NA,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(0, 300),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta8 = list(
                            hyperid = 57208,
                            name = "beta6",
                            short.name = "beta6",
                            output.name = "MUST BE FIXED",
                            output.name.intern = "MUST BE FIXED",
                            initial = NA,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(0, 300),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta9 = list(
                            hyperid = 57209,
                            name = "beta7",
                            short.name = "beta7",
                            output.name = "MUST BE FIXED",
                            output.name.intern = "MUST BE FIXED",
                            initial = NA,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(0, 300),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta10 = list(
                            hyperid = 57210,
                            name = "beta8",
                            short.name = "beta8",
                            output.name = "MUST BE FIXED",
                            output.name.intern = "MUST BE FIXED",
                            initial = NA,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(0, 300),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta11 = list(
                            hyperid = 57211,
                            name = "beta9",
                            short.name = "beta9",
                            output.name = "MUST BE FIXED",
                            output.name.intern = "MUST BE FIXED",
                            initial = NA,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(0, 300),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta12 = list(
                            hyperid = 57212,
                            name = "beta10",
                            short.name = "beta",
                            output.name = "MUST BE FIXED",
                            output.name.intern = "MUST BE FIXED",
                            initial = NA,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(0, 300),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        )
                    ),
                    survival = FALSE,
                    discrete = FALSE,
                    link = c("default", "identity", "log"),
                    pdf = "bgev"
                ),

                gamma = list(
                    doc = "The Gamma likelihood",
                    hyper = list(
                        theta = list(
                            hyperid = 58001,
                            name = "precision parameter",
                            short.name = "prec",
                            output.name = "Precision-parameter for the Gamma observations",
                            output.name.intern = "Intern precision-parameter for the Gamma observations",
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

                mgamma = list(
                    doc = "The modal Gamma likelihood",
                    hyper = list(
                        theta = list(
                            hyperid = 58002,
                            name = "precision parameter",
                            short.name = "prec",
                            output.name = "Precision-parameter for the modal Gamma observations",
                            output.name.intern = "Intern precision-parameter for the modal Gamma observations",
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
                    pdf = "mgamma"
                ),

                gammasurv = list(
                    doc = "The Gamma likelihood (survival)",
                    hyper = list(
                        theta1 = list(
                            hyperid = 58101,
                            name = "precision parameter",
                            short.name = "prec",
                            output.name = "Precision-parameter for the Gamma surv observations",
                            output.name.intern = "Intern precision-parameter for the Gamma surv observations",
                            initial = log(1),
                            fixed = FALSE,
                            prior = "loggamma",
                            param = c(1, 0.01),
                            to.theta = function(x) log(x),
                            from.theta = function(x) exp(x)
                        ),
                        theta2 = list(
                            hyperid = 58102,
                            name = "beta1",
                            short.name = "beta1",
                            output.name = "beta1 for Gamma-Cure",
                            output.name.intern = "beta1 for Gamma-Cure",
                            initial = -7,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(-4, 100),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta3 = list(
                            hyperid = 58103,
                            name = "beta2",
                            short.name = "beta2",
                            output.name = "beta2 for Gamma-Cure",
                            output.name.intern = "beta2 for Gamma-Cure",
                            initial = 0,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(0, 100),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta4 = list(
                            hyperid = 58104,
                            name = "beta3",
                            short.name = "beta3",
                            output.name = "beta3 for Gamma-Cure",
                            output.name.intern = "beta3 for Gamma-Cure",
                            initial = 0,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(0, 100),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta5 = list(
                            hyperid = 58105,
                            name = "beta4",
                            short.name = "beta4",
                            output.name = "beta4 for Ga mma-Cure",
                            output.name.intern = "beta4 for Gamma-Cure",
                            initial = 0,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(0, 100),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta6 = list(
                            hyperid = 58106,
                            name = "beta5",
                            short.name = "beta5",
                            output.name = "beta5 for Gamma-Cure",
                            output.name.intern = "beta5 for Gamma-Cure",
                            initial = 0,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(0, 100),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta7 = list(
                            hyperid = 58107,
                            name = "beta6",
                            short.name = "beta6",
                            output.name = "beta6 for Gamma-Cure",
                            output.name.intern = "beta6 for Gamma-Cure",
                            initial = 0,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(0, 100),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta8 = list(
                            hyperid = 58108,
                            name = "beta7",
                            short.name = "beta7",
                            output.name = "beta7 for Gamma-Cure",
                            output.name.intern = "beta7 for Gamma-Cure",
                            initial = 0,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(0, 100),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta9 = list(
                            hyperid = 58109,
                            name = "beta8",
                            short.name = "beta8",
                            output.name = "beta8 for Gamma-Cure",
                            output.name.intern = "beta8 for Gamma-Cure",
                            initial = 0,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(0, 100),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta10 = list(
                            hyperid = 58110,
                            name = "beta9",
                            short.name = "beta9",
                            output.name = "beta9 for Gamma-Cure",
                            output.name.intern = "beta9 for Gamma-Cure",
                            initial = 0,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(0, 100),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta11 = list(
                            hyperid = 58111,
                            name = "beta10",
                            short.name = "beta10",
                            output.name = "beta10 for Gamma-Cure",
                            output.name.intern = "beta10 for Gamma-Cure",
                            initial = 0,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(0, 100),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        )
                    ), 
                    survival = TRUE,
                    discrete = FALSE,
                    link = c("default", "log", "neglog", "quantile"),
                    pdf = "gammasurv"
                ),
                
                mgammasurv = list(
                    doc = "The modal Gamma likelihood (survival)",
                    hyper = list(
                        theta1 = list(
                            hyperid = 58121,
                            name = "precision parameter",
                            short.name = "prec",
                            output.name = "Precision-parameter for the modal Gamma surv observations",
                            output.name.intern = "Intern precision-parameter for the modal Gamma surv observations",
                            initial = log(1),
                            fixed = FALSE,
                            prior = "loggamma",
                            param = c(1, 0.01),
                            to.theta = function(x) log(x),
                            from.theta = function(x) exp(x)
                        ),
                        theta2 = list(
                            hyperid = 58122,
                            name = "beta1",
                            short.name = "beta1",
                            output.name = "beta1 for modal Gamma-Cure",
                            output.name.intern = "beta1 for modal Gamma-Cure",
                            initial = -7,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(-4, 100),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta3 = list(
                            hyperid = 58123,
                            name = "beta2",
                            short.name = "beta2",
                            output.name = "beta2 for modal Gamma-Cure",
                            output.name.intern = "beta2 for modal Gamma-Cure",
                            initial = 0,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(0, 100),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta4 = list(
                            hyperid = 58124,
                            name = "beta3",
                            short.name = "beta3",
                            output.name = "beta3 for modal Gamma-Cure",
                            output.name.intern = "beta3 for modal Gamma-Cure",
                            initial = 0,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(0, 100),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta5 = list(
                            hyperid = 58125,
                            name = "beta4",
                            short.name = "beta4",
                            output.name = "beta4 for Ga mma-Cure",
                            output.name.intern = "beta4 for modal Gamma-Cure",
                            initial = 0,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(0, 100),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta6 = list(
                            hyperid = 58126,
                            name = "beta5",
                            short.name = "beta5",
                            output.name = "beta5 for modal Gamma-Cure",
                            output.name.intern = "beta5 for modal Gamma-Cure",
                            initial = 0,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(0, 100),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta7 = list(
                            hyperid = 58127,
                            name = "beta6",
                            short.name = "beta6",
                            output.name = "beta6 for modal Gamma-Cure",
                            output.name.intern = "beta6 for modal Gamma-Cure",
                            initial = 0,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(0, 100),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta8 = list(
                            hyperid = 58128,
                            name = "beta7",
                            short.name = "beta7",
                            output.name = "beta7 for modal Gamma-Cure",
                            output.name.intern = "beta7 for modal Gamma-Cure",
                            initial = 0,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(0, 100),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta9 = list(
                            hyperid = 58129,
                            name = "beta8",
                            short.name = "beta8",
                            output.name = "beta8 for modal Gamma-Cure",
                            output.name.intern = "beta8 for modal Gamma-Cure",
                            initial = 0,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(0, 100),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta10 = list(
                            hyperid = 58130,
                            name = "beta9",
                            short.name = "beta9",
                            output.name = "beta9 for modal Gamma-Cure",
                            output.name.intern = "beta9 for modal Gamma-Cure",
                            initial = 0,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(0, 100),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta11 = list(
                            hyperid = 58131,
                            name = "beta10",
                            short.name = "beta10",
                            output.name = "beta10 for modal Gamma-Cure",
                            output.name.intern = "beta10 for modal Gamma-Cure",
                            initial = 0,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(0, 100),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        )
                    ), 
                    survival = TRUE,
                    discrete = FALSE,
                    link = c("default", "log", "neglog"),
                    pdf = "agamma"
                ),
                
                gammajw = list(
                    doc = "A special case of the Gamma likelihood",
                    hyper =  list(), 
                    survival = FALSE,
                    discrete = FALSE,
                    link = c("default", "log", "neglog"),
                    pdf = "gammajw"
                ),

                gammajwsurv = list(
                    doc = "A special case of the Gamma likelihood (survival)",
                    hyper = list(
                        theta1 = list(
                            hyperid = 58200,
                            name = "beta1",
                            short.name = "beta1",
                            output.name = "beta1 for GammaJW-Cure", 
                            output.name.intern = "beta1 for GammaJW-Cure", 
                            initial = -7,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(-4, 100),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta2 = list(
                            hyperid = 58201,
                            name = "beta2",
                            short.name = "beta2",
                            output.name = "beta1 for GammaJW-Cure", 
                            output.name.intern = "beta1 for GammaJW-Cure", 
                            initial = 0,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(0, 100),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                       theta3 = list(
                            hyperid = 58202,
                            name = "beta3",
                            short.name = "beta3",
                            output.name = "beta3 for GammaJW-Cure", 
                            output.name.intern = "beta3 for GammaJW-Cure", 
                            initial = 0,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(0, 100),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta4 = list(
                            hyperid = 58203,
                            name = "beta4",
                            short.name = "beta4",
                            output.name = "beta4 for GammaJW-Cure", 
                            output.name.intern = "beta4 for GammaJW-Cure", 
                            initial = 0,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(0, 100),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta5 = list(
                            hyperid = 58204,
                            name = "beta5",
                            short.name = "beta5",
                            output.name = "beta5 for GammaJW-Cure", 
                            output.name.intern = "beta5 for GammaJW-Cure", 
                            initial = 0,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(0, 100),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta6 = list(
                            hyperid = 58205,
                            name = "beta6",
                            short.name = "beta6",
                            output.name = "beta6 for GammaJW-Cure", 
                            output.name.intern = "beta6 for GammaJW-Cure", 
                            initial = 0,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(0, 100), 
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta7 = list(
                            hyperid = 58206,
                            name = "beta7",
                            short.name = "beta7",
                            output.name = "beta7 for GammaJW-Cure", 
                            output.name.intern = "beta7 for GammaJW-Cure", 
                            initial = 0,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(0, 100),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta8 = list(
                            hyperid = 58207,
                            name = "beta8",
                            short.name = "beta8",
                            output.name = "beta8 for GammaJW-Cure", 
                            output.name.intern = "beta8 for GammaJW-Cure", 
                            initial = 0,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(0, 100),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta9 = list(
                            hyperid = 58208,
                            name = "beta9",
                            short.name = "beta9",
                            output.name = "beta9 for GammaJW-Cure", 
                            output.name.intern = "beta9 for GammaJW-Cure", 
                            initial = 0,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(0, 100),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta10 = list(
                            hyperid = 58209,
                            name = "beta10",
                            short.name = "beta10",
                            output.name = "beta10 for GammaJW-Cure", 
                            output.name.intern = "beta10 for GammaJW-Cure", 
                            initial = 0,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(0, 100),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        )
                    ), 
                    survival = TRUE,
                    discrete = FALSE,
                    link = c("default", "log"),
                    pdf = "gammajw"
                ),

                gammacount = list(
                    doc = "A Gamma generalisation of the Poisson likelihood",
                    hyper = list(
                        theta = list(
                            hyperid = 59001,
                            name = "log alpha",
                            short.name = "alpha",
                            output.name = "Log-alpha parameter for Gammacount observations",
                            output.name.intern = "Alpha parameter for Gammacount observations",
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
                    pdf = "gammacount"
                ),

                qkumar = list(
                    doc = "A quantile version of the Kumar likelihood",
                    hyper = list(
                        theta = list(
                            hyperid = 60001,
                            name = "precision parameter",
                            short.name = "prec",
                            output.name = "precision for qkumar observations",
                            output.name.intern = "log precision for qkumar observations", 
                            initial = 1,
                            fixed = FALSE,
                            prior = "loggamma",
                            param = c(1, 0.1),
                            ## the 'sc' constant is defined in inla.h, and must be the same.
                            ## I know, this is hard-coded for the moment. Should be a generic
                            ## way of doing this...
                            to.theta = function(x, sc = 0.1) log(x) / sc,
                            from.theta = function(x, sc = 0.1) exp(sc * x)
                        )
                    ),
                    survival = FALSE,
                    discrete = FALSE,
                    link = c("default", "logit", "loga", "cauchit"),
                    pdf = "qkumar"
                ),

                qloglogistic = list(
                    doc = "A quantile loglogistic likelihood",
                    hyper = list(
                        theta = list(
                            hyperid = 60011,
                            name = "log alpha",
                            short.name = "alpha",
                            output.name = "alpha for qloglogistic observations", 
                            output.name.intern = "log alpha for qloglogistic observations",
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
                        theta1 = list(
                            hyperid = 60021,
                            name = "log alpha",
                            short.name = "alpha",
                            output.name = "alpha for qloglogisticsurv observations", 
                            output.name.intern = "log alpha for qloglogisticsurv observations",
                            initial = 1,
                            fixed = FALSE,
                            prior = "loggamma",
                            param = c(25, 25),
                            to.theta = function(x) log(x),
                            from.theta = function(x) exp(x)
                        ),
                        theta2 = list(
                            hyperid = 60022,
                            name = "beta1",
                            short.name = "beta1",
                            output.name = "beta1 for qlogLogistic-Cure", 
                            output.name.intern = "beta1 for logLogistic-Cure", 
                            initial = -5,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(-4, 100),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta3 = list(
                            hyperid = 60023,
                            name = "beta2",
                            short.name = "beta2",
                            output.name = "beta2 for qlogLogistic-Cure", 
                            output.name.intern = "beta2 for logLogistic-Cure", 
                            initial = 0,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(0, 100),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta4 = list(
                            hyperid = 60024,
                            name = "beta3",
                            short.name = "beta3",
                            output.name = "beta3 for qlogLogistic-Cure", 
                            output.name.intern = "beta3 for qlogLogistic-Cure", 
                            initial = 0,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(0, 100),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta5 = list(
                            hyperid = 60025,
                            name = "beta4",
                            short.name = "beta4",
                            output.name = "beta4 for qlogLogistic-Cure", 
                            output.name.intern = "beta4 for qlogLogistic-Cure", 
                            initial = 0,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(0, 100),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta6 = list(
                            hyperid = 60026,
                            name = "beta5",
                            short.name = "beta5",
                            output.name = "beta5 for qlogLogistic-Cure", 
                            output.name.intern = "beta5 for qlogLogistic-Cure", 
                            initial = 0,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(0, 100),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta7 = list(
                            hyperid = 60027,
                            name = "beta6",
                            short.name = "beta6",
                            output.name = "beta6 for qlogLogistic-Cure", 
                            output.name.intern = "beta6 for qlogLogistic-Cure", 
                            initial = 0,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(0, 100),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta8 = list(
                            hyperid = 60028,
                            name = "beta7",
                            short.name = "beta7",
                            output.name = "beta7 for qlogLogistic-Cure", 
                            output.name.intern = "beta7 for qlogLogistic-Cure", 
                            initial = 0,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(0, 100),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta9 = list(
                            hyperid = 60029,
                            name = "beta8",
                            short.name = "beta8",
                            output.name = "beta8 for qlogLogistic-Cure", 
                            output.name.intern = "beta8 for qlogLogistic-Cure", 
                            initial = 0,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(0, 100),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta10 = list(
                            hyperid = 60030,
                            name = "beta9",
                            short.name = "beta9",
                            output.name = "beta9 for qlogLogistic-Cure", 
                            output.name.intern = "beta9 for qlogLogistic-Cure", 
                            initial = 0,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(0, 100),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta11 = list(
                            hyperid = 60031,
                            name = "beta10",
                            short.name = "beta10",
                            output.name = "beta10 for qlogLogistic-Cure", 
                            output.name.intern = "beta10 for qlogLogistic-Cure", 
                            initial = 0,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(0, 100),
                            to.theta = function(x) x,
                            from.theta = function(x) x
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
                            hyperid = 61001,
                            name = "precision parameter",
                            short.name = "phi",
                            output.name = "precision parameter for the beta observations",
                            output.name.intern = "intern precision-parameter for the beta observations",
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
                    link = c("default", "logit", "loga", "cauchit", "probit", "cloglog", "ccloglog", "loglog"),
                    pdf = "beta"
                ),

                betabinomial = list(
                    doc = "The Beta-Binomial likelihood",
                    hyper = list(
                        theta = list(
                            hyperid = 62001,
                            name = "overdispersion",
                            short.name = "rho",
                            output.name = "overdispersion for the betabinomial observations",
                            output.name.intern = "intern overdispersion for the betabinomial observations",
                            initial = 0,
                            fixed = FALSE,
                            prior = "gaussian",
                            param = c(0.0, 0.4),
                            to.theta = function(x) log(x / (1 - x)),
                            from.theta = function(x) exp(x) / (1 + exp(x))
                        )
                    ),
                    survival = FALSE,
                    discrete = TRUE,
                    link = c("default", "logit", "loga", "cauchit", "probit", "cloglog", "ccloglog", "loglog", "robit", "sn"),
                    pdf = "betabinomial"
                ),

                betabinomialna = list(
                    doc = "The Beta-Binomial Normal approximation likelihood",
                    hyper = list(
                        theta = list(
                            hyperid = 62101,
                            name = "overdispersion",
                            short.name = "rho",
                            output.name = "overdispersion for the betabinomialna observations",
                            output.name.intern = "intern overdispersion for the betabinomialna observations",
                            initial = 0,
                            fixed = FALSE,
                            prior = "gaussian",
                            param = c(0.0, 0.4),
                            to.theta = function(x) log(x / (1 - x)),
                            from.theta = function(x) exp(x) / (1 + exp(x))
                        )
                    ),
                    survival = FALSE,
                    discrete = TRUE,
                    link = c("default", "logit", "loga", "cauchit", "probit", "cloglog", "ccloglog", "loglog", "robit", "sn"),
                    pdf = "betabinomialna"
                ),

                cbinomial = list(
                    doc = "The clustered Binomial likelihood",
                    hyper = list(),
                    survival = FALSE,
                    discrete = TRUE,
                    link = c("default", "logit", "loga", "cauchit", "probit", "cloglog", "ccloglog", "loglog", "robit", "sn"),
                    pdf = "cbinomial"
                ),

                nbinomial = list(
                    doc = "The negBinomial likelihood",
                    hyper = list(
                        theta = list(
                            hyperid = 63001,
                            name = "size",
                            short.name = "size",
                            output.name = "size for the nbinomial observations (1/overdispersion)",
                            output.name.intern = "log size for the nbinomial observations (1/overdispersion)",
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
                    link = c("default", "logit", "loga", "cauchit", "probit", "cloglog", "ccloglog", "loglog"),
                    pdf = "nbinomial"
                ),

                cennbinomial2 = list(
                    doc = "The CenNegBinomial2 likelihood (similar to cenpoisson2)",
                    hyper = list(
                        theta = list(
                            hyperid = 63101,
                            name = "size",
                            short.name = "size",
                            output.name = "size for the cennbinomial2 observations (1/overdispersion)",
                            output.name.intern = "log size for the cennbinomial2 observations (1/overdispersion)",
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
                    pdf = "cennbinomial2"
                ),

                simplex = list(
                    doc = "The simplex likelihood",
                    hyper = list(
                        theta = list(
                            hyperid = 64001,
                            name = "log precision",
                            short.name = "prec",
                            output.name = "Precision for the Simplex observations", 
                            output.name.intern = "Log precision for the Simplex observations",
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
                    link = c("default", "logit", "loga", "cauchit", "probit", "cloglog", "ccloglog", "loglog"),
                    pdf = "simplex"
                ),

                gaussian = list(
                    doc = "The Gaussian likelihoood",
                    hyper = list(
                        theta1 = list(
                            hyperid = 65001,
                            name = "log precision",
                            short.name = "prec",
                            output.name = "Precision for the Gaussian observations", 
                            output.name.intern = "Log precision for the Gaussian observations", 
                            initial = 4,
                            fixed = FALSE,
                            prior = "loggamma",
                            param = c(1, 0.00005),
                            to.theta = function(x) log(x),
                            from.theta = function(x) exp(x)
                        ),
                        theta2 = list(
                            hyperid = 65002,
                            name = "log precision offset",
                            short.name = "precoffset",
                            output.name = "NOT IN USE", 
                            output.name.intern = "NOT IN USE", 
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
                    link = c("default", "identity", "logit", "loga", "cauchit", "log", "logoffset"),
                    pdf = "gaussian"
                ),

                stdgaussian = list(
                    doc = "The stdGaussian likelihoood",
                    hyper = list(), 
                    survival = FALSE,
                    discrete = FALSE,
                    link = c("default", "identity", "logit", "loga", "cauchit", "log", "logoffset"),
                    pdf = "gaussian"
                ),

                gaussianjw = list(
                    doc = "The GaussianJW likelihoood",
                    hyper = list(
                        theta1 = list(
                            hyperid = 65101,
                            name = "beta1",
                            short.name = "beta1",
                            output.name = "beta1 for GaussianJW observations", 
                            output.name.intern = "beta1 for GaussianJW observations", 
                            initial = 0,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(0, 100),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta2 = list(
                            hyperid = 65102,
                            name = "beta2",
                            short.name = "beta2",
                            output.name = "beta2 for GaussianJW observations", 
                            output.name.intern = "beta2 for GaussianJW observations", 
                            initial = 1,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(1, 100),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ), 
                        theta3 = list(
                            hyperid = 65103,
                            name = "beta3",
                            short.name = "beta3",
                            output.name = "beta3 for GaussianJW observations", 
                            output.name.intern = "beta3 for GaussianJW observations", 
                            initial = -1,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(-1, 100),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        )
                    ),
                    survival = FALSE,
                    discrete = FALSE,
                    link = c("default", "logit", "probit"), 
                    pdf = "gaussianjw"
                ),

                agaussian = list(
                    doc = "The aggregated Gaussian likelihoood",
                    hyper = list(
                        theta = list(
                            hyperid = 66001,
                            name = "log precision",
                            short.name = "prec",
                            output.name = "Precision for the AggGaussian observations",
                            output.name.intern = "Log precision for the AggGaussian observations",
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
                    link = c("default", "identity", "logit", "loga", "cauchit", "log", "logoffset"),
                    pdf = "agaussian"
                ),

                ggaussian = list(
                    doc = "Generalized Gaussian",
                    hyper = list(
                        theta1 = list(
                            hyperid = 66501,
                            name = "beta1",
                            short.name = "beta1",
                            output.name = "beta1 for ggaussian observations",
                            output.name.intern = "beta1 for ggaussian observations",
                            initial = 4,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(9.33, 0.61), ## as for loggamma(1, 0.00005)
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta2 = list(
                            hyperid = 66502,
                            name = "beta2",
                            short.name = "beta2",
                            output.name = "beta2 for ggaussian observations",
                            output.name.intern = "beta2 for ggaussian observations",
                            initial = 0,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(0, 10),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta3 = list(
                            hyperid = 66503,
                            name = "beta3",
                            short.name = "beta3",
                            output.name = "beta3 for ggaussian observations",
                            output.name.intern = "beta3 for ggaussian observations",
                            initial = 0,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(0, 10),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta4 = list(
                            hyperid = 66504,
                            name = "beta4",
                            short.name = "beta4",
                            output.name = "beta4 for ggaussian observations",
                            output.name.intern = "beta4 for ggaussian observations",
                            initial = 0,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(0, 10),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta5 = list(
                            hyperid = 66505,
                            name = "beta5",
                            short.name = "beta5",
                            output.name = "beta5 for ggaussian observations",
                            output.name.intern = "beta5 for ggaussian observations",
                            initial = 0,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(0, 10),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta6 = list(
                            hyperid = 66506,
                            name = "beta6",
                            short.name = "beta6",
                            output.name = "beta6 for ggaussian observations",
                            output.name.intern = "beta6 for ggaussian observations",
                            initial = 0,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(0, 10),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta7 = list(
                            hyperid = 66507,
                            name = "beta7",
                            short.name = "beta7",
                            output.name = "beta7 for ggaussian observations",
                            output.name.intern = "beta7 for ggaussian observations",
                            initial = 0,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(0, 10),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta8 = list(
                            hyperid = 66508,
                            name = "beta8",
                            short.name = "beta8",
                            output.name = "beta8 for ggaussian observations",
                            output.name.intern = "beta8 for ggaussian observations",
                            initial = 0,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(0, 10),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta9 = list(
                            hyperid = 66509,
                            name = "beta9",
                            short.name = "beta9",
                            output.name = "beta9 for ggaussian observations",
                            output.name.intern = "beta9 for ggaussian observations",
                            initial = 0,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(0, 10),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta10 = list(
                            hyperid = 66510,
                            name = "beta10",
                            short.name = "beta10",
                            output.name = "beta10 for ggaussian observations",
                            output.name.intern = "beta10 for ggaussian observations",
                            initial = 0,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(0, 10),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        )
                    ),
                    survival = FALSE,
                    discrete = FALSE,
                    link = c("default", "identity"), 
                    link.simple = c("default", "log"), 
                    pdf = "ggaussian"
                ), 

                ggaussianS = list(
                    doc = "Generalized GaussianS",
                    hyper = list(
                        theta1 = list(
                            hyperid = 66601,
                            name = "beta1",
                            short.name = "beta1",
                            output.name = "beta1 for ggaussianS observations",
                            output.name.intern = "beta1 for ggaussianS observations",
                            initial = 0,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(0, 0.001),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta2 = list(
                            hyperid = 66602,
                            name = "beta2",
                            short.name = "beta2",
                            output.name = "beta2 for ggaussianS observations",
                            output.name.intern = "beta2 for ggaussianS observations",
                            initial = 0,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(0, 0.001),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta3 = list(
                            hyperid = 66603,
                            name = "beta3",
                            short.name = "beta3",
                            output.name = "beta3 for ggaussianS observations",
                            output.name.intern = "beta3 for ggaussianS observations",
                            initial = 0,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(0, 0.001),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta4 = list(
                            hyperid = 66604,
                            name = "beta4",
                            short.name = "beta4",
                            output.name = "beta4 for ggaussianS observations",
                            output.name.intern = "beta4 for ggaussianS observations",
                            initial = 0,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(0, 0.001),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta5 = list(
                            hyperid = 66605,
                            name = "beta5",
                            short.name = "beta5",
                            output.name = "beta5 for ggaussianS observations",
                            output.name.intern = "beta5 for ggaussianS observations",
                            initial = 0,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(0, 0.001),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta6 = list(
                            hyperid = 66606,
                            name = "beta6",
                            short.name = "beta6",
                            output.name = "beta6 for ggaussianS observations",
                            output.name.intern = "beta6 for ggaussianS observations",
                            initial = 0,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(0, 0.001),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta7 = list(
                            hyperid = 66607,
                            name = "beta7",
                            short.name = "beta7",
                            output.name = "beta7 for ggaussianS observations",
                            output.name.intern = "beta7 for ggaussianS observations",
                            initial = 0,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(0, 0.001),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta8 = list(
                            hyperid = 66608,
                            name = "beta8",
                            short.name = "beta8",
                            output.name = "beta8 for ggaussianS observations",
                            output.name.intern = "beta8 for ggaussianS observations",
                            initial = 0,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(0, 0.001),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta9 = list(
                            hyperid = 66609,
                            name = "beta9",
                            short.name = "beta9",
                            output.name = "beta9 for ggaussianS observations",
                            output.name.intern = "beta9 for ggaussianS observations",
                            initial = 0,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(0, 0.001),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta10 = list(
                            hyperid = 66610,
                            name = "beta10",
                            short.name = "beta10",
                            output.name = "beta10 for ggaussianS observations",
                            output.name.intern = "beta10 for ggaussianS observations",
                            initial = 0,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(0, 0.001),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        )
                    ),
                    survival = FALSE,
                    discrete = FALSE,
                    link = c("default", "log"), 
                    link.simple = c("default", "identity"), 
                    pdf = "ggaussian"
                ), 

                bcgaussian = list(
                    doc = "The Box-Cox Gaussian likelihoood",
                    hyper = list(
                        theta1 = list(
                            hyperid = 65010,
                            name = "log precision",
                            short.name = "prec",
                            output.name = "Precision for the Box-Cox Gaussian observations", 
                            output.name.intern = "Log precision for the Box-Cox Gaussian observations", 
                            initial = 4,
                            fixed = FALSE,
                            prior = "loggamma",
                            param = c(1, 0.00005),
                            to.theta = function(x) log(x),
                            from.theta = function(x) exp(x)
                        ),
                        theta2 = list(
                            hyperid = 65011,
                            name = "Box-Cox transformation parameter",
                            short.name = "lambda",
                            output.name = "NOT IN USE", 
                            output.name.intern = "NOT IN USE", 
                            initial = 1, 
                            fixed = FALSE,
                            prior = "gaussian",
                            param = c(1, 8),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        )
                    ),
                    status = "disabled", 
                    survival = FALSE,
                    discrete = FALSE,
                    link = c("default", "identity"), 
                    pdf = "bcgaussian"
                ),

                rcpoisson = list(
                    doc = "Randomly censored Poisson",
                    hyper = list(
                        theta1 = list(
                            hyperid = 66701,
                            name = "beta1",
                            short.name = "beta1",
                            output.name = "beta1 rcpoisson observations",
                            output.name.intern = "beta1 rcpoisson observations",
                            initial = 0,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(0, 100),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta2 = list(
                            hyperid = 66702,
                            name = "beta2",
                            short.name = "beta2",
                            output.name = "beta2 rcpoisson observations",
                            output.name.intern = "beta2 rcpoisson observations",
                            initial = 0,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(0, 100),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta3 = list(
                            hyperid = 66703,
                            name = "beta3",
                            short.name = "beta3",
                            output.name = "beta3 rcpoisson observations",
                            output.name.intern = "beta3 rcpoisson observations",
                            initial = 0,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(0, 100),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta4 = list(
                            hyperid = 66704,
                            name = "beta4",
                            short.name = "beta4",
                            output.name = "beta4 rcpoisson observations",
                            output.name.intern = "beta4 rcpoisson observations",
                            initial = 0,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(0, 100),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta5 = list(
                            hyperid = 66705,
                            name = "beta5",
                            short.name = "beta5",
                            output.name = "beta5 rcpoisson observations",
                            output.name.intern = "beta5 rcpoisson observations",
                            initial = 0,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(0, 100),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta6 = list(
                            hyperid = 66706,
                            name = "beta6",
                            short.name = "beta6",
                            output.name = "beta6 rcpoisson observations",
                            output.name.intern = "beta6 rcpoisson observations",
                            initial = 0,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(0, 100),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta7 = list(
                            hyperid = 66707,
                            name = "beta7",
                            short.name = "beta7",
                            output.name = "beta7 rcpoisson observations",
                            output.name.intern = "beta7 rcpoisson observations",
                            initial = 0,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(0, 100),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta8 = list(
                            hyperid = 66708,
                            name = "beta8",
                            short.name = "beta8",
                            output.name = "beta8 rcpoisson observations",
                            output.name.intern = "beta8 rcpoisson observations",
                            initial = 0,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(0, 100),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta9 = list(
                            hyperid = 66709,
                            name = "beta9",
                            short.name = "beta9",
                            output.name = "beta9 rcpoisson observations",
                            output.name.intern = "beta9 rcpoisson observations",
                            initial = 0,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(0, 100),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta10 = list(
                            hyperid = 66710,
                            name = "beta10",
                            short.name = "beta10",
                            output.name = "beta10 rcpoisson observations",
                            output.name.intern = "beta10 rcpoisson observations",
                            initial = 0,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(0, 100),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        )
                    ),
                    status = "experimental", 
                    survival = FALSE,
                    discrete = TRUE,
                    link = c("default", "log"), 
                    pdf = "rcpoisson"
                ), 

                tpoisson = list(
                    doc = "Thinned Poisson",
                    hyper = list(
                        theta1 = list(
                            hyperid = 66721,
                            name = "beta1",
                            short.name = "beta1",
                            output.name = "beta1 tpoisson observations",
                            output.name.intern = "beta1 tpoisson observations",
                            initial = 0,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(0, 100),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta2 = list(
                            hyperid = 66722,
                            name = "beta2",
                            short.name = "beta2",
                            output.name = "beta2 tpoisson observations",
                            output.name.intern = "beta2 tpoisson observations",
                            initial = 0,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(0, 100),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta3 = list(
                            hyperid = 66723,
                            name = "beta3",
                            short.name = "beta3",
                            output.name = "beta3 tpoisson observations",
                            output.name.intern = "beta3 tpoisson observations",
                            initial = 0,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(0, 100),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta4 = list(
                            hyperid = 66724,
                            name = "beta4",
                            short.name = "beta4",
                            output.name = "beta4 tpoisson observations",
                            output.name.intern = "beta4 tpoisson observations",
                            initial = 0,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(0, 100),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta5 = list(
                            hyperid = 66725,
                            name = "beta5",
                            short.name = "beta5",
                            output.name = "beta5 tpoisson observations",
                            output.name.intern = "beta5 tpoisson observations",
                            initial = 0,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(0, 100),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta6 = list(
                            hyperid = 66726,
                            name = "beta6",
                            short.name = "beta6",
                            output.name = "beta6 tpoisson observations",
                            output.name.intern = "beta6 tpoisson observations",
                            initial = 0,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(0, 100),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta7 = list(
                            hyperid = 66727,
                            name = "beta7",
                            short.name = "beta7",
                            output.name = "beta7 tpoisson observations",
                            output.name.intern = "beta7 tpoisson observations",
                            initial = 0,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(0, 100),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta8 = list(
                            hyperid = 66728,
                            name = "beta8",
                            short.name = "beta8",
                            output.name = "beta8 tpoisson observations",
                            output.name.intern = "beta8 tpoisson observations",
                            initial = 0,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(0, 100),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta9 = list(
                            hyperid = 66729,
                            name = "beta9",
                            short.name = "beta9",
                            output.name = "beta9 tpoisson observations",
                            output.name.intern = "beta9 tpoisson observations",
                            initial = 0,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(0, 100),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta10 = list(
                            hyperid = 66730,
                            name = "beta10",
                            short.name = "beta10",
                            output.name = "beta10 tpoisson observations",
                            output.name.intern = "beta10 tpoisson observations",
                            initial = 0,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(0, 100),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        )
                    ),
                    status = "experimental", 
                    survival = FALSE,
                    discrete = TRUE,
                    link = c("default", "log"), 
                    pdf = "tpoisson"
                ), 

                circularnormal = list(
                    doc = "The circular Gaussian likelihoood",
                    hyper = list(
                        theta = list(
                            hyperid = 67001,
                            name = "log precision parameter",
                            short.name = "prec",
                            output.name = "Precision parameter for the Circular Normal observations",
                            output.name.intern = "Log precision parameter for the Circular Normal observations",
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
                    status = "disabled"
                ),

                wrappedcauchy = list(
                    doc = "The wrapped Cauchy likelihoood",
                    hyper = list(
                        theta = list(
                            hyperid = 68001,
                            name = "log precision parameter",
                            short.name = "prec",
                            output.name = "Precision parameter for the Wrapped Cauchy observations",
                            output.name.intern = "Log precision parameter for the Wrapped Cauchy observations",
                            initial = 2,
                            fixed = FALSE,
                            prior = "loggamma",
                            param = c(1, 0.005),
                            to.theta = function(x) log(x / (1 - x)),
                            from.theta = function(x) exp(x) / (1 + exp(x))
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
                            hyperid = 69001,
                            name = "logshape",
                            short.name = "shape",
                            output.name = "Shape parameter for iid-gamma",
                            output.name.intern = "Log shape parameter for iid-gamma",
                            initial = 0,
                            fixed = FALSE,
                            prior = "loggamma",
                            param = c(100, 100),
                            to.theta = function(x) log(x),
                            from.theta = function(x) exp(x)
                        ),
                        theta2 = list(
                            hyperid = 69002,
                            name = "lograte",
                            short.name = "rate",
                            output.name = "Rate parameter for iid-gamma",
                            output.name.intern = "Log rate parameter for iid-gamma",
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
                            hyperid = 70001,
                            name = "log.a",
                            short.name = "a",
                            output.name = "a parameter for iid-beta",
                            output.name.intern = "Log a parameter for iid-beta", 
                            initial = 1,
                            fixed = FALSE,
                            prior = "loggamma",
                            param = c(1, 1),
                            to.theta = function(x) log(x),
                            from.theta = function(x) exp(x)
                        ),
                        theta2 = list(
                            hyperid = 70002,
                            name = "log.b",
                            short.name = "b",
                            output.name = "Rate parameter for iid-gamma",
                            output.name.intern = "Log rate parameter for iid-gamma",
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
                    link = c("default", "logit", "loga"),
                    pdf = "iidlogitbeta",
                    status = "experimental"
                ),

                loggammafrailty = list(
                    doc = "(experimental)",
                    hyper = list(
                        theta = list(
                            hyperid = 71001,
                            name = "log precision",
                            short.name = "prec",
                            output.name = "precision for the gamma frailty", 
                            output.name.intern = "log precision for the gamma frailty",
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
                            hyperid = 72001,
                            name = "log precision",
                            short.name = "prec",
                            output.name = "precision for the logistic observations",
                            output.name.intern = "log precision for the logistic observations", 
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

                sn = list(
                    doc = "The Skew-Normal likelihoood",
                    hyper = list(
                        theta1 = list(
                            hyperid = 74001,
                            name = "log precision",
                            short.name = "prec",
                            output.name = "precision for skew-normal observations", 
                            output.name.intern = "log precision for skew-normal observations",
                             initial = 4,
                            fixed = FALSE,
                            prior = "loggamma",
                            param = c(1, 0.00005),
                            to.theta = function(x) log(x),
                            from.theta = function(x) exp(x)
                        ),
                        theta2 = list(
                            hyperid = 74002,
                            name = "logit skew",
                            short.name = "skew",
                            output.name = "Skewness for skew-normal observations", 
                            output.name.intern = "Intern skewness for skew-normal observations",
                            initial = 0.00123456789,
                            fixed = FALSE,
                            prior = "pc.sn",
                            param = 10,
                            ## This value defined by GMRFLib_SN_SKEWMAX
                            to.theta = function(x, skew.max = 0.988) log((1 + x / skew.max) / (1 - x / skew.max)),
                            from.theta = function(x, skew.max = 0.988) skew.max * (2 * exp(x) / (1 + exp(x)) - 1)
                        )
                    ),
                    survival = FALSE,
                    discrete = FALSE,
                    link = c("default", "identity"),
                    pdf = "sn"
                ),

                gev = list(
                    doc = "The Generalized Extreme Value likelihood",
                    hyper = list(
                        theta1 = list(
                            hyperid = 76001,
                            name = "log precision",
                            short.name = "prec",
                            output.name = "precision for GEV observations", 
                            output.name.intern = "log precision for GEV observations", 
                             initial = 4,
                            fixed = FALSE,
                            prior = "loggamma",
                            param = c(1, 0.00005),
                            to.theta = function(x) log(x),
                            from.theta = function(x) exp(x)
                        ),
                        theta2 = list(
                            hyperid = 76002,
                            name = "tail parameter",
                            short.name = "tail",
                            output.name = "tail parameter for GEV observations",
                            output.name.intern = "tail parameter for GEV observations",
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
                    status = "disabled: Use likelihood model 'bgev' instead; see inla.doc('bgev')",
                    pdf = "gev"
                ),

                lognormal = list(
                    doc = "The log-Normal likelihood",
                    hyper = list(
                        theta = list(
                            hyperid = 77101,
                            name = "log precision",
                            short.name = "prec",
                            output.name = "Precision for the lognormal observations",
                            output.name.intern = "Log precision for the lognormal observations",
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
                        theta1 = list(
                            hyperid = 78001,
                            name = "log precision",
                            short.name = "prec",
                            output.name = "Precision for the lognormalsurv observations", 
                            output.name.intern = "Log precision for the lognormalsurv observations", 
                            initial = 0,
                            fixed = FALSE,
                            prior = "loggamma",
                            param = c(1, 0.00005),
                            to.theta = function(x) log(x),
                            from.theta = function(x) exp(x)
                        ),
                        theta2 = list(
                            hyperid = 78002,
                            name = "beta1",
                            short.name = "beta1",
                            output.name = "beta1 for logNormal-Cure", 
                            output.name.intern = "beta1 for logNormal-Cure", 
                            initial = -7,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(-4, 100),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta3 = list(
                            hyperid = 78003,
                            name = "beta2",
                            short.name = "beta2",
                            output.name = "beta2 for logNormal-Cure", 
                            output.name.intern = "beta2 for logNormal-Cure", 
                            initial = 0,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(0, 100),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta4 = list(
                            hyperid = 78004,
                            name = "beta3",
                            short.name = "beta3",
                            output.name = "beta3 for logNormal-Cure", 
                            output.name.intern = "beta3 for logNormal-Cure", 
                            initial = 0,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(0, 100),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta5 = list(
                            hyperid = 78005,
                            name = "beta4",
                            short.name = "beta4",
                            output.name = "beta4 for logNormal-Cure", 
                            output.name.intern = "beta4 for logNormal-Cure", 
                            initial = 0,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(0, 100),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta6 = list(
                            hyperid = 78006,
                            name = "beta5",
                            short.name = "beta5",
                            output.name = "beta5 for logNormal-Cure", 
                            output.name.intern = "beta5 for logNormal-Cure", 
                            initial = 0,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(0, 100),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta7 = list(
                            hyperid = 78007,
                            name = "beta6",
                            short.name = "beta6",
                            output.name = "beta6 for logNormal-Cure", 
                            output.name.intern = "beta6 for logNormal-Cure", 
                            initial = 0,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(0, 100),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta8 = list(
                            hyperid = 78008,
                            name = "beta7",
                            short.name = "beta7",
                            output.name = "beta7 for logNormal-Cure", 
                            output.name.intern = "beta7 for logNormal-Cure", 
                            initial = 0,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(0, 100),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta9 = list(
                            hyperid = 78009,
                            name = "beta8",
                            short.name = "beta8",
                            output.name = "beta8 for logNormal-Cure", 
                            output.name.intern = "beta8 for logNormal-Cure", 
                            initial = 0,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(0, 100),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta10 = list(
                            hyperid = 78010,
                            name = "beta9",
                            short.name = "beta9",
                            output.name = "beta9 for logNormal-Cure", 
                            output.name.intern = "beta9 for logNormal-Cure", 
                            initial = 0,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(0, 100),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta11 = list(
                            hyperid = 78011,
                            name = "beta10",
                            short.name = "beta10",
                            output.name = "beta10 for logNormal-Cure", 
                            output.name.intern = "beta10 for logNormal-Cure", 
                            initial = 0,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(0, 100),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        )
                    ),
                    survival = TRUE,
                    discrete = FALSE,
                    link = c("default", "identity"),
                    pdf = "lognormal"
                ),

                exponential = list(
                    doc = "The Exponential likelihood",
                    hyper = list(),
                    survival = FALSE,
                    discrete = FALSE,
                    link = c("default", "log"),
                    pdf = "exponential"
                ),

                exponentialsurv = list(
                    doc = "The Exponential likelihood (survival)",
                    hyper = list(
                        theta1 = list(
                            hyperid = 78020,
                            name = "beta1",
                            short.name = "beta1",
                            output.name = "beta1 for Exp-Cure", 
                            output.name.intern = "beta1 for Exp-Cure", 
                            initial = -4,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(-1, 100),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta2 = list(
                            hyperid = 78021,
                            name = "beta2",
                            short.name = "beta2",
                            output.name = "beta2 for Exp-Cure", 
                            output.name.intern = "beta2 for Exp-Cure", 
                            initial = 0,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(0, 100),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta3 = list(
                            hyperid = 78022,
                            name = "beta3",
                            short.name = "beta3",
                            output.name = "beta3 for Exp-Cure", 
                            output.name.intern = "beta3 for Exp-Cure", 
                            initial = 0,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(0, 100),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta4 = list(
                            hyperid = 78023,
                            name = "beta4",
                            short.name = "beta4",
                            output.name = "beta4 for Exp-Cure", 
                            output.name.intern = "beta4 for Exp-Cure", 
                            initial = 0,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(0, 100),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta5 = list(
                            hyperid = 78024,
                            name = "beta5",
                            short.name = "beta5",
                            output.name = "beta5 for Exp-Cure", 
                            output.name.intern = "beta5 for Exp-Cure", 
                            initial = 0,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(0, 100),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta6 = list(
                            hyperid = 78025,
                            name = "beta6",
                            short.name = "beta6",
                            output.name = "beta6 for Exp-Cure", 
                            output.name.intern = "beta6 for Exp-Cure", 
                            initial = 0,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(0, 100),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta7 = list(
                            hyperid = 78026,
                            name = "beta7",
                            short.name = "beta7",
                            output.name = "beta7 for Exp-Cure", 
                            output.name.intern = "beta7 for Exp-Cure", 
                            initial = 0,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(0, 100),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta8 = list(
                            hyperid = 78027,
                            name = "beta8",
                            short.name = "beta8",
                            output.name = "beta8 for Exp-Cure", 
                            output.name.intern = "beta8 for Exp-Cure", 
                            initial = 0,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(0, 100),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta9 = list(
                            hyperid = 78028,
                            name = "beta9",
                            short.name = "beta9",
                            output.name = "beta9 for Exp-Cure", 
                            output.name.intern = "beta9 for Exp-Cure", 
                            initial = 0,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(0, 100),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta10 = list(
                            hyperid = 78029,
                            name = "beta10",
                            short.name = "beta10",
                            output.name = "beta10 for Exp-Cure", 
                            output.name.intern = "beta10 for Exp-Cure", 
                            initial = 0,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(0, 100),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        )
                    ), 
                    survival = TRUE,
                    discrete = FALSE,
                    link = c("default", "log", "neglog"),
                    pdf = "exponential"
                ),

                coxph = list(
                    doc = "Cox-proportional hazard likelihood",
                    hyper = list(),
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
                            hyperid = 79001,
                            name = "log alpha",
                            short.name = "alpha",
                            output.name = "alpha parameter for weibull", 
                            output.name.intern = "alpha_intern for weibull", 
                            initial = -2,
                            fixed = FALSE,
                            prior = "pc.alphaw",
                            param = c(5),
                            ## the 'sc' constant is defined in inla.h, and must be the same.
                            ## I know, this is hard-coded for the moment. Should be a generic
                            ## way of doing this...
                            to.theta = function(x, sc = 0.1) log(x) / sc,
                            from.theta = function(x, sc = 0.1) exp(sc * x)
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
                            hyperid = 79101,
                            name = "log alpha",
                            short.name = "alpha",
                            output.name = "alpha parameter for weibullsurv", 
                            output.name.intern = "alpha_intern for weibullsurv", 
                            initial = -2,
                            fixed = FALSE,
                            prior = "pc.alphaw",
                            param = c(5),
                            ## the 'sc' constant is defined in inla.h, and must be the same.
                            ## I know, this is hard-coded for the moment. Should be a generic
                            ## way of doing this...
                            to.theta = function(x, sc = 0.1) log(x) / sc,
                            from.theta = function(x, sc = 0.1) exp(sc * x)
                        ),
                        theta2 = list(
                            hyperid = 79102,
                            name = "beta1",
                            short.name = "beta1",
                            output.name = "beta1 for Weibull-Cure", 
                            output.name.intern = "beta1 for Weibull-Cure", 
                            initial = -7,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(-4, 100),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta3 = list(
                            hyperid = 79103,
                            name = "beta2",
                            short.name = "beta2",
                            output.name = "beta2 for Weibull-Cure", 
                            output.name.intern = "beta2 for Weibull-Cure", 
                            initial = 0,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(0, 100),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta4 = list(
                            hyperid = 79104,
                            name = "beta3",
                            short.name = "beta3",
                            output.name = "beta3 for Weibull-Cure", 
                            output.name.intern = "beta3 for Weibull-Cure", 
                            initial = 0,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(0, 100),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta5 = list(
                            hyperid = 79105,
                            name = "beta4",
                            short.name = "beta4",
                            output.name = "beta4 for Weibull-Cure", 
                            output.name.intern = "beta4 for Weibull-Cure", 
                            initial = 0,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(0, 100),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta6 = list(
                            hyperid = 79106,
                            name = "beta5",
                            short.name = "beta5",
                            output.name = "beta5 for Weibull-Cure", 
                            output.name.intern = "beta5 for Weibull-Cure", 
                            initial = 0,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(0, 100),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta7 = list(
                            hyperid = 79107,
                            name = "beta6",
                            short.name = "beta6",
                            output.name = "beta6 for Weibull-Cure", 
                            output.name.intern = "beta6 for Weibull-Cure", 
                            initial = 0,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(0, 100),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta8 = list(
                            hyperid = 79108,
                            name = "beta7",
                            short.name = "beta7",
                            output.name = "beta7 for Weibull-Cure", 
                            output.name.intern = "beta7 for Weibull-Cure", 
                            initial = 0,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(0, 100),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta9 = list(
                            hyperid = 79109,
                            name = "beta8",
                            short.name = "beta8",
                            output.name = "beta8 for Weibull-Cure", 
                            output.name.intern = "beta8 for Weibull-Cure", 
                            initial = 0,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(0, 100),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta10 = list(
                            hyperid = 79110,
                            name = "beta9",
                            short.name = "beta9",
                            output.name = "beta9 for Weibull-Cure", 
                            output.name.intern = "beta9 for Weibull-Cure", 
                            initial = 0,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(0, 100),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta11 = list(
                            hyperid = 79111,
                            name = "beta10",
                            short.name = "beta10",
                            output.name = "beta10 for Weibull-Cure", 
                            output.name.intern = "beta10 for Weibull-Cure", 
                            initial = 0,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(0, 100),
                            to.theta = function(x) x,
                            from.theta = function(x) x
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
                            hyperid = 80001,
                            name = "log alpha",
                            short.name = "alpha",
                            output.name = "alpha for loglogistic observations", 
                            output.name.intern = "log alpha for loglogistic observations", 
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
                        theta1 = list(
                            hyperid = 80011,
                            name = "log alpha",
                            short.name = "alpha",
                            output.name = "alpha for loglogisticsurv observations", 
                            output.name.intern = "log alpha for loglogisticsurv observations", 
                            initial = 1,
                            fixed = FALSE,
                            prior = "loggamma",
                            param = c(25, 25),
                            to.theta = function(x) log(x),
                            from.theta = function(x) exp(x)
                        ),
                        theta2 = list(
                            hyperid = 80012,
                            name = "beta1",
                            short.name = "beta1",
                            output.name = "beta1 for logLogistic-Cure", 
                            output.name.intern = "beta1 for logLogistic-Cure", 
                            initial = -5,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(-4, 100),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta3 = list(
                            hyperid = 80013,
                            name = "beta2",
                            short.name = "beta2",
                            output.name = "beta2 for logLogistic-Cure", 
                            output.name.intern = "beta2 for logLogistic-Cure", 
                            initial = 0,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(0, 100),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta4 = list(
                            hyperid = 80014,
                            name = "beta3",
                            short.name = "beta3",
                            output.name = "beta3 for logLogistic-Cure", 
                            output.name.intern = "beta3 for logLogistic-Cure", 
                            initial = 0,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(0, 100),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta5 = list(
                            hyperid = 80015,
                            name = "beta4",
                            short.name = "beta4",
                            output.name = "beta4 for logLogistic-Cure", 
                            output.name.intern = "beta4 for logLogistic-Cure", 
                            initial = 0,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(0, 100),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta6 = list(
                            hyperid = 80016,
                            name = "beta5",
                            short.name = "beta5",
                            output.name = "beta5 for logLogistic-Cure", 
                            output.name.intern = "beta5 for logLogistic-Cure", 
                            initial = 0,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(0, 100),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta7 = list(
                            hyperid = 80017,
                            name = "beta6",
                            short.name = "beta6",
                            output.name = "beta6 for logLogistic-Cure", 
                            output.name.intern = "beta6 for logLogistic-Cure", 
                            initial = 0,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(0, 100),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta8 = list(
                            hyperid = 80018,
                            name = "beta7",
                            short.name = "beta7",
                            output.name = "beta7 for logLogistic-Cure", 
                            output.name.intern = "beta7 for logLogistic-Cure", 
                            initial = 0,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(0, 100),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta9 = list(
                            hyperid = 80019,
                            name = "beta8",
                            short.name = "beta8",
                            output.name = "beta8 for logLogistic-Cure", 
                            output.name.intern = "beta8 for logLogistic-Cure", 
                            initial = 0,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(0, 100),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta10 = list(
                            hyperid = 80020,
                            name = "beta9",
                            short.name = "beta9",
                            output.name = "beta9 for logLogistic-Cure", 
                            output.name.intern = "beta9 for logLogistic-Cure", 
                            initial = 0,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(0, 100),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta11 = list(
                            hyperid = 80021,
                            name = "beta10",
                            short.name = "beta10",
                            output.name = "beta10 for logLogistic-Cure", 
                            output.name.intern = "beta10 for logLogistic-Cure", 
                            initial = 0,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(0, 100),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        )
                    ),
                    survival = TRUE,
                    discrete = FALSE,
                    link = c("default", "log", "neglog"),
                    pdf = "loglogistic"
                ),

                stochvol = list(
                    doc = "The Gaussian stochvol likelihood",
                    hyper = list(
                        theta = list(
                            hyperid = 82001,
                            name = "log precision",
                            short.name = "prec",
                            output.name = "Offset precision for stochvol",
                            output.name.intern = "Log offset precision for stochvol", 
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

                stochvolln = list(
                    doc = "The Log-Normal stochvol likelihood",
                    hyper = list(
                        theta = list(
                            hyperid = 82011,
                            name = "offset",
                            short.name = "c",
                            output.name = "Mean offset for stochvolln",
                            output.name.intern = "Mean offset for stochvolln", 
                            initial = 0, 
                            fixed = FALSE,
                            prior = "normal",
                            param = c(0, 10),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        )
                    ),
                    survival = FALSE,
                    discrete = FALSE,
                    link = c("default", "log"),
                    pdf = "stochvolln"
                ),

                stochvolsn = list(
                    doc = "The SkewNormal stochvol likelihood",
                    hyper = list(
                        theta1 = list(
                            hyperid = 82101,
                            name = "logit skew",
                            short.name = "skew",
                            output.name = "Skewness for stochvol_sn observations", 
                            output.name.intern = "Intern skewness for stochvol_sn observations", 
                            initial = 0.00123456789,
                            fixed = FALSE,
                            prior = "pc.sn",
                            param = 10,
                            ## This value defined by GMRFLib_SN_SKEWMAX
                            to.theta = function(x, skew.max = 0.988) log((1 + x / skew.max) / (1 - x / skew.max)),
                            from.theta = function(x, skew.max = 0.988) skew.max * (2 * exp(x) / (1 + exp(x)) - 1)
                        ), 
                        theta2 = list(
                            hyperid = 82102,
                            name = "log precision",
                            short.name = "prec",
                            output.name = "Offset precision for stochvol_sn",
                            output.name.intern = "Log offset precision for stochvol_sn",
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
                    pdf = "stochvolsn"
                ),

                stochvolt = list(
                    doc = "The Student-t stochvol likelihood",
                    hyper = list(
                        theta = list(
                            hyperid = 83001,
                            name = "log degrees of freedom",
                            short.name = "dof",
                            output.name = "degrees of freedom for stochvol student-t",
                            output.name.intern = "dof_intern for stochvol student-t",
                            initial = 4,
                            fixed = FALSE,
                            prior = "pc.dof",
                            param = c(15, 0.5),
                            to.theta = function(x) log(x - 2),
                            from.theta = function(x) 2 + exp(x)
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
                            hyperid = 84001,
                            name = "skewness",
                            short.name = "skew",
                            output.name.intern = "skewness_param_intern for stochvol-nig",
                            output.name = "skewness parameter for stochvol-nig", 
                            initial = 0,
                            fixed = FALSE,
                            prior = "gaussian",
                            param = c(0, 10),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta2 = list(
                            hyperid = 84002,
                            name = "shape",
                            short.name = "shape",
                            output.name = "shape parameter for stochvol-nig",
                            output.name.intern = "shape_param_intern for stochvol-nig",
                            initial = 0,
                            fixed = FALSE,
                            prior = "loggamma",
                            param = c(1, 0.5),
                            to.theta = function(x) log(x - 1),
                            from.theta = function(x) 1 + exp(x)
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
                            hyperid = 85001,
                            name = "logit probability",
                            short.name = "prob",
                            output.name = "zero-probability parameter for zero-inflated poisson_0", 
                            output.name.intern = "intern zero-probability parameter for zero-inflated poisson_0",
                            initial = -1,
                            fixed = FALSE,
                            prior = "gaussian",
                            param = c(-1, 0.2),
                            to.theta = function(x) log(x / (1 - x)),
                            from.theta = function(x) exp(x) / (1 + exp(x))
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
                            hyperid = 86001,
                            name = "logit probability",
                            short.name = "prob",
                            output.name = "zero-probability parameter for zero-inflated poisson_1", 
                            output.name.intern = "intern zero-probability parameter for zero-inflated poisson_1",
                            initial = -1,
                            fixed = FALSE,
                            prior = "gaussian",
                            param = c(-1, 0.2),
                            to.theta = function(x) log(x / (1 - x)),
                            from.theta = function(x) exp(x) / (1 + exp(x))
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
                            hyperid = 87001,
                            name = "log alpha",
                            short.name = "a",
                            output.name = "zero-probability parameter for zero-inflated poisson_2", 
                            output.name.intern = "intern zero-probability parameter for zero-inflated poisson_2",
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

                zeroinflatedcenpoisson0 = list(
                    doc = "Zero-inflated censored Poisson, type 0",
                    hyper = list(
                        theta = list(
                            hyperid = 87101,
                            name = "logit probability",
                            short.name = "prob",
                            output.name = "zero-probability parameter for zero-inflated poisson_0",
                            output.name.intern = "intern zero-probability parameter for zero-inflated poisson_0",
                             initial = -1,
                            fixed = FALSE,
                            prior = "gaussian",
                            param = c(-1, 0.2),
                            to.theta = function(x) log(x / (1 - x)),
                            from.theta = function(x) exp(x) / (1 + exp(x))
                        )
                    ),
                    survival = FALSE,
                    discrete = FALSE,
                    link = c("default", "log"),
                    pdf = "zeroinflated"
                ),

                zeroinflatedcenpoisson1 = list(
                    doc = "Zero-inflated censored Poisson, type 1",
                    hyper = list(
                        theta = list(
                            hyperid = 87201,
                            name = "logit probability",
                            short.name = "prob",
                            output.name = "zero-probability parameter for zero-inflated poisson_1",
                            output.name.intern = "intern zero-probability parameter for zero-inflated poisson_1",
                            initial = -1,
                            fixed = FALSE,
                            prior = "gaussian",
                            param = c(-1, 0.2),
                            to.theta = function(x) log(x / (1 - x)),
                            from.theta = function(x) exp(x) / (1 + exp(x))
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
                            hyperid = 88001,
                            name = "overdispersion",
                            short.name = "rho",
                            output.name = "rho for zero-inflated betabinomial_0", 
                            output.name.intern = "rho_intern for zero-inflated betabinomial_0",
                            initial = 0,
                            fixed = FALSE,
                            prior = "gaussian",
                            param = c(0.0, 0.4),
                            to.theta = function(x) log(x / (1 - x)),
                            from.theta = function(x) exp(x) / (1 + exp(x))
                        ),
                        theta2 = list(
                            hyperid = 88002,
                            name = "logit probability",
                            short.name = "prob",
                            output.name = "zero-probability parameter for zero-inflated betabinomial_0", 
                            output.name.intern = "intern zero-probability parameter for zero-inflated betabinomial_0", 
                            initial = -1,
                            fixed = FALSE,
                            prior = "gaussian",
                            param = c(-1, 0.2),
                            to.theta = function(x) log(x / (1 - x)),
                            from.theta = function(x) exp(x) / (1 + exp(x))
                        )
                    ),
                    survival = FALSE,
                    discrete = TRUE,
                    link = c("default", "logit", "loga", "cauchit", "probit", "cloglog", "ccloglog", "loglog", "robit", "sn"),
                    pdf = "zeroinflated"
                ),

                zeroinflatedbetabinomial1 = list(
                    doc = "Zero-inflated Beta-Binomial, type 1",
                    hyper = list(
                        theta1 = list(
                            hyperid = 89001,
                            name = "overdispersion",
                            short.name = "rho",
                            output.name = "rho for zero-inflated betabinomial_1", 
                            output.name.intern = "rho_intern for zero-inflated betabinomial_1",
                            initial = 0,
                            fixed = FALSE,
                            prior = "gaussian",
                            param = c(0.0, 0.4),
                            to.theta = function(x) log(x / (1 - x)),
                            from.theta = function(x) exp(x) / (1 + exp(x))
                        ),
                        theta2 = list(
                            hyperid = 89002,
                            name = "logit probability",
                            short.name = "prob",
                            output.name = "zero-probability parameter for zero-inflated betabinomial_1", 
                            output.name.intern = "intern zero-probability parameter for zero-inflated betabinomial_1", 
                            initial = -1,
                            fixed = FALSE,
                            prior = "gaussian",
                            param = c(-1, 0.2),
                            to.theta = function(x) log(x / (1 - x)),
                            from.theta = function(x) exp(x) / (1 + exp(x))
                        )
                    ),
                    survival = FALSE,
                    discrete = TRUE,
                    link = c("default", "logit", "loga", "cauchit", "probit", "cloglog", "ccloglog", "loglog", "robit", "sn"),
                    pdf = "zeroinflated"
                ),

                zeroinflatedbinomial0 = list(
                    doc = "Zero-inflated Binomial, type 0",
                    hyper = list(
                        theta = list(
                            hyperid = 90001,
                            name = "logit probability",
                            short.name = "prob",
                            output.name = "zero-probability parameter for zero-inflated binomial_0", 
                            output.name.intern = "intern zero-probability parameter for zero-inflated binomial_0", 
                            initial = -1,
                            fixed = FALSE,
                            prior = "gaussian",
                            param = c(-1, 0.2),
                            to.theta = function(x) log(x / (1 - x)),
                            from.theta = function(x) exp(x) / (1 + exp(x))
                        )
                    ),
                    survival = FALSE,
                    discrete = FALSE,
                    link = c("default", "logit", "loga", "cauchit", "probit", "cloglog", "ccloglog", "loglog", "robit", "sn"),
                    pdf = "zeroinflated"
                ),

                zeroinflatedbinomial1 = list(
                    doc = "Zero-inflated Binomial, type 1",
                    hyper = list(
                        theta = list(
                            hyperid = 91001,
                            name = "logit probability",
                            short.name = "prob",
                            output.name = "zero-probability parameter for zero-inflated binomial_1", 
                            output.name.intern = "intern zero-probability parameter for zero-inflated binomial_1", 
                            initial = -1,
                            fixed = FALSE,
                            prior = "gaussian",
                            param = c(-1, 0.2),
                            to.theta = function(x) log(x / (1 - x)),
                            from.theta = function(x) exp(x) / (1 + exp(x))
                        )
                    ),
                    survival = FALSE,
                    discrete = FALSE,
                    link = c("default", "logit", "loga", "cauchit", "probit", "cloglog", "ccloglog", "loglog", "robit", "sn"),
                    pdf = "zeroinflated"
                ),

                zeroinflatedbinomial2 = list(
                    doc = "Zero-inflated Binomial, type 2",
                    hyper = list(
                        theta = list(
                            hyperid = 92001,
                            name = "alpha",
                            short.name = "alpha",
                            output.name = "zero-probability parameter for zero-inflated binomial_2", 
                            output.name.intern = "intern zero-probability parameter for zero-inflated binomial_2",
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
                    link = c("default", "logit", "loga", "cauchit", "probit", "cloglog", "ccloglog", "loglog", "robit", "sn"),
                    pdf = "zeroinflated"
                ),

                zeroninflatedbinomial2 = list(
                    doc = "Zero and N inflated binomial, type 2",
                    hyper = list(
                        theta1 = list(
                            hyperid = 93001,
                            name = "alpha1",
                            short.name = "alpha1",
                            output.name = "alpha1 parameter for zero-n-inflated binomial_2",
                            output.name.intern = "intern alpha1 parameter for zero-n-inflated binomial_2",
                            initial = -1,
                            fixed = FALSE,
                            prior = "gaussian",
                            param = c(-1, 0.2),
                            to.theta = function(x) log(x),
                            from.theta = function(x) exp(x)
                        ),
                        theta2 = list(
                            hyperid = 93002,
                            name = "alpha2",
                            short.name = "alpha2",
                            output.name = "alpha2 parameter for zero-n-inflated binomial_2", 
                            output.name.intern = "intern alpha2 parameter for zero-n-inflated binomial_2", 
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
                    link = c("default", "logit", "loga", "cauchit", "probit", "cloglog", "ccloglog", "loglog", "robit", "sn"),
                    pdf = NA
                ),

                zeroninflatedbinomial3 = list(
                    doc = "Zero and N inflated binomial, type 3",
                    hyper = list(
                        theta1 = list(
                            hyperid = 93101,
                            name = "alpha0",
                            short.name = "alpha0",
                            output.name = "alpha0 parameter for zero-n-inflated binomial_3",
                            output.name.intern = "intern alpha0 parameter for zero-n-inflated binomial_3",
                            initial = 1,
                            fixed = FALSE,
                            prior = "loggamma",
                            param = c(1, 1),
                            to.theta = function(x) log(x),
                            from.theta = function(x) exp(x)
                        ),
                        theta2 = list(
                            hyperid = 93102,
                            name = "alphaN",
                            short.name = "alphaN",
                            output.name.intern = "intern alphaN parameter for zero-n-inflated binomial_3",
                            output.name = "alphaN parameter for zero-n-inflated binomial_3", 
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
                    link = c("default", "logit", "loga", "cauchit", "probit", "cloglog", "ccloglog", "loglog", "robit", "sn"),
                    pdf = "zeroinflated"
                ),

                zeroinflatedbetabinomial2 = list(
                    doc = "Zero inflated Beta-Binomial, type 2",
                    hyper = list(
                        theta1 = list(
                            hyperid = 94001,
                            name = "log alpha",
                            short.name = "a",
                            output.name = "zero-probability parameter for zero-inflated betabinomial_2", 
                            output.name.intern = "intern zero-probability parameter for zero-inflated betabinomial_2",
                            initial = log(2),
                            fixed = FALSE,
                            prior = "gaussian",
                            param = c(log(2), 1),
                            to.theta = function(x) log(x),
                            from.theta = function(x) exp(x)
                        ),
                        theta2 = list(
                            hyperid = 94002,
                            name = "beta",
                            short.name = "b",
                            output.name = "overdispersion parameter for zero-inflated betabinomial_2",
                            output.name.intern = "intern overdispersion parameter for zero-inflated betabinomial_2",
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
                    link = c("default", "logit", "loga", "cauchit", "probit", "cloglog", "ccloglog", "loglog", "robit", "sn"),
                    pdf = "zeroinflated"
                ),

                zeroinflatednbinomial0 = list(
                    doc = "Zero inflated negBinomial, type 0",
                    hyper = list(
                        theta1 = list(
                            hyperid = 95001,
                            name = "log size",
                            short.name = "size",
                            output.name = "size for nbinomial_0 zero-inflated observations",
                            output.name.intern = "log size for nbinomial_0 zero-inflated observations",
                            initial = log(10),
                            fixed = FALSE,
                            prior = "pc.mgamma",
                            param = 7,
                            to.theta = function(x) log(x),
                            from.theta = function(x) exp(x)
                        ),
                        theta2 = list(
                            hyperid = 95002,
                            name = "logit probability",
                            short.name = "prob",
                            output.name = "zero-probability parameter for zero-inflated nbinomial_0", 
                            output.name.intern = "intern zero-probability parameter for zero-inflated nbinomial_0",
                            initial = -1,
                            fixed = FALSE,
                            prior = "gaussian",
                            param = c(-1, 0.2),
                            to.theta = function(x) log(x / (1 - x)),
                            from.theta = function(x) exp(x) / (1 + exp(x))
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
                            hyperid = 96001,
                            name = "log size",
                            short.name = "size",
                            output.name = "size for nbinomial_1 zero-inflated observations",
                            output.name.intern = "log size for nbinomial_1 zero-inflated observations",
                            initial = log(10),
                            fixed = FALSE,
                            prior = "pc.mgamma",
                            param = 7,
                            to.theta = function(x) log(x),
                            from.theta = function(x) exp(x)
                        ),
                        theta2 = list(
                            hyperid = 96002,
                            name = "logit probability",
                            short.name = "prob",
                            output.name = "zero-probability parameter for zero-inflated nbinomial_1", 
                            output.name.intern = "intern zero-probability parameter for zero-inflated nbinomial_1",
                            initial = -1,
                            fixed = FALSE,
                            prior = "gaussian",
                            param = c(-1, 0.2),
                            to.theta = function(x) log(x / (1 - x)),
                            from.theta = function(x) exp(x) / (1 + exp(x))
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
                            hyperid = 97001,
                            name = "log size",
                            short.name = "size",
                            output.name = "size for zero-inflated nbinomial_1_strata2", 
                            output.name.intern = "log size for zero-inflated nbinomial_1_strata2",
                            initial = log(10),
                            fixed = FALSE,
                            prior = "pc.mgamma",
                            param = 7,
                            to.theta = function(x) log(x),
                            from.theta = function(x) exp(x)
                        ),
                        theta2 = list(
                            hyperid = 97002,
                            name = "logit probability 1",
                            short.name = "prob1",
                            output.name = "zero-probability1 for zero-inflated nbinomial_1_strata2", 
                            output.name.intern = "intern zero-probability1 for zero-inflated nbinomial_1_strata2",
                            initial = -1,
                            fixed = FALSE,
                            prior = "gaussian",
                            param = c(-1, 0.2),
                            to.theta = function(x) log(x / (1 - x)),
                            from.theta = function(x) exp(x) / (1 + exp(x))
                        ),
                        theta3 = list(
                            hyperid = 97003,
                            name = "logit probability 2",
                            short.name = "prob2",
                            output.name = "zero-probability2 for zero-inflated nbinomial_1_strata2", 
                            output.name.intern = "intern zero-probability2 for zero-inflated nbinomial_1_strata2",
                            initial = -1,
                            fixed = FALSE,
                            prior = "gaussian",
                            param = c(-1, 0.2),
                            to.theta = function(x) log(x / (1 - x)),
                            from.theta = function(x) exp(x) / (1 + exp(x))
                        ),
                        theta4 = list(
                            hyperid = 97004,
                            name = "logit probability 3",
                            short.name = "prob3",
                            output.name = "zero-probability3 for zero-inflated nbinomial_1_strata2", 
                            output.name.intern = "intern zero-probability3 for zero-inflated nbinomial_1_strata2",
                            initial = -1,
                            fixed = TRUE,
                            prior = "gaussian",
                            param = c(-1, 0.2),
                            to.theta = function(x) log(x / (1 - x)),
                            from.theta = function(x) exp(x) / (1 + exp(x))
                        ),
                        theta5 = list(
                            hyperid = 97005,
                            name = "logit probability 4",
                            short.name = "prob4",
                            output.name = "zero-probability4 for zero-inflated nbinomial_1_strata2", 
                            output.name.intern = "intern zero-probability4 for zero-inflated nbinomial_1_strata2",
                            initial = -1,
                            fixed = TRUE,
                            prior = "gaussian",
                            param = c(-1, 0.2),
                            to.theta = function(x) log(x / (1 - x)),
                            from.theta = function(x) exp(x) / (1 + exp(x))
                        ),
                        theta6 = list(
                            hyperid = 97006,
                            name = "logit probability 5",
                            short.name = "prob5",
                            output.name = "zero-probability5 for zero-inflated nbinomial_1_strata2", 
                            output.name.intern = "intern zero-probability5 for zero-inflated nbinomial_1_strata2",
                            initial = -1,
                            fixed = TRUE,
                            prior = "gaussian",
                            param = c(-1, 0.2),
                            to.theta = function(x) log(x / (1 - x)),
                            from.theta = function(x) exp(x) / (1 + exp(x))
                        ),
                        theta7 = list(
                            hyperid = 97007,
                            name = "logit probability 6",
                            short.name = "prob6",
                            output.name = "zero-probability6 for zero-inflated nbinomial_1_strata2", 
                            output.name.intern = "intern zero-probability6 for zero-inflated nbinomial_1_strata2",
                            initial = -1,
                            fixed = TRUE,
                            prior = "gaussian",
                            param = c(-1, 0.2),
                            to.theta = function(x) log(x / (1 - x)),
                            from.theta = function(x) exp(x) / (1 + exp(x))
                        ),
                        theta8 = list(
                            hyperid = 97008,
                            name = "logit probability 7",
                            short.name = "prob7",
                            output.name = "zero-probability7 for zero-inflated nbinomial_1_strata2", 
                            output.name.intern = "intern zero-probability7 for zero-inflated nbinomial_1_strata2",
                            initial = -1,
                            fixed = TRUE,
                            prior = "gaussian",
                            param = c(-1, 0.2),
                            to.theta = function(x) log(x / (1 - x)),
                            from.theta = function(x) exp(x) / (1 + exp(x))
                        ),
                        theta9 = list(
                            hyperid = 97009,
                            name = "logit probability 8",
                            short.name = "prob8",
                            output.name = "zero-probability8 for zero-inflated nbinomial_1_strata2", 
                            output.name.intern = "intern zero-probability8 for zero-inflated nbinomial_1_strata2",
                            initial = -1,
                            fixed = TRUE,
                            prior = "gaussian",
                            param = c(-1, 0.2),
                            to.theta = function(x) log(x / (1 - x)),
                            from.theta = function(x) exp(x) / (1 + exp(x))
                        ),
                        theta10 = list(
                            hyperid = 97010,
                            name = "logit probability 9",
                            short.name = "prob9",
                            output.name = "zero-probability9 for zero-inflated nbinomial_1_strata2", 
                            output.name.intern = "intern zero-probability9 for zero-inflated nbinomial_1_strata2",
                            initial = -1,
                            fixed = TRUE,
                            prior = "gaussian",
                            param = c(-1, 0.2),
                            to.theta = function(x) log(x / (1 - x)),
                            from.theta = function(x) exp(x) / (1 + exp(x))
                        ),
                        theta11 = list(
                            hyperid = 97011,
                            name = "logit probability 10",
                            short.name = "prob10",
                            output.name = "zero-probability10 for zero-inflated nbinomial_1_strata2", 
                            output.name.intern = "intern zero-probability10 for zero-inflated nbinomial_1_strata2",
                            initial = -1,
                            fixed = TRUE,
                            prior = "gaussian",
                            param = c(-1, 0.2),
                            to.theta = function(x) log(x / (1 - x)),
                            from.theta = function(x) exp(x) / (1 + exp(x))
                        )
                    ),
                    survival = FALSE,
                    discrete = FALSE,
                    link = c("default", "log"),
                    pdf = "zeroinflated"
                ),

                zeroinflatednbinomial1strata3 = list(
                    doc = "Zero inflated negBinomial, type 1, strata 3",
                    hyper = list(
                        theta1 = list(
                            hyperid = 98001,
                            name = "logit probability",
                            short.name = "prob",
                            output.name = "zero-probability for zero-inflated nbinomial_1_strata3", 
                            output.name.intern = "intern zero-probability for zero-inflated nbinomial_1_strata3",
                            initial = -1,
                            fixed = FALSE,
                            prior = "gaussian",
                            param = c(-1, 0.2),
                            to.theta = function(x) log(x / (1 - x)),
                            from.theta = function(x) exp(x) / (1 + exp(x))
                        ),
                        theta2 = list(
                            hyperid = 98002,
                            name = "log size 1",
                            short.name = "size1",
                            output.name = "size1 for zero-inflated nbinomial_1_strata3",
                            output.name.intern = "log_size1 for zero-inflated nbinomial_1_strata3",
                            initial = log(10),
                            fixed = FALSE,
                            prior = "pc.mgamma",
                            param = 7,
                            to.theta = function(x) log(x),
                            from.theta = function(x) exp(x)
                        ),
                        theta3 = list(
                            hyperid = 98003,
                            name = "log size 2",
                            short.name = "size2",
                            output.name = "size2 for zero-inflated nbinomial_1_strata3",
                            output.name.intern = "log_size2 for zero-inflated nbinomial_1_strata3",
                            initial = log(10),
                            fixed = FALSE,
                            prior = "pc.mgamma",
                            param = 7,
                            to.theta = function(x) log(x),
                            from.theta = function(x) exp(x)
                        ),
                        theta4 = list(
                            hyperid = 98004,
                            name = "log size 3",
                            short.name = "size3",
                            output.name = "size3 for zero-inflated nbinomial_1_strata3",
                            output.name.intern = "log_size3 for zero-inflated nbinomial_1_strata3",
                            initial = log(10),
                            fixed = TRUE,
                            prior = "pc.mgamma",
                            param = 7,
                            to.theta = function(x) log(x),
                            from.theta = function(x) exp(x)
                        ),
                        theta5 = list(
                            hyperid = 98005,
                            name = "log size 4",
                            short.name = "size4",
                            output.name = "size4 for zero-inflated nbinomial_1_strata3",
                            output.name.intern = "log_size4 for zero-inflated nbinomial_1_strata3",
                            initial = log(10),
                            fixed = TRUE,
                            prior = "pc.mgamma",
                            param = 7,
                            to.theta = function(x) log(x),
                            from.theta = function(x) exp(x)
                        ),
                        theta6 = list(
                            hyperid = 98006,
                            name = "log size 5",
                            short.name = "size5",
                            output.name = "size5 for zero-inflated nbinomial_1_strata3",
                            output.name.intern = "log_size5 for zero-inflated nbinomial_1_strata3",
                            initial = log(10),
                            fixed = TRUE,
                            prior = "pc.mgamma",
                            param = 7,
                            to.theta = function(x) log(x),
                            from.theta = function(x) exp(x)
                        ),
                        theta7 = list(
                            hyperid = 98007,
                            name = "log size 6",
                            short.name = "size6",
                            output.name = "size6 for zero-inflated nbinomial_1_strata3",
                            output.name.intern = "log_size6 for zero-inflated nbinomial_1_strata3",
                            initial = log(10),
                            fixed = TRUE,
                            prior = "pc.mgamma",
                            param = 7,
                            to.theta = function(x) log(x),
                            from.theta = function(x) exp(x)
                        ),
                        theta8 = list(
                            hyperid = 98008,
                            name = "log size 7",
                            short.name = "size7",
                            output.name = "size7 for zero-inflated nbinomial_1_strata3",
                            output.name.intern = "log_size7 for zero-inflated nbinomial_1_strata3",
                            initial = log(10),
                            fixed = TRUE,
                            prior = "pc.mgamma",
                            param = 7,
                            to.theta = function(x) log(x),
                            from.theta = function(x) exp(x)
                        ),
                        theta9 = list(
                            hyperid = 98009,
                            name = "log size 8",
                            short.name = "size8",
                            output.name = "size8 for zero-inflated nbinomial_1_strata3",
                            output.name.intern = "log_size8 for zero-inflated nbinomial_1_strata3",
                            initial = log(10),
                            fixed = TRUE,
                            prior = "pc.mgamma",
                            param = 7,
                            to.theta = function(x) log(x),
                            from.theta = function(x) exp(x)
                        ),
                        theta10 = list(
                            hyperid = 98010,
                            name = "log size 9",
                            short.name = "size9",
                            output.name = "size9 for zero-inflated nbinomial_1_strata3",
                            output.name.intern = "log_size9 for zero-inflated nbinomial_1_strata3",
                            initial = log(10),
                            fixed = TRUE,
                            prior = "pc.mgamma",
                            param = 7,
                            to.theta = function(x) log(x),
                            from.theta = function(x) exp(x)
                        ),
                        theta11 = list(
                            hyperid = 98011,
                            name = "log size 10",
                            short.name = "size10",
                            output.name = "size10 for zero-inflated nbinomial_1_strata3",
                            output.name.intern = "log_size10 for zero-inflated nbinomial_1_strata3",
                            initial = log(10),
                            fixed = TRUE,
                            prior = "pc.mgamma",
                            param = 7,
                            to.theta = function(x) log(x),
                            from.theta = function(x) exp(x)
                        )
                    ),
                    survival = FALSE,
                    discrete = FALSE,
                    link = c("default", "log"),
                    pdf = "zeroinflated"
                ),

                zeroinflatednbinomial2 = list(
                    doc = "Zero inflated negBinomial, type 2",
                    hyper = list(
                        theta1 = list(
                            hyperid = 99001,
                            name = "log size",
                            short.name = "size",
                            output.name = "size for nbinomial zero-inflated observations",
                            output.name.inter = "log size for nbinomial zero-inflated observations",
                            initial = log(10),
                            fixed = FALSE,
                            prior = "pc.mgamma",
                            param = 7,
                            to.theta = function(x) log(x),
                            from.theta = function(x) exp(x)
                        ),
                        theta2 = list(
                            hyperid = 99002,
                            name = "log alpha",
                            short.name = "a",
                            output.name = "parameter alpha for zero-inflated nbinomial2",
                            output.name.intern = "parameter alpha.intern for zero-inflated nbinomial2",
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
                            hyperid = 100001,
                            name = "log precision",
                            short.name = "prec",
                            output.name = "precision for the student-t observations", 
                            output.name.intern = "log precision for the student-t observations", 
                            initial = 0,
                            fixed = FALSE,
                            prior = "loggamma",
                            param = c(1, 0.00005),
                            to.theta = function(x) log(x),
                            from.theta = function(x) exp(x)
                        ),
                        theta2 = list(
                            hyperid = 100002,
                            name = "log degrees of freedom",
                            short.name = "dof",
                            output.name = "degrees of freedom for student-t",
                            output.name.intern = "dof_intern for student-t", 
                            initial = 5,
                            fixed = FALSE,
                            prior = "pc.dof",
                            param = c(15, 0.5),
                            to.theta = function(x) log(x - 2),
                            from.theta = function(x) 2 + exp(x)
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
                            hyperid = 101001,
                            name = "log degrees of freedom",
                            short.name = "dof",
                            output.name.intern = "dof_intern for tstrata",
                            output.name = "degrees of freedom for tstrata",
                            initial = 4,
                            fixed = FALSE,
                            prior = "pc.dof",
                            param = c(15, 0.5),
                            to.theta = function(x) log(x - 5),
                            from.theta = function(x) 5 + exp(x)
                        ),
                        theta2 = list(
                            hyperid = 101002,
                            name = "log precision1",
                            short.name = "prec1",
                            output.name = "Prec for tstrata strata",
                            output.name.intern = "Log prec for tstrata strata",
                            initial = 2,
                            fixed = FALSE,
                            prior = "loggamma",
                            param = c(1, 0.00005),
                            to.theta = function(x) log(x),
                            from.theta = function(x) exp(x)
                        ),
                        theta3 = list(
                            hyperid = 101003,
                            name = "log precision2",
                            short.name = "prec2",
                            output.name = "Prec for tstrata strata[2]",
                            output.name.intern = "Log prec for tstrata strata[2]",
                            initial = 2,
                            fixed = FALSE,
                            prior = "loggamma",
                            param = c(1, 0.00005),
                            to.theta = function(x) log(x),
                            from.theta = function(x) exp(x)
                        ),
                        theta4 = list(
                            hyperid = 101004,
                            name = "log precision3",
                            short.name = "prec3",
                            output.name = "Prec for tstrata strata[3]",
                            output.name.intern = "Log prec for tstrata strata[3]",
                            initial = 2,
                            fixed = FALSE,
                            prior = "loggamma",
                            param = c(1, 0.00005),
                            to.theta = function(x) log(x),
                            from.theta = function(x) exp(x)
                        ),
                        theta5 = list(
                            hyperid = 101005,
                            name = "log precision4",
                            short.name = "prec4",
                            output.name = "Prec for tstrata strata[4]",
                            output.name.intern = "Log prec for tstrata strata[4]",
                            initial = 2,
                            fixed = FALSE,
                            prior = "loggamma",
                            param = c(1, 0.00005),
                            to.theta = function(x) log(x),
                            from.theta = function(x) exp(x)
                        ),
                        theta6 = list(
                            hyperid = 101006,
                            name = "log precision5",
                            short.name = "prec5",
                            output.name = "Prec for tstrata strata[5]",
                            output.name.intern = "Log prec for tstrata strata[5]",
                            initial = 2,
                            fixed = FALSE,
                            prior = "loggamma",
                            param = c(1, 0.00005),
                            to.theta = function(x) log(x),
                            from.theta = function(x) exp(x)
                        ),
                        theta7 = list(
                            hyperid = 101007,
                            name = "log precision6",
                            short.name = "prec6",
                            output.name = "Prec for tstrata strata[6]",
                            output.name.intern = "Log prec for tstrata strata[6]",
                            initial = 2,
                            fixed = FALSE,
                            prior = "loggamma",
                            param = c(1, 0.00005),
                            to.theta = function(x) log(x),
                            from.theta = function(x) exp(x)
                        ),
                        theta8 = list(
                            hyperid = 101008,
                            name = "log precision7",
                            short.name = "prec7",
                            output.name = "Prec for tstrata strata[7]",
                            output.name.intern = "Log prec for tstrata strata[7]",
                            initial = 2,
                            fixed = FALSE,
                            prior = "loggamma",
                            param = c(1, 0.00005),
                            to.theta = function(x) log(x),
                            from.theta = function(x) exp(x)
                        ),
                        theta9 = list(
                            hyperid = 101009,
                            name = "log precision8",
                            short.name = "prec8",
                            output.name = "Prec for tstrata strata[8]",
                            output.name.intern = "Log prec for tstrata strata[8]",
                            initial = 2,
                            fixed = FALSE,
                            prior = "loggamma",
                            param = c(1, 0.00005),
                            to.theta = function(x) log(x),
                            from.theta = function(x) exp(x)
                        ),
                        theta10 = list(
                            hyperid = 101010,
                            name = "log precision9",
                            short.name = "prec9",
                            output.name = "Prec for tstrata strata[9]",
                            output.name.intern = "Log prec for tstrata strata[9]",
                            initial = 2,
                            fixed = FALSE,
                            prior = "loggamma",
                            param = c(1, 0.00005),
                            to.theta = function(x) log(x),
                            from.theta = function(x) exp(x)
                        ),
                        theta11 = list(
                            hyperid = 101011,
                            name = "log precision10",
                            short.name = "prec10",
                            output.name = "Prec for tstrata strata[10]",
                            output.name.intern = "Log prec for tstrata strata[10]",
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
                            hyperid = 101101,
                            name = "beta1",
                            short.name = "beta1",
                            output.name = "beta[1] for NMix observations",
                            output.name.intern = "beta[1] for NMix observations",
                            initial = log(10),
                            fixed = FALSE,
                            prior = "normal",
                            param = c(0, 0.5),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta2 = list(
                            hyperid = 101102,
                            name = "beta2",
                            short.name = "beta2",
                            output.name = "beta[2] for NMix observations",
                            output.name.intern = "beta[2] for NMix observations",
                            initial = 0,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(0, 1),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta3 = list(
                            hyperid = 101103,
                            name = "beta3",
                            short.name = "beta3",
                            output.name = "beta[3] for NMix observations",
                            output.name.intern = "beta[3] for NMix observations",
                            initial = 0,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(0, 1),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta4 = list(
                            hyperid = 101104,
                            name = "beta4",
                            short.name = "beta4",
                            output.name = "beta[4] for NMix observations",
                            output.name.intern = "beta[4] for NMix observations",
                            initial = 0,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(0, 1),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta5 = list(
                            hyperid = 101105,
                            name = "beta5",
                            short.name = "beta5",
                            output.name = "beta[5] for NMix observations",
                            output.name.intern = "beta[5] for NMix observations",
                            initial = 0,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(0, 1),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta6 = list(
                            hyperid = 101106,
                            name = "beta6",
                            short.name = "beta6",
                            output.name = "beta[6] for NMix observations",
                            output.name.intern = "beta[6] for NMix observations",
                            initial = 0,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(0, 1),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta7 = list(
                            hyperid = 101107,
                            name = "beta7",
                            short.name = "beta7",
                            output.name = "beta[7] for NMix observations",
                            output.name.intern = "beta[7] for NMix observations",
                            initial = 0,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(0, 1),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta8 = list(
                            hyperid = 101108,
                            name = "beta8",
                            short.name = "beta8",
                            output.name = "beta[8] for NMix observations",
                            output.name.intern = "beta[8] for NMix observations",
                            initial = 0,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(0, 1),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta9 = list(
                            hyperid = 101109,
                            name = "beta9",
                            short.name = "beta9",
                            output.name = "beta[9] for NMix observations",
                            output.name.intern = "beta[9] for NMix observations",
                            initial = 0,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(0, 1),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta10 = list(
                            hyperid = 101110,
                            name = "beta10",
                            short.name = "beta10",
                            output.name = "beta[10] for NMix observations",
                            output.name.intern = "beta[10] for NMix observations",
                            initial = 0,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(0, 1),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta11 = list(
                            hyperid = 101111,
                            name = "beta11",
                            short.name = "beta11",
                            output.name = "beta[11] for NMix observations",
                            output.name.intern = "beta[11] for NMix observations",
                            initial = 0,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(0, 1),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta12 = list(
                            hyperid = 101112,
                            name = "beta12",
                            short.name = "beta12",
                            output.name = "beta[12] for NMix observations",
                            output.name.intern = "beta[12] for NMix observations",
                            initial = 0,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(0, 1),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta13 = list(
                            hyperid = 101113,
                            name = "beta13",
                            short.name = "beta13",
                            output.name = "beta[13] for NMix observations",
                            output.name.intern = "beta[13] for NMix observations",
                            initial = 0,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(0, 1),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta14 = list(
                            hyperid = 101114,
                            name = "beta14",
                            short.name = "beta14",
                            output.name = "beta[14] for NMix observations",
                            output.name.intern = "beta[14] for NMix observations",
                            initial = 0,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(0, 1),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta15 = list(
                            hyperid = 101115,
                            name = "beta15",
                            short.name = "beta15",
                            output.name = "beta[15] for NMix observations",
                            output.name.intern = "beta[15] for NMix observations",
                            initial = 0,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(0, 1),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        )
                    ),
                    survival = FALSE,
                    discrete = TRUE,
                    link = c("default", "logit", "loga", "probit"),
                    pdf = "nmix"
                ),

                nmixnb = list(
                    doc = "NegBinomial-Poisson mixture",
                    hyper = list(
                        theta1 = list(
                            hyperid = 101121,
                            name = "beta1",
                            short.name = "beta1",
                            output.name = "beta[1] for NMixNB observations",
                            output.name.intern = "beta[1] for NMixNB observations",
                            initial = log(10),
                            fixed = FALSE,
                            prior = "normal",
                            param = c(0, 0.5),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta2 = list(
                            hyperid = 101122,
                            name = "beta2",
                            short.name = "beta2",
                            output.name = "beta[2] for NMixNB observations",
                            output.name.intern = "beta[2] for NMixNB observations",
                            initial = 0,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(0, 1),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta3 = list(
                            hyperid = 101123,
                            name = "beta3",
                            short.name = "beta3",
                            output.name = "beta[3] for NMixNB observations",
                            output.name.intern = "beta[3] for NMixNB observations",
                            initial = 0,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(0, 1),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta4 = list(
                            hyperid = 101124,
                            name = "beta4",
                            short.name = "beta4",
                            output.name = "beta[4] for NMixNB observations",
                            output.name.intern = "beta[4] for NMixNB observations",
                            initial = 0,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(0, 1),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta5 = list(
                            hyperid = 101125,
                            name = "beta5",
                            short.name = "beta5",
                            output.name = "beta[5] for NMixNB observations",
                            output.name.intern = "beta[5] for NMixNB observations",
                            initial = 0,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(0, 1),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta6 = list(
                            hyperid = 101126,
                            name = "beta6",
                            short.name = "beta6",
                            output.name = "beta[6] for NMixNB observations",
                            output.name.intern = "beta[6] for NMixNB observations",
                            initial = 0,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(0, 1),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta7 = list(
                            hyperid = 101127,
                            name = "beta7",
                            short.name = "beta7",
                            output.name = "beta[7] for NMixNB observations",
                            output.name.intern = "beta[7] for NMixNB observations",
                            initial = 0,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(0, 1),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta8 = list(
                            hyperid = 101128,
                            name = "beta8",
                            short.name = "beta8",
                            output.name = "beta[8] for NMixNB observations",
                            output.name.intern = "beta[8] for NMixNB observations",
                            initial = 0,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(0, 1),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta9 = list(
                            hyperid = 101129,
                            name = "beta9",
                            short.name = "beta9",
                            output.name = "beta[9] for NMixNB observations",
                            output.name.intern = "beta[9] for NMixNB observations",
                            initial = 0,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(0, 1),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta10 = list(
                            hyperid = 101130,
                            name = "beta10",
                            short.name = "beta10",
                            output.name = "beta[10] for NMixNB observations",
                            output.name.intern = "beta[10] for NMixNB observations",
                            initial = 0,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(0, 1),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta11 = list(
                            hyperid = 101131,
                            name = "beta11",
                            short.name = "beta11",
                            output.name = "beta[11] for NMixNB observations",
                            output.name.intern = "beta[11] for NMixNB observations",
                            initial = 0,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(0, 1),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta12 = list(
                            hyperid = 101132,
                            name = "beta12",
                            short.name = "beta12",
                            output.name = "beta[12] for NMixNB observations",
                            output.name.intern = "beta[12] for NMixNB observations",
                            initial = 0,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(0, 1),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta13 = list(
                            hyperid = 101133,
                            name = "beta13",
                            short.name = "beta13",
                            output.name = "beta[13] for NMixNB observations",
                            output.name.intern = "beta[13] for NMixNB observations",
                            initial = 0,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(0, 1),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta14 = list(
                            hyperid = 101134,
                            name = "beta14",
                            short.name = "beta14",
                            output.name = "beta[14] for NMixNB observations",
                            output.name.intern = "beta[14] for NMixNB observations",
                            initial = 0,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(0, 1),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta15 = list(
                            hyperid = 101135,
                            name = "beta15",
                            short.name = "beta15",
                            output.name = "beta[15] for NMixNB observations",
                            output.name.intern = "beta[15] for NMixNB observations",
                            initial = 0,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(0, 1),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta16 = list(
                            hyperid = 101136,
                            name = "overdispersion",
                            short.name = "overdispersion",
                            output.name = "overdispersion for NMixNB observations",
                            output.name.intern = "log_overdispersion for NMixNB observations",
                            initial = 0,
                            fixed = FALSE,
                            prior = "pc.gamma",
                            param = 7,
                            to.theta = function(x) log(x),
                            from.theta = function(x) exp(x)
                        )
                    ),
                    survival = FALSE,
                    discrete = TRUE,
                    link = c("default", "logit", "loga", "probit"),
                    pdf = "nmixnb"
                ),

                gp = list(
                    doc = "Generalized Pareto likelihood",
                    hyper = list(
                        theta = list(
                            hyperid = 101201,
                            name = "tail",
                            short.name = "xi",
                            output.name = "Tail parameter for the gp observations",
                            output.name.intern = "Intern tail parameter for the gp observations",
                            initial = -4,
                            fixed = FALSE,
                            prior = "pc.gevtail",
                            param = c(7, 0.0, 0.5),
                            to.theta = function(x, interval = c(REPLACE.ME.low, REPLACE.ME.high)) log(-(interval[1] - x) / (interval[2] - x)),
                            from.theta = function(x, interval = c(REPLACE.ME.low, REPLACE.ME.high)) interval[1] + (interval[2] - interval[1]) * exp(x) / (1.0 + exp(x))
                        )
                    ),
                    survival = FALSE,
                    discrete = TRUE,
                    link = c("default", "quantile"),
                    pdf = "genPareto"
                ),

                dgp = list(
                    doc = "Discrete generalized Pareto likelihood",
                    hyper = list(
                        theta = list(
                            hyperid = 101301,
                            name = "tail",
                            short.name = "xi",
                            output.name = "Tail parameter for the dgp observations",
                            output.name.intern = "Intern tail parameter for the dgp observations",
                            initial = 2,
                            fixed = FALSE,
                            prior = "pc.gevtail",
                            param = c(7, 0.0, 0.5),
                            to.theta = function(x, interval = c(REPLACE.ME.low, REPLACE.ME.high)) log(-(interval[1] - x) / (interval[2] - x)),
                            from.theta = function(x, interval = c(REPLACE.ME.low, REPLACE.ME.high)) interval[1] + (interval[2] - interval[1]) * exp(x) / (1.0 + exp(x))
                        )
                    ),
                    survival = FALSE,
                    discrete = TRUE,
                    link = c("default", "quantile"),
                    pdf = "dgp"
                ),

                logperiodogram = list(
                    doc = "Likelihood for the log-periodogram",
                    hyper = list(),
                    survival = FALSE,
                    discrete = FALSE,
                    link = c("default", "identity"),
                    pdf = NA
                ),

                tweedie = list(
                    doc = "Tweedie distribution",
                    hyper = list(
                        theta1 = list(
                            hyperid = 102101,
                            name = "p",
                            short.name = "p",
                            output.name = "p parameter for Tweedie",
                            output.name.intern = "p_intern parameter for Tweedie",
                            initial = 0,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(0, 100),
                            to.theta = function(x, interval = c(1.0, 2.0)) log(-(interval[1] - x) / (interval[2] - x)),
                            from.theta = function(x, interval = c(1.0, 2.0)) interval[1] + (interval[2] - interval[1]) * exp(x) / (1.0 + exp(x))
                        ),
                        theta2 = list(
                            hyperid = 102201,
                            name = "dispersion",
                            short.name = "phi",
                            output.name = "Dispersion parameter for Tweedie",
                            output.name.intern = "Log dispersion parameter for Tweedie",
                            initial = -4,
                            fixed = FALSE,
                            prior = "loggamma",
                            param = c(100, 100),
                            to.theta = function(x) log(x),
                            from.theta = function(x) exp(x)
                        )
                    ),
                    survival = FALSE,
                    discrete = FALSE,
                    link = c("default", "log"),
                    pdf = "tweedie"
                ),

                fmri = list(
                    doc = "fmri distribution (special nc-chi)",
                    hyper = list(
                        theta1 = list(
                            hyperid = 103101,
                            name = "precision",
                            short.name = "prec",
                            output.name = "Precision for fmri", 
                            output.name.intern = "Log precision for fmri", 
                            initial = 0,
                            fixed = FALSE,
                            prior = "loggamma",
                            param = c(10, 10),
                            to.theta = function(x) log(x), 
                            from.theta = function(x) exp(x)
                        ),
                        theta2 = list(
                            ## this parameter is ment to be fixed. this is why the identity
                            ## mapping is used. an error is thown if fixed=FALSE is set.
                            hyperid = 103202,
                            name = "dof",
                            short.name = "df",
                            output.name = "NOT IN USE",
                            output.name.intern = "NOT IN USE",
                            initial = 4,
                            fixed = TRUE,
                            prior = "normal",
                            param = c(0, 1),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        )
                    ),
                    survival = FALSE,
                    discrete = FALSE,
                    link = c("default", "log"),
                    pdf = "fmri"
                ), 

                fmrisurv = list(
                    doc = "fmri distribution (special nc-chi)",
                    hyper = list(
                        theta1 = list(
                            hyperid = 104101,
                            name = "precision",
                            short.name = "prec",
                            output.name = "Precision for fmrisurv", 
                            output.name.intern = "Log precision for fmrisurv", 
                            initial = 0,
                            fixed = FALSE,
                            prior = "loggamma",
                            param = c(10, 10),
                            to.theta = function(x) log(x), 
                            from.theta = function(x) exp(x)
                        ),
                        theta2 = list(
                            ## this parameter is ment to be fixed. this is why the identity
                            ## mapping is used. an error is thown if fixed=FALSE is set.
                            hyperid = 104201,
                            name = "dof",
                            short.name = "df",
                            output.name = "NOT IN USE",
                            output.name.intern = "NOT IN USE",
                            initial = 4,
                            fixed = TRUE,
                            prior = "normal",
                            param = c(0, 1),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        )
                    ),
                    survival = TRUE,
                    discrete = FALSE,
                    link = c("default", "log"),
                    pdf = "fmri"
                ),

                gompertz = list(
                    doc = "gompertz distribution",
                    hyper = list(
                        theta = list(
                            hyperid = 105101,
                            name = "shape",
                            short.name = "alpha",
                            output.name.intern = "alpha_intern for Gompertz", 
                            output.name = "alpha parameter for Gompertz", 
                            initial = -1,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(0, 1),
                            ## the 'sc' constant is defined in inla.h, and must be the same.
                            ## I know, this is hard-coded for the moment. Should be a generic
                            ## way of doing this...
                            to.theta = function(x, sc = 0.1) log(x) / sc,
                            from.theta = function(x, sc = 0.1) exp(sc * x)
                        )
                    ),
                    survival = FALSE,
                    discrete = FALSE,
                    link = c("default", "log", "neglog"),
                    pdf = "gompertz"
                ), 

                gompertzsurv = list(
                    doc = "gompertz distribution",
                    hyper = list(
                        theta1 = list(
                            hyperid = 106101,
                            name = "shape",
                            short.name = "alpha",
                            output.name.intern = "alpha_intern for Gompertz-surv", 
                            output.name = "alpha parameter for Gompertz-surv", 
                            initial = -10,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(0, 1),
                            ## the 'sc' constant is defined in inla.h, and must be the same.
                            ## I know, this is hard-coded for the moment. Should be a generic
                            ## way of doing this...
                            to.theta = function(x, sc = 0.1) log(x) / sc,
                            from.theta = function(x, sc = 0.1) exp(sc * x)
                        ),
                        theta2 = list(
                            hyperid = 106102,
                            name = "beta1",
                            short.name = "beta1",
                            output.name = "beta1 for Gompertz-Cure", 
                            output.name.intern = "beta1 for Gompertz-Cure", 
                            initial = -5,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(-4, 100),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta3 = list(
                            hyperid = 106103,
                            name = "beta2",
                            short.name = "beta2",
                            output.name = "beta2 for Gompertz-Cure", 
                            output.name.intern = "beta2 for Gompertz-Cure", 
                            initial = 0,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(0, 100),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta4 = list(
                            hyperid = 106104,
                            name = "beta3",
                            short.name = "beta3",
                            output.name = "beta3 for Gompertz-Cure", 
                            output.name.intern = "beta3 for Gompertz-Cure", 
                            initial = 0,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(0, 100),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta5 = list(
                            hyperid = 106105,
                            name = "beta4",
                            short.name = "beta4",
                            output.name = "beta4 for Gompertz-Cure", 
                            output.name.intern = "beta4 for Gompertz-Cure", 
                            initial = 0,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(0, 100),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta6 = list(
                            hyperid = 106106,
                            name = "beta5",
                            short.name = "beta5",
                            output.name = "beta5 for Gompertz-Cure", 
                            output.name.intern = "beta5 for Gompertz-Cure", 
                            initial = 0,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(0, 100),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta7 = list(
                            hyperid = 106107,
                            name = "beta6",
                            short.name = "beta6",
                            output.name = "beta6 for Gompertz-Cure", 
                            output.name.intern = "beta6 for Gompertz-Cure", 
                            initial = 0,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(0, 100),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta8 = list(
                            hyperid = 106108,
                            name = "beta7",
                            short.name = "beta7",
                            output.name = "beta7 for Gompertz-Cure", 
                            output.name.intern = "beta7 for Gompertz-Cure", 
                            initial = 0,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(0, 100),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta9 = list(
                            hyperid = 106109,
                            name = "beta8",
                            short.name = "beta8",
                            output.name = "beta8 for Gompertz-Cure", 
                            output.name.intern = "beta8 for Gompertz-Cure", 
                            initial = 0,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(0, 100),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta10 = list(
                            hyperid = 106110,
                            name = "beta9",
                            short.name = "beta9",
                            output.name = "beta9 for Gompertz-Cure", 
                            output.name.intern = "beta9 for Gompertz-Cure", 
                            initial = 0,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(0, 100),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta11 = list(
                            hyperid = 106111,
                            name = "beta10",
                            short.name = "beta10",
                            output.name = "beta10 for Gompertz-Cure", 
                            output.name.intern = "beta10 for Gompertz-Cure", 
                            initial = 0,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(0, 100),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        )
                    ),
                    survival = TRUE,
                    discrete = FALSE,
                    link = c("default", "log", "neglog"),
                    pdf = "gompertz"
                ) 
            )
    )
}

`inla.models.section.prior` <- function() {
    list(
        prior =
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
                laplace = list(
                    doc = "Laplace prior",
                    nparameters = 2L,
                    pdf = "laplace"
                ),
                linksnintercept = list(
                    doc = "Skew-normal-link intercept-prior",
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
                flat = list(
                    doc = "A constant prior",
                    nparameters = 0L,
                    pdf = "various-flat"
                ),
                logflat = list(
                    doc = "A constant prior for log(theta)",
                    nparameters = 0L,
                    pdf = "various-flat"
                ),
                logiflat = list(
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
                    doc = "Dirichlet prior",
                    nparameters = 1L,
                    pdf = "dirichlet"
                ),

                ## this is the 'no prior needed' prior
                none = list(
                    doc = "No prior",
                    nparameters = 0L
                ),

                ## this is the 'flag an error if used' prior
                invalid = list(
                    doc = "Void prior",
                    nparameters = 0L
                ),

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

                ## this is the generic one, which is case-specific and possibly adaptive
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

                wishartkd =  list(
                    doc = "Wishart prior",
                    nparameters = 324L,
                    pdf = NULL
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
                ),
                
                "rprior:" = list(
                    doc = "A R-function defining the prior",
                    status = "experimental", 
                    nparameters = 0L,
                    pdf = "rprior"
                )
            )
    )
}

`inla.models.section.wrapper` <- function() {
    list(
        wrapper =
            list(
                joint = list(
                    doc = "(experimental)",
                    hyper = list(
                        theta = list(
                            hyperid = 102001,
                            name = "log precision",
                            short.name = "prec",
                            output.name = "NOT IN USE",
                            output.name.intern = "NOT IN USE",
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

`inla.models.section.lp.scale` <- function() {
    list(
        lp.scale =
            list(
                lp.scale = list(
                    hyper = list(
                        theta1 = list(
                            hyperid = 103001,
                            name = "beta1",
                            short.name = "b1",
                            output.name = "beta[1] for lp_scale", 
                            output.name.intern = "beta[1] for lp_scale", 
                            initial = 1,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(1, 10),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ), 
                        theta2 = list(
                            hyperid = 103002,
                            name = "beta2",
                            short.name = "b2",
                            output.name = "beta[2] for lp_scale", 
                            output.name.intern = "beta[2] for lp_scale", 
                            initial = 1,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(1, 10),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ), 
                        theta3 = list(
                            hyperid = 103003,
                            name = "beta3",
                            short.name = "b3",
                            output.name = "beta[3] for lp_scale", 
                            output.name.intern = "beta[3] for lp_scale", 
                            initial = 1,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(1, 10),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ), 
                        theta4 = list(
                            hyperid = 103004,
                            name = "beta4",
                            short.name = "b4",
                            output.name = "beta[4] for lp_scale", 
                            output.name.intern = "beta[4] for lp_scale", 
                            initial = 1,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(1, 10),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ), 
                        theta5 = list(
                            hyperid = 103005,
                            name = "beta5",
                            short.name = "b5",
                            output.name = "beta[5] for lp_scale", 
                            output.name.intern = "beta[5] for lp_scale", 
                            initial = 1,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(1, 10),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ), 
                        theta6 = list(
                            hyperid = 103006,
                            name = "beta6",
                            short.name = "b6",
                            output.name = "beta[6] for lp_scale", 
                            output.name.intern = "beta[6] for lp_scale", 
                            initial = 1,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(1, 10),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ), 
                        theta7 = list(
                            hyperid = 103007,
                            name = "beta7",
                            short.name = "b7",
                            output.name = "beta[7] for lp_scale", 
                            output.name.intern = "beta[7] for lp_scale", 
                            initial = 1,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(1, 10),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ), 
                        theta8 = list(
                            hyperid = 103008,
                            name = "beta8",
                            short.name = "b8",
                            output.name = "beta[8] for lp_scale", 
                            output.name.intern = "beta[8] for lp_scale", 
                            initial = 1,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(1, 10),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ), 
                        theta9 = list(
                            hyperid = 103009,
                            name = "beta9",
                            short.name = "b9",
                            output.name = "beta[9] for lp_scale", 
                            output.name.intern = "beta[9] for lp_scale", 
                            initial = 1,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(1, 10),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ), 
                        theta10 = list(
                            hyperid = 103010,
                            name = "beta10",
                            short.name = "b10",
                            output.name = "beta[10] for lp_scale", 
                            output.name.intern = "beta[10] for lp_scale", 
                            initial = 1,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(1, 10),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ), 
                        theta11 = list(
                            hyperid = 103011,
                            name = "beta11",
                            short.name = "b11",
                            output.name = "beta[11] for lp_scale", 
                            output.name.intern = "beta[11] for lp_scale", 
                            initial = 1,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(1, 10),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ), 
                        theta12 = list(
                            hyperid = 103012,
                            name = "beta12",
                            short.name = "b12",
                            output.name = "beta[12] for lp_scale", 
                            output.name.intern = "beta[12] for lp_scale", 
                            initial = 1,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(1, 10),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ), 
                        theta13 = list(
                            hyperid = 103013,
                            name = "beta13",
                            short.name = "b13",
                            output.name = "beta[13] for lp_scale", 
                            output.name.intern = "beta[13] for lp_scale", 
                            initial = 1,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(1, 10),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ), 
                        theta14 = list(
                            hyperid = 103014,
                            name = "beta14",
                            short.name = "b14",
                            output.name = "beta[14] for lp_scale", 
                            output.name.intern = "beta[14] for lp_scale", 
                            initial = 1,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(1, 10),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ), 
                        theta15 = list(
                            hyperid = 103015,
                            name = "beta15",
                            short.name = "b15",
                            output.name = "beta[15] for lp_scale", 
                            output.name.intern = "beta[15] for lp_scale", 
                            initial = 1,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(1, 10),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ), 
                        theta16 = list(
                            hyperid = 103016,
                            name = "beta16",
                            short.name = "b16",
                            output.name = "beta[16] for lp_scale", 
                            output.name.intern = "beta[16] for lp_scale", 
                            initial = 1,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(1, 10),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ), 
                        theta17 = list(
                            hyperid = 103017,
                            name = "beta17",
                            short.name = "b17",
                            output.name = "beta[17] for lp_scale", 
                            output.name.intern = "beta[17] for lp_scale", 
                            initial = 1,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(1, 10),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ), 
                        theta18 = list(
                            hyperid = 103018,
                            name = "beta18",
                            short.name = "b18",
                            output.name = "beta[18] for lp_scale", 
                            output.name.intern = "beta[18] for lp_scale", 
                            initial = 1,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(1, 10),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ), 
                        theta19 = list(
                            hyperid = 103019,
                            name = "beta19",
                            short.name = "b19",
                            output.name = "beta[19] for lp_scale", 
                            output.name.intern = "beta[19] for lp_scale", 
                            initial = 1,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(1, 10),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ), 
                        theta20 = list(
                            hyperid = 103020,
                            name = "beta20",
                            short.name = "b20",
                            output.name = "beta[20] for lp_scale", 
                            output.name.intern = "beta[20] for lp_scale", 
                            initial = 1,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(1, 10),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ), 
                        theta21 = list(
                            hyperid = 103021,
                            name = "beta21",
                            short.name = "b21",
                             output.name = "beta[21] for lp_scale", 
                            output.name.intern = "beta[21] for lp_scale", 
                           initial = 1,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(1, 10),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ), 
                        theta22 = list(
                            hyperid = 103022,
                            name = "beta22",
                            short.name = "b22",
                            output.name = "beta[22] for lp_scale", 
                            output.name.intern = "beta[22] for lp_scale", 
                            initial = 1,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(1, 10),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ), 
                        theta23 = list(
                            hyperid = 103023,
                            name = "beta23",
                            short.name = "b23",
                            output.name = "beta[23] for lp_scale", 
                            output.name.intern = "beta[23] for lp_scale", 
                            initial = 1,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(1, 10),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ), 
                        theta24 = list(
                            hyperid = 103024,
                            name = "beta24",
                            short.name = "b24",
                            output.name = "beta[24] for lp_scale", 
                            output.name.intern = "beta[24] for lp_scale", 
                            initial = 1,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(1, 10),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ), 
                        theta25 = list(
                            hyperid = 103025,
                            name = "beta25",
                            short.name = "b25",
                            output.name = "beta[25] for lp_scale", 
                            output.name.intern = "beta[25] for lp_scale", 
                            initial = 1,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(1, 10),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ), 
                        theta26 = list(
                            hyperid = 103026,
                            name = "beta26",
                            short.name = "b26",
                             output.name = "beta[26] for lp_scale", 
                            output.name.intern = "beta[26] for lp_scale", 
                           initial = 1,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(1, 10),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ), 
                        theta27 = list(
                            hyperid = 103027,
                            name = "beta27",
                            short.name = "b27",
                            output.name = "beta[27] for lp_scale", 
                            output.name.intern = "beta[27] for lp_scale", 
                            initial = 1,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(1, 10),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ), 
                        theta28 = list(
                            hyperid = 103028,
                            name = "beta28",
                            short.name = "b28",
                            output.name = "beta[28] for lp_scale", 
                            output.name.intern = "beta[28] for lp_scale", 
                            initial = 1,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(1, 10),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ), 
                        theta29 = list(
                            hyperid = 103029,
                            name = "beta29",
                            short.name = "b29",
                            output.name = "beta[29] for lp_scale", 
                            output.name.intern = "beta[29] for lp_scale", 
                            initial = 1,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(1, 10),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ), 
                        theta30 = list(
                            hyperid = 103030,
                            name = "beta30",
                            short.name = "b30",
                            output.name = "beta[30] for lp_scale", 
                            output.name.intern = "beta[30] for lp_scale", 
                            initial = 1,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(1, 10),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ), 
                        theta31 = list(
                            hyperid = 103031,
                            name = "beta31",
                            short.name = "b31",
                            output.name = "beta[31] for lp_scale", 
                            output.name.intern = "beta[31] for lp_scale", 
                            initial = 1,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(1, 10),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ), 
                        theta32 = list(
                            hyperid = 103032,
                            name = "beta32",
                            short.name = "b32",
                            output.name = "beta[32] for lp_scale", 
                            output.name.intern = "beta[32] for lp_scale", 
                            initial = 1,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(1, 10),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ), 
                        theta33 = list(
                            hyperid = 103033,
                            name = "beta33",
                            short.name = "b33",
                            output.name = "beta[33] for lp_scale", 
                            output.name.intern = "beta[33] for lp_scale", 
                            initial = 1,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(1, 10),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ), 
                        theta34 = list(
                            hyperid = 103034,
                            name = "beta34",
                            short.name = "b34",
                            output.name = "beta[34] for lp_scale", 
                            output.name.intern = "beta[34] for lp_scale", 
                            initial = 1,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(1, 10),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ), 
                        theta35 = list(
                            hyperid = 103035,
                            name = "beta35",
                            short.name = "b35",
                            output.name = "beta[35] for lp_scale", 
                            output.name.intern = "beta[35] for lp_scale", 
                            initial = 1,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(1, 10),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ), 
                        theta36 = list(
                            hyperid = 103036,
                            name = "beta36",
                            short.name = "b36",
                            output.name = "beta[36] for lp_scale", 
                            output.name.intern = "beta[36] for lp_scale", 
                            initial = 1,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(1, 10),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ), 
                        theta37 = list(
                            hyperid = 103037,
                            name = "beta37",
                            short.name = "b37",
                            output.name = "beta[37] for lp_scale", 
                            output.name.intern = "beta[37] for lp_scale", 
                            initial = 1,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(1, 10),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ), 
                        theta38 = list(
                            hyperid = 103038,
                            name = "beta38",
                            short.name = "b38",
                            output.name = "beta[38] for lp_scale", 
                            output.name.intern = "beta[38] for lp_scale", 
                            initial = 1,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(1, 10),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ), 
                        theta39 = list(
                            hyperid = 103039,
                            name = "beta39",
                            short.name = "b39",
                            output.name = "beta[39] for lp_scale", 
                            output.name.intern = "beta[39] for lp_scale", 
                            initial = 1,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(1, 10),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ), 
                        theta40 = list(
                            hyperid = 103040,
                            name = "beta40",
                            short.name = "b40",
                            output.name = "beta[40] for lp_scale", 
                            output.name.intern = "beta[40] for lp_scale", 
                            initial = 1,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(1, 10),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ), 
                        theta41 = list(
                            hyperid = 103041,
                            name = "beta41",
                            short.name = "b41",
                            output.name = "beta[41] for lp_scale", 
                            output.name.intern = "beta[41] for lp_scale", 
                            initial = 1,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(1, 10),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ), 
                        theta42 = list(
                            hyperid = 103042,
                            name = "beta42",
                            short.name = "b42",
                            output.name = "beta[42] for lp_scale", 
                            output.name.intern = "beta[42] for lp_scale", 
                            initial = 1,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(1, 10),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ), 
                        theta43 = list(
                            hyperid = 103043,
                            name = "beta43",
                            short.name = "b43",
                            output.name = "beta[43] for lp_scale", 
                            output.name.intern = "beta[43] for lp_scale", 
                            initial = 1,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(1, 10),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ), 
                        theta44 = list(
                            hyperid = 103044,
                            name = "beta44",
                            short.name = "b44",
                            output.name = "beta[44] for lp_scale", 
                            output.name.intern = "beta[44] for lp_scale", 
                            initial = 1,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(1, 10),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ), 
                        theta45 = list(
                            hyperid = 103045,
                            name = "beta45",
                            short.name = "b45",
                            output.name = "beta[45] for lp_scale", 
                            output.name.intern = "beta[45] for lp_scale", 
                            initial = 1,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(1, 10),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ), 
                        theta46 = list(
                            hyperid = 103046,
                            name = "beta46",
                            short.name = "b46",
                            output.name = "beta[46] for lp_scale", 
                            output.name.intern = "beta[46] for lp_scale", 
                            initial = 1,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(1, 10),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ), 
                        theta47 = list(
                            hyperid = 103047,
                            name = "beta47",
                            short.name = "b47",
                            output.name = "beta[47] for lp_scale", 
                            output.name.intern = "beta[47] for lp_scale", 
                            initial = 1,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(1, 10),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ), 
                        theta48 = list(
                            hyperid = 103048,
                            name = "beta48",
                            short.name = "b48",
                            output.name = "beta[48] for lp_scale", 
                            output.name.intern = "beta[48] for lp_scale", 
                            initial = 1,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(1, 10),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ), 
                        theta49 = list(
                            hyperid = 103049,
                            name = "beta49",
                            short.name = "b49",
                            output.name = "beta[49] for lp_scale", 
                            output.name.intern = "beta[49] for lp_scale", 
                            initial = 1,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(1, 10),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ), 
                        theta50 = list(
                            hyperid = 103050,
                            name = "beta50",
                            short.name = "b50",
                            output.name = "beta[50] for lp_scale", 
                            output.name.intern = "beta[50] for lp_scale", 
                            initial = 1,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(1, 10),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ), 
                        theta51 = list(
                            hyperid = 103051,
                            name = "beta51",
                            short.name = "b51",
                            output.name = "beta[51] for lp_scale", 
                            output.name.intern = "beta[51] for lp_scale", 
                            initial = 1,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(1, 10),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ), 
                        theta52 = list(
                            hyperid = 103052,
                            name = "beta52",
                            short.name = "b52",
                            output.name = "beta[52] for lp_scale", 
                            output.name.intern = "beta[52] for lp_scale", 
                            initial = 1,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(1, 10),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ),
                        theta53 = list(
                            hyperid = 103053,
                            name = "beta53",
                            short.name = "b53",
                            output.name = "beta[53] for lp_scale", 
                            output.name.intern = "beta[53] for lp_scale", 
                            initial = 1,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(1, 10),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ), 
                        theta54 = list(
                            hyperid = 103054,
                            name = "beta54",
                            short.name = "b54",
                            output.name = "beta[54] for lp_scale", 
                            output.name.intern = "beta[54] for lp_scale", 
                            initial = 1,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(1, 10),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ), 
                        theta55 = list(
                            hyperid = 103055,
                            name = "beta55",
                            short.name = "b55",
                            output.name = "beta[55] for lp_scale", 
                            output.name.intern = "beta[55] for lp_scale", 
                            initial = 1,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(1, 10),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ), 
                        theta56 = list(
                            hyperid = 103056,
                            name = "beta56",
                            short.name = "b56",
                            output.name = "beta[56] for lp_scale", 
                            output.name.intern = "beta[56] for lp_scale", 
                            initial = 1,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(1, 10),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ), 
                        theta57 = list(
                            hyperid = 103057,
                            name = "beta57",
                            short.name = "b57",
                            output.name = "beta[57] for lp_scale", 
                            output.name.intern = "beta[57] for lp_scale", 
                            initial = 1,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(1, 10),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ), 
                        theta58 = list(
                            hyperid = 103058,
                            name = "beta58",
                            short.name = "b58",
                            output.name = "beta[58] for lp_scale", 
                            output.name.intern = "beta[58] for lp_scale", 
                            initial = 1,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(1, 10),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ), 
                        theta59 = list(
                            hyperid = 103059,
                            name = "beta59",
                            short.name = "b59",
                            output.name = "beta[59] for lp_scale", 
                            output.name.intern = "beta[59] for lp_scale", 
                            initial = 1,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(1, 10),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ), 
                        theta60 = list(
                            hyperid = 103060,
                            name = "beta60",
                            short.name = "b60",
                            output.name = "beta[60] for lp_scale", 
                            output.name.intern = "beta[60] for lp_scale", 
                            initial = 1,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(1, 10),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ), 
                        theta61 = list(
                            hyperid = 103061,
                            name = "beta61",
                            short.name = "b61",
                            output.name = "beta[61] for lp_scale", 
                            output.name.intern = "beta[61] for lp_scale", 
                            initial = 1,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(1, 10),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ), 
                        theta62 = list(
                            hyperid = 103062,
                            name = "beta62",
                            short.name = "b62",
                            output.name = "beta[62] for lp_scale", 
                            output.name.intern = "beta[62] for lp_scale", 
                            initial = 1,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(1, 10),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ), 
                        theta63 = list(
                            hyperid = 103063,
                            name = "beta63",
                            short.name = "b63",
                            output.name = "beta[63] for lp_scale", 
                            output.name.intern = "beta[63] for lp_scale", 
                            initial = 1,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(1, 10),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ), 
                        theta64 = list(
                            hyperid = 103064,
                            name = "beta64",
                            short.name = "b64",
                            output.name = "beta[64] for lp_scale", 
                            output.name.intern = "beta[64] for lp_scale", 
                            initial = 1,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(1, 10),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ), 
                        theta65 = list(
                            hyperid = 103065,
                            name = "beta65",
                            short.name = "b65",
                            output.name = "beta[65] for lp_scale", 
                            output.name.intern = "beta[65] for lp_scale", 
                            initial = 1,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(1, 10),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ), 
                        theta66 = list(
                            hyperid = 103066,
                            name = "beta66",
                            short.name = "b66",
                            output.name = "beta[66] for lp_scale", 
                            output.name.intern = "beta[66] for lp_scale", 
                            initial = 1,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(1, 10),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ), 
                        theta67 = list(
                            hyperid = 103067,
                            name = "beta67",
                            short.name = "b67",
                            output.name = "beta[67] for lp_scale", 
                            output.name.intern = "beta[67] for lp_scale", 
                            initial = 1,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(1, 10),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ), 
                        theta68 = list(
                            hyperid = 103068,
                            name = "beta68",
                            short.name = "b68",
                            output.name = "beta[68] for lp_scale", 
                            output.name.intern = "beta[68] for lp_scale", 
                            initial = 1,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(1, 10),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ), 
                        theta69 = list(
                            hyperid = 103069,
                            name = "beta69",
                            short.name = "b69",
                            output.name = "beta[69] for lp_scale", 
                            output.name.intern = "beta[69] for lp_scale", 
                            initial = 1,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(1, 10),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ), 
                        theta70 = list(
                            hyperid = 103070,
                            name = "beta70",
                            short.name = "b70",
                            output.name = "beta[70] for lp_scale", 
                            output.name.intern = "beta[70] for lp_scale", 
                            initial = 1,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(1, 10),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ), 
                        theta71 = list(
                            hyperid = 103071,
                            name = "beta71",
                            short.name = "b71",
                            output.name = "beta[71] for lp_scale", 
                            output.name.intern = "beta[71] for lp_scale", 
                            initial = 1,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(1, 10),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ), 
                        theta72 = list(
                            hyperid = 103072,
                            name = "beta72",
                            short.name = "b72",
                            output.name = "beta[72] for lp_scale", 
                            output.name.intern = "beta[72] for lp_scale", 
                            initial = 1,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(1, 10),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ), 
                        theta73 = list(
                            hyperid = 103073,
                            name = "beta73",
                            short.name = "b73",
                            output.name = "beta[73] for lp_scale", 
                            output.name.intern = "beta[73] for lp_scale", 
                            initial = 1,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(1, 10),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ), 
                        theta74 = list(
                            hyperid = 103074,
                            name = "beta74",
                            short.name = "b74",
                            output.name = "beta[74] for lp_scale", 
                            output.name.intern = "beta[74] for lp_scale", 
                            initial = 1,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(1, 10),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ), 
                        theta75 = list(
                            hyperid = 103075,
                            name = "beta75",
                            short.name = "b75",
                            output.name = "beta[75] for lp_scale", 
                            output.name.intern = "beta[75] for lp_scale", 
                            initial = 1,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(1, 10),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ), 
                        theta76 = list(
                            hyperid = 103076,
                            name = "beta76",
                            short.name = "b76",
                            output.name = "beta[76] for lp_scale", 
                            output.name.intern = "beta[76] for lp_scale", 
                            initial = 1,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(1, 10),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ), 
                        theta77 = list(
                            hyperid = 103077,
                            name = "beta77",
                            short.name = "b77",
                            output.name = "beta[77] for lp_scale", 
                            output.name.intern = "beta[77] for lp_scale", 
                            initial = 1,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(1, 10),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ), 
                        theta78 = list(
                            hyperid = 103078,
                            name = "beta78",
                            short.name = "b78",
                            output.name = "beta[78] for lp_scale", 
                            output.name.intern = "beta[78] for lp_scale", 
                            initial = 1,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(1, 10),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ), 
                        theta79 = list(
                            hyperid = 103079,
                            name = "beta79",
                            short.name = "b79",
                            output.name = "beta[79] for lp_scale", 
                            output.name.intern = "beta[79] for lp_scale", 
                            initial = 1,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(1, 10),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ), 
                        theta80 = list(
                            hyperid = 103080,
                            name = "beta80",
                            short.name = "b80",
                            output.name = "beta[80] for lp_scale", 
                            output.name.intern = "beta[80] for lp_scale", 
                            initial = 1,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(1, 10),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ), 
                        theta81 = list(
                            hyperid = 103081,
                            name = "beta81",
                            short.name = "b81",
                            output.name = "beta[81] for lp_scale", 
                            output.name.intern = "beta[81] for lp_scale", 
                            initial = 1,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(1, 10),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ), 
                        theta82 = list(
                            hyperid = 103082,
                            name = "beta82",
                            short.name = "b82",
                            output.name = "beta[82] for lp_scale", 
                            output.name.intern = "beta[82] for lp_scale", 
                            initial = 1,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(1, 10),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ), 
                        theta83 = list(
                            hyperid = 103083,
                            name = "beta83",
                            short.name = "b83",
                            output.name = "beta[83] for lp_scale", 
                            output.name.intern = "beta[83] for lp_scale", 
                            initial = 1,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(1, 10),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ), 
                        theta84 = list(
                            hyperid = 103084,
                            name = "beta84",
                            short.name = "b84",
                            output.name = "beta[84] for lp_scale", 
                            output.name.intern = "beta[84] for lp_scale", 
                            initial = 1,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(1, 10),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ), 
                        theta85 = list(
                            hyperid = 103085,
                            name = "beta85",
                            short.name = "b85",
                            output.name = "beta[85] for lp_scale", 
                            output.name.intern = "beta[85] for lp_scale", 
                            initial = 1,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(1, 10),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ), 
                        theta86 = list(
                            hyperid = 103086,
                            name = "beta86",
                            short.name = "b86",
                            output.name = "beta[86] for lp_scale", 
                            output.name.intern = "beta[86] for lp_scale", 
                            initial = 1,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(1, 10),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ), 
                        theta87 = list(
                            hyperid = 103087,
                            name = "beta87",
                            short.name = "b87",
                            output.name = "beta[87] for lp_scale", 
                            output.name.intern = "beta[87] for lp_scale", 
                            initial = 1,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(1, 10),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ), 
                        theta88 = list(
                            hyperid = 103088,
                            name = "beta88",
                            short.name = "b88",
                            output.name = "beta[88] for lp_scale", 
                            output.name.intern = "beta[88] for lp_scale", 
                            initial = 1,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(1, 10),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ), 
                        theta89 = list(
                            hyperid = 103089,
                            name = "beta89",
                            short.name = "b89",
                            output.name = "beta[89] for lp_scale", 
                            output.name.intern = "beta[89] for lp_scale", 
                            initial = 1,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(1, 10),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ), 
                        theta90 = list(
                            hyperid = 103090,
                            name = "beta90",
                            short.name = "b90",
                            output.name = "beta[90] for lp_scale", 
                            output.name.intern = "beta[90] for lp_scale", 
                            initial = 1,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(1, 10),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ), 
                        theta91 = list(
                            hyperid = 103091,
                            name = "beta91",
                            short.name = "b91",
                            output.name = "beta[91] for lp_scale", 
                            output.name.intern = "beta[91] for lp_scale", 
                            initial = 1,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(1, 10),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ), 
                        theta92 = list(
                            hyperid = 103092,
                            name = "beta92",
                            short.name = "b92",
                            output.name = "beta[92] for lp_scale", 
                            output.name.intern = "beta[92] for lp_scale", 
                            initial = 1,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(1, 10),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ), 
                        theta93 = list(
                            hyperid = 103093,
                            name = "beta93",
                            short.name = "b93",
                            output.name = "beta[93] for lp_scale", 
                            output.name.intern = "beta[93] for lp_scale", 
                            initial = 1,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(1, 10),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ), 
                        theta94 = list(
                            hyperid = 103094,
                            name = "beta94",
                            short.name = "b94",
                            output.name = "beta[94] for lp_scale", 
                            output.name.intern = "beta[94] for lp_scale", 
                            initial = 1,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(1, 10),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ), 
                        theta95 = list(
                            hyperid = 103095,
                            name = "beta95",
                            short.name = "b95",
                            output.name = "beta[95] for lp_scale", 
                            output.name.intern = "beta[95] for lp_scale", 
                            initial = 1,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(1, 10),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ), 
                        theta96 = list(
                            hyperid = 103096,
                            name = "beta96",
                            short.name = "b96",
                            output.name = "beta[96] for lp_scale", 
                            output.name.intern = "beta[96] for lp_scale", 
                            initial = 1,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(1, 10),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ), 
                        theta97 = list(
                            hyperid = 103097,
                            name = "beta97",
                            short.name = "b97",
                            output.name = "beta[97] for lp_scale", 
                            output.name.intern = "beta[97] for lp_scale", 
                            initial = 1,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(1, 10),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ), 
                        theta98 = list(
                            hyperid = 103098,
                            name = "beta98",
                            short.name = "b98",
                            output.name = "beta[98] for lp_scale", 
                            output.name.intern = "beta[98] for lp_scale", 
                            initial = 1,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(1, 10),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ), 
                        theta99 = list(
                            hyperid = 103099,
                            name = "beta99",
                            short.name = "b99",
                            output.name = "beta[99] for lp_scale", 
                            output.name.intern = "beta[99] for lp_scale", 
                            initial = 1,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(1, 10),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        ), 
                        theta100 = list(
                            hyperid = 103100,
                            name = "beta100",
                            short.name = "b100",
                            output.name = "beta[100] for lp_scale", 
                            output.name.intern = "beta[100] for lp_scale", 
                            initial = 1,
                            fixed = FALSE,
                            prior = "normal",
                            param = c(1, 10),
                            to.theta = function(x) x,
                            from.theta = function(x) x
                        )
                    ), 
                    pdf = "lp.scale"
                )
            )
    )
}

## Documentation is generated by inla.models.generate.roxygen()
## and (for r-inla.org/doc/) inla.models.generate.tex()
`inla.models` <- function() {
    ## this is not very clean solution, but for the moment is ok. the
    ## inla.models() function takes just to much time!!!

    envir <- inla.get.inlaEnv()

    if (exists("inla.models", envir = envir) &&
        !is.null(get("inla.models", envir = envir)) && 
        exists("rinla.version", envir = envir) &&
        get("rinla.version", envir = envir) == inla.version("version")) {
        return(get("inla.models", envir = envir))
    } else {
        ## have to split it, as option keep.source has an upper limit...
        models <- c(
            inla.models.section.latent(),
            inla.models.section.group(),
            inla.models.section.scopy(),
            inla.models.section.mix(),
            inla.models.section.link(),
            inla.models.section.predictor(),
            inla.models.section.hazard(),
            inla.models.section.likelihood(),
            inla.models.section.prior(),
            inla.models.section.wrapper(), 
            inla.models.section.lp.scale()
        )
        ## set "read.only" attribute for the `hyper' at those elements
        ## that cannot be changed.
        for (section in names(models)) {
            for (model in names(models[[section]])) {
                a <- models[[section]][[model]]
                if (!is.null(a$hyper) && length(a$hyper) > 0) {
                    for (theta in names(a$hyper)) {
                        for (elm in names(a$hyper[[theta]])) {
                            val <- FALSE
                            if (elm == "prior") {
                                if (a$hyper[[theta]][[elm]] == "none" ||
                                    a$hyper[[theta]][[elm]] == "wishart1d" ||
                                    a$hyper[[theta]][[elm]] == "wishart2d" ||
                                    a$hyper[[theta]][[elm]] == "wishart3d" ||
                                    a$hyper[[theta]][[elm]] == "wishart4d" ||
                                    a$hyper[[theta]][[elm]] == "wishart5d") {
                                    val <- TRUE
                                }
                            }
                            if (elm == "to.theta" || elm == "from.theta") {
                                  val <- TRUE
                              }
                            ## append the new attribute, so the code
                            ## depends on if the object has an attribute
                            ## from before or not.
                            if (is.null(attributes(models[[section]][[model]]$hyper[[theta]][[elm]]))) {
                                attributes(models[[section]][[model]]$hyper[[theta]][[elm]]) <-
                                    list("inla.read.only" = val)
                            } else {
                                attributes(models[[section]][[model]]$hyper[[theta]][[elm]]) <-
                                    c(
                                        "inla.read.only" = val,
                                        attributes(models[[section]][[model]]$hyper[[theta]][[elm]])
                                    )
                            }
                        }
                    }
                }
            }
        }

        assign("inla.models", models, envir = envir)
        assign("rinla.version", inla.version("version"), envir = envir)

        return(models)
    }

    stop("This should not happen")
}

`inla.is.model` <- function(model, section = NULL,
                            stop.on.error = TRUE, ignore.case = FALSE) {
    mm <- inla.models()
    if (is.null(section)) {
        stop("No section given; please fix...")
    }
    section <- match.arg(section, names(mm))
    models <- names((mm[names(mm) == section])[[1]])

    if (is.character(model) && length(model) > 0) {
        if (ignore.case) {
            m <- tolower(model)
            ms <- tolower(models)
        } else {
            m <- model
            ms <- models
        }

        ret <- c()
        for (i in 1L:length(m)) {
            if (is.element(m[i], ms)) {
                ret[i] <- TRUE
            } else {
                if (stop.on.error) {
                    print("Valid models are:")
                    print(models)
                    stop(paste(c(
                        "\n\tUnknown name [", model[i], "]\n", "\tValid choices are: ",
                        models
                    ), sep = " ", collapse = " "))
                } else {
                    ret[i] <- FALSE
                }
            }
        }
    } else {
        if (stop.on.error) {
              stop(paste("\n\tModel [", model, "] is not of type character()\n", sep = ""))
          }
    }
    return(ret)
}

`inla.model.properties` <- function(
                                    model,
                                    section = NULL,
                                    stop.on.error = TRUE,
                                    ignore.case = FALSE) {
    if (is.null(section)) {
        stop("No 'section' given; please fix...")
    }
    if (is.null(model)) {
        stop("No 'model' given; please fix...")
    }
    mm <- inla.models()
    section <- match.arg(section, names(mm))
    m <- inla.model.properties.generic(
        inla.trim.family(model),
        mm[names(mm) == section][[1]],
        stop.on.error, ignore.case,
        ## would like to know the section for a possible warning/error if the status says so.
        section = section
    )

    if (is.null(m)) {
        return(NULL)
    }

    return(m)
}

`inla.model.properties.generic` <- function(model, models, stop.on.error = TRUE, ignore.case = TRUE, section = "UKNOWN") {
    ## argument 'section' is for 'status' only...

    m <- ifelse(ignore.case, tolower(model), model)
    if (ignore.case) {
        ms <- tolower(names(models))
    } else {
        ms <- names(models)
    }
    ms <- inla.trim.family(ms)
    k <- grep(paste("^", m, "$", sep = ""), ms)
    if (length(k) == 0L) {
        if (stop.on.error) {
            stop(paste("Name '", model, "' is unknown. Valid choices are: ", inla.paste(ms), ".", sep = ""))
        }
        return(NULL)
    } else {
        ## if 'status' is set, then issue a warning/error depending on
        ## 'status'. do this the first time only if status is
        ## 'experimental'.
        status <- models[[k]]$status
        if (is.null(status)) {
            ## do nothing; all ok.
        } else {
            status.core <- strsplit(status, ":")[[1]][1]
            stopifnot(any(inla.strcasecmp(status.core, c("experimental", "disabled", "changed"))))
            envir <- inla.get.inlaEnv()
            var <- paste("processed.status.for.model.", model, ".in.section.", section, sep = "")
            if (inla.strcasecmp(status.core, "experimental")) {
                if (!(exists(var, envir = envir) && get(var, envir = envir))) {
                    assign(var, TRUE, envir = envir)
                    msg <- paste0(
                        "Model '", model, "' in section '", section, "' is marked as '", status,
                        "'; changes may appear at any time.",
                        "\n  ",
                        "Further warnings are disabled."
                    )
                    warning(msg)
                } else {
                    ## the warning is already given; do nothing
                }
            } else if (inla.strcasecmp(status.core, "disabled")) {
                assign(var, TRUE, envir = envir)
                msg <- paste("Model '", model, "' in section '", section, "' is marked as '",
                    status, ".\n",
                    "Usage is either not recommended and/or unsupported.\n",
                    "Email <help@r-inla.org> for requests.\n",
                    "\n",
                    sep = ""
                )
                var <- paste("enable.model.", section, ".", model, sep = "")
                if (!(exists(var, envir = envir) && get(var, envir = envir))) {
                    msg <- paste0(c(msg, paste(
                        "  You can enable this model setting variable '", var,
                        "'\n  to 'TRUE' in environment inla.get.inlaEnv().\n",
                        "  If you chose to do so, you are on your own."
                    )))
                    stop(msg)
                }
            } else if (inla.strcasecmp(status.core, "changed")) {
                assign(var, TRUE, envir = envir)
                msg <- paste0(
                    "Model '", model, "' in section '", section, "' is marked as '",
                    status, ".\n",
                    "  There have been a change in the model definition, which is not backward compatible.\n",
                    "  Please refer to the documentation before proceeeding."
                )
                var <- paste("enable.model.", section, ".", model, sep = "")
                if (!(exists(var, envir = envir) && get(var, envir = envir))) {
                    msg <- paste0(c(msg, paste0(
                        "\n  You can bypass this check setting variable '", var,
                        "'\n  to 'TRUE' in environment inla.get.inlaEnv().\n"
                    )))
                    stop(msg)
                }
            }
        }

        return(models[[k]])
    }
}

`inla.model.validate.link.function` <- function(model, link) {
    valid.links <- inla.model.properties(model, "likelihood")$link

    stopifnot(!is.null(valid.links))
    stopifnot(length(valid.links) >= 2L)
    stopifnot(valid.links[1L] == "default")

    link <- tolower(link)
    if (is.element(link, valid.links)) {
        ## this is the convention: the default link is the second
        ## entry in the list. the first entry is always "default"
        if (link == "default") {
            link <- valid.links[2L]
        }
    } else {
        stop(inla.paste(c(
            "Link function `", link, "' is not valid or yet implemented.",
            "\n",
            "Valid ones are: ", inla.paste(valid.links), "\n"
        ), sep = ""))
    }

    return(link)
}

`inla.model.validate.link.simple.function` <- function(model, link) {
    valid.links <- inla.model.properties(model, "likelihood")$link.simple

    if (is.null(valid.links)) {
        return (NULL )
    }

    stopifnot(!is.null(valid.links))
    stopifnot(length(valid.links) >= 2L)
    stopifnot(valid.links[1L] == "default")

    link <- tolower(link)
    if (is.element(link, valid.links)) {
        ## this is the convention: the default link is the second
        ## entry in the list. the first entry is always "default"
        if (link == "default") {
            link <- valid.links[2L]
        }
    } else {
        stop(inla.paste(c(
            "Link function `", link, "' is not valid or yet implemented.",
            "\n",
            "Valid ones are: ", inla.paste(valid.links), "\n"
        ), sep = ""))
    }

    return(link)
}
