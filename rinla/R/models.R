## models and its hyperparameters are defined here

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

              me = list(
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
                                      name = "prec.obs",
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
                      pdf = "me"
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
                      pdf = "besagproper"
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
                                      param = c(0, 0.15),  ## same as for AR1
                                      to.theta = function(x) log((1+x)/(1-x)),
                                      from.theta = function(x) 2*exp(x)/(1+exp(x))-1
                                      ),
                              theta3 = list(
                                      name = "pacf2",
                                      short.name = "pacf2",
                                      initial = 0,
                                      fixed = TRUE,
                                      prior = "none",
                                      param = numeric(0),
                                      to.theta = function(x) log((1+x)/(1-x)),
                                      from.theta = function(x) 2*exp(x)/(1+exp(x))-1
                                      ),
                              theta4 = list(
                                      name = "pacf3",
                                      short.name = "pacf3",
                                      initial = 0,
                                      fixed = TRUE,
                                      prior = "none",
                                      param = numeric(0),
                                      to.theta = function(x) log((1+x)/(1-x)),
                                      from.theta = function(x) 2*exp(x)/(1+exp(x))-1
                                      ),
                              theta5 = list(
                                      name = "pacf4",
                                      short.name = "pacf4",
                                      initial = 0,
                                      fixed = TRUE,
                                      prior = "none",
                                      param = numeric(0),
                                      to.theta = function(x) log((1+x)/(1-x)),
                                      from.theta = function(x) 2*exp(x)/(1+exp(x))-1
                                      ),
                              theta6 = list(
                                      name = "pacf5",
                                      short.name = "pacf5",
                                      initial = 0,
                                      fixed = TRUE,
                                      prior = "none",
                                      param = numeric(0),
                                      to.theta = function(x) log((1+x)/(1-x)),
                                      from.theta = function(x) 2*exp(x)/(1+exp(x))-1
                                      ),
                              theta7 = list(
                                      name = "pacf6",
                                      short.name = "pacf6",
                                      initial = 0,
                                      fixed = TRUE,
                                      prior = "none",
                                      param = numeric(0),
                                      to.theta = function(x) log((1+x)/(1-x)),
                                      from.theta = function(x) 2*exp(x)/(1+exp(x))-1
                                      ),
                              theta8 = list(
                                      name = "pacf7",
                                      short.name = "pacf7",
                                      initial = 0,
                                      fixed = TRUE,
                                      prior = "none",
                                      param = numeric(0),
                                      to.theta = function(x) log((1+x)/(1-x)),
                                      from.theta = function(x) 2*exp(x)/(1+exp(x))-1
                                      ),
                              theta9 = list(
                                      name = "pacf8",
                                      short.name = "pacf8",
                                      initial = 0,
                                      fixed = TRUE,
                                      prior = "none",
                                      param = numeric(0),
                                      to.theta = function(x) log((1+x)/(1-x)),
                                      from.theta = function(x) 2*exp(x)/(1+exp(x))-1
                                      ),
                              theta10 = list(
                                      name = "pacf9",
                                      short.name = "pacf9",
                                      initial = 0,
                                      fixed = TRUE,
                                      prior = "none",
                                      param = numeric(0),
                                      to.theta = function(x) log((1+x)/(1-x)),
                                      from.theta = function(x) 2*exp(x)/(1+exp(x))-1
                                      ),
                              theta11 = list(
                                      name = "pacf10",
                                      short.name = "pacf10",
                                      initial = 0,
                                      fixed = TRUE,
                                      prior = "none",
                                      param = numeric(0),
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
                                      param = c(-2, 1),
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
                                      param = c(0, 0.001),
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
                                      from.theta = function(x) exp(x)                                                                                                                          )
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
                      aug.factor = 1L,
                      aug.constr = NULL,
                      n.div.by = NULL,
                      n.required = FALSE,
                      set.default.values = FALSE,
                      pdf = NA
                      )
              ##
              )
         )
}

`inla.models.section.group` = function()
{
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
                                      to.theta = function(x, ngroup) log((1+x*(ngroup-1))/(1-x)),
                                      from.theta = function(x, ngroup) (exp(x)-1)/(exp(x) + ngroup -1)
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
                              )
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
                              )
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
                      link = c("default", "log"),
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
                      link = c("default", "log"),
                      pdf = "gpoisson",
                      status = "experimental"
                      ),

              binomial = list(
                      hyper = list(
                              ),
                      survival = FALSE,
                      discrete = TRUE,
                      link = c("default", "logit", "probit", "cloglog", "log", "h"),
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

              binomialtest = list(
                      hyper = list(
                              theta = list(
                                      name = "psi",
                                      short.name = "psi",
                                      initial = 0,
                                      fixed = FALSE,
                                      prior = "normal",
                                      param = c(0, 100),
                                      to.theta = function(x) x,
                                      from.theta = function(x) x
                                      )
                              ),
                      survival = FALSE,
                      discrete = TRUE,
                      link = c("default", "logit", "probit", "cloglog", "log"),
                      pdf = "binomialtest",
                      status = "experimental"
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
                      link = c("default", "logit", "probit", "cloglog", "log"),
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
                      link = c("default", "log"),
                      pdf = "nbinomial"
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
                      link = c("default", "identity", "logit"),
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

              sas = list(
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
                                      name = "skewness",
                                      short.name = "skew",
                                      initial = 0,
                                      fixed = FALSE,
                                      prior = "normal",
                                      param = c(0, 25),
                                      to.theta = function(x) x,
                                      from.theta = function(x) x
                                      ),
                              theta3 = list(
                                      name = "kurtosis",
                                      short.name = "kurt",
                                      initial = 3,
                                      fixed = FALSE,
                                      prior = "normal",
                                      param = c(3, 100),
                                      to.theta = function(x) (x),
                                      from.theta = function(x) (x)
                                      )
                              ),
                      survival = FALSE,
                      discrete = FALSE,
                      link = c("default", "identity"),
                      pdf = "sas",
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
                      discrete = FALSE,
                      link = c("default", "identity"),
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
                      pdf = "add-file-here"
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
                                      initial = 3,
                                      fixed = FALSE,
                                      prior = "loggamma",
                                      param = c(1, 0.00005),
                                      to.theta = function(x) log(x),
                                      from.theta = function(x) exp(x)
                                      ),
                              theta2 = list(
                                      name = "log degrees of freedom",
                                      short.name = "dof",
                                      initial = 0,
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

              ## this is the 'no prior needed' prior
              none = list(nparameters = 0L),

              betacorrelation = list(
                      nparameters = 2L,
                      pdf = "betacorrelation"
                      ),

              logitbeta = list(
                      nparameters = 2L,
                      pdf = "logitbeta"
                      ),

              jeffreystdf = list(
                      nparameters = 0L,
                      pdf = "jeffreystdf"
                      ),

              "expression:" = list(
                  nparameters = -1L,
                  pdf = "expression"
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
        ## have to split it,  as option keep.source has an upper limit; what a shame...
        models = c(
                inla.models.section.latent(),
                inla.models.section.group(),
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
                        "'; do not use.", sep="")
                stop(msg)
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
