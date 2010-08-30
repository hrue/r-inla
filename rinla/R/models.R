
inla.models = function()
{
    return (
            list(
                 models=list(
                         ##
                         iid = list(ntheta = 1, theta = "precision", nparameters = 2, npriors = 1,
                                 constr = FALSE, nrow.ncol = FALSE, augmented = FALSE, aug.factor = 1, aug.constr = NULL,
                                 n.div.by = NULL, n.required = FALSE, set.default.values = FALSE),
                         positive = list(ntheta = 1, theta = "positive fixed-effect", nparameters = 2, npriors = 1,
                                 constr = FALSE, nrow.ncol = FALSE, augmented = FALSE, aug.factor = 1, aug.constr = NULL,
                                 n.div.by = NULL, n.required = FALSE, set.default.values = FALSE),
                         rw1 = list(ntheta = 1, theta = "precision", nparameters = 2, npriors = 1,
                                 constr = TRUE, nrow.ncol = FALSE, augmented = FALSE, aug.factor = 1, aug.constr = NULL,
                                 n.div.by = NULL, n.required = FALSE, set.default.values = FALSE),
                         rw2 = list(ntheta = 1, theta = "precision", nparameters = 2, npriors = 1,
                                 constr = TRUE, nrow.ncol = FALSE, augmented = FALSE, aug.factor = 1, aug.constr = NULL,
                                 n.div.by = NULL, n.required = FALSE, set.default.values = FALSE),
                         crw2 = list(ntheta = 1, theta = "precision", nparameters = 2, npriors = 1,
                                 constr = TRUE, nrow.ncol = FALSE, augmented = TRUE, aug.factor = 2, aug.constr = 1,
                                 n.div.by = NULL, n.required = FALSE, set.default.values = FALSE),
                         seasonal = list(ntheta = 1, theta = "precision", nparameters = 2, npriors = 1,
                                 constr = FALSE, nrow.ncol = FALSE, augmented = FALSE, aug.factor = 1, aug.constr = NULL,
                                 n.div.by = NULL, n.required = FALSE, set.default.values = FALSE),
                         besag = list(ntheta = 1, theta = "precision", nparameters = 2, npriors = 1,
                                 constr = TRUE, nrow.ncol = FALSE, augmented = FALSE, aug.factor = 1, aug.constr = NULL,
                                 n.div.by = NULL, n.required = TRUE, set.default.values = TRUE),
                         besag2 = list(ntheta = 2, theta = c("precision", "a"), nparameters = 4, npriors = 2,
                                 constr = FALSE, nrow.ncol = FALSE, augmented = FALSE, aug.factor = 1, aug.constr = c(1,2),
                                 n.div.by = 2, n.required = TRUE, set.default.values = TRUE),
                         bym = list(ntheta = 2, theta = c("precision iid", "precision spatial"), nparameters = 4, npriors = 2,
                                 constr = TRUE, nrow.ncol = FALSE, augmented = TRUE, aug.factor = 2, aug.constr = 2,
                                 n.div.by = NULL, n.required = TRUE, set.default.values = TRUE),
                         ar1 = list(ntheta = 2, theta = c("precision", "lag one correlation"), nparameters = 4, npriors = 2,
                                 constr = FALSE, nrow.ncol = FALSE, augmented = FALSE, aug.factor = 1, aug.constr = NULL,
                                 n.div.by = NULL, n.required = FALSE, set.default.values = FALSE),
                         generic = list(ntheta = 1, theta = "precision", nparameters = 2, npriors =1,
                                 constr = FALSE, nrow.ncol = FALSE, augmented = FALSE, aug.factor = 1, aug.constr = NULL,
                                 n.div.by = NULL, n.required = TRUE, set.default.values = TRUE),
                         generic0 = list(ntheta = 1, theta = "precision", nparameters = 2, npriors =1,
                                 constr = FALSE, nrow.ncol = FALSE, augmented = FALSE, aug.factor = 1, aug.constr = NULL,
                                 n.div.by = NULL, n.required = TRUE, set.default.values = TRUE),
                         generic1 = list(ntheta = 2, theta = c("precision", "beta"), nparameters = 4, npriors =2,
                                 constr = FALSE, nrow.ncol = FALSE, augmented = FALSE, aug.factor = 1, aug.constr = NULL,
                                 n.div.by = NULL, n.required = TRUE, set.default.values = TRUE),
                         generic2 = list(ntheta = 2, theta = c("precision-cmatrix", "h2"), nparameters = 4, npriors =2,
                                 constr = FALSE, nrow.ncol = FALSE, augmented = FALSE, aug.factor = 2, aug.constr = 2,
                                 n.div.by = NULL, n.required = TRUE, set.default.values = TRUE),
                         sphere = list(ntheta = 4, theta = c("theta.T", "theta.K", "theta.TK", "theta.OC"), nparameters = 8, npriors = 4,
                                 constr = FALSE, nrow.ncol = FALSE, augmented = FALSE, aug.factor = 1, aug.constr = NULL,
                                 n.div.by = NULL, n.required = FALSE, set.default.values = FALSE),
                         spde = list(ntheta = 4, theta = c("theta.T", "theta.K", "theta.TK", "theta.OC"), nparameters = 8, npriors = 4,
                                 constr = FALSE, nrow.ncol = FALSE, augmented = FALSE, aug.factor = 1, aug.constr = NULL,
                                 n.div.by = NULL, n.required = TRUE, set.default.values = TRUE),
                         iid1d = list(ntheta = 1, theta = "precision", nparameters = 2, npriors = 1,
                                 constr = FALSE, nrow.ncol = FALSE, augmented = FALSE, aug.factor = 1, aug.constr = NULL,
                                 n.div.by = NULL, n.required = FALSE, set.default.values = FALSE),
                         iid2d = list(ntheta = 3, theta = c("precision1", "precision2", "correlation"),
                                 nparameters = 4, npriors = 1, constr = FALSE, nrow.ncol = FALSE, augmented = TRUE,
                                 aug.factor = 1, aug.constr = c(1,2), n.div.by = 2, n.required = TRUE, set.default.values = TRUE),
                         iid3d = list(ntheta = 6,
                                 theta = c("precision1", "precision2", "precision3", "correlation12", "correlation13", "correlation23"),
                                 nparameters = 7, npriors = 1, constr = FALSE, nrow.ncol = FALSE, augmented = TRUE,
                                 aug.factor = 1, aug.constr = c(1,2,3), n.div.by = 3, n.required = TRUE, set.default.values = TRUE),
                         "2diid" = list(ntheta = 3, theta = c("precision1", "precision2", "correlation"),
                                 nparameters = 6, npriors = 3, constr = FALSE, nrow.ncol = FALSE, augmented = FALSE,
                                 aug.factor = 1, aug.constr = NULL, n.div.by = NULL, n.required = FALSE, set.default.values = FALSE),
                         z = list(ntheta = 1, theta = "precision", nparameters = 2, npriors = 0,
                                 constr = FALSE, nrow.ncol = FALSE, augmented = FALSE, aug.factor = 1, aug.constr = NULL,
                                 n.div.by = NULL, n.required = FALSE, set.default.values = FALSE),
                         rw2d = list(ntheta = 1, theta = "precision", nparameters = 2, npriors = 1,
                                 constr = TRUE, nrow.ncol = TRUE, augmented = FALSE, aug.factor = 1, aug.constr = NULL,
                                 n.div.by = NULL, n.required = FALSE, set.default.values = FALSE),
                         matern2d = list(ntheta = 2, theta = c("precision", "range"), nparameters = 4, npriors = 2,
                                 constr = FALSE, nrow.ncol = TRUE, augmented = FALSE, aug.factor = 1, aug.constr = NULL,
                                 n.div.by = NULL, n.required = FALSE, set.default.values = FALSE),
                         copy = list(ntheta = 1, theta = "beta", nparameters = 2, npriors = 1,
                                 constr = FALSE, nrow.ncol = FALSE, augmented = FALSE, aug.factor = 1, aug.constr = NULL,
                                 n.div.by = NULL, n.required = FALSE, set.default.values = FALSE)
                         ##
                         ),
                 lmodels = list(
                         ##
                         poisson = list(ntheta = 0, theta = NULL, nparameters = 0, npriors = 0, survival = FALSE, discrete = TRUE),
                         poissonext = list(ntheta = 2, theta = c("ThetaE1", "ThetaE2"), nparameters = 4, npriors = 2, survival = FALSE, discrete = TRUE),
                         binomial = list(ntheta = 0, theta = NULL, nparameters = 0, npriors = 0, survival = FALSE, discrete = TRUE),
                         nbinomial = list(ntheta = 1, theta = "size", nparameters = 2, npriors = 1, survival = FALSE, discrete = TRUE),
                         exponential = list(ntheta = 0, theta = NULL, nparameters = 0, npriors = 0, survival = TRUE, discrete = FALSE),
                         coxph = list(ntheta = 0, theta = NULL, nparameters = 0, npriors = 0, survival = TRUE, discrete = TRUE),
                         gaussian = list(ntheta = 1, theta = "precision", nparameters = 2, npriors = 1, survival = FALSE, discrete = FALSE),
                         normal =   list(ntheta = 1, theta = "precision", nparameters = 2, npriors = 1, survival = FALSE, discrete = FALSE),
                         ## the same
                         skewnormal = list(ntheta = 2, theta = c("inverse.scale", "skewness"), nparameters = 4, npriors = 2, survival = FALSE, discrete = FALSE),
                         sn         = list(ntheta = 2, theta = c("inverse.scale", "skewness"), nparameters = 4, npriors = 2, survival = FALSE, discrete = FALSE),
                         ##
                         laplace = list(ntheta = 1, theta = "precision", nparameters = 2, npriors = 1, survival = FALSE, discrete = FALSE),
                         weibull = list(ntheta = 1, theta = "alpha", nparameters = 2, npriors = 1, survival = TRUE, survival = TRUE, discrete = FALSE),
                         weibullcure = list(ntheta = 2, theta = c("alpha", "prob"), nparameters = 4, npriors = 2, survival = TRUE, discrete = FALSE),
                         stochvol = list(ntheta = 0, theta = NULL, nparameters = 0, npriors = 0, survival = FALSE, discrete = FALSE),
                         stochvolt = list(ntheta = 1, theta = "degrees of freedom", nparameters = 2, npriors = 1, survival = FALSE, discrete = FALSE),
                         zeroinflatedpoisson0 = list(ntheta = 1, theta = "prob", nparameters = 2, npriors = 1, survival = FALSE, discrete = TRUE),
                         zeroinflatedpoisson1 = list(ntheta = 1, theta = "prob", nparameters = 2, npriors = 1, survival = FALSE, discrete = TRUE),
                         zeroinflatedpoisson2 = list(ntheta = 1, theta = "alpha", nparameters = 2, npriors = 1, survival = FALSE, discrete = TRUE),
                         zeroinflatedbinomial0 = list(ntheta = 1, theta = "prob", nparameters = 2, npriors = 1, survival = FALSE, discrete = TRUE),
                         zeroinflatedbinomial1 = list(ntheta = 1, theta = "prob", nparameters = 2, npriors = 1, survival = FALSE, discrete = TRUE),
                         zeroinflatedbinomial2 = list(ntheta = 1, theta = "alpha", nparameters = 2, npriors = 1, survival = FALSE, discrete = TRUE),
                         zeroinflatednbinomial0 = list(ntheta = 2, theta = c("size", "prob"), nparameters = 4, npriors = 2, survival = FALSE, discrete = TRUE),
                         zeroinflatednbinomial1 = list(ntheta = 2, theta = c("size", "prob"), nparameters = 4, npriors = 2, survival = FALSE, discrete = TRUE),
                         zeroinflatednbinomial2 = list(ntheta = 2, theta = c("size", "alpha"), nparameters = 4, npriors = 2, survival = FALSE, discrete = TRUE),
                         t = list(ntheta = 2, theta =c("precision", "degrees of freedom"), nparameters = 4, npriors = 2, survival = FALSE, discrete = FALSE),
                         stochvolnig = list(ntheta = 2, theta = c("skewness", "shape"), nparameters = 4, npriors = 2, survival = FALSE, discrete = FALSE),
                         logperiodogram = list(ntheta = 0, theta = NULL, nparameters = 0, npriors = 0, survival = FALSE, discrete = FALSE)
                         ##
                         ),
                 
                 priors = list(
                         ##
                         ## these are the same
                         normal = list(),
                         gaussian = list(),
                         ##
                         ##
                         wishart = list(),
                         loggamma = list(),
                         ##
                         ## these are the same
                         minuslogsqrtruncnormal = list(),
                         logtnormal = list(),
                         logtgaussian = list(),
                         ##
                         ##
                         bymjoint=list(),
                         flat=list(),
                         logflat=list()
                         ##
                         )
                 )
            )
}
inla.is.generic = function(model, stop.on.error, models, ignore.case)
{
    if (is.character(model) && length(model) > 0) {
        if (ignore.case) {
            m = tolower(model)
            ms = tolower(names(models))
        } else {
            m = model
            ms = names(models)
        }

        ret = c()
        for(i in 1:length(m)) {
        
            if (is.element(m[i], ms)) {
                ret[i] = TRUE
            } else {
                if (stop.on.error)
                    stop(paste("\n\tUnknown model [", model[i], "]\n",sep=""))
                else
                    ret[i] = FALSE
            }
        }
    } else {
        if (stop.on.error)
            stop(paste("\n\tModel [", model, "] is not of type character()\n",sep=""))
    }
    return (ret)
}
inla.is.model = function(model, stop.on.error = TRUE, ignore.case = FALSE)
{
    inla.is.generic(model, stop.on.error, inla.models()$models, ignore.case)
}
inla.is.lmodel = function(lmodel, stop.on.error = TRUE, ignore.case = FALSE)
{
    inla.is.generic(lmodel, stop.on.error, inla.models()$lmodels, ignore.case)
}
inla.is.prior = function(prior, stop.on.error = TRUE, ignore.case = FALSE)
{
    inla.is.generic(prior, stop.on.error, inla.models()$priors, ignore.case)
}

inla.model.properties.generic = function(model, stop.on.error, models, ignore.case)
{
    if (length(model) > 1)
        stop("...only implemented for length(model)=1; FIXME")
    
    m = ifelse(ignore.case, tolower(model), model)
    if (ignore.case)
        ms = tolower(names(models))
    else
        ms = names(models)

    if (inla.is.generic(model, stop.on.error, models, ignore.case)) {
        k = grep(paste("^", m, "$", sep=""), ms)
        return (models[[k]])
    }
    else
        return (FALSE)
}
inla.model.properties = function(model, stop.on.error = TRUE, ignore.case = FALSE)
{
    inla.model.properties.generic(inla.trim.family(model), stop.on.error, inla.models()$models, ignore.case)
}
inla.lmodel.properties = function(lmodel, stop.on.error = TRUE, ignore.case = FALSE)
{
    inla.model.properties.generic(inla.trim.family(lmodel), stop.on.error, inla.models()$lmodels, ignore.case)
}
inla.prior.properties = function(prior, stop.on.error = TRUE, ignore.case = FALSE)
{
    inla.model.properties.generic(inla.trim.family(prior), stop.on.error, inla.models()$priors, ignore.case)
}


