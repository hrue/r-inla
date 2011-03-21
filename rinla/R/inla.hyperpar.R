##!\name{inla.hyperpar}
##!\alias{inla.hyperpar}
##!\alias{hyperpar.inla}
##!
##!\title{Improved estimates for the hyperparameters}
##!\description{
##! Improve the estimates of the  posterior marginals
##! for the hyperparameters of
##! the model using the grid integration strategy.
##!}
##!\usage{
##!`inla.hyperpar` =
##!    function(object, skip.configurations, verbose,  dz, diff.logdens, h, restart, quantiles,  keep)
##!}
##!
##!\arguments{

`inla.hyperpar` =
    function(

             ##!\item{object}{An object of class
             ##!\code{inla},  ie a result of a call to
             ##!\code{inla()}}
             object,

             ##!\item{skip.configurations}{ A boolean variable; skip
             ##!configurations if the values at the main axis are to
             ##!small. (Default TRUE)}
             skip.configurations = TRUE,

             ##!\item{verbose}{Boolean indicating wheather the inla
             ##!program should run in a verbose mode.}
             verbose = FALSE,

             ##!\item{dz}{Step length in the standardized scale used
             ##!in the construction of the grid, default 0.75.}
             dz = 0.75,

             ##!\item{diff.logdens}{The difference of the
             ##!log.density for the hyperpameters to stop numerical
             ##!integration using int.strategy='grid'.  Default 7}
             diff.logdens = 7,

             ##!\item{h}{The step-length for the gradient
             ##!calculations for the hyperparameters. Default 0.01.}
             h = NULL,

             ##!\item{restart}{A boolean defining wheather the
             ##!optimizer should start again to ind the mode or if it
             ##!should use the mode contained in the \code{object}}
             restart = FALSE,

             ##!\item{quantiles}{A vector of quantiles,
             ##!to compute for each posterior marginal.}
             quantiles = c(0.025, 0.5, 0.975),

             ##!\item{keep}{A boolean variable indicating the
             ##!working files (ini file, data files and results files)
             ##!should be kept}
             keep = FALSE
             )
{
    ##!}
    ##!\value{
    ##!The object returned is the same as  \code{object} but the
    ##!estimates of the hyperparameters are replaced by improved estimates.}
    ##!\references{
    ##!See the references in \code{inla}}
    ##!\author{Havard Rue \email{hrue@math.ntnu.no}}
    ##!\note{This function might take a long time or if the number of
    ##!hyperparameters in the model is large. If it complains and says \code{I cannot get enough memory},
    ##!try to increase the value of the argument \code{dz} or decrease \code{diff.logdens}.}
    ##!\seealso{\code{\link{inla}}}

    stopifnot(class(object) == "inla")
    
    verbose.save = verbose
    dz.save = dz
    quantiles.save = quantiles
    diff.logdens.save = diff.logdens
    keep.save = keep

    for (nm in names(object))
        inla.eval(paste(nm, " = object$", nm, sep=""))

    verbose = verbose.save
    dz = dz.save
    quantiles = quantiles.save
    diff.logdens = diff.logdens.save
    keep = keep.save

    ## if these are not in the `object', then they are NULL
    if (!is.element("E", names(object)))
        E=NULL
    if (!is.element("offset", names(object)))
        offset=NULL
    if (!is.element("Ntrials", names(object)))
        Ntrials=NULL
    if (!is.element("scale", names(object)))
        scale=NULL
    if (!is.element("inla.arg", names(object)))
        inla.arg=NULL
    if (!is.element("lincomb", names(object)))
        lincomb=NULL
    if (!is.element("lincomb", names(object)))
        lincomb=NULL
    if (!is.element("control.lincomb", names(object)))
        control.lincomb=list()

    ## call inla() with the attached variables. Note that variables in
    ## the call of this function, has precedance over the attached'
    ## ones.

    result = inla(
            formula = as.formula(formula),
            family = family, 
            data = data,
            quantiles = quantiles,
            E = E,
            offset = offset,
            Ntrials = Ntrials,
            scale = scale,
            verbose = verbose,      
            lincomb = lincomb,
            control.lincomb = control.lincomb,
            control.compute=list(hyperpar=TRUE, dic=FALSE,
                    mlik=TRUE, cpo=FALSE,
                    smtp=control.compute$smtp,
                    strategy = control.compute$strategy),
            control.predictor = list(compute=FALSE,
                    cdf=NULL, quantiles=NULL,
                    hyper = control.predictor$hyper, 
                    cross = control.predictor$cross,
                    A = control.predictor$A, 
                    precision = control.predictor$precision), 
            control.data = control.data,
            control.inla = list(strategy=NULL,
                    int.strategy="grid",
                    h=inla.ifelse(is.null(h), control.inla$h, h),
                    dz=dz,
                    diff.logdens= diff.logdens,
                    print.joint.hyper=TRUE,
                    force.diagonal=TRUE,
                    adjust.weights = FALSE,
                    skip.configurations=skip.configurations,
                    tolerance = control.inla$tolerance,
                    lincomb.derived.only = control.inla$lincomb.derived.only,
                    step.len = control.inla$step.len),
            control.fixed = control.fixed,
            control.results = list(return.marginals.random=FALSE,
                    return.marginals.predictor=FALSE),
            control.mode = list(result = object, restart=restart),
            inla.call = inla.call,
            inla.arg = inla.arg,
            only.hyperparam = TRUE,
            keep = keep,
            working.directory = NULL,
            silent = silent,
            user.hook = inla.ifelse(exists("user.hook"), user.hook, NULL),
            user.hook.arg = inla.ifelse(exists("user.hook.arg"), user.hook.arg, NULL),
            ##
            .internal = .internal
            )

    ## these are the entries that we want to replace
    replace.names = c("summary.hyperpar", "marginals.hyperpar", "internal.marginals.hyperpar",
            "joint.hyper", "mlik", "version", "control.inla")

    for (nm in replace.names) {
        idx.result = which(names(result) == nm)
        idx.object = which(names(object) == nm)
        object[[idx.object]] = result[[idx.result]]
    }

    ## these are old names that was used before as the output of
    ## inla.hyperpar(). add them with an informative value to ease the
    ## transition...
    old.names = c("summary", "marginals", "internal.marginals", "log.joint")
    for (nm in old.names) {
        inla.eval(paste("object$", nm, "=", "'", 
                        "Variable ...$", nm, " does not longer exists. The output of inla.hyperpar(...) is now the same as the output of inla(...).",
                        "'",  sep=""))
    }

    return(object)
}
