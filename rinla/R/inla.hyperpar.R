
`hyperpar.inla` =
function(...)
{
    cat("\n\n*** WARNING ***  hyperpar.inla() is now called inla.hyperpar()")
    return (inla.hyperpar(...))
}

`inla.hyperpar` =
    function(object,
             skip.configurations = TRUE,
             verbose = FALSE,
             dz = 0.75,
             diff.logdens = 7,
             h = NULL, 
             restart = FALSE,
             quantiles = c(0.025, 0.5, 0.975),
             keep = FALSE
             )
{
    if(class(object) == "inla")
        rr = inla.hyperpar.inla(object, skip.configurations = skip.configurations,
                verbose = verbose, dz = dz, diff.logdens = diff.logdens,
                h = h, restart = restart, quantiles = quantiles, keep = keep)
    else
        stop("Object class not identified")
    
    class(rr) = "hyperpar.inla"
    return(rr)
}


### wrapper for the default models
`inla.hyperpar.inla` =
    function(object,
             skip.configurations,
             verbose,
             dz,
             h,
             restart,
             quantiles,
             diff.logdens,
             keep
             )
{
    ## this function is just a wrapper for inla() as it just change some
    ## parameters to better compute the hyperparameters, and sets the
    ## computed mode to be known (using the type = mode section)

    ## ... we need to evaluate object within here, to prevent to pass
    ## object along with data into inla(). this creates problems in
    ## the eval() part in inla(); so we do this the rough-way:

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

    rr <- inla(formula = as.formula(formula),
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
               control.compute=list(hyperpar=TRUE,dic=FALSE,
                       mlik=TRUE,cpo=FALSE,
                       smtp=control.compute$smtp),
               control.predictor = list(compute=FALSE,
                       cdf=NULL, quantiles=NULL,
                       fixed=control.predictor$fixed,
                       param=control.predictor$param,
                       initial = control.predictor$initial),
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
                       derived.only = control.inla$derived.only,
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

    ret = list(summary = rr$summary.hyperpar, marginals = rr$marginals.hyperpar, 
            internal.marginals = rr$internal.marginals.hyperpar, log.joint = rr$joint.hyper, mlik = rr$mlik,
            version = rr$version)
    return(ret)
}


