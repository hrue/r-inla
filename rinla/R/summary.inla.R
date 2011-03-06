`summary.inla` = function(object, ...)
{
    ## provides a summary for a inla object
    ret = list()
    ret = c(ret, list(call=object$call))
    ret = c(ret, list(cpu.used = object$cpu.used))

    if(!is.null(object$summary.fixed)) {
        fixed = object$summary.fixed
        if(!is.null(object$kld.fixed))
            fixed = cbind(fixed, kld=object$kld.fixed)
        
        ret = c(ret, list(fixed=fixed))
    }

    if(!is.null(object$summary.lincomb)) {
        lincomb = object$summary.lincomb
        if(!is.null(object$kld.lincomb))
            lincomb = cbind(lincomb, kld=object$kld.lincomb)
        
        ret = c(ret, list(lincomb=lincomb))
    }

    if(!is.null(object$summary.hyperpar))
        ret = c(ret, list(hyperpar=object$summary.hyperpar))
    

    if(!is.null(object$summary.random)) {
        random.name = names(object$summary.random)
        random.mod =  object$model.random
        random.max.kld = numeric(length(random.name))
        for (i in 1:length(random.name))
            random.max.kld[i] = max(object$summary.random[[i]]$kld)
        ret = c(ret, list(random.names=random.name, random.model=random.mod, random.max.kld=random.max.kld))
    }
    
    
    neffp = object$neffp
    ret = c(ret, list(neffp=neffp))
    
    if (!is.null(object$dic)) {
        ret = c(ret, list(dic=object$dic))
    }

    if(!is.null(object$mlik))
        ret = c(ret, list(mlik=object$mlik))
    
    if(!is.null(object$pit))
        ret = c(ret, list(pit=object$pit, cpo=object$cpo))

    if(!is.null(object$summary.linear.predictor))
        ret = c(ret, list(linear.predictor= object$summary.linear.predictor))
    
    ret = c(ret, list(family=object$family))
    class(ret) = "summary.inla"

    return (ret)
}

`summary.hyperpar.inla` =
    function(object,...)
{
    ##provides a summary for a hyperpar inla object
    ret = object$summary
    return (ret)
}

