`summary.inla` = function(object, digits = 4, ...)
{
    ## provides a summary for a inla object
    ret = list()

    maxlen = 2048L
    if (sum(nchar(object$call)) > maxlen) {
        ret = c(ret,
                list(call=paste(substr(inla.paste(object$call), 1L, maxlen),
                             "<<<the rest is not shown>>>")))
    } else {
        ret = c(ret, list(call=object$call))
    }

    ## might not be if using collect directly
    if (inla.is.element("cpu.used", object)) {
        ret = c(ret, list(cpu.used = round(object$cpu.used, digits))) 
    } else {
        ret = c(ret,  list(cpu.used = NA))
    }

    if(!is.null(object$summary.fixed)) {
        fixed = round(object$summary.fixed, digits)
        ret = c(ret, list(fixed=fixed))
    }

    if(!is.null(object$summary.lincomb) && any(names(object) == "summary.lincomb")) {
        lincomb = round(as.matrix(object$summary.lincomb), digits)
        ret = c(ret, list(lincomb=lincomb))
    }

    if(!is.null(object$summary.lincomb.derived)) {
        lincomb = round(as.matrix(object$summary.lincomb.derived), digits)
        ret = c(ret, list(lincomb.derived=lincomb))
    }

    if(!is.null(object$summary.hyperpar))
        ret = c(ret, list(hyperpar=round(object$summary.hyperpar, digits)))
    

    if(!is.null(object$summary.random)) {
        random.name = names(object$summary.random)
        random.mod =  object$model.random
        ret = c(ret, list(random.names=random.name, random.model=random.mod))
    }
    
    neffp = object$neffp
    ret = c(ret, list(neffp = round(neffp, digits)))
    
    if (!is.null(object$dic)) {
        ret = c(ret, list(dic = lapply(object$dic, round, digits = digits)))
    }

    if(!is.null(object$mlik))
        ret = c(ret, list(mlik = round(object$mlik, digits)))
    
    if(!is.null(object$cpo$cpo) && length(object$cpo$cpo) > 0L)
        ret = c(ret, list(cpo = lapply(object$cpo, round, digits = digits)))

    if(!is.null(object$summary.linear.predictor))
        ret = c(ret, list(linear.predictor= round(object$summary.linear.predictor, digits)))
    
    ret = c(ret, list(family=object$family))
    class(ret) = "summary.inla"

    return (ret)
}
