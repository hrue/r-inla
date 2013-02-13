`summary.inla` = function(object, ...)
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
    ret = c(ret, list(cpu.used = object$cpu.used))

    if(!is.null(object$summary.fixed)) {
        fixed = object$summary.fixed
        ret = c(ret, list(fixed=fixed))
    }

    if(!is.null(object$summary.lincomb) && any(names(object) == "summary.lincomb")) {
        lincomb = as.matrix(object$summary.lincomb)
        ret = c(ret, list(lincomb=lincomb))
    }

    if(!is.null(object$summary.lincomb.derived)) {
        lincomb = as.matrix(object$summary.lincomb.derived)
        ret = c(ret, list(lincomb.derived=lincomb))
    }

    if(!is.null(object$summary.hyperpar))
        ret = c(ret, list(hyperpar=object$summary.hyperpar))
    

    if(!is.null(object$summary.random)) {
        random.name = names(object$summary.random)
        random.mod =  object$model.random
        ret = c(ret, list(random.names=random.name, random.model=random.mod))
    }
    
    neffp = object$neffp
    ret = c(ret, list(neffp=neffp))
    
    if (!is.null(object$dic)) {
        ret = c(ret, list(dic=object$dic))
    }

    if(!is.null(object$mlik))
        ret = c(ret, list(mlik=object$mlik))
    
    if(!is.null(object$cpo$cpo))
        ret = c(ret, list(cpo=object$cpo))

    if(!is.null(object$summary.linear.predictor))
        ret = c(ret, list(linear.predictor= object$summary.linear.predictor))
    
    ret = c(ret, list(family=object$family))
    class(ret) = "summary.inla"

    return (ret)
}
