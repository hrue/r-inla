## Nothing to export

inla.extract.Q = function(what, formula,  data, ...)
{
    ## extract the Q for 'what'. this is an interface so that we can extract a Q-matrix of any
    ## f() model
    stopifnot(!missing(what))
    stopifnot(!missing(formula))
    stopifnot(!missing(data))
    
    y.name = as.character(formula)[2]
    i = which(names(data) == y.name)
    stopifnot(length(i) == 1L)
    if (is.matrix(data[[i]])) {
        data[[i]] = rep(0, nrow = nrow(data[[i]]))
    } else if (is.vector(data[[i]])) {
        data[[i]] = rep(0, length(data[[i]]))
    } else if (is.list(data[[i]])) {
        stop("data as list, is not yet implemented")
    }

    log.prec = 0  ## the fixed noise level
    args = list(...)
    args$control.compute = list(config=TRUE)
    args$control.fixed = list(prec = exp(log.prec), prec.intercept = exp(log.prec))
    args$family = "gaussian"
    args$control.family = list(hyper = list(
                                   prec = list(
                                       initial = log.prec,
                                       fixed=TRUE)),
                               variant = 1, 
                               control.link = list(
                                   model = "identity"))
    args$formula = formula
    args$data = data

    result = do.call("inla", args)
    conf = result$misc$configs
    cont = conf$contents
    i = which(cont$tag == what)
    if (length(i) != 1L)
        stop(paste0("Cannot find '", what, "'"))
    p.len = cont$length[1]
    i.start = cont$start[i]
    i.len = cont$length[i]

    k = NA
    for (j in seq_along(conf$config)) {
        if (conf$config[[j]]$log.posterior.orig == 0.0) {
            k = j
            break
        }
    }
    stopifnot(!is.na(k))
    Q = conf$config[[k]]$Q
    d = diag(Q)
    for(i in i.start:(i.start + i.len - 1L)) {
        Q[i, i] = Q[i, i] - sum(Q[1:p.len, i]^2 / (d[1:p.len] - exp(log.prec)))
    }
    Q = Q[i.start:(i.start + i.len -1), i.start:(i.start + i.len -1)]
    Q = Q + t(Q)
    diag(Q) = diag(Q) / 2.0

    return (Q)
}
