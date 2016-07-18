## 'spde' model functions
## Export: inla.spde.precision.inla.spde1 inla.spde.result.inla.spde1
## Export: inla.spde.precision!inla.spde1 inla.spde.result!inla.spde1
## Export: inla.spde1.create inla.spde.create
## Export: inla.spde1.imatern inla.spde1.matern
## Export: inla.spde1.matern.osc inla.spde1.models inla.spde1.precision
## Export: inla.spde1.result

inla.spde1.create =
    function(mesh,
             model=c("matern", "imatern", "matern.osc"),
             param=NULL,
             ...)
{
    inla.require.inherits(mesh, "inla.mesh", "'mesh'")

    model = match.arg(model)
    if (is.null(param))
        param = list()

    spde = (list(model = model,
                 mesh = mesh,
                 internal = list(),
                 f = (list(model="spde",
                           n=nrow(mesh$loc)
                           ))
                 ))
    class(spde) = c("inla.spde1", "inla.spde", "inla.model.class")

    if (identical(model, "matern") ||
        identical(model, "imatern") ||
        identical(model, "matern.osc")
        ) {

        if (is.null(param))
            param = list()
        if (is.null(param$alpha))
            param$alpha = 2
        if (is.null(param$basis.T))
            param$basis.T = matrix(1, mesh$n, 1)
        else if (!is.matrix(param$basis.T)) {
            len = length(as.vector(param$basis.T))
            if (len == 1L)
                param$basis.K = matrix(as.vector(param$basis.T), mesh$n, 1)
            else
                param$basis.T = as.matrix(param$basis.T)
        }
        if (nrow(param$basis.T) != mesh$n)
            stop(paste("'basis.T' has ", nrow(param$basis.T),
                       " rows; expected ", mesh$n, ".", sep=""))
        if (identical(model, "matern") ||
            identical(model, "matern.osc")
            ) {
            if (is.null(param$basis.K)) {
                param$basis.K = matrix(1, mesh$n, 1)
            } else if (!is.matrix(param$basis.K)) {
                len = length(as.vector(param$basis.K))
                if (len == 1L)
                    param$basis.K = matrix(as.vector(param$basis.K), mesh$n, 1)
                else
                    param$basis.K = as.matrix(param$basis.K)
            }
        } else {
            param$basis.K = matrix(0, mesh$n, 1)
        }
        if (nrow(param$basis.K) != mesh$n)
            stop(paste("'basis.K' has ", nrow(param$basis.K),
                       " rows; expected ", mesh$n, ".", sep=""))
        spde$internal = (c(spde$internal,
                           list(alpha = param$alpha,
                                basis.T = param$basis.T,
                                basis.K = param$basis.K)
                           ))

        mesh.range = (max(c(diff(range(mesh$loc[,1])),
                            diff(range(mesh$loc[,2])),
                            diff(range(mesh$loc[,3]))
                            )))

        spde$internal = (c(spde$internal,
                           inla.fmesher.smorg(mesh$loc,
                                              mesh$graph$tv,
                                              fem=2,
                                              output=list("c0", "g1", "g2"))))

        if (param$alpha==2) {
            kappa0 = sqrt(8)/(mesh.range*0.2)
            tau0 = 1/sqrt(4*pi*kappa0^2)/1.0
        } else if (param$alpha==1) {
            spde$internal$g2 = spde$internal$g1
            spde$internal$g1 = spde$internal$g1*0.0

            kappa0 = sqrt(sqrt(8)/(mesh.range*0.2))
            tau0 = 1/sqrt(4*pi)/1.0
        }

        if (identical(model, "matern")) {
            spde$f$hyper.default = (list(theta1=(list(initial=log(tau0),
                                                      param=c(log(tau0), 0.1))),
                                         theta2=(list(initial=log(kappa0^2),
                                                      param=c(log(kappa0^2), 0.1))),
                                         theta3=(list(initial=(log(tau0)+log(kappa0^2))/2,
                                                      param=c((log(tau0)+log(kappa0^2))/2, 0.1)))
                                         ))
        } else if (identical(model, "imatern")) {
            spde$f$hyper.default = (list(theta1=(list(initial=log(tau0),
                                                      param=c(log(tau0), 0.1))),
                                         theta2=(list(initial=-20,
                                                      fixed=TRUE)),
                                         theta3=(list(initial=log(tau0),
                                                      param=c(log(tau0), 0.1)))
                                         ))
        } else if (identical(model, "matern.osc")) {
            spde$f$hyper.default = (list(theta1=(list(initial=log(tau0),
                                                      param=c(log(tau0), 0.1))),
                                         theta2=(list(initial=log(kappa0^2),
                                                      param=c(log(kappa0^2), 0.1))),
                                         theta3=(list(initial=(log(tau0)+log(kappa0^2))/2,
                                                      param=c((log(tau0)+log(kappa0^2))/2, 0.1))),
                                         theta4=(list(fixed=FALSE,
                                                      initial=0,
                                                      param=c(0, 0.01)))
                                         ))
        }
    } else {
        stop(paste("Model '", model, "' unknown or not implemented.", sep=""))
    }

    return(invisible(spde))
}








inla.spde1.query = function(spde, ...)
{
    inla.require.inherits(spde, "inla.spde1", "'spde'")

    not.known = function (spde, queryname)
    {
        stop(paste("Query '", queryname,
                   "' unknown.", sep=""))
    }
    not.implemented = function (spde, queryname)
    {
        stop(paste("Query '", queryname,
                   "' not implemented for inla.spde1 model '",
                   spde$model, "'.", sep=""))
    }
    param.to.fcn =
        function(basis, theta, values, n, name)
        {
            if (!is.null(values)) {
                if ((is.vector(values) &&
                     (length(values)==n)) ||
                    (is.matrix(values) &&
                     (nrow(values)==n))) {
                    fcn = as.vector(values)
                } else {
                    fcn = as.vector(values)
                }
                if (length(fcn) == 1L)
                    fcn = rep(fcn, n)
                else if (length(fcn) != n)
                    stop(paste("Length of '", name, "' is ", length(fcn),
                               ", should be ", n, ".", sep=""))
            } else {
                fcn = exp(basis %*% theta)
            }
            return(fcn)
        }

    result = list()
    queries = inla.parse.queries(...)
    if (length(queries)==0L)
        return(result)

    for (query.idx in 1:length(queries)) {
        query = names(queries)[query.idx]
        param = queries[[query.idx]]
        answer = NULL
        query = (match.arg(query, c("precision",
                                    "sample"
                                    )))
        if (identical(query, "precision")) {
            if (identical(spde$model, "matern")) {
                tau = (param.to.fcn(spde$internal$basis.T,
                                    param$theta.T, param$tau,
                                    spde$mesh$n, "tau"))
                kappa2 = (param.to.fcn(spde$internal$basis.K,
                                       param$theta.K, param$kappa2,
                                       spde$mesh$n, "kappa2"))
                dT = Diagonal(spde$mesh$n, tau)
                dK2 = Diagonal(spde$mesh$n, kappa2)
                tmp = dK2 %*% spde$internal$g1
                answer = (dT %*% (dK2 %*% spde$internal$c0 %*% dK2+
                                  tmp + t(tmp) +
                                  spde$internal$g2) %*% dT)
            } else if (identical(spde$model, "matern.osc")) {
                tau = (param.to.fcn(spde$internal$basis.T,
                                    param$theta.T, param$tau,
                                    spde$mesh$n, "tau"))
                kappa2 = (param.to.fcn(spde$internal$basis.K,
                                       param$theta.K, param$kappa2,
                                       spde$mesh$n, "kappa2"))
                osc = param$osc
                dT = Diagonal(spde$mesh$n, tau)
                dK2 = Diagonal(spde$mesh$n, kappa2)
                tmp = dK2 %*% spde$internal$g1
                answer = (dT %*% (dK2 %*% spde$internal$c0 %*% dK2+
                                  cos(pi*osc)*(tmp + t(tmp)) +
                                  spde$internal$g2) %*% dT)
            } else if (identical(spde$model, "imatern")) {
                tau = (param.to.fcn(spde$internal$basis.T,
                                    param$theta.T, param$tau,
                                    spde$mesh$n, "tau"))
                dT = Diagonal(spde$mesh$n, tau)
                answer = (dT %*% ((1e-10)*spde$internal$c0+
                                  spde$internal$g2) %*% dT)
            } else {
                not.implemented(spde,query)
            }
        } else if (identical(query, "sample")) {
            Q = inla.spde1.query(spde, precision=param)$precision
            finn = (inla.finn(Q, seed=(inla.ifelse(is.null(param$seed),
                                                   0L,
                                                   param$seed))))
            answer = finn$sample
        } else if (!identical(query, "")) {
            not.known(spde,query)
        }
        ## Expand the result list:
        result[query.idx] = list(NULL)
        names(result)[query.idx] = query
        ## Set the answer:
        if (!is.null(answer))
            result[[query.idx]] = answer
    }

    return(result)
}





inla.spde1.result = function(inla, name, spde, do.transform=TRUE, ...)
{
    warning("'inla.spde1.result' is not fully implemented yet.")
    inla.require.inherits(inla, "inla", "'inla'")
    inla.require.inherits(spde, "inla.spde1", "'spde'")
    if (!spde$f$model=="spde") {
        stop("'inla.spde1.result' only supports internal inla models 'spde'")
    }

    result = list()
    ## Setup rownames for UserFunction1
    UF.names =
        c(rownames(rep(list(NULL), spde$mesh$n),
                   do.NULL=FALSE, prefix="tau."),
          rownames(rep(list(NULL), spde$mesh$n),
                   do.NULL=FALSE, prefix="kappa2."))
    rownames(inla$summary.random[["UserFunction1"]]) <- UF.names

    if (!is.null(inla$marginals.random[["UserFunction1"]])) {
        names(inla$marginals.random[["UserFunction1"]]) <-
            UF.names
    }

    ## Values
    result$summary.values = inla$summary.random[[name]]

    ## Marginals for values
    if (!is.null(inla$marginals.random[[name]])) {
        result$marginals.values = inla$marginals.random[[name]]
    }

    ## Theta
    result$summary.theta.tau =
        inla.extract.el(inla$summary.hyperpar,
                        paste("^T\\.[^ ]+ for ", name, "-basisT$",
                              sep=""))
    result$summary.theta.kappa2 =
        inla.extract.el(inla$summary.hyperpar,
                        paste("^K\\.[^ ]+ for ", name, "-basisK$",
                              sep=""))
    rownames(result$summary.theta.tau) <-
        rownames(rep(list(NULL), nrow(result$summary.theta.tau)),
                 do.NULL=FALSE, prefix="theta.")
    rownames(result$summary.theta.kappa2) <-
        rownames(rep(list(NULL), nrow(result$summary.theta.kappa2)),
                 do.NULL=FALSE, prefix="theta.")

    ## Marginals for theta
    if (!is.null(inla$marginals.hyperpar)) {
        result$marginals.theta.tau =
            inla.extract.el(inla$marginals.hyperpar,
                            paste("^T\\.[^ ]+ for ", name, "-basisT$",
                                  sep=""))
        result$marginals.theta.kappa2 =
            inla.extract.el(inla$marginals.hyperpar,
                            paste("^K\\.[^ ]+ for ", name, "-basisK$",
                                  sep=""))
        names(result$marginals.theta.tau) <-
            rownames(rep(list(NULL), length(result$marginals.theta.tau)),
                     do.NULL=FALSE, prefix="theta.")
        names(result$marginals.theta.kappa2) <-
            rownames(rep(list(NULL), length(result$marginals.theta.kappa2)),
                     do.NULL=FALSE, prefix="theta.")
    }

    ## log-tau/kappa2
    result$summary.log.tau =
        inla.extract.el(inla$summary.random[["UserFunction1"]],
                        paste("^tau\\.[^ ]+$",
                              sep=""))
    result$summary.log.kappa2 =
        inla.extract.el(inla$summary.random[["UserFunction1"]],
                        paste("^kappa2\\.[^ ]+$",
                              sep=""))

    ## Approximate log-range:
    forward.names =
        c("mean", "0.025quant", "0.5quant", "0.975quant")
    reverse.names =
        c("mean", "0.975quant", "0.5quant", "0.025quant")
    result$summary.log.range.nominal = result$summary.log.kappa2
    result$summary.log.range.nominal[,forward.names] =
        (log(sqrt(8)) -
         result$summary.log.kappa2[,reverse.names]/2)
    result$summary.log.range.nominal[,"sd"] =
        result$summary.log.kappa2[,"sd"]/4

    ## Approximate log-variance:
    result$summary.log.variance.nominal = result$summary.log.kappa2
    result$summary.log.variance.nominal[,forward.names] =
        -(log(4*pi) +
          result$summary.log.tau[,reverse.names]*2 +
          result$summary.log.kappa2[,reverse.names])
    result$summary.log.range.nominal[,"sd"] =
        sqrt(result$summary.log.tau[,"sd"]^2*4 +
             result$summary.log.kappa2[,"sd"]^2 +
             2*2*result$summary.log.tau[,"sd"]*
             result$summary.log.kappa2[,"sd"])

    ## Marginals for log-tau/kappa2
    if (!is.null(inla$marginals.random[["UserFunction1"]])) {
        result$marginals.log.tau =
            inla.extract.el(inla$marginals.random[["UserFunction1"]],
                            paste("^tau\\.[^ ]+$",
                                  sep=""))
        result$marginals.log.kappa2 =
            inla.extract.el(inla$marginals.random[["UserFunction1"]],
                            paste("^kappa2\\.[^ ]+$",
                                  sep=""))

        if (do.transform) {
            result$marginals.tau =
                lapply(result$marginals.log.tau[1],
                       function(x) inla.tmarginal(function(y) exp(y), x))
            result$marginals.kappa =
                lapply(result$marginals.log.kappa2[1],
                       function(x) inla.tmarginal(function(y) exp(y/2), x))
##            result$marginals.variance.nominal =
##                lapply(result$marginals.log.variance.nominal,
##                       function(x) inla.tmarginal(function(y) exp(y), x))
            result$marginals.range.nominal =
                lapply(result$marginals.log.kappa2[1],
                       function(x) inla.tmarginal(function(y) sqrt(8)*exp(-y/2), x))
        }

    }

    return(result)
}


inla.spde1.precision = function(spde, ...)
{
    inla.require.inherits(spde, "inla.spde1", "'spde'")

    return(inla.spde1.query(spde, precision=list(...))$precision)
}


inla.spde1.models = function()
{
    return(c("matern", "imatern", "matern.osc"))
}

inla.spde1.matern = function(mesh, ...)
{
    return(inla.spde1.create(mesh, model="matern", ...))
}
inla.spde1.imatern = function(mesh, ...)
{
    return(inla.spde1.create(mesh, model="imatern", ...))
}
inla.spde1.matern.osc = function(mesh, ...)
{
    return(inla.spde1.create(mesh, model="matern.osc", ...))
}

## spde.common-connections:
inla.spde.precision.inla.spde1 = inla.spde1.precision
inla.spde.result.inla.spde1 = inla.spde1.result

## Backwards compatibility:
inla.spde.create = inla.spde1.create
