## 'spde' model functions


inla.spde.create =
    function(mesh,
             model=c("matern", "imatern", "matern.osc"),
             param=NULL,
             ...)
{
    inla.require.inherits(mesh, "inla.mesh", "'mesh'")

    model = match.arg(model)
    if (is.null(param))
        param = list()

    spde.prefix = inla.fmesher.make.prefix(NULL, NULL)

    spde = (list(model = model,
                 mesh = mesh,
                 internal = list(),
                 f = (list(model="spde",
                           spde.prefix=spde.prefix,
                           n=nrow(mesh$loc)
                           ))
                 ))
    class(spde) = "inla.spde"

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
            stop(paste("'basis.T' has ", nrow(basis.T),
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
            stop(paste("'basis.K' has ", nrow(basis.K),
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
        ## inla checks PREFIX valididy by looking for "s":
        fmesher.write(spde$mesh$loc, spde.prefix, "s")
        ## Write the precision building blocks:
        fmesher.write(spde$internal$c0, spde.prefix, "c0")
        fmesher.write(spde$internal$g1, spde.prefix, "g1")
        fmesher.write(spde$internal$g2, spde.prefix, "g2")
        fmesher.write(spde$internal$basis.T, spde.prefix, "basisT")
        fmesher.write(spde$internal$basis.K, spde.prefix, "basisK")

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








inla.spde.query = function(spde, ...)
{
    inla.require.inherits(spde, "inla.spde", "'spde'")

    not.known = function (spde, queryname)
    {
        stop(paste("Query '", queryname,
                   "' unknown.", sep=""))
    }
    not.implemented = function (spde, queryname)
    {
        stop(paste("Query '", queryname,
                   "' not implemented for inla.spde model '",
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
            Q = inla.spde.query(spde, precision=param)$precision
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




## Deprecated functions:

inla.spde = function(...)
{
    args = list(...)
    if (length(args)>0) {
        if (inherits(args[[1]],"inla.spde")) {
            warning("'inla.spde(spde, ...)' is deprecated.  Use 'inla.spde.query(spde, ...)' instead.")
            return(inla.spde.query(...))
        }
    }
    warning("'inla.spde(...)' is deprecated.  Use 'inla.spde.create(...)' instead.")
    return(inla.spde.create(...))
}

