## 'spde2' model functions




homogenise_B_matrix = function(B, n.spde, n.theta)
{
    if (!is.numeric(B))
        stop("B matrix must be numeric.")
    if (is.matrix(B)) {
        if ((nrow(B) != 1) && (nrow(B) != n.spde)) {
            stop(inla.paste(list("B matrix has",
                                 as.character(nrow(B)),
                                 "rows but should have 1 or",
                                 as.character(n.spde),
                                 sep=" ")))
        }
        if ((ncol(B) != 1) && (ncol(B) != 1+n.theta)) {
            stop(inla.paste(list("B matrix has",
                                 as.character(ncol(B)),
                                 "columns but should have 1 or",
                                 as.character(1+n.theta),
                                 sep=" ")))
        }
        if (ncol(B) == 1) {
            return(cbind(as.vector(B), matrix(0.0, n.spde, n.theta)))
        } else if (ncol(B) == 1+n.theta) {
            if (nrow(B) == 1) {
                return(matrix(as.vector(B), n.spde, 1+n.theta, byrow=TRUE))
            } else if (nrow(B) == n.spde) {
                return(B)
            }
        }
    } else { ## !is.matrix(B)
        if ((length(B) == 1) || (length(B) == n.spde)) {
            return(cbind(B, matrix(0.0, n.spde, n.theta)))
        } else if (length(B) == 1+n.theta) {
            return(matrix(B, n.spde, 1+n.theta, byrow=TRUE))
        } else {
            stop(inla.paste(list("Length of B vector is",
                                 as.character(length(B)),
                                 "but should be 1,",
                                 as.character(1+n.theta), "or",
                                 as.character(n.spde)),
                            sep=" "))
        }
    }
    stop(inla.paste(list("Unrecognised structure for B matrix"),
                    sep=" "))
}


inla.spde.generic2 =
    function(M0, M1, M2, B0, B1, B2, theta.mu, theta.Q,
             transform = c("logit","log","identity"),
             theta.initial = theta.mu,
             fixed = rep(FALSE, length(theta.mu)),
             theta.fixed = theta.initial[fixed],
             BLC = cbind(0.0, diag(nrow=length(theta.mu))),
             ...)
{
    transform = match.arg(transform)

    n.spde=nrow(M0)
    n.theta = length(theta.mu)

    spde = (list(model = "spde2.generic",
                 internal = list(),
                 f = list()
                 ))
    class(spde) = "inla.spde"

    B0 = homogenise_B_matrix(B0, n.spde, n.theta)
    B1 = homogenise_B_matrix(B1, n.spde, n.theta)
    B2 = homogenise_B_matrix(B2, n.spde, n.theta)
    BLC = homogenise_B_matrix(BLC, nrow(BLC), n.theta)

    param.generic =
        list(M0=M0, M1=M1, M2=M2,
             B0=B0, B1=B1, B2=B2, BLC=BLC,
             theta.mu=theta.mu, theta.Q=theta.Q,
             transform=transform,
             theta.initial=theta.initial,
             fixed=fixed,
             theta.fixed=theta.fixed)
    ## Copy full generic parameters to internal inla-representation.
    param.inla = param.generic
    ## TODO: Check dimension consistency of inputs

    ## If any "fixed", move information from B[,1+(1:n.theta)] into
    ## B[,1], theta.mu, and theta.Q
    ## NOTE: Should this be in the f function instead?
    ##       No, the information is also needed elsewhere.
    if (any(fixed)) {
        param.inla$theta.initial = theta.initial[!fixed]
        param.inla$B0[,1] =
            B0[,c(TRUE, fixed)] %*% c(1.0, theta.fixed)
        param.inla$B1[,1] =
            B1[,c(TRUE, fixed)] %*% c(1.0, theta.fixed)
        param.inla$B2[,1] =
            B2[,c(TRUE, fixed)] %*% c(1.0, theta.fixed)
        param.inla$BLC[,1] =
            BLC[,c(TRUE, fixed)] %*% c(1.0, theta.fixed)
        param.inla$B0 = param.inla$B0[,c(TRUE, !fixed), drop=FALSE]
        param.inla$B1 = param.inla$B1[,c(TRUE, !fixed), drop=FALSE]
        param.inla$B2 = param.inla$B2[,c(TRUE, !fixed), drop=FALSE]
        param.inla$BLC = param.inla$BLC[,c(TRUE, !fixed), drop=FALSE]
        param.inla$theta.Q = param.inla$theta.Q[!fixed, !fixed, drop=FALSE]
        if (!all(fixed)) {
            param.inla$theta.mu =
                theta.mu - solve(param.inla$theta.Q,
                                 theta.Q[!fixed, fixed, drop=FALSE]
                                 %*% (theta.fixed[fixed]-theta.mu[fixed]))
        } else { ## All fixed.
            param.inla$theta.mu = as.vector(matrix(0.0, 0, 1))
        }
    }

    spde$internal$param.generic = param.generic
    spde$internal$param.inla = param.inla

    ## NOTE: Should the prefix be set inside the f function instead?
    ##       Yes, it should.
    spde.prefix = inla.fmesher.make.prefix(NULL, NULL)

    spde$f = (list(model="spde2",
                   spde2.prefix=spde.prefix,
                   n=n.spde,
                   spde2.transform=transform,
                   hyper.default =
                   (list(theta1 =
                         list(prior="mvnorm",
                              param=(c(param.inla$theta.mu,
                                       as.matrix(param.inla$theta.Q)))
                              )
                         ))
                   ))

    ## NOTE: Should this be in the f function instead?
    ##       Yes, it should.
    fmesher.write(spde$internal$param.inla$M0, spde.prefix, "M0")
    fmesher.write(spde$internal$param.inla$M1, spde.prefix, "M1")
    fmesher.write(spde$internal$param.inla$M2, spde.prefix, "M2")
    fmesher.write(spde$internal$param.inla$B0, spde.prefix, "B0")
    fmesher.write(spde$internal$param.inla$B1, spde.prefix, "B1")
    fmesher.write(spde$internal$param.inla$B2, spde.prefix, "B2")
    fmesher.write(spde$internal$param.inla$BLC, spde.prefix, "BLC")

    return(spde)
}


inla.spde.matern =
    function(mesh,
             alpha=2,
             B.tau = matrix(c(0,1,0),1,3),
             B.kappa2 = matrix(c(0,0,1),1,3),
             B.prec = NULL,
             sigma2.approx = 1,
             theta.prior.mean = NULL,
             theta.prior.prec = NULL)
{
    if (is.null(B.kappa2))
        stop("B.kappa2 must be specified.")
    is.stationary = (nrow(B.kappa2)==1)

    d = 2
    nu = alpha-d/2

    n.spde = mesh$n
    n.theta = ncol(B.kappa2)-1L

    if (!is.null(B.tau)) {
        is.stationary = is.stationary && (nrow(B.tau)==1)
    }
    if (!is.null(B.prec)) {
        if (!is.null(B.tau))
            stop("Only one of B.tau and B.prec may be null.")
        is.stationary = is.stationary && (nrow(B.prec)==1)
        ## prec = 1/gamma(nu) * gamma(alpha)*(4*pi)^(d/2)*kappa2^nu * tau^2
        ## 2*log(tau) = - log(gamma(alpha)/gamma(nu)) - (d/2)*log(4*pi)
        ##             - nu*log(kappa2) + log(prec)
        ##           = - lgamma(alpha) + lgamma(nu) - (d/2)*log(4*pi)
        ##             + (- nu * B.kappa2 + B.prec) * c(1,theta)
        B.kappa2 = homogenise_B_matrix(B.kappa2, n.spde, n.theta)
        B.prec = homogenise_B_matrix(B.prec, n.spde, n.theta)

        B.tau =
            cbind(- lgamma(max(1.5,alpha)) + lgamma(max(0.5,nu))
                  - (d/2)*log(4*pi) - max(0.5,nu)*B.kappa2[,1] + B.prec[,1],
                  - max(0.5,nu)*B.kappa2[,-1,drop=FALSE]
                  + B.prec[,-1,drop=FALSE] )
    } else {
        if (is.null(B.tau))
            stop("At most one of B.tau and B.prec may be null.")
        B.kappa2 = homogenise_B_matrix(B.kappa2, n.spde, n.theta)
        B.tau = homogenise_B_matrix(B.tau, n.spde, n.theta)
        B.prec =
            cbind(+ lgamma(max(1.5,alpha)) - lgamma(max(0.5,nu))
                  + (d/2)*log(4*pi) + max(0.5,nu)*B.kappa2[,1] + B.tau[,1],
                  max(0.5,nu)*B.kappa2[,-1,drop=FALSE]
                  + B.tau[,-1,drop=FALSE] )
    }
    if (is.stationary) {
        B.tau = B.tau[1,,drop=FALSE]
        B.kappa2 = B.kappa2[1,,drop=FALSE]
        B.prec = B.prec[1,,drop=FALSE]
    }
    B.range =
        cbind(0.5*log(8*max(0.5,nu))-B.kappa2[,1]/2,
              -B.kappa2[,-1,drop=FALSE]/2)

    B.theta = cbind(0,diag(1, n.theta))
    rownames(B.theta) <- rownames(B.theta, do.NULL=FALSE, prefix="theta.")
    rownames(B.tau) <- rownames(B.tau, do.NULL=FALSE, prefix="tau.")
    rownames(B.kappa2) <- rownames(B.kappa2, do.NULL=FALSE, prefix="kappa2.")
    rownames(B.prec) <- rownames(B.prec, do.NULL=FALSE, prefix="prec.")
    rownames(B.range) <- rownames(B.range, do.NULL=FALSE, prefix="range.")
    BLC = rbind(B.theta, B.tau, B.kappa2, B.prec, B.range)

    fem =
        inla.fmesher.smorg(mesh$loc,
                           mesh$graph$tv,
                           fem=alpha,
                           output=list("c0", "g1", "g2"))

    ## TODO: Adjust fem matrices and B* matrices for fractional alpha:s.
    ## TODO: Construct priors.

    spde =
        inla.spde.generic2(M0=fem$c0, M1=fem$g1, M2=fem$g2,
                           B0=B.tau, B1=B.kappa2, B2=1,
                           theta.mu = rep(0,n.theta),
                           theta.Q = diag(1,n.theta)*0.1,
                           transform = "identity",
                           BLC = BLC)
    spde$model = "matern"

    spde$B =
        list(B.tau = B.tau, B.kappa2 = B.kappa2,
             B.prec = B.prec, B.range = B.range, BLC = BLC)

    return(spde)

    mesh.range = (max(c(diff(range(mesh$loc[,1])),
                        diff(range(mesh$loc[,2])),
                        diff(range(mesh$loc[,3]))
                        )))

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
    } else {
        stop(paste("Model '", model, "' unknown or not implemented.", sep=""))
    }

    return(invisible(spde))
}




inla.extract.el = function(M, ...)
{
    if (is.null(M))
        return(NULL)
    UseMethod("inla.extract.el", M)
}

inla.regex.match =  function(x, match) {
    return(strsplit(x, match)[[1]][1]=="")
}

inla.extract.el.matrix = function(M, match, by.row=TRUE)
{
    if (by.row) {
        return(M[sapply(rownames(M), inla.regex.match, match=match),,drop=FALSE])
    } else {
        return(M[,sapply(colnames(M), inla.regex.match, match=match),drop=FALSE])
    }
}

inla.extract.el.data.frame = function(M, match, by.row=TRUE)
{
    if (by.row) {
        return(M[sapply(rownames(M), inla.regex.match, match=match),,drop=FALSE])
    } else {
        return(M[,sapply(colnames(M), inla.regex.match, match=match),drop=FALSE])
    }
}

inla.extract.el.list = function(M, match)
{
    return(M[sapply(names(M), inla.regex.match, match=match)])
}


inla.spde.result = function(inla, name, spde, do.transform=TRUE, ...)
{
    warning("'inla.spde.result' is not fully implemented yet.")
    inla.require.inherits(inla, "inla", "'inla'")
    inla.require.inherits(spde, "inla.spde", "'spde'")

    result = list()
    if (spde$f$model=="spde") {
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
                            paste("^tau´.[^ ]+$",
                                  sep=""))
        result$summary.log.tau =
            inla.extract.el(inla$summary.random[["UserFunction1"]],
                            paste("^kappa2´.[^ ]+$",
                                  sep=""))

        ## Marginals for log-tau/kappa2
        if (!is.null(inla$marginals.random[["UserFunction1"]])) {
            result$marginals.log.tau =
                inla.extract.el(inla$marginals.random[["UserFunction1"]],
                                paste("^tau´.[^ ]+$",
                                      sep=""))
            result$marginals.log.tau =
                inla.extract.el(inla$marginals.random[["UserFunction1"]],
                                paste("^kappa2´.[^ ]+$",
                                      sep=""))
        }

    } else if (spde$f$model=="spde2") {
        ## Setup rownames
        BLC.names = rownames(spde$internal$param.inla$BLC)

        ## Values
        result$summary.values = inla$summary.random[[name]]

        ## Marginals for values
        if (!is.null(inla$marginals.random[[name]])) {
            result$marginals.values = inla$marginals.random[[name]]
        }

        ## Theta
        result$summary.hyperpar =
            inla.extract.el(inla$summary.hyperpar,
                            paste("Theta[^ ]+ for ", name, "$", sep=""))

        ## Marginals for theta
        if (!is.null(inla$marginals$gyperpar)) {
            result$marginals.hyperpar =
                inla.extract.el(inla$marginals.hyperpar,
                                paste("Theta[^ ]+ for ", name, "$", sep=""))
        }

        ## log-tau/kappa2/prec/range
        if (!is.null(inla$summary.spde2.blc[[name]])) {
            rownames(inla$summary.spde2.blc[[name]]) <- BLC.names

            result$summary.theta =
                inla.extract.el(inla$summary.spde2.blc[[name]],
                                "theta\\.[^ ]+$")
            result$summary.log.tau =
                inla.extract.el(inla$summary.spde2.blc[[name]],
                                "tau\\.[^ ]+$")
            result$summary.log.kappa2 =
                inla.extract.el(inla$summary.spde2.blc[[name]],
                                "kappa2\\.[^ ]+$")
            result$summary.log.prec =
                inla.extract.el(inla$summary.spde2.blc[[name]],
                                "prec\\.[^ ]+$")
            result$summary.log.range =
                inla.extract.el(inla$summary.spde2.blc[[name]],
                                "range\\.[^ ]+$")
        }

        ## Marginals for log-tau/kappa2/prec
        if (!is.null(inla$marginals.spde2.blc[[name]])) {
            names(inla$marginals.spde2.blc[[name]]) <- BLC.names

            result$marginals.theta =
                inla.extract.el(inla$marginals.spde2.blc[[name]],
                                "theta\\.[^ ]+$")
            result$marginals.log.tau =
                inla.extract.el(inla$marginals.spde2.blc[[name]],
                                "tau\\.[^ ]+$")
            result$marginals.log.kappa2 =
                inla.extract.el(inla$marginals.spde2.blc[[name]],
                                "kappa2\\.[^ ]+$")
            result$marginals.log.prec =
                inla.extract.el(inla$marginals.spde2.blc[[name]],
                                "prec\\.[^ ]+$")
            result$marginals.log.range =
                inla.extract.el(inla$marginals.spde2.blc[[name]],
                                "range\\.[^ ]+$")

            if (do.transform) {
                result$marginals.tau =
                    lapply(result$marginals.log.tau,
                           function(x) inla.tmarginal(function(y) exp(y), x))
                result$marginals.kappa =
                    lapply(result$marginals.log.kappa,
                           function(x) inla.tmarginal(function(y) exp(y/2), x))
                result$marginals.prec =
                    lapply(result$marginals.log.prec,
                           function(x) inla.tmarginal(function(y) exp(y), x))
                result$marginals.range =
                    lapply(result$marginals.log.range,
                           function(x) inla.tmarginal(function(y) exp(y), x))
            }
        }
    } else {
        stop(paste("inla.spde model type '", spde$model,
                   "' is not supported by inla.spde.result", sep=""))
    }

    return(result)
}









##inla.spde.create(mesh, model=list("matern"), ...)
##inla.spde.create(mesh, model=list("heat", Qw=..., t=...), ...)
##inla.spde.create(mesh, model=list("imatern"), ...)
