## 'spde2' model functions






inla.spde2.generic =
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

    spde = (list(model = "generic",
                 n.spde = n.spde,
                 n.theta = n.theta,
                 param.inla = list(),
                 f = list()
                 ))
    class(spde) = c("inla.spde2", "inla.model.class")

    B0 = inla.spde.homogenise_B_matrix(B0, n.spde, n.theta)
    B1 = inla.spde.homogenise_B_matrix(B1, n.spde, n.theta)
    B2 = inla.spde.homogenise_B_matrix(B2, n.spde, n.theta)
    BLC = inla.spde.homogenise_B_matrix(BLC, nrow(BLC), n.theta)

    param.inla =
        list(n = n.spde,
             n.theta = n.theta,
             M0=M0, M1=M1, M2=M2,
             B0=B0, B1=B1, B2=B2, BLC=BLC,
             theta.mu=theta.mu, theta.Q=theta.Q,
             transform=transform,
             theta.initial=theta.initial,
             fixed=fixed,
             theta.fixed=theta.fixed)
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

    spde$param.inla = param.inla

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
    fmesher.write(spde$param.inla$M0, spde.prefix, "M0")
    fmesher.write(spde$param.inla$M1, spde.prefix, "M1")
    fmesher.write(spde$param.inla$M2, spde.prefix, "M2")
    fmesher.write(spde$param.inla$B0, spde.prefix, "B0")
    fmesher.write(spde$param.inla$B1, spde.prefix, "B1")
    fmesher.write(spde$param.inla$B2, spde.prefix, "B2")
    fmesher.write(spde$param.inla$BLC, spde.prefix, "BLC")

    return(spde)
}


inla.spde2.matern =
    function(mesh,
             alpha=2,
             B.tau = matrix(c(0,1,0),1,3),
             B.kappa = matrix(c(0,0,1),1,3),
             sigma2.approx = 1,
             theta.prior.mean = NULL,
             theta.prior.prec = NULL,
             fractional.method = c("parsimonious", "null"))
{
    inla.require.inherits(mesh, "inla.mesh", "'mesh'")
    fractional.method = match.arg(fractional.method)
    if (is.null(B.tau))
        stop("B.tau must not be NULL.")
    if (is.null(B.kappa))
        stop("B.kappa must not be NULL.")
    is.stationary = (nrow(B.kappa)==1) && (nrow(B.tau)==1)

    d = 2
    nu = alpha-d/2

    n.spde = mesh$n
    n.theta = ncol(B.kappa)-1L

    B.kappa = inla.spde.homogenise_B_matrix(B.kappa, n.spde, n.theta)
    B.tau = inla.spde.homogenise_B_matrix(B.tau, n.spde, n.theta)
    B.prec =
        cbind(+ lgamma(max(1.5,alpha)) - lgamma(max(0.5,nu))
              + (d/2)*log(4*pi) + 2*max(0.5,nu)*B.kappa[,1]
              + 2*B.tau[,1],
              2*max(0.5,nu)*B.kappa[,-1,drop=FALSE]
              + 2*B.tau[,-1,drop=FALSE] )
    if (is.stationary) {
        B.tau = B.tau[1,,drop=FALSE]
        B.kappa = B.kappa[1,,drop=FALSE]
        B.prec = B.prec[1,,drop=FALSE]
    }
    B.range =
        cbind(0.5*log(8*max(0.5,nu))-B.kappa[,1],
              -B.kappa[,-1,drop=FALSE])

    B.theta = cbind(0,diag(1, n.theta))
    rownames(B.theta) <- rownames(B.theta, do.NULL=FALSE, prefix="theta.")
    rownames(B.tau) <- rownames(B.tau, do.NULL=FALSE, prefix="tau.")
    rownames(B.kappa) <- rownames(B.kappa, do.NULL=FALSE, prefix="kappa2.")
    rownames(B.prec) <- rownames(B.prec, do.NULL=FALSE, prefix="prec.")
    rownames(B.range) <- rownames(B.range, do.NULL=FALSE, prefix="range.")
    BLC = rbind(B.theta, B.tau, B.kappa, B.prec, B.range)

    fem =
        inla.fmesher.smorg(mesh$loc,
                           mesh$graph$tv,
                           fem=alpha,
                           output=list("c0", "g1", "g2"))

    if (alpha==2) {
        B.phi0 = B.tau
        B.phi1 = 2*B.kappa
        M0 = fem$c0
        M1 = fem$g1
        M2 = fem$g2
    } else if (alpha==1) {
        B.phi0 = B.tau
        B.phi1 = B.kappa
        M0 = fem$c0
        M1 = fem$g1*0
        M2 = fem$g1
    } else if (!is.stationary) {
        stop("Non-stationary Matern with fractional alpha is not implemented.")
    } else if (alpha>1) {
        if (fractional.method == "parsimonious") {
            lambda = alpha-floor(alpha)
            b = matrix(c(1,0,0, 1,1,0, 1,2,1),3,3) %*%
                solve(matrix(1/(c(4:2, 3:1, 2:0)+lambda), 3, 3),
                      1/(c(4:2)+lambda-alpha))
        } else if (fractional.method == "null") {
            b = c(1,alpha,alpha*(alpha-1)/2)
        } else {
            stop(paste("Unknown fractional.method '", fractional.method,
                       "'.", sep=""))
        }
        B.phi0 = B.tau + (alpha-2)*B.kappa
        B.phi1 = 2*B.kappa
        M0 = fem$c0*b[1]
        M1 = fem$g1*b[2]/2
        M2 = fem$g2*b[3]
    } else if (alpha>0) {
        if (fractional.method == "parsimonious") {
            lambda = alpha-floor(alpha)
            b = matrix(c(1,0,1,1),2,2) %*%
                solve(matrix(1/(c(2:1, 1:0)+lambda), 2, 2),
                      1/(c(2:1)+lambda-alpha))
        } else if (fractional.method == "null") {
            b = c(1,alpha)
        } else {
            stop(paste("Unknown fractional.method '", fractional.method,
                       "'.", sep=""))
        }
        B.phi0 = B.tau + (alpha-1)*B.kappa
        B.phi1 = B.kappa
        M0 = fem$c0*b[1]
        M1 = fem$g1*0
        M2 = fem$g1*b[2]
    } else {
        stop("Unsupported alpha value")
    }

    ## TODO: Construct priors.

    spde =
        inla.spde2.generic(M0=M0, M1=M1, M2=M2,
                           B0=B.phi0, B1=B.phi1, B2=1,
                           theta.mu = rep(0,n.theta),
                           theta.Q = diag(1,n.theta)*0.1,
                           transform = "identity",
                           BLC = BLC)
    spde$model = "matern"
    spde$BLC = BLC


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


inla.spde2.theta2phi0 = function(spde, theta)
{
    inla.require.inherits(spde, "inla.spde2", "'spde'")

    return(exp(spde$param.inla$B0[,1] +
               spde$param.inla$B0[,-1] %*% theta))
}

inla.spde2.theta2phi1 = function(spde, theta)
{
    inla.require.inherits(spde, "inla.spde2", "'spde'")

    return(exp(spde$param.inla$B1[,1] +
               spde$param.inla$B1[,-1] %*% theta))
}

inla.spde2.theta2phi2 = function(spde, theta)
{
    inla.require.inherits(spde, "inla.spde2", "'spde'")

##    warning("TODO: support link functions for phi2")
    return((spde$param.inla$B2[,1] +
            spde$param.inla$B2[,-1] %*% theta))
}

inla.spde2.precision =
    function(spde, theta=NULL,
             phi0=inla.spde2.theta2phi0(spde, theta),
             phi1=inla.spde2.theta2phi1(spde, theta),
             phi2=inla.spde2.theta2phi2(spde, theta))
{
    inla.require.inherits(spde, "inla.spde2", "'spde'")
    if (spde$f$model!="spde2") {
        stop("'inla.spde2.precision' only supports internal inla models 'spde2'")
    }

    D0 = Diagonal(spde$n.spde, phi0)
    D1 = Diagonal(spde$n.spde, phi1)
    D12 = Diagonal(spde$n.spde, phi1*phi2)
    Q = (D0 %*% (D1 %*% spde$param.inla$M0 %*% D1 +
                 D12 %*% spde$param.inla$M1 +
                 t(spde$param.inla$M1) %*% D12 +
                 spde$param.inla$M2) %*% D0)

    return(Q)
}



inla.spde2.result = function(inla, name, spde, do.transform=TRUE, ...)
{
    inla.require.inherits(inla, "inla", "'inla'")
    inla.require.inherits(spde, "inla.spde2", "'spde'")
    if (!spde$f$model=="spde2") {
        stop("'inla.spde2.result' only supports internal inla models 'spde2'")
    }

    result = list()
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

    ## log-tau/kappa/prec/range
    if (!is.null(inla$summary.spde2.blc[[name]])) {
        rownames(inla$summary.spde2.blc[[name]]) <- BLC.names

        result$summary.theta =
            inla.extract.el(inla$summary.spde2.blc[[name]],
                            "theta\\.[^ ]+$")
        result$summary.log.tau =
            inla.extract.el(inla$summary.spde2.blc[[name]],
                            "tau\\.[^ ]+$")
        result$summary.log.kappa =
            inla.extract.el(inla$summary.spde2.blc[[name]],
                            "kappa\\.[^ ]+$")
        result$summary.log.prec =
            inla.extract.el(inla$summary.spde2.blc[[name]],
                            "prec\\.[^ ]+$")
        result$summary.log.range =
            inla.extract.el(inla$summary.spde2.blc[[name]],
                            "range\\.[^ ]+$")
    }

    ## Marginals for log-tau/kappa/prec
    if (!is.null(inla$marginals.spde2.blc[[name]])) {
        names(inla$marginals.spde2.blc[[name]]) <- BLC.names

        result$marginals.theta =
            inla.extract.el(inla$marginals.spde2.blc[[name]],
                            "theta\\.[^ ]+$")
        result$marginals.log.tau =
            inla.extract.el(inla$marginals.spde2.blc[[name]],
                            "tau\\.[^ ]+$")
        result$marginals.log.kappa =
            inla.extract.el(inla$marginals.spde2.blc[[name]],
                            "kappa\\.[^ ]+$")
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

    return(result)
}



inla.spde2.models = function()
{
    return(c("generic","matern"))
}


## spde.common-connections:
inla.spde.precision.inla.spde2 = inla.spde2.precision
inla.spde.result.inla.spde2 = inla.spde2.result


##inla.spde.create(mesh, model=list("matern"), ...)
##inla.spde.create(mesh, model=list("heat", Qw=..., t=...), ...)
##inla.spde.create(mesh, model=list("imatern"), ...)


## Deprecated:
##inla.spde.generic2 = inla.spde2.generic
