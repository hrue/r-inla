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
    class(spde) = c("inla.spde2", "inla.spde", "inla.model.class")

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
            B0[,c(TRUE, fixed), drop=FALSE] %*% c(1.0, theta.fixed)
        param.inla$B1[,1] =
            B1[,c(TRUE, fixed), drop=FALSE] %*% c(1.0, theta.fixed)
        param.inla$B2[,1] =
            B2[,c(TRUE, fixed), drop=FALSE] %*% c(1.0, theta.fixed)
        param.inla$BLC[,1] =
            BLC[,c(TRUE, fixed), drop=FALSE] %*% c(1.0, theta.fixed)
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
    for (k in 1:n.theta) {
        eval(parse(text=
                   paste("spde$f$hyper.default$theta", k,
                         "$initial = param.inla$theta.mu[k]", sep="")))
    }

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
             prior.variance.nominal = 1,
             prior.range.nominal = NULL,
             prior.tau = NULL,
             prior.kappa = NULL,
             theta.prior.mean = NULL,
             theta.prior.prec = 0.1,
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
    alpha.nominal = max(1.5, alpha)
    nu.nominal = max(0.5, nu)

    n.spde = mesh$n
    n.theta = ncol(B.kappa)-1L

    B.kappa = inla.spde.homogenise_B_matrix(B.kappa, n.spde, n.theta)
    B.tau = inla.spde.homogenise_B_matrix(B.tau, n.spde, n.theta)
    B.prec =
        cbind(+ lgamma(alpha.nominal) - lgamma(nu.nominal)
              + (d/2)*log(4*pi) + 2*nu.nominal*B.kappa[,1]
              + 2*B.tau[,1],
              2*nu.nominal*B.kappa[,-1,drop=FALSE]
              + 2*B.tau[,-1,drop=FALSE] )
    if (is.stationary) {
        B.tau = B.tau[1,,drop=FALSE]
        B.kappa = B.kappa[1,,drop=FALSE]
        B.prec = B.prec[1,,drop=FALSE]
    }
    B.variance = -B.prec
    B.range =
        cbind(0.5*log(8*nu.nominal)-B.kappa[,1],
              -B.kappa[,-1,drop=FALSE])

    B.theta = cbind(0,diag(1, n.theta))
    rownames(B.theta) <- rownames(B.theta, do.NULL=FALSE, prefix="theta.")
    rownames(B.tau) <- rownames(B.tau, do.NULL=FALSE, prefix="tau.")
    rownames(B.kappa) <- rownames(B.kappa, do.NULL=FALSE, prefix="kappa.")
    rownames(B.variance) <-
        rownames(B.variance, do.NULL=FALSE, prefix="variance.nominal.")
    rownames(B.range) <-
        rownames(B.range, do.NULL=FALSE, prefix="range.nominal.")
    BLC = rbind(B.theta, B.tau, B.kappa, B.variance, B.range)

    fem =
        inla.fmesher.smorg(mesh$loc,
                           mesh$graph$tv,
                           fem=2,
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
        stop(paste("Unsupported alpha value (", alpha,
                   "). Supported values are 0<alpha <= 2", sep=""))
    }

    ## Construct priors.
    if (is.null(theta.prior.prec)) {
        theta.prior.prec = diag(0.1, n.theta, n.theta)
    } else {
        theta.prior.prec = as.matrix(theta.prior.prec)
        if (ncol(theta.prior.prec) == 1) {
            theta.prior.prec =
                diag(as.vector(theta.prior.prec), n.theta, n.theta)
        }
        if ((nrow(theta.prior.prec) != n.theta) ||
            (ncol(theta.prior.prec) != n.theta)) {
            stop(paste("Size of theta.prior.prec is (",
                       paste(dim(theta.prior.prec), collapse=",", sep=""),
                       ") but should be (",
                       paste(c(n.theta, n.theta), collapse=",", sep=""),
                       ")."))
        }
    }

    if (is.null(theta.prior.mean)) {
        if (is.null(prior.range.nominal)) {
            mesh.range = (max(c(diff(range(mesh$loc[,1])),
                                diff(range(mesh$loc[,2])),
                                diff(range(mesh$loc[,3]))
                                )))
            prior.range.nominal = mesh.range*0.2
        }

        if (is.null(prior.kappa)) {
            prior.kappa = sqrt(8*nu.nominal)/prior.range.nominal
        }
        if (is.null(prior.tau)) {
            prior.tau =
                sqrt(gamma(nu.nominal)/gamma(alpha.nominal)/
                     (4*pi*prior.kappa^(2*nu.nominal)*prior.variance.nominal))
        }

        if (n.theta>0) {
            theta.prior.mean =
                solve(rbind(B.tau[,-1,drop=FALSE], B.kappa[,-1,drop=FALSE]),
                      c(log(prior.tau) - B.tau[,1],
                        log(prior.kappa) - B.kappa[,1]))
        } else {
            theta.prior.mean = rep(0, n.theta) ## Empty vector
        }
    }

    spde =
        inla.spde2.generic(M0=M0, M1=M1, M2=M2,
                           B0=B.phi0, B1=B.phi1, B2=1,
                           theta.mu = theta.prior.mean,
                           theta.Q = theta.prior.prec,
                           transform = "identity",
                           BLC = BLC)
    spde$model = "matern"
    spde$BLC = BLC


    return(invisible(spde))
}


inla.spde2.theta2phi0 = function(spde, theta)
{
    inla.require.inherits(spde, "inla.spde2", "'spde'")

    return(exp(spde$param.inla$B0[,1, drop=TRUE] +
               spde$param.inla$B0[,-1, drop=FALSE] %*% theta))
}

inla.spde2.theta2phi1 = function(spde, theta)
{
    inla.require.inherits(spde, "inla.spde2", "'spde'")

    return(exp(spde$param.inla$B1[,1, drop=TRUE] +
               spde$param.inla$B1[,-1, drop=FALSE] %*% theta))
}

inla.spde2.theta2phi2 = function(spde, theta)
{
    inla.require.inherits(spde, "inla.spde2", "'spde'")

##    warning("TODO: support link functions for phi2")
    return((spde$param.inla$B2[,1, drop=TRUE] +
            spde$param.inla$B2[,-1, drop=FALSE] %*% theta))
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
    BLC.names = rownames(spde$param.inla$BLC)

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

    ## log-tau/kappa/variance/range
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
        result$summary.log.variance.nominal =
            inla.extract.el(inla$summary.spde2.blc[[name]],
                            "variance.nominal\\.[^ ]+$")
        result$summary.log.range.nominal =
            inla.extract.el(inla$summary.spde2.blc[[name]],
                            "range.nominal\\.[^ ]+$")
    }

    ## Marginals for log-tau/kappa/variance/range
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
        result$marginals.log.variance.nominal =
            inla.extract.el(inla$marginals.spde2.blc[[name]],
                            "variance.nominal\\.[^ ]+$")
        result$marginals.log.range.nominal =
            inla.extract.el(inla$marginals.spde2.blc[[name]],
                            "range.nominal\\.[^ ]+$")

        if (do.transform) {
            result$marginals.tau =
                lapply(result$marginals.log.tau,
                       function(x) inla.tmarginal(function(y) exp(y), x))
            result$marginals.kappa =
                lapply(result$marginals.log.kappa,
                       function(x) inla.tmarginal(function(y) exp(y), x))
            result$marginals.variance.nominal =
                lapply(result$marginals.log.variance.nominal,
                       function(x) inla.tmarginal(function(y) exp(y), x))
            result$marginals.range.nominal =
                lapply(result$marginals.log.range.nominal,
                       function(x) inla.tmarginal(function(y) exp(y), x))
        }
    }

    return(result)
}



inla.spde2.models = function()
{
    return(c("generic", "matern"))
}


## spde.common-connections:
inla.spde.precision.inla.spde2 = inla.spde2.precision
inla.spde.result.inla.spde2 = inla.spde2.result


##inla.spde.create(mesh, model=list("matern"), ...)
##inla.spde.create(mesh, model=list("heat", Qw=..., t=...), ...)
##inla.spde.create(mesh, model=list("imatern"), ...)


## Deprecated:
##inla.spde.generic2 = inla.spde2.generic
