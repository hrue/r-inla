#' @include spde.common.R
NULL

#' Generic spde2 model creation.
#' 
#' Creates and inla.spde2 object describing the internal structure of an
#' `'spde2'` model.
#' 
#' 
#' @aliases inla.spde2.generic inla.spde2
#' @param M0 The symmetric `M0` matrix.
#' @param M1 The square `M1` matrix.
#' @param M2 The symmetric `M2` matrix.
#' @param B0 Basis definition matrix for \eqn{\phi_0}{phi0}.
#' @param B1 Basis definition matrix for \eqn{\phi_2}{phi1}.
#' @param B2 Basis definition matrix for \eqn{\phi_2}{phi2}.
#' @param theta.mu Prior expectation for the \eqn{\theta}{theta} vector
#' @param theta.Q Prior precision for the \eqn{\theta}{theta} vector
#' @param transform Transformation link for \eqn{\phi_2}{phi2}.  Valid settings
#' are `"logit"`, `"log"`, and `"identity"`
#' @param theta.initial Initial value for the \eqn{\theta}{theta} vector.
#' Default `theta.mu`
#' @param fixed Logical vector.  For every `TRUE` value, treat the
#' corresponding `theta` value as known.
#' @param theta.fixed Vector holding the values of fixed `theta` values.
#' Default `=theta.initial[fixed]`
#' @param BLC Basis definition matrix for linear combinations of `theta`.
#' @param \dots Additional parameters, currently unused.
#' @param theta parameter values to be mapped.
#' @return For `inla.spde2.generic`, an [inla.spde2()] object.
#' 
#' For `inla.spde2.theta2phi0/1/2`, a vector of \eqn{\phi}{phi} values.
#' @author Finn Lindgren \email{finn.lindgren@@gmail.com}
#' @seealso [inla.spde2.models()], [inla.spde2.matern()]
#' @export inla.spde2.generic
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

    if (!inla.is.matrix(M0)) stop("M0 must be of a matrix type.")
    if (!inla.is.matrix(M1)) stop("M1 must be of a matrix type.")
    if (!inla.is.matrix(M2)) stop("M2 must be of a matrix type.")
    M0 = inla.as.dgTMatrix(M0)
    M1 = inla.as.dgTMatrix(M1)
    M2 = inla.as.dgTMatrix(M2)

    n.spde = nrow(M0)
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

    ## Remove all BLC rows that would produce point masses, which is
    ## not supported by inla:
    if (n.theta==0) {
        param.inla$BLC = matrix(0,0,1)
    } else {
        param.inla$BLC =
            param.inla$BLC[rowSums(abs(param.inla$BLC[,-1,drop=FALSE]))>0,,
                           drop=FALSE]
    }

    spde$param.inla = param.inla

    spde$f = (list(model="spde2",
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
    if (n.theta>0)
        for (k in 1:n.theta) {
            eval(parse(text=
                       paste("spde$f$hyper.default$theta", k,
                             "$initial = param.inla$theta.mu[k]", sep="")))
        }

    return(spde)
}









inla.internal.spde2.matern.B.tau =
    function(precision.fcn, n.theta, delta = 1e-3, ...)
{
    theta.0 = rep(0, n.theta)
    Q.0 = precision.fcn(theta.0, ...)
    logv.0 = log(diag(inla.qinv(Q.0)))

    B.tau = matrix(0, nrow(Q.0), 1L+n.theta)
    for (k in seq_len(n.theta)) {
        theta = theta.0
        theta[k] = theta[k] + delta
        Q = precision.fcn(theta, ...)
        logv = log(diag(inla.qinv(Q)))

        dlogv = (logv-logv.0)/delta

        B.tau[, k+1] = dlogv/2
    }
    B.tau[,1L] = logv.0/2

    return(B.tau)
}





#' Approximate variance-compensating basis functions
#' 
#' Calculates an approximate basis for `tau` and `kappa` for an
#' `inla.spde2.matern` model where `tau` is a rescaling parameter.
#' 
#' 
#' @param mesh An [inla.mesh()] object.
#' @param B.sd Desired basis for log-standard deviations.
#' @param B.range Desired basis for spatial range.
#' @param method Construction method selector.  Expert option only.
#' @param local.offset.compensation If `FALSE`, only compensate in the
#' average for the `tau` offset.
#' @param alpha The model `alpha` parameter.
#' @param \dots Additional parameters passed on to internal
#' `inla.spde2.matern` calls.
#' @return List of basis specifications \item{B.tau }{Basis for
#' `log(tau)`} \item{B.kappa }{Basis for `log(kappa)`} Intended for
#' passing on to [inla.spde2.matern()].
#' @author Finn Lindgren \email{finn.lindgren@@gmail.com}
#' @seealso [inla.spde2.matern()]
#' @export inla.spde2.matern.sd.basis
inla.spde2.matern.sd.basis =
    function(mesh, B.sd, B.range, method=1,
             local.offset.compensation=FALSE,
             alpha=2, ...)
{
    p.sd = ncol(B.sd)-1L
    p.kappa = ncol(B.range)-1L

    d = inla.ifelse(inherits(mesh, "inla.mesh"), 2, 1)
    nu = alpha-d/2
    nu.nominal = max(0.5, nu)
    alpha.nominal = max(nu.nominal+d/2, alpha)

    n.spde = inla.ifelse(d==2, mesh$n, mesh$m)

    ## log(kappa) = 0.5*log(8*nu) - log(range)
    B.kappa = (cbind(log(8*nu.nominal)/2 - B.range[,1],
                     -B.range[,-1,drop=FALSE] ))

    if (method==1) {
        spde =
            inla.spde2.matern(mesh, B.tau=as.matrix(0), B.kappa=B.kappa, ...)

        B.tau =
            inla.internal.spde2.matern.B.tau(
                precision.fcn=
                function(theta, spde) {
                    return(inla.spde2.precision(spde, theta=theta))
                },
                n.theta = p.kappa,
                spde = spde)
        if (!local.offset.compensation) {
            B.tau[,1] = rep(mean(B.tau[,1]), spde$n.spde)
        }
    } else {
        spde =
            inla.spde2.matern(mesh,
                              B.tau=as.matrix(0),
                              B.kappa=B.kappa[,1,drop=FALSE],
                              alpha=alpha, ...)
        lvar.0 =
            log(diag(inla.qinv(inla.spde2.precision(spde, theta=rep(0,0)))))
        if (!local.offset.compensation) {
            lvar.0 = rep(mean(lvar.0), spde$n.spde)
        }

        B.tau = -B.kappa*nu
        B.tau[,1] = lvar.0/2
    }

    log.var.scaling = lgamma(nu.nominal)-lgamma(alpha.nominal)-log(4*pi)*d/2

    B.tau =
        cbind(B.tau[,1]-B.sd[,1],
              -B.sd[,-1,drop=FALSE],
              B.tau[,-1,drop=FALSE])
    B.kappa =
        cbind(B.kappa[,1],
              matrix(0, n.spde, p.sd),
              B.kappa[,-1,drop=FALSE])

    return(list(B.tau=B.tau, B.kappa=B.kappa))

}



inla.internal.test.spde2.sd.basis = function (k=1, dth=0.1, r=1000, globe=25, compensate=FALSE) {

    alpha=2
    mesh = inla.mesh.create(globe=globe)

    B.common = inla.mesh.basis(mesh, "b.spline", n=3)
    B.common = cbind(B.common, B.common[,1]*(mesh$loc[,3]>0))
    B.common[,1] = B.common[,1]*(mesh$loc[,3]<=0)

    B.sd = cbind(0, B.common)
    B.range = cbind(log(r/6370), B.common)

    sd.basis0 = inla.spde2.matern.sd.basis(mesh, B.sd=B.sd, B.range=B.range, method=0, alpha=alpha, local.offset.compensation=compensate)
    sd.basis1 = inla.spde2.matern.sd.basis(mesh, B.sd=B.sd, B.range=B.range, method=1, alpha=alpha, local.offset.compensation=compensate)

    spde0 =
        inla.spde2.matern(mesh, alpha=alpha,
                          B.tau=sd.basis0$B.tau,
                          B.kappa=sd.basis0$B.kappa)
    spde1 =
        inla.spde2.matern(mesh, alpha=alpha,
                          B.tau=sd.basis1$B.tau,
                          B.kappa=sd.basis1$B.kappa)

    th0 = c(0,0,0,0,0,0,0,0)
    th1 = th0
    th1[k] = th0[k]+dth
    v0.0=diag(inla.qinv(inla.spde2.precision(spde0,theta=th0)))
    v0.1=diag(inla.qinv(inla.spde2.precision(spde0,theta=th1)))
    v1.0=diag(inla.qinv(inla.spde2.precision(spde1,theta=th0)))
    v1.1=diag(inla.qinv(inla.spde2.precision(spde1,theta=th1)))

    op=par(mfrow=c(2,2))
    on.exit(par(op))
    plot(asin(mesh$loc[,3]),rowSums(sd.basis0$B.kappa[,1+k,drop=FALSE]), pch=20,
         ylim=(range(c(0, range(rowSums(sd.basis0$B.kappa[,1+k,drop=FALSE])),
                       range(-rowSums(sd.basis0$B.tau[,1+k,drop=FALSE])),
                       range(-rowSums(sd.basis1$B.tau[,1+k,drop=FALSE]))))),
         main="Basis functions")
    points(asin(mesh$loc[,3]),-rowSums(sd.basis0$B.tau[,1+k,drop=FALSE]),col=2,pch=20)
    points(asin(mesh$loc[,3]),-rowSums(sd.basis1$B.tau[,1+k,drop=FALSE]),col=4,pch=20)
    legend("bottomright",
           legend=paste("Method", 0:1),
           col=c(2,4), pch=20)

    plot(asin(mesh$loc[,3]),v0.1, col=2,pch=20,
         ylim=range(c(0,range(v0.0),range(v1.0),range(v0.1),range(v1.1))),
         main="Absolute variances")
    points(asin(mesh$loc[,3]),v1.1, col=4,pch=20)
    points(asin(mesh$loc[,3]),v0.0, col=2,pch=1)
    points(asin(mesh$loc[,3]),v1.0, col=4,pch=1)
    legend("bottomright",
           legend=paste("Method", 0:1),
           col=c(2,4), pch=20)

    plot(asin(mesh$loc[,3]),v0.1-v0.0, col=2,pch=20,
         ylim=range(c(0,range(v0.1-v0.0),range(v1.1-v1.0))),
         main="Absolute variance error")
    points(asin(mesh$loc[,3]),v1.1-v1.0, col=4,pch=20)
    legend("bottomright",
           legend=paste("Method", 0:1),
           col=c(2,4), pch=20)

    plot(asin(mesh$loc[,3]),v0.1/v0.0-1, col=2,pch=20,
         ylim=range(c(0,range(v0.1/v0.0-1),range(v1.1/v1.0-1))),
         main="Relative variance error")
    points(asin(mesh$loc[,3]),v1.1/v1.0-1, col=4,pch=20)
    legend("bottomright",
           legend=paste("Method", 0:1),
           col=c(2,4), pch=20)

    return(environment())
}






#' Parameter settings for `inla.spde2.matern` models.
#' 
#' Construct parameter settings for `inla.spde2.matern` models.
#' 
#' 
#' @param mesh The mesh to build the model on, as an [inla.mesh()]
#' object.
#' @param alpha Fractional operator order, \eqn{0<\alpha\leq 2} supported.
#' (\eqn{\nu=\alpha-d/2})
#' @param B.tau Matrix with specification of log-linear model for \eqn{\tau}.
#' @param B.kappa Matrix with specification of log-linear model for
#' \eqn{\kappa}.
#' @param prior.variance.nominal Nominal prior mean for the field variance
#' @param prior.range.nominal Nominal prior mean for the spatial range
#' @param prior.tau Prior mean for tau (overrides
#' `prior.variance.nominal`)
#' @param prior.kappa Prior mean for kappa (overrides
#' `prior.range.nominal`)
#' @param theta.prior.mean (overrides `prior.*`)
#' @param theta.prior.prec Scalar, vector or matrix, specifying the joint prior
#' precision for \eqn{theta}.
#' @author Finn Lindgren \email{finn.lindgren@@gmail.com}
#' @seealso [inla.spde2.matern()]
#' @export param2.matern.orig
param2.matern.orig =
    function(mesh,
             alpha=2,
             B.tau = matrix(c(0,1,0),1,3),
             B.kappa = matrix(c(0,0,1),1,3),
             prior.variance.nominal = 1,
             prior.range.nominal = NULL,
             prior.tau = NULL,
             prior.kappa = NULL,
             theta.prior.mean = NULL,
             theta.prior.prec = 0.1)
{
    ## NOTE: For d==1, degree==2, the B.* must be given per basis
    ## function, not per knot.

    inla.require.inherits(mesh, c("inla.mesh", "inla.mesh.1d"), "'mesh'")
    if (is.null(B.tau))
        stop("B.tau must not be NULL.")
    if (is.null(B.kappa))
        stop("B.kappa must not be NULL.")
    is.stationary = (nrow(B.kappa)==1) && (nrow(B.tau)==1)

    d = inla.ifelse(inherits(mesh, "inla.mesh"), 2, 1)
    nu = alpha-d/2
    nu.nominal = max(0.5, nu)
    alpha.nominal = max(nu.nominal+d/2, alpha)

    n.spde = inla.ifelse(d==2, mesh$n, mesh$m)
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

    if (n.theta>0) {
        B.theta = cbind(0,diag(1, n.theta))
        rownames(B.theta) <- rownames(B.theta, do.NULL=FALSE, prefix="theta.")
    } else {
        B.theta = NULL
    }
    rownames(B.tau) <- rownames(B.tau, do.NULL=FALSE, prefix="tau.")
    rownames(B.kappa) <- rownames(B.kappa, do.NULL=FALSE, prefix="kappa.")
    rownames(B.variance) <-
        rownames(B.variance, do.NULL=FALSE, prefix="variance.nominal.")
    rownames(B.range) <-
        rownames(B.range, do.NULL=FALSE, prefix="range.nominal.")
    BLC = rbind(B.theta, B.tau, B.kappa, B.variance, B.range)


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
            mesh.range =
                inla.ifelse(d==2,
                            (max(c(diff(range(mesh$loc[,1])),
                                   diff(range(mesh$loc[,2])),
                                   diff(range(mesh$loc[,3]))
                                   ))),
                            diff(mesh$interval))
            prior.range.nominal = mesh.range*0.2
        } else {
          if (!is.numeric(prior.range.nominal) ||
              (length(prior.range.nominal) != 1)) {
            stop(paste0("'prior.range.nominal' must be NULL or a single scalar value.\n",
                        "Did you intend to supply 'prior.range' to inla.spde2.pcmatern instead?"))
          }
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
                qr.solve(rbind(B.tau[,-1,drop=FALSE], B.kappa[,-1,drop=FALSE]),
                         c(log(prior.tau) - B.tau[,1],
                           log(prior.kappa) - B.kappa[,1]))
        } else {
            theta.prior.mean = rep(0, n.theta) ## Empty vector
        }
    }

    param =
        list(is.stationary=is.stationary,
             B.tau=B.tau, B.kappa=B.kappa, BLC=BLC,
             theta.prior.mean=theta.prior.mean,
             theta.prior.prec=theta.prior.prec)
    return(param)
}



#' Matern SPDE model object for INLA
#' 
#' Create an `inla.spde2` model object for a Matern model. Use
#' `inla.spde2.pcmatern` instead for a PC prior for the parameters.
#' 
#' This method constructs a Matern SPDE model, with spatial scale parameter
#' \eqn{\kappa(u)} and variance rescaling parameter \eqn{\tau(u)}.
#' 
#' \deqn{(\kappa^2(u)-\Delta)^{\alpha/2}(\tau(u) }{(kappa^2(u)-Delta)^(alpha/2)
#' (tau(u) x(u)) = W(u)}\deqn{ x(u))=W(u)}{(kappa^2(u)-Delta)^(alpha/2) (tau(u)
#' x(u)) = W(u)}
#' 
#' Stationary models are supported for \eqn{0 < \alpha \leq 2}, with spectral
#' approximation methods used for non-integer \eqn{\alpha}, with approximation
#' method determined by `fractional.method`.
#' 
#' Non-stationary models are supported for \eqn{\alpha=2} only, with \itemize{
#' \item \eqn{\log\tau(u) = B^\tau_0(u) + \sum_{k=1}^p B^\tau_k(u) }{log tau(u)
#' = B.tau_0(u) + sum_{k=1}^p B.tau_k(u) theta_k}\eqn{ \theta_k}{log tau(u) =
#' B.tau_0(u) + sum_{k=1}^p B.tau_k(u) theta_k}
#' 
#' \item \eqn{\log\kappa(u) = B^{\kappa}_0(u) + \sum_{k=1}^p B^{\kappa}_k(u)
#' }{log kappa(u) = B.kappa_0(u) + sum_{k=1}^p B.kappa_k(u) theta_k}\eqn{
#' \theta_k}{log kappa(u) = B.kappa_0(u) + sum_{k=1}^p B.kappa_k(u) theta_k}
#' 
#' }
#' 
#' The same parameterisation is used in the stationary cases, but with
#' \eqn{B^\tau_0}, \eqn{B^\tau_k}, \eqn{B^\kappa_0}, and \eqn{B^\tau_k}
#' constant across \eqn{u}.
#' 
#' Integration and other general linear constraints are supported via the
#' `constr`, `extraconstr.int`, and `extraconstr` parameters,
#' which also interact with `n.iid.group`.
#' 
#' @param mesh The mesh to build the model on, as an [inla.mesh()] or
#' [inla.mesh.1d()] object.
#' @param alpha Fractional operator order, \eqn{0<\alpha\leq 2} supported.
#' (\eqn{\nu=\alpha-d/2})
#' @param param Parameter, e.g. generated by `param2.matern.orig`
#' @param constr If `TRUE`, apply an integrate-to-zero constraint.
#' Default `FALSE`.
#' @param extraconstr.int Field integral constraints.
#' @param extraconstr Direct linear combination constraints on the basis
#' weights.
#' @param fractional.method Specifies the approximation method to use for
#' fractional (non-integer) `alpha` values. `'parsimonious'` gives an
#' overall approximate minimal covariance error, `'null'` uses
#' approximates low-order properties.
#' @param B.tau Matrix with specification of log-linear model for \eqn{\tau}.
#' @param B.kappa Matrix with specification of log-linear model for
#' \eqn{\kappa}.
#' @param prior.variance.nominal Nominal prior mean for the field variance
#' @param prior.range.nominal Nominal prior mean for the spatial range
#' @param prior.tau Prior mean for tau (overrides
#' `prior.variance.nominal`)
#' @param prior.kappa Prior mean for kappa (overrides
#' `prior.range.nominal`)
#' @param theta.prior.mean (overrides `prior.*`)
#' @param theta.prior.prec Scalar, vector or matrix, specifying the joint prior
#' precision for \eqn{theta}.
#' @param n.iid.group If greater than 1, build an explicitly iid replicated
#' model, to support constraints applied to the combined replicates, for
#' example in a time-replicated spatial model. Constraints can either be
#' specified for a single mesh, in which case it's applied to the average of
#' the replicates (`ncol(A)` should be `mesh$n` for 2D meshes,
#' `mesh$m` for 1D), or as general constraints on the collection of
#' replicates (`ncol(A)` should be `mesh$n * n.iid.group` for 2D
#' meshes, `mesh$m * n.iid.group` for 1D).
#' @param \dots Additional parameters for special uses.
#' @return An `inla.spde2` object.
#' @author Finn Lindgren \email{finn.lindgren@@gmail.com}
#' @seealso [inla.mesh.2d()], [inla.mesh.create()],
#' [inla.mesh.1d()], [inla.mesh.basis()],
#' [inla.spde2.pcmatern()], [inla.spde2.generic()]
#' @examples
#' 
#' n = 100
#' field.fcn = function(loc) (10*cos(2*pi*2*(loc[,1]+loc[,2])))
#' loc = matrix(runif(n*2),n,2)
#' ## One field, 2 observations per location
#' idx.y = rep(1:n,2)
#' y = field.fcn(loc[idx.y,]) + rnorm(length(idx.y))
#' 
#' mesh = inla.mesh.create(loc, refine=list(max.edge=0.05))
#' spde = inla.spde2.matern(mesh)
#' data = list(y=y, field=mesh$idx$loc[idx.y])
#' formula = y ~ -1 + f(field, model=spde)
#' result = inla(formula, data=data, family="normal")
#' 
#' ## Plot the mesh structure:
#' plot(mesh)
#' \donttest{
#' if (require(rgl)) {
#'   col.pal = colorRampPalette(c("blue","cyan","green","yellow","red"))
#'   ## Plot the posterior mean:
#'   plot(mesh, rgl=TRUE,
#'        result$summary.random$field[,"mean"],
#'        color.palette = col.pal)
#'   ## Plot residual field:
#'   plot(mesh, rgl=TRUE,
#'        result$summary.random$field[,"mean"]-field.fcn(mesh$loc),
#'        color.palette = col.pal)
#' }
#' }
#' result.field = inla.spde.result(result, "field", spde)
#' plot(result.field$marginals.range.nominal[[1]])
#' 
#' @export inla.spde2.matern
inla.spde2.matern =
    function(mesh,
             alpha=2,
             param = NULL,
             constr = FALSE,
             extraconstr.int = NULL,
             extraconstr = NULL,
             fractional.method = c("parsimonious", "null"),
             B.tau = matrix(c(0,1,0),1,3),
             B.kappa = matrix(c(0,0,1),1,3),
             prior.variance.nominal = 1,
             prior.range.nominal = NULL,
             prior.tau = NULL,
             prior.kappa = NULL,
             theta.prior.mean = NULL,
             theta.prior.prec = 0.1,
             n.iid.group = 1,
             ...)
{
  if(!is.null(list(...)[["prior.pc.rho"]]) &&
     !is.null(list(...)[["prior.pc.sig"]])){
      ## Temporary implementation of PC prior for standard deviation and range
      ##    - Changes parametrization to range and standard deviation
      ##    - Sets prior according to hyperparameters for range   : prior.pc.rho
      ##                                              and std.dev.: prior.pc.sig
      warning("You're using a deprecated experimental PC prior matern model that will be removed in a future version of the package. Use 'inla.spde2.pcmatern' instead.")
      prior.pc.rho <- list(...)[["prior.pc.rho"]]
      prior.pc.sig <- list(...)[["prior.pc.sig"]]

      ## Call inla.spde2.matern with range and standard deviation parametrization
      d = inla.ifelse(inherits(mesh, "inla.mesh"), 2, 1)
      nu = alpha-d/2
      kappa0 = log(8*nu)/2
      tau0   = 0.5*(lgamma(nu)-lgamma(nu+d/2)-d/2*log(4*pi))-nu*kappa0
      spde   = inla.spde2.matern(mesh = mesh,
                                 B.tau   = cbind(tau0,   nu,  -1),
                                 B.kappa = cbind(kappa0, -1, 0))

      ## Change prior information
      param = c(prior.pc.rho, prior.pc.sig)
      spde$f$hyper.default$theta1$prior = "pcspdega"
      spde$f$hyper.default$theta1$param = param
      spde$f$hyper.default$theta1$initial = log(prior.pc.rho[1])+1
      spde$f$hyper.default$theta2$initial = log(prior.pc.sig[1])-1

      ## End and return
      return(invisible(spde))
    }

    ## Standard code
    inla.require.inherits(mesh, c("inla.mesh", "inla.mesh.1d"), "'mesh'")
    fractional.method = match.arg(fractional.method)

    if (is.null(param)) {
        param =
            param2.matern.orig(
                mesh, alpha,
                B.tau, B.kappa,
                prior.variance.nominal,
                prior.range.nominal,
                prior.tau,
                prior.kappa,
                theta.prior.mean,
                theta.prior.prec)
    } else {
        deprecated =
            !c(missing(B.tau), missing(B.kappa),
               missing(prior.variance.nominal),
               missing(prior.range.nominal),
               missing(prior.tau), missing(prior.kappa),
               missing(theta.prior.mean), missing(theta.prior.prec))
        deprecated =
            c("B.tau", "B.kappa",
              "prior.variance.nominal",
              "prior.range.nominal",
              "prior.tau", "prior.kappa",
              "theta.prior.mean", "theta.prior.prec")[deprecated]
        if (length(deprecated) > 0) {
            warning(paste("'param' specified;  ",
                          "Ignoring deprecated parameter(s) ",
                          paste(deprecated, collapse=", "), ".", sep=""))
        }
    }

    d <- fmesher::fm_manifold_dim(fmesher::fm_as_fm(mesh))
    nu = alpha-d/2
    nu.nominal = max(0.5, nu)
    alpha.nominal = max(nu.nominal+d/2, alpha)

    n.spde = inla.ifelse(d==2, mesh$n, mesh$m)
    n.theta = ncol(B.kappa)-1L

    if (fmesher_deprecate_allow(2L)) {
        fem <- fmesher::fm_fem(fmesher::fm_as_fm(mesh), order = 2)
        if ((d == 1) && (mesh$degree == 2)) {
            fem$c0 <- fem$c1 ## Use higher order matrix.
        }
    } else {
        if (d == 2) {
            fem <-
                inla.fmesher.smorg(mesh$loc,
                                   mesh$graph$tv,
                                   fem=2,
                                   output=list("c0", "c1", "g1", "g2"))
        } else {
            fem <- inla.mesh.1d.fem(mesh)
            if (mesh$degree == 2) {
                fem$c0 <- fem$c1 ## Use higher order matrix.
            }
        }
    }
    
    if (alpha==2) {
        B.phi0 = param$B.tau
        B.phi1 = 2*param$B.kappa
        M0 = fem$c0
        M1 = fem$g1
        M2 = fem$g2
    } else if (alpha==1) {
        B.phi0 = param$B.tau
        B.phi1 = param$B.kappa
        M0 = fem$c0
        M1 = fem$g1*0
        M2 = fem$g1
    } else if (!param$is.stationary) {
        stop("Non-stationary Matern with fractional alpha is not implemented.")
    } else if ((alpha<2) && (alpha>1)) {
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
        B.phi0 = param$B.tau + (alpha-2)*param$B.kappa
        B.phi1 = 2*param$B.kappa
        M0 = fem$c0*b[1]
        M1 = fem$g1*b[2]/2
        M2 = fem$g2*b[3]
    } else if ((alpha<1) && (alpha>0)) {
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
        B.phi0 = param$B.tau + (alpha-1)*param$B.kappa
        B.phi1 = param$B.kappa
        M0 = fem$c0*b[1]
        M1 = fem$g1*0
        M2 = fem$g1*b[2]
    } else {
        stop(paste("Unsupported alpha value (", alpha,
                   "). Supported values are 0 < alpha <= 2", sep=""))
    }

    if (n.iid.group == 1) {
      spde =
        inla.spde2.generic(M0=M0, M1=M1, M2=M2,
                           B0=B.phi0, B1=B.phi1, B2=1,
                           theta.mu = param$theta.prior.mean,
                           theta.Q = param$theta.prior.prec,
                           transform = "identity",
                           BLC = param$BLC)
    } else {
      if (nrow(B.phi0) > 1) {
        B.phi0 <- kronecker(matrix(1, n.iid.group, 1), B.phi0)
      }
      if (nrow(B.phi1) > 1) {
        B.phi1 <- kronecker(matrix(1, n.iid.group, 1), B.phi1)
      }
      spde =
        inla.spde2.generic(M0=kronecker(Diagonal(n.iid.group), M0),
                           M1=kronecker(Diagonal(n.iid.group), M1),
                           M2=kronecker(Diagonal(n.iid.group), M2),
                           B0=B.phi0, B1=B.phi1, B2=1,
                           theta.mu = param$theta.prior.mean,
                           theta.Q = param$theta.prior.prec,
                           transform = "identity",
                           BLC = param$BLC)
    }
    spde$model = "matern"
    spde$BLC = param$BLC

    if (constr || !is.null(extraconstr.int) || !is.null(extraconstr)) {
      A.constr = matrix(numeric(0), 0, n.spde*n.iid.group)
      e.constr = matrix(numeric(0), 0, 1)
      if (constr) {
        A.constr <- rbind(A.constr,
                          matrix(colSums(fem$c1)/n.iid.group,
                                 1, n.spde*n.iid.group))
        e.constr <- rbind(e.constr, 0)
      }
      if (!is.null(extraconstr.int)) {
        if (ncol(extraconstr.int$A) == n.spde) {
          A.constr <-
            rbind(A.constr,
                  kronecker(matrix(1/n.iid.group, 1, n.iid.group),
                            as.matrix(extraconstr.int$A %*% fem$c1)))
        } else {
          A.constr <-
            rbind(A.constr,
                  as.matrix(extraconstr.int$A %*%
                              kronecker(Diagonal(n.iid.group),
                                        fem$c1)))
        }
        e.constr <- rbind(e.constr, as.matrix(extraconstr.int$e))
      }
      if (!is.null(extraconstr)) {
        if (ncol(extraconstr$A) == n.spde) {
          A.constr <-
            rbind(A.constr,
                  kronecker(matrix(1/n.iid.group, 1, n.iid.group),
                            as.matrix(extraconstr$A)))
        } else {
          A.constr <- rbind(A.constr, as.matrix(extraconstr$A))
        }
        e.constr <- rbind(e.constr, as.matrix(extraconstr$e))
      }

      spde$f$constr = FALSE
      spde$f$extraconstr = list(A=A.constr, e=e.constr)
    }

    ## Attach the mesh, so downstream code can have access
    spde$mesh <- mesh

    return(invisible(spde))
}







#' Matern SPDE model object with PC prior for INLA
#' 
#' Create an `inla.spde2` model object for a Matern model, using a PC
#' prior for the parameters.
#' 
#' This method constructs a Matern SPDE model, with spatial range \eqn{\rho}
#' and standard deviation parameter \eqn{\sigma}.  In the parameterisation
#' 
#' \deqn{(\kappa^2-\Delta)^{\alpha/2}(\tau }{(kappa^2-Delta)^(alpha/2) (tau
#' x(u)) = W(u)}\deqn{ x(u))=W(u)}{(kappa^2-Delta)^(alpha/2) (tau x(u)) = W(u)}
#' 
#' the spatial scale parameter \eqn{\kappa=\sqrt{8\nu}/\rho}, where
#' \eqn{\nu=\alpha-d/2}, and \eqn{\tau} is proportional to \eqn{1/\sigma}.
#' 
#' Stationary models are supported for \eqn{0 < \alpha \leq 2}{0 < alpha <= 2},
#' with spectral approximation methods used for non-integer \eqn{\alpha}, with
#' approximation method determined by `fractional.method`.
#' 
#' Integration and other general linear constraints are supported via the
#' `constr`, `extraconstr.int`, and `extraconstr` parameters,
#' which also interact with `n.iid.group`.
#' 
#' The joint PC prior density for the spatial range, \eqn{\rho}, and the
#' marginal standard deviation, \eqn{\sigma}, and is \deqn{ }{p(rho, sigma) =
#' (d R)/2 rho^(-1-d/2) exp(-R rho^(-d/2)) S exp(-S sigma) }\deqn{ \pi(\rho,
#' \sigma) = }{p(rho, sigma) = (d R)/2 rho^(-1-d/2) exp(-R rho^(-d/2)) S exp(-S
#' sigma) }\deqn{ \frac{d \lambda_\rho}{2} \rho^{-1-d/2} \exp(-\lambda_\rho
#' \rho^{-d/2}) }{p(rho, sigma) = (d R)/2 rho^(-1-d/2) exp(-R rho^(-d/2)) S
#' exp(-S sigma) }\deqn{ \lambda_\sigma\exp(-\lambda_\sigma \sigma) }{p(rho,
#' sigma) = (d R)/2 rho^(-1-d/2) exp(-R rho^(-d/2)) S exp(-S sigma) } where
#' \eqn{\lambda_\rho}{R} and \eqn{\lambda_\sigma}{S} are hyperparameters that
#' must be determined by the analyst. The practical approach for this in INLA
#' is to require the user to indirectly specify these hyperparameters through
#' \deqn{P(\rho < \rho_0) = p_\rho} and \deqn{P(\sigma > \sigma_0) = p_\sigma}
#' where the user specifies the lower tail quantile and probability for the
#' range (\eqn{\rho_0} and \eqn{p_\rho}) and the upper tail quantile and
#' probability for the standard deviation (\eqn{\sigma_0} and
#' \eqn{\alpha_\sigma}).
#' 
#' This allows the user to control the priors of the parameters by supplying
#' knowledge of the scale of the problem. What is a reasonable upper magnitude
#' for the spatial effect and what is a reasonable lower scale at which the
#' spatial effect can operate? The shape of the prior was derived through a
#' construction that shrinks the spatial effect towards a base model of no
#' spatial effect in the sense of distance measured by Kullback-Leibler
#' divergence.
#' 
#' The prior is constructed in two steps, under the idea that having a spatial
#' field is an extension of not having a spatial field. First, a spatially
#' constant random effect (\eqn{\rho = \infty}) with finite variance is more
#' complex than not having a random effect (\eqn{\sigma = 0}). Second, a
#' spatial field with spatial variation (\eqn{\rho < \infty}) is more complex
#' than the random effect with no spatial variation. Each of these extensions
#' are shrunk towards the simpler model and, as a result, we shrink the spatial
#' field towards the base model of no spatial variation and zero variance
#' (\eqn{\rho = \infty} and \eqn{\sigma = 0}).
#' 
#' The details behind the construction of the prior is presented in Fuglstad,
#' et al. (2016) and is based on the PC prior framework (Simpson, et al.,
#' 2015).
#' 
#' @param mesh The mesh to build the model on, as an [inla.mesh()] or
#' [inla.mesh.1d()] object.
#' @param alpha Fractional operator order, \eqn{0<\alpha\leq 2}{0 < alpha <= 2}
#' supported, for \eqn{\nu=\alpha-d/2>0}.
#' @param param Further model parameters. Not currently used.
#' @param constr If `TRUE`, apply an integrate-to-zero constraint.
#' Default `FALSE`.
#' @param extraconstr.int Field integral constraints.
#' @param extraconstr Direct linear combination constraints on the basis
#' weights.
#' @param fractional.method Specifies the approximation method to use for
#' fractional (non-integer) `alpha` values. `'parsimonious'` gives an
#' overall approximate minimal covariance error, `'null'` uses
#' approximates low-order properties.
#' @param n.iid.group If greater than 1, build an explicitly iid replicated
#' model, to support constraints applied to the combined replicates, for
#' example in a time-replicated spatial model. Constraints can either be
#' specified for a single mesh, in which case it's applied to the average of
#' the replicates (`ncol(A)` should be `mesh$n` for 2D meshes,
#' `mesh$m` for 1D), or as general constraints on the collection of
#' replicates (`ncol(A)` should be `mesh$n * n.iid.group` for 2D
#' meshes, `mesh$m * n.iid.group` for 1D).
#' @param prior.range A length 2 vector, with `(range0,Prange)` specifying
#' that \eqn{P(\rho < \rho_0)=p_\rho}, where \eqn{\rho} is the spatial range of
#' the random field. If `Prange` is `NA`, then `range0` is used
#' as a fixed range value.
#' @param prior.sigma A length 2 vector, with `(sigma0,Psigma)` specifying
#' that \eqn{P(\sigma > \sigma_0)=p_\sigma}, where \eqn{\sigma} is the marginal
#' standard deviation of the field.  If `Psigma` is `NA`, then
#' `sigma0` is used as a fixed range value.
#' @return An `inla.spde2` object.
#' @author Finn Lindgren \email{finn.lindgren@@gmail.com}
#' @seealso [inla.mesh.2d()], [inla.mesh.create()],
#' [inla.mesh.1d()], [inla.mesh.basis()],
#' [inla.spde2.matern()], [inla.spde2.generic()]
#' @references Fuglstad, G.-A., Simpson, D., Lindgren, F., and Rue, H. (2016)
#' Constructing Priors that Penalize the Complexity of Gaussian Random Fields.
#' arXiv:1503.00256
#' 
#' Simpson, D., Rue, H., Martins, T., Riebler, A., and Sørbye, S. (2015)
#' Penalising model component complexity: A principled, practical approach to
#' constructing priors. arXiv:1403.4630
#' @examples
#' 
#' ## Spatial interpolation
#'   n = 100
#'   field.fcn = function(loc) (10*cos(2*pi*2*(loc[,1]+loc[,2])))
#'   loc = matrix(runif(n*2),n,2)
#'   ## One field, 2 observations per location
#'   idx.y = rep(1:n,2)
#'   y = field.fcn(loc[idx.y,]) + rnorm(length(idx.y))
#' 
#'   mesh = inla.mesh.2d(loc, max.edge=0.05, cutoff=0.01)
#'   spde = inla.spde2.pcmatern(mesh,
#'            prior.range=c(0.01,0.1), prior.sigma=c(100,0.1))
#'   data = list(y=y, field=mesh$idx$loc[idx.y])
#'   formula = y ~ -1 + f(field, model=spde)
#'   result = inla(formula, data=data, family="normal")
#' 
#'   ## Plot the mesh structure:
#'   plot(mesh)
#'   \donttest{
#'   if (require(rgl)) {
#'     col.pal = colorRampPalette(c("blue","cyan","green","yellow","red"))
#'     ## Plot the posterior mean:
#'     plot(mesh, rgl=TRUE,
#'          result$summary.random$field[,"mean"],
#'          color.palette = col.pal)
#'     ## Plot residual field:
#'     plot(mesh, rgl=TRUE,
#'          result$summary.random$field[,"mean"]-field.fcn(mesh$loc),
#'          color.palette = col.pal)
#'   }
#'   }
#' 
#'   result.field = inla.spde.result(result, "field", spde)
#'   par(mfrow=c(2,1))
#'   plot(result.field$marginals.range.nominal[[1]],
#'        type="l", main="Posterior density for range")
#'   plot(inla.tmarginal(sqrt, result.field$marginals.variance.nominal[[1]]),
#'        type="l", main="Posterior density for std.dev.")
#'   par(mfrow=c(1,1))
#' 
#' ## Spatial model
#'   set.seed(1234234)
#'   
#'   ## Generate spatial locations
#'   nObs = 200
#'   loc = matrix(runif(nObs*2), nrow = nObs, ncol = 2)
#' 
#'   ## Generate observation of spatial field
#'   nu = 1.0
#'   rhoT = 0.2
#'   kappaT = sqrt(8*nu)/rhoT
#'   sigT = 1.0
#'   Sig = sigT^2*inla.matern.cov(nu = nu,
#'                                kappa = kappaT,
#'                                x = as.matrix(dist(loc)),
#'                                d = 2,
#'                                corr = TRUE)
#'   L = t(chol(Sig))
#'   u = L %*% rnorm(nObs)
#' 
#'   ## Construct observation with nugget
#'   sigN = 0.1
#'   y = u + sigN*rnorm(nObs)
#' 
#'   ## Create the mesh and spde object
#'   mesh = inla.mesh.2d(loc,
#'                       max.edge = 0.05,
#'                       cutoff = 0.01)
#'   spde = inla.spde2.pcmatern(mesh,
#'                              prior.range = c(0.01, 0.05),
#'                              prior.sigma = c(10, 0.05))
#' 
#'   ## Create projection matrix for observations
#'   A = inla.spde.make.A(mesh = mesh,
#'                        loc = loc)
#' 
#'   ## Run model without any covariates
#'   idx = 1:spde$n.spde
#'   res = inla(y ~ f(idx, model = spde) - 1,
#'              data = list(y = y, idx = idx, spde = spde),
#'              control.predictor = list(A = A))
#' 
#' ## Re-run model with fixed range
#'   spde.fixed = inla.spde2.pcmatern(mesh,
#'                                    prior.range = c(0.2, NA),
#'                                    prior.sigma = c(10, 0.05))
#' 
#'   res.fixed = inla(y ~ f(idx, model = spde) - 1,
#'                    data = list(y = y, idx = idx, spde = spde.fixed),
#'                    control.predictor = list(A = A))
#' 
#' 
#' @export inla.spde2.pcmatern
inla.spde2.pcmatern =
    function(mesh,
             alpha = 2,
             param = NULL,
             constr = FALSE,
             extraconstr.int = NULL,
             extraconstr = NULL,
             fractional.method = c("parsimonious", "null"),
             n.iid.group = 1,
             prior.range = NULL,
             prior.sigma = NULL)
{
  ## Implementation of PC prior for standard deviation and range
  ##    - Sets the parametrization to range and standard deviation
  ##    - Sets prior according to hyperparameters for range   : prior.range
  ##                                              and std.dev.: prior.sigma
  ## Calls inla.spde2.matern to construct the object, then changes the prior
  if (inherits(mesh, "inla.mesh")) {
    d <- 2
  } else if (inherits(mesh, "inla.mesh.1d")) {
    d <- 1
  } else {
    stop(paste("Unknown mesh class '",
               paste(class(mesh), collapse=",", sep=""),
               "'.", sep=""))
  }

  if (missing(prior.range) || is.null(prior.range) ||
      !is.vector(prior.range) || (length(prior.range) != 2)) {
    stop("'prior.range' should be a length 2 vector 'c(range0,tailprob)' or a fixed range specified with 'c(range,NA)'.")
  }
  if (missing(prior.sigma) || is.null(prior.sigma) ||
      !is.vector(prior.sigma) || (length(prior.sigma) != 2)) {
    stop("'prior.sigma' should be a length 2 vector 'c(sigma0,tailprob)' or a fixed sigma specified with 'c(sigma,NA)'.")
  }
  if (prior.range[1] <= 0){
    stop("'prior.range[1]' must be a number greater than 0 specifying a spatial range")
  }
  if (prior.sigma[1] <= 0){
    stop("'prior.sigma[1]' must be a number greater than 0 specifying a standard deviation")
  }
  if (!is.na(prior.range[2]) &&
      ((prior.range[2] <= 0) || (prior.range[2] >= 1))) {
    stop("'prior.range[2]' must be a probaility strictly between 0 and 1 (or NA to specify a fixed range)")
  }
  if (!is.na(prior.sigma[2]) &&
      ((prior.sigma[2] <= 0) || (prior.sigma[2] >= 1))) {
    stop("'prior.sigma[2]' must be a probaility strictly between 0 and 1 (or NA to specify a fixed sigma)")
  }

  nu <- alpha-d/2
  if (nu <= 0) {
    stop(paste("Smoothness nu = alpha-dim/2 = ", nu,
               ", but must be > 0.", sep=""))
  }

  kappa0 <- log(8*nu)/2
  tau0   <- 0.5*(lgamma(nu)-lgamma(nu+d/2)-d/2*log(4*pi))-nu*kappa0
  spde   <- inla.spde2.matern(mesh = mesh,
                              B.tau   = cbind(tau0,   nu,  -1),
                              B.kappa = cbind(kappa0, -1, 0),
                              alpha = alpha,
                              param = NULL,
                              constr = constr,
                              extraconstr.int = extraconstr.int,
                              extraconstr = extraconstr,
                              fractional.method = fractional.method,
                              n.iid.group = n.iid.group)

  ## Calculate hyperparameters
  is.fixed.range <- is.na(prior.range[2])
  if (is.fixed.range) {
    lam1 <- 0
    initial.range <- log(prior.range[1])
  } else {
    lam1 <- -log(prior.range[2])*prior.range[1]^(d/2)
    initial.range <- log(prior.range[1]) + 1
  }

  is.fixed.sigma <- is.na(prior.sigma[2])
  if (is.fixed.sigma){
    lam2 <- 0
    initial.sigma <- log(prior.sigma[1])
  } else{
    lam2 <- -log(prior.sigma[2])/prior.sigma[1]
    initial.sigma <- log(prior.sigma[1]) - 1
  }

  pcmatern.param = c(lam1, lam2, d)

  ## Change prior information
  spde$f$hyper.default <-
    list(theta1=list(prior="pcmatern",
                     param=pcmatern.param,
                     initial=initial.range,
                     fixed=is.fixed.range),
         theta2=list(initial=initial.sigma,
                     fixed=is.fixed.sigma))

  ## Change the model descriptor
  spde$model = "pcmatern"

  invisible(spde)
}





param2.iheat <- function(mesh.space,
                         mesh.time,
                         gamma.E = NULL,
                         alpha.E = NULL,
                         theta.prior.mean = NULL,
                         theta.prior.prec = NULL)
{
    return(list(gamma.E=gamma.E,
                alpha.E=alpha.E,
                theta.prior.mean=theta.prior.mean,
                theta.prior.prec=theta.prior.prec))
}

inla.spde2.iheat =
    function(mesh.space,
             mesh.time,
             order=2,
             param = NULL,
             constr = FALSE,
             extraconstr.int = NULL,
             extraconstr = NULL)
{
    ## order=1:
    ## gamma.t dot(u) - Laplacian (u) = E
    ##
    ## order=2:
    ## gamma.t dot(u) - Laplacian (u) = v
    ## gamma.t dot(v) - Laplacian (v) = E
    ##
    ## (1 - gamma.E Laplacian)^(alpha.E/2) E = W/sqrt(gamma.s)
    ## gamma.E and alpha.E known

    inla.require.inherits(mesh.space, c("inla.mesh", "inla.mesh.1d"),
                          "'mesh.space'")
    inla.require.inherits(mesh.time, c("inla.mesh.1d"),
                          "'mesh.time'")

  if (is.null(param)) {
    stop("Use param2.iheat() to construct the prior parameter settings")
    ##        param =
    ##            param2.iheat(
    ##                mesh.space, mesh.time,
    ##                theta.prior.mean,
    ##                theta.prior.prec)
  }

    d.space = inla.ifelse(inherits(mesh.space, "inla.mesh.1d"), 1, 2)
    d.time = 1
    n.space = inla.ifelse(d.space==2, mesh.space$n, mesh.space$m)
    n.time = mesh.time$m
    n.spde = n.space*n.time
    n.theta = 2L ## gamma.s, gamma.t

    if (d.space==2) {
        fem.space = inla.mesh.fem(mesh.space$loc, 2)
    } else {
        fem.space = inla.mesh.1d.fem(mesh.space)
        if (mesh.space$degree==2) {
            fem.space$c0 = fem.space$c1 ## Use higher order matrix.
        }
    }
    fem.time = inla.mesh.1d.fem(mesh.time)
    if (mesh.time$degree==2) {
        fem.time$c0 = fem.time$c1 ## Use higher order matrix.
    }

  ## TODO: the rest
  if (FALSE) {

    if (alpha==2) {
        B.phi0 = param$B.tau
        B.phi1 = 2*param$B.kappa
        M0 = fem$c0
        M1 = fem$g1
        M2 = fem$g2
    } else if (alpha==1) {
        B.phi0 = param$B.tau
        B.phi1 = param$B.kappa
        M0 = fem$c0
        M1 = fem$g1*0
        M2 = fem$g1
    } else if (!param$is.stationary) {
        stop("Non-stationary Matern with fractional alpha is not implemented.")
    } else if ((alpha<2) && (alpha>1)) {
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
        B.phi0 = param$B.tau + (alpha-2)*param$B.kappa
        B.phi1 = 2*param$B.kappa
        M0 = fem$c0*b[1]
        M1 = fem$g1*b[2]/2
        M2 = fem$g2*b[3]
    } else if ((alpha<1) && (alpha>0)) {
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
        B.phi0 = param$B.tau + (alpha-1)*param$B.kappa
        B.phi1 = param$B.kappa
        M0 = fem$c0*b[1]
        M1 = fem$g1*0
        M2 = fem$g1*b[2]
    } else {
        stop(paste("Unsupported alpha value (", alpha,
                   "). Supported values are 0 < alpha <= 2", sep=""))
    }

    spde =
        inla.spde2.generic(M0=M0, M1=M1, M2=M2,
                           B0=B.phi0, B1=B.phi1, B2=1,
                           theta.mu = param$theta.prior.mean,
                           theta.Q = param$theta.prior.prec,
                           transform = "identity",
                           BLC = param$BLC)
    spde$model = "matern"
    spde$BLC = param$BLC

    if (constr || !is.null(extraconstr.int)) {
        if (constr) {
            A.constr = matrix(colSums(fem$c1), 1, n.spde)
            e.constr = 0
        } else {
            A.constr = matrix(numeric(0), 0, n.spde)
            e.constr = c()
        }
        if (!is.null(extraconstr.int)) {
            A.constr =
                rbind(A.constr,
                      matrix(extraconstr.int$A %*% fem$c1,
                             nrow(extraconstr.int$A), n.spde))
            e.constr = c(e.constr, extraconstr.int$e)
        }
        if (!is.null(extraconstr)) {
            A.constr = rbind(A.constr, extraconstr$A)
            e.constr = c(e.constr, extraconstr$e)
        }

        spde$f$constr = FALSE
        spde$f$extraconstr = list(A=A.constr, e=e.constr)
    } else if (!is.null(extraconstr)) {
        spde$f$constr = FALSE
        spde$f$extraconstr = extraconstr
    }

    return(invisible(spde))
  }

  invisible(NULL)
}








#' @export
#' @param spde An spde model object
#' @param theta Parameters in the model's internal scale
#' @describeIn inla.spde2.matern Convert from theta vector to phi0 values in
#' the internal spde2 model representation
inla.spde2.theta2phi0 = function(spde, theta)
{
    inla.require.inherits(spde, "inla.spde2", "'spde'")

    if (spde$n.theta>0)
        return(exp(spde$param.inla$B0[,1, drop=TRUE] +
                   spde$param.inla$B0[,-1, drop=FALSE] %*% theta))
    else
        return(exp(spde$param.inla$B0[,1, drop=TRUE]))
}

#' @export
#' @describeIn inla.spde2.matern Convert from theta vector to phi1 values in
#' the internal spde2 model representation
inla.spde2.theta2phi1 = function(spde, theta)
{
    inla.require.inherits(spde, "inla.spde2", "'spde'")

    if (spde$n.theta>0)
        return(exp(spde$param.inla$B1[,1, drop=TRUE] +
                   spde$param.inla$B1[,-1, drop=FALSE] %*% theta))
    else
        return(exp(spde$param.inla$B1[,1, drop=TRUE]))
}

#' @export
#' @describeIn inla.spde2.matern Convert from theta vector to phi2 values in
#' the internal spde2 model representation
inla.spde2.theta2phi2 = function(spde, theta)
{
    inla.require.inherits(spde, "inla.spde2", "'spde'")

    if (spde$n.theta>0) {
        phi = ((spde$param.inla$B2[,1, drop=TRUE] +
                spde$param.inla$B2[,-1, drop=FALSE] %*% theta))
    } else {
        phi = ((spde$param.inla$B2[,1, drop=TRUE]))
    }
    if (spde$param.inla$transform == "identity") {
        return(phi)
    } else if (spde$param.inla$transform == "logit") {
        return(cos(pi/(1+exp(-phi))))
    } else if (spde$param.inla$transform == "log") {
        return(2*exp(phi)-1)
    } else {
        warning(paste("Unknown link function '",
                      spde$param.inla$transform,
                      "' phi2.  Using identity link instead.", sep=""))
        return(phi)
    }
}

#' @export
#' @rdname inla.spde.precision
inla.spde2.precision =
    function(spde, theta=NULL,
             phi0=inla.spde2.theta2phi0(spde, theta),
             phi1=inla.spde2.theta2phi1(spde, theta),
             phi2=inla.spde2.theta2phi2(spde, theta),
             ...)
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


#' @export
#' @rdname inla.spde.result
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
    if (!is.null(inla$marginals$hyperpar)) {
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


#' @export
#' @rdname inla.spde.models
inla.spde2.models = function()
{
    return(c("generic", "matern", "pcmatern"))
}


## spde.common-connections:
#' @export
#' @method inla.spde.precision inla.spde2
#' @rdname inla.spde.precision
inla.spde.precision.inla.spde2 = inla.spde2.precision
#' @export
#' @method inla.spde.result inla.spde2
#' @rdname inla.spde.result
inla.spde.result.inla.spde2 = inla.spde2.result


## Deprecated:
##inla.spde.generic2 = inla.spde2.generic
