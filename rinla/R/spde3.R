## 'spde3' model functions
## External: inla.spde.precision!inla.spde3
## External: inla.spde.result!inla.spde3 inla.spde3.generic
## External: inla.spde3.matern param2.matern.orig
## External: inla.spde3.matern.sd.basis inla.spde3.models
## External: inla.spde3.precision inla.spde3.result
## External: inla.spde3.theta2phi0 inla.spde3.theta2phi1
## External: inla.spde3.theta2phi2 inla.spde3.theta2phi3
##
## Internal: inla.internal.spde3.matern.B.tau
## Internal: inla.internal.test.spde3.sd.basis





inla.spde3.generic <-
    function(M0, M1, M2, M3, B0, B1, B2, B3, theta.mu, theta.Q,
             transform = c("identity","log","logit","shiftedlog","oldlogit"),
             theta.initial = theta.mu,
             fixed = rep(FALSE, length(theta.mu)),
             theta.fixed = theta.initial[fixed],
             BLC = cbind(0.0, diag(nrow=length(theta.mu))),
             ...)
{
    warning("The spde3 interface is still experimental and may be changed without notice.")
    transform = match.arg(transform)

    if (!inla.is.matrix(M0)) stop("M0 must be of a matrix type.")
    if (!inla.is.matrix(M1)) stop("M1 must be of a matrix type.")
    if (!inla.is.matrix(M2)) stop("M2 must be of a matrix type.")
    if (!inla.is.matrix(M3)) stop("M3 must be of a matrix type.")
    M0 = inla.as.dgTMatrix(M0)
    M1 = inla.as.dgTMatrix(M1)
    M2 = inla.as.dgTMatrix(M2)
    M3 = inla.as.dgTMatrix(M3)

    n.spde = nrow(M0)
    n.theta = length(theta.mu)

    spde = (list(model = "generic",
                 n.spde = n.spde,
                 n.theta = n.theta,
                 param.inla = list(),
                 f = list()
                 ))
    class(spde) = c("inla.spde3", "inla.spde", "inla.model.class")

    B0 = inla.spde.homogenise_B_matrix(B0, n.spde, n.theta)
    B1 = inla.spde.homogenise_B_matrix(B1, n.spde, n.theta)
    B2 = inla.spde.homogenise_B_matrix(B2, n.spde, n.theta)
    B3 = inla.spde.homogenise_B_matrix(B3, n.spde, n.theta)
    BLC = inla.spde.homogenise_B_matrix(BLC, nrow(BLC), n.theta)

    param.inla =
        list(n = n.spde,
             n.theta = n.theta,
             M0=M0, M1=M1, M2=M2, M3=M3,
             B0=B0, B1=B1, B2=B2, B3=B3, BLC=BLC,
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
        param.inla$B3[,1] =
            B3[,c(TRUE, fixed), drop=FALSE] %*% c(1.0, theta.fixed)
        param.inla$BLC[,1] =
            BLC[,c(TRUE, fixed), drop=FALSE] %*% c(1.0, theta.fixed)
        param.inla$B0 = param.inla$B0[,c(TRUE, !fixed), drop=FALSE]
        param.inla$B1 = param.inla$B1[,c(TRUE, !fixed), drop=FALSE]
        param.inla$B2 = param.inla$B2[,c(TRUE, !fixed), drop=FALSE]
        param.inla$B3 = param.inla$B3[,c(TRUE, !fixed), drop=FALSE]
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

    spde$f = (list(model="spde3",
                   n=n.spde,
                   spde3.transform=transform,
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









inla.internal.spde3.matern.B.tau =
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



inla.spde3.matern.sd.basis =
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
            inla.spde3.matern(mesh, B.tau=as.matrix(0), B.kappa=B.kappa, ...)

        B.tau =
            inla.internal.spde3.matern.B.tau(
                precision.fcn=
                function(theta, spde) {
                    return(inla.spde3.precision(spde, theta=theta))
                },
                n.theta = p.kappa,
                spde = spde)
        if (!local.offset.compensation) {
            B.tau[,1] = rep(mean(B.tau[,1]), spde$n.spde)
        }
    } else {
        spde =
            inla.spde3.matern(mesh,
                              B.tau=as.matrix(0),
                              B.kappa=B.kappa[,1,drop=FALSE],
                              alpha=alpha, ...)
        lvar.0 =
            log(diag(inla.qinv(inla.spde3.precision(spde, theta=rep(0,0)))))
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



inla.internal.test.spde3.sd.basis = function (k=1, dth=0.1, r=1000, globe=25, compensate=FALSE) {

    alpha=2
    mesh = inla.mesh.create(globe=globe)

    B.common = inla.mesh.basis(mesh, "b.spline", n=3)
    B.common = cbind(B.common, B.common[,1]*(mesh$loc[,3]>0))
    B.common[,1] = B.common[,1]*(mesh$loc[,3]<=0)

    B.sd = cbind(0, B.common)
    B.range = cbind(log(r/6370), B.common)

    sd.basis0 = inla.spde3.matern.sd.basis(mesh, B.sd=B.sd, B.range=B.range, method=0, alpha=alpha, local.offset.compensation=compensate)
    sd.basis1 = inla.spde3.matern.sd.basis(mesh, B.sd=B.sd, B.range=B.range, method=1, alpha=alpha, local.offset.compensation=compensate)

    spde0 =
        inla.spde3.matern(mesh, alpha=alpha,
                          B.tau=sd.basis0$B.tau,
                          B.kappa=sd.basis0$B.kappa)
    spde1 =
        inla.spde3.matern(mesh, alpha=alpha,
                          B.tau=sd.basis1$B.tau,
                          B.kappa=sd.basis1$B.kappa)

    th0 = c(0,0,0,0,0,0,0,0)
    th1 = th0
    th1[k] = th0[k]+dth
    v0.0=diag(inla.qinv(inla.spde3.precision(spde0,theta=th0)))
    v0.1=diag(inla.qinv(inla.spde3.precision(spde0,theta=th1)))
    v1.0=diag(inla.qinv(inla.spde3.precision(spde1,theta=th0)))
    v1.1=diag(inla.qinv(inla.spde3.precision(spde1,theta=th1)))

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

inla.spde3.matern =
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
             theta.prior.prec = 0.1)
{
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

    d = inla.ifelse(inherits(mesh, "inla.mesh"), 2, 1)
    nu = alpha-d/2
    nu.nominal = max(0.5, nu)
    alpha.nominal = max(nu.nominal+d/2, alpha)

    n.spde = inla.ifelse(d==2, mesh$n, mesh$m)
    n.theta = ncol(B.kappa)-1L

    if (d==2) {
        fem =
            inla.fmesher.smorg(mesh$loc,
                               mesh$graph$tv,
                               fem=2,
                               output=list("c0", "c1", "g1", "g2"))
    } else {
        fem = inla.mesh.1d.fem(mesh)
        if (mesh$degree==2) {
            fem$c0 = fem$c1 ## Use higher order matrix.
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

    spde =
        inla.spde3.generic(M0=M0, M1=M1, M2=M2,
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








param3.iheat <- function(mesh.space,
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

inla.spde3.iheat =
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
    stop("Use param3.iheat() to construct the prior parameter settings")
    ##        param =
    ##            param3.iheat(
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
        inla.spde3.generic(M0=M0, M1=M1, M2=M2,
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










inla.spde3.theta2phi0 = function(spde, theta)
{
    inla.require.inherits(spde, "inla.spde3", "'spde'")

    if (spde$n.theta>0)
        return(exp(spde$param.inla$B0[,1, drop=TRUE] +
                   spde$param.inla$B0[,-1, drop=FALSE] %*% theta))
    else
        return(exp(spde$param.inla$B0[,1, drop=TRUE]))
}

inla.spde3.theta2phi1 = function(spde, theta)
{
    inla.require.inherits(spde, "inla.spde3", "'spde'")

    if (spde$n.theta>0)
        return(exp(spde$param.inla$B1[,1, drop=TRUE] +
                   spde$param.inla$B1[,-1, drop=FALSE] %*% theta))
    else
        return(exp(spde$param.inla$B1[,1, drop=TRUE]))
}

inla.spde3.theta2phi2 = function(spde, theta)
{
    inla.require.inherits(spde, "inla.spde3", "'spde'")

    if (spde$n.theta>0) {
        phi = ((spde$param.inla$B2[,1, drop=TRUE] +
                spde$param.inla$B2[,-1, drop=FALSE] %*% theta))
    } else {
        phi = ((spde$param.inla$B2[,1, drop=TRUE]))
    }
    if (spde$param.inla$transform == "identity") {
        return(phi)
    } else if (spde$param.inla$transform == "log") {
        return(exp(phi))
    } else if (spde$param.inla$transform == "logit") {
        return((1-exp(phi))/(1+exp(phi)))
    } else if (spde$param.inla$transform == "shiftedlog") {
        return(2*exp(phi)-1)
    } else if (spde$param.inla$transform == "oldlogit") {
        return(cos(pi/(1+exp(-phi))))
    } else {
        warning(paste("Unknown link function '",
                      spde$param.inla$transform,
                      "' phi2.  Using identity link instead.", sep=""))
        return(phi)
    }
}

inla.spde3.theta2phi3 = function(spde, theta)
{
    inla.require.inherits(spde, "inla.spde3", "'spde'")

    if (spde$n.theta>0)
        return(exp(spde$param.inla$B3[,1, drop=TRUE] +
                   spde$param.inla$B3[,-1, drop=FALSE] %*% theta))
    else
        return(exp(spde$param.inla$B3[,1, drop=TRUE]))
}

inla.spde3.precision =
    function(spde, theta=NULL,
             phi0=inla.spde3.theta2phi0(spde, theta),
             phi1=inla.spde3.theta2phi1(spde, theta),
             phi2=inla.spde3.theta2phi2(spde, theta),
             phi3=inla.spde3.theta2phi3(spde, theta),
             ...)
{
    inla.require.inherits(spde, "inla.spde3", "'spde'")
    if (spde$f$model!="spde3") {
        stop("'inla.spde3.precision' only supports internal inla models 'spde3'")
    }

    D0 <- Diagonal(spde$n.spde, phi0)
    D1 <- Diagonal(spde$n.spde, phi1)
    D12 <- Diagonal(spde$n.spde, phi1*phi2)
    D3 <- Diagonal(nrow(spde$param.inla$M3), phi3)
    Q <- (D0 %*% (D1 %*% spde$param.inla$M0 %*% D1 +
                  D12 %*% spde$param.inla$M1 +
                  t(spde$param.inla$M1) %*% D12 +
                  spde$param.inla$M2 +
                  t(spde$param.inla$M3) %*% D3 %*% spde$param.inla$M3)
          %*% D0)

    return(Q)
}



inla.spde3.result = function(inla, name, spde, do.transform=TRUE, ...)
{
    inla.require.inherits(inla, "inla", "'inla'")
    inla.require.inherits(spde, "inla.spde3", "'spde'")
    if (!spde$f$model=="spde3") {
        stop("'inla.spde3.result' only supports internal inla models 'spde3'")
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
    if (!is.null(inla$summary.spde3.blc[[name]])) {
        rownames(inla$summary.spde3.blc[[name]]) <- BLC.names

        result$summary.theta =
            inla.extract.el(inla$summary.spde3.blc[[name]],
                            "theta\\.[^ ]+$")
        result$summary.log.tau =
            inla.extract.el(inla$summary.spde3.blc[[name]],
                            "tau\\.[^ ]+$")
        result$summary.log.kappa =
            inla.extract.el(inla$summary.spde3.blc[[name]],
                            "kappa\\.[^ ]+$")
        result$summary.log.variance.nominal =
            inla.extract.el(inla$summary.spde3.blc[[name]],
                            "variance.nominal\\.[^ ]+$")
        result$summary.log.range.nominal =
            inla.extract.el(inla$summary.spde3.blc[[name]],
                            "range.nominal\\.[^ ]+$")
    }

    ## Marginals for log-tau/kappa/variance/range
    if (!is.null(inla$marginals.spde3.blc[[name]])) {
        names(inla$marginals.spde3.blc[[name]]) <- BLC.names

        result$marginals.theta =
            inla.extract.el(inla$marginals.spde3.blc[[name]],
                            "theta\\.[^ ]+$")
        result$marginals.log.tau =
            inla.extract.el(inla$marginals.spde3.blc[[name]],
                            "tau\\.[^ ]+$")
        result$marginals.log.kappa =
            inla.extract.el(inla$marginals.spde3.blc[[name]],
                            "kappa\\.[^ ]+$")
        result$marginals.log.variance.nominal =
            inla.extract.el(inla$marginals.spde3.blc[[name]],
                            "variance.nominal\\.[^ ]+$")
        result$marginals.log.range.nominal =
            inla.extract.el(inla$marginals.spde3.blc[[name]],
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

inla.spde3.result = function(inla, name, spde, do.transform=TRUE, ...)
{
    inla.require.inherits(inla, "inla", "'inla'")
    inla.require.inherits(spde, "inla.spde3", "'spde'")
    if (!spde$f$model=="spde3") {
        stop("'inla.spde3.result' only supports internal inla models 'spde3'")
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
    if (!is.null(inla$summary.spde3.blc[[name]])) {
        rownames(inla$summary.spde3.blc[[name]]) <- BLC.names

        result$summary.theta =
            inla.extract.el(inla$summary.spde3.blc[[name]],
                            "theta\\.[^ ]+$")
        result$summary.log.tau =
            inla.extract.el(inla$summary.spde3.blc[[name]],
                            "tau\\.[^ ]+$")
        result$summary.log.kappa =
            inla.extract.el(inla$summary.spde3.blc[[name]],
                            "kappa\\.[^ ]+$")
        result$summary.log.variance.nominal =
            inla.extract.el(inla$summary.spde3.blc[[name]],
                            "variance.nominal\\.[^ ]+$")
        result$summary.log.range.nominal =
            inla.extract.el(inla$summary.spde3.blc[[name]],
                            "range.nominal\\.[^ ]+$")
    }

    ## Marginals for log-tau/kappa/variance/range
    if (!is.null(inla$marginals.spde3.blc[[name]])) {
        names(inla$marginals.spde3.blc[[name]]) <- BLC.names

        result$marginals.theta =
            inla.extract.el(inla$marginals.spde3.blc[[name]],
                            "theta\\.[^ ]+$")
        result$marginals.log.tau =
            inla.extract.el(inla$marginals.spde3.blc[[name]],
                            "tau\\.[^ ]+$")
        result$marginals.log.kappa =
            inla.extract.el(inla$marginals.spde3.blc[[name]],
                            "kappa\\.[^ ]+$")
        result$marginals.log.variance.nominal =
            inla.extract.el(inla$marginals.spde3.blc[[name]],
                            "variance.nominal\\.[^ ]+$")
        result$marginals.log.range.nominal =
            inla.extract.el(inla$marginals.spde3.blc[[name]],
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



inla.spde3.models = function()
{
    return(c("generic", "matern"))
}


## spde.common-connections:
inla.spde.precision.inla.spde3 = inla.spde3.precision
inla.spde.result.inla.spde3 = inla.spde3.result
