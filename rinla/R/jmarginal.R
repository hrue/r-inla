#' Sample, transform and evaluate from a joint marginal approximation
#' 
#' Sample, transform and evalue from from a joint marginal approximation as
#' returned using argument \code{selection} in \code{inla}.
#' 
#' 
#' @aliases inla.joint.marginal inla.joint.marginal.eval joint.marginal
#' joint.marginal.eval rjmarginal rjmarginal.eval inla.rjmarginal
#' inla.rjmarginal.eval inla.tjmarginal tjmarginal inla.1djmarginal 1djmarginal
#' @param n The number of samples
#' @param jmarginal A marginal object given either by a \code{inla} object or
#' \code{result$selection}
#' @param constr Optional linear constraints; see \code{?INLA::f} and argument
#' \code{extraconstr}
#' @param fun A function which is evaluated for each sample, similar to
#' \code{inla.posterior.sample.eval}: please see the documentation for this
#' functions for details.
#' @param samples The samples, as in the form of the output from
#' \code{inla.rjmarginal}
#' @param A A matrix used for the linear combination
#' @param ... Arguments passed on to other methods (printing and summarising)
#' @returns THESE FUNCTIONS ARE EXPERIMENTAL FOR THE MOMENT (JULY 2020)
#' 
#' \code{inla.rjmarginal} returns a list with the samples in \code{samples}
#' (matrix) and the corresponding log-densities in \code{log.density} (vector).
#' Each column in \code{samples} contains one sample.
#' 
#' \code{inla.rjmarginal.eval} returns a matrix, where each row is the (vector)
#' function evaluated at each sample.
#' 
#' \code{inla.tjmarginal} returns a \code{inla.jmarginal}-object of the linear
#' combination defined by the matrix \code{A}.
#' 
#' \code{inla.1djmarginal} return the marginal densities from a joint
#' approximation.
#' @author Cristian Chiuchiolo and Havard Rue \email{hrue@@r-inla.org}
#' @seealso \code{\link{inla}}
#' @examples
#' 
#'  n = 10
#'  x = 1+rnorm(n)
#'  xx = 3 + rnorm(n)
#'  y = 1 + x + xx + rnorm(n)
#'  selection = list(xx=1, Predictor = 3:4, x=1)
#'  r = inla(y ~ 1 + x + xx,
#'           data = data.frame(y, x, xx),
#'           selection = selection)
#'  ns = 100
#'  xx = inla.rjmarginal(ns, r)
#' 
#'  print(cbind(mean = r$selection$mean, sample.mean = rowMeans(xx$samples)))
#'  print("cov matrix")
#'  print(round(r$selection$cov.matrix, dig=3))
#'  print("sample cov matrix")
#'  print(round(cov(t(xx$samples)), dig=3))
#' 
#'  skew = function(z) mean((z-mean(z))^3)/var(z)^1.5
#'  print(round(cbind(skew = r$selection$skewness,
#'                    sample.skew = apply(xx$sample, 1, skew)), dig=3))
#' 
#'  ## illustrating the eval function
#'  n = 10
#'  x = rnorm(n)
#'  eta = 1 + x
#'  y = eta + rnorm(n, sd=0.1)
#'  selection = list(x = 1, Predictor = c(1, 2, 4, 5),  '(Intercept)' = 1)
#'  r = inla(y ~ 1 + x,
#'           data = data.frame(y, x),
#'           selection = selection)
#'  xx = inla.rjmarginal(100,  r)
#'  xx.eval = inla.rjmarginal.eval(function() c(x, Predictor, Intercept),  xx)
#'  print(cbind(xx$samples[, 1]))
#'  print(cbind(xx.eval[, 1]))
#' 
#'  constr <- list(A = matrix(1, ncol = n, nrow = 1), e = 1)
#'  x <- inla.rjmarginal(10, r, constr = constr)
#' 
#'  A <- matrix(rnorm(n^2), n, n)
#'  b <- inla.tjmarginal(r, A)
#'  b.marg <- inla.1djmarginal(b)
#' 
#' @name joint.marginal
#' @rdname jmarginal
#' @export

`inla.rjmarginal` <- function(n, jmarginal, constr) {
    if (missing(jmarginal) || missing(n) || n <= 0) {
        return(list(samples = matrix(ncol = 0, nrow = 0), log.density = numeric(0)))
    }

    if (inherits(jmarginal, "inla")) {
        jmarginal <- jmarginal$selection
    } else if (inherits(jmarginal, "inla.jmarginal")) {
        ## ok
    } else {
        stop("Unknown object: argument 'jmarginal'")
    }

    ## the setup depends on 'constr' or not
    cov.matrix <- inla.ensure.spd(jmarginal$cov.matrix)
    if (missing(constr)) {
        constr <- NULL
        Q <- solve(cov.matrix)
        Qi <- cov.matrix
        mm <- jmarginal$mean
    } else {
        constr$nc <- nrow(constr$A) ## this is required in the 'configs'
        Q <- solve(cov.matrix)
        Qi <- cov.matrix
        QiA <- Qi %*% t(constr$A)
        SS <- solve(constr$A %*% QiA)
        Qi <- Qi - QiA %*% SS %*% t(QiA)
        mm <- jmarginal$mean - QiA %*% SS %*% (constr$A %*% jmarginal$mean - constr$e)
    }

    ## build a fake 'inla'-object, so we can feed that into 'inla.posterior.sample'
    m <- length(jmarginal$mean)
    r <- list(
        misc = list(
            configs = list(
                n = m,
                nconfig = 1,
                max.log.posterior = 0,
                nz = m * (m + 1) / 2,
                contents = list(
                    tag = jmarginal$names,
                    start = 1:m,
                    length = rep(1, m)
                ),
                ntheta = 0,
                constr = constr,
                config = list(list(
                    theta = NULL,
                    Q = inla.as.sparse(Q),
                    Qinv = inla.as.sparse(Qi),
                    mean = mm,
                    improved.mean = mm,
                    skewness = jmarginal$skewness,
                    log.posterior = 0,
                    log.posterior.orig = 0
                ))
            )
        )
    )
    class(r) <- "inla"

    x <- inla.posterior.sample(n, r,
                               use.improved.mean = TRUE, skew.corr = TRUE,
                               add.names = FALSE
                               )
    xx <- matrix(unlist(lapply(x, function(z) z$latent)), ncol = n)
    log.dens <- unlist(lapply(x, function(z) z$logdens$joint))
    rownames(xx) <- jmarginal$names
    colnames(xx) <- paste0("sample:", 1:n)
    names(log.dens) <- colnames(xx)

    return(list(samples = xx, log.density = log.dens))
}

#' @rdname jmarginal
#' @export
`inla.rjmarginal.eval` <- function(fun, samples, ...) {
    stopifnot(all(names(samples) == c("samples", "log.density")))

    var <- "inla.jmarginal.eval.warning.given"
    if (!(exists(var, envir = inla.get.inlaEnv()) &&
          get(var, envir = inla.get.inlaEnv()) == TRUE)) {
        warning("Function 'inla.rjmarginal.eval()' is experimental.")
        assign(var, TRUE, envir = inla.get.inlaEnv())
    }

    nm.split <- strsplit(rownames(samples$samples), ":")
    ## determine the names
    nms <- unique(unlist(lapply(nm.split, function(x) x[1])))
    contents <- rep(list(list()), length(nms))
    for (i in seq_along(nms)) {
        contents[[i]]$tag <- nms[i]
        contents[[i]]$start <- min(which(nms[i] == lapply(nm.split, function(x) x[1])))
        contents[[i]]$sample.idx <- which(nms[i] == lapply(nm.split, function(x) x[1]))
    }

    my.fun <- function(a.sample, .contents, .fun, ...) {
        env <- new.env()
        for (i in seq_along(.contents)) {
            .tmp <- a.sample[contents[[i]]$sample.idx]
            assign(.contents[[i]]$tag, .tmp, envir = env)
        }
        if (exists("(Intercept)", envir = env)) {
            assign("Intercept", get("(Intercept)", envir = env),
                   envir = env
                   )
        }
        parent.env(env) <- .GlobalEnv
        environment(.fun) <- env
        return(.fun(...))
    }

    ret <- apply(samples$samples, 2, my.fun, .fun = fun, .contents = contents, ...)
    ret <- matrix(ret, ncol = ncol(samples$samples))
    colnames(ret) <- paste0("sample:", 1:ncol(ret))
    rownames(ret) <- paste0("fun[", 1:nrow(ret), "]")

    return(ret)
}

#' @rdname jmarginal
#' @method print inla.jmarginal
#' @export
`print.inla.jmarginal` <- function(x, ...) {
    x$.private <- NULL
    class(x) <- class(list())
    print(x)
}

#' @rdname jmarginal
#' @param object Object to be summarised
#' @method summary inla.jmarginal
#' @export
`summary.inla.jmarginal` <- function(object, ...) {
    inla.require("sn", stop.on.error = TRUE)
    mode.sn <- function(xi, omega, alpha) {
        med <- sn::qsn(0.5, xi, omega, alpha)
        res <- optimize(
            f = dsn, interval = c(med - omega, med + omega), maximum = TRUE,
            ## arguments to 'dsn'
            log = TRUE, xi = xi, omega = omega, alpha = alpha
        )
        return(res$maximum)
    }

    n.sel <- object$names
    mu <- object$mean
    std <- sqrt(diag(object$cov.matrix))
    dsn.xi <- object$marginal.sn.par$xi
    dsn.omega <- object$marginal.sn.par$omega
    dsn.alpha <- object$marginal.sn.par$alpha
    prob <- c(0.025, 0.5, 0.975)
    mode <- c()
    qsn.eval <- matrix(NA, nrow = length(n.sel), ncol = length(prob))
    for (i in seq_along(n.sel)) {
        qsn.eval[i, ] <- sn::qsn(prob, xi = dsn.xi[i], omega = dsn.omega[i], alpha = dsn.alpha[i])
        mode[i] <- mode.sn(xi = dsn.xi[i], omega = dsn.omega[i], alpha = dsn.alpha[i])
    }
    obj <- cbind(mu, std, qsn.eval, mode)
    rownames(obj) <- n.sel
    colnames(obj) <- c("mean", "sd", "0.025quant", "0.5quant", "0.975quant", "mode")
    ret <- list(msg = "Joint marginal is computed for: ", matrix = obj)
    class(ret) <- "summary.inla.jmarginal"

    return(ret)
}

#' @rdname jmarginal
#' @param x Object to be printed
#' @method print summary.inla.jmarginal
#' @export
`print.summary.inla.jmarginal` <- function(x, ...) {
    print(x$msg)
    print(as.matrix(x$matrix))
    return(invisible())
}

#' @rdname jmarginal
#' @export
`inla.tjmarginal` <- function(jmarginal, A) {
    if (inherits(jmarginal, "inla")) {
        jmarginal <- jmarginal$selection
    } else if (inherits(jmarginal, "inla.jmarginal")) {
        ## ok
    } else {
        stop("Unknown object: argument 'jmarginal'")
    }
    stopifnot(is.matrix(A))

    sel.len <- length(jmarginal$names)
    if (ncol(A) != sel.len){
        stop("Incorrect dimension for the selected components")
    }

    if (is.null(rownames(A))) {
        names.sel <- sapply(1:nrow(A), function(x) paste0("Lin:", x))
    } else {
        names.sel <- rownames(A)
    }
    
    skew.max <- 0.99
    moments <- jmarginal$.private$moments
    skewness.sel <- jmarginal$skewness
    mom1 <- moments[[1]]
    mom2 <- moments[[2]]
    mom3 <- moments[[3]]
    S <- jmarginal$cov.matrix
    var.sel <- diag(S)
    mu.tjoint <- A %*% mom1
    S.tjoint <- A %*% S %*% t(A)
    m2.tjoint <- diag(S.tjoint) + mu.tjoint^2
    m3.tjoint <- skew.tjoint <- c()
    for (lc in 1:nrow(A)) {
        x <- A[lc, ]
        lc.ind <- which(x != 0 & !is.na(x))
        coef <- A[lc, lc.ind]
        skewness.lin <- skewness.sel[lc.ind]
        var.lin <- var.sel[lc.ind]
        mom1.m <- mu.tjoint[lc]
        mom2.m <- m2.tjoint[lc]
        mom3.cm <- sum(coef^3*skewness.lin*var.lin^(1.5)) 
        skew.m <- (mom3.cm)*((mom2.m - mom1.m^2)^(-1.5)) # global skewnesses
        m3.tjoint[lc] <- mom3.cm+3*mom1.m*mom2.m-2*mom1.m^3
        if (any(abs(skew.m) > skew.max)) {
            skew.m <- pmax(-skew.max, pmin(skew.max, skew.m))
            warning(paste0("One or more marginal skewness are too high. Coerced to be ", skew.max))
        }
        skew.tjoint[lc] <- skew.m
    }
    sn.par <- inla.sn.reparam(moments = list(mean = as.numeric(mu.tjoint),
                                             variance = diag(S.tjoint),
                                             skewness = skew.tjoint))
    output <- list()
    output$names <- names.sel
    output$mean <- mu.tjoint
    output$cov.matrix <- S.tjoint
    output$skewness <- skew.tjoint
    output$marginal.sn.par$xi <- sn.par$xi
    output$marginal.sn.par$omega <- sn.par$omega
    output$marginal.sn.par$alpha <- sn.par$alpha
    output$.private$moments <- list(as.numeric(mu.tjoint), as.numeric(m2.tjoint), m3.tjoint)
    names(output$.private$moments[[1]]) <- names.sel
    names(output$.private$moments[[2]]) <- names.sel
    names(output$.private$moments[[3]]) <- names.sel
    class(output) <- "inla.jmarginal"

    return(output)
}

#' @rdname jmarginal
#' @export
`inla.1djmarginal` <- function(jmarginal) {
    inla.require("sn", stop.on.error = TRUE)
    if (inherits(jmarginal, "inla")) {
        jmarginal <- jmarginal$selection
    } else if (inherits(jmarginal, "inla.jmarginal")) {
        ## ok
    } else {
        stop("Unknown object: argument 'jmarginal'")
    }

    n.sel <- jmarginal$names
    mu <- jmarginal$mean
    std <- sqrt(diag(jmarginal$cov.matrix))
    dsn.xi <- jmarginal$marginal.sn.par$xi
    dsn.omega <- jmarginal$marginal.sn.par$omega
    dsn.alpha <- jmarginal$marginal.sn.par$alpha
    obj <- vector("list", length(n.sel))
    names(obj) <- n.sel

    ## copy from density.c
    q.many <- c(
        0.0000001, 0.000001, 0.00001, 0.0001, 0.0005, 0.001, 0.005, 0.01, 0.025,
        0.05, 0.075, 0.10, 0.125, 0.15, 0.175, 0.2, 0.225, 0.25, 0.275, 0.30,
        0.325, 0.35, 0.375, 0.40, 0.425, 0.45, 0.46, 0.47, 0.475, 0.48, 0.49,
        0.50, 0.51, 0.52, 0.525, 0.53, 0.54, 0.55, 0.575, 0.60, 0.625, 0.65,
        0.675, 0.70, 0.725, 0.75, 0.775, 0.80, 0.825, 0.85, 0.875, 0.9, 0.925,
        0.95, 0.975, 0.99, 0.995, 0.999, 0.9995, 0.9999, 0.99999, 0.999999, 0.9999999
    )

    for (i in seq_along(n.sel)) {
        val <- sn::qsn(q.many, xi = dsn.xi[i], omega = dsn.omega[i], alpha = dsn.alpha[i])
        dsn.eval <- sn::dsn(val, xi = dsn.xi[i], omega = dsn.omega[i], alpha = dsn.alpha[i])
        obj[[i]] <- cbind(val, dsn.eval)
        colnames(obj[[i]]) <- c("x", "y")
    }

    return(obj)
}
