#' @title Utility functions for the PC prior for a correlation matrix
#' 
#' @description
#' Functions to evaluate and sample from the PC prior for a correlation matrix.
#' 
#' The parameterisation of a correlation matrix of dimension `p` has
#' `dim` parameters: `theta` which are in the interval -pi to pi.
#' The alternative parameterisation is through the off-diagonal elements
#' `r` of the correlation matrix `R`. The functions
#' `inla.pc.cormat.<A>2<B>` convert between parameterisations `<A>`
#' to parameterisations `<B>`, where both `<A>` and `<B>` are
#' one of `theta`, `r` and `R`, and `p` and `dim`.
#' 
#' @aliases inla.pc.cormat inla.pc.cormat.dim2p cormat.dim2p
#' inla.pc.cormat.p2dim cormat.p2dim inla.pc.cormat.theta2R cormat.theta2R
#' inla.pc.cormat.R2theta cormat.R2theta inla.pc.cormat.r2R cormat.r2R
#' inla.pc.cormat.R2r cormat.R2r inla.pc.cormat.r2theta cormat.r2theta
#' inla.pc.cormat.theta2r cormat.theta2r inla.pc.cormat.permute cormat.permute
#' inla.pc.cormat.rtheta cormat.rtheta inla.pc.cormat.dtheta cormat.dtheta
#' @param dim The dimension of `theta`, the parameterisatin of the
#' correlation matrix
#' @param p The dimension the correlation matrix
#' @param theta A vector of parameters for the correlation matrix
#' @param r The off diagonal elements of a correlation matrix
#' @param R A correlation matrix
#' @param n Number of observations
#' @param lambda The rate parameter in the prior
#' @param log Logical. Return the density in natural or log-scale.
#' @returns `inla.pc.cormat.rtheta` generate samples from the prior,
#' returning a matrix where each row is a sample of `theta`.
#' `inla.pc.cormat.dtheta` evaluates the density of `theta`.
#' `inla.pc.cormat.permute` randomly permutes a correlation matrix, which
#' is useful if an exchangable sample of a correlation matrix is required.
#' @author Havard Rue \email{hrue@@r-inla.org}
#' @examples
#' 
#'   p = 4
#'   print(paste("theta has length", inla.pc.cormat.p2dim(p)))
#'   theta = inla.pc.cormat.rtheta(n=1, p=4, lambda = 1)
#'   print("sample theta:")
#'   print(theta)
#'   print(paste("log.dens", inla.pc.cormat.dtheta(theta, log=TRUE)))
#'   print("r:")
#'   r = inla.pc.cormat.theta2r(theta)
#'   print(r)
#'   print("A sample from the non-exchangable prior, R:")
#'   R = inla.pc.cormat.r2R(r)
#'   print(R)
#'   print("A sample from the exchangable prior, R:")
#'   R = inla.pc.cormat.permute(R)
#'   print(R)
#'  
#' @name pc.cormat
#' @rdname pc-cormat
NULL



#' @rdname pc-cormat
#' @export
inla.pc.cormat.dim2p <- function(dim) {
    p <- round(1 / 2 + 1 / 2 * sqrt(1 + 4L * dim * 2L))
    stopifnot(abs(dim - inla.pc.cormat.p2dim(p)) < sqrt(.Machine$double.eps))
    return(p)
}

#' @rdname pc-cormat
#' @export
inla.pc.cormat.p2dim <- function(p) {
    return(p * (p - 1L) / 2L)
}

#' @rdname pc-cormat
#' @export
inla.pc.cormat.theta2R <- function(theta) {
    p <- inla.pc.cormat.dim2p(length(theta))
    theta.m <- matrix(NA, p, p)
    theta.m[lower.tri(theta.m)] <- theta

    B <- matrix(0, p, p)
    B[1, 1] <- 1.0
    for (i in 2L:p) {
        for (j in 1L:i) {
            if (j == 1) {
                B[i, j] <- cos(theta.m[i, j])
            } else if (j >= 2L && j <= i - 1) {
                B[i, j] <- cos(theta.m[i, j]) * prod(sin(theta.m[i, 1:(j - 1)]))
            } else if (j == i) {
                B[i, j] <- prod(sin(theta.m[i, 1:(j - 1)]))
            } else {
                stop("This should not happen.")
            }
        }
    }
    R <- B %*% t(B)
    diag(R) <- 1.0 ## so that its exactly 1, not just numerically 1

    return(R)
}

#' @rdname pc-cormat
#' @export
inla.pc.cormat.R2theta <- function(R) {
    L <- t(chol(R))
    p <- dim(R)[1]

    ## overwrite L with theta
    for (i in 2:p) {
        if (i == 2) {
            L[2, 1] <- acos(L[2, 1])
        } else {
            for (j in 1:(i - 1)) {
                if (j == 1) {
                    L[i, j] <- acos(L[i, j])
                } else {
                    L[i, j] <- acos(L[i, j] / prod(sin(L[i, 1:(j - 1)])))
                }
            }
        }
    }
    theta <- L[lower.tri(L)]

    return(theta)
}

#' @rdname pc-cormat
#' @export
inla.pc.cormat.r2R <- function(r) {
    p <- inla.pc.cormat.dim2p(length(r))
    R <- matrix(1, p, p)
    R[lower.tri(R)] <- r
    R <- t(R)
    R[lower.tri(R)] <- r
    diag(R) <- 1.0

    return(R)
}

#' @rdname pc-cormat
#' @export
inla.pc.cormat.R2r <- function(R) {
    return(R[lower.tri(R)])
}

#' @rdname pc-cormat
#' @export
inla.pc.cormat.r2theta <- function(r) {
    R <- inla.pc.cormat.r2R(r)
    theta <- inla.pc.cormat.R2theta(R)

    return(theta)
}

#' @rdname pc-cormat
#' @export
inla.pc.cormat.theta2r <- function(theta) {
    R <- inla.pc.cormat.theta2R(theta)
    r <- inla.pc.cormat.R2r(R)

    return(r)
}

#' @rdname pc-cormat
#' @export
inla.pc.cormat.permute <- function(R) {
    r <- inla.pc.cormat.R2r(R)
    r <- r[order(runif(length(r)))]
    R <- inla.pc.cormat.r2R(r)

    return(R)
}

#' @rdname pc-cormat
#' @export
inla.pc.cormat.rtheta <- function(n = 1, p, lambda = 1) {
    stopifnot(!missing(p) && p > 1)
    stopifnot(lambda > 0)
    stopifnot(n >= 1)

    rcormat <- function(n, p, lambda) {
        q <- inla.pc.cormat.p2dim(p)
        ## sample a point on the simplex
        r <- rexp(1, rate = lambda)
        gamma <- rexp(q)
        gamma <- (gamma / sum(gamma)) * r^2 / 2
        theta <- asin(exp(-gamma))
        ## sample branch by random
        branch <- sample(c(0, 1), size = q, replace = TRUE)
        theta <- branch * theta + (1 - branch) * (pi - theta)
        return(theta)
    }

    x <- unlist(lapply(rep(1, n), rcormat, p = p, lambda = lambda))
    x <- matrix(x, n, inla.pc.cormat.p2dim(p), byrow = TRUE)
    return(x)
}

#' @rdname pc-cormat
#' @export
inla.pc.cormat.dtheta <- function(theta, lambda = 1, log = FALSE) {
    ## reimplementation using the simplex function
    p <- length(theta)
    gamma <- -log(sin(theta))
    ldens <- (inla.pc.multvar.simplex.d(gamma, lambda = lambda, log = TRUE, b = rep(1, p))
    + sum(log(abs(1 / tan(theta)))) - p * log(2))

    return(if (log) ldens else exp(ldens))
}
