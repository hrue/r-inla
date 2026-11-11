#' Generate samples from a GMRF using the GMRFLib implementation
#' 
#' This function generate samples, or evaluate the log-densities of samples, from a GMRF
#' 
#' 
#' @aliases inla.qsample qsample
#' @param n Number of samples. Only used if `sample` is not `NULL`
#' @param Q The precision matrix or a filename containing it.
#' @param b The linear term
#' @param mu The mu term
#' @param sample A matrix of optional samples where each column is a sample. If set, then use
#'     this instead of generating new samples
#' @param constr Optional linear constraints; see `?INLA::f` and argument `extraconstr`
#' @param logdens If `TRUE`, compute also the log-density of each sample, and the (constrained)
#'     mean. Note that the output format then change.
#' @param num.threads Number of threads to use. If in format "A:B:C", then "A" is used.
#' @param seed Seed for RNG, if non-zero, set seed using [set.seed()]. If `seed !=0` then
#'     argument `num.threads` is set to 1. In that case, the results are almost surely
#'     reproducible.
#' @param selection This argument is no longer in use and will be removed in the future. Will
#'     issue a warning if used
#' @param verbose This argument is no longer in use and will be removed in the future.. Will
#'     issue a waring if used
#' @param ... Additional arguments no longer supported
#' 
#' @return The log-density has form -1/2(x-mu)^T Q (x-mu) + b^T x. If `logdens` is `FALSE`. If
#'     `samples` are not `NULL`, then `samples` are used instead instead of generating new ones.
#'     `inla.qsample` returns the samples in a matrix if `logdens` is `FALSE`, where each column
#'     is a sample. If `logdens` is `TRUE`, then a list with names `sample`, `logdens` and
#'     `mean` is returned. The (constrained) mean is stored in the vector `mean`.
#' 
#' @author Havard Rue \email{hrue@@r-inla.org}
#' @examples
#' 
#' @name qsample
#' @rdname qsample
#' @export
`inla.qsample` <- function(n = 1L, Q, b = NULL, mu = NULL, sample = NULL,
                           constr = NULL, logdens = !is.null(sample), num.threads = 2L,
                           seed = 0L, selection = list(), verbose = FALSE, ...)
{
    Q <- inla.as.sparse(Q)
    if (all(Q@i >= Q@j) || all(Q@i <= Q@j)) {
        Q <- Q + t(Q)
        diag(Q) <- diag(Q) / 2
    }

    if (is.list(selection) && length(selection) == 0) {
        warning("Argument 'selection' in 'inla.qsample()' is no longer supported")
    }
    if (verbose) {
        warning("Argument 'verbose' in 'inla.qsample()' is no longer supported")
    }
    if (seed != 0) {
        set.seed(seed)
        cores <- 1
    } else {
        if (is.null(num.threads)) num.threads <- inla.getOption("num.threads")
        cores <- as.integer(max(1, strsplit(inla.parse.num.threads(num.threads), ":")[[1]][1]))
    }
    
    m <- nrow(Q)
    s <- sTiles::sTiles(Q, cores = cores)
    new.samples <- if (is.null(sample)) TRUE else FALSE
    if (is.null(sample)) {
        X <- matrix(rnorm(n*m), m, n)
        X <- sTiles::sTiles_solve(s, X, system = "Lt") 
    }

    lin.term <- rep(0, m)
    if (!is.null(b)) {
        lin.term <- lin.term + b
    }
    if (!is.null(mu)) {
        lin.term <- lin.term + Q %*% mu
    }
    u.mean <- sTiles::sTiles_solve(s, lin.term, system = "A")
    if (new.samples) {
        for(i in seq_len(n)) {
            X[, i] <- X[, i] + u.mean
        }
    }
    
    if (!is.null(constr)) {
        A <- constr$A
        e <- constr$e
        nc <- length(e)
        Qinv.tA <- sTiles::sTiles_solve(s, t(A), system = "A")
        A.Qinv.tA <- A %*% Qinv.tA
        SS <- solve(A.Qinv.tA)
        H <- Qinv.tA %*% SS
        if (new.samples) {
            for(i in seq_len(nc)) {
                X[, i] <- X[, i] - H %*% (A %*% X[, i] -e)
            }
        }
        c.mean <- u.mean - H %*% (A %*% u.mean - e)
    } else {
        c.mean <- u.mean
    }

    if (!new.samples) {
        X <- sample
    }
    
    if (logdens) {
        ld <- inla.dmvnorm(X, Q = Q, mu = u.mean, constr = constr, s = s, cores = cores)
        r <- list(sample = X, logdens = ld, mean = c.mean)
    } else {
        r <- X
    }

    if (!is.null(s)) {
        sTiles::sTiles_close(s)
    }

    return (r)
}

`inla.dmvnorm` <- function(X, Q, mu = NULL, constr = NULL, s = NULL, cores = 2L)
{
    ## return the log-density for each column of X, when the GMRF is defined with precision 'Q'
    ## and (unconstrained) mean 'mu'. reuse 's' if provided. if there are constraints, 'mu' is
    ## the unconstrained mean

    no.sTiles.object <- is.null(s)
    if (no.sTiles.object) {
        s <- sTiles::sTiles(Q, cores = cores)
    }
    log.det <- sTiles::sTiles_logdet(s)

    n <- ncol(Q)
    X <- matrix(X, nrow = n)
    m <- ncol(X)
    r <- numeric(m)
    if (is.null(mu)) {
        mu <- rep(0, n)
    }
    for(i in seq_len(m)) {
        XX <- X[, i] - mu
        r[i] <- - n / 2.0 * log(2*pi) + 0.5 * log.det - 0.5 * (t(XX) %*% Q %*% XX)
    }

    if (!is.null(constr)) {
        A <- constr$A
        S.tA <- sTiles::sTiles_solve(s, t(A), system = "A")
        A.S.tA <- A %*% S.tA
        b <- mvtnorm::dmvnorm(t(A %*% X), sigma = A.S.tA, mean = A %*% mu, log = TRUE)
        d <- -0.5 * determinant(A %*% t(A))$modulus
        r <- r + d - b
    }
    
    if (no.sTiles.object) {
        sTiles::sTiles_close(s)
    }
    
    return (r)
}
