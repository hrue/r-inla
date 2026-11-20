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

`inla.qsample` <- function(
                           n = 1L,
                           Q,
                           b,
                           mu,
                           sample,
                           constr,
                           reordering = INLA::inla.reorderings(),
                           seed = 0L,
                           logdens = !(missing(sample)),
                           compute.mean = !(missing(sample)),
                           num.threads = NULL, 
                           selection = NULL,
                           verbose = inla.getOption("verbose"),
                           .debug = FALSE) 
{
    t.dir <- inla.tempdir()
    smtp <- match.arg(inla.getOption("smtp"), c("stiles", "taucs", "band", "default"))
    if (smtp == "stiles") smtp <- "taucs"
    stopifnot(!missing(Q))
    stopifnot(n >= 1L)

    if (is.null(num.threads)) {
        num.threads <- inla.getOption("num.threads")
    }
    if (seed == 0) {
        num.threads <- inla.parse.num.threads(num.threads)
    }

    if (is.list(reordering)) {
        ## argument is the output from inla.qreordering()
        reordering <- reordering$name
    }
    reordering <- match.arg(reordering)

    Q <- inla.sparse.check(Q)
    if (is(Q, "dgTMatrix")) {
        Q.file <- inla.write.fmesher.file(Q, filename = inla.tempfile(tmpdir = t.dir))
    } else if (is.character(Q)) {
        Q.file <- Q
    } else {
        stop("This should not happen.")
    }

    b.file <- inla.tempfile(tmpdir = t.dir)
    mu.file <- inla.tempfile(tmpdir = t.dir)
    constr.file <- inla.tempfile(tmpdir = t.dir)
    x.file <- inla.tempfile(tmpdir = t.dir)
    sample.file <- inla.tempfile(tmpdir = t.dir)
    rng.file <- inla.tempfile(tmpdir = t.dir)
    cmean.file <- inla.tempfile(tmpdir = t.dir)
    selection.file <- inla.tempfile(tmpdir = t.dir)

    if (!missing(b)) {
        stopifnot(length(b) == nrow(Q))
        b <- matrix(b, nrow(Q), 1)
        inla.write.fmesher.file(b, filename = b.file)
    }

    if (!missing(mu)) {
        stopifnot(length(mu) == nrow(Q))
        mu <- matrix(mu, nrow(Q), 1)
        inla.write.fmesher.file(mu, filename = mu.file)
    }

    if (!missing(constr) && !is.null(constr)) {
        stopifnot(is.list(constr))
        A <- as.matrix(constr$A)
        e <- as.numeric(constr$e)
        stopifnot(ncol(A) == ncol(Q))
        stopifnot(nrow(A) == length(e))
        xx <- matrix(c(nrow(A), c(A), c(e)), ncol = 1)
        inla.write.fmesher.file(xx, filename = constr.file)
    }

    if (!missing(sample) && !is.null(sample)) {
        sample <- as.matrix(sample)
        stopifnot(nrow(sample) == nrow(Q))
        stopifnot(ncol(sample) > 0L)
        inla.write.fmesher.file(sample, filename = sample.file)
        n <- ncol(sample) ## redefine n here
    }

    if (!missing(selection) && !is.null(selection)) {
        if (!missing(sample)) {
            stop("Cannot use 'selection' and 'sample' at the same time")
        }
        selection <- as.matrix(selection, ncol = 1)
        stopifnot(nrow(selection) <= nrow(Q))
        inla.write.fmesher.file(selection - 1, filename = selection.file)
    } else {
        ## make the code easier below
        selection <- seq_len(nrow(Q))
    }

    envir <- inla.get.inlaEnv()
    if (seed < 0L) {
        if (!exists("GMRFLib.rng.state", envir = envir)) {
            seed <- 0L
        } else {
            rng.state <- get("GMRFLib.rng.state", envir = envir)
            fp <- file(rng.file, "wb")
            writeBin(as.raw(rng.state), fp)
            close(fp)
        }
    }

    inla.set.environment()
    env <- inla.run.environment.set()
    if (inla.os("linux") || inla.os("mac") || inla.os("mac.arm64")) {
        if (.debug) {
            print(paste(
                shQuote(inla.call.no.remote()), "-s -m qsample",
                paste0("-t", num.threads), "-r", reordering, "-z", seed, "-S", smtp,
                if (verbose) "-v" else "",
                Q.file, x.file, as.integer(n), rng.file,
                sample.file, b.file, mu.file, constr.file, cmean.file, selection.file
            ))
        }
        s <- system(paste(
            shQuote(inla.call.no.remote()), "-s -m qsample",
            paste0("-t", num.threads), "-r", reordering, "-z", seed, "-S", smtp,
            if (verbose) "-v" else "",
            Q.file, x.file, as.integer(n), rng.file,
            sample.file, b.file, mu.file, constr.file, cmean.file, selection.file
        ), intern = TRUE)
    } else if (inla.os("windows")) {
        s <- system(paste(
            shQuote(inla.call.no.remote()), "-s -m qsample",
            paste0("-t", num.threads), "-r", reordering, "-z", seed, "-S", smtp,
            if (verbose) "-v" else "",
            Q.file, x.file, as.integer(n), rng.file,
            sample.file, b.file, mu.file, constr.file, cmean.file, selection.file
        ), intern = TRUE)
    } else {
        stop("\n\tNot supported architecture.")
    }

    inla.run.environment.unset(env)
    fp <- file(rng.file, "rb")
    siz <- file.info(rng.file)$size
    rng.state <- readBin(fp, raw(), siz)
    close(fp)
    assign("GMRFLib.rng.state", rng.state, envir = envir)

    x <- inla.read.fmesher.file(x.file)
    cmean <- inla.read.fmesher.file(cmean.file)

    nx <- dim(x)[1L] - 1L
    samples <- matrix(x[-(nx + 1L), , drop = FALSE], nx, n)
    colnames(samples) <- paste0("sample:", 1L:n)
    rownames(samples) <- paste0("x:", selection)
    ld <- c(x[nx + 1L, ])
    names(ld) <- paste("logdens", 1L:n, sep = "")

    unlink(t.dir, recursive = TRUE)

    if (logdens || compute.mean) {
        return(list(sample = samples, logdens = ld, mean = cmean))
    } else {
        return(samples)
    }
}

`inla.qsample.stiles` <- function(n = 1L, Q, b = NULL, mu = NULL, sample = NULL,
                           constr = NULL, logdens = !is.null(sample), num.threads = NULL,
                           seed = 0L, selection = list(), verbose = FALSE, ...)
{
    Q <- inla.as.sparse(Q)
    if (all(Q@i >= Q@j) || all(Q@i <= Q@j)) {
        Q <- Q + t(Q)
        diag(Q) <- diag(Q) / 2
    }

    if (!(is.list(selection) && length(selection) == 0)) {
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
        cores <- numeric(2)
        cores[1] <- max(1, as.integer(strsplit(inla.parse.num.threads(num.threads), ":")[[1]][1]))
        cores[2] <- max(1, as.integer(strsplit(inla.parse.num.threads(num.threads), ":")[[1]][2]))
        cores <- prod(cores)
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

