#' Generate samples from a GMRF using the GMRFLib implementation
#' 
#' This function generate samples from a GMRF using the GMRFLib implementation
#' 
#' 
#' @aliases inla.qsample qsample
#' @param n Number of samples. Only used if `missing(sample)`
#' @param Q The precision matrix or a filename containing it.
#' @param b The linear term
#' @param mu The mu term
#' @param sample A matrix of optional samples where each column is a sample. If
#' set, then evaluate the log-density for each sample only.
#' @param constr Optional linear constraints; see `?INLA::f` and argument
#' `extraconstr`
#' @param reordering The type of reordering algorithm to be used for
#' `TAUCS`; either one of the names listed in `inla.reorderings()` or
#' the output from `inla.qreordering(Q)`.  The default is "auto" which try
#' several reordering algorithm and use the best one for this particular
#' matrix.
#' @param seed Control the RNG. If `seed=0L` then GMRFLib will set the
#' seed intelligently/at 'random', and this is and should be the default
#' behaviour.  If `seed < 0L` then the saved state of the RNG will be
#' reused if possible, otherwise, GMRFLib will set the seed intelligently/at
#' 'random'.  If `seed > 0L` then this value is used as the seed for the
#' RNG.
#' 
#' PLEASE NOTE1: If `seed!=0` then the computations will run in serial
#' mode, over-riding whatever is set in `num.threads` (a warning might be
#' issued).
#' 
#' PLEASE NOTE2: If the PARDISO sparse matrix library is used, continuity of
#' the samples with respect to small changes in the precision matrix, can be
#' expected but is not guaranteed. If this feature is required, please use the
#' TAUCS sparse matrix library.
#' @param logdens If `TRUE`, compute also the log-density of each sample.
#' Note that the output format then change.
#' @param compute.mean If `TRUE`, compute also the (constrained) mean.
#' Note that the output format then change.
#' @param num.threads Maximum number of threads the `inla`-program will
#' use, or as 'A:B' defining the number threads in the outer (A) and inner (B)
#' layer for nested parallelism.  `seed!=0` requires serial comptuations.
#' @param selection A vector of indices of each sample to return. `NULL`
#' means return the whole sample.  (Note that the log-density retured, is for
#' the whole sample.)  The use of `selection` cannot be combined with the
#' use of `sample`.
#' @param verbose Logical. Run in verbose mode or not.
#' @param .debug Logical. Internal debug-mode.
#' @return The log-density has form -1/2(x-mu)^T Q (x-mu) + b^T x
#' 
#' If `logdens` is `FALSE`, then `inla.qsample` returns the
#' samples in a matrix, where each column is a sample.  If `logdens` or
#' `compute.mean` is `TRUE`, then a list with names `sample`,
#' `logdens` and `mean` is returned. The samples are stored in the
#' matrix `sample` where each column is a sample, and the log densities of
#' each sample are stored as the vector `logdens`.  The mean (include
#' corrections for the constraints, if any) is store in the vector `mean`.
#' @author Havard Rue \email{hrue@@r-inla.org}
#' @examples
#' 
#'  g = system.file("demodata/germany.graph", package="INLA")
#'  Q = inla.graph2matrix(g)
#'  diag(Q) = dim(Q)[1L]
#'  x = inla.qsample(10, Q)
#'  \dontrun{matplot(x)}
#'  x = inla.qsample(10, Q, logdens=TRUE)
#'  \dontrun{matplot(x$sample)}
#' 
#'  n = 3
#'  Q = diag(n)
#'  ns = 2
#' 
#'  ## sample and evaluate a sample
#'  x = inla.qsample(n, Q=Q, logdens=TRUE)
#'  xx = inla.qsample(Q=Q,  sample = x$sample)
#'  print(x$logdens - xx$logdens)
#' 
#'  ## the use of a constraint
#'  constr = list(A = matrix(rep(1, n), 1, n), e = 0)
#'  x = inla.qsample(n, Q=Q, constr=constr)
#'  print(constr$A %*% x)
#' 
#'  ## control the RNG (require serial mode)
#'  x = inla.qsample(n, Q=Q, seed = 123, num.threads="1:1")
#'  ## restart from same seed,  only sample 1
#'  xx = inla.qsample(n=1, Q=Q, seed = 123, num.threads="1:1")
#'  ## continue from the save state, sample the remaining 2
#'  xxx = inla.qsample(n=n-1, Q=Q, seed = -1, num.threads="1:1")
#'  ## should be 0
#'  print(x - cbind(xx, xxx))
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
                           logdens = ifelse(missing(sample), FALSE, TRUE),
                           compute.mean = ifelse(missing(sample), FALSE, TRUE),
                           num.threads = if (seed == 0L) "0:0" else NULL,
                           selection = NULL,
                           verbose = inla.getOption("verbose"),
                           .debug = FALSE) 
{
    t.dir <- inla.tempdir()
    smtp <- match.arg(inla.getOption("smtp"), c("taucs", "band", "default", "pardiso"))
    stopifnot(!missing(Q))
    stopifnot(n >= 1L)

    if (is.null(num.threads)) {
        num.threads <- inla.getOption("num.threads")
    }
    if (seed != 0L) {
        num.threads.user <- inla.parse.num.threads(num.threads)
        num.threads <- inla.parse.num.threads("1:1")
        if (num.threads != num.threads.user) {
            warning("Since 'seed!=0', parallel model is disabled and serial model is selected")
        }
    } else {
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
        selection <- 1:nrow(Q)
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
    if (inla.os("linux") || inla.os("mac") || inla.os("mac.arm64")) {
        if (.debug) {
            print(paste(
                shQuote(inla.call.no.remote()), "-s -m qsample",
                paste0("-t", num.threads), "-r", reordering, "-z", seed, "-S", smtp,
                if (verbose) "-v" else "",
                Q.file, x.file, as.integer(n), rng.file,
                sample.file, b.file, mu.file, constr.file, cmean.file, selection.file
            ))
            browser()
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
        return(list(sample = samples, logdens = ld, mean = c(cmean)))
    } else {
        return(samples)
    }
}
