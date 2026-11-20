#' Solves linear SPD systems
#' 
#' This function is a wrapper to 'sTiles_solve' with a slightly different interface
#' 
#' @aliases inla.qsolve qsolve
#' @param Q A SPD (sparse) matrix
#' @param B The right hand side matrix
#' @param method The system to solve, one of `"solve"`, `"forward"`
#' or `"backward"`. Let `Q = L L^T`, where `L` is lower
#' triangular (the Cholesky triangle), then `method="solve"` solves
#' `L L^T X = B` or equivalently `Q X = B`, `method="forward"`
#' solves `L X = B`, and `method="backward"` solves `L^T X = B`.
#' @param num.threads Number of threads to use. If in format "A:B:C", then "A" is used.
#'
#' @return `inla.qsolve` returns a matrix `X`, which is the solution
#' of `Q X = B`, `L X = B` or `L^T X = B` depending on the value
#' of `method`.
#'
#' @note Note that reordering of `Q`, and then `L`, complicates matters for forward and backward
#'     solve: see [sTiles_solve()] for details.
#' 
#' @examples
#' n <- 10
#' A <- matrix(rnorm(n^2) * rbinom(n^2, size = 1, prob = 1/n), n, n)
#' diag(A) <- rnorm(n)
#' Q <- inla.as.sparse(A %*% t(A))
#' m <- 3
#' B <- matrix(rnorm(n*m), n, m)
#' X.solve <- inla.qsolve(Q, B, method = "solve")
#' X.forward <- inla.qsolve(Q, B, method = "forward")
#' X.backward <- inla.qsolve(Q, B, method = "backward")
#' print(max(abs(Q %*% X.solve - B)))
#'
#' @author Havard Rue \email{hrue@@r-inla.org}
#' @name qsolve
#' @rdname qsolve
#' @export

`inla.qsolve` <- function(Q, B, reordering = inla.reorderings(),
                               method = c("solve", "forward", "backward"))
{
    t.dir <- inla.tempdir()
    smtp <- match.arg(inla.getOption("smtp"), c("stiles", "taucs", "band", "default"))
    if (smtp == "stiles") smtp <- "taucs"
    Q <- inla.sparse.check(Q)
    if (is(Q, "dgTMatrix")) {
        Qfile <- inla.write.fmesher.file(Q, filename = inla.tempfile(tmpdir = t.dir))
    } else if (is.character(Q)) {
        Qfile <- Q
    } else {
        stop("This should not happen.")
    }

    ## ensure B is a dense matrix if its a sparse
    if (is(B, "Matrix")) {
        B <- as.matrix(B)
    }
    if (is.matrix(B)) {
        Bfile <- inla.write.fmesher.file(B, filename = inla.tempfile(tmpdir = t.dir))
    } else if (is.character(B)) {
        Bfile <- B
    } else {
        stop("This should not happen.")
    }

    if (is.list(reordering)) {
        ## argument is the output from inla.qreordering()
        reordering <- reordering$name
    }
    reordering <- match.arg(reordering)
    method <- match.arg(method)

    inla.set.environment()
    env <- inla.run.environment.set()
    Xfile <- inla.tempfile(tmpdir = t.dir)
    if (inla.os("linux") || inla.os("mac") || inla.os("mac.arm64")) {
        s <- system(paste(
            shQuote(inla.call.no.remote()), "-s -m qsolve", "-r",
            reordering, "-S", smtp, Qfile, Xfile, Bfile, method
        ), intern = TRUE)
    } else if (inla.os("windows")) {
        s <- system(paste(
            shQuote(inla.call.no.remote()), "-s -m qsolve", "-r",
            reordering, "-S", smtp, Qfile, Xfile, Bfile, method
        ), intern = TRUE)
    } else {
        stop("\n\tNot supported architecture.")
    }

    inla.run.environment.unset(env)
    X <- inla.read.fmesher.file(Xfile)
    unlink(t.dir, recursive = TRUE)

    return(X)
}

`inla.qsolve.stiles` <- function(Q, B, method = c("solve", "forward", "backward"),
                                 num.threads = NULL)
{
    Q <- inla.as.sparse(Q)
    if (all(Q@i >= Q@j) || all(Q@i <= Q@j)) {
        Q <- Q + t(Q)
        diag(Q) <- diag(Q) / 2
    }

    if (is.null(num.threads)) num.threads <- inla.getOption("num.threads")
    cores <- numeric(2)
    cores[1] <- max(1, as.integer(strsplit(inla.parse.num.threads(num.threads), ":")[[1]][1]))
    cores[2] <- max(1, as.integer(strsplit(inla.parse.num.threads(num.threads), ":")[[1]][2]))
    cores <- prod(cores)

    method <- match.arg(method)
    if (method == "solve") {
        system <- "A"
    } else if (method == "forward") {
        system <- "L"
    } else if (method == "backward") {
        system <- "Lt"
    } else {
        stop("This should not happen")
    }

    Q <- inla.as.sparse(Q)
    s <- sTiles::sTiles(Q, cores = cores)
    X <- sTiles::sTiles_solve(s, B, system = system)
    sTiles::sTiles_close(s)

    return (X)
}
