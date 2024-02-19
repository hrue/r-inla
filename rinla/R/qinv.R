#' Computes (parts of) the inverse of a SPD sparse matrix
#' 
#' This routine use the GMRFLib implementation which compute parts of the
#' inverse of a SPD sparse matrix.  The diagonal and values for the neighbours
#' in the inverse, are provided.
#' 
#' 
#' @aliases inla.qinv qinv
#' @param Q A SPD matrix, either as a (dense) matrix or sparseMatrix.
#' @param constr Optional linear constraints; see `?INLA::f` and argument
#' `extraconstr`
#' @param reordering The type of reordering algorithm to be used for
#' `TAUCS`; either one of the names listed in `inla.reorderings()` or
#' the output from `inla.qreordering(Q)`.  The default is "auto" which try
#' several reordering algorithm and use the best one for this particular
#' matrix.
#' @param num.threads Maximum number of threads the `inla`-program will
#' use, or as 'A:B' defining the number threads in the outer (A) and inner (B)
#' layer for nested parallelism.
#' @return `inla.qinv` returns a `sparseMatrix` of type
#' `dgTMatrix` with the diagonal and values for the neigbours in the
#' inverse. Note that the full inverse is NOT provided!
#' @author Havard Rue \email{hrue@@r-inla.org}
#' @examples
#' 
#' 
#'  ## dense matrix example
#'  n = 10
#'  A = matrix(runif(n^2), n, n)
#'  Q = A %*% t(A)
#'  print(mean(abs(inla.qinv(Q) - solve(Q))))
#' 
#'  ## sparse matrix example
#'  rho = 0.9
#'  Q = toeplitz(c(1+rho^2, -rho,  rep(0, n-3), -rho)) / (1-rho^2)
#'  Q = inla.as.dgTMatrix(Q)
#'  Q.inv = inla.qinv(Q)
#' 
#'  ## compute the marginal variances as a vector from a precision matrix
#'  marginal.variances = diag(inla.qinv(Q))
#' 
#'  ## read the sparse matrix from a file in the 'i, j, value' format
#'  filename = tempfile()
#'  write(t(cbind(Q@i+1L,  Q@j+1L,  Q@x)), ncol=3, file=filename)
#'  Qinv = inla.qinv(filename)
#'  unlink(filename)
#'  
#' @name qinv
#' @rdname qinv
#' @export
`inla.qinv` <- function(Q, constr, reordering = INLA::inla.reorderings(),
                        num.threads = NULL) {
    t.dir <- inla.tempdir()
    smtp <- match.arg(inla.getOption("smtp"), c("taucs", "band", "default", "pardiso"))
    if (is.null(num.threads)) {
        num.threads <- inla.getOption("num.threads")
    }
    num.threads <- inla.parse.num.threads(num.threads)

    Q <- inla.sparse.check(Q)
    if (is(Q, "dgTMatrix")) {
        qinv.file <- inla.write.fmesher.file(Q, filename = inla.tempfile(tmpdir = t.dir))
    } else if (is.character(Q)) {
        qinv.file <- Q
    } else {
        stop("This chould not happen.")
    }

    constr.file <- inla.tempfile(tmpdir = t.dir)
    if (!missing(constr) && !is.null(constr)) {
        stopifnot(is.list(constr))
        A <- as.matrix(constr$A)
        e <- as.numeric(constr$e)
        stopifnot(ncol(A) == ncol(Q))
        stopifnot(nrow(A) == length(e))
        xx <- matrix(c(nrow(A), c(A), c(e)), ncol = 1)
        inla.write.fmesher.file(xx, filename = constr.file)
    }

    if (is.list(reordering)) {
        ## argument is the output from inla.qreordering()
        reordering <- reordering$name
    }
    reordering <- match.arg(reordering)

    inla.set.environment()
    out.file <- inla.tempfile(tmpdir = t.dir)
    if (inla.os("linux") || inla.os("mac") || inla.os("mac.arm64")) {
        s <- system(paste(
            shQuote(inla.call.no.remote()), "-s -m qinv",
            "-r", reordering, "-S", smtp, paste0("-t", num.threads),
            qinv.file, constr.file, out.file
        ), intern = TRUE)
    } else if (inla.os("windows")) {
        s <- system(paste(
            shQuote(inla.call.no.remote()), "-s -m qinv",
            "-r", reordering, "-S", smtp, paste0("-t", num.threads),
            qinv.file, constr.file, out.file
        ), intern = TRUE)
    } else {
        stop("\n\tNot supported architecture.")
    }

    Qinv <- inla.read.fmesher.file(out.file)
    unlink(t.dir, recursive = TRUE)

    return(Qinv)
}
