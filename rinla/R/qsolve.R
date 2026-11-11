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

`inla.qsolve` <- function(Q, B, method = c("solve", "forward", "backward"), num.threads = "2L")
{
    Q <- inla.as.sparse(Q)
    if (all(Q@i >= Q@j) || all(Q@i <= Q@j)) {
        Q <- Q + t(Q)
        diag(Q) <- diag(Q) / 2
    }

    if (is.null(num.threads)) num.threads <- inla.getOption("num.threads")
    cores <- as.integer(max(1, strsplit(inla.parse.num.threads(num.threads), ":")[[1]][1]))

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
