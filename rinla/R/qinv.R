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
#' @param num.threads Number of threads to be used
#' @return `inla.qinv` returns a `sparseMatrix` of type
#' `dgTMatrix` with the inverse of `Q`,  where `Q` itself is non-zero.
#' @author Havard Rue \email{hrue@@r-inla.org}
#' @examples
#' 
#'  ## dense matrix example
#'  n = 10
#'  A = matrix(runif(n^2), n, n)
#'  Q = A %*% t(A)
#'  print(mean(abs(as.matrix(inla.qinv(Q)) - solve(Q))))
#' 
#'  ## sparse matrix example
#'  rho = 0.9
#'  Q = toeplitz(c(1+rho^2, -rho,  rep(0, n-3), -rho)) / (1-rho^2)
#'  Q = inla.as.dgTMatrix(Q)
#'  Q.inv = inla.qinv(Q)
#'  marginal.variances = diag(inla.qinv(Q))
#' 
#' @name qinv
#' @rdname qinv
#' @export
`inla.qinv` <- function(Q, constr = NULL, num.threads = NULL)
{
    if (!requireNamespace("sTiles", quietly = TRUE)) {
        stop("package 'sTiles' is required")
    }
    stopifnot(!missing(Q))
    if (is.null(num.threads)) {
        num.threads <- inla.getOption("num.threads")
    }
    num.threads <- strsplit(inla.parse.num.threads(num.threads), ":")[[1]][1]

    s <- sTiles::sTiles(Q, inverse = TRUE, cores = num.threads)
    Q.inv <- sTiles::sTiles_selinv_matrix(s)
    Q.inv <- inla.as.sparse(Q.inv)
    Q.graph <- Q
    Q.graph[Q.graph != 0.0] <- 1
    Q.inv <- inla.as.sparse(Q.inv * Q.graph)

    if (!is.null(constr)) {
        stopifnot(is.matrix(constr$A))
        Qinv_At <- as.matrix(sTiles::sTiles_solve(s, t(constr$A)))
        LL <- t(chol(solve(constr$A %*% Qinv_At)))
        BB <- Qinv_At %*% LL
        ii <- Q.inv@i + 1L
        jj <- Q.inv@j + 1L
        d  <- numeric(length(ii))
        for (l in seq_len(ncol(BB))) {
            d <- d + BB[ii, l] * BB[jj, l]
        }
        Q.inv@x <- Q.inv@x - d
    }

    Q.graph <- NULL
    sTiles::sTiles_close(s)

    return(Q.inv)
}
