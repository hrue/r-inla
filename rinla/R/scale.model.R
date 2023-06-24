#' Scale an intrinsic GMRF model
#' 
#' This function scales an intrinsic GMRF model so the geometric mean of the
#' marginal variances is one
#' 
#' 
#' @aliases inla.scale.model scale.model
#' @param Q A SPD matrix, either as a (dense) matrix or `sparseMatrix`
#' @param constr Linear constraints spanning the null-space of `Q`; see
#' `?INLA::f` and argument `extraconstr`
#' @param eps A small constant added to the diagonal of `Q` if
#' `constr`
#' @return `inla.scale.model` returns a `sparseMatrix` of type
#' `dgTMatrix` scaled so the geometric mean of the marginal variances (of
#' the possible non-singular part of `Q`) is one, for each connected
#' component of the matrix.
#' @author Havard Rue \email{hrue@@r-inla.org}
#' @examples
#' 
#'  ## Q is singular
#'  data(Germany)
#'  g = system.file("demodata/germany.graph", package="INLA")
#'  Q = -inla.graph2matrix(g)
#'  diag(Q) = 0
#'  diag(Q) = -rowSums(Q)
#'  n = dim(Q)[1]
#'  Q.scaled = inla.scale.model(Q, constr = list(A = matrix(1, 1, n), e=0))
#'  print(diag(MASS::ginv(as.matrix(Q.scaled))))
#' 
#'  ## Q is singular with 3 connected components
#'  g = inla.read.graph("6 1 2 2 3 2 2 1 3 3 2 1 2 4 1 5 5 1 4 6 0")
#'  print(paste("Number of connected components", g$cc$n))
#'  Q = -inla.graph2matrix(g)
#'  diag(Q) = 0
#'  diag(Q) = -rowSums(Q)
#'  n = dim(Q)[1]
#'  Q.scaled = inla.scale.model(Q, constr = list(A = matrix(1, 1, n), e=0))
#'  print(diag(MASS::ginv(as.matrix(Q.scaled))))
#' 
#'  ## Q is non-singular with 3 connected components. no constraints needed
#'  diag(Q) = diag(Q) + 1
#'  Q.scaled = inla.scale.model(Q)
#'  print(diag(MASS::ginv(as.matrix(Q.scaled))))
#'  
#' @name scale.model
#' @rdname scale.model
#' @export
inla.scale.model.internal <- function(Q, constr = NULL, eps = sqrt(.Machine$double.eps)) {
    ## return also the scaled marginal variances

    marg.var <- rep(0, nrow(Q))
    Q <- inla.as.sparse(Q)
    g <- inla.read.graph(Q)
    for (k in seq_len(g$cc$n)) {
        i <- g$cc$nodes[[k]]
        n <- length(i)
        QQ <- Q[i, i, drop = FALSE]
        if (n == 1) {
            QQ[1, 1] <- 1
            marg.var[i] <- 1
        } else {
            cconstr <- constr
            if (!is.null(constr)) {
                ## the GMRFLib will automatically drop duplicated constraints; how convenient...
                cconstr$A <- constr$A[, i, drop = FALSE]
                eeps <- eps
            } else {
                eeps <- 0
            }
            idx.zero <- which(rowSums(abs(cconstr$A)) == 0)
            if (length(idx.zero) > 0) {
                cconstr$A <- cconstr$A[-idx.zero,, drop = FALSE]
                cconstr$e <- cconstr$e[-idx.zero]
            }
            res <- inla.qinv(QQ + Diagonal(n) * max(diag(QQ)) * eeps, constr = cconstr)
            fac <- exp(mean(log(diag(res))))
            QQ <- fac * QQ
            marg.var[i] <- diag(res) / fac
        }
        Q[i, i] <- QQ
    }
    return(list(Q = Q, var = marg.var))
}

#' @rdname scale.model
#' @export
inla.scale.model <- function(Q, constr = NULL, eps = sqrt(.Machine$double.eps)) {
    res <- inla.scale.model.internal(Q = Q, constr = constr, eps = eps)
    return(res$Q)
}
