#' Get the parameterisation of the `scopy` model
#' 
#' This function provide access to the parameterisation of the `scopy` model
#' 
#' @aliases inla.scopy.define scopy.define
#' @param n Number of parameters, either 2 (intercept and slope)
#' or >= 5 (intercept, slope and a spline representing the deviation from it)
#' @param w A vector a positive weights defining the weighted sum-to-zero constraint.
#' By default, all weights are set to 1. Only in effect if n >= 5.
#' @return A list with the number of parameters and matrix defining basis vectors for the scopy model
#' @author Havard Rue \email{hrue@@r-inla.org}
#' @name scopy.define
#' @rdname scopy.define
#' @export

`inla.scopy.define` <- function(n = 5L, w = NULL)
{
    stopifnot(n == 2 || n >= 5)
    if (n == 2) {
        VV = cbind(1, c(-0.5, 0.5))
    } else {
        if (is.null(w)) w <- rep(1, n)
        stopifnot(all(w > 0))
        w <- w / sum(w)
        w <- matrix(w, 1, n)
        Q <- inla.rw2(n, scale.model = TRUE)
        QQ <- Q + 1E-8 * diag(n)
        SS <- solve(QQ)
        SS <- SS - SS %*% t(w) %*% solve(w %*% SS %*% t(w)) %*% w %*% SS
        SS <- 0.5 * (SS + t(SS))
        e <- eigen(SS)
        ## linear is the first one, but the sign might be reversed
        if (all(diff(e$vectors[, 1]) < 0)) e$vectors[, 1] <- -e$vectors[, 1]
        stopifnot(all(diff(e$vectors[, 1]) > 0))
        ## intercept is last one
        VV <- matrix(0, n, n)
        VV[, 1:2] <- cbind(1, e$vectors[, 1] / (e$vectors[n, 1] - e$vectors[1, 1]))
        for (i in 2:(n-1)) {
            VV[, i + 1] <- e$vectors[, i] * sqrt(e$values[i])
        }
    }
    return (list(n = n, W = VV))
}
