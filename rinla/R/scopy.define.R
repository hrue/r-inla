#' Get the parameterisation of the `scopy` model
#' 
#' This function provide access to the parameterisation of the `scopy` model
#' 
#' @aliases inla.scopy.define scopy.define
#' @param n Number of parameters, either 2 (intercept and slope)
#' or >= 5 (intercept, slope and a spline representing the deviation from it)
#' @return A list with the number of parameters and matrix defining basis vectors for the scopy model
#' @author Havard Rue \email{hrue@@r-inla.org}
#' @name scopy.define
#' @rdname scopy.define
#' @export

`inla.scopy.define` <- function(n = 5L)
{
    stopifnot(n == 2 || n >= 5)
    if (n == 2) {
        VV = cbind(1, c(-0.5, 0.5))
    } else {
        Q <- inla.rw2(n, scale.model = TRUE)
        e <- eigen(Q)
        idx.remove <- (n-1):n
        V <- e$vectors[, -idx.remove, drop = FALSE]
        lambda <- e$values[ -idx.remove]
        ## this is the sqrt(covariance.matrix) which makes the remaining easy. put them in reverse
        ## order, which make more sense
        VV <- matrix(0, n, n)
        VV[, 1:2] <- cbind(rep(1, n), seq(-0.5, 0.5, len = n))
        m <- n - 2
        for (i in seq_len(m)) {
            j <- m - i + 3
            VV[, j] <- V[, i] * sqrt(1 / lambda[i])
        }
    }
    return (list(n = n, W = VV))
}
