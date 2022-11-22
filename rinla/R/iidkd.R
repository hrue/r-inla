## Export: inla.iidkd.sample

## !\name{inla.iidkd.sample}
## !\alias{inla.iidkd.sample}
## !\alias{iidkd.sample}
## !
## !\title{Provide samples from the \code{iidkd} component (experimental)}
## !\description{
## ! This function provide samples of the \code{iidkd} component 
## ! using more interpretable parameters
## !}
## !\usage{
## !inla.iidkd.sample(n = 10^4, result, name, return.cov = FALSE)
## !}
## !
## !\arguments{
## !\item{n}{Integer Number of samples to use}
## !\item{result}{inla-object An object of class \code{inla}, ie a result of a call to \code{inla()}}
## !\item{name}{Character The name of the \code{iidkd} component}
## !\item{return.cov}{Logical Return samples of the covariance matrix instead of
## !      stdev/correlation matrix described below?}
## !}
## !\value{A list of sampled matrices, with (default) correlations on the off-diagonal and
## !standard-deviations on the diagonal}
## !\seealso{\code{inla.doc("iidkd")}}
## !\author{Havard Rue \email{hrue@r-inla.org}}

`inla.iidkd.sample` <- function(n = 10^4, result, name, return.cov = FALSE) 
{
    stopifnot(!missing(result))
    stopifnot(inherits(result, "inla"))
    stopifnot(!missing(name))
    stopifnot(n > 0)

    k.max <- 10
    theta.max <- (k.max * (k.max + 1L)) / 2L

    names <- paste0("Theta", 1:theta.max, " for ", name)
    idx <- which(rownames(result$summary.hyperpar) %in% names)
    if (length(idx) == 0) {
        return (list())
    }
    len <- length(idx)
    k <- 0
    for(kk in 2:k.max) {
        if (len == kk*(kk+1)/2) {
            k <- kk
            break
        }
    }
    if (k == 0) {
        stop(paste0("Cannot find correct dimension for: ", name))
    }

    theta2matrix <- function(i, thetas, dim, dim2, ret.cov) {
        theta <- thetas[i, ]
        L <- matrix(0, dim, dim)
        diag(L) <- exp(theta[1:dim])
        L[lower.tri(L)] <- theta[dim + 1:dim2]
        ## rescue Inf values
        L[is.infinite(L)] <- 1/.Machine$double.eps
        S <- try(solve(L %*% t(L)), silent = TRUE)

        ## if singular, then invert what is invertable
        if (inherits(S, "try-error")) {
            S <- inla.ginv(L %*% t(L))
            ## do an extra check
            diag(S) <- pmax(0, diag(S))
            e <- eigen(S)
            S <- e$vectors %*% diag(pmax(0, e$values)) %*% t(e$vectors)
        } 
        S <- (S + t(S))/2.0 ## force matrix to be symmetric

        if (ret.cov) {
            return (S)
        } else {
            sd <- sqrt(diag(S))
            iSigma <- diag(1/sd)
            iSigma[is.infinite(iSigma)] <- 0
            Cor <- iSigma %*% S %*% iSigma
            diag(Cor) <- sd
            return(Cor)
        }
    }
    xx <- inla.hyperpar.sample(n, result, intern = TRUE)[, idx, drop = FALSE]
    res <- lapply(1:n, theta2matrix, dim = k, dim2 = k*(k-1)/2, thetas = xx,
                  ret.cov = return.cov)

    return (res)
}
