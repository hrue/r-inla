## Export: inla.pc.rcormat

##! \name{pc.rcormat}
##! \alias{inla.pc.rcormat}
##! \alias{pc.rcormat}
##! \alias{pc.dcormat}
##! \alias{inla.pc.dcormat}
##! 
##! \title{Utility functions for the PC prior for a correlation matrix}
##! 
##! \description{Functions to evaluate and sample from the 
##!              PC prior for a correlation matrix.}
##! \usage{
##! inla.pc.rcormat(n=1, p, lambda=1, permute = TRUE)
##! inla.pc.dcormat(cormat, lambda=1, permute = TRUE, log=FALSE)
##! }
##! \arguments{
##!   \item{n}{Number of observations}
##!   \item{p}{The dimension of the matrix}
##!   \item{lambda}{The rate parameter}
##!   \item{permute}{Make the prior exchangable.}
##!   \item{cormat}{A correlation matrix}
##!   \item{log}{Logical. Return the density in natural or log-scale.}
##! }
##! \details{
##!     The permute-option can have a high cost as the number of terms in the
##!     sum is \code{factorial(p)}.
##! }
##!\value{%%
##!  \code{inla.pc.dcormat} gives the density and 
##!  \code{inla.pc.rcormat} generates random deviates (a
##!        list of correlation matrices).
##! }
##! \author{Havard Rue \email{hrue@math.ntnu.no}}
##! \examples{
##!  A = inla.pc.rcormat(n=1, p=3, lambda = 1, permute=FALSE)[[1]]
##!  print(inla.pc.dcormat(A))
##!  A = inla.pc.rcormat(n=10, p=3, lambda = 1, permute=TRUE)
##!  lapply(A, inla.pc.dcormat)
##! }

inla.pc.rcormat = function(n=1, p, lambda=1, permute = TRUE)
{
    rcormat = function(n, p, lambda, permute)
    {
        stopifnot(all(n == 1))
        stopifnot(!missing(p) && p > 1)
        stopifnot(!missing(lambda) && lambda > 0)
    
        ## number of parameters
        m = p*(p-1)/2L

        ## sample a point on the simplex
        r = rexp(1, rate = lambda)
        gamma = rexp(m)
        gamma = (gamma/sum(gamma)) * r^2/2 
        theta.vec = asin(exp(-gamma))
        ## sample branch by random
        branch = sample(c(0, 1), size = m, replace=TRUE)
        theta.vec = branch * theta.vec + (1-branch) * (pi - theta.vec)

        ## put theta.vec into a matrix, as it makes the code below easier.
        theta = matrix(NA, p, p)
        theta[lower.tri(theta)] = theta.vec
        
        B = matrix(0, p, p)
        B[1, 1] = 1
        for(i in 2L:p) {
            for(j in 1L:i) {
                if (j == 1) {
                    B[i, j] = cos(theta[i, j])
                } else if (j >= 2L && j <= i-1) {
                    B[i, j] = cos(theta[i, j])*prod(sin(theta[i, 1:(j-1)]))
                } else if (j == i) {
                    B[i, j] = prod(sin(theta[i, 1:(j-1)]))
                } else {
                    stop("This should not happen.")
                }
            }
        }

        mat = B %*% t(B)
        if (permute) {
            r = mat[upper.tri(mat)]
            r = r[order(runif(m))]
            mat[upper.tri(mat)] = r
            ## to make the same order....
            mat = t(mat)
            mat[upper.tri(mat)] = r
        } 
        ## so its truly 1,  not just numerically 1
        diag(mat) = 1.0
        return (mat)
    }

    stopifnot(n >= 1)
    return (lapply(rep(1, n), rcormat,
                   p=p, lambda = lambda, permute = permute))
}

inla.pc.dcormat = function(cormat, lambda=1, permute = TRUE, log=FALSE)
{
    r2theta = function(R)
    {
        ## from a correlation matrix, compute the corresponding
        ## theta-vector.
        L = t(chol(R))
        p = dim(R)[1]
        ## overwrite L with theta
        for(i in 2:p) {
            if (i == 2) {
                L[2, 1] = acos(L[2, 1])
            } else {
                for(j in 1:(i-1)) {
                    if (j == 1) {
                        L[i, j] = acos(L[i, j])
                    } else {
                        L[i, j] = acos(L[i, j]/prod(sin(L[i, 1:(j-1)])))
                    }
                }
            }
        }
        return (L[lower.tri(L)])
    }

    log.dens = function(R)
    {
        p = dim(R)[1]
        theta = r2theta(R)
        d = sqrt(-log(det(R)))
        log.dens = (log(lambda) 
                    + lfactorial(p) 
                    + (3*p/2 + 1) * log(2.0)
                    - 0.5 * log(p+1) 
                    - (2*p+1) * log(d)
                    - lambda * d
                    + sum(log(abs(1/tan(theta)))))
        return (log.dens)
    }

    p = dim(cormat)[1]
    R = cormat
    stopifnot(all(diag(R) == 1.0))
    
    if (!permute) {
        ld = log.dens(R)
    } else {
        stopifnot(INLA:::inla.require("combinat"))
        r = R[lower.tri(R)]
        m = length(r)
        perm = permn(m)
        if (FALSE) {
            ## old code
            ld = numeric(m)
            for(k in 1:length(perm)) {
                rr = r[perm[[k]]]
                R[lower.tri(R)] = rr
                R = t(R)
                R[lower.tri(R)] = rr
                ld[k] = log.dens(R)
            }
        } else {
            ld = lapply(1:length(perm),
                    function(k, RR) {
                        rr = r[perm[[k]]]
                        RR[lower.tri(RR)] = rr
                        RR = t(RR)
                        RR[lower.tri(RR)] = rr
                        return(log.dens(RR))
                    }, RR=R)
            ld = unlist(ld)
        }
        ld.max = max(ld)
        ld = ld.max + log(mean(exp(ld - ld.max)))
    }
    
    return (if (log) ld else exp(ld))
}    
