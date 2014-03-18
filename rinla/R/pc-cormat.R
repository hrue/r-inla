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
##! }
##!\value{%%
##!  \code{inla.pc.dcormat} gives the density and 
##!  \code{inla.pc.rcormat} generates random deviates (a
##!        list of correlation matrices).
##! }
##! \author{Havard Rue \email{hrue@math.ntnu.no}}
##! \examples{
##! }

inla.pc.cormat.dim2p = function(dim)
{
    p = round(1/2 + 1/2 * sqrt(1+4L*dim*2L))
    stopifnot(abs(dim - inla.pc.cormat.p2dim(p)) < sqrt(.Machine$double.eps))
    return (p)
}

inla.pc.cormat.p2dim = function(p)
{
    return (p*(p-1L)/2L)
}

inla.pc.cormat.theta2R = function(theta)
{
    ## convert from theta to a correlation matrix

    p = inla.pc.cormat.dim2p(length(theta))
    theta.m = matrix(NA, p, p)
    theta.m[lower.tri(theta.m)] = theta

    B = matrix(0, p, p)
    B[1, 1] = 1.0
    for(i in 2L:p) {
        for(j in 1L:i) {
            if (j == 1) {
                B[i, j] = cos(theta.m[i, j])
            } else if (j >= 2L && j <= i-1) {
                B[i, j] = cos(theta.m[i, j])*prod(sin(theta.m[i, 1:(j-1)]))
            } else if (j == i) {
                B[i, j] = prod(sin(theta.m[i, 1:(j-1)]))
            } else {
                stop("This should not happen.")
            }
        }
    }
    R = B %*% t(B)

    return(R)
}

inla.pc.cormat.R2theta = function(R)
{
    ## convert from correlation matrix to theta vector
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
    theta = L[lower.tri(L)]

    return (theta)
}

inla.pc.cormat.r2R = function(r)
{
    p = inla.pc.cormat.dim2p(length(r))
    R = matrix(1, p, p)
    R[lower.tri(R)] = r
    R = t(R)
    R[lower.tri(R)] = r
    return (R)
}

inla.pc.cormat.R2r = function(R)
{
    return (R[lower.tri(R)])
}

inla.pc.cormat.permute = function(R)
{
    r = inla.pc.cormat.R2r(R)
    r = r[order(runif(length(r)))]
    R = inla.pc.cormat.r2R(r)

    return (R)
    
}

inla.pc.cormat.rtheta = function(n=1, p, lambda = 1)
{
    stopifnot(!missing(p) && p > 1)
    stopifnot(lambda > 0)
    stopifnot(n >= 1)
    
    rcormat = function(n, p, lambda)
    {
        m = inla.pc.cormat.p2dim(p)

        ## sample a point on the simplex
        r = rexp(1, rate = lambda)
        gamma = rexp(m)
        gamma = (gamma/sum(gamma)) * r^2/2 
        theta = asin(exp(-gamma))
        ## sample branch by random
        branch = sample(c(0, 1), size = m, replace=TRUE)
        theta = branch * theta + (1-branch) * (pi - theta)
        return (theta)
    }
    
    x = unlist(lapply(rep(1, n), rcormat, p=p, lambda = lambda))
    x = matrix(x, n, inla.pc.cormat.p2dim(p), byrow=TRUE)
    return (x)
}

inla.pc.cormat.dtheta = function(theta, lambda = 1, log = TRUE)
{
    ## q is the dimension of theta. p is the dimension of R. this is
    ## not the notation used in the paper.
    q = length(theta)
    R = inla.pc.cormat.theta2R(theta)
    d = sqrt(-log(det(R)))
    ld = (log(lambda) 
          + lfactorial(q) 
          + (3*q/2 + 1) * log(2.0)
          - 0.5 * log(q+1) 
          - (2*q+1) * log(d)
          - lambda * d
          + sum(log(abs(1/tan(theta)))))
    
    return (if (log) ld else exp(ld))
}    

inla.pc.cormat.internal.test.1 = function(n, p, lim = 0.99)
{
    m = inla.pc.cormat.p2dim(p)
    ppi = lim*pi
    X = matrix(runif(n*m, min=-ppi, max=ppi), n, m)
    INLA:::inla.require("parallel")
    d = mclapply(1:n,
            function(k, lambda, log) {
                return (inla.pc.cormat.dtheta(X[k, ], lambda = lambda, log=log))
            }, lambda = 1, log = FALSE,
            mc.cores = detectCores())
    d = unlist(d)
    mu = mean(d)
    s = sd(d)/sqrt(n)
    z = mu/s

    return(list(mean = mu, sd = s, z = z))
}

    
    
