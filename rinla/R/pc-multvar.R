## Export: inla.pc.multvar.simplex.d
## Export: inla.pc.multvar.simplex.r
## Export: inla.pc.multvar.sphere.d
## Export: inla.pc.multvar.sphere.r
## Export: inla.pc.multvar.simplex.general.d
## Export: inla.pc.multvar.simplex.general.r
## Export: inla.pc.multvar.sphere.general.d
## Export: inla.pc.multvar.sphere.general.r

##! \name{pc.multvar}
##! \alias{inla.pc.multvar}
##! \alias{inla.pc.multvar.simplex}
##! \alias{inla.pc.multvar.simplex.d}
##! \alias{inla.pc.multvar.simplex.r}
##! \alias{inla.pc.multvar.simplex.general}
##! \alias{inla.pc.multvar.simplex.general.d}
##! \alias{inla.pc.multvar.simplex.general.r}
##! \alias{inla.pc.multvar.sphere}
##! \alias{inla.pc.multvar.sphere.d}
##! \alias{inla.pc.multvar.sphere.r}
##! \alias{inla.pc.multvar.sphere.general}
##! \alias{inla.pc.multvar.sphere.general.d}
##! \alias{inla.pc.multvar.sphere.general.r}
##!
##! \title{Multivariate PC priors}
##! 
##! \description{Functions to evaluate and simulate from 
##!              multivariate PC priors: The simplex and sphere case}
##! \usage{
##! inla.pc.multvar.h.default(x, inverse = FALSE, derivative = FALSE)
##! inla.pc.multvar.simplex.r(n, p, lambda = 1, h = inla.pc.multvar.h.default)
##! inla.pc.multvar.simplex.d(x, lambda = 1, log = FALSE, h = inla.pc.multvar.h.default)
##! inla.pc.multvar.simplex.general.r(n = NULL, lambda = 1, h = inla.pc.multvar.h.default, b = NULL)
##! inla.pc.multvar.simplex.general.d(x = NULL, lambda = 1, log = FALSE, h = inla.pc.multvar.h.default, b = NULL)
##! inla.pc.multvar.sphere.r(n, p, lambda = 1, h = inla.pc.multvar.h.default)
##! inla.pc.multvar.sphere.d(x, lambda = 1, log = FALSE, h = inla.pc.multvar.h.default)
##! inla.pc.multvar.sphere.general.r(n = NULL, lambda = 1, h = inla.pc.multvar.h.default, H = NULL)
##! inla.pc.multvar.sphere.general.d(x = NULL, lambda = 1, log = FALSE, h = inla.pc.multvar.h.default, H = NULL)
##! }
##! \arguments{
##!   \item{x}{Samples to evaluate. If input is a matrix then each row is a sample. If input is
##!            a vector then this is the sample.}
##!   \item{inverse}{Compute the inverse of the h()-function.}
##!   \item{derivative}{Compute the derivative of the h()-function. (derivative of the inverse
##!                     function is not used).}
##!   \item{n}{Number of samples to generate.}
##!   \item{p}{The dimension of the PC-prior.}
##!   \item{lambda}{The lambda-parameter in the PC-prior.}
##!   \item{log}{Evaluate the density in log-scale or ordinary scale.}
##!   \item{h}{The h()-function,  defaults to \code{inla.pc.multvar.h.default}. See that code
##!            for an example of how to write a user-spesific function.}
##!   \item{b}{The b-vector (gradient) in the general expression for the simplex option,  \code{d(xi) = h(b^T xi)}}
##!   \item{b}{The H(essian)-matrix in the general expression for the sphere option,  \code{d(xi) =
##!            h(1/2 *xi^T H xi)}}
##!}
##! \details{
##!   These functions implements general multivariate PC-priors of the simplex and sphere type.
##! }
##! \value{%%
##!     \code{inla.pc.multvar.simplex.r} generate samples from the standard simplex case,  and
##!     \code{inla.pc.multvar.simplex.d} evaluate the density.
##!     \code{inla.pc.multvar.simplex.general.r} generate samples from the general simplex case,  and
##!     \code{inla.pc.multvar.simplex.general.d} evaluate the density.
##!     \code{inla.pc.multvar.sphere.r} generate samples from the standard sphere case,  and
##!     \code{inla.pc.multvar.sphere.d} evaluate the density.
##!     \code{inla.pc.multvar.sphere.general.r} generate samples from the general sphere case,  and
##!     \code{inla.pc.multvar.sphere.general.d} evaluate the density.
##! }
##! \author{Havard Rue \email{hrue@math.ntnu.no}}
##! \examples{
##! }

inla.pc.multvar.h.default = function(x, inverse = FALSE, derivative = FALSE)
{
    ## default: h(x) = sqrt(2*x)
    if (derivative) {
        if (inverse) {
            stopifnot("Options inverse=TRUE and derivative=TRUE is not allowed.")
        } else {
            return (1/sqrt(2*x))
        }
    } else {
        if (inverse) {
            return (1/2 * x^2) 
        } else {
            return (sqrt(2*x))
        }
    }
}

inla.pc.multvar.simplex.r = function(n, p, lambda = 1, h = inla.pc.multvar.h.default)
{
    ## simulate from the pc prior where d = h(\sum x_i), x_i >=0
    stopifnot(n > 0)
    stopifnot(p > 0)
    stopifnot(lambda > 0)

    x = matrix(NA, n, p)
    for(i in 1:n) {
        ## simulate z on the unit simplex
        z = rexp(p)
        z = z/sum(z)
        ## sample the size of the simplex
        d = rexp(1, rate = lambda)
        r = h(d, inverse=TRUE)
        ## and scale it
        x[i, ] = z*r
    }
    colnames(x) = paste("x", INLA:::inla.num(1:p), sep="")
    rownames(x) = paste("sample", INLA:::inla.num(1:n), sep="")

    return (x)
}

inla.pc.multvar.simplex.d = function(x, lambda = 1, log = FALSE, 
    h = inla.pc.multvar.h.default)
{
    ## evaluate the density from the pc prior where d = h(\sum x_i), x_i >=0.

    simplex.log.volume = function(n)
    {
        ## volumne of the n-simplex with n+1 vertices, see
        ## https://en.wikipedia.org/wiki/Simplex#Volume
        return (-lfactorial(n))
    }

    ## each row is a sample
    if (!is.matrix(x)) {
        x = matrix(c(x), ncol = length(x))
    }
    p = ncol(x)
    stopifnot(p > 0)
    stopifnot(lambda > 0)

    ## compute the density
    ldens = apply(x, 1,
        function(z) {
            r = sum(z)
            d = h(r)
            ld = (log(lambda) - lambda * d + log(abs(h(r, derivative = TRUE)))
                  - (p-1)*log(r) - simplex.log.volume(p-1))
            return (ld)
        })

    return (if (log) ldens else exp(ldens))
}

inla.pc.multvar.simplex.general.r = function(
    n = NULL, lambda = 1, h = inla.pc.multvar.h.default, b = NULL)
{
    return (inla.pc.multvar.simplex.general.core(
        n = n, lambda = lambda, h = h, b = b, mode = "r"))
}

inla.pc.multvar.simplex.general.d = function(
    x = NULL, lambda = 1, log = FALSE, h = inla.pc.multvar.h.default, b = NULL)
{
    return (inla.pc.multvar.simplex.general.core(
        x = x, lambda = lambda, log = log, h = h, b = b, mode = "d"))
}

inla.pc.multvar.simplex.general.core = function(
    x = NULL, n = NULL, lambda = 1, log = FALSE,
    h = inla.pc.multvar.h.default, b = NULL, mode = c("r", "d"))
{
    ## this is the case where d = h( b^T x)
    ## for vector b>0

    mode = match.arg(mode)
    ## need either b or H or both
    stopifnot(!is.null(b))

    ## ensure that b is a matrix
    if (!is.null(b) && !is.matrix(b)) {
        b = matrix(c(b), ncol = 1)
    }
    p = nrow(b)

    ## ensure x, if given, is a matrix. this is relevant if x is just one sample which can be
    ## given as a vector
    if (!is.null(x) && !is.matrix(x)) {
        x = matrix(c(x), ncol = length(x))
    }

    if (mode == "r") {
        X = inla.pc.multvar.simplex.r(n, p = p, lambda = lambda,  h = h)
        for(i in 1:n) {
            X[i, ] = t(1/b) %*% X[i,, drop = FALSE]
        }
        return (X)
    } else if (mode == "d") {
        n = nrow(x)
        ldens.x = numeric(n)
        for(i in 1:n) {
            theta = c(b) * c(x[i,])
            ldens = inla.pc.multvar.simplex.d(theta, lambda = lambda, log = TRUE, h = h)
            ldens.x[i] = ldens + sum(log(abs(b)))
        }
        return (if (log) ldens.x else exp(ldens.x))
    }

    stop("Should not happen")
    return()
}


inla.pc.multvar.sphere.r = function(n, p, lambda = 1, h = inla.pc.multvar.h.default)
{
    ## simulate from the pc prior where d = h(1/2 * \sum x_i^2)
    stopifnot(n > 0)
    stopifnot(p > 0)
    stopifnot(lambda > 0)

    x = matrix(NA, n, p)
    for(i in 1:n) {
        ## simulate z on the unit sphere
        z = rnorm(p)
        z = z/sqrt(sum(z^2))
        ## sample the radii of the sphere
        d = rexp(1, rate = lambda)
        r = sqrt(2*h(d, inverse=TRUE))
        ## and scale it
        x[i, ] = z*r
    }
    colnames(x) = paste("x", INLA:::inla.num(1:p), sep="")
    rownames(x) = paste("sample", INLA:::inla.num(1:n), sep="")

    return (x)
}

inla.pc.multvar.sphere.d = function(x, lambda = 1, log = FALSE, 
    h = inla.pc.multvar.h.default)
{
    ## evaluate the density from the pc prior where d = h(1/2 * \sum x_i^2), x_i >=0.

    sphere.log.volume = function(n)
    {
        ## volumne of the sphere in n-dimensions with radii=1
        ## https://en.wikipedia.org/wiki/Volume_of_an_n-ball
        return (n/2*log(pi) - lgamma(n/2 +1))
    }

    ## each row is a sample
    if (!is.matrix(x)) {
        x = matrix(c(x), ncol = length(x))
    }
    p = ncol(x)
    stopifnot(p > 0)
    stopifnot(lambda > 0)

    ldens = apply(x, 1,
        function(z) {
            r = sqrt(sum(z^2))
            d = h(0.5*r^2)
            ldens = (- log(2) + log(lambda) - lambda * d + log(abs(h(0.5*r^2, derivative = TRUE)))
                     - (p-1)*log(r) - sphere.log.volume(p-1) + log(mean(abs(z))))
            return (ldens)
        })

    return (if (log) ldens else exp(ldens))
}

inla.pc.multvar.sphere.general.r = function(
    n = NULL, lambda = 1, h = inla.pc.multvar.h.default, H = NULL)
{
    return (inla.pc.multvar.sphere.general.core(
        n = n, lambda = lambda, h = h, H = H, mode = "r"))
}

inla.pc.multvar.sphere.general.d = function(
    x = NULL, lambda = 1, log = FALSE, h = inla.pc.multvar.h.default, H = NULL)
{
    return (inla.pc.multvar.sphere.general.core(
        x = x, lambda = lambda, log = log, h = h, H = H, mode = "d"))
}


inla.pc.multvar.sphere.general.core = function(
    x = NULL, n = NULL, lambda = 1, log = FALSE,
    h = inla.pc.multvar.h.default, H = NULL, mode = c("r", "d"))
{
    ## this is the case where d = h(1/2 x^{T} H x)

    mode = match.arg(mode)
    stopifnot(!is.null(H)) 
    p = nrow(H)

    ## ensure x, if given, is a matrix. this is relevant if x is just one sample which can be
    ## given as a vector
    if (!is.null(x) && !is.matrix(x)) {
        x = matrix(c(x), ncol = length(x))
    }

    ## this is the general case involving rescaling and rotation
    eig = eigen(H)
    V = eig$vectors
    Lam.sqrt = diag(sqrt(eig$values))
    rep.Lam.sqrt = diag(1/sqrt(eig$values))

    ## x = Lam.sqrt %*% V^T z
    ## and z is such that  x^T H x = z^T z
    x.to.z = Lam.sqrt %*% t(V) 
    z.to.x = V %*% rep.Lam.sqrt

    if (mode == "r") {
        X = inla.pc.multvar.sphere.r(n, p = p, lambda = lambda,  h = h)
        for(i in 1:n) {
            z = t(X[i,, drop=FALSE])
            X[i, ] = z.to.x %*% z
        }
        return (X)
    } else if (mode == "d") {
        ljac = 0.5*sum(log(eig$values))
        n = nrow(x)
        ldens.x = numeric(n)
        for(i in 1:n) {
            z = x.to.z %*% t(x[i,, drop=FALSE])
            ldens = inla.pc.multvar.sphere.d(t(z), lambda = lambda, log = TRUE, h = h)
            ldens.x[i] = ldens + ljac
        }
        return (if (log) ldens.x else exp(ldens.x))
    }
    stop("Should not happen")
    return()
}

inla.pc.multvar.testing = function(p, lambda=1)
{
    ## test some of the functions
    library(cubature)

    upper = inla.pc.multvar.h.default(1/lambda + 6 * 1/lambda, inverse=TRUE)
    r = (adaptIntegrate(inla.pc.multvar.simplex.d,
                        lowerLimit = rep(0, p),
                        upperLimit = rep(upper, p),
                        lambda = lambda))
    print(paste("TEST 1: should be 1. result = ", r$integral))

    b = matrix(runif(p), ncol=1)
    b = b/max(b)
    r = (adaptIntegrate(inla.pc.multvar.simplex.general.d,
                        lowerLimit = rep(0, p),
                        upperLimit = rep(upper, p),
                        lambda = lambda, b=b))
    print(paste("TEST 2: should be 1. result = ", r$integral))

    r = (adaptIntegrate(inla.pc.multvar.sphere.d,
                        lowerLimit = rep(0, p),
                        upperLimit = rep(upper, p),
                        lambda = lambda))
    print(paste("TEST 3: should be ", 1/2^p, " result = ", r$integral))
    
    H = diag(runif(p))
    H = matrix(runif(p^2), p, p)
    H = H %*% t(H)
    H = H/max(H)
    print(H)
    library(R2Cuba)
    r = cuhre(p, 1, inla.pc.multvar.sphere.general.d,
        lower = rep(-upper+0.01, p),
        upper = rep(upper, p),
        lambda = lambda, H=H,
        rel.tol =1e-6, 
        flags = list(verbose=1))
    
    print(paste("TEST 4: should be ", 1, " result = ", r$value))
}
