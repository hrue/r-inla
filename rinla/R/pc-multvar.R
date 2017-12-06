## Export: inla.pc.multvar.h.default
## Export: inla.pc.multvar.simplex.d
## Export: inla.pc.multvar.simplex.r
## Export: inla.pc.multvar.sphere.d
## Export: inla.pc.multvar.sphere.r

##! \name{pc.multvar}
##! \alias{pc.multvar}
##! \alias{inla.pc.multvar}
##! \alias{inla.pc.multvar.simplex}
##! \alias{inla.pc.multvar.simplex.d}
##! \alias{inla.pc.multvar.simplex.r}
##! \alias{inla.pc.multvar.sphere}
##! \alias{inla.pc.multvar.sphere.d}
##! \alias{inla.pc.multvar.sphere.r}
##! \alias{inla.pc.multvar.h.default}
##! \alias{pc.multvar}
##! \alias{pc.multvar.simplex}
##! \alias{pc.multvar.simplex.d}
##! \alias{pc.multvar.simplex.r}
##! \alias{pc.multvar.sphere}
##! \alias{pc.multvar.sphere.d}
##! \alias{pc.multvar.sphere.r}
##! \alias{pc.multvar.h.default}
##!
##! \title{Multivariate PC priors}
##! 
##! \description{Functions to evaluate and simulate from 
##!              multivariate PC priors: The simplex and sphere case}
##! \usage{
##! inla.pc.multvar.h.default(x, inverse = FALSE, derivative = FALSE)
##! inla.pc.multvar.simplex.r(n = NULL, lambda = 1, h = inla.pc.multvar.h.default, b = NULL)
##! inla.pc.multvar.simplex.d(x = NULL, lambda = 1, log = FALSE, h = inla.pc.multvar.h.default, b = NULL)
##! inla.pc.multvar.sphere.r(n = NULL, lambda = 1, h = inla.pc.multvar.h.default, H = NULL)
##! inla.pc.multvar.sphere.d(x = NULL, lambda = 1, log = FALSE, h = inla.pc.multvar.h.default, H = NULL)
##! }
##! \arguments{
##!   \item{x}{Samples to evaluate. If input is a matrix then each row is a sample. If input is
##!            a vector then this is the sample.}
##!   \item{inverse}{Compute the inverse of the h()-function.}
##!   \item{derivative}{Compute the derivative of the h()-function. (derivative of the inverse
##!                     function is not used).}
##!   \item{n}{Number of samples to generate.}
##!   \item{lambda}{The lambda-parameter in the PC-prior.}
##!   \item{log}{Evaluate the density in log-scale or ordinary scale.}
##!   \item{h}{The h()-function,  defaults to \code{inla.pc.multvar.h.default}. See that code
##!            for an example of how to write a user-spesific function.}
##!   \item{b}{The b-vector (gradient) in the expression for the simplex option,  \code{d(xi) = h(b^T xi)}}
##!   \item{H}{The H(essian)-matrix in the expression for the sphere option,  \code{d(xi) =
##!            h(1/2 *xi^T H xi)}. If \code{H} is a vector, then it is interpreted as the
##!            diagonal of a (sparse) diagonal matrix.}
##!}
##! \details{
##!   These functions implements multivariate PC-priors of the simplex and sphere type.
##! }
##! \value{%%
##!     \code{inla.pc.multvar.simplex.r} generate samples from the simplex case, and
##!     \code{inla.pc.multvar.simplex.d} evaluate the density.
##!     \code{inla.pc.multvar.sphere.r} generate samples from the sphere case,  and
##!     \code{inla.pc.multvar.sphere.d} evaluate the density.
##!     \code{inla.pc.multvar.h.default} implements the default h()-function and illustrate how
##!     to code your own spesific one, if needed. 
##! }
##! \author{Havard Rue \email{hrue@r-inla.org}}
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

inla.pc.multvar.simplex.r.core = function(n, p, lambda = 1, h = inla.pc.multvar.h.default)
{
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
    colnames(x) = paste("x", inla.num(1:p), sep="")
    rownames(x) = paste("sample", inla.num(1:n), sep="")

    return (x)
}

inla.pc.multvar.simplex.d.core = function(x, lambda = 1, log = FALSE, 
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

inla.pc.multvar.simplex.r = function(
    n = NULL, lambda = 1, h = inla.pc.multvar.h.default, b = NULL)
{
    return (inla.pc.multvar.simplex.core(
        n = n, lambda = lambda, h = h, b = b, mode = "r"))
}

inla.pc.multvar.simplex.d = function(
    x = NULL, lambda = 1, log = FALSE, h = inla.pc.multvar.h.default, b = NULL)
{
    return (inla.pc.multvar.simplex.core(
        x = x, lambda = lambda, log = log, h = h, b = b, mode = "d"))
}

inla.pc.multvar.simplex.core = function(
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
        X = inla.pc.multvar.simplex.r.core(n, p = p, lambda = lambda,  h = h)
        for(i in 1:n) {
            X[i, ] = t(1/b) * X[i,, drop = FALSE]
        }
        return (X)
    } else if (mode == "d") {
        n = nrow(x)
        ldens.x = numeric(n)
        for(i in 1:n) {
            theta = c(b) * c(x[i,])
            ldens = inla.pc.multvar.simplex.d.core(theta, lambda = lambda, log = TRUE, h = h)
            ldens.x[i] = ldens + sum(log(abs(b)))
        }
        return (if (log) ldens.x else exp(ldens.x))
    }

    stop("Should not happen")
    return()
}

inla.pc.multvar.sphere.r.core = function(n, p, lambda = 1, h = inla.pc.multvar.h.default)
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
    colnames(x) = paste("x", inla.num(1:p), sep="")
    rownames(x) = paste("sample", inla.num(1:n), sep="")

    return (x)
}

inla.pc.multvar.sphere.d.core = function(x, lambda = 1, log = FALSE, 
    h = inla.pc.multvar.h.default)
{
    ## evaluate the density from the pc prior where d = h(1/2 * \sum x_i^2), x_i >=0.

    sphere.log.surface = function(n)
    {
        ## surface of the sphere in n-dimensions with radii=1
        ## n=0 line, n=1 circle, n=2 sphere etc
        return (log(n+1) +  (n+1)/2*log(pi) - lgamma((n+1)/2 +1))
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
            ldens = (log(lambda) - lambda * d + log(abs(h(0.5*r^2, derivative = TRUE)))
                     - (p-2)*log(r) - sphere.log.surface(p-1))
            return (ldens)
        })

    return (if (log) ldens else exp(ldens))
}

inla.pc.multvar.sphere.r = function(
    n = NULL, lambda = 1, h = inla.pc.multvar.h.default, H = NULL)
{
    return (inla.pc.multvar.sphere.core(
        n = n, lambda = lambda, h = h, H = H, mode = "r"))
}

inla.pc.multvar.sphere.d = function(
    x = NULL, lambda = 1, log = FALSE, h = inla.pc.multvar.h.default, H = NULL)
{
    return (inla.pc.multvar.sphere.core(
        x = x, lambda = lambda, log = log, h = h, H = H, mode = "d"))
}


inla.pc.multvar.sphere.core = function(
    x = NULL, n = NULL, lambda = 1, log = FALSE,
    h = inla.pc.multvar.h.default, H = NULL, mode = c("r", "d"))
{
    ## this is the case where d = h(1/2 x^{T} H x)

    mode = match.arg(mode)
    stopifnot(!is.null(H)) 

    ## ensure x, if given, is a matrix. this is relevant if x is just one sample which can be
    ## given as a vector
    if (!is.null(x) && !is.matrix(x)) {
        x = matrix(c(x), ncol = length(x))
    }
    ## if the class is ddiMatrix, convert it
    if (is(H, class(Diagonal(1)))) {
        H = diag(H)
    }

    if (is.vector(H)) {
        ## in case H is a diagonal matrix
        p = length(H)
        d = sort(H, decreasing = TRUE)
        ## create the reverse-diagonal, as this is how R compute the eigenvectors. In this case
        ## we will get exactly the same results as with H as a dense diagonal matrix. Copy this
        ## trick from lava::revdiag
        V = Diagonal(p, x = rep(0, p))
        V[cbind(rev(seq(p)), seq(p))] = rep(1, p)
        ##
        Lam.sqrt = Diagonal(p, x = sqrt(d))
        rep.Lam.sqrt = Diagonal(p, x = 1.0/sqrt(d))
        eig = list(values = d)
    } else {
        p = nrow(H)
        ## this case involving rescaling and rotation
        eig = eigen(H)
        V = eig$vectors
        Lam.sqrt = diag(sqrt(eig$values), ncol = p)
        rep.Lam.sqrt = diag(1/sqrt(eig$values), ncol = p)
    }

    ## x = Lam.sqrt %*% V^T z
    ## and z is such that  x^T H x = z^T z
    x.to.z = Lam.sqrt %*% t(V) 
    z.to.x = V %*% rep.Lam.sqrt

    if (mode == "r") {
        X = inla.pc.multvar.sphere.r.core(n, p = p, lambda = lambda,  h = h)
        for(i in 1:n) {
            z = t(X[i,, drop=FALSE])
            X[i, ] = as.numeric(z.to.x %*% z)
        }
        return (X)
    } else if (mode == "d") {
        ljac = 0.5*sum(log(eig$values))
        n = nrow(x)
        ldens.x = numeric(n)
        for(i in 1:n) {
            z = x.to.z %*% t(x[i,, drop=FALSE])
            ldens = inla.pc.multvar.sphere.d.core(t(z), lambda = lambda, log = TRUE, h = h)
            ldens.x[i] = ldens + ljac
        }
        return (if (log) ldens.x else exp(ldens.x))
    }
    stop("Should not happen")
    return()
}
