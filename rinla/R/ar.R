## Export: inla.ar.pacf2phi
## Export: inla.ar.phi2pacf
## Export: inla.ar.pacf2acf
## Export: inla.ar.phi2acf

##!\name{inla.ar}
##!\alias{inla.ar.pacf2phi}
##!\alias{ar.pacf2phi}
##!\alias{pacf2phi}
##!\alias{inla.ar.phi2pacf}
##!\alias{ar.phi2pacf}
##!\alias{phi2pacf}
##!\alias{inla.ar.phi2acf}
##!\alias{ar.phi2acf}
##!\alias{phi2acf}
##!\alias{inla.ar.pacf2acf}
##!\alias{ar.pacf2acf}
##!\alias{pacf2acf}
##!
##!\title{Convert between parameterizations for the AR(p) model}
##!
##!\description{These functions convert between the AR(p) coefficients \code{phi}, 
##!             the partial autorcorrelation coefficients \code{pacf} and the
##!             autocorrelation function \code{acf}.
##!             The \code{phi}-parameterization is the same as used for \code{arima}-models in \code{R}; see \code{?arima}
##!             and the parameter-vector \code{a} in \code{Details}.}
##!\usage{
##!   inla.ar.pacf2phi(pac)
##!   inla.ar.phi2pacf(phi)
##!   inla.ar.pacf2acf(pac, lag.max = length(pac))
##!   inla.ar.phi2acf(phi, lag.max = length(phi))
##!}
##!
##!\arguments{
##!  \item{pac}{The partial autorcorrelation coefficients}
##!  \item{phi}{The AR(p) parameters \code{phi}}
##!  \item{lag.max}{The maximum lag to compute the ACF for}
##!}
##!\value{
##!  \code{inla.ar.pacf2phi}  returns \code{phi} for given \code{pacf}.
##!  \code{inla.ar.phi2pacf}  returns \code{pac} for given \code{phi}.
##!  \code{inla.ar.phi2acf}  returns \code{acf} for given \code{phi}.
##!  \code{inla.ar.pacf2acf}  returns \code{acf} for given \code{pacf}.
##!}
##!\author{Havard Rue \email{hrue@r-inla.org}}
##!\examples{
##! pac = runif(5)
##! phi = inla.ar.pacf2phi(pac)
##! pac2 = inla.ar.phi2pacf(phi)
##! print(paste("Error:", max(abs(pac2-pac))))
##! print("Correlation matrix (from pac)")
##! print(toeplitz(inla.ar.pacf2acf(pac)))
##! print("Correlation matrix (from phi)")
##! print(toeplitz(inla.ar.phi2acf(phi)))
##!}


## functions for the AR model, same as its C-versions in ar.c
inla.ar.pacf2phi = function(pac)
{
    ## I know, the R-coding is a bit weird as these are translated
    ## from arima.c in R...
    
    ## run the Durbin-Levinson recursions to find phi_{j.}, ( j = 2,
    ## ..., p and phi_{p.} are the autoregression coefficients
    p = length(pac)
    stopifnot(p > 0)
    phi = work = pac
    if (p > 1L) {
        for (j in 1L:(p-1L)) {
            a = phi[j+1L];
            phi[1L:j] = work[1L:j] = work[1L:j] - a * phi[j:1L]
        }
    }
    
    return(phi)
}
inla.ar.phi2pacf = function(phi)
{
    ## I know, the R-coding is a bit weird as these are translated
    ## from arima.c in R...

    ## Run the Durbin-Levinson recursions backwards to find the PACF
    ## phi_{j.} from the autoregression coefficients
    p = length(phi)
    stopifnot(p > 0)
    work = pac = phi
    
    if (p > 1L) {
        for(j in (p-1L):1L) {
            a = pac[j+1L];
            pac[1L:j] = work[1L:j]  = (pac[1L:j] + a * pac[j:1L]) / (1.0 - a^2);
        }
    }
    
    return (pac)
}
inla.ar.phi2acf = function(phi, lag.max = length(phi))
{
    ## return acf for given phi
    p = length(phi)
    stopifnot(p > 0)

    A = matrix(0, p, p)
    b = numeric(p)

    for(i in 1:p) {
        for(j in 1:p) {
            if (i == j) {
                A[i, j] = -1.0
            } else {
                lag = abs(i-j)
                A[i, lag] = phi[j] + A[i, lag]
            }
        }
        b[i] = -phi[i]
    }
    r = try(solve(A, b), silent = TRUE)
    ## if the model is singular,  then return nothing
    if (inherits(r, "try-error")) {
        return (numeric(0))
    }
    r = pmax(-1, pmin(1, r)) ## known to be true
    r = c(1, r)
    if (lag.max > p) {
        r = c(r, rep(0, lag.max-p))
        for(i in (p+1):(lag.max+1)) {
            r[i] = sum(phi * r[(i-1):(i-1-p+1)])
        }
        r = pmax(-1, pmin(1, r)) ## known to be true
    }

    return (r)
}
inla.ar.pacf2acf = function(pac, lag.max = length(pac))
{
    return (inla.ar.phi2acf(inla.ar.pacf2phi(pac), lag.max = lag.max))
}
