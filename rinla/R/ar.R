##!\name{inla.ar}
##!\alias{inla.ar.pacf2phi}
##!\alias{ar.pacf2phi}
##!\alias{pacf2phi}
##!\alias{inla.ar.phi2pacf}
##!\alias{ar.phi2pacf}
##!\alias{phi2pacf}
##!
##!\title{Convert between parameterizations for the AR(p) model}
##!
##!\description{These functions convert between the AR(p) coefficients \code{phi} and
##!             the partial autorcorrelation coefficients \code{pacf}.
##!             The \code{phi}-parameterization is the same as used for \code{arima}-models in \code{R}; see \code{?arima}
##!             and the parameter-vector \code{a} in \code{Details}.}
##!\usage{
##!   phi = inla.ar.pacf2phi(pac)
##!   pac = inla.ar.phi2pacf(phi)
##!}
##!
##!\arguments{
##!  \item{pac}{The partial autorcorrelation coefficients}
##!  \item{phi}{The AR(p) parameters \code{phi}}
##!}
##!\value{
##!  \code{inla.ar.pacf2phi}  returns \code{phi} for given \code {pacf}.
##!  \code{inla.ar.phi2pacf}  returns \code{pac} for given \code {phi}.
##!}
##!\author{Havard Rue \email{hrue@math.ntnu.no}}
##!\examples{
##! pac = runif(5)
##! phi = inla.ar.pacf2phi(pac)
##! pac2 = inla.ar.phi2pac(phi)
##! print(paste("Error:", max(abs(pac2-pac))))
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
    for (j in 1L:(p-1L)) {
        a = phi[j+1L];
        for (k in 0:(j-1L)) {
            work[k+1L] = work[k+1L] - a * phi[j - k];
        }
        phi[1:j] = work[1:j]
    }

    return(phi)
}
inla.ar.phi2pacf = function(phi)
{
    ## I know, the R-coding is a bit weird as these are translated
    ## from arima.c in R...

    ## Run the Durbin-Levinson recursions backwards to find the PACF
    ## phi_{j.} from the autoregression coefficients
    p = length(pac)
    stopifnot(p > 0)
    work = pac = phi

    for(j in (p-1L):1L) {
        a = pac[j+1L];
        for(k in 0L:(j-1L)) {
            work[k+1L]  = (pac[k+1L] + a * pac[j - k]) / (1.0 - a^2);
        }
        pac[1L:j] = work[1L:j];
    }
    
    return (pac)
}
