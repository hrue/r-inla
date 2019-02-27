## Export: inla.fgn

##! \name{fgn}
##! \alias{inla.fgn}
##! \alias{fgn}
##!
##! \title{Return the coefficients in the 3-component AR(1) mixture representing FGN(H)}
##!
##! \description{This function will return the coefficients in the 3-component AR(1)
##!              mixture representing FGN(H)}
##!
##! \usage{
##!     inla.fgn(H, K=4L, lag.max = NULL, approx = TRUE)
##! } 
##!
##! \arguments{
##!   \item{H}{The Hurst coeffcient (0<H<1),  or a vector of those}
##!   \item{K}{The number of components in representation,  must be 3L or 4L}
##!   \item{lag.max}{Integer. If positive integer, return the coeffcients implicitely as the ACF
##!                  from 0 to \code{lag.max}}
##!   \item{approx}{Logical. If \code{lag.max} is an positive integer and \code{approx} is \code{FALSE},
##!                   then return the true ACF instead of the approximated one.}
##!  }
##! \value{
##!  \code{inla.fgn} returns a named matrix. 
##!  If \code{is.null(lag.max)},  then 
##!  first column is \code{H},
##!  columns \code{1+1:K} are lag one correlations (or phi's),
##!  and columns \code{1+K+1:K} are the weights.
##!  If \code{lag.max > 0},  then return the ACFs in columns \code{2+(0:lag.max)},
##!  for the H in column 1,  either the approximated ones or the the true ones.
##!
##!  This function is EXPERIMENTAL!!!
##! }
##! \author{Havard Rue \email{hrue@r-inla.org}}
##!
##! \examples{
##!     r = c(inla.fgn(0.7))
##!     r_m = inla.fgn(seq(0.6, 0.8, by=0.01))
##! } 

`inla.fgn` = function(H, K=4L, lag.max = NULL, approx = TRUE)
{
    if (!any(K == c(3L, 4L))) {
        stop(paste0("Number of components 'K' must be 3 or 4,  not ",  K))
    }
    
    t.dir = inla.tempdir()
    in.file = inla.tempfile(tmpdir = t.dir)
    out.file = inla.tempfile(tmpdir = t.dir)
    inla.write.fmesher.file(matrix(c(K, as.numeric(H)), ncol = 1), file = in.file)
    if (inla.os("linux") || inla.os("mac")) {
        s = system(paste(shQuote(inla.getOption("inla.call")), "-s -m fgn", in.file, out.file), intern=TRUE)
    } else if(inla.os("windows")) {
        s = system(paste(shQuote(inla.getOption("inla.call")), "-s -m fgn", in.file, out.file), intern=TRUE)
    } else {
        stop("\n\tNot supported architecture.")
    }

    res = inla.read.fmesher.file(out.file)
    stopifnot(ncol(res) == 2*K+1)
    colnames(res) = c("H",  paste0("phi", 1:K), paste0("weight", 1:K))

    unlink(t.dir, recursive = TRUE)

    if (!is.null(lag.max)) {
        lag.max = as.integer(lag.max)
        stopifnot(lag.max > 0)
        m = length(H)
        n = lag.max + 1
        ACF = matrix(NA, m, n)
        for(i in 1:m) {
            if (approx) {
                phi = res[i, 1+1:K]
                w = res[i, 1+K+1:K]
                a = rep(0, n)
                for(j in 1:K) {
                    a = a + w[j] * phi[j]^(0:(n-1))
                }
            } else {
                a = inla.acvfFGN(H[i], n-1)
            }
            ACF[i, ] = a
        }
        ACF[, 1] = 1.0 ## just to make sure
        colnames(ACF) = paste0("acf", (1:ncol(ACF))-1)
        return (cbind(H=H, ACF))
    } else {
        return (res)
    }
}

`inla.acvfFGN` = function(H, maxlag) 
{
    ## a copy of function FGN::acvfFGN as library FGN in not there for R-3.5
    h2 <- 2 * H
    k <- 1:maxlag
    return (c(1, 0.5 * ((k + 1)^h2 - 2 * k^h2 + (k - 1)^h2)))
}
