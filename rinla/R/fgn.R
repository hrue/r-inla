## Export: inla.fgn

##! \name{fgn}
##! \alias{inla.fgn}
##! \alias{fgn}
##!
##! \title{Lookup coefficients in the 3-component AR(1) mixture representing FGN(H)}
##!
##! \description{This function will lookup the coefficients in the 3-component AR(1)
##!              mixture representing FGN(H)}
##!
##! \usage{
##!     inla.fgn(H)
##! } 
##!
##! \arguments{
##!   \item{H}{The Hurst coeffcient (0<H<1)}
##!  }
##! \value{
##!  A named vector of length 7, where the first element is \code{H},
##!  the next three are the lag one correlations (or phi's),
##!  and the last three are the weights.
##!
##!  This function is EXPERIMENTAL!!!
##! }
##! \author{Havard Rue \email{hrue@r-inla.org}}
##!
##! \examples{
##!     inla.fgn(0.7)
##! } 

`inla.fgn` = function(H)
{
    stopifnot(is.numeric(H) && length(H) == 1)
    out.file = inla.tempfile()
    if (inla.os("linux") || inla.os("mac")) {
        s = system(paste(shQuote(inla.getOption("inla.call")), "-s -m fgn", as.character(H), out.file), intern=TRUE)
    } else if(inla.os("windows")) {
        s = system(paste(shQuote(inla.getOption("inla.call")), "-s -m fgn", as.character(H), out.file), intern=TRUE)
    } else {
        stop("\n\tNot supported architecture.")
    }

    res = c(inla.read.fmesher.file(out.file))
    k = (length(res)-1) %/% 2L
    names(res) = c("H",  paste0("phi", 1:k), paste0("weight", 1:k))
    unlink(out.file)

    return (res)
}
