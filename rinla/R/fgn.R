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
##!   \item{H}{The Hurst coeffcient (0<H<1),  or a vector of those}
##!  }
##! \value{
##!  A named matrix with 7 columns, where the first column is \code{H},
##!  column 2,  3 and 4 are lag one correlations (or phi's),
##!  and column 5,  6 and 7 are the weights.
##!
##!  This function is EXPERIMENTAL!!!
##! }
##! \author{Havard Rue \email{hrue@r-inla.org}}
##!
##! \examples{
##!     r = c(inla.fgn(0.7))
##!     r_m = inla.fgn(seq(0.6, 0.8, by=0.01))
##! } 

`inla.fgn` = function(H)
{
    in.file = inla.tempfile()
    out.file = inla.tempfile()
    inla.write.fmesher.file(matrix(as.numeric(H), ncol = 1), file = in.file)
    if (inla.os("linux") || inla.os("mac")) {
        s = system(paste(shQuote(inla.getOption("inla.call")), "-s -m fgn", in.file, out.file), intern=TRUE)
    } else if(inla.os("windows")) {
        s = system(paste(shQuote(inla.getOption("inla.call")), "-s -m fgn", in.file, out.file), intern=TRUE)
    } else {
        stop("\n\tNot supported architecture.")
    }

    res = inla.read.fmesher.file(out.file)
    k = (ncol(res)-1L) %/% 2L
    colnames(res) = c("H",  paste0("phi", 1:k), paste0("weight", 1:k))
    unlink(in.file)
    unlink(out.file)

    return (res)
}
