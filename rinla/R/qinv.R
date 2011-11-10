##! \name{qinv}
##! \alias{inla.qinv}
##! \alias{qinv}
##! 
##! \title{Computes (parts of) the inverse of a SPD sparse matrix}
##! 
##! \description{This routine use the GMRFLib implementation
##!              which compute parts of the inverse of a SPD sparse matrix.
##!              The diagonal and values for the neighbours in the inverse, are provided.}
##! 
##! \usage{
##!     inla.qinv(Q)
##! }
##! 
##! \arguments{
##! 
##!   \item{Q}{A SPD matrix,  either as a (dense) matrix,  sparseMatrix,  or
##!           a (ascii-)filename with entries in the following format \code{i j Qij}.}
##!  }
##! \value{
##!   \code{inla.qinv} returns a \code{sparseMatrix} of type \code{dgTMatrix} with the
##!   diagonal and values for the neigbours in the inverse. Note that the full inverse is NOT provided!
##! }
##! \author{Havard Rue \email{hrue@math.ntnu.no}}
##! 
##! \examples{
##! 
##! ## dense matrix example
##! n = 10
##! A = matrix(runif(n^2), n, n)
##! Q = A %*% t(A)
##! print(mean(abs(inla.qinv(Q) - solve(Q))))
##! 
##! ## sparse matrix example
##! rho = 0.9
##! Q = toeplitz(c(1+rho^2, -rho,  rep(0, n-3), -rho)) / (1-rho^2)
##! Q = inla.as.dgTMatrix(Q)
##! Q.inv = inla.qinv(Q)
##! 
##! ## compute the marginal variances as a vector from a precision matrix
##! marginal.variances = diag(inla.qinv(Q))
##! 
##! ## read the sparse matrix from a file in the 'i, j, value' format
##! filename = inla.tempfile()
##! write(t(cbind(Q@i+1L,  Q@j+1L,  Q@x)), ncol=3, file=filename)
##! Qinv = inla.qinv(filename)
##! unlink(filename)
##! }

`inla.qinv` = function(C)
{
    if (!is(C, "dgTMatrix")) {
        C = inla.sparse.check(C)
    }
    
    if (is(C, "dgTMatrix")) {
        qinv.file = inla.sparse2file(C, c.indexing = TRUE, symmetric = TRUE)
        remove = TRUE
    } else if (is.character(C)) {
        qinv.file = C
        remove = FALSE
    } else {
        stop("This chould not happen.")
    }
        
    out.file = inla.tempfile()
    if (inla.os("linux") || inla.os("mac")) {
        s = system(paste(shQuote(inla.getOption("inla.call")), "-s -m qinv", qinv.file, out.file), intern=TRUE)
    } else if(inla.os("windows")) {
        s = system(paste(shQuote(inla.getOption("inla.call")), "-s -m qinv", qinv.file, out.file), intern=TRUE)
    } else {
        stop("\n\tNot supported architecture.")
    }

    Qinv = inla.read.fmesher.file(out.file)

    if (remove) {
        unlink(qinv.file)
    }
    unlink(out.file)
    
    return (Qinv)
}
