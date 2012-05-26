##! \name{qsolve}
##! \alias{inla.qsolve}
##! \alias{qsolve}
##! 
##! \title{Solves linear SPD systems}
##! 
##! \description{This routine use the GMRFLib implementation
##!              to solve linear systems with a SPD matrix.}
##! \usage{
##!     inla.qsolve(Q, B, reordering)
##! }
##! 
##! \arguments{
##!   \item{Q}{A SPD matrix,  either as a (dense) matrix,  sparseMatrix  or a filename
##!            containing the matrix.}
##!   \item{B}{The right hand side matrix. (Must be a matrix even if ncol is 1.)}
##!   \item{reordering}{The type of reordering algorithm to be used,  one of
##!        "auto", "default", "identity", "band", "metis", "genmmd", "amd", "amdbar", "md", "mmd", "amdc" and "amdbarc", or
##!        the output from \code{inla.qreordering(Q)}.
##!        The default is "auto" which try several reordering algorithm and use the best one for this particular matrix.}
##!}
##!\value{
##!  \code{inla.qsolve} returns a matrix A,  which is the solution of QA=B.
##!}
##!\author{Havard Rue \email{hrue@math.ntnu.no}}
##! 
##!\examples{
##! n = 5
##! Q = matrix(runif(n^2), n, n)
##! Q = inla.as.dgTMatrix(Q \%*\% t(Q))
##! B = matrix(runif(n^2 -n), n, n-1)
##! A = inla.qsolve(Q, B)
##!}

`inla.qsolve` = function(Q, B, reordering = inla.reorderings())
{
    Q = inla.sparse.check(Q)
    if (is(Q, "dgTMatrix")) {
        Qfile = inla.write.fmesher.file(Q)
        Qremove = TRUE
    } else if (is.character(Q)) {
        Qfile = Q
        Qremove = FALSE
    } else {
        stop("This should not happen.")
    }
    
    if (is.matrix(B)) {
        Bfile = inla.write.fmesher.file(B)
        Bremove = TRUE
    } else if (is.character(B)) {
        Bfile = B
        Bremove = FALSE
    } else {
        stop("This should not happen.")
    }
        
    if (is.list(reordering)) {
        ## argument is the output from inla.qreordering()
        reordering = reordering$name
    }
    reordering = match.arg(reordering)

    Afile = inla.tempfile()
    if (inla.os("linux") || inla.os("mac")) {
        s = system(paste(shQuote(inla.getOption("inla.call")), "-s -m qsolve", "-r",  reordering, Qfile, Afile, Bfile), intern=TRUE)
    } else if(inla.os("windows")) {
        s = system(paste(shQuote(inla.getOption("inla.call")), "-s -m qsolve", "-r",  reordering, Qfile, Afile, Bfile), intern=TRUE)
    } else {
        stop("\n\tNot supported architecture.")
    }

    A = inla.read.fmesher.file(Afile)

    if (Qremove) {
        unlink(Qfile)
    }
    if (Bremove) {
        unlink(Bfile)
    }
    
    return (A)
}
