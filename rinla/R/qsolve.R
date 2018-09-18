## Export: inla.qsolve

##! \name{qsolve}
##! \alias{inla.qsolve}
##! \alias{qsolve}
##! 
##! \title{Solves linear SPD systems}
##! 
##! \description{This routine use the GMRFLib implementation
##!              to solve linear systems with a SPD matrix.}
##! \usage{
##!     inla.qsolve(Q, B, reordering, method = c("solve", "forward", "backward"))
##! }
##! 
##! \arguments{
##!   \item{Q}{A SPD matrix,  either as a (dense) matrix,  sparse-matrix  or a filename
##!            containing the matrix (in the fmesher-format).}
##!   \item{B}{The right hand side matrix, either as a (dense) matrix,  sparse-matrix  or a filename
##!            containing the matrix (in the fmesher-format).
##!            (Must be a matrix or sparse-matrix even if \code{ncol(B)} is 1.)}
##!   \item{reordering}{The type of reordering algorithm to be used for \code{TAUCS};
##!        either one of the names listed in \code{inla.reorderings()} 
##!        or the output from \code{inla.qreordering(Q)}.
##!        The default is "auto" which try several reordering algorithm and use the best one for this particular matrix.}
##!   \item{method}{The system to solve, one of \code{"solve"},  \code{"forward"} or \code{"backward"}. Let \code{Q = L L^T},
##!                 where \code{L} is lower triangular
##!                 (the Cholesky triangle),  then \code{method="solve"} solves \code{L L^T X = B} or
##!                  equivalently \code{Q X = B},  \code{method="forward"} solves \code{L X = B},  and
##!                 \code{method="backward"} solves \code{L^T X  = B}. }
##!}
##!\value{
##!  \code{inla.qsolve} returns a matrix \code{X},  which is the solution of \code{Q X = B},  \code{L X = B} or \code{L^T X = B} 
##!  depending on the value of \code{method}.                      
##!}
##!\author{Havard Rue \email{hrue@r-inla.org}}
##! 
##!\examples{
##!n = 10
##!QQ = matrix(runif(n^2), n, n)
##!Q = inla.as.dgTMatrix(QQ \%*\% t(QQ))
##!B = matrix(runif(n^2-n), n, n-1)
##!
##!X = inla.qsolve(Q, B, method = "solve")
##!print(paste("err", sum(abs( Q \%*\% X - B))))
##!
##!L = t(chol(Q))
##!X = inla.qsolve(Q, B, method = "forward")
##!print(paste("err", sum(abs( L \%*\% X - B))))
##!
##!X = inla.qsolve(Q, B, method = "backward")
##!print(paste("err", sum(abs( t(L) \%*\% X - B))))
##!
##!Q.file = INLA:::inla.write.fmesher.file(Q)
##!B.file = INLA:::inla.write.fmesher.file(B)
##!X = inla.qsolve(Q.file, B.file, method = "backward")
##!print(paste("err", sum(abs( t(L) \%*\% X - B))))
##!unlink(Q.file)
##!unlink(B.file)
##!}

`inla.qsolve` = function(Q, B, reordering = inla.reorderings(), method = c("solve", "forward", "backward"))
{
    smtp = match.arg(inla.getOption("smtp"), c("taucs", "band", "default", "pardiso"))
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
    
    ## ensure B is a dense matrix if its a sparse
    if (is(B, "Matrix")) {
        B = as.matrix(B)
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
    method = match.arg(method)

    Xfile = inla.tempfile()
    inla.set.sparselib.env(NULL)
    if (inla.os("linux") || inla.os("mac")) {
        s = system(paste(shQuote(inla.getOption("inla.call")), "-s -m qsolve", "-r", 
                         reordering, "-S", smtp, Qfile, Xfile, Bfile, method), intern=TRUE)
    } else if(inla.os("windows")) {
        s = system(paste(shQuote(inla.getOption("inla.call")), "-s -m qsolve", "-r",
                         reordering, "-S", smtp, Qfile, Xfile, Bfile, method), intern=TRUE)
    } else {
        stop("\n\tNot supported architecture.")
    }

    X = inla.read.fmesher.file(Xfile)

    if (Qremove) {
        unlink(Qfile)
    }
    if (Bremove) {
        unlink(Bfile)
    }
    
    return (X)
}
