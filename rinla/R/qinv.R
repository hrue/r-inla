## Export: inla.qinv

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
##!     inla.qinv(Q, constr, reordering = INLA::inla.reorderings())
##! }
##!
##! \arguments{
##!
##!   \item{Q}{A SPD matrix,  either as a (dense) matrix or sparseMatrix.}
##!   \item{constr}{Optional linear constraints; 
##!                 see \code{?INLA::f} and argument \code{extraconstr}}
##!   \item{reordering}{The type of reordering algorithm to be used for \code{TAUCS};
##!        either one of the names listed in \code{inla.reorderings()} 
##!        or the output from \code{inla.qreordering(Q)}.
##!        The default is "auto" which try several reordering algorithm and use the best one for this particular matrix.}
##!  }
##! \value{
##!   \code{inla.qinv} returns a \code{sparseMatrix} of type \code{dgTMatrix} with the
##!   diagonal and values for the neigbours in the inverse. Note that the full inverse is NOT provided!
##! }
##! \author{Havard Rue \email{hrue@r-inla.org}}
##!
##! \examples{
##!
##! ## dense matrix example
##! n = 10
##! A = matrix(runif(n^2), n, n)
##! Q = A \%*\% t(A)
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
##! filename = INLA:::inla.tempfile()
##! write(t(cbind(Q@i+1L,  Q@j+1L,  Q@x)), ncol=3, file=filename)
##! Qinv = inla.qinv(filename)
##! unlink(filename)
##! }

`inla.qinv` = function(Q, constr, reordering = INLA::inla.reorderings())
{
    smtp = match.arg(inla.getOption("smtp"), c("taucs", "band", "default", "pardiso"))
    num.threads = inla.getOption("num.threads")

    Q = inla.sparse.check(Q)
    if (is(Q, "dgTMatrix")) {
        qinv.file = inla.write.fmesher.file(Q)
        remove = TRUE
    } else if (is.character(Q)) {
        qinv.file = Q
        remove = FALSE
    } else {
        stop("This chould not happen.")
    }

    constr.file = inla.tempfile()
    if (!missing(constr) && !is.null(constr)) {
        stopifnot(is.list(constr))
        A = as.matrix(constr$A)
        e = as.numeric(constr$e)
        stopifnot(ncol(A) == ncol(Q))
        stopifnot(nrow(A) == length(e))
        xx = matrix(c(nrow(A), c(A), c(e)), ncol = 1)
        inla.write.fmesher.file(xx, filename = constr.file)
    } 

    if (is.list(reordering)) {
        ## argument is the output from inla.qreordering()
        reordering = reordering$name
    }
    reordering = match.arg(reordering)

    out.file = inla.tempfile()
    inla.set.sparselib.env(NULL)
    if (inla.os("linux") || inla.os("mac")) {
        s = system(paste(shQuote(inla.getOption("inla.call")), "-s -m qinv",
                         "-r",  reordering, "-t", num.threads, 
                         "-S", smtp, qinv.file, constr.file, out.file), intern=TRUE)
    } else if(inla.os("windows")) {
        s = system(paste(shQuote(inla.getOption("inla.call")), "-s -m qinv",
                         "-r",  reordering, "-t", num.threads, 
                         "-S", smtp, qinv.file, constr.file, out.file), intern=TRUE)
    } else {
        stop("\n\tNot supported architecture.")
    }

    Qinv = inla.read.fmesher.file(out.file)

    if (remove) {
        unlink(qinv.file)
    }
    unlink(out.file)
    unlink(constr.file)

    return (Qinv)
}
