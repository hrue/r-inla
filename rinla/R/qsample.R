##! \name{qsample}
##! \alias{inla.qsample}
##! \alias{qsample}
##! 
##! \title{Generate samples from a GMRF using the GMRFLib implementation}
##! 
##! \description{This function generate samples from a GMRF using the GMRFLib implementation}
##! \usage{
##!     inla.qsample(n, Q, reordering = "auto",  seed = 0L)
##! }
##! 
##! \arguments{
##!   \item{n}{Number of samples}
##!   \item{Q}{The precision matrix or a filename containing it.}
##!   \item{reordering}{The type of reordering algorithm to be used; either one of the names listed in \code{inla.reorderings()}
##!        or the output from \code{inla.qreordering(Q)}.
##!        The default is "auto" which try several reordering algorithm and use the best one for this particular matrix.}
##!   \item{seed}{The seed to be used,  where \code{seed=0L} means that GMRFLib should decide the seed.}
##!}
##!\value{
##!  \code{inla.qsample} returns a list with names \code{sample} and
##!  \code{logdens}. The samples are stored in the matrix
##!  \code{sample} where each column is a sample, and the log
##!  densities of each sample are stored the vector \code{logdens}.
##!}
##!\author{Havard Rue \email{hrue@math.ntnu.no}}
##! 
##!\examples{
##! g = system.file("demodata/germany.graph", package="INLA")
##! G = inla.graph2matrix(g)
##! diag(G) = dim(G)[1L]
##! x = inla.qsample(10, G)
##! matplot(x$sample)
##!}

`inla.qsample` = function(n = 1L, Q, reordering = inla.reorderings(), seed = 0L)
{
    stopifnot(!missing(Q))
    stopifnot(n >= 1L)

    if (is.list(reordering)) {
        ## argument is the output from inla.qreordering()
        reordering = reordering$name
    }
    reordering = match.arg(reordering)

    Q = inla.sparse.check(Q)
    if (is(Q, "dgTMatrix")) {
        Q.file = inla.write.fmesher.file(Q)
        remove = TRUE
    } else if (is.character(Q)) {
        Q.file = Q
        remove = FALSE
    } else {
        stop("This should not happen.")
    }

    x.file = inla.tempfile()
    if (inla.os("linux") || inla.os("mac")) {
        s = system(paste(shQuote(inla.getOption("inla.call")), "-s -m qsample", 
                "-r", reordering, "-z", seed, Q.file, x.file, n), intern=TRUE)
    } else if(inla.os("windows")) {
        s = system(paste(shQuote(inla.getOption("inla.call")), "-s -m qsample",
                "-r", reordering, "-z", seed, Q.file, x.file, n), intern=TRUE)
    } else {
        stop("\n\tNot supported architecture.")
    }

    if (remove) {
        unlink(Q.file)
    }

    x = inla.read.fmesher.file(x.file)
    unlink(x.file)

    nx = dim(Q)[1L]
    samples = matrix(x[-(nx + 1L),, drop=FALSE], nx, n)
    colnames(samples) = paste("sample", 1L:n, sep="")
    rownames(samples) = paste("x", 1L:nx, sep="")
    logdens = x[nx+1L,, drop=TRUE]
    names(logdens) = paste("logdens", 1L:n, sep="")
    result = list(sample=samples, logdens = logdens)

    return (result)
}
