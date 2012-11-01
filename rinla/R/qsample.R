##! \name{qsample}
##! \alias{inla.qsample}
##! \alias{qsample}
##! 
##! \title{Generate samples from a GMRF using the GMRFLib implementation}
##! 
##! \description{This function generate samples from a GMRF using the GMRFLib implementation}
##! \usage{
##!     inla.qsample(n, Q, reordering = "auto",  seed = 0L, logdens = FALSE)
##! }
##! 
##! \arguments{
##!   \item{n}{Number of samples}
##!   \item{Q}{The precision matrix or a filename containing it.}
##!   \item{reordering}{The type of reordering algorithm to be used; either one of the names listed in \code{inla.reorderings()}
##!        or the output from \code{inla.qreordering(Q)}.
##!        The default is "auto" which try several reordering algorithm and use the best one for this particular matrix.}
##!   \item{seed}{The seed-option to control the RNG. If \code{seed=0L} then GMRFLib will set the seed intelligently.
##!               If \code{seed < 0L}  then the saved state of the RNG should be reused if possible,  if not
##!               then the behaviour for \code{seed=0L} is used. If \code{seed > 0L} then this is used as the seed for the RNG.}
##!   \item{logdens}{If \code{TRUE}, compute also the log-density of each sample.}
##!}
##!\value{
##!  If \code{logdens} is \code{FALSE} (default),  then \code{inla.qsample} returns 
##!  the samples in a matrix,  where each column is a sample. 
##!  If \code{logdens} is \code{TRUE}, then a list 
##!   with names \code{sample} and
##!  \code{logdens} is returned. The samples are stored in the matrix
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
##! matplot(x)
##! x = inla.qsample(10, G, logdens=TRUE)
##! matplot(x$sample)
##!}

`inla.qsample` = function(n = 1L, Q, reordering = inla.reorderings(), seed = 0L, logdens = FALSE)
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
    rng.file = inla.tempfile()

    envir = inla.get.inlaEnv()
    if (seed < 0L) {
        if (!exists("GMRFLib.rng.state", envir = envir)) {
            seed = 0L
        } else {
            rng.state = get("GMRFLib.rng.state", envir = envir)
            fp = file(rng.file, "wb")
            writeBin(as.raw(rng.state), fp)
            close(fp)
        }
    }

    if (inla.os("linux") || inla.os("mac")) {
        s = system(paste(shQuote(inla.getOption("inla.call")), "-s -m qsample", 
                "-r", reordering, "-z", seed, Q.file, x.file, n, rng.file), intern=TRUE)
    } else if(inla.os("windows")) {
        s = system(paste(shQuote(inla.getOption("inla.call")), "-s -m qsample",
                "-r", reordering, "-z", seed, Q.file, x.file, n, rng.file), intern=TRUE)
    } else {
        stop("\n\tNot supported architecture.")
    }

    if (remove) {
        unlink(Q.file)
    }

    fp = file(rng.file, "rb")
    siz = file.info(rng.file)$size
    rng.state = readBin(fp, raw(), siz)
    close(fp)    
    assign("GMRFLib.rng.state", rng.state, envir = envir)
    unlink(rng.file)

    x = inla.read.fmesher.file(x.file)
    unlink(x.file)

    nx = dim(Q)[1L]
    samples = matrix(x[-(nx + 1L),, drop=FALSE], nx, n)
    colnames(samples) = paste("sample", 1L:n, sep="")
    rownames(samples) = paste("x", 1L:nx, sep="")
    ld = c(x[nx+1L, ])
    names(ld) = paste("logdens", 1L:n, sep="")

    if (logdens) {
        return (list(sample=samples, logdens = ld))
    } else {
        return (samples)
    }
}
