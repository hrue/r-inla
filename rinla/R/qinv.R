## function for computing the marginals variances of a sparse matrix

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
