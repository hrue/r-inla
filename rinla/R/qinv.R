## function for computing the marginals variances of a sparse matrix

`inla.qinv` = function(C)
{
    if (is.matrix(C)) {
        qinv.file = inla.sparse2file(inla.matrix2sparse(C), c.indexing = TRUE)
    } else {
        qinv.file = inla.sparse2file(C, c.indexing = TRUE)
    }
        
    if (inla.os("linux") || inla.os("mac")) {
        s = system(paste(shQuote(inla.getOption("inla.call")), "-s -m qinv", qinv.file), intern=TRUE)
    } else if(inla.os("windows")) {
        s = system(paste(shQuote(inla.getOption("inla.call")), "-s -m qinv", qinv.file), intern=TRUE)
    } else {
        stop("\n\tNot supported architecture.")
    }

    unlink(qinv.file)

    return (s)
}

