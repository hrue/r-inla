## function for computing the marginals variances of a sparse matrix

`inla.qinv` = function(C)
{
    if (is.matrix(C) || is(C, "dgTMatrix")) {
        qinv.file = inla.sparse2file(inla.matrix2sparse(C), c.indexing = TRUE)
        remove = TRUE
    } else if (is.character(C)) {
        qinv.file = C
        remove = FALSE
    } else {
        qinv.file = inla.sparse2file(C, c.indexing = TRUE)
        remove = TRUE
    }
        
    if (inla.os("linux") || inla.os("mac")) {
        s = system(paste(shQuote(inla.getOption("inla.call")), "-s -m qinv", qinv.file), intern=TRUE)
    } else if(inla.os("windows")) {
        s = system(paste(shQuote(inla.getOption("inla.call")), "-s -m qinv", qinv.file), intern=TRUE)
    } else {
        stop("\n\tNot supported architecture.")
    }

    if (remove) {
        unlink(qinv.file)
    }

    return (as.numeric(s))
}

