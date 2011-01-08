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























