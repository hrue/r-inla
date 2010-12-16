## this function interface with mode = finn in the inlaprogram, which
## computes whatever Finn wants...

`inla.finn` = function(C,
        reordering = c("auto", "default", "identity", "band", "metis", "genmmd", "amd", "md", "mmd"))
{
    reordering = match.arg(reordering)

    if (is.matrix(C) || is(C, "dgTMatrix")) {
        finn.file = inla.sparse2file(inla.matrix2sparse(C), c.indexing = TRUE)
        remove = TRUE
    } else if (is.character(C)) {
        finn.file = C
        remove = FALSE
    } else {
        finn.file = inla.sparse2file(C, c.indexing = TRUE)
        remove = TRUE
    }
        
    if (inla.os("linux") || inla.os("mac")) {
        s = system(paste(shQuote(inla.getOption("inla.call")), "-s -m finn",
                "-r", reordering, finn.file), intern=TRUE)
    } else if(inla.os("windows")) {
        s = system(paste(shQuote(inla.getOption("inla.call")), "-s -m finn",
                "-r", reordering, finn.file), intern=TRUE)
    } else {
        stop("\n\tNot supported architecture.")
    }

    if (remove) {
        unlink(finn.file)
    }

    s = as.numeric(s)
    n = length(s) %/% 2L
    return ( list(sample = s[1L:n], reordering = s[n + (1L:n)] + 1L))
}
