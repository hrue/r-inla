## this function interface with mode = finn in the inlaprogram, which
## computes whatever Finn wants...

`inla.finn` = function(Q, reordering = inla.reorderings(), seed = 0L)
{
    reordering = match.arg(reordering)

    Q = inla.sparse.check(Q)

    if (is(Q, "dgTMatrix")) {
        finn.file = inla.sparse2file(Q, c.indexing = TRUE)
        remove = TRUE
    } else if (is.character(Q)) {
        finn.file = Q
        remove = FALSE
    } else {
        stop("This should not happen.")
    }

    if (inla.os("linux") || inla.os("mac")) {
        s = system(paste(shQuote(inla.getOption("inla.call")), "-s -m finn", 
                "-r", reordering, "-z", seed, finn.file), intern=TRUE)
    } else if(inla.os("windows")) {
        s = system(paste(shQuote(inla.getOption("inla.call")), "-s -m finn",
                "-r", reordering, "-z", seed, finn.file), intern=TRUE)
    } else {
        stop("\n\tNot supported architecture.")
    }

    if (remove) {
        unlink(finn.file)
    }

    s = as.numeric(s)
    n = length(s) %/% 2L

    r = s[n + (1L:n)] + 1L
    ir = r
    ir[r] = 1:length(r)

    return ( list(sample = s[1L:n], reordering = r, ireordering = ir) )
}
