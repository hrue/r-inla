## Nothing to Export

`inla.fmesher.call` <- function(...) {
    stop("This feature is disabled and 'fmesher' is moved to a separate package.")
}

`inla.fmesher.make.prefix` <- function(dir = tempdir(), prefix = NULL) {
    if (is.null(prefix)) {
        if (is.null(dir)) {
            prefix <- paste(tempfile(), ".", sep = "")
        } else {
            inla.dir.create(dir)
            prefix <- paste(tempfile(tmpdir = dir), ".", sep = "")
        }
    } else {
        inla.dir.create(inla.dirname(prefix))
    }
    return(prefix)
}
