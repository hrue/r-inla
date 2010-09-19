### Binary I/O and interpretation of vectors

`inla.interpret.vector` = function(xx, debug=FALSE)
{
    if (is.null(xx))
        return (NULL)
    
    if (debug)
        cat("\tEnter inla.interpret.vector ... length(xx) = ", length(xx), "\n")

    len = length(xx)
    ii = 1L
    res = numeric()
    m = 0L
    np = NA
    no.missing = FALSE

    np = as.integer(xx[2L])
    stopifnot((length(xx) %% (2L*np+2)) == 0L)
    m = 2L * length(xx) %/% (2L*np+2)

    res.ncol = m
    res = matrix(0, np, res.ncol)

    m = 1L
    ii = 1L
    np2 = 2L*np + 2L
    for(j in 1:(res.ncol %/% 2)) {
        idx = seq.int(ii + 2L, by=2L, len=np)
        res[, c(m, m+1L)] = xx[c(idx, idx+1L)]
        m = m + 2L
        ii = ii + np2
    }

    stopifnot(length(xx) == ii-1L)
    stopifnot(res.ncol == m-1L) ## since we start from 1L

    if (debug)
        cat("\tLeave inla.interpret.vector ...\n")

    return(res)
}

`inla.interpret.vector.idx` = function(xx, debug=FALSE)
{
    if (is.null(xx))
        return (NULL)
    
    if (debug)
        cat("\tEnter inla.interpret.vector.idx ... length(xx) = ", length(xx), "\n")

    len = length(xx)
    ii = 1L
    res = numeric()
    m = 0L
    np = NA
    no.missing = FALSE
    
    np = as.integer(xx[2L])
    stopifnot((length(xx) %% (2L*np + 2L)) == 0L)
    m = 3L * length(xx) %/% (2L*np+2L)
    res.ncol = m
    res = matrix(0, np, res.ncol)

    m = 1L
    ii = 1L
    np2 = 2L * np + 2L
    for(j in 1:(m %/% 3L)) {
        idx = as.integer(xx[ii])
        jj = seq.int(ii + 2L, by=2L, len = np)
        res[, m] = idx
        res[, c(m+1L, m+2L)] = xx[c(jj, jj+1L)]
        m = m + 3L
        ii = ii + np2
    }
    stopifnot(length(xx) == ii-1L)
    stopifnot(res.ncol == m-1L) ## since we start from 1L

    if (debug)
        cat("\tLeave inla.interpret.vector.idx ...\n")

    return(res)
}

`inla.interpret.vector.list` = function(xx, debug=FALSE)
{
    if (is.null(xx))
        return (NULL)
    
    if (debug)
        cat("\tEnter inla.interpret.vector.list ... length(xx) = ", length(xx), "\n")

    len = length(xx)
    ind = ii = 1L
    np = as.integer(xx[2L])
    m = length(xx) %/% (2L*np + 2L)
    stopifnot(length(xx) %% (2L*np + 2L) == 0L)
    
    res = lapply(1L:m,
            function(mm, xx, np) {
                i = seq.int((mm-1L)*(2L*np+2L)+3L, len=np, by=2L)
                return (cbind(xx[i],xx[i+1L]))
            }, xx, np)
    
    if (debug)
        cat("\tLeave inla.interpret.vector.list ...\n")

    return(res)
}

`inla.read.binary.file` = function(file)
{
    if (file.exists(file)) {
        fp = file(file,"rb")
        len = file.info(file)$size
        xx = readBin(fp, numeric(), n=len)
        close(fp)
    } else {
        warning(paste("File does not exits: ", file))
        xx = NULL
    } 

    return (xx)
}
