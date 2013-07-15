## Nothing to Export

### Binary I/O and interpretation of vectors

`inla.interpret.vector` = function(xx, debug=FALSE)
{
    if (is.null(xx))
        return (NULL)
    
    if (debug)
        cat("\tEnter inla.interpret.vector ... length(xx) = ", length(xx), "\n")

    len = length(xx)
    np = as.integer(xx[2L])
    stopifnot((len %% (2L*np+2L)) == 0L)
    res.ncol = (2L * len) %/% (2L*np+2L)
    res = matrix(0.0, np, res.ncol)
    m = 1L
    ii = 1L
    np2 = 2L*np + 2L
    for (j in 1L:(res.ncol %/% 2L)) {
        idx = seq.int(ii + 2L, by=2L, length.out=np)
        res[, c(m, m+1L)] = xx[c(idx, idx+1L)]
        m = m + 2L
        ii = ii + np2
    }
    stopifnot(len == ii-1L)
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
    res = numeric()
    np = as.integer(xx[2L])
    stopifnot((len %% (2L*np + 2L)) == 0L)
    m = (3L * len) %/% (2L*np+2L)
    res.ncol = m
    res = matrix(0.0, np, res.ncol)

    m = 1L
    ii = 1L
    np2 = 2L * np + 2L
    for(j in 1L:(m %/% 3L)) {
        idx = as.integer(xx[ii])
        jj = seq.int(ii + 2L, by=2L, length.out = np)
        res[, m] = idx
        res[, c(m+1L, m+2L)] = xx[c(jj, jj+1L)]
        m = m + 3L
        ii = ii + np2
    }
    stopifnot(len == ii-1L)
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

    ## this gives approximately the correct length, most likely, an
    ## upper bound.
    m.approx = length(xx) %/% 54L
    nps = integer(m.approx)  ## all lengths
    inds = integer(m.approx) ## all offsets where the density starts

    ## need to get the lengths and the offsets as they vary.
    count = 1L
    ind = 2L
    while (ind < length(xx)) {
        np = as.integer(xx[ind])
        inds[count]= ind+1L
        nps[count] = np
        count = count + 1L
        ind = ind + 2L*np + 2L
    }

    count = count - 1L
    res = lapply(1L:count,
            function(m, xx, inds, nps) {
                i = seq.int(inds[m], length.out = nps[m], by=2L)
                return (cbind(xx[i], xx[i + 1L]))
            }, xx, inds, nps)
    
    if (debug)
        cat("\tLeave inla.interpret.vector.list ...\n")

    return(res)
}

`inla.read.binary.file` = function(file, silent = TRUE)
{
    if (file.exists(file)) {
        fp = file(file,"rb")
        len = file.info(file)$size
        xx = readBin(fp, numeric(), n=len)
        close(fp)
    } else {
        if (!silent) {
            warning(paste("File does not exits: ", file))
        }
        xx = NULL
    } 

    return (xx)
}
