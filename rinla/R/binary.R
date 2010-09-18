### Binary I/O and interpretation of vectors

## keep old version in the source for a while...

`inla.interpret.vector` = function(xx)
{
    if (is.null(xx))
        return (NULL)
    
    ## This new version first determine the size of the matrix, then
    ## just fills in the matrix. This is much more efficient than
    ## doing this dynamically.

    len = length(xx)
    ii = 1L
    res = numeric()
    by = 2L
    m = 0L
    np = NA
    ## first do the counting
    while(ii <= len) {
        np = as.integer(xx[ii + 1L])
        if (np) {
            m = m + 2L
        }
        ii = ii + 2L*np + 2L
    }

    ## then allocate the storage
    res.ncol = m
    res = matrix(NA, np, res.ncol)

    m = 1L
    ii = 1L
    while(ii <= len) {
        np = as.integer(xx[ii+1L])
        if (np > 0) {
            start = ii + 2L
            end = ii + 2L*np + 1L
            idx = seq.int(start, end, by=by)
            res[, m] = xx[idx]
            res[, m+1] = xx[idx+1]
            m = m + 2L
            ii = end + 1L
        } else {
            ii = ii + 2L
        }
    }
    stopifnot(res.ncol == m-1L) ## since we start m at 1L
    return(res)
}

`inla.interpret.vector.idx` = function(xx)
{
    ## same here, count first...

    if (is.null(xx))
        return (NULL)
    
    len = length(xx)
    ii = 1L
    res = numeric()
    by = 2L
    m = 0L
    np = NA
    while(ii <= len) {
        np = as.integer(xx[ii + 1L])
        if (np) {
            m = m + 3L
        }
        ii = ii + 2L*np + 2L
    }
    res.ncol = m
    res = matrix(NA, np, res.ncol)

    m = 1L
    ii = 1L
    while(ii <= len) {
        idx = as.integer(xx[ii])
        np = as.integer(xx[ii+1L])
        if (np > 0) {
            start = ii + 2L
            end = ii + 2L*np + 1L
            jj = seq.int(start, end, by=by)
            res[, m] = idx
            res[, m+1] = xx[jj]
            res[, m+2] = xx[jj+1L]
            m = m + 3L
            ii = end + 1L
        } else {
            ii = ii + 2L;
        }
    }
    stopifnot(res.ncol == m-1L) ## since we start from 1L
    return(res)
}

`inla.interpret.vector.list` = function(xx)
{
    if (is.null(xx))
        return (NULL)
    
    len = length(xx)
    ind = ii = 1L
    res = list()
    by=2L
    while(ii <= len) {
        np = as.integer(xx[ii+1L])
        if (np>0) {
            start = ii + 2L
            end = ii + 2L*np + 1L
            jj = seq.int(start, end, by=by)
            res[[ind]] = cbind(xx[jj], xx[jj+1L])
            ind = ind + 1L
            ii = end + 1L
        } else {
            res[[ind]] = NA
            ii = ii + 2L
        }
    }
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
