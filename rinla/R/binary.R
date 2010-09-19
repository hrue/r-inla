### Binary I/O and interpretation of vectors

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
    m = 0L
    np = NA
    no.missing = FALSE
    ## first do the counting

    np = as.integer(xx[2L])
    if ((length(xx) %% (2L*np+2)) == 0L) {
        ## no missing, all np > 0
        m = 2L * length(xx) %/% (2L*np+2)
        no.missing = TRUE
    } else {
        ## have to check carefully
        while(ii <= len) {
            npp = as.integer(xx[ii + 1L])
            if (npp) {
                np = npp
                m = m + 2L
            }
            ii = ii + 2L*npp + 2L
        }
    }
    ## then allocate the storage
    res.ncol = m
    res = matrix(0, np, res.ncol)

    m = 1L
    ii = 1L
    if (no.missing) {
        ## we can do this a bit more cleanly
        np = as.integer(xx[2L])
        np2 = 2L*np + 2L
        for(j in 1:(res.ncol %/% 2)) {
            idx = seq.int(ii + 2L, by=2L, len=np)
            res[, c(m, m+1L)] = xx[c(idx, idx+1L)]
            m = m + 2L
            ii = ii + np2
        }
    } else {
        while(ii <= len) {
            np = as.integer(xx[ii+1L])
            if (np > 0) {
                start = ii + 2L
                end = ii + 2L*np + 1L
                idx = seq.int(start, end, by=2L)
                res[, c(m,m+1L)] = xx[c(idx,idx+1L)]
                m = m + 2L
                ii = end + 1L
            } else {
                ii = ii + 2L
            }
        }
    }

    stopifnot(length(xx) == ii-1L)
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
    m = 0L
    np = NA
    no.missing = FALSE
    
    np = as.integer(xx[2L])
    if ((length(xx) %% (2L*np + 2L)) == 0L) {
        m = 3L * length(xx) %/% (2L*np+2L)
        no.missing = TRUE
    } else {
        while(ii <= len) {
            np = as.integer(xx[ii + 1L])
            if (np) {
                m = m + 3L
            }
            ii = ii + 2L*np + 2L
        }
    }
    res.ncol = m
    res = matrix(NA, np, res.ncol)

    m = 1L
    ii = 1L
    if (no.missing) {
        ## cleaner code
        np = as.integer(xx[2L])
        np2 = 2L * np + 2L
        for(j in 1:(m %/% 3L)) {
            idx = as.integer(xx[ii])
            jj = seq.int(ii + 2L, by=2L, len = np)
            res[, m] = idx
            res[, c(m+1L, m+2L)] = xx[c(jj, jj+1L)]
            m = m + 3L
            ii = ii + np2
        }
    } else {
        while(ii <= len) {
            idx = as.integer(xx[ii])
            np = as.integer(xx[ii+1L])
            if (np > 0) {
                start = ii + 2L
                end = ii + 2L*np + 1L
                jj = seq.int(start, end, by=2L)
                res[, m] = idx
                res[, c(m+1L, m+2L)] = xx[c(jj, jj+1L)]
                m = m + 3L
                ii = end + 1L
            } else {
                ii = ii + 2L;
            }
        }
    }

    stopifnot(length(xx) == ii-1L)
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
    while(ii <= len) {
        np = as.integer(xx[ii+1L])
        if (np>0) {
            jj = seq.int(ii + 2L, len=np, by=2L)
            res[[ind]] = cbind(xx[jj], xx[jj+1L])
            ii = ii + 2L*np + 2L
        } else {
            res[[ind]] = NA
            ii = ii + 2L
        }
        ind = ind + 1L
    }

    stopifnot(length(xx) == ii-1L)

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
