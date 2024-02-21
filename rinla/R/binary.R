## Nothing to Export

### Binary I/O and interpretation of vectors

`inla.interpret.vector` <- function(xx, debug = FALSE)
{
    if (is.null(xx)) {
        return(NULL)
    }

    if (debug) {
        cat("\tEnter inla.interpret.vector ... length(xx) = ", length(xx), "\n")
    }

    len <- length(xx)
    np <- as.integer(xx[2L])
    stopifnot((len %% (2L * np + 2L)) == 0L)
    res.ncol <- (2L * len) %/% (2L * np + 2L)
    res.ncol2 <- res.ncol %/% 2L
    m <- 1L
    ii <- 3L
    np2 <- 2L * np + 2L
    res <- matrix(NA, np, res.ncol)

    if (FALSE) {
        ## slow: this loops over columns
        for (j in 1L:(res.ncol %/% 2L)) {
            idx <- seq.int(ii, by = 2L, length.out = np)
            res[, m] <- xx[idx]
            res[, m + 1L] <- xx[idx + 1L]
            m <- m + 2L
            ii <- ii + np2
        }
        stopifnot(res.ncol == m - 1L) ## since we start from 1L
        stopifnot(len == ii - 3L)
    }
    if (FALSE) {
        ## better: this loops over rows
        idx <- seq.int(from = 3L, by = np2, length.out = res.ncol2)
        idx2 <- seq.int(1L, by = 2L, length.out = res.ncol2)
        for(i in 1L:np) {
            res[i, idx2] <- xx[idx]
            res[i, idx2 + 1L] <- xx[idx + 1L]
            idx <- idx + 2L
        }
    }
    if (TRUE) {
        ## best 
        idx <- seq.int(from = 3L, by = np2, length.out = res.ncol2)
        idx2 <- seq.int(1L, by = 2L, length.out = res.ncol2)
        for(i in 1L:np) {
            res[i, idx2] <- xx[idx + 2L * (i-1L)]
            res[i, idx2 + 1L] <- xx[idx + 2L * (i-1L) + 1L]
        }
    }

    if (debug) {
        cat("\tLeave inla.interpret.vector ...\n")
    }

    return(res)
}

`inla.interpret.vector.list` <- function(xx, debug = FALSE)
{
    if (is.null(xx)) {
        return(NULL)
    }

    if (debug) {
        cat("\tEnter inla.interpret.vector.list ... length(xx) = ", length(xx), "\n")
    }

    ## this gives approximately the correct length, most likely, an
    ## upper bound.
    m.approx <- length(xx) %/% 54L
    nps <- integer(m.approx) ## all lengths
    inds <- integer(m.approx) ## all offsets where the density starts

    ## need to get the lengths and the offsets as they vary.
    count <- 0L
    ind <- 2L
    while (ind < length(xx)) {
        count <- count + 1L
        np <- as.integer(xx[ind])
        inds[count] <- ind + 1L
        nps[count] <- np
        ind <- ind + 2L * np + 2L
    }

    fun <- function(m, xx, inds, nps) {
        ii <- seq.int(inds[m], length.out = nps[m], by = 2L)
        return(cbind(xx[ii], xx[ii + 1L]))
    }
    res <- lapply(1L:count, fun, xx, inds, nps)

    if (debug) {
        cat("\tLeave inla.interpret.vector.list ...\n")
    }

    return(res)
}

`inla.read.binary.file` <- function(file, silent = TRUE)
{
    if (file.exists(file)) {
        fp <- file(file, "rb")
        len <- file.info(file)$size
        xx <- readBin(fp, numeric(), n = len)
        xx[is.nan(xx)] <- NA
        close(fp)
    } else {
        if (!silent) {
            warning(paste("File does not exits: ", file))
        }
        xx <- NULL
    }

    return(xx)
}
