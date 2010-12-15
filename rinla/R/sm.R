### Some utilities for sparse matrices using the `Matrix' library

`inla.sparse.check` = function(C, dims)
{
    msg = paste("Call:", deparse(match.call()))
    if (!is.list(C) && !is.data.frame(C))
        stop(paste(msg, "\n",
                   "INLA-format for sparse-matrix is now `C = list(i=.., j=.., values=...)'. Please fix."))
    if (!(!is.null(C$i) && !is.null(C$j) && !is.null(C$values)))
        stop(paste(msg, "\n",
                   "INLA-format for sparse-matrix is now `C = list(i=.., j=.., values=...)'. Please fix."))
    if ((length(C$i) != length(C$j)) || (length(C$i) != length(C$values)))
        stop(paste(msg, "\n",
                   paste("Length of C$i, C$j and C$values differ:", length(C$i), length(C$j), length(C$values))))

    if (missing(dims)) {
        ni = max(C$i)
        nj = max(C$j)
    } else {
        if (length(dims) == 1) {
            ni = nj = dims
        } else {
            ni = dims[1]
            nj = dims[2]
        }
    }
    if (!(all(C$i <= ni) && all(C$i >= 1)))
        stop(paste(msg, "\n", "C$i is not a subset of the numbers 1, 2, ...,", ni))
    
    if (!(all(C$j <= nj) && all(C$j >= 1)))
        stop(paste(msg, "\n", "C$j is not a subset of the numbers 1, 2, ...,", nj))
    
    if ((length(C$i) != length(C$values)) || (length(C$j) != length(C$i))) {
        stop(paste(msg, "\n",
                   "Wrong dimensions:",
                   ", length(C$i)", length(C$i),
                   ", length(C$j)", length(C$j),
                   ", length(C$values)", length(C$values)))
    }
    
    return (TRUE)
}

`inla.as.dgTMatrix` = function(A)
{
    ## convert into dgTMatrix format of Matrix. Argument A is any matrix
    if (is(A, "dgTMatrix")) {
        return (A)
    } else {
        return (as(A, "dgTMatrix"))
    }
}

`inla.sparse2dgTMatrix` = function(A, dims)
{
    ## convert a inla-sparse-matrix into a dgTMatrix-format of library
    ## 'Matrix'

    require(Matrix)
    if (!missing(dims) && length(dims) == 1) {
        dims = rep(dims, 2)
    }
    B = inla.as.dgTMatrix( sparseMatrix(i = A$i, j = A$j, x = A$values, dims = dims) )
    return (B)
}

`inla.sparse.get` = function(A, row, col)
{
    ## extract a list of the a spesific row or col of a sparseMatrix
    ## (preferable of class dgTMatrix) A. the list returned is of type
    ## list(i=..., j=..., values=...)

    ## need this library
    require(Matrix)

    if (missing(row) && missing(col))
        return (list(i=numeric(0), j=numeric(0), values=numeric(0)))
    
    if (!missing(row) && !missing(col))
        stop("Only one of 'row' and 'col' can be specified.")

    ## need this particular format. I think this is faster for
    ## repeated use, as we do not need to duplicate the A matrix in
    ## memory. but perhaps I'm wrong...
    if (!is(A, "dgTMatrix"))
        A = inla.as.dgTMatrix(A)

    if (!missing(row)) {
        stopifnot(row >=1 && row <= A@Dim[1])
        stopifnot(length(row) == 1)

        idx = which(A@i == row-1)
        return (list(i = row, j = A@j[idx]+1, values = A@x[idx]))
    } else {
        stopifnot(col >=1 && col <= A@Dim[2])
        stopifnot(length(col) == 1)

        idx = which(A@j == col-1)
        return (list(i = A@i[idx]+1, j = col, values = A@x[idx]))
    }
}

`inla.sparse2matrix` = function(C, symmetric = TRUE)
{
    ## convert a sparse-matrix into a matrix.
    
    inla.sparse.check(C)
    stopifnot(length(C$i) == length(C$j))
    stopifnot(length(C$i) == length(C$values))

    n = max(c(C$i, C$j))
    A = matrix(0,n,n)
    for(k in 1:length(C$i)) {
        ## this could be vectorised with more work...
        i = C$i[k]
        j = C$j[k]
        v = C$values[k]
        A[i,j] = v
        if (symmetric)
            A[j,i] = v
    }
    return (A)
}

`inla.matrix2sparse` = function(Q, symmetric = TRUE)
{
    ## convert a possibly symmetric martix matrix into the form:
    ## list(i=,j=,values=). Q is either a matrix or a dgTMatrix

    require(Matrix)
    
    if (is(Q, "dgTMatrix")) {

        ii = Q@i + 1L
        jj = Q@j + 1L
        values = Q@x

        if (symmetric) {
            idx = (ii >= jj)
            ii = ii[idx]
            jj = jj[idx]
            values = values[idx]
        }

    } else {

        stopifnot(is.matrix(Q))

        n = dim(Q)
        if (n[1] != n[2])
            stop(paste("Matrix must be a square matrix, dim(Q) =", dim(Q)))
    
        n = n[1]
        ii = c()
        jj = c()
        values = c()

        for(i in 1:n) {
            if (symmetric){
                idx = which(Q[i, i:n] != 0)
                offset = i-1
            } else {
                idx = which(Q[i, 1:n] != 0)
                offset = 0
            }
            if (length(idx) > 0) {
                idx = idx + offset
                ii = c(ii, rep(i, length(idx)))
                jj = c(jj, idx)
                values = c(values, Q[i, idx])
            }
        }
    }

    return (list(i=ii, j=jj, values = values))
}

`inla.sparse2file` = function(C, filename = NULL, c.indexing = FALSE, binary = FALSE)
{
    if (c.indexing)
        off = 1L
    else
        off = 0L
    
    if (is.null(filename)) {
        filename = tempfile()
    }
    if (!binary) {
        opt = options()
        options(digits=16)
        write(t(cbind(as.integer(C$i-off), as.integer(C$j-off), C$values)), ncolumns=3, file = filename)
        options(opt)
    } else {
        Q  =  cbind(as.integer(C$i), as.integer(C$j), C$values) ## 1-based indexing
        inla.write.fmesher.file(Q, filename = filename)
    }

    return (filename)
}
