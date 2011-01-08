### Some utilities for sparse matrices using the `Matrix' library

`inla.sparse.check` = function(A, must.be.squared = TRUE)
{
    if (is.character(A)) {
        if (!file.exists(A)) {
            stop(paste("File not found:", A))
        }
    }
    
    if (is.list(A))
        stop("Define matrix using Matrix::sparseMatrix() instead!!! The list(i=,j=,values=)-format is obsolete!")

    if (!is(A, "dgTMatrix")) {
        A = inla.as.dgTMatrix(A)
    }

    if (must.be.squared) {
        if (dim(A)[1] != dim(A)[2]) {
            stop(paste(c("Matrix is not a square matrix:", dim(A))))
        }
    }

    return (A)
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

    warning("THIS FUNCTION IS OBSOLETE!")
    
    if (!missing(dims) && length(dims) == 1) {
        dims = rep(dims, 2)
    }
    B = inla.as.dgTMatrix( sparseMatrix(i = A$i, j = A$j, x = A$values, dims = dims) )
    return (B)
}

`inla.sparse.get` = function(A, row, col)
{
    ## extract a list of the a spesific row or col of a dgTMatrix
    ## A. the list returned is of type list(i=..., j=..., values=...)

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

`inla.sparse2matrix` = function(A, symmetric = TRUE)
{
    warning("THIS FUNCTION IS OBSOLETE!!!")

    ## convert a sparse-matrix into a matrix.
    
    inla.sparse.check(A)
    stopifnot(length(A$i) == length(A$j))
    stopifnot(length(A$i) == length(A$values))

    n = max(c(A$i, A$j))
    A = matrix(0,n,n)
    for(k in 1:length(A$i)) {
        ## this could be vectorised with more work...
        i = A$i[k]
        j = A$j[k]
        v = A$values[k]
        A[i,j] = v
        if (symmetric)
            A[j,i] = v
    }
    return (A)
}

`inla.matrix2sparse` = function(Q, symmetric = TRUE)
{
    warning("THIS FUNCTION IS OBSOLETE!!!")

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

`inla.sparse2file` = function(A, filename = NULL, c.indexing = FALSE, binary = FALSE, symmetric = FALSE)
{
    if (is.null(filename)) {
        filename = tempfile()
    }

    if (!is(A, "dgTMatrix"))
        A = inla.as.dgTMatrix(A)

    if (!binary) {
        if (c.indexing)
            off = 0L
        else
            off = 1L
    
        if (symmetric) {
            idx = (A@i >= A@j)
        } else {
            idx = 1:length(A@i)
        }

        opt = options()
        options(digits=16)
        write(t(cbind(as.integer(A@i[idx]+off), as.integer(A@j[idx]+off), A@x[idx])), ncolumns=3, file = filename)
        options(opt)
    } else {
        off = 1L
        Q  =  cbind(as.integer(A@i+off), as.integer(A@j+off), A@x)  ## 1-based indexing
        inla.write.fmesher.file(Q, filename = filename)
    }

    return (filename)
}
