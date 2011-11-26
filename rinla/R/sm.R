### Some utilities for sparse matrices using the `Matrix' library

`inla.sparse.dim` = function(A) {
    ## return the dimension of the matrix A
    A = inla.sparse.check(A, must.be.squared = FALSE)

    if (is.character(A)) {
        Am = read.table(A, col.names = c("i", "j", "Aij"))
        if (min(Am$i) == 0 || min(Am$j) == 0) {
            cindex = 1
        } else {
            cindex = 0
        }
        return (c(max(Am$i) + cindex, max(Am$j) + cindex))
    } else {
        return (dim(A))
    }
}


`inla.sparse.check` = function(A, must.be.squared = TRUE)
{
    ## just check if matrix A exists and return either a filename or a
    ## dgTMatrix.
    if (is.character(A)) {
        if (!file.exists(A)) {
            stop(paste("File not found:", A))
        }
        return (A)
    }

    if (is.list(A))
        stop("Define matrix using Matrix::sparseMatrix() instead!!! The list(i=, j=, values=)-format is obsolete!")

    if (!is(A, "dgTMatrix")) {
        A = inla.as.dgTMatrix(A)
    }

    if (must.be.squared) {
        if (dim(A)[1] != dim(A)[2]) {
            stop(paste(c("Matrix is not a square matrix:", dim(A)[1], "x", dim(A)[2])))
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
        ## Convert via virtual class TsparseMatrix;
        ## this allows more general conversions than direct conversion.
        return (as(as(A, "TsparseMatrix"), "dgTMatrix"))
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
    ## extract a list of the a specific row or col of a dgTMatrix
    ## A. the list returned is of type list(i=..., j=..., values=...)

    if (missing(row) && missing(col))
        return (list(i=numeric(0), j=numeric(0), values=numeric(0)))

    if (!missing(row) && !missing(col))
        stop("Only one of 'row' and 'col' can be specified.")

    ## need this particular format. I think this is faster for
    ## repeated use, as we do not need to duplicate the A matrix in
    ## memory. but perhaps I'm wrong...
    if (!is(A, "dgTMatrix")) {
        ## This can be slow, so its better to stop and say that the
        ## matrix has to be converted upfront.
        stop("Matrix is not of type 'dgTMatrix'; please convert it with inla.as.dgTMatrix().")
        A = inla.as.dgTMatrix(A)
    }

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
    A = matrix(0, n, n)
    for(k in 1:length(A$i)) {
        ## this could be vectorised with more work...
        i = A$i[k]
        j = A$j[k]
        v = A$values[k]
        A[i, j] = v
        if (symmetric)
            A[j, i] = v
    }
    return (A)
}

`inla.matrix2sparse` = function(Q, symmetric = TRUE)
{
    warning("THIS FUNCTION IS OBSOLETE!!!")

    ## convert a possibly symmetric martix matrix into the form:
    ## list(i=, j=, values=). Q is either a matrix or a dgTMatrix

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

`inla.sparse2file` = function(A, filename = NULL, c.indexing = FALSE,
        binary = inla.getOption("internal.binary.mode"), symmetric = FALSE)
{
    if (is.null(filename)) {
        filename = tempfile()
    }

    if (!is(A, "dgTMatrix")) {
        A = inla.as.dgTMatrix(A)
    }
    dims = dim(A)

    if (!binary) {
        if (c.indexing) {
            off = 0L
        } else {
            off = 1L
        }

        if (symmetric) {
            idx = which(A@i >= A@j)
        } else {
            idx = 1:length(A@i)
        }

        opt = options()
        options(digits=16)

        ## make sure that at the elements (1, 1) and (n, m) are
        ## written out (just set them to zero), so the dimension of
        ## the matrix will be correct. In case they are there in any
        ## case, the entries will be added up.
        ##
        ## FIXME: switch to the fmesher-binary format later, but we
        ## need a TAG to identify such files uniquely.
        ##

        write(t(cbind(as.integer(c(0L + off, dims[1]-1L + off, A@i[idx] + off)),
                      as.integer(c(0L + off, dims[2]-1L + off, A@j[idx] + off)),
                      c(0.0, 0.0, A@x[idx]))),
              ncolumns=3, file = filename)
        options(opt)

    } else {

        ## make sure that at the elements (1, 1) and (n, m) are
        ## written out (just set them to zero), so the dimension of
        ## the matrix will be correct. In case they are there in any
        ## case, the entries will be added up.
        ##
        ## FIXME: switch to the fmesher-binary format ONLY later, but
        ## we need a TAG to identify such files uniquely.
        ##

        off = 1L
        Q = list(i = as.integer(c(0L + off, dims[1]-1L + off, A@i + off)),
                j = as.integer(c(0L + off, dims[2]-1L + off, A@j + off)),
                values = c(0.0, 0.0, A@x))
        inla.write.fmesher.file(Q, filename = filename)
    }

    return (filename)
}
