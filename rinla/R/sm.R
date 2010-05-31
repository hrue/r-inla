### Some utilities for sparse matrices using the `Matrix' library

`inla.sm.check` = function(C, dims)
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
    return (as(A, "dgTMatrix"))
}

`inla.sm2dgTMatrix` = function(A, dims)
{
    ## convert a inla-sparse-matrix into a dgTMatrix-format of library
    ## 'Matrix'

    require(Matrix)
    B = inla.as.dgTMatrix( sparseMatrix(i = A$i, j = A$j, x = A$values, dims = dims) )
    return (B)
}

`inla.sm.get` = function(A, irow, icol)
{
    ## extract a list of the a spesific row or col of a sparseMatrix
    ## (preferable of class dgTMatrix) A. the list returned is of type
    ## list(i=..., j=..., values=...)

    ## need this library
    require(Matrix)

    if (missing(irow) && missing(icol))
        return (list(i=numeric(0), j=numeric(0), values=numeric(0)))
    
    if (!missing(irow) && !missing(icol))
        stop("Only one of 'irow' and 'icol' can be specified.")

    ## need this particular format
    A = inla.as.dgTMatrix(A)

    if (!missing(irow)) {
        stopifnot(irow >=1 && irow <= A@Dim[1])
        stopifnot(length(irow) == 1)

        idx = which(A@i == irow-1)
        return (list(i = irow, j = A@j[idx]+1, values = A@x[idx]))
    } else {
        stopifnot(icol >=1 && icol <= A@Dim[2])
        stopifnot(length(icol) == 1)

        idx = which(A@j == icol-1)
        return (list(i = A@i[idx]+1, j = icol, values = A@x[idx]))
    }
}
