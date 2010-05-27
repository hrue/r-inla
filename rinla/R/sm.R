### Some utilities for sparse matrices using the `Matrix' library

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
    if (!is(A, "dgTMatrix")) {
        A = as(A, "dgTMatrix")
    }

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
