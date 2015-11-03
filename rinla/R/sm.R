## Export: inla.as.sparse inla.as.dgTMatrix

##!\name{inla.as.sparse}
##!\alias{inla.as.sparse}
##!\alias{inla.as.dgTMatrix}
##!
##!\title{Convert a matrix or sparse matrix into the sparse formate used by INLA}
##!
##!\description{Convert a matrix or sparse matrix into the sparse format used by INLA (dgTMatrix)}
##!
##!\usage{
##! inla.as.sparse(...)
##! inla.as.dgTMatrix(A, unique = TRUE)
##!}
##!
##!\arguments{
##!  \item{...}{The arguments. The matrix or sparse matrix,  and the additonal arguments}
##!  \item{A}{The matrix}
##!  \item{unique}{If the internal representation needs to be unique or can have duplicated entries.
##!                Do not use this option unless you know what you are doing.}
##!}
##!
##!\value{%%
##!  \code{inla.as.sparse} and \code{inla.as.dgTMatrix} is the same function.
##!  The returned value is a sparse matrix in the sparse-format used by INLA
##!}
##!%%
##!
##!\author{Havard Rue \email{hrue@math.ntnu.no}}
##!
##!\examples{
##! A = matrix(1:9, 3, 3)
##! inla.as.sparse(A)
##!}

`inla.as.sparse` = function(...)
{
    return (inla.as.dgTMatrix(...))
}



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

    if (must.be.squared) {
        if (dim(A)[1] != dim(A)[2]) {
            stop(paste(c("Matrix is not a square matrix:", dim(A)[1], "x", dim(A)[2])))
        }
    }

    return (inla.as.dgTMatrix(A))
}

`inla.as.dgTMatrix` = function(A, unique = TRUE)
{
    ## convert into dgTMatrix format of Matrix. Argument A is any
    ## matrix.  make sure the representation is unique if the UNIQUE
    ## flag it TRUE. (ie no double triplets etc)

    if (unique) {
        ## convert through the 'dgCMatrix'-class to make it unique;
        ## (there is no method or default for coercing "dtTMatrix" to
        ## "dgCMatrix"). Yes: its 'CsparseMatrix'.
        return (as(as(as(A, "CsparseMatrix"), "dgCMatrix"), "dgTMatrix"))
    } else {
        if (is(A, "dgTMatrix")) {
            return (A)
        } else {
            ## Convert via virtual class TsparseMatrix;
            ## this allows more general conversions than direct conversion.
            return (as(as(A, "TsparseMatrix"), "dgTMatrix"))
        }
    }
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
