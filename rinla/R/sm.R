#' Convert a matrix or sparse matrix into the sparse formate used by INLA
#' 
#' Convert a matrix or sparse matrix into the sparse format used by INLA
#' (dgTMatrix)
#' 
#' 
#' @aliases inla.as.sparse inla.as.dgTMatrix
#' @param ... The arguments. The matrix or sparse matrix, and the additonal
#' arguments
#' @param A The matrix
#' @param unique Logical. If `TRUE`, then ensure that the internal
#' representation is unique and there are no duplicated entries.  (Do not
#' change this unless you know what you are doing.)
#' @param na.rm Replace `NA`'s in the matrix with zeros.
#' @param zeros.rm Remove zeros in the matrix.
#' @returns `inla.as.sparse` and `inla.as.dgTMatrix` is the same
#' function.  The returned value is a sparse matrix in the
#' `dgTMatrix`-format.
#' @author Havard Rue \email{hrue@@r-inla.org}
#' @examples
#' 
#'  A = matrix(1:9, 3, 3)
#'  inla.as.sparse(A)
#' 
#' @rdname sm
#' @export inla.as.sparse
`inla.as.sparse` <- function(A, unique = TRUE, na.rm = FALSE, zeros.rm = FALSE) {
    ## convert into dgTMatrix format of Matrix. Argument A is any
    ## matrix.  make sure the representation is unique if the UNIQUE
    ## flag it TRUE. (ie no double triplets etc).
    ## if na.rm then replace NA's with zeros.
    ## if zeros.rm then remove zeros.
    
    if (!inherits(A, "Matrix")) {
        A <- as(A, "Matrix")
    }

    if (unique) {
        ## Convert via double, general, Tsparse Matrix classes, following
        ## the recommendations of the Matrix package developer vignette.
        ## Convert to 'CsparseMatrix'-class first to ensure the T-representation
        ## is unique, i.e. has no duplicate (i,j) pairs.
        A <- as(as(as(as(A,
                         "dMatrix"),
                      "generalMatrix"),
                   "CsparseMatrix"),
                "TsparseMatrix")
    } else {
        if (!is(A, "dgTMatrix")) {
            ## Convert via double, general, Tsparse Matrix classes, following
            ## the recommendations of the Matrix package developer vignette.
            A <- as(as(as(A,
                          "dMatrix"),
                       "generalMatrix"),
                    "TsparseMatrix")
        }
    }

    if (na.rm) {
        idx.na <- is.na(A@x)
        if (any(idx.na)) {
            A@x[idx.na] <- 0.0
        }
    }

    if (zeros.rm) {
        x.zero <- (A@x == 0.0)
        if (any(x.zero)) {
            idx.zero <- which(x.zero)
            A@x <- A@x[-idx.zero]
            A@i <- A@i[-idx.zero]
            A@j <- A@j[-idx.zero]
        }
    }

    return(A)
}

#' @rdname sm
#' @export
`inla.as.dgTMatrix` <- function(...) {
    return(inla.as.sparse(...))
}

### Some utilities for sparse matrices using the `Matrix' library

`inla.sparse.dim` <- function(A) {
    ## return the dimension of the matrix A
    A <- inla.sparse.check(A, must.be.squared = FALSE)

    if (is.character(A)) {
        Am <- read.table(A, col.names = c("i", "j", "Aij"))
        if (min(Am$i) == 0 || min(Am$j) == 0) {
            cindex <- 1
        } else {
            cindex <- 0
        }
        return(c(max(Am$i) + cindex, max(Am$j) + cindex))
    } else {
        return(dim(A))
    }
}


`inla.sparse.check` <- function(A, must.be.squared = TRUE) {
    ## just check if matrix A exists and return either a filename or a
    ## dgTMatrix.
    if (is.character(A)) {
        if (!file.exists(A)) {
            stop(paste("File not found:", A))
        }
        return(A)
    }

    if (is.list(A)) {
          stop("Define matrix using Matrix::sparseMatrix() instead!!! The list(i=, j=, values=)-format is obsolete!")
      }

    if (must.be.squared) {
        if (dim(A)[1] != dim(A)[2]) {
            stop(paste(c("Matrix is not a square matrix:", dim(A)[1], "x", dim(A)[2])))
        }
    }

    return(inla.as.dgTMatrix(A))
}

`inla.sparse.get` <- function(A, row, col) {
    ## extract a list of the a specific row or col of a dgTMatrix
    ## A. the list returned is of type list(i=..., j=..., values=...)

    if (missing(row) && missing(col)) {
          return(list(i = numeric(0), j = numeric(0), values = numeric(0)))
      }

    if (!missing(row) && !missing(col)) {
          stop("Only one of 'row' and 'col' can be specified.")
      }

    ## need this particular format. I think this is faster for
    ## repeated use, as we do not need to duplicate the A matrix in
    ## memory. but perhaps I'm wrong...
    if (!is(A, "dgTMatrix")) {
        ## This can be slow, so its better to stop and say that the
        ## matrix has to be converted upfront.
        stop("Matrix is not of type 'dgTMatrix'; please convert it with inla.as.dgTMatrix().")
        A <- inla.as.dgTMatrix(A)
    }

    if (!missing(row)) {
        stopifnot(row >= 1 && row <= A@Dim[1])
        stopifnot(length(row) == 1)

        idx <- which(A@i == row - 1)
        return(list(i = row, j = A@j[idx] + 1, values = A@x[idx]))
    } else {
        stopifnot(col >= 1 && col <= A@Dim[2])
        stopifnot(length(col) == 1)

        idx <- which(A@j == col - 1)
        return(list(i = A@i[idx] + 1, j = col, values = A@x[idx]))
    }
}

inla.sm.write <- function(A, filename = "SparseMatrix.dat") {
    stopifnot(!missing(A))
    A <- inla.as.sparse(A)

    x.int <- as.integer(c(dim(A)[1], dim(A)[2], length(A@i), A@i + 1, A@j + 1))
    fp <- file(filename, "wb")
    writeBin(x.int, fp)
    writeBin(A@x, fp)
    close(fp)

    return(filename)
}

inla.sm.read <- function(filename = "SparseMatrix.dat") {
    fp <- file(filename, "rb")
    nn <- readBin(fp, integer(), n = 3L)
    nx <- nn[3]
    ii <- readBin(fp, integer(), n = nx)
    jj <- readBin(fp, integer(), n = nx)
    xx <- readBin(fp, double(), n = nx)
    close(fp)

    M <- inla.as.sparse(sparseMatrix(
        i = ii - 1L,
        j = jj - 1L,
        x = xx,
        dims = nn[1:2],
        index1 = FALSE,
        repr = "T"
    ))
    return(M)
}
