## Nothing to export

## obsolete functions goes here.

`inla.obsolete` = function(old, new) {
    stop(inla.paste(c("INLA-function `", as.character(old)[1], "' is obsolete. Use function `",
                      as.character(new), "' instead."), sep=""))
}

`inla.Cmatrix2matrix` = function(C, symmetric = TRUE)
{
    inla.obsolete(match.call(), "inla.sparse2matrix")
    return (invisible())
}

`inla.matrix2Cmatrix` = function(Q, symmetric = TRUE)
{
    inla.obsolete(match.call(), "inla.matrix2sparse")
    return (invisible())
}
    
`inla.Cmatrix2file` = function(Cmatrix, filename = NULL, c.indexing = FALSE)
{
    inla.obsolete(match.call(), "inla.write.fmesher.file")
    return (invisible())
}

`inla.sparse2dgTMatrix` = function(A, dims)
{
    ## convert a inla-sparse-matrix into a dgTMatrix-format of library
    ## 'Matrix'
    stop("THIS FUNCTION IS OBSOLETE: inla.sparse2dgTMatrix")
    if (!missing(dims) && length(dims) == 1) {
        dims = rep(dims, 2)
    }
    B = inla.as.dgTMatrix( sparseMatrix(i = A$i, j = A$j, x = A$values, dims = dims) )
    return (B)
}

`inla.sparse2matrix` = function(A, symmetric = TRUE)
{
    stop("THIS FUNCTION IS OBSOLETE: inla.sparse2matrix")

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
    stop("THIS FUNCTION IS OBSOLETE: inla.matrix2sparse")

    ## convert a possibly symmetric martix matrix into the form:
    ## list(i=, j=, values=). Q is either a matrix or a dgTMatrix

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
    ## FIXME: switch to the fmesher-binary format ONLY later, but we
    ## need a TAG to identify such files uniquely.

    stop("THIS FUNCTION IS OBSOLETE: inla.sparse2file")

    if (is.null(filename)) {
        filename = tempfile()
    }

    A = inla.as.dgTMatrix(A)
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
        ## the matrix will be correct. 

        off = 1L
        Q = list(i = as.integer(c(0L + off, dims[1]-1L + off, A@i + off)),
                j = as.integer(c(0L + off, dims[2]-1L + off, A@j + off)),
                values = c(0.0, 0.0, A@x))
        inla.write.fmesher.file(Q, filename = filename)
    }

    return (filename)
}

inla.hyperpar.sampler = function(...) {
    warning("New name is 'inla.hyperpar.sample'")
    return (inla.hyperpar.sample(...))
}
