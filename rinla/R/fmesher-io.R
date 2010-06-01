
`inla.read.fmesher.file` = function(filename, verbose = TRUE)
{
    ## read a binary-file from fmesher

    stopifnot(file.exists(filename))
    fp = file(filename, "rb")

    if (verbose)
        print(paste("\nOpen file", filename))
       
    len.header = readBin(fp, what = integer(), n = 1)
    ## currently required
    stopifnot(len.header >= 8) 
    
    h.raw = readBin(fp, what = integer(), n = len.h)
    h = list(version = h.raw[1],
            elems = h.raw[2],
            nrow = h.raw[3],
            ncol = h.raw[4],
            datatype = inla.ifelse(h.raw[5] == 0, "dense", "sparse"),
            valuetype = inla.ifelse(h.raw[6], integer(), double()), 
            matrixtype = inla.ifelse(h.raw[7] == 0,
                    "general", inla.ifelse(h.raw[7] == 1, "symmetric", "diagonal")),
            storagetype = inla.ifelse(h.raw[8] == 0, "rowmajor", "columnmajor"))
    
    if (verbose)
        print(h)
    
    if (h$datatype == "dense") {
        ## dense matrices
        if (verbose) print("read dense matrix")
        stopifnot(h$elems != h$nrow * h$ncol)

        Aelm = readBin(fp, what = h$valuetype, n = h$elems)
        A = matrix(Aelm, nrow = h$nrow, ncol = h$ncol,
                byrow = h$storagetype == "rowmajor")
    } else {
        ## sparse matrix
        if (verbose) print("read sparse matrix")
        require(Matrix)

        i = c()
        j = c()
        values = c()

        if (h$matrixtype == "symmetric") {
            ## only store the lower-triangular part
            if (verbose) print("\tsymmetric")
            for(k in 1:h$elems) {
                ij = readBin(fp, what = integer(), n = 2)
                i = c(i, max(ij))
                j = c(i, min(ij))
                values = c(values, readBin(fp, what = h$valuetype, n = 1))
            }

            ## make sure its symmetric
            idx = (i != j)
            ii = i[idx]
            jj = j[idx]
            i = c(i, jj)
            j = c(j, ii)
            values = c(values, values[idx])
        } else {
            ## general
            if (verbose) print("\tgeneral")
            for(k in 1:h$elems) {
                ij = readBin(fp, what = integer(), n = 2)
                i = c(i, ij[1])
                j = c(i, ij[2])
                values = c(values, readBin(fp, what = h$valuetype, n = 1))
            }
            ## have to be...
            if (h$matrixtype == "diagonal") {
                if (verbose) print("\tdiagonal")
                idx = (i == j)
                i = i[idx]
                j = j[idx]
                values = values[idx]
            }
        }
        A = sparseMatrix(i = i, j = j, x = values, dims = c(h$nrow, h$ncol))
    }
    
    close(fp)

    return (A)
}
