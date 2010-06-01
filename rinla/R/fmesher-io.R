
`inla.read.fmesher.file` = function(filename, verbose = FALSE, debug = FALSE)
{
    ##
    ## read a binary-file from fmesher in format specified by FL.
    ##

    if (debug)
        verbose=TRUE
    
    ## internal simple checking routine, which also do debug output.
    read.check = function(x, h)
    {
        msg = deparse(match.call())
        if (length(x) != h$elems)
            stop(paste("Reading file", filename, ". Fail to read", h$elems, ", got", length(x), "."))
        if (h$debug)
            print(inla.paste(c(msg, x)))

        return (invisible())
    }
               
    stopifnot(file.exists(filename))
    fp = file(filename, "rb")

    if (verbose)
        print(paste("\nOpen file", filename))
       
    len.h = readBin(fp, what = integer(), n = 1)
    ## currently required
    stopifnot(len.h >= 8)
    if (verbose)
        print(paste("header is", len.h, "integers."))
    
    h.raw = readBin(fp, what = integer(), n = len.h)
    ## negative entries specify something different not yet defined.
    for(i in 1:len.h) {
        if (h.raw[i] < 0)
            stop(paste("Entry ", i, " in the header of file ", filename, " is negative: ", h.raw[i],
                       ". Do not know what to do.", sep=""))
    }

    h = list(filename = filename, verbose = verbose, debug = debug,
            version = h.raw[1],
            elems = h.raw[2],
            nrow = h.raw[3],
            ncol = h.raw[4],
            datatype = inla.ifelse(h.raw[5] == 0, "dense", "sparse"),
            valuetype = inla.ifelse(h.raw[6] == 0, integer(), double()), 
            matrixtype = inla.ifelse(h.raw[7] == 0, "general",
                    inla.ifelse(h.raw[7] == 1, "symmetric", "diagonal")),
            storagetype = inla.ifelse(h.raw[8] == 0, "rowmajor", "columnmajor"))
    
    if (verbose)
        print(h)
    
    if (h$datatype == "dense") {
        ##
        ## dense matrix
        ##
        if (h$matrixtype != "general")
            stop(paste("Read", filename, ". Type (`dense' && !`general') is not yet implemented."))

        stopifnot(h$elems == h$nrow * h$ncol)
        Aelm = readBin(fp, what = h$valuetype, n = h$elems)
        read.check(Aelm, h)
        A = matrix(Aelm, nrow = h$nrow, ncol = h$ncol, byrow = (h$storagetype == "rowmajor"))
    } else if (h$datatype == "sparse") {
        ##
        ## sparse matrix
        ##
        require(Matrix)

        if (h$storagetype == "rowmajor") {
            ##
            ## rowmajor format
            ##
            i = c()
            j = c()
            values = c()
            if (h$matrixtype == "symmetric") {
                ##
                ## symmetric
                ##
                for(k in 1:h$elems) {
                    ij = readBin(fp, what = integer(), n = 2)
                    i = c(i, max(ij))
                    j = c(j, min(ij))
                    values = c(values, readBin(fp, what = h$valuetype, n = 1))
                }
                read.check(i, h)
                read.check(j, h)
                read.check(values, h)
                
                ## oops. Matrix adds replicated elements!!! 
                if (!(all(i >= j) || all(i <= j)))
                    stop(paste("Reading file", filename,
                               ". Both upper and lower part of symmetric sparse matrix",
                               "is specified. Do not know what to do."))
                idx = (i != j)
                ii = i[idx]
                jj = j[idx]
                i = c(i, jj)
                j = c(j, ii)
                values = c(values, values[idx])
            } else if (h$matrixtype == "general" || h$matrixtype == "diagonal") {
                ##
                ## general/diagonal
                ##
                for(k in 1:h$elems) {
                    ij = readBin(fp, what = integer(), n = 2)
                    i = c(i, ij[1])
                    j = c(j, ij[2])
                    values = c(values, readBin(fp, what = h$valuetype, n = 1))
                }
                read.check(i, h)
                read.check(j, h)
                read.check(values, h)

                if (h$matrixtype == "diagonal") {
                    idx = (i == j)
                    i = i[idx]
                    j = j[idx]
                    values = values[idx]
                }
            } else {
                stop("This should not happen.")
            }
        } else {
            ##
            ## columnmajor format
            ##

            ##
            ## other format: (i,j,values)
            ##
            i = readBin(fp, what = integer(0), n = h$elems)
            j = readBin(fp, what = integer(0), n = h$elems)
            values = readBin(fp, what = h$valuetype, n = h$elems)
            read.check(i, h)
            read.check(j, h)
            read.check(values, h)

            if (h$matrixtype == "symmetric") {
                ##
                ## symmetric: lower or upper triangular part is given
                ##

                ## oops. Matrix adds replicated elements!!! 
                if (!(all(i >= j) || all(i <= j)))
                    stop(paste("Reading file", filename,
                               ". Both upper and lower part of symmetric sparse matrix",
                               "is specified. Do not know what to do..."))

                idx = (i != j)
                ii = i[idx]
                jj = j[idx]
                ## yes, this is correct...
                i = c(i, jj)
                j = c(j, ii)
                values = c(values, values[idx])
            } else if (h$matrixtype == "diagonal") {
                idx = (i == j)
                i = i[idx]
                j = j[idx]
                values = values[idx]
            }
        }
        A = sparseMatrix(i = i, j = j, x = values, dims = c(h$nrow, h$ncol), index1=FALSE)
    } else {
        stop("This should not happen.")
    }
    
    close(fp)

    return (A)
}
