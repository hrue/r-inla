## Nothing to export

## The 'fmesher'-file format is for internal use only
## but it is used in some low-level examples.

`inla.is.fmesher.file` = function(filename)
{
    ## return TRUE if file exists and is a fmesher file
    if (!is.character(filename) || !file.exists(filename))
        return (FALSE)

    fp = file(filename, "rb")
    len.h = readBin(fp, what = integer(), n = 1)
    close(fp)

    ## the only test we can make now,  is if the length of the header is 8.
    if (len.h == 8L) {
        return (TRUE)
    } else {
        return (FALSE)
    }
}

`inla.read.fmesher.file` = function(filename, verbose = FALSE, debug = FALSE)
{
    ##
    ## read a binary-file from fmesher in format specified by FL.
    ##

    if (debug) {
        verbose=TRUE
    }

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
    if (verbose) {
        print(paste("header is", len.h, "integers."))
    }

    h.raw = readBin(fp, what = integer(), n = len.h)
    ## negative entries specify something different not yet defined.
    for(i in 1:len.h) {
        if (h.raw[i] < 0) {
            stop(paste("Entry ", i, " in the header of file ", filename, " is negative: ", h.raw[i],
                       ". Do not know what to do.", sep=""))
        }
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

    if (verbose) {
        print(h)
    }

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
            ## other format: (i, j, values)
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
        A = sparseMatrix(i = i, j = j, x = values, dims = c(h$nrow, h$ncol),
                index1 = FALSE)
        A = inla.as.dgTMatrix(A)
    } else {
        stop("This should not happen.")
    }

    close(fp)

    return (A)
}

`inla.write.fmesher.file` = function(A, filename = tempfile(), verbose = FALSE, debug = FALSE, auto.convert = FALSE)
{
    ##
    ## write a binary-file from fmesher in format specified by FL.
    ##

    if (debug) {
        verbose=TRUE
    }

    if (verbose)
        print(paste("\nOpen file to write", filename))

    version = 0

    ## decide how to treat with almost integers...
    if (auto.convert) {
        if (is.matrix(A)) {
            A = inla.affirm.integer(A)
        } else if (is.list(A)) {
            A$values = inla.affirm.integer(A$values)
        } else if (is(A, "Matrix")) {
            ## Cannot store Matrix matrices as integers; do nothing.
        } else if (is.vector(A)) {
            A = inla.affirm.integer(A)
        } else {
            stop(inla.paste(c("Unknown type of matrix:",
                              deparse(match.call()),
                              "\n class(A) =",
                              class(A)
                              )))
        }
    }

    if (is.matrix(A)) {
        ##
        nrow = dim(A)[1]
        ncol = dim(A)[2]
        elems = nrow*ncol
        datatype = 0
        valuetype = inla.ifelse(is.integer(A), integer(), double())
        matrixtype = 0  ## general
        storagetype = 1 ## columnmajor
    } else if (is.list(A)) {
        nrow = max(A$i)
        ncol = max(A$j)
        datatype = 1 ## sparse
        valuetype = inla.ifelse(is.integer(A$values), integer(), double())
        matrixtype = 0  ## general
        storagetype = 1 ## columnmajor

        i = A$i-1L
        j = A$j-1L
        values = A$values
        elems = length(values)
    } else if (is(A, "Matrix")) {
        ##
        A = inla.as.dgTMatrix(A)
        nrow = dim(A)[1]
        ncol = dim(A)[2]
        datatype = 1 ## sparse
        valuetype = inla.ifelse(is.integer(A), integer(), double())
        matrixtype = 0  ## general
        storagetype = 1 ## columnmajor
        i = A@i
        j = A@j
        values = A@x
        elems = length(i)
    } else if (is.vector(A)) {
        ## diagonal
        nrow = length(A)
        ncol = length(A)
        elems = length(A)
        datatype = 1 ## sparse
        valuetype = inla.ifelse(is.integer(A), integer(), double())
        matrixtype = 2  ## diagonal
        storagetype = 1 ## columnmajor

        i = j = 0:(nrow-1)
        values = A
    } else {
            stop(inla.paste(c("Unknown type of matrix:",
                              deparse(match.call()),
                              "\n class(A) =",
                              class(A)
                              )))
    }

    h = integer(8)
    valuetp = inla.ifelse(identical(valuetype, integer()), 0, 1)
    h = c(version, elems, nrow, ncol, datatype, valuetp, matrixtype, storagetype)
    if (verbose)
        print(h)
    if (debug) {
        if (!is.matrix(A)) {
            print("i (zero-based indexing)")
            print(i)
            print("j (zero-based indexing)")
            print(j)
            print("values")
            print(values)
        }
    }

    fp = file(filename, "wb")
    writeBin(as.integer(length(h)), fp)
    writeBin(as.integer(h), fp)

    if (datatype == 0) {
        ## dense
        if (identical(valuetype, integer())) {
            writeBin(as.integer(as.vector(A)), fp)
        } else {
            writeBin(as.double(as.vector(A)), fp)
        }
    } else {
        ## sparse: columnorder
        writeBin(as.integer(i), fp)
        writeBin(as.integer(j), fp)
        if (identical(valuetype, integer())) {
            writeBin(as.integer(values), fp)
        } else {
            writeBin(as.double(values), fp)
        }
    }
    close(fp)

    return (filename)
}


inla.fmesher.make.dir = function(dir)
{
    dir.start = dir

    ## if this already exists then create one more
    k = 1
    while (file.exists(dir)) {
        dir = paste(dir.start, "-", k, sep="")
        k=k+1
    }

    return(dir)
}


fmesher.write = function(m, prefix, matrixname)
{
    filename = paste(prefix, matrixname, sep="")
    return(inla.write.fmesher.file(m, filename))
}

fmesher.read = function(prefix, matrixname)
{
    filename = paste(prefix, matrixname, sep="")
    if (!file.exists(filename))
        stop(paste("File '", filename, "' does not exist.", sep=""))
    return(inla.read.fmesher.file(filename))
}
