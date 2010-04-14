### RCSId = "$Id: binary.R,v 1.14 2010/01/25 08:49:08 hrue Exp $"
### Binary I/O and interpretation of vectors

`inla.interpret.vector` =
    function(xx)
{
    len = length(xx)
    ii = 1
    res = numeric()
    by = 2
    while(ii <= len) {
        np = xx[ii+1]
        if (np > 0) {
            start = ii+2
            end = ii+2*np+1
            res = cbind(res,cbind(xx[seq(start, end, by=by)], xx[seq(start+1, end, by=by)]))
            ii = end + 1
        } else {
            ii = ii + 2
        }
    }
    rm(xx)
    return(res)
}

`inla.interpret.vector.idx` =
    function(xx)
{
    len = length(xx)
    ii = 1
    res = numeric()
    by = 2
    while(ii <= len) {
        idx = xx[ii]
        np = xx[ii+1]
        if (np > 0) {
            start = ii+2
            end = ii+2*np+1
            res = cbind(res,cbind(rep(idx, length(seq(start,end,by=by))),
                xx[seq(start, end, by=by)], xx[seq(start+1, end, by=by)]))
            ii = end + 1
        } else {
            ii = ii + 2;
        }
    }
    rm(xx)
    return(res)
}

 `inla.interpret.vector.list` =
    function(xx)
{
    len = length(xx)
    ind = ii = 1
    res = list()
    by=2
    while(ii <= len) {
        np = xx[ii+1]
        if (np>0) {
            start = ii+2
            end = ii+2*np+1
            res[[ind]] = cbind(xx[seq(start, end, by=by)], xx[seq(start+1, end, by=by)])
            ind = ind+1
            ii = end + 1
        } else {
            res[[ind]] = NA
            ii = ii + 2
        }
    }
    rm(xx)
    return(res)
}

`inla.read.binary.file` =
    function(file)
{
    len = file.info(file)$size

    if (is.na(len))
        stop(paste("File does not exits: ", file))

    fp = file(file,"rb")
    xx = readBin(fp, numeric(), n=len)
    close(fp)

    return (xx)
}
