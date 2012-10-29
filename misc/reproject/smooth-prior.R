smooth.sas.prior = function(file = "sas-prior-table.dat", multicore = TRUE, h = 0L, only.logjac = FALSE, dif.lim = 0.25, ...)
{
    library(quantreg)
    library(MASS)

    do.smooth.matrix = function(A, edge.check = FALSE, h = 0, dif.lim = 1) 
    {
        ##
        A.smooth = A
        nr = nrow(A)
        nc = ncol(A)

        for(i in 1L:nr) {
            print(i)
            neig.i = unique(pmax(1L, pmin(nr, (i-h):(i+h))))
            for(j in 1L:nc) {
                neig.j = unique(pmax(1L, pmin(nc, (j-h):(j+h))))
                x = rep(neig.i,  length(neig.j))
                y = rep(neig.j, each = length(neig.i))
                z = A[ x + (y-1)*nr ]

                remove.idx = which(is.na(z))
                if (length(remove.idx) > 0) {
                    x = x[-remove.idx] -i
                    y = y[-remove.idx] -j
                    z = z[-remove.idx]
                } else {
                    x = x - i
                    y = y - j
                }

                if (length(z) >= 6 && length(unique(x)) > 3 && length(unique(y)) > 3 && !is.na(A[i, j])) {
                    if (edge.check) {
                        dif = abs(A[i, j]- z)
                        weights = numeric(length(z))
                        weights[dif > dif.lim] = 0.0
                        weights[dif <= dif.lim] = 1.0
                    } else {
                        weights = rep(1, length(z))
                    }

                    ## or rq/lqs/rlm
                    fit = lm(z ~ 1 + x + y + xy + x2 + y2,
                            data = list(
                                    x = x,
                                    y = y, 
                                    xy = x*y, 
                                    x2 = x^2,
                                    y2 = y^2, 
                                    z = z),
                            weights = weights)
                    A.smooth[i, j] = as.numeric(predict(fit, data.frame(x=0, y=0, x2=0, y2=0, xy=0)))
                } 
            }
        }
        return (A.smooth)
    }


    fd = file(file, "rb")
    x = readBin(fd, integer(), 3L)
    nx = x[1]
    ny = x[2]
    nz = x[3]

    skew = readBin(fd, double(), nx)
    kurt = readBin(fd, double(), ny)
    level = matrix(readBin(fd, double(), nz), nx, ny)
    len = matrix(readBin(fd, double(), nz), nx, ny)
    point = matrix(readBin(fd, double(), nz), nx, ny)
    logjac = readBin(fd, double(), nz)
    if (length(logjac) == nz) {
        logjac = matrix(logjac, nx, ny)
    } else {
        logjac = NULL
    }
    close(fd)
    
    if (!only.logjac) {
        if (multicore) {
            library(multicore)
            cat("smooth...\n")
            llevel = mcparallel(do.smooth.matrix(level, h=h, dif.lim = dif.lim, ...), name = "level")
            llen = mcparallel(do.smooth.matrix(len, h=h, dif.lim = dif.lim, ...), name = "len")
            ppoint = mcparallel(do.smooth.matrix(point/len, h=h, edge.check=TRUE, dif.lim = dif.lim, ...), name = "point")
            ans = collect(list(llevel, llen, ppoint))
            level = ans$level
            len = ans$len
            point = len*ans$point
        } else {
            level = do.smooth.matrix(level, h=h, dif.lim = dif.lim, ...)
            len = do.smooth.matrix(len, h=h, dif.lim = dif.lim, ...)
            point = len * do.smooth.matrix(point/len, h=h, edge.check=TRUE, dif.lim = dif.lim, ...)
        }
    } else {
        stopifnot(!is.null(logjac))
        logjac = do.smooth.matrix(logjac, h=h, edge.check = edge.check, dif.lim=dif.lim, ...)
    }

    fd = file(paste(file, "-SMOOTH-h-", as.integer(h), sep=""), "wb")
    writeBin(as.integer(c(nx, ny, nz)), fd)
    writeBin(skew, fd)
    writeBin(kurt, fd)
    writeBin(c(level), fd)
    writeBin(c(len), fd)
    writeBin(c(point), fd)
    if (!is.null(logjac)) {
        writeBin(c(logjac), fd)
    }
    close(fd)
}
