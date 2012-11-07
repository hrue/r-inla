sas.smooth.prior = function(
        file = "sas-prior-table.dat",
        h = 0L,
        only.logjac = FALSE,
        dif.lim = 0.25,
        ...)
{
    do.smooth.matrix = function(A, edge.check = FALSE, h = 0, dif.lim,  method = c("lm", "rlm"))
    {
        stopifnot(!missing(dif.lim))
        method = match.arg(method)

        A[is.infinite(A)] = NA
        A.smooth = A
        nr = nrow(A)
        nc = ncol(A)

        library(multicore)
        A.smooth = mclapply(1:(nr*nc),
                function(idx, nr, nc, A, edge.check, method) { 
                    i = (idx-1) %% nr + 1
                    j = ((idx-1) %/% nr) + 1
                    A.smooth = A[i, j]
                    neig.i = unique(pmax(1L, pmin(nr, (i-h):(i+h))))
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
                    
                    n.x = length(unique(x))
                    n.y = length(unique(y))
                    
                    if (length(z) >= 3 && n.x > 1 && n.y > 1 && !is.na(A[i, j])) {
                        if (edge.check) {
                            dif = abs(A[i, j]- z)
                            weights = numeric(length(z))
                            weights[dif > dif.lim] = 0.0
                            weights[dif <= dif.lim] = 1.0
                        } else {
                            weights = rep(1, length(z))
                        }
                        nobs = length(weights[weights > 0])
                        
                        if (nobs >= 6 && n.x > 2 && n.y > 2) {
                            formula = z ~ 1 + x + y + xy + x2 + y2
                        } else if (nobs > 2 && n.x > 1 && n.y > 1) {
                            formula = z ~ 1 + x + y
                        } else if (nobs > 0) {
                            formula = z ~ 1
                        } else {
                            formula = NULL
                        }
                        
                        if (!is.null(formula)) {
                            fit = do.call(method,
                                    args = list(
                                            formula = formula, 
                                            data = list(
                                                    x = x,
                                                    y = y, 
                                                    xy = x*y, 
                                                    x2 = x^2,
                                                    y2 = y^2, 
                                                    z = z),
                                            weights = weights))
                            A.smooth = as.numeric(predict(fit, data.frame(x=0, y=0, x2=0, y2=0, xy=0)))
                        }
                    } 
                    return (A.smooth)
                },
                nr=nr, nc=nc, A=A, edge.check=edge.check, method=method)
        
        A.smooth = matrix(unlist(A.smooth), nr, nc)
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
        method = "lm"
        level = do.smooth.matrix(level, h=h, dif.lim = dif.lim, method=method, ...)
        len = do.smooth.matrix(len, h=h, dif.lim = dif.lim, method=method, ...)
        point = len * do.smooth.matrix(point/len, h=h, edge.check=TRUE, dif.lim = dif.lim, method=method, ...)
    } else {
        library(MASS)
        stopifnot(!is.null(logjac))
        method = "rlm"  ## rlm is in MASS
        logjac = do.smooth.matrix(logjac, h=h, edge.check = FALSE, dif.lim=dif.lim, method=method, ...)
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
sas.remove.logjac = function(file = "sas-prior-table.dat")
{
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
    close(fd)
    
    fd = file(paste(file, "-no-logjac", sep=""), "wb")
    writeBin(as.integer(c(nx, ny, nz)), fd)
    writeBin(skew, fd)
    writeBin(kurt, fd)
    writeBin(c(level), fd)
    writeBin(c(len), fd)
    writeBin(c(point), fd)
    close(fd)
}
