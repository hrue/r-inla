### Various utility functions

`inla.graph.convert.1` = function(in.file, graph.file = "graph.txt", c.indexing = FALSE)
{
    ## convert a file with each line consist of (1-based) node and all
    ## its neigbours. nodes with zero neighbours need not to be
    ## spesified.

    if (c.indexing)
        off = 1
    else
        off = 0
    
    xx = max(scan(in.file))
    n = max(xx)
    rm(xx)
    
    lines = readLines(in.file)
    nlines = length(lines)
    nbs = list(list())
    nnbs = numeric(n)
    
    for(i in 1:n)
        nbs[[i]] = NA

    for(i in 1:nlines) {
        xx = as.integer(unlist(strsplit(readLines("tt")[i], "  *")))
        xx = xx[!is.na(xx)]
        m = length(xx)
        if (m > 0) {
            nnbs[ xx[1] ] = m -1
            nbs[[ xx[1] ]] =  xx[-1]
        }
    }

    sink(graph.file)
    cat(n, '\n')
    for(i in 1:n) {
        if (nnbs[i] > 0) {
            cat(i-off, nnbs[i], nbs[[i]] - off, '\n')
        } else {
            cat(i-off, nnbs[i], '\n')
        }
    }        
    sink()
    rm(xx)
    return (graph.file)
}

`geobugs2inla` = function(adj, num, graph.file="graph.txt", c.indexing = FALSE)
{
    inla.graph.convert.2(adj, num, graph.file, c.indexing = c.indexing)
}

`inla.graph.convert.2` = function(adj, num, graph.file="graph.txt", c.indexing = FALSE)
{
    ## A function for converting GeoBUGS adjacency data into the INLA
    ## graph format. Kindly provided by Aki Havunlinna tkk.fi; thanks.

    if (c.indexing)
        off = 1
    else
        off = 0

    sink(graph.file)
    len <- length(num)
    cat(len, '\n')
    spl <- split(adj-1, rep(1:len, num))
    for(i in 1:len)  {
        ni <- num[i]
        cat(i-off, ni, spl[[i]] + (1-off),'\n')
    }
    sink()
    return (graph.file)
}

`inla.Cmatrix2matrix` = function(C, symmetric = TRUE)
{
    ## convert a Cij matrix into a matrix.
    
    stopifnot(length(C$i) == length(C$j))
    stopifnot(length(C$i) == length(C$Cij))

    n = max(c(C$i, C$j))
    A = matrix(0,n,n)
    for(k in 1:length(C$i)) {
        ## this could be vectorised with more work...
        i = C$i[k]
        j = C$j[k]
        v = C$Cij[k]
        A[i,j] = v
        if (symmetric)
            A[j,i] = v
    }
    return (A)
}

`inla.matrix2Cmatrix` = function(Q, symmetric = TRUE)
{
    ## convert a symmetric(if TRUE) matrix into the form:
    ## list(i=,j=,Cij=)

    if (!is.matrix(Q))
        stop("Argument must be a matrix")
    
    n = dim(Q)
    if (n[1] != n[2])
        stop(paste("Matrix must be a square matrix, dim(Q) =", dim(Q)))
    
    n = n[1]
    ii = c()
    jj = c()
    Cij = c()

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
            Cij = c(Cij, Q[i, idx])
        }
    }

    return (list(i=ii, j=jj, Cij = Cij))
}

`inla.Cmatrix2file` = function(Cmatrix, filename = NULL, c.indexing = FALSE)
{
    if (c.indexing)
        off = 1
    else
        off = 0
    
    if (is.null(filename)) {
        filename = tempfile()
    }
    write(t(cbind(Cmatrix$i-off, Cmatrix$j-off, Cmatrix$Cij)), ncolumns=3, file = filename)

    return (filename)
}


`inla.matrix2graph` = function(Q, graph.file = "graph.txt", c.indexing = FALSE)
{
    if (c.indexing)
        off = 1
    else
        off = 0

    if (!is.matrix(Q))
        stop("Argument must be a matrix")
    
    n = dim(Q)
    if (n[1] != n[2])
        stop(paste("Matrix must be a square matrix, dim(Q) =", dim(Q)))
    
    n = dim(Q)[1]
    diag(Q) = 0
    sink(graph.file)

    cat(n, "\n")
    for(i in 1:n) {
        idx = which(Q[i,] != 0)
        nb = length(idx)
        if (nb > 0)
            cat(i-off, nb, idx-off, "\n")
        else
            cat(i-off, nb, "\n")
    }
    sink()
    return (graph.file)
}

`inla.num` = function(i, width = 8, digits=4)
{
    ##return (gsub(" ", "0", format(i,width=width)))
    return(formatC(i, format="g", width = width, flag = "0", digits=digits))
}

`inla.trim` = function(string)
{
    ## trim leading and trailing whitespaces and dots. there is a
    ## function in R.oo called `trim' that do this, but I don't want
    ## INLA to be dependent on R.oo. This function also works the
    ## string is a list of strings.

    string = gsub("^[ .]+", "", string)
    string = gsub("[ .]+$", "", string)
    return (string)
}

`inla.namefix` = function(string)
{
    ## makes inla-name from string
    if (FALSE) {
        string = inla.trim(string)
        string = gsub("[^A-Za-z0-9.*]+", ".", string)
        string = gsub("[.]+", ".", string)
    }
    return (string)
}

`inla.nameunfix` = function(string)
{
    ## makes a nice printable version of STRING whic is a `name'
    if (FALSE) {
        string = gsub("[ .]+", " ", string)
        string = gsub("([0-9]) ([0-9])", "\\1.\\2", string)
    }
    return (string)
}

`inla.strcmp` = function(s,ss)
{
    ## compare two strings
    return (s == ss)
}

`inla.strcasecmp` = function(s,ss)
{
    ## compare two strings, ignore case
    return (tolower(s) == tolower(ss))
}

`inla.pause` = function(msg = NULL) 
{
    ## just print a msg and wait for the next return
    if (!is.null(msg))
        cat(msg, "\n")
    scan(quiet=TRUE, multi.line=TRUE,what=character(0))
}
`inla.call.builtin` = function()
{
    if (inla.os("mac"))
        fnm = system.file("bin/mac/inla", package="INLA")
    else if (inla.os("linux"))
        fnm = system.file("bin/linux/inla", package="INLA")
    else if (inla.os("windows"))
        fnm = system.file("bin/windows/inla.exe", package="INLA")
    else
        stop("Unknown OS")

    if (file.exists(fnm))
        return (fnm)
    else
        stop(paste("INLA installation error; no such file",fnm))
}

`inla.paste` = function(strings, sep = " ")
{
    ## merge strings
    res = NULL
    for(s in strings) {
        if (is.null(res))
            res = s
        else
            res = paste(res, s, sep=sep)
    }
    return(res)
}

`inla.compare.results` = function(dir.inla = NULL, dir.mcmc = NULL)
{
    ## This function is useful for comparing the results of INLA with
    ## those obtained from a MCMC run (ie running `inla -m mcmc ...')

    if (is.null(dir.inla)) {
        print("Enter directory for INLA results ")
        dir.inla = scan(nmax=1,what = character(0))
    }
    if (!file.exists(dir.inla))
         stop(paste(" *** ERROR *** No such directory ", dir.inla))
    if (file.exists(paste(dir.inla, "/", "results.files", sep="")))
        dir.inla = paste(dir.inla, "/", "results.files", sep="")

    if (is.null(dir.mcmc)) {
        print("Enter directory for MCMC (from ``inla -m mcmc'') results ")
        dir.mcmc = scan(nmax=1,what = character(0))
    }
    if (!file.exists(dir.mcmc))
         stop(paste(" *** ERROR *** No such directory ", dir.mcmc))
    if (file.exists(paste(dir.mcmc, "/", "results.files", sep="")))
        dir.mcmc = paste(dir.mcmc, "/", "results.files", sep="")
 
    options(device.ask.default=TRUE)
    dev.new()

    dd = c()
    for(d in dir(dir.inla))
        if (file.exists(paste(dir.inla, "/", d, "/marginal-densities.dat", sep="")))
            dd = c(dd, d)

    quit = FALSE
    while(TRUE) {

        cat("\nChose directory to view\n")
        i=1
        for(a in dd) {
            cat("[",i,"] ", a, "\n", sep="")
            i = i+1
        }
        cat("\n(<0: histogram, >0: density: q:quit) : ")

        ok = FALSE
        while(!ok) {
            id = readline("")
            if (id == "q" || id == "quit") {
                ok = TRUE
                quit = TRUE
            } else {
                id = as.integer(id)
                if (!is.na(id))
                    ok = TRUE
                else
                    cat("***Input error.\n")
            }
        }

        if (quit)
            return (invisible())
            
        if (id < 0) {
            hist = TRUE
            id = abs(id)
        } else {
            hist = FALSE
        }
    
        d.inla = paste(dir.inla, "/", dd[id], sep="")
        d.mcmc = paste(dir.mcmc, "/", dd[id], sep="")

        fnm = paste(d.inla,"/","marginal-densities.dat",sep="")
        f = paste(d.inla,"/","TAG",sep="")
        TAG = inla.ifelse(file.exists(f), inla.paste(scan(paste(d.inla,"/","TAG",sep=""),
            what=as.character())), "TAG")
        if (file.exists(fnm)) {
            ##print(paste("read marginals from ",fnm))
            marg = read.table(fnm)

            n.marg = dim(marg)[1]
            n.points = dim(marg)[2]

            fnm = paste(d.mcmc,"/","trace.dat",sep="")
            if (file.exists(fnm)) {
                samples = read.table(fnm)
                di = dim(samples)
                samples = samples[-dim(samples)[1],]
                ## with ncol=1, the previous line remove dimension attributes
                samples = as.matrix(samples,nrow=di[1]-1,ncol=di[2])
                n.samples = dim(samples)[1]

                if (dim(samples)[2] == n.marg) {
                    if (n.marg <= 4)
                        par(mfrow=c(n.marg,1))
                    else
                        par(mfrow=c(4,3))
        
                    for(i in 1:n.marg) {
                        x.values = as.double(marg[i, seq(2, n.points, by=2)])
                        f.values = as.double(marg[i, seq(3, n.points, by=2)])
                        plot(x.values, f.values, type="l", lty=1, lwd=2,
                             main=paste("Marginal for ", marg[i,1]," ", TAG, sep=""),
                             sub = fnm, frame.plot = FALSE)
                        if (hist)
                            hist(samples[,i], n = 50, add = TRUE, probability = TRUE)
                        else
                            lines(density(samples[,i]), lty=3, lwd=2)
                    }
                }
                else
                    print("n.marg != columns(samples) ", n.marg,dim(samples)[2])
            }
            else
                print(paste(" *** ERROR *** file does not exists: ",fnm))

            rm(marg)
            rm(samples)
        }
        else
            print(paste(" *** ERROR *** file does not exits: ",fnm))
    }
}
`inla.only.for.developers` = function()
{
    require(R.utils)
    if (getUsername.System() != "hrue")
        stop("This function is for developers only.")
    return (invisible())
}
`inla.my.update` = function(dir = "~/p/inla/google-code/inla/rinla/R")
{
    inla.only.for.developers()

    files = dir(dir, pattern = "[.][Rr]$")
    for (f in files)
        source(paste(dir, "/", f, sep=""))
    cat("Source in ", dir, " loaded (", length(files), " files)\n", sep="")

    return (invisible())
}
`inla.remove` = function(name, from)
{
    ## remove NAME FROM. Works for both lists and data.frames

    if (is.list(from) || is.data.frame(from))
        for(nm in name)
            if (length(grep(nm, names(from))) > 0)
                inla.eval(paste("from$", nm, "=NULL", sep=""))

    return (from)
}

`inla.ifelse` = function(test, yes, no)
{
    if (test)
        return (yes)
    else
        return (no)
}

`inla.lattice2node.mapping` = function(nrow, ncol)
{
    ## return a matrix with the mapping

    stopifnot( nrow > 0 && ncol > 0 )

    mapping = matrix(NA, nrow=nrow, ncol=ncol)
    for(i in 1:nrow) {
        j = 1:ncol
        mapping[i,j] = inla.lattice2node(i, j, nrow, ncol)
    }
    return (mapping)
}
`inla.node2lattice.mapping` = function(nrow, ncol)
{
    stopifnot( nrow > 0 && ncol > 0 )

    return (inla.node2lattice(1:(nrow*ncol), nrow, ncol))
}

`inla.lattice2node` = function(irow, icol, nrow, ncol)
{
    ## convert from a lattice point (irow, icol) to a node-number in
    ## the graph; similar to the GMRFLib_lattice2node()-function in
    ## GMRFLib.  Indices goes from irow=1...nrow, to icol=1...ncol and
    ## node=1.....nrow*ncol.

    stopifnot( nrow > 0 && ncol > 0 )

    ##return ((icol-1) + (irow-1)*ncol + 1)
    return ((irow-1) + (icol-1)*nrow + 1)
}

`inla.node2lattice` = function(node, nrow, ncol)
{
    ## convert from a node-number in the graph to a lattice point
    ## (irow, icol); similar to the GMRFLib_node2lattice()-function in
    ## GMRFLib.  Indices goes from irow=1...nrow, to icol=1...ncol and
    ## node=1.....nrow*ncol.

    stopifnot( nrow > 0 && ncol > 0 )

    ##irow = (node - 1) %/% ncol
    ##icol = (node - 1) %% ncol
    icol = (node - 1) %/% nrow
    irow = (node - 1) %% nrow

    return (list(irow = irow+1, icol = icol+1))
}
            
`inla.matrix2vector` = function(a.matrix)
{
    ## utility function for mapping a matrix to inla's internal `node'
    ## representation by inla.lattice2node() and inla.node2lattice()

    if (!is.matrix(a.matrix))
        stop("Argument must be a matrix")

    return (as.vector(a.matrix))
}
`inla.vector2matrix` = function(a.vector, nrow, ncol)
{
    ## utility function for mapping from inla's internal `node'
    ## representation, inla.lattice2node() and inla.node2lattice(),
    ## and to a matrix

    n = length(a.vector)
    if (missing(nrow) && !missing(ncol))
        nrow = n %/% ncol
    if (!missing(nrow) && missing(ncol))
        ncol = n %/% nrow
    if (n != nrow*ncol)
        stop(paste("Length of vector", n, "does not equal to nrow*ncol", nrow*ncol))
    
    return (matrix(a.vector, nrow, ncol))
}
`inla.squishplot` = function (xlim, ylim, asp = 1, newplot = TRUE) 
{
    ## This function is a copy from package TeachingDemos

    if (length(xlim) < 2) 
        stop("xlim must be a vector of length 2")
    if (length(ylim) < 2) 
        stop("ylim must be a vector of length 2")
    if (newplot) 
        plot.new()
    tmp <- par(c("plt", "pin", "xaxs", "yaxs"))
    if (tmp$xaxs == "i") {
        xlim <- range(xlim)
    } else {
        tmp.r <- diff(range(xlim))
        xlim <- range(xlim) + c(-1, 1) * 0.04 * tmp.r
    }
    if (tmp$yaxs == "i") {
        ylim <- range(ylim)
    } else {
        tmp.r <- diff(range(ylim))
        ylim <- range(ylim) + c(-1, 1) * 0.04 * tmp.r
    }
    tmp2 <- (ylim[2] - ylim[1])/(xlim[2] - xlim[1])
    tmp.y <- tmp$pin[1] * tmp2 * asp
    if (tmp.y < tmp$pin[2]) {
        par(pin = c(tmp$pin[1], tmp.y))
        par(plt = c(tmp$plt[1:2], par("plt")[3:4]))
    } else {
        tmp.x <- tmp$pin[2]/tmp2/asp
        par(pin = c(tmp.x, tmp$pin[2]))
        par(plt = c(par("plt")[1:2], tmp$plt[3:4]))
    }
    return(invisible(tmp["plt"]))
}

`inla.display.matrix` = function(x, wrap=TRUE, xaxt=FALSE, yaxt=FALSE, ...)
{
    ## display a matrix as an image with correct layout and autoscaling
    
    if (!is.matrix(x)) stop("First argument must be a matrix")
    n = dim(x)
    y = x
    if (wrap) {
        for(j in 1:n[2]) 
            y[1:n[1],j] = x[n[1]:1,j]
    }

    ## use the image.plot-function in package fields; its much better...
    inla.squishplot(c(0,1), c(0,1), n[1]/n[2])
    if (require("fields", quietly = TRUE)) {
        image.plot(t(y), col=gray(seq(0,1,len=256)), bty="n", xaxt="n", yaxt="n", ...)
    } else {
        warning("Please install package `fields'")
        image(t(y), col=gray(seq(0,1,len=256)), bty="n", xaxt="n", yaxt="n", ...) 
    }

    box()
    if (xaxt) {
        title(xlab="")
        nn = (n[2]+1)%/%2
        axis(1, at=c(0,(nn-1)/(n[2]-1), 1), labels=as.character(c(1, nn, n[2])))
    }
    if (yaxt) {
        title(ylab="")
        nn = (n[1]+1)%/%2
        axis(2, at=c(0,(nn-1)/(n[1]-1), 1), labels=as.character(c(1, nn, n[1])))
    }
}

`inla.group` = function(x, n = 25, method = c("cut", "quantile"))
{
    ## group covariates into N groups using method "quantile" or
    ## "cut", i.e., the functions quantile() or cut().  the cut use
    ## 'even length' wheras the 'quantile' use even length quantiles.

    ## I make the "cut" default, as then we have control over the
    ## minimum distance between each cell, whereas the "quantile"
    ## approach, we do not have a such control.

    if (n < 1)
        stop("Number of groups must be > 0")
    
    if (n == 1)
        return (rep(median(x), length(x)))
    
    method = match.arg(method)
    if (method == "cut") {
        ## use equal length
        a = cut(x, n)
    } else {
        ## use break-points corresponding to the quantiles
        aq = quantile(x, probs = c(0, ppoints(n-1),1))
        a = cut(x, breaks = as.numeric(aq), include.lowest=TRUE)
    }
    ## the rest is then the same
    nlev = nlevels(a)
    xx = list()
    for(i in 1:nlev)
        xx[[i]] = list()
    
    for(i in 1:length(x))
        xx[[as.numeric(a[i])]] = c(unlist(xx[[as.numeric(a[i])]]),x[i])
    values = numeric(nlev)
    
    ff.local = function(xx) {
        if (length(xx) > 0)
            return (median(xx))
        else
            return (NA)
    }
    values = unlist(sapply(xx, ff.local))
    x.group = as.numeric(values[as.numeric(a)])

    return (x.group)
}

`inla.group.old` = function(x,n)
{
    ### old version
    cutpoints = seq(range(x)[1],range(x)[2], length.out=(n+1))
    lev = (cutpoints[1:(n)]+cutpoints[2:(n+1)])/2
    int = cut(x, breaks=cutpoints, include.lowest=TRUE, labels=FALSE)
    return(lev[int])
}

`inla.scale` = function(x)
{
    ## return [(x - mean(x))/sd(x)]
    return ((x-mean(x))/sd(x))
}

`inla.trim.family` = function(family)
{
    ## remove `_' and `.' and space and tabs... and then convert to lowercase
    return (tolower(gsub("[_. \t]+", "", family)))
}

`inla.is.list.of.lists` = function(a.list)
{
    ## return TRUE if `a.list' is a list of lists, otherwise FALSE
    if (length(a.list) == 0)
        return (FALSE)
    else
        return (all(sapply(a.list, is.list)))
}

`inla.as.list.of.lists` = function(a)
{
    if (is.matrix(a) || is.data.frame(a)) {
        return(as.list(as.data.frame(a)))
    } else if (inla.is.list.of.lists(a) || TRUE) {
        return (a)
    } else {
        stop("Argument if of unknown type; do not know what to do...")
    }
}

`inla.replicate.list` = function(a.list, nrep)
{
    ## make a list() into a list of lists.  example: list(a=1) =>
    ## list(list(a=1), list(a=1), ...)

    if (nrep <= 0)
        stop("nrep must be > 0")

    return (rep(list(a.list), nrep))
}

`inla.tictac` = function(num, elms = "|/-\\|/-\\")
{
    len = nchar(elms)
    i = (num %% len) + 1
    return(substr(elms,i,i))
}

`inla.2list` = function(x = NULL)
{
    ## convert a vector `x' or string 'x' into the string ``c(x[1],
    ## x[2], x[3])'' so it can be evaluated using
    ## inla.eval()

    ## example:
    ## > x=1:3
    ## > inla.2list(x)
    ## [1] "c(1,2,3)"
    
    if (is.null(x))
        return (NULL)

    if (length(x) == 0)
        return ("numeric(0)")
    s = inla.paste(as.character(x))
    s = gsub("^[ ]+", "", s)
    s = gsub("[ ]+$", "", s)
    s = gsub("c[ ]*[(][ ]*", "", s)
    s = gsub("[ ]*[)][ ]*", "", s)
    return (paste("c(", gsub("[, ]+", ",", s, ")"), ")", sep=""))
}

`inla.reload` = function(lib.loc = NULL)
{
    ## detach the INLA-library if its already loaded and then reload
    ## load the INLA-library

    path = dirname(searchpaths()[grep("/INLA$", searchpaths())])
    pkg = "package:INLA"
    if (any(search() == pkg))
        detach(pos = match(pkg, search()))
    if (is.null(lib.loc) && length(path) > 0)
        lib.loc = path
    library(INLA, lib.loc = lib.loc)
}

`inla.even` = function(n)
{
    return (inla.divisible(n,by=2))
}

`inla.odd` = function(n)
{
    return (inla.divisible(n,by=-2))
}

`inla.divisible` = function(n, by=2)
{
    ### if by>0, return TRUE if `n' is divisible by `by', and if by<0,
    ### return TRUE if `n' is not divisible by `-by'. if by==0, return
    ### TRUE.
    if (by == 0)
        return (rep(TRUE, length(n)))
    if (by > 0)
        return ((n%%by) == 0)
    else 
        return ((n%%(-by)) != 0)
}

`inla.one.of` = function(family, candidates)
{
    ## check if family is one of the canidates
    return (any(inla.trim.family(family) == inla.trim.family(candidates)))
}

`inla.idx` = function(idx, n = max(idx), group = rep(1, n), ngroup = max(group),
        replicate = rep(1, n), nrep = max(replicate))
{
    ## this function might be useful to convert from (idx,group,rep)
    ## to idx, in the same way as done internally in inla.R
    
    stopifnot(n >= 1)
    stopifnot(ngroup >= 1)
    stopifnot(nrep >= 1)
    stopifnot(all(group >= 1))
    stopifnot(all(replicate >= 1))
    stopifnot(all(idx >= 1))
    stopifnot(all(idx <= n))
    stopifnot(ngroup >= max(group))
    stopifnot(nrep >= max(replicate))
    
    return (idx + (group-1)*n + (replicate-1)*n*ngroup)
}

`inla.get.HOME` = function()
{
    return (as.character(inla.ifelse(inla.os("windows"), gsub("\\\\", "/", Sys.getenv("USERPROFILE")), Sys.getenv("HOME"))))
}

`inla.get.USER` = function()
{
    u = ""
    for (U in c("USER", "USERNAME", "LOGNAME")) {
        u = Sys.getenv(U)
        if (u != "")
            break;
    }
    if (u == "")
        u = "UnknownUserName"

    return (as.character(u))
}

`inla.remote` = function()
{
    if (inla.os("windows") && !inla.cygwin.check.path()) {
        cat("\n\n\tRemote computing in INLA from Windows require CYGWIN.\n")
        cat("\tThe CYGWIN-installation is not found. Please install it from www.cygwin.com\n")
        cat("\tMake sure to install also all `ssh' components; search for it, select them and install.\n")
        cat("\tIts easiest to use the standard installation path; c:/cygwin, otherwise give the new path using\n")
        cat("\t\tinla.setOption(\"cygwin\", \"NEWPATH\")\n")
        cat("\n")
        return (NULL)
    }
            
    inla=system.file("bin/remote/inla.remote", package="INLA")
    inla.setOption("inla.call", inla)
    
    f = paste(inla.get.HOME(), "/.inlarc", sep="")
    if (!file.exists(f)) {
        if (inla.os("windows"))
            ini=system.file("bin/remote/dot.inlarc-example-windows", package="INLA")
        else
            ini=system.file("bin/remote/dot.inlarc-example", package="INLA")
        
        if (!file.exists(ini))
            stop(paste("Missing file in installation:", ini, "This should not happen..."))
        file.copy(ini, f)
        cat(paste("Remote computing: Please edit remote-host details in\n\t\t", f, "\n"))
    } else {
        cat(paste("Remote computing is defined in\n\t\t", f, "\n"))
    }
    if (inla.os("windows")) {
        cat("\n\tIf you need to copy ssh-keys, you can do that using this script from within a cygwin-terminal\n")
        cat(paste("\t\t", inla.cygwin.map.filename(system.file("bin/remote/ssh-copy-id", package="INLA")), "\n\n"))
    }

    return (invisible())
}

`inla.eval` = function(command)
{
    eval.parent(parse(text=command))
}

`inla.tempfile` = function(pattern = "file", tmpdir = tempdir())
{
    ## just replace \ in Windows with /

    return (gsub("\\\\", "/", tempfile(pattern, tmpdir)))
}

`inla.tempdir` = function()
{
    ## just replace \ in Windows with /

    return (gsub("\\\\", "/", tempdir()))
}

`inla.ssh.copy.id` = function () 
{
    ## print the path to the copy-id script

    f = system.file("bin/remote/ssh-copy-id", package="INLA")
    if (inla.os("windows"))
        return (inla.cygwin.map.filename(f))
    else
        return (f)
}

`inla.formula2character` = function(formula)
{
    ## convert a formula to characters without the 500 character
    ## constraint. Remove the ``()'' at the end.
    return (gsub("[(][)]$","", inla.paste(deparse(formula))))
}

`inla.qinv` = function(Cmatrix) {

    if (is.matrix(Cmatrix)) {
        qinv.file = inla.Cmatrix2file(inla.matrix2Cmatrix(Cmatrix), c.indexing = TRUE)
    } else {
        qinv.file = inla.Cmatrix2file(Cmatrix, c.indexing = TRUE)
    }
        
    if (inla.os("linux") || inla.os("mac")) {
        s = system(paste(shQuote(inla.getOption("inla.call")), "-s -m qinv", qinv.file), intern=TRUE)
    } else if(inla.os("windows")) {
        s = system(paste(shQuote(inla.getOption("inla.call")), "-s -m qinv", qinv.file), intern=TRUE)
    } else {
        stop("\n\tNot supported architecture.")
    }

    unlink(qinv.file)

    return (s)
}
