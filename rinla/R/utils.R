## Various utility functions

`inla.numlen` = function(n)
{
    ## number of digits required to represent a specific (integer)
    ## numb{e{r
    return (floor(log10(max(abs(n))))+1)
}

`inla.num` = function(x, width = if (length(x) > 1) inla.numlen(x) else 8, digits=max(4, width))
{
    ## format numbers using preceeding zeros.
    ## > inla.num(sqrt(2))
    ##  [1] "0001.414"
    ## > inla.num(2)
    ##  [1] "00000002"

    ## for a sequence of numbers, the width is calculated
    ## automatically.> inla.num(1:10)
    ## > inla.num(1:10)
    ##  [1] "01" "02" "03" "04" "05" "06" "07" "08" "09" "10"

    return(formatC(x, format="g", width = width, flag = "0", digits=digits))
}


`inla.trim` = function(string)
{
    ## trim leading and trailing whitespaces and dots. there is a
    ## function in R.oo called `trim' that do this, but I don't want
    ## INLA to be dependent on R.oo. This function also works the
    ## string is a list of strings.

    string = gsub("^[ \t.]+", "", string)
    string = gsub("[ \t.]+$", "", string)
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

`inla.strcmp` = function(s, ss)
{
    ## compare two strings
    return (s == ss)
}

`inla.strncmp` = function(s, ss)
{
    ## compare two strings
    if (length(s) == 1 && length(ss) > 1) {
        ans = c()
        for(i in 1:length(ss)) {
            ans = c(ans,  inla.strncmp(s, ss[i]))
        }
        return (ans)
    } else if (length(s) > 1 && length(ss) > 1) {
        stop("length(s) > 1 && length(ss) > 1: not allowed.")
    } else {
        return (substr(s, 1, nchar(ss)) == ss)
    }
}

`inla.strcasecmp` = function(s, ss)
{
    ## compare two strings, ignore case
    return (tolower(s) == tolower(ss))
}

`inla.strncasecmp` = function(s, ss)
{
    ## compare two strings, ignore case
    if (length(s) == 1 && length(ss) > 1) {
        ans = c()
        for(i in 1:length(ss)) {
            ans = c(ans,  inla.strncasecmp(s, ss[i]))
        }
        return (ans)
    } else if (length(s) > 1 && length(ss) > 1) {
        stop("length(s) > 1 && length(ss) > 1: not allowed.")
    } else {
        return (substr(tolower(s), 1, nchar(ss)) == tolower(ss))
    }
}

`inla.pause` = function(msg = NULL)
{
    ## just print a msg and wait for the next return
    if (!is.null(msg))
        cat(msg, "\n")
    scan(quiet=TRUE, multi.line=TRUE, what=character(0))
}
`inla.call.builtin` = function()
{
    if (inla.os("mac")) {
        fnm = system.file(paste("bin/mac/", inla.os.32or64bit(), "bit/inla", sep=""), package="INLA")
    } else if (inla.os("linux")) {
        fnm = system.file(paste("bin/linux/inla", inla.os.32or64bit(), sep=""), package="INLA")
    } else if (inla.os("windows")) {
        fnm = system.file(paste("bin/windows/inla", inla.os.32or64bit(), ".exe", sep=""), package="INLA")
    } else {
        stop("Unknown OS")
    }

    if (file.exists(fnm)) {
        return (fnm)
    } else {
        stop(paste("INLA installation error; no such file", fnm))
    }
}
`inla.fmesher.call.builtin` = function()
{
    if (inla.os("mac")) {
        fnm = system.file(paste("bin/mac/", inla.os.32or64bit(), "bit/fmesher", sep=""), package="INLA")
    } else if (inla.os("linux")) {
        fnm = system.file(paste("bin/linux/fmesher", inla.os.32or64bit(), sep=""), package="INLA")
    } else if (inla.os("windows")) {
        fnm = system.file(paste("bin/windows/fmesher", inla.os.32or64bit(), ".exe", sep=""), package="INLA")
    } else {
        stop("Unknown OS")
    }

    if (file.exists(fnm)) {
        return (fnm)
    } else {
        stop(paste("INLA installation error; no such file", fnm))
    }
}

`inla.paste` = function(strings, sep = " ")
{
    return (paste(strings, collapse = sep, sep=""))
}

`inla.only.for.developers` = function()
{
    if (!(Sys.getenv("USER") %in%
                    c("hrue", "martino", "finnkrl", "finn",
                      "danielsimpson", "rieblera", "ariebler")))
        warning("This function is for developers only...")
    return (invisible())
}

`inla.my.update` = function(dir, binaries=FALSE, ignore.regexp= NULL)
{
    ##  Set binaries=TRUE to set the inla.call and fmesher.call options
    ##  To override the default binaries path, set binaries="/the/path/bin"

    inla.only.for.developers()

    if (Sys.getenv("USER") %in% c("hrue")) {
        dir.default = "~/p/inla/google-code/inla/rinla/R"
        bin.default = "~/p/inla/work/local/bin"
    } else if (Sys.getenv("USER") %in% c("elias")) {
        dir.default = "~/inla-project/source/inla/rinla/R"
        bin.default = "~/inla-project/compile/local/bin"
    } else if (Sys.getenv("USER") %in% c("rieblera", "ariebler")) {
        dir.default = "~/inla/rinla/R"
        bin.default = "~/local/bin"
    } else {
        dir.default = "~/hg/inla/rinla/R"
        bin.default = "~/local/bin"
    }
    if (!missing(binaries)) {
        if (is.character(binaries)) {
            bin.default = binaries
        }
        binaries = TRUE
    } else {
        binaries = FALSE
    }
    if (missing(dir)) {
        dir = path.expand(dir.default)
    }
    if (binaries) {
        bin.path = path.expand(bin.default)
    }

    files = dir(dir, pattern = "[.][Rr]$")
    ## remove files matching 'ignore.regexp'
    if (!is.null(ignore.regexp)) {
        idx = grep(ignore.regexp, files)
        files = files[-idx]
    }

    ## source the files in a temporary environment and not the globalenv
    tmp.env = new.env()
    for (ff in files) {
        fff = paste(dir, "/", ff, sep="")
        local({source(fff, local=TRUE)}, envir = tmp.env)
    }

    ## replace the ones in the INLA-namespace
    funcs = ls(tmp.env)
    env = as.environment("package:INLA")
    nfuncs = 0
    for (func in funcs) {
        if (existsFunction(f=func, where = tmp.env)) {
            locked = try(bindingIsLocked(func, env), silent = TRUE)
            if (class(locked) %in% "try-error") {
                ## then this function does not exists in package:INLA,
                ## so we assign it in the globalenv()
                assign(func, get(func, envir = tmp.env), envir = globalenv())
            } else {
                ## otherwise, we change it in the environment and
                ## namespace of INLA.
                if (locked) {
                    unlockBinding(func, env)
                }
                ## 'assignInNamespace' might not be 'allowed' in the future.
                try(assignInNamespace(func, get(func, envir = tmp.env), ns = "INLA", envir = env), silent=TRUE)
                assign(func, get(func, envir = tmp.env), envir = env)
                if (locked) {
                    lockBinding(func, env)
                }
            }
            nfuncs = nfuncs + 1
        }
    }

    cat("Source files in ", dir, ". Loaded ", length(files), " files and replaced ", nfuncs, " functions.\n", sep="")

    if (binaries) {
        inla.setOption("inla.call", paste(bin.path, "/", "inla", sep=""))
        inla.setOption("fmesher.call", paste(bin.path, "/", "fmesher", sep=""))
        cat("Define new values for 'inla.call' and 'fmesher.call': ", bin.path, "/{inla,fmesher}\n", sep="")
    }

    ## hash the models again
    assign("hgid", "hash it again, please!", envir = inla.get.inlaEnv())
    assign("inla.models", NULL, envir = inla.get.inlaEnv())
    cat("Reset stored 'inla.models()' in .inlaEnv\n")

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
    if (length(test) > 1)
        stop("oops: len(test) > 1")

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
        mapping[i, j] = inla.lattice2node(i, j, nrow, ncol)
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

    if (length(irow) == length(icol) && length(irow) > 1) {
        ## this makes it kind of 'vectorize' for two arguments...
        n = length(irow)
        k = 1:n
        return (sapply(k,
                       function(irow, icol, nrow, ncol, k) {
                           return (inla.lattice2node(irow[k], icol[k], nrow, ncol))
                       }, irow = irow, icol = icol, nrow = nrow, ncol = ncol))
    } else {
        return ((irow-1) + (icol-1)*nrow + 1)
    }
}

`inla.node2lattice` = function(node, nrow, ncol)
{
    ## convert from a node-number in the graph to a lattice point
    ## (irow, icol); similar to the GMRFLib_node2lattice()-function in
    ## GMRFLib.  Indices goes from irow=1...nrow, to icol=1...ncol and
    ## node=1.....nrow*ncol.

    stopifnot( nrow > 0 && ncol > 0 )

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

`inla.sparse.matrix.pattern` = function(A, factor=1.0, size=NULL, reordering = NULL,
        binary.pattern = TRUE)
{
    ## Calculate sparse matrix pattern, with optional resolution reduction
    A = inla.as.dgTMatrix(A)
    n = dim(A)
    if (is.null(reordering)) {
        AA = list(i=A@i+1L, j=A@j+1L, values=A@x)
    } else {
        if (is.numeric(reordering)) {
            AA = list(i=reordering[A@i+1L], j=reordering[A@j+1L], values=A@x)
        }
        else if (inla.is.element("reordering", reordering)) {
            AA = list(i=reordering$reordering[A@i+1L], j=reordering$reordering[A@j+1L], values=A@x)
        } else {
            stop("This should not happen.")
        }
    }

    if (is.null(size)) {
        ## Resize by factor:
        size = ceiling(n*factor)
    } else if (length(size)==1L) {
        ## Keep aspect ratio:
        size = c(size, size/n[1]*n[2])
    }
    fac = size/n

    ## duplicated entries will simply add up, so we need to truncate
    M = inla.as.dgTMatrix(sparseMatrix(i=pmin(size[1], ceiling(AA$i*fac[1])), j=pmin(size[2], ceiling(AA$j*fac[2])), x=1,
            dims = size))
    if (binary.pattern)
        M[ M != 0 ] = 1

    return (M)
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

`inla.display.matrix` = function(x, wrap=TRUE, xaxt=FALSE, yaxt=FALSE, col=gray(seq(0, 1, len=256)), ...)
{
    ## display a matrix as an image with correct layout and autoscaling

    x = as.matrix(x)
    n = dim(x)
    y = x
    if (wrap) {
        ii = 1L:n[1L]
        jj = n[1L]:1L
        for(j in 1L:n[2L])
            y[ii, j] = x[jj, j]
    }

    ## use the image.plot-function in package fields; its much better...
    inla.squishplot(c(0, 1), c(0, 1), n[1]/n[2])
    if (inla.require("fields")) {
        image.plot(t(y), col=col, bty="n", xaxt="n", yaxt="n", ...)
    } else {
        warning("Please install package `fields'")
        image(t(y), col=col, bty="n", xaxt="n", yaxt="n", ...)
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







`inla.group` = function(x, n = 25, method = c("cut", "quantile"), idx.only = FALSE)
{
    `inla.group.core` = function(x, n = 25, method = c("cut", "quantile"), idx.only)
    {
        ## group covariates into N groups using method "quantile" or
        ## "cut", i.e., the functions quantile() or cut().  the cut
        ## use 'even length' wheras the 'quantile' use even length
        ## quantiles.

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
            aq = unique(quantile(x, probs = c(0, ppoints(n-1), 1)))
            a = cut(x, breaks = as.numeric(aq), include.lowest=TRUE)
        }
        ## the rest is then the same
        nlev = nlevels(a)
        xx = list()
        for(i in 1:nlev)
            xx[[i]] = list()

        for(i in 1:length(x))
            xx[[as.numeric(a[i])]] = c(unlist(xx[[as.numeric(a[i])]]), x[i])
        values = numeric(nlev)

        ff.local = function(xx) {
            if (length(xx) > 0)
                return (median(xx))
            else
                return (NA)
        }
        if (!idx.only) {
            values = unlist(sapply(xx, ff.local))
            return (as.numeric(values[as.numeric(a)]))
        } else {
            return (as.numeric(a))
        }
    }

    if (missing(x))
        return (NULL)

    if (any(is.na(x))) {
        idx.ok = !is.na(x)
        x[idx.ok] = inla.group.core(x[idx.ok], n, method, idx.only)

        return (x)
    } else {
        return (inla.group.core(x, n, method, idx.only))
    }
}

`inla.group.old` = function(x, n)
{
    ### old version
    cutpoints = seq(range(x)[1], range(x)[2], length.out=(n+1))
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
    ## remove `_' and `.' and space and tabs, and then convert to
    ## lowercase. remove also everything after a `:', so that "name:
    ## a" is "name:".
    family = gsub(":.*$", ":", family)
    family = tolower(gsub("[_. \t]+", "", family))
    return (family)
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
    if (is.matrix(a) || is.data.frame(a) || length(dim(a))==2) {
        return(as.list(as.data.frame(as.matrix(a))))
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
    return(substr(elms, i, i))
}

`inla.2list` = function(x = NULL)
{
    ## convert a vector `x' or string 'x' into the string ``c(x[1],
    ## x[2], x[3])'' so it can be evaluated using
    ## inla.eval()

    if (is.numeric(x)) {
        return (as.character(enquote(as.numeric(x)))[2])
    } else {
        return (as.character(enquote(x))[2])
    }

    ## example:
    ## > x=1:3
    ## > inla.2list(x)
    ## [1] "c(1, 2, 3)"

    if (is.null(x))
        return (NULL)

    if (!is.character(x)) {
        return (paste("c(", paste(x, collapse=","), ")"))
    } else {
        if (length(x) == 0)
            return ("numeric(0)")
        s = inla.paste(as.character(x))
        s = gsub("^[ ]+", "", s)
        s = gsub("[ ]+$", "", s)
        s = gsub("c[ ]*[(][ ]*", "", s)
        s = gsub("[ ]*[)][ ]*", "", s)
        return (paste("c(", gsub("[, ]+", ",", s, ")"), ")", sep=""))
    }
}

`inla.even` = function(n)
{
    return (inla.divisible(n, by=2L))
}

`inla.odd` = function(n)
{
    return (inla.divisible(n, by=-2L))
}

`inla.divisible` = function(n, by=2L)
{
    ### if by>0, return TRUE if `n' is divisible by `by', and if by<0,
    ### return TRUE if `n' is not divisible by `-by'. if by==0, return
    ### TRUE.
    if (by == 0L)
        return (rep(TRUE, length(n)))
    if (by > 0L)
        return ((n%%by) == 0L)
    else
        return ((n%%(-by)) != 0L)
}

`inla.one.of` = function(family, candidates)
{
    if (is.null(candidates) || length(candidates) == 0)
        return (FALSE)

    ## check if family is one of the canidates
    return (any(inla.trim.family(family) == inla.trim.family(candidates)))
}

`inla.idx` = function(idx, n = max(idx), group = rep(1, n), ngroup = max(group),
        replicate = rep(1, n), nrep = max(replicate))
{
    ## this function might be useful to convert from (idx, group, rep)
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


`inla.eval` = function(command,
        envir = parent.frame(),
        enclos = if (is.list(envir) || is.pairlist(envir))
        parent.frame() else baseenv())
{
    return (eval(parse(text=command), envir, enclos))
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


`inla.formula2character` = function(formula)
{
    ## convert a formula to characters without the 500 character
    ## constraint. Remove the ``()'' at the end.
    return (gsub("[(][)]$","", inla.paste(deparse(formula))))
}

`inla.unique.rows` = function(A)
{
    ## return the unique rows and index-list how to map rows of A.

    ## Example:
    ## > A
    ##       [, 1] [, 2]
    ##  [1,]    1    1
    ##  [2,]    2    2
    ##  [3,]    3    3
    ##  [4,]    1    1
    ##  [5,]    2    2
    ##  [6,]    3    3
    ##  [7,]    4    4
    ##  [8,]    4    4
    ##  [9,]    1    1

    ## > inla.unique.rows(A)
    ## $rows
    ##     [, 1] [, 2]
    ## [1,]    1    1
    ## [2,]    2    2
    ## [3,]    3    3
    ## [4,]    4    4
    ##
    ## $idx
    ## [1] 1 2 3 1 2 3 4 4 1


    stopifnot(is.matrix(A))

    ## we will use the builtin function 'duplicated', but since we
    ## have a matrix, we form a vector where each element is the
    ## character string of that column.

    n = dim(A)[1]
    ncol = dim(A)[2]
    a = apply(A, 1, function(x) inla.paste(c(as.character(x)), sep="<|>"))

    ## then we know which columns that are duplicated. the next step
    ## is to make an index-table over unique rows, and if duplicated,
    ## the index to the first equal unique row.

    dup = as.integer(duplicated(a))
    nu = sum(!dup)
    uA = matrix(NA, nu, ncol)
    k = 1
    for(i in 1:n) {
        if (dup[i] == 1) {
            idx = which(a[i] == a[1:(i-1)])[1]
            dup[i] = dup[ idx ]
        } else {
            dup[i] = k
            uA[k, ] = A[i, ]
            k = k+1
        }
    }

    ## we need both the unique ones and the mapping

    return (list(rows = uA, idx = dup))
}

`inla.is.dir` = function(dir)
{
    return (!is.na(file.info(dir)$isdir) && file.info(dir)$isdir)
}

`inla.dirname` = function(path)
{
    if (identical(substr(path, nchar(path), nchar(path)), .Platform$file.sep)) {
        return (substr(path, 1, nchar(path)-1))
    } else {
        return (dirname(path))
    }
}

`inla.affirm.integer` = function(A,...)
{
    is.wholenumber = function(x, tol = .Machine$double.eps*2)
    {
        return (abs(x - round(x)) <= tol)
    }
    is.integer.values = function(A,...)
    {
        return (is.integer(A) || all(is.wholenumber(A,...)))
    }
    if (is.integer.values(A,...)) {
        A = round(A)
        storage.mode(A) = "integer"
    }
    return(A)
}

`inla.affirm.double` = function(A,...)
{
    if (!is(A, "Matrix")) {
        storage.mode(A) = "double"
    }
    return (A)
}

`inla.matrix2list` = function(A, byrow = FALSE)
{
    ## convert a matrix to a list of list, with columns as the first
    ## index (byrow=FALSE) with with rows as the first index
    ## (byrow=TRUE).

    if (byrow) {
        return (inla.matrix2list(t(A)))
    }

    stopifnot(is.matrix(A))
    a = c()
    for(i in 1:nrow(A)) {
        b = list(as.list(A[i,]))
        names(b) = rownames(A)[i]
        a = c(a, b)
    }
    return(a)
}

`inla.dir.create` = function(dir, showWarnings = TRUE, recursive = TRUE, mode = "0777", StopOnError = TRUE)
{
    ## change the default behaviour, so that if DIR exists an error is
    ## NOT flagged, and recursive is default TRUE.
    if (inla.is.dir(dir))
        return (TRUE)

    result = try(dir.create(dir, showWarnings = showWarnings, recursive = recursive, mode = mode))
    if ((inherits(result, "try-error") || !result) && StopOnError) {
        stop(paste("Fail to create directory [", dir, "]", sep=""))
    }

    return (result)
}

`inla.is.element` = function(name, alist)
{
    ## return TRUE if element with name NAME is a member of LIST and
    ## the value is non null and not NA.
    if (any(names(alist) == name)) {
        idx = which(names(alist) == name)
        if (!is.null(alist[[idx]]) && !is.na(alist[[idx]])) {
            return (TRUE)
        } else {
            return (FALSE)
        }
    } else {
        return (FALSE)
    }
}
`inla.get.element` = function(name, alist)
{
    if (inla.is.element(name, alist)) {
        return (alist[[which(names(alist) == name)]])
    } else {
        return (NULL)
    }
}
`inla.require.inherits` = function(x, what, name="Object")
{
    if (!inherits(x, what)) {
        n.what = length(as.list(what))
        stop(paste(name, " must inherit from class ",
                   inla.ifelse(n.what==1,
                               paste("\"", what, "\".", sep=""),
                               paste("\"",
                                     inla.paste(what[1:(n.what-1)],
                                                "\", \""),
                                     "\"",
                                     inla.ifelse(n.what==2, "", ","),
                                     " or \"",
                                     what[n.what],
                                     "\".",
                                     sep="")),
                   sep=""))
    }
    return(invisible())
}


`inla.function2source` = function(the.function,  newline = "<<NEWLINE>>")
{
    ## take a function and return the souce with 'newline' as the
    ## newline
    attributes(the.function) = NULL
    return (paste(deparse(the.function), collapse=newline))
}

`inla.source2function` = function(the.source, newline = "<<NEWLINE>>")
{
    ## take function source, output from inla.function2source(), and
    ## return the function, using 'newline' as newline.
    return (inla.eval(strsplit(the.source, split = newline)[[1]]))
}

`inla.writeLines` = function(filename, lines)
{
    ## write a sequence of lines to file in binary format. need to do
    ## it like this in order to solve the unix->windows encoding
    ## issue, reading and writing this same file on different
    ## platforms.

    fp = file(filename, "wb")
    len = length(lines)
    writeBin(len, fp)

    for(i in 1L:len) {
        nc = nchar(lines[i])
        writeBin(nc, fp)
        writeChar(lines[i], fp, nchars = nchar(lines[i]), eos=NULL)
    }
    close(fp)
}

`inla.readLines` = function(filename)
{
    ## read a sequence of lines to file in binary format. need to do
    ## it like this in order to solve the unix->windows encoding
    ## issue, reading and writing this same file on different
    ## platforms.

    if (!file.exists(filename)) {
        return (NULL)
    }

    fp = file(filename, "rb")
    len = readBin(fp, integer(), n = 1L)
    lines = character(len)

    for(i in 1L:len) {
        nc = readBin(fp, integer(), n = 1L)
        lines[i] = readChar(fp, nc)
    }
    close(fp)

    return(lines)
}

`inla.is.matrix` = function(A)
{
    ## return TRUE if A is a 'matrix' neglecting possible formats
    ## (dense, sparse, etc...)
    if (is.matrix(A) || is(A,  "Matrix")) {
        ## the common cases
        return (TRUE)
    } else {
        ## then its something else
        d = dim(A)
        if (is.null(d)) {
            return (FALSE)
        } else {
            return (length(d) == 2)
        }
    }
}

`inla.dev.new` = function(...)
{
    ## If running in RStudio then don't open a new device,  otherwise,  do.
    dev = getOption("device")
    if (is.character(dev) && inla.strncasecmp(dev, "RStudioGD")) {
        return (NULL)
    } else {
        return (dev.new(...))
    }
}

`inla.factor2matrix` = function(f, sparse=FALSE)
{
    ok = !is.na(f)
    levels = levels(f)
    if (sparse) {
        factor.matrix = sparseMatrix(i=which(ok), j=as.integer(f[ok]), x=1,
                dims=c(length(f), nlevels(f)))
    } else {
        factor.matrix = matrix(0, length(f), nlevels(f))
        for (k in 1:nlevels(f)) {
            factor.matrix[ok, k] = (f[ok]==levels[k])
        }
    }
    colnames(factor.matrix) = levels(f)
    rownames(factor.matrix) = names(f)
    return(factor.matrix)
}

`inla.is.installed` = function(pkg) 
{
    ## return TRUE if PKG is installed and FALSE if not
    return (is.element(pkg, installed.packages()[,1L]))
}

`inla.require` = function(pkg)
{
    ## load PKG if it exists, but be silent. return status
    w = getOption("warn")
    options(warn = -1L)
    value = (inla.is.installed(pkg) && require(pkg, quietly = TRUE, character.only = TRUE))
    options(warn = w)

    return (value)
}
