## Nothing to export

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
    ## trim leading and trailing whitespaces. there is a function in R.oo called `trim' that do
    ## this, but I don't want INLA to be dependent on R.oo. This function also works the string
    ## is a list of strings.

    string = gsub("^[ \t]+", "", string)
    string = gsub("[ \t]+$", "", string)
    return (string)
}

`inla.namefix` = function(string)
{
    ## must be the same as in iniparser.h
    re = "[$]" 
    re.to = "|S|"
    old.string = inla.trim(string)
    ## special characters, need to do something
    while(TRUE) {
        string = gsub(re, re.to, old.string)
        if (string == old.string) break
        old.string = string
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
    ## cannot call inla.getOption() here as it leads to an infinite recursive call. do this
    ## manually instead.
    if (exists("inla.options", env = inla.get.inlaEnv())) {
        opt = get("inla.options", env = inla.get.inlaEnv())
        mkl = if (!is.null(opt$mkl) && opt$mkl) "mkl." else ""
    } else {
        mkl = ""
    }

    if (inla.os("mac")) {
        fnm = system.file(paste("bin/mac/", inla.os.32or64bit(), "bit/inla.", mkl, "run", sep=""), package="INLA")
    } else if (inla.os("linux")) {
        fnm = system.file(paste("bin/linux/", inla.os.32or64bit(), "bit/inla.", mkl, "run", sep=""), package="INLA")
    } else if (inla.os("windows")) {
        fnm = system.file(paste("bin/windows/", inla.os.32or64bit(), "bit/inla.exe", sep=""), package="INLA")
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
        fnm = system.file(paste("bin/mac/", inla.os.32or64bit(), "bit/fmesher.run", sep=""), package="INLA")
    } else if (inla.os("linux")) {
        fnm = system.file(paste("bin/linux/", inla.os.32or64bit(), "bit/fmesher.run", sep=""), package="INLA")
    } else if (inla.os("windows")) {
        fnm = system.file(paste("bin/windows/", inla.os.32or64bit(), "bit/fmesher.exe", sep=""), package="INLA")
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
                      "danielsimpson", "rieblera", "ariebler", "fuglstad")))
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
    } else if (Sys.getenv("USER") %in% c("fuglstad")) {
        dir.default = "~/build64/rinla/R"
        bin.default = "~/build64/local/bin"
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
        res = try(local({source(fff, local=TRUE)}, envir = tmp.env))
        if (class(res) %in% "try-error")
            warning(paste0("Got an error while sourcing file ",  fff))
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
        inla.setOption("inla.call", path.expand(paste(bin.path, "/", "inla", sep="")))
        inla.setOption("fmesher.call", path.expand(paste(bin.path, "/", "fmesher", sep="")))
        cat("Define new values for 'inla.call' and 'fmesher.call'\n", sep="")
    }

    ## hash the models again
    assign("hgid", "(Undefined)", envir = inla.get.inlaEnv())
    assign("inla.models", NULL, envir = inla.get.inlaEnv())
    cat("Reset stored 'inla.models()' in .inlaEnv\n")

    return (invisible())
}

`inla.remove` = function(name, from)
{
    ## remove NAME FROM. Works for both lists and data.frames

    if (FALSE) {
        ## OLD CODE
        if (is.list(from) || is.data.frame(from))
            for(nm in name)
                if (length(grep(nm, names(from))) > 0)
                    inla.eval(paste("from$", nm, "=NULL", sep=""))
    } else {
        ## NEW CODE
        if (is.list(from) || is.data.frame(from))
            for(nm in name)
                from[which(nm == names(from))] = NULL
    }

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



`inla.scale` = function(x)
{
    ## didn't knew about this one...
    return (scale(x))
}

`inla.trim.family` = function(family)
{
    ## remove `_' and `.' and space and tabs, and then convert to
    ## lowercase. remove also everything after a `:', so that "name:
    ## a" is "name:".
    family = gsub(":.*$", ":", family)
    family = tolower(gsub("[_ \t.]+", "", family))
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
    if (inla.is.dir(dir)) {
        return (dir)
    }

    result = try(dir.create(dir, showWarnings = showWarnings, recursive = recursive, mode = mode))
    if ((inherits(result, "try-error") || !result)) {
        if (StopOnError) {
            stop(paste("Fail to create directory [", dir, "]. Stop.", sep=""))
        }
        result = NULL
    }

    return (result)
}

`inla.is.element` = function(name, alist)
{
    ## return TRUE if element with name NAME is a member of LIST and
    ## the value is non null and not NA.
    if (any(names(alist) == name)) {
        idx = which(names(alist) == name)
        if (length(alist[[idx]]) > 0 && !is.null(alist[[idx]]) && !is.na(alist[[idx]])) {
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

`inla.require` = function(pkg, ...)
{
    ## follow the `new' standard...
    return (requireNamespace(pkg, quietly = TRUE, ...))
}

`inla.require.old` = function(pkg)
{
    ## load PKG if it exists, but be silent. return status
    w = getOption("warn")
    options(warn = -1L)
    value = (inla.is.installed(pkg) && require(pkg, quietly = TRUE, character.only = TRUE))
    options(warn = w)

    return (value)
}

`inla.inlaprogram.has.crashed` = function()
{
    stop("The inla-program exited with an error. Unless you interupted it yourself, please rerun with verbose=TRUE and check the output carefully.\n  If this does not help, please contact the developers at <help@r-inla.org>.")
}

`inla.eval.dots` = function(..., stop.if.no.name = TRUE, allowed.names = NULL) {
    ## evaluate named argument in the parent frame. allowed.names can
    ## be a list of allowed names, or if NULL then all names are
    ## allowed. this function will give an error if one argument has
    ## no name and stop.of.no.name is TRUE
    if (!length(list(...))) {
        return (invisible())
    }
    dots = lapply(match.call(), eval, envir = parent.frame())[-1L]
    for(i in seq_along(dots)) {
        nm = names(dots)[i]
        if (!is.null(nm)) {
            if (!is.null(allowed.names)) {
                if (!(nm %in% allowed.names)) {
                    stop(paste("This argument is not allowed:", nm))
                }
            }
            assign(nm, dots[[i]], envir = parent.frame())
        } else {
            if (stop.if.no.name) {
                stop(paste("The", i, "th argument has no name."))
            }
        }
    }
    return (invisible())
}


`match.arg.vector` <- function(arg=NULL,
                               choices,
                               length=NULL) {
    ## Like match.arg, but for a vector of options 'arg'
    if (is.null(length)) {
        length = inla.ifelse(is.null(arg), 1, length(arg))
    }
    if (is.null(arg)) {
        arg = match.arg(arg, choices)
    } else {
        for (k in seq_along(arg)) {
            arg[k] = match.arg(arg[k], choices)
        }
    }
    if (length(arg) < length) {
        arg = c(arg, rep(arg, length-length(arg)))
    } else if (length(arg) > length) {
        stop('Option list too long.')
    }
    return(arg)
}
`inla.get.var` = function(var, data = NULL)
{
    if (is.null(data)) {
        return (get(var))
    } else {
        return (get(var, envir = as.environment(data)))
    }
}
`inla.ginv` = function(x, tol = sqrt(.Machine$double.eps), rankdef = NULL)
{
    ## from MASS:::ginv, but with added option 'rankdef'.
    if (!is.matrix(x)) {
        x <- as.matrix(x)
    }
    if (length(dim(x)) > 2 || !(is.numeric(x) || is.complex(x))) {
        stop("'x' must be a numeric or complex matrix")
    }

    xsvd <- svd(x)
    if (is.complex(x)) {
        xsvd$u <- Conj(xsvd$u)
    }

    if (is.null(rankdef) || rankdef == 0) {
        Positive <- xsvd$d > max(tol * xsvd$d[1], 0)
    } else {
        n = length(xsvd$d)
        stopifnot(rankdef >= 1 && rankdef <= n)
        Positive <- c(rep(TRUE, n - rankdef), rep(FALSE, rankdef))
    }

    if (all(Positive)) {
        xsvd$v %*% (1/xsvd$d * t(xsvd$u))
    } else if (!any(Positive)) {
        array(0, dim(x)[2:1])
    } else {
        xsvd$v[, Positive, drop = FALSE] %*% ((1/xsvd$d[Positive]) *
                                              t(xsvd$u[, Positive, drop = FALSE]))
    }
}
`inla.gdet` = function(x, tol = sqrt(.Machine$double.eps), rankdef = NULL, log=TRUE)
{
    x = as.matrix(x)
    lambda = eigen(x, only.values = TRUE)$values
    if (is.null(rankdef) || rankdef == 0) {
        non.zero = (lambda > max(tol * max(lambda), 0))
        lambda = lambda[non.zero]
    } else {
        if (rankdef > 0) {
            lambda = sort(lambda)
            lambda = lambda[-(1:rankdef)]
        }
    }

    if (log) {
        return (sum(log(lambda)))
    } else {
        return (prod(lambda))
    }
}
## nice to have these around
`inla.rw1` = function(n, ...)
{
    inla.rw(n, order=1L, ...)
}
`inla.rw2` = function(n, ...)
{
    inla.rw(n, order=2L, ...)
}
`inla.rw` = function(n, order=1L, sparse=TRUE, scale.model = FALSE)
{
    if (scale.model) {
        Q = inla.rw(n, order = order,  sparse=sparse,  scale.model=FALSE)
        fac = exp(mean(log(diag(inla.ginv(as.matrix(Q), rankdef = order)))))
        Q = fac * Q
        if (sparse) Q = inla.as.sparse(Q)
        return (Q)
    } else {
        U = diff(diag(n), diff=order)
        if (sparse) {
            return (inla.as.sparse(t(U) %*% U))
        } else {
            return (t(U) %*% U)
        }
    }
}

##
`inla.mclapply` = function(..., mc.cores = NULL, parallel = TRUE)
{
    if (parallel && !inla.os("windows")) {
        if (is.null(mc.cores)) {
            mc.cores = inla.getOption("num.threads")
        }
        return (parallel::mclapply(..., mc.cores = mc.cores))
    } else {
        return (lapply(...))
    }
}
`inla.cmpfun` = function(fun, options = list(optimize = 3L))
{
    if (inla.require("compiler")) {
        return (compiler::cmpfun(fun, options = options))
    } else {
        return(fun)
    }
}

`inla.sn.reparam` = function(moments, param)
{
    stopifnot((!missing(moments) && missing(param)) ||
              (missing(moments) && !missing(param)))

    if (!missing(param)) {
        if (is.list(param)) {
            xi = param$xi
            omega = param$omega
            alpha = param$alpha
        } else {
            xi = param[1]
            omega = param[2]
            alpha = param[3]
        }
        delta = alpha/sqrt(1 + alpha^2)
        mean = xi + omega * delta * sqrt(2/pi)
        variance = omega^2 * (1 - 2*delta^2/pi)
        skewness = (4-pi)/2 * (delta * sqrt(2/pi))^3 / (1-2*delta^2/pi)^(3/2)

        return (list(mean = mean, variance = variance, skewness = skewness))
    } else {
        if (is.list(moments)) {
            mean = moments$mean
            variance = moments$variance
            skewness = moments$skewness
        } else {
            mean = moments[1]
            variance = moments[2]
            skewness = moments[3]
        }
        delta = sqrt(pi/2 *abs(skewness)^(2/3) /
                (abs(skewness)^(2/3) + ((4-pi)/2)^(2/3)))
        delta = delta * sign(skewness)
        alpha = delta/sqrt(1-delta^2)
        omega = sqrt(variance / (1-2*delta^2/pi))
        xi = mean - omega * delta * sqrt(2/pi)

        return (list(xi = xi, omega = omega, alpha = alpha))
    }
}

`inla.runjags2dataframe` = function(runjags.object)
{
    ## convert from runjags-output to a data.frame
    stopifnot(inla.require("runjags"))
    return (as.data.frame(runjags::combine.mcmc(runjags.object, collapse.chains=TRUE)))
}

`inla.check.location` = function(loc, term, model, section = "latent")
{
    if (is.null(loc) || length(loc) <= 1) {
        return (invisible())
    }
    lim = inla.model.properties(model, section)$min.diff
    if (is.null(lim))
        return (invisible())

    min.diff = min(diff(sort(loc)))/diff(range(loc))
    if (min.diff < lim) {
        stop(paste(sep = "",
                   "Locations are too close for f(",
                   term, ", model=\"",
                   model, "\", ...): ",
                   " min(diff(sort(x)))/diff(range(x)) = ", 
                   format(min.diff, scientific=TRUE, digits=4),
                   " < ", format(lim, scientific=TRUE, digits=4), "\n",
                   "  You can fix this by some kind of binning, see ?inla.group", "\n",
                   "  If you want/need to bypass this check at your own risk, do", "\n",
                   "\t> m = get(\"inla.models\", INLA:::inla.get.inlaEnv())\n",
                   "\t> m$", section, "$", model, "$min.diff = NULL\n",
                   "\t> assign(\"inla.models\", m, INLA:::inla.get.inlaEnv())"))

    }
    return (invisible())
}

`inla.dynload.workaround` = function()
{
    if (inla.os("linux")) {
        ## setup the static builds instead
        mkl = inla.getOption("mkl")
        if (is.null(mkl)) mkl = FALSE
        d = dirname(inla.call.builtin())
        inla.setOption(inla.call = paste(d, "/inla.", if (mkl) "mkl." else "", "static", sep=""))
        inla.setOption(fmesher.call = paste(d,"/fmesher.static", sep=""))
    }
    return (invisible())
}

`inla.matern.cf` = function(dist, range = 1, nu = 0.5) 
{
    ## the matern correlation function, with parameter 'range' as defined for the SPDE models.
    ## this function can vectorize over 'dist'
    kappa = sqrt(8.0 * nu) / range  ## this the definition used
    d = kappa * dist
    res = numeric(length(dist))
    is.zero = (dist == 0.0)
    d = d[!is.zero]
    res[is.zero] = 1.0
    res[!is.zero] = 1.0 / 2.0^(nu - 1.0) / gamma(nu) * d^nu * besselK(d, nu)
    return (res)
}

