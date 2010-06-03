
### Functions to deal with `inla.options'

`inla.getOption` = function(option = c("inla.call", "inla.arg", "fmesher.call", "fmesher.arg", "num.threads", "keep",
                                    "working.directory", "silent", "debug",
                                    "cygwin", "ssh.auth.sock", "cygwin.home"))
{
    if (missing(option))
        stop("argument is required.")
    
    ## get default options from the global list `inla.options'
    option = match.arg(option, several.ok = TRUE)
    if (exists("inla.options", envir = .GlobalEnv))
        opt = get("inla.options", envir = .GlobalEnv)
    else
        opt = list()

    if (is.null(opt$inla.call))
        inla.call = inla.call.builtin()
    else if (inla.strcasecmp(opt$inla.call, "remote") || inla.strcasecmp(opt$inla.call, "inla.remote"))
        inla.call = gsub("\\\\", "/", system.file("bin/remote/inla.remote", package="INLA"))
    else
        inla.call = opt$inla.call

    if (is.null(opt$fmesher.call))
        fmesher.call = inla.fmesher.call.builtin()
    else
        fmesher.call = opt$fmesher.call

    all.opt = list(
            inla.call = inla.call,
            fmesher.call = fmesher.call,
            inla.arg = NULL,
            fmesher.arg = NULL,
            num.threads = inla.ifelse(!is.null(opt$num.threads), opt$num.threads, NULL),
            keep = inla.ifelse(!is.null(opt$keep), opt$keep, FALSE),
            working.directory = inla.ifelse(!is.null(opt$working.directory), opt$working.directory, NULL),
            silent = inla.ifelse(!is.null(opt$silent), opt$silent, TRUE),
            debug = inla.ifelse(!is.null(opt$debug), opt$debug, FALSE),
            cygwin = inla.ifelse(!is.null(opt$cygwin), gsub("\\\\", "/", opt$cygwin), "C:/cygwin"),
            cygwin.home = inla.ifelse(!is.null(opt$cygwin.home), opt$cygwin.home, paste("/home/", inla.get.USER(), sep="")),
            ssh.auth.sock = inla.ifelse(!is.null(opt$ssh.auth.sock), opt$ssh.auth.sock,
                    paste("/tmp/ssh-auth-sock-", inla.get.USER(), sep=""))
            )

    res = c()
    if (length(option) > 0) {
        for (i in 1:length(option))
            res = c(res, inla.eval(paste("all.opt$", option[i], sep="")))
    } else {
        res = all.opt
    }
    return (res)
}

`inla.setOption` = function(option = c("inla.call", "inla.arg", "fmesher.call", "fmesher.arg", "num.threads", "keep",
                                    "working.directory", "silent", "debug",
                                    "cygwin", "ssh.auth.sock", "cygwin.home"), value)
{
    option = match.arg(option, several.ok = FALSE)
    if (!exists("inla.options", envir = .GlobalEnv))
        assign("inla.options", list(), envir = .GlobalEnv)
    if (is.character(value)) {
        eval(parse(text=paste("inla.options$", option, "=", shQuote(value), sep="")),
             envir = .GlobalEnv)
    } else {
        eval(parse(text=paste("inla.options$", option, "=", inla.ifelse(is.null(value), "NULL", value), sep="")),
             envir = .GlobalEnv)
    }
}
