### Functions to deal with options in INLA

`inla.getOption` = function(option = c("inla.call", "inla.arg", "fmesher.call", "fmesher.arg", "num.threads", "keep",
                                    "working.directory", "silent", "debug", "internal.binary.mode", "internal.experimental.mode", 
                                    "cygwin", "ssh.auth.sock", "cygwin.home",
                                    "enable.inla.argument.weights",
                                    "show.warning.graph.file",
                                    "show.warning.control.data"))
{
    if (missing(option))
        stop("argument is required.")

    envir = inla.get.inlaEnv()

    option = match.arg(option, several.ok = TRUE)
    if (exists("inla.options", envir = envir))
        opt = get("inla.options", envir = envir)
    else
        opt = list()

    if (is.null(opt$inla.call))
        inla.call = inla.call.builtin()
    else if (inla.strcasecmp(opt$inla.call, "remote") || inla.strcasecmp(opt$inla.call, "inla.remote"))
        inla.call = gsub("\\\\", "/", system.file("bin/remote/inla.remote", package="INLA"))
    else
        inla.call = opt$inla.call

    all.opt = list(
            inla.call = inla.call,
            fmesher.call = inla.ifelse(is.null(opt$fmesher.call), inla.fmesher.call.builtin(), opt$fmesher.call),
            inla.arg = NULL,
            fmesher.arg = inla.ifelse(is.null(opt$fmesher.arg), "", opt$fmesher.arg),
            num.threads = inla.ifelse(!is.null(opt$num.threads), opt$num.threads, NULL),
            keep = inla.ifelse(!is.null(opt$keep), opt$keep, FALSE),
            working.directory = inla.ifelse(!is.null(opt$working.directory), opt$working.directory, NULL),
            silent = inla.ifelse(!is.null(opt$silent), opt$silent, TRUE),
            debug = inla.ifelse(!is.null(opt$debug), opt$debug, FALSE),
            internal.binary.mode = inla.ifelse(!is.null(opt$internal.binary.mode),  opt$internal.binary.mode, TRUE), 
            internal.experimental.mode = inla.ifelse(!is.null(opt$internal.experimental.mode),  opt$internal.experimental.mode, FALSE), 
            cygwin = inla.ifelse(!is.null(opt$cygwin), gsub("\\\\", "/", opt$cygwin), "C:/cygwin"),
            cygwin.home = inla.ifelse(!is.null(opt$cygwin.home), opt$cygwin.home, paste("/home/", inla.get.USER(), sep="")),
            ssh.auth.sock = inla.ifelse(!is.null(opt$ssh.auth.sock), opt$ssh.auth.sock,
                    paste("/tmp/ssh-auth-sock-", inla.get.USER(), sep="")),
            enable.inla.argument.weights = inla.ifelse(is.null(opt$enable.inla.argument.weights), FALSE, opt$enable.inla.argument.weights), 
            show.warning.graph.file = inla.ifelse(is.null(opt$show.warning.graph.file), TRUE, opt$show.warning.graph.file), 
            show.warning.control.data = inla.ifelse(is.null(opt$show.warning.control.data), TRUE, opt$show.warning.control.data)
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

`inla.setOption` = function(...)
{
    ## now supports more formats, and also the common one
    ##     inla.setOption("keep", TRUE)
    ## and
    ##     inla.setOption(keep=TRUE)
    ## and
    ##     inla.setOption(keep=TRUE, num.threads=10)

    `inla.setOption.core` = function(option = c("inla.call",
                                             "inla.arg",
                                             "fmesher.call",
                                             "fmesher.arg",
                                             "num.threads",
                                             "keep",
                                             "working.directory",
                                             "silent",
                                             "debug",
                                             "internal.binary.mode",
                                             "internal.experimental.mode", 
                                             "cygwin",
                                             "ssh.auth.sock",
                                             "cygwin.home",
                                             "enable.inla.argument.weights",
                                             "show.warning.graph.file", 
                                             "show.warning.control.data"), value)
    {
        envir = inla.get.inlaEnv()

        option = match.arg(option, several.ok = FALSE)
        if (!exists("inla.options", envir = envir))
            assign("inla.options", list(), envir = envir)
        if (is.character(value)) {
            eval(parse(text=paste("inla.options$", option, "=", shQuote(value), sep="")),
                 envir = envir)
        } else {
            eval(parse(text=paste("inla.options$", option, "=", inla.ifelse(is.null(value), "NULL", value), sep="")),
                 envir = envir)
        }
        return (invisible())
    }

    called = list(...)
    len = length(names(called))
    if (len > 0L) {
        for(i in 1L:len) {
            do.call(inla.setOption.core, args = list(names(called)[i], called[[i]]))
        }
    } else {
        inla.setOption.core(...)
    }
    return (invisible())
}
