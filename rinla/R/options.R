## Export: inla.setOption inla.getOption

##! \name{inla.option}
##! \alias{inla.option}
##! \alias{inla.options}
##! \alias{inla.setOption}
##! \alias{inla.getOption}
##! 
##! \title{Set and get global options for INLA}
##! 
##! \description{Set and get global options for INLA}
##! \usage{
##! inla.setOption(...)
##! inla.getOption(option)
##! }
##! 
##! \arguments{
##! 
##!   \item{...}{Option and value,  like \code{option=value} or \code{option, value}; see the Examples}
##!   \item{option}{The option to get. If \code{option = NULL} then
##!     \code{inla.getOption} then \code{inla.getOption} will return a named list of
##!     current values, otherwise, \code{option} must be one of
##! 
##!     inla.call: The path to the inla-program.
##! 
##!     inla.arg: Additional arguments to \code{inla.call}
##! 
##!     fmesher.call: The path to the fmesher-program
##!
##!     fmesher.arg: Additional arguments to \code{fmesher.call}
##!     
##!     num.threads: Number of threads to use.
##!
##!     blas.num.threads: Number of threads to use for openblas and mklblas (see \code{inla} for details)
##!     
##!     smtp: Sparse matrix library to use,  one of \code{band}, \code{taucs} (\code{default}) or \code{pardiso}
##!
##!     mkl: Use binaries buildt with Intel MKL?  (If possible)
##!
##!     pardiso.license: The full path to the PARDISO license file
##!     
##!     keep: Keep temporary files?
##! 
##!     working.directory: The name of the working directory.
##!
##!     silent: Run the inla-program in a silent mode?
##! 
##!     debug : Run the inla-program in a debug mode?
##!
##!     internal.binary.mode : if \code{FALSE} the (some) output are in ascii format instead of binary format.
##!                            Using this option,  then \code{inla.collect.results} will fail (Expert mode)
##!
##!     internal.experimental.mode :  Expert option
##!
##!     cygwin : The home of the Cygwin installation (default "C:/cygwin") [Remote computing for Windows only]
##! 
##!     ssh.auth.sock: The ssh bind-adress (value of $SSH_AUTH_SOCK int the
##!     Cygwin-shell). [Remote computing for Windows only]
##!
##!     enable.inla.argument.weights : if \code{TRUE} the \code{inla} accepts argument \code{weights} 
##!
##!     show.warning.graph.file : Give a warning for using the obsolete argument
##!                               \code{graph.file} instead of \code{graph} 
##!
##!     scale.model.default : The default value of argument \code{scale.model} which
##!                           optionally scale intrinisic models to have generalized
##!                           unit average variance 
##!
##! The options are stored in the variable \code{inla.options} in the
##!     \code{.GlobalEnv}-environment.
##!   }
##! }
##! 
##! \author{Havard Rue \email{hrue@r-inla.org}}
##! 
##! \examples{
##! ## set number of threads
##! inla.setOption("num.threads", 2)
##! ## alternative format
##! inla.setOption(num.threads=2)
##! ## check it
##! inla.getOption("num.threads")
##!}


`inla.getOption` = function(
    option = c("inla.call",
        "inla.arg",
        "fmesher.call",
        "fmesher.arg",
        "num.threads",
        "blas.num.threads",
        "smtp", 
        "mkl", 
        "pardiso.license", 
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
        "scale.model.default"))
{
    default.opt = list(
        inla.call = inla.call.builtin(), 
        fmesher.call = inla.fmesher.call.builtin(), 
        inla.arg = NULL,
        fmesher.arg = "", 
        num.threads = parallel::detectCores(), 
        blas.num.threads = 1L, 
        smtp = "default", 
        mkl = FALSE, 
        pardiso.license = NULL, 
        keep = FALSE, 
        working.directory = NULL, 
        silent = TRUE, 
        debug = FALSE, 
        internal.binary.mode = TRUE, 
        internal.experimental.mode = FALSE, 
        cygwin = "C:/cygwin",
        cygwin.home = paste("/home/", inla.get.USER(), sep=""), 
        ssh.auth.sock = paste("/tmp/ssh-auth-sock-", inla.get.USER(), sep=""),
        enable.inla.argument.weights = FALSE, 
        show.warning.graph.file = TRUE, 
        scale.model.default = FALSE
    )

    ## with no argument, return a named list of current values
    if (missing(option)) {
        opt.names = names(default.opt)
        option = opt.names
    } else {
        opt.names = NULL
    }

    envir = inla.get.inlaEnv()
    option = match.arg(option, several.ok = TRUE)
    if (exists("inla.options", envir = envir)) {
        opt = get("inla.options", envir = envir)
    } else {
        opt = list()
    }

    if (is.null(opt$inla.call)) {
        inla.call = inla.call.builtin()
    } else if (inla.strcasecmp(opt$inla.call, "remote") ||
               inla.strcasecmp(opt$inla.call, "inla.remote")) {
        inla.call = gsub("\\\\", "/", system.file("bin/remote/inla.remote", package="INLA"))
    } else {
        inla.call = opt$inla.call
    }

    if (is.null(opt$fmesher.call)) {
        fmesher.call = inla.fmesher.call.builtin()
    } else {
        fmesher.call = opt$fmesher.call
    }

    res = c()
    for (i in 1:length(option)) {
        if (inla.is.element(option[i], opt)) {
            val = list(inla.get.element(option[i], opt))
        } else {
            val = list(inla.get.element(option[i], default.opt))
        }
        if (!is.null(opt.names)) {
            names(val) = opt.names[i]
        }
        res = c(res, val)
    }

    if (is.null(opt.names)) {
        res = unlist(res)
    }

    return (res)
}

`inla.setOption` = function(...)
{
    ## supports formats:
    ##     inla.setOption("keep", TRUE)
    ## and
    ##     inla.setOption(keep=TRUE)
    ## and
    ##     inla.setOption(keep=TRUE, num.threads=10)

    `inla.setOption.core` = function(
        option = c("inla.call",
            "inla.arg",
            "fmesher.call",
            "fmesher.arg",
            "num.threads",
            "blas.num.threads",
            "smtp",
            "mkl", 
            "pardiso.license", 
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
            "scale.model.default"), value)
    {
        envir = inla.get.inlaEnv()

        option = match.arg(option, several.ok = FALSE)
        if (!exists("inla.options", envir = envir)) {
            assign("inla.options", list(), envir = envir)
        }
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
