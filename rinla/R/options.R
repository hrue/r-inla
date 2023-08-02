#' Set and get global options for INLA
#' 
#' Set and get global options for INLA
#' 
#' 
#' @aliases inla.option inla.options inla.setOption inla.getOption
#' @param ... Option and value, like `option=value` or `option,
#' value`; see the Examples
#' @param option The option to get. If `option = NULL` then
#' `inla.getOption` then `inla.getOption` will return a named list of
#' current values, otherwise, `option` must be one of
#' \describe{
#' \item{inla.call}{The path to the inla-program.}
#' 
#' \item{inla.arg}{Additional arguments to `inla.call`}
#' 
#' \item{fmesher.call}{The path to the fmesher-program}
#' 
#' \item{fmesher.arg}{Additional arguments to `fmesher.call`}
#' 
#' \item{num.threads}{Character string with the number of threads to use as
#' `A:B`, see `?inla`}
#' 
#' \item{smtp}{Sparse matrix library to use, one of `band`, `taucs`
#' (`default`) or `pardiso`}
#' 
#' \item{safe}{Run in safe-mode (ie try to automatically fix convergence errors)
#' (default `TRUE`)}
#' 
#' \item{pardiso.license}{The full path to the PARDISO license file or a
#' newline-separated string with license key(s)}
#' 
#' \item{keep}{Keep temporary files?}
#' 
#' \item{verbose}{Verbose output?}
#' 
#' \item{save.memory}{Save memory at the cost of (minor) accuracy and computing time?}
#' 
#' \item{working.directory}{The name of the working directory.}
#' 
#' \item{silent}{Run the inla-program in a silent mode?}
#' 
#' \item{debug}{Run the inla-program in a debug mode?}
#' 
#' \item{cygwin}{The home of the Cygwin installation (default "C:/cygwin") (Remote
#' computing for Windows only) (No longer in use!)}
#' 
#' \item{ssh.auth.sock}{The ssh bind-adress (value of $SSH_AUTH_SOCK int the
#' Cygwin-shell). (Remote computing for Windows only)}
#' 
#' \item{show.warning.graph.file}{Give a warning for using the obsolete argument
#' `graph.file` instead of `graph`}
#' 
#' \item{scale.model.default}{The default value of argument `scale.model` which
#' optionally scale intrinisic models to have generalized unit average variance}
#' 
#' \item{short.summary}{Use a less verbose output for `summary`. Useful for
#' Markdown documents.}
#' 
#' \item{inla.timeout}{The timeout limit, in whole seconds, for calls to the inla
#' binary. Default is 0, meaning no timeout limit.  Set to a positive integer
#' to terminate inla calls if they run to long.  Fractional seconds are rounded
#' up to the nearest integer. This feature is EXPERIMENTAL and might change at
#' a later stage.}
#' 
#' \item{fmesher.timeout}{The timeout limit, in whole seconds, for calls to the
#' fmesher binary. Default is 0, meaning no timeout limit.  Set to a positive
#' integer to terminate fmesher calls that may enter infinite loops due to
#' special geometry regularity. Fractional seconds are rounded up to the
#' nearest integer.}
#' 
#' \item{inla.mode}{Which mode to use in INLA? Default is `"compact"`. Other
#' options are `"classic"` and `"twostage"`.}
#' 
#' \item{fmesher.evolution}{Control use of fmesher methods during the transition
#' to a separate fmesher package. Levels of
#' `fmesher.evolution`:
#' \describe{
#' \item{`1L`}{
#' (current default) uses the intermediate `fm_*` methods in `fmesher` that were already
#' available via `inlabru` from 2.8.0.}
#' \item{`2L`}{
#' uses the full range of `fmesher` package methods.}
#' }
#' Further levels may be added as the package development progresses.}
#' 
#' \item{fmesher.evolution.warn}{logical; whether to show warnings about deprecated
#' use of legacy INLA methods with fmesher package replacements. When `TRUE`,
#' shows deprecation messages for many CRS and mesh
#' related methods, pointing to their `fm_*` replacements. Default is currently `FALSE`.}
#' }
#' @author Havard Rue \email{hrue@@r-inla.org}
#' @examples
#' 
#'  ## set number of threads
#'  inla.setOption("num.threads", "4:1")
#'  ## alternative format
#'  inla.setOption(num.threads="4:1")
#'  ## check it
#'  inla.getOption("num.threads")
#' 
#' @name inla.option
#' @rdname options
NULL

`inla.getOption.default` <- function() {
    ## this function is not exported. it need to be separate, to avoid infinite recursion
    return(
        list(
            inla.arg = NULL,
            fmesher.arg = "",
            num.threads = paste0(parallel::detectCores(all.tests = TRUE, logical = FALSE), ":1"),
            smtp = "default",
            safe = TRUE, 
            pardiso.license = NULL,
            keep = FALSE,
            verbose = FALSE,
            save.memory = FALSE,
            working.directory = NULL,
            silent = TRUE,
            debug = FALSE,
            show.warning.graph.file = TRUE,
            scale.model.default = FALSE,
            short.summary = FALSE,
            inla.timeout = 0, 
            fmesher.timeout = 0,
            inla.mode = "compact",
            fmesher.evolution = 2L,
            fmesher.evolution.warn = TRUE
        )
    )
}

#' @rdname options
#' @export
`inla.getOption` <- function(
                             option = c(
                                 "inla.call",
                                 "inla.arg",
                                 "fmesher.call",
                                 "fmesher.arg",
                                 "num.threads",
                                 "smtp",
                                 "safe", 
                                 "pardiso.license",
                                 "keep",
                                 "verbose",
                                 "save.memory",
                                 "working.directory",
                                 "silent",
                                 "debug",
                                 "show.warning.graph.file",
                                 "scale.model.default",
                                 "short.summary",
                                 "inla.timeout", 
                                 "fmesher.timeout",
                                 "inla.mode",
                                 "fmesher.evolution",
                                 "fmesher.evolution.warn"
                             )) {
    ## we 'inla.call' and 'fmesher.call' separately to avoid infinite recursion
    default.opt <- inla.getOption.default()
    default.opt$inla.call <- inla.call.builtin()
    default.opt$fmesher.call <- inla.fmesher.call.builtin()

    ## with no argument, return a named list of current values
    if (missing(option)) {
        opt.names <- names(default.opt)
        option <- opt.names
    } else {
        opt.names <- NULL
    }

    envir <- inla.get.inlaEnv()
    option <- match.arg(option, several.ok = TRUE)
    if (exists("inla.options", envir = envir)) {
        opt <- get("inla.options", envir = envir)
    } else {
        opt <- list()
    }

    if (is.null(opt$inla.call)) {
        inla.call <- inla.call.builtin()
    } else if (inla.strcasecmp(opt$inla.call, "remote") ||
        inla.strcasecmp(opt$inla.call, "inla.remote")) {
        inla.call <- gsub("\\\\", "/", system.file("bin/remote/inla.remote", package = "INLA"))
    } else {
        inla.call <- opt$inla.call
    }

    if (is.null(opt$fmesher.call)) {
        fmesher.call <- inla.fmesher.call.builtin()
    } else {
        fmesher.call <- opt$fmesher.call
    }

    res <- c()
    for (i in 1:length(option)) {
        if (inla.is.element(option[i], opt)) {
            val <- list(inla.get.element(option[i], opt))
        } else {
            val <- list(inla.get.element(option[i], default.opt))
        }
        if (!is.null(opt.names)) {
            names(val) <- opt.names[i]
        }
        res <- c(res, val)
    }

    if (is.null(opt.names)) {
        res <- unlist(res)
    }

    return(res)
}

#' @rdname options
#' @export
`inla.setOption` <- function(...) {
    ## supports formats:
    ## inla.setOption("keep", TRUE)
    ## and
    ## inla.setOption(keep=TRUE)
    ## and
    ## inla.setOption(keep=TRUE, num.threads=10)

    `inla.setOption.core` <- function(
                                      option = c(
                                          "inla.call",
                                          "inla.arg",
                                          "fmesher.call",
                                          "fmesher.arg",
                                          "num.threads",
                                          "smtp",
                                          "safe", 
                                          "pardiso.license",
                                          "keep",
                                          "verbose",
                                          "save.memory",
                                          "working.directory",
                                          "silent",
                                          "debug",
                                          "show.warning.graph.file",
                                          "scale.model.default",
                                          "short.summary",
                                          "inla.timeout", 
                                          "fmesher.timeout",
                                          "inla.mode",
                                          "fmesher.evolution",
                                          "fmesher.evolution.warn"
                                      ), value) {
        envir <- inla.get.inlaEnv()
        option <- match.arg(option, several.ok = FALSE)
        if (!exists("inla.options", envir = envir)) {
            assign("inla.options", list(), envir = envir)
        }
        if (is.character(value)) {
            eval(parse(text = paste("inla.options$", option, "=", shQuote(value), sep = "")),
                envir = envir
            )
        } else {
            eval(parse(text = paste("inla.options$", option, "=", inla.ifelse(is.null(value), "NULL", value), sep = "")),
                envir = envir
            )
        }
        return(invisible())
    }

    called <- list(...)
    len <- length(names(called))
    if (len > 0L) {
        for (i in 1L:len) {
            do.call(inla.setOption.core, args = list(names(called)[i], called[[i]]))
        }
    } else {
        inla.setOption.core(...)
    }

    ## add checks that nothing very wrong is set
    dummy <- match.arg(inla.getOption("inla.mode"),
                       c("compact", "classic", "twostage", "experimental"),
                       several.ok = FALSE)
    if (dummy == "experimental") {
        inla.setOption(inla.mode = "compact")
    }

    return(invisible())
}
