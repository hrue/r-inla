#' Set and get global options for INLA
#' 
#' Set and get global options for INLA
#' 
#' 
#' @aliases inla.option inla.options inla.setOption inla.getOption
#' @param ... Option and value, like `option=value` or `option,
#' value`; see the Examples
#' @param option The option to get. If `option = NULL` then
#' `inla.getOption` will return a named list of
#' current values, otherwise, `option` must be one (or a vector of several) of
#' \describe{
#' \item{inla.call}{The path to the inla-program.}
#' 
#' \item{inla.arg}{Additional arguments to `inla.call`}
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
#' \item{keep}{Keep temporary files?}
#' 
#' \item{verbose}{Verbose output?}
#' 
#' \item{save.memory}{Save memory at the cost of (minor) accuracy and computing time?}
#' 
#' \item{internal.opt}{Do internal online optimisations or not}
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
#' \item{inla.mode}{Which mode to use in INLA? Default is `"compact"`. Other
#' options are `"classic"` and `"twostage"`.}
#' 
#' \item{malloc.lib}{Which malloc library to use: `"je"`, `"tc"`, `"mi"`, `"compiler"` or `"default"`.
#' Option `"compiler"` use the compiler's implementation. The library
#' is loaded using `LD_PRELOAD` and similar functionality. Loosely, `jemalloc` is from Facebook,
#' `tcmalloc` is from Google and `mimalloc` is from Microsoft. This option is not available for
#' Windows and not all options might be available for every arch. 
#' If `malloc.lib` is a complete path to an external library, that file will be used
#' instead of one of the supported ones.}
#'
#' \item{fmesher.evolution}{Control use of fmesher methods during the transition
#' to a separate fmesher package. Only supported value is `2L`, for using the
#' full range of `fmesher` package methods. A higher value may be used in the future
#' to force use of newer methods, but this is not currently implemented.
#' Default is `2L`.}
#' 
#' \item{fmesher.evolution.warn}{logical; whether to show warnings about deprecated
#' use of legacy INLA methods with fmesher package replacements. When `TRUE`,
#' shows deprecation messages for many CRS and mesh
#' related methods, pointing to their `fm_*` replacements. Default
#' since July 2025 is `TRUE`.}
#' 
#' \item{fmesher.evolution.verbosity}{logical or character; at what minimum
#' severity to show warnings about deprecated
#' use of legacy INLA methods with fmesher package replacements.
#' When set to "default" (default), "soft", "warn", or "stop", indicates the
#' minimum warning level used when `fmesher.evolution.warn` is `TRUE`.}
#' 
#' \item{INLAjoint.features}{logical Do not use. By purpose left undocumented}
#'
#' \item{numa}{logical Enable NUMA features (Linux only)}
#'
#' \item{disable.values.warning}{logical Disable warning about missing 'values'}
#'}
#' @param list.out logical; If `TRUE`, then `inla.getOption` will return a named
#' list of values, even if only one option is requested.  Default is `NULL`,
#' which means that a named list will be returned if more than one option is
#' requested, and a single value will be returned if only one option is
#' requested. Using `list.out = FALSE` will return a plain vector of values,
#' without names, even if more than one option is requested; the data type may
#' depend on which options are requested, and may be unpredictable, so use only
#' if you know what you're doing.
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
            num.threads = paste0(max(1, min(16, parallel::detectCores(all.tests = TRUE, logical = FALSE))), ":1"),
            smtp = "default",
            safe = TRUE, 
            keep = FALSE,
            verbose = FALSE,
            save.memory = FALSE,
            internal.opt = TRUE,
            working.directory = NULL,
            silent = TRUE,
            debug = FALSE,
            show.warning.graph.file = TRUE,
            scale.model.default = FALSE,
            short.summary = FALSE,
            inla.timeout = 0, 
            inla.mode = "compact",
            malloc.lib = "mi", 
            fmesher.evolution = 2L,
            fmesher.evolution.warn = TRUE,
            fmesher.evolution.verbosity = "default",
            INLAjoint.features = FALSE,
            numa = FALSE,
            disable.values.warning = FALSE
        )
    )
}

`inla.validOptions` <- function() {
    c("inla.call", names(inla.getOption.default()))
} 

# Non-exported function, used internally in place of `match.arg` to match
# one or more INLA options, with a custom warning message for unknown or
# ambiguous options.
# This allows for partial matching of options while providing informative
# feedback to the user, with logic similar to match.arg.
`inla.matchOption` <- function(option, valid.options = NULL, several.ok = FALSE) {
    if ((length(option) > 1L) && !several.ok) {
        stop("Only one option can be selected.")
    }
    if (is.null(valid.options)) {
        valid.options <- inla.validOptions()
    }
    match.opt <- pmatch(option,
                        valid.options,
                        nomatch = NA_integer_,
                        duplicates.ok = FALSE)
    if (any(is.na(match.opt))) {
        warning(
            paste0("Unknown or ambiguous INLA option '",
                   option[is.na(match.opt)],
                   "' will be ignored.",
                   " Valid options are: ",
                   paste0(valid.options, collapse = ", ")),
            immediate. = TRUE)
    }
    # Return the matched options, expanding any partial matches,
    # and removing any that were not matched
    valid.options[match.opt[!is.na(match.opt)]]
}

#' @rdname options
#' @export
`inla.getOption` <- function(option = NULL, list.out = NULL) {
    ## we get 'inla.call' separately to avoid infinite recursion
    default.opt <- inla.getOption.default()
    default.opt$inla.call <- inla.call.builtin()
    valid.opt <- names(default.opt)

    ## with no argument, return a named list of current values
    if (missing(option) || is.null(option)) {
        option <- valid.opt
    }
    if (is.null(list.out)) {
        list.out <- length(option) > 1L
    }

    option <- inla.matchOption(option, valid.opt, several.ok = TRUE)

    envir <- inla.get.inlaEnv()
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

    res <- list()
    for (i in seq_along(option)) {
        if (inla.is.element(option[i], opt)) {
            val <- list(inla.get.element(option[i], opt))
        } else {
            val <- list(inla.get.element(option[i], default.opt))
        }
        res[i] <- val
        names(res)[i] <- option[i]
    }

    if (!list.out) {
        names(res) <- NULL
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

    valid.opt <- inla.validOptions()

    `inla.setOption.core` <- function(
option, value) {
        envir <- inla.get.inlaEnv()
        option <- inla.matchOption(option, valid.opt, several.ok = FALSE)
        if (length(option) == 0L) {
            return(invisible())
        }
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
        inla.setOption.core("inla.mode", "compact")
    }

    arg <- inla.getOption("malloc.lib")
    if (is.null(arg) || arg == "default") {
        inla.setOption.core("malloc.lib", inla.getOption.default()$malloc.lib)
    }
    if (length(grep("^/", arg)) == 1 && !file.exists(arg)) {
        warning(paste0("User-defined library for option 'malloc.lib, ", arg,
                       ", does not exists. Use malloc.lib='default'"))
        inla.setOption.core("malloc.lib", inla.getOption.default()$malloc.lib)
    } else {
        arg <- match.arg(arg, c("default", "compiler", "je", "tc", "mi"), several.ok = FALSE)
        if (arg != "default" && arg != "compiler") {
            avail <- dir(paste0(dirname(inla.call.builtin()),"/malloc"), full.names = TRUE)
            idx <- grep(paste0("lib", arg, "malloc"), avail)
            if (length(idx) == 0 || length(idx) > 1) {
                if (FALSE) {
                    warning(paste0("Value for option 'malloc.lib, ", arg, ", is not availble. ",
                                   "Use malloc.lib='default'"))
                }
                inla.setOption.core("malloc.lib", "default")
            }
        }
    }
    
    arg <- inla.getOption("working.directory")
    if (!is.null(arg) && inla.anyMultibyteUTF8Characters(arg)) {
        warning(paste0("*** working.directory=[", arg,
                       "] contains multibyte characters. Normally, this will fail.\n"))
    }

    return(invisible())
}


## useful function
`inla.enabled.INLAjoint.features` <- function() {
    return (as.logical(inla.getOption('INLAjoint.features')))
}
