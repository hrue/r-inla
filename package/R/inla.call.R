`inla.remote.script` <- function(what = NULL) {
    if (is.null(what) || inla.os("windows")) {
        return (NULL)
    }

    f <- system.file(file.path("remote", paste0("inla.", what)), package = "INLA")
    if (nzchar(f)) {
        return(f)
    } 

    f <- system.file(file.path("remote", what), package = "INLA")
    if (nzchar(f)) {
        return(f)
    } 

    return(what)
}

`inla.call.no.remote` <- function() {
    ## return what is defined in options$inla.call except for 'remote', for
    ## which we revert back to the resolved binary. must = TRUE because the
    ## caller is about to RUN this: an NA here would reach system() as the
    ## literal string "NA".
    inla.call <- inla.getOption("inla.call")
    if (is.null(inla.call) || any(inla.strcasecmp(inla.call, c("remote", "inla.remote")))) {
        inla.call <- inla.binary.path(check = TRUE)
    }
    return(inla.call)
}

`inla.fmesher.call.builtin` <- function() {
    return ("fmesher.call.is.no.longer.in.use")
}

#' @title Path of the inla binary in use
#'
#' @description
#' The absolute path of the `inla` program this session would run. That is the
#' binary an `inla()` call uses, and the one to point at a `Model.ini` when
#' running a saved model without R:
#'
#' ```
#' <path> -v -t4:1 Model.ini
#' ```
#'
#' A model directory containing `Model.ini` comes from
#' `inla(..., keep = TRUE, working.directory = <dir>)`.
#'
#' @param check Check that the binary exists and return `NULL` if not.
#'        
#' @returns The path to the inla-binary. `NULL` when no 
#'        binary can be resolved at all.
#'
#' @examples
#' \dontrun{
#' inla.binary.path()
#' system2(inla.binary.path(), c("-v", "-t4:1", "Model.ini"))
#' }
#'
#' @seealso [inla.stiles.install()], [inla.stiles.status()]
#' @export inla.binary.path

`inla.binary.path` <- function(check = TRUE) {
    ## return the path to the inla-binary or an internal script. environment variable
    ## 'INLA_CALL' is checked first (expert use only, not documented), then 'inla.call' in
    ## inla-options, then local cache, and at last the bin-directory in the installed package.
    ## and error is issued if the binary is not found and 'check' is TRUE

    path <- Sys.getenv("INLA_CALL")
    if (nchar(path) == 0) {
        path <- NULL
    }
    if (is.null(path)) {
        path <- inla.getOption("inla.call")
    }

    if (!is.null(path)) {
        ## this will translate 'remote', 'inla.remote', 'keep', ....
        path <- inla.remote.script(path)
        path <- normalizePath(path, mustWork = FALSE)
    }

    user.defined <- if (is.null(path)) FALSE else TRUE

    ## Nothing set? Look where inla.stiles.install() puts it.
    if (is.null(path) || !is.character(path) || !nzchar(path)) {
        exe <- if (.Platform$OS.type == "windows") "inla.exe" else "inla.run"
        want <- file.path(paste0("Version_", inla.version("version")), "bin", exe)
        for (b in inla.binary.bases()) {
            full.path <- file.path(b, want)
            if (file.access(full.path) == 0) {
                path <- full.path
                break
            }
        }
    }

    if ((is.null(path) || file.access(path) != 0) && isTRUE(check)) {
        if (user.defined) {
            stop("User-defined inla-binary not found, try 'inla.setOption(inla.call=NULL)' or remove 'inla.call'-argument",
                 call. = FALSE)
        } else {
            stop("No inla-binary found. Run 'inla.stiles.install()' to install one.", call. = FALSE)
        }
    }

    if (!is.null(path)) {
        abs <- suppressWarnings(normalizePath(path, mustWork = check))
        if (nzchar(abs)) path <- abs
        if (isTRUE(check)) {
            stopifnot(file.access(path) == 0)
        }
    }

    return(path)
}

`inla.call.check` <- function(inla.call = NULL) {
    ## that that 'inla.call' exists
    stopifnot(!is.null(inla.call))

    if (file.access(inla.call) != 0) {
        stop(paste0("inla-binary not found: [", inla.call, "]"))
    }
}
