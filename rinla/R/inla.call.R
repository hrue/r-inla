`inla.call.builtin` <- function() {
    ## cannot call inla.getOption() here as it leads to an infinite recursive call. do this
    ## manually instead.
    if (exists("inla.options", envir = inla.get.inlaEnv())) {
        opt <- get("inla.options", envir = inla.get.inlaEnv())
        if (!is.null(opt$inla.call)) {
            return (opt$inla.call)
        }
    }
    
    ## OBSOLETE: the package no longer ships a binary. Those builds moved to
    ## inst/obsolete/, so the lookups under bin/ that used to live here found
    ## nothing and failed with an empty filename in the message
    ## ("no such file " and then a blank), which said nothing about what to do
    ## next. A binary now arrives through inla.stiles.install(), which fetches
    ## one from the releases and sets inla.call to it, so that is what to say.
    ## This still raises, exactly as the old lookup did once the file was
    ## missing: there is genuinely no binary to return.
    stop("No inla binary is installed. Run 'inla.stiles.install()' to install one.",
         call. = FALSE)
}

`inla.call.no.remote` <- function() {
    ## return what is defined in options$inla.call except for 'remote', for which we revert back
    ## to the builtin one
    inla.call <- inla.getOption("inla.call")
    if (is.null(inla.call) || any(inla.strcasecmp(inla.call, c("remote", "inla.remote")))) {
        inla.call <- inla.call.builtin()
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
#' @param check Verify the file exists and answers `-V`. On failure the reason
#'        is returned as the `"error"` attribute rather than thrown, so a
#'        script can test the path without a `tryCatch`.
#'
#' @returns The path, invisibly when `check` fails. `NA_character_` when no
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
    ## The active setting first: inla.stiles.install() and any inla.setOption()
    ## write here, so this is what an inla() call would actually execute. Fall
    ## back to the launcher shipped inside the package, which is what a fresh
    ## install runs before any binary has been fetched.
    path <- tryCatch(inla.getOption("inla.call"), error = function(e) NULL)
    if (is.null(path) || !is.character(path) || !nzchar(path[1])) {
        path <- tryCatch(inla.call.builtin(), error = function(e) NA_character_)
    }
    if (is.null(path) || !is.character(path) || !nzchar(path[1])) {
        return(NA_character_)
    }
    path <- path[1]

    ## Absolute, so the value stays correct after a setwd() -- a Model.ini run
    ## is usually done from inside the model directory.
    abs <- suppressWarnings(normalizePath(path, mustWork = FALSE))
    if (nzchar(abs)) path <- abs

    if (!isTRUE(check)) {
        return(path)
    }
    if (!file.exists(path)) {
        attr(path, "error") <- "no such file"
        return(invisible(path))
    }
    out <- suppressWarnings(tryCatch(
        system2(path, "-V", stdout = TRUE, stderr = TRUE),
        error = function(e) NULL))
    if (is.null(out) || !any(grepl("version", out, ignore.case = TRUE))) {
        attr(path, "error") <- "did not answer -V"
        return(invisible(path))
    }
    attr(path, "version") <- sub(".*version:[[:space:]]*", "",
                                 grep("version", out, ignore.case = TRUE, value = TRUE)[1])
    path
}
