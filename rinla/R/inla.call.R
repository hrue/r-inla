`inla.call.no.remote` <- function() {
    ## return what is defined in options$inla.call except for 'remote', for
    ## which we revert back to the resolved binary. must = TRUE because the
    ## caller is about to RUN this: an NA here would reach system() as the
    ## literal string "NA".
    inla.call <- inla.getOption("inla.call")
    if (is.null(inla.call) || any(inla.strcasecmp(inla.call, c("remote", "inla.remote")))) {
        inla.call <- inla.binary.path(check = FALSE, must = TRUE)
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
#' @param must Raise instead of returning `NA_character_` when no binary can be
#'        resolved. Use it where the path is about to be executed, so that a
#'        missing binary is reported as such rather than reaching `system()` as
#'        the string `"NA"`.
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

`inla.binary.path` <- function(check = TRUE, must = FALSE) {
    ## The active setting, which inla.stiles.install() and any inla.setOption()
    ## write to: this is what an inla() call would actually execute. There is no
    ## second place to look any more. The package used to ship a binary and
    ## inla.call.builtin() found it; those builds moved to inst/obsolete/ and
    ## this function replaced it, so an unset inla.call means no binary at all.
    path <- tryCatch(inla.getOption("inla.call"), error = function(e) NULL)
    if (is.null(path) || !is.character(path) || !nzchar(path[1])) {
        if (isTRUE(must)) {
            stop("No inla binary is installed. Run 'inla.stiles.install()' to install one.",
                 call. = FALSE)
        }
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
