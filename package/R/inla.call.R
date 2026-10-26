`inla.remote.script` <- function(what = c("inla.remote", "inla.submit")) {
    ## What inla.call="remote" and "submit" translate into: the scripts in
    ## inst/remote/, not inst/bin/remote/, which does not exist and so returned
    ## "" and ran nothing. They are bash over ssh, so NULL on Windows.
    what <- match.arg(what)
    if (inla.os("windows")) {
        return(NULL)
    }
    f <- system.file(file.path("remote", what), package = "INLA")
    if (!nzchar(f)) {
        stop("Missing from the installation: inst/remote/", what, call. = FALSE)
    }
    return(f)
}

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

    ## "remote" and "submit" are names inla() translates into its own scripts,
    ## not binaries, so they are neither warned about nor returned as a path.
    special <- c("remote", "inla.remote", "submit", "inla.submit")
    if (!is.null(path) && is.character(path) && tolower(path[1]) %in% special) {
        return(path[1])
    }

    ## A binary set here is used as given, with no version check. Say so once.
    if (!is.null(path) && is.character(path) && nzchar(path[1])) {
        key <- paste0(".binary.set.warned:", path[1])
        e <- inla.get.inlaEnv()
        if (!(exists(key, envir = e) && isTRUE(get(key, envir = e)))) {
            assign(key, TRUE, envir = e)
            warning("using the inla binary set in inla.call: ", path[1],
                    "\n  Its version is not checked against the R package.",
                    call. = FALSE)
        }
    }

    ## Nothing set? Look where inla.stiles.install() puts it.
    ##
    ## The option is normally set at attach time, or from the block the
    ## installer writes into ~/.Rprofile. Neither happens in a VANILLA session,
    ## and R CMD build runs vignettes with --vanilla: no profile is read at all,
    ## and .onAttach returns early because R_CMD is set. So every vignette that
    ## fits a model failed with "No inla binary is installed" even though one
    ## was installed and working, which is why releases shipped without them.
    ##
    ## The binary lives at a known path, so look there rather than depend on a
    ## profile having run. Same layouts the installer searches: <root>/bin/inla*
    ## on unix, and the flat Windows bundle where inla.exe sits at the top.
    ## Same order the installer writes in: inside the package first, so the
    ## binary found is the one that shipped with this version, then the cache.
    if (is.null(path) || !is.character(path) || !nzchar(path[1])) {
        exe <- if (.Platform$OS.type == "windows") "inla.exe" else "inla*"
        ## The release this package asks for comes FIRST. Globbing the install
        ## directories and taking the first hit sorts alphabetically, so with
        ## Version_26.08.31 and Version_26.09.19 both installed the older one
        ## won and became the default silently.
        need <- tryCatch(utils::packageDescription("INLA")[["Config/INLA/BinaryVersion"]],
                         error = function(e) NULL)
        want <- if (!is.null(need) && nzchar(need)) {
            c(paste0("Version_", need), paste0("v", need))
        } else character(0)
        ## The declared version in ANY location beats a mismatch in any other,
        ## so the exact match is searched everywhere before any fallback.
        bases <- inla.binary.bases()
        exact <- as.vector(t(outer(bases, want, file.path)))
        ## Newest first: taking the first hit gave the OLDEST. Tags are
        ## zero-padded, so decreasing name order is version order.
        other <- sort(unlist(lapply(bases, inla.binary.dirs)), decreasing = TRUE)
        ## Then "latest", where the release build unpacks the binary it renders
        ## the vignettes with, and where older installs put it.
        other <- c(other, file.path(bases, "latest"))
        pick <- function(dirs) {
            ## as.character: with no dirs unlist() gives NULL, and
            ## file.exists(NULL) is an error, not FALSE.
            cand <- as.character(unlist(lapply(dirs, function(d)
                c(Sys.glob(file.path(d, "bin", exe)), Sys.glob(file.path(d, exe))))))
            cand <- cand[file.exists(cand) & !dir.exists(cand)]
            if (length(cand)) cand[1] else NULL
        }
        path <- pick(exact)
        if (is.null(path)) {
            ## No binary for the declared version: warn once. The attach
            ## message is not enough, it is skipped when R_CMD is set.
            path <- pick(other)
            if (!is.null(path) && !is.null(need) && nzchar(need)) {
                key <- paste0(".binary.mismatch.warned:", path)
                e <- inla.get.inlaEnv()
                if (!(exists(key, envir = e) && isTRUE(get(key, envir = e)))) {
                    assign(key, TRUE, envir = e)
                    warning("no inla binary for INLA ", need, "; using ", path,
                            ".\n  Run inla.stiles.install() to get the matching one.",
                            call. = FALSE)
                }
            }
        }
    }


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
