### The upgrade utility

#' Upgrade the INLA-package
#'
#' Upgrades the `INLA` R package and then the matching inla-program.
#'
#' @aliases inla.upgrade inla.update
#' @param ... `testing=TRUE` installs from the `testing` branch of the install
#'   repository instead of `master`.
#' @return `inla.upgrade` returns nothing
#' @author Havard Rue \email{hrue@@r-inla.org}
#' @seealso `update.packages`
#' @rdname upgrade
#' @export
`inla.update` <- function(...) {
    ## Same two channels as inla.stiles.install(), so the R package and the
    ## binary can be kept on the same one. They are separate downloads: this
    ## chooses which RELEASE's R tarball to install, the tag there chooses
    ## which release's binary. Mixing them is allowed (the version check then
    ## reports the mismatch) but rarely wanted.
    args <- list(...)
    ## Reject unknown names rather than drop them. A value read out of `...`
    ## is invisible when misspelled, so inla.update(testng = TRUE) would
    ## quietly install stable and look like it had worked. That already
    ## happened once here, with an argument that was accepted and ignored.
    nms <- if (is.null(names(args))) rep("", length(args)) else names(args)
    if (length(args) && any(!nzchar(nms))) {
        stop("arguments must be named: inla.update(testing = TRUE).", call. = FALSE)
    }
    unknown <- setdiff(nms, "testing")
    if (length(unknown)) {
        stop("unknown argument", if (length(unknown) > 1L) "s" else "", ": ",
             paste(unknown, collapse = ", "), ". Only `testing` is used.",
             call. = FALSE)
    }
    testing <- if (is.null(args$testing)) FALSE else args$testing
    if (!is.logical(testing) || length(testing) != 1L || is.na(testing)) {
        stop("`testing` must be TRUE or FALSE.", call. = FALSE)
    }
    ## The R package comes from the install repository, whose root IS the
    ## package: no subdir, and none of the symlinks that make Windows tar.exe
    ## abort. `testing` selects the branch there.
    repo <- inla.install.repo()
    ref  <- if (testing) "testing" else "master"

    if (!requireNamespace("remotes", quietly = TRUE)) {
        stop("package 'remotes' is needed here. install.packages(\"remotes\")",
             call. = FALSE)
    }
    cat('\nInstalling the R package:\n  remotes::install_github("', repo,
        '", ref = "', ref, '")\n', sep = "")
    cat("If you have 'library(INLA)' in your '~/.Rprofile',  this will fail...\n\n")
    remotes::install_github(repo, ref = ref)

    ## The binary belongs to the R package version, so it is refreshed here too.
    ## In a separate R process: this session still holds the namespace that was
    ## loaded before the upgrade, and would install the binary for that version.
    cat("\nInstalling the matching inla-program:\n\n")
    rbin <- file.path(R.home("bin"), if (.Platform$OS.type == "windows") "R.exe" else "R")
    st <- tryCatch(
        system2(rbin, c("--vanilla", "-q", "-e", shQuote("INLA::inla.stiles.install()")),
                env = paste0("R_LIBS=",
                             paste(.libPaths(), collapse = .Platform$path.sep))),
        error = function(e) 1L)
    if (!identical(as.integer(st), 0L)) {
        cat("\nThe inla-program was NOT installed. Restart R and run\n",
            "  INLA::inla.stiles.install()\n", sep = "")
    }
    return (invisible())
}


#' @rdname upgrade
#' @export inla.upgrade
`inla.upgrade` <- function(...) {
    return (inla.update(...))
}
