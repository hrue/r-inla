### The upgrade utility

#' Upgrade the INLA-package
#' 
#' Function to upgrade the `INLA`-package to the most recent version
#' 
#' @aliases inla.upgrade inla.update
#' @param ... `testing=TRUE` installs the newest release of any kind,
#'   prereleases included, instead of the one marked "Latest" on GitHub.
#'   The default installs the stable release.
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
    ## Install the RELEASE TARBALL, not the repository.
    ##
    ## remotes::install_github() unpacks the raw repository, which carries ~211
    ## symlinks (the documentation PDFs, inst/include, build-configs). Windows
    ## tar.exe cannot create a symlink without Developer Mode or admin rights,
    ## so on Windows the extraction aborts with "Invalid argument" once per
    ## symlink and the install never starts. R CMD build follows symlinks and
    ## stores real files, so the tarball attached to each release has none of
    ## them and installs anywhere with no special privileges.
    ##
    ## Falls back to install_github when no tarball is published for the
    ## newest release (older releases predate it), which is the old behaviour.
    repo <- "hrue/r-inla"
    url <- tryCatch({
        ## stable  = the release carrying GitHub's "Latest" badge, which skips
        ##           prereleases and can be pointed at any older build.
        ## testing = the newest release of any kind. /releases is ordered
        ##           newest first and does not skip prereleases.
        api <- if (testing) {
            paste0("https://api.github.com/repos/", repo, "/releases?per_page=1")
        } else {
            paste0("https://api.github.com/repos/", repo, "/releases/latest")
        }
        js <- paste(readLines(api, warn = FALSE), collapse = "")
        m <- regmatches(js, gregexpr('"browser_download_url"[^"]*"[^"]+"', js))[[1]]
        u <- sub('.*"([^"]+)"$', "\\1", m)
        u <- grep("/INLA_.*\\.tar\\.gz$", u, value = TRUE)
        if (length(u)) u[1] else NULL
    }, error = function(e) NULL, warning = function(w) NULL)

    if (!is.null(url)) {
        cat("\nInstalling the R package from the newest release:\n  ", url, "\n\n")
        f <- file.path(tempdir(), basename(url))
        utils::download.file(url, f, mode = "wb", quiet = FALSE)
        utils::install.packages(f, repos = NULL, type = "source")
        return (invisible())
    }

    cat('\nNo release tarball found; falling back to GitHub.\n')
    ## subdir = "package", NOT "rinla": the generated source archives no longer
    ## carry rinla/ (it holds the symlinks that make Windows tar.exe abort).
    ## package/ is the same sources with the symlinks resolved, so this works on
    ## every platform and no longer needs Developer Mode.
    cat('Run remotes::install_github("hrue/r-inla", subdir = "package", ref = "master")\n')
    cat("If you have 'library(INLA)' in your '~/.Rprofile',  this will fail...\n\n")
    remotes::install_github("hrue/r-inla", subdir = "package", ref = "master")
    return (invisible())
}


#' @rdname upgrade
#' @export inla.upgrade
`inla.upgrade` <- function(...) {
    return (inla.update(...))
}
