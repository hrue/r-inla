### The upgrade utility

#' Upgrade the INLA-package
#' 
#' Function to upgrade the `INLA`-package to the most recent version
#' 
#' @aliases inla.upgrade inla.update
#' @param ... Arguments not used
#' @return `inla.upgrade` returns nothing
#' @author Havard Rue \email{hrue@@r-inla.org}
#' @seealso `update.packages`
#' @rdname upgrade
#' @export
`inla.update` <- function(channel = c("stable", "testing"), ...) {
    ## Same two channels as inla.stiles.install(), so the R package and the
    ## binary can be kept on the same one. They are separate downloads: the
    ## channel here chooses which RELEASE's R tarball to install, and the
    ## channel there chooses which release's binary. Mixing them is allowed
    ## (Config/INLA/BinaryVersion then reports the mismatch) but rarely wanted.
    channel <- match.arg(channel)
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
        api <- if (channel == "stable") {
            paste0("https://api.github.com/repos/", repo, "/releases/latest")
        } else {
            paste0("https://api.github.com/repos/", repo, "/releases?per_page=1")
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
