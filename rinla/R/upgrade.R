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
`inla.update` <- function(...) {
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
        js <- paste(readLines(paste0("https://api.github.com/repos/", repo,
                                     "/releases/latest"), warn = FALSE),
                    collapse = "")
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
    cat('Run remotes::install_github("hrue/r-inla", subdir = "rinla", ref = "master")\n')
    cat("If you have 'library(INLA)' in your '~/.Rprofile',  this will fail...\n")
    cat("On Windows this also needs Developer Mode (the repository has symlinks).\n\n")
    remotes::install_github("hrue/r-inla", subdir = "rinla", ref = "master")
    return (invisible())
}


#' @rdname upgrade
#' @export inla.upgrade
`inla.upgrade` <- function(...) {
    return (inla.update(...))
}
