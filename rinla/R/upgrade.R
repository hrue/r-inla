
### The upgrade utility

#' Upgrade the INLA-package
#' 
#' Functions to upgrade the `INLA`-package to the current version.
#' 
#' 
#' @aliases inla.upgrade inla.update
#' @param lib Location to install the library.
#' @param testing If `TRUE`, then look for a test-version if the
#' INLA-package.
#' @param ask same argument as in `update.packages`
#' @return `inla.upgrade` will update the INLA package to the current
#' version, and `inla.update` do the same for backward compatibility. This
#' function is simple wrapper for `update.packages` using the INLA
#' repository.
#' @author Havard Rue \email{hrue@@r-inla.org}
#' @seealso `update.packages`
#' @rdname upgrade
#' @export
`inla.update` <- function(lib = NULL, testing = FALSE, ask = TRUE) {
    inla.upgrade(lib = lib, testing = testing, ask = ask)
}


#' @rdname upgrade
#' @export inla.upgrade
`inla.upgrade` <- function(lib = NULL, testing = FALSE, ask = TRUE) {

    repo <- c(CRAN = "https://cran.rstudio.com",
              INLA = paste0("https://inla.r-inla-download.org/R/",
                            if (testing) "testing" else "stable"))

    ## default timeout (60 is sometimes to low
    min.timeout <- 300
    opt <- options()
    if (is.null(opt$timeout) || opt$timeout < min.timeout) {
        options(timeout = min.timeout)
    } else {
        min.timeout <- 0
    }

    if (require("INLA", quietly = TRUE, lib.loc = lib, character.only = TRUE, warn.conflicts = FALSE)) {
        suppressWarnings(new.pack <- any(old.packages(repos = repo)[, 1] == "INLA"))
        if (new.pack) {
            if (.Platform$OS.type == "windows") {
                cat(
                    sep = "",
                    "\n *** Windows locks the INLA-package's DLL when its loaded, see",
                    "\n ***     https://cran.r-project.org/bin/windows/base/rw-FAQ.html",
                    "\n *** Section 4.8,  so you cannot update a package that is in use.",
                    "\n *** We recommend to remove the INLA-package and then reinstall, like",
                    "\n         remove.packages(\"INLA\")"
                )
                cat(paste0("\n         install.packages(\"INLA\", repos=\"", repo["INLA"], "\")"))
                cat("\n *** and then restart R.", "\n")
            } else {
                suppressWarnings(update.packages(repos = repo, oldPkgs = "INLA", ask = ask))
                cat("\n *** If the INLA-package was updated, you need to restart R and load it again.\n\n")
            }
        } else {
            cat("\n *** You already have the latest version.\n\n")
        }
    } else {
        install.packages(pkgs = "INLA", lib = lib, repos = repo, dependencies = TRUE)
        library("INLA")
    }

    if (min.timeout) {
        options(opt)
    }

    return(invisible())
}
