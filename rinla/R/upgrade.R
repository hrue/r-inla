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
    cat('\nRun remotes::install_github("hrue/r-inla", subdir = "rinla", ref = "master")\n')
    cat("If you have 'library(INLA)' in your '~/.Rprofile',  this will fail...\n\n")
    remotes::install_github("hrue/r-inla", subdir = "rinla", ref = "master")
    return (invisible())
}


#' @rdname upgrade
#' @export inla.upgrade
`inla.upgrade` <- function(...) {
    return (inla.update(...))
}
