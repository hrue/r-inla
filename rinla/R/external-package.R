#' Return the path to the cgeneric-library for a pre-compiled external package
#' 
#' Return the path to the cgeneric-library for a pre-compiled external package
#' 
#' 
#' @param package the name of a package, given as a name or literal character
#' string
#' @return This function returns the complete path or \code{NULL} if file does
#' not exists
#' @author Havard Rue \email{hrue@@r-inla.org}
#' @rdname external-package
#' @export inla.external.lib
`inla.external.lib` <- function(package) {
    package <- as.character(substitute(package))
    fnm <- normalizePath(paste0(dirname(inla.call.builtin()), "/external/",
                                package, "/lib", package, ".so"))
    return (if (file.exists(fnm)) fnm else NULL)
}
