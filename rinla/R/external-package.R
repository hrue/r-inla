#' Return the path to the cgeneric-library for a pre-compiled external package
#' 
#' Return the path to the cgeneric-library for a pre-compiled external package
#' 
#' 
#' @param package the name of a package, given as a name or literal character
#' string
#' @return This function returns the complete path or `NULL` if file does
#' not exists
#' @author Havard Rue \email{hrue@@r-inla.org}
#' @rdname external-package
#' @export inla.external.lib
`inla.external.lib` <- function(package) {
    package <- as.character(substitute(package))
    if (inla.os("mac")) {
        fnm <- system.file(paste0("bin/mac/", inla.os.32or64bit(), "bit"), package = "INLA")
    } else if (inla.os("mac.arm64")) {
        fnm <- system.file(paste0("bin/mac.arm64/", inla.os.32or64bit(), "bit"), package = "INLA")
    } else if (inla.os("linux")) {
        fnm <- system.file(paste0("bin/linux/", inla.os.32or64bit(), "bit"), package = "INLA")
    } else if (inla.os("windows")) {
        fnm <- system.file(paste0("bin/windows/", inla.os.32or64bit(), "bit"), package = "INLA")
    } else {
        stop("Unknown OS")
    }

    fnm <- normalizePath(paste0(fnm, "/external/", package, "/lib", package, ".so"))
    return (if (file.exists(fnm)) fnm else NULL)
}
