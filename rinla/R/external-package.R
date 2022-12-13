## Export: inla.external.lib

## !\name{inla.external.lib}
## !\alias{inla.external.lib}
## !
## !\title{Return the path to the cgeneric-library for a pre-compiled external package}
## !\description{
## ! Return the path to the cgeneric-library for a pre-compiled external package
## !}
## !\usage{
## !inla.external.lib(package)
## !}
## !\arguments{
## !\item{package}{the name of a package, given as a name or literal character string}
## !}
## !\value{This function returns the complete path or \code{NULL} if file does not exists}
## !\author{Havard Rue \email{hrue@r-inla.org}}

`inla.external.lib` <- function(package) {
    package <- as.character(substitute(package))
    fnm <- normalizePath(paste0(dirname(INLA:::inla.call.builtin()), "/external/", package, "/lib", package, ".so"))
    return (if (file.exists(fnm)) fnm else NULL)
}
