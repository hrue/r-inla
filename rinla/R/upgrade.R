## Export: inla.upgrade inla.update

##!\name{inla.upgrade}
##!\alias{inla.upgrade}
##!\alias{inla.update}
##!
##!\title{Upgrade the INLA-package}
##!
##!\description{Functions to upgrade the \code{INLA}-package to the
##!  current version.}
##!
##!\usage{
##!inla.upgrade(lib = NULL, testing=FALSE)
##!inla.update(lib = NULL, testing=FALSE)
##!}
##!
##!\arguments{
##!  \item{lib}{Location to install the library.}
##!  \item{testing}{If \code{TRUE}, then look for a test-version if the
##!    INLA-package.}
##!}
##!\value{\code{inla.upgrade} will update the INLA package to the current version, 
##! and \code{inla.update} do the same 
##!  for backward compatibility. This function is simple wrapper for
##!  \code{update.packages} using the INLA repository.}
##!\author{Havard Rue \email{hrue@math.ntnu.no}}
##!\seealso{\code{update.packages}}

### The upgrade utility

`inla.update` = function(lib = NULL, testing = FALSE)
{
    inla.upgrade(lib=lib, testing=testing)
}

`inla.upgrade` = function(lib = NULL, testing = FALSE)
{
    repo=c(CRAN = "http://cran.rstudio.com",
        INLA = paste("http://www.math.ntnu.no/inla/R/",
            (if (testing) "testing" else "stable"),  sep=""))
    if (require("INLA", quietly = TRUE, lib.loc = lib, character.only=TRUE)) {
        update.packages(repos = repo, oldPkgs = "INLA")
        cat("\n *** Please restart R to load package 'INLA'\n\n")
    } else {
        install.packages(pkgs = "INLA", lib = lib, repos = repo)
        library("INLA")
    }

    return (invisible())
}
