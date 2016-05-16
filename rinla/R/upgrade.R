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
##!inla.upgrade(lib = NULL, testing=FALSE, ask = TRUE)
##!inla.update(lib = NULL, testing=FALSE, ask = TRUE)
##!}
##!
##!\arguments{
##!  \item{lib}{Location to install the library.}
##!  \item{testing}{If \code{TRUE}, then look for a test-version if the
##!    INLA-package.}
##!  \item{ask}{same argument as in \code{update.packages}}
##!}
##!\value{\code{inla.upgrade} will update the INLA package to the current version, 
##! and \code{inla.update} do the same 
##!  for backward compatibility. This function is simple wrapper for
##!  \code{update.packages} using the INLA repository.}
##!\author{Havard Rue \email{hrue@math.ntnu.no}}
##!\seealso{\code{update.packages}}

### The upgrade utility

`inla.update` = function(lib = NULL, testing = FALSE, ask = TRUE)
{
    inla.upgrade(lib=lib, testing=testing, ask = ask)
}

`inla.upgrade` = function(lib = NULL, testing = FALSE, ask = TRUE)
{
    repo=c(CRAN = "https://cran.rstudio.com",
        INLA = paste("https://www.math.ntnu.no/inla/R/",
            (if (testing) "testing" else "stable"),  sep=""))
    if (require("INLA", quietly = TRUE, lib.loc = lib, character.only=TRUE)) {
        update.packages(repos = repo, oldPkgs = "INLA", ask = ask)
        cat("\n *** Please restart R to load package 'INLA'\n\n")
    } else {
        install.packages(pkgs = "INLA", lib = lib, repos = repo,
                         dependencies = TRUE)
        library("INLA")
    }

    return (invisible())
}
