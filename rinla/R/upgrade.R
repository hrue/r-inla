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
##!\author{Havard Rue \email{hrue@r-inla.org}}
##!\seealso{\code{update.packages}}

### The upgrade utility

`inla.update` = function(lib = NULL, testing = FALSE, ask = TRUE)
{
    inla.upgrade(lib=lib, testing=testing, ask = ask)
}

`inla.upgrade` = function(lib = NULL, testing = FALSE, ask = TRUE)
{
    repo=c(CRAN = "https://cran.rstudio.com",
        INLA = paste("https://inla.r-inla-download.org/R/",
            (if (testing) "testing" else "stable"),  sep=""))
    if (require("INLA", quietly = TRUE, lib.loc = lib,
                character.only=TRUE, warn.conflicts=FALSE)) {
        suppressWarnings(new.pack <- any(old.packages(repos = repo)[,1] == "INLA"))
        if (new.pack) {
            if (.Platform$OS.type == "windows") {
                cat(sep="", 
                    "\n *** Windows locks the INLA-package's DLL when its loaded, see",
                    "\n ***     https://cran.r-project.org/bin/windows/base/rw-FAQ.html",
                    "\n *** Section 4.8,  so you cannot update a package that is in use.",
                    "\n *** We recommend to remove the INLA-package and then reinstall, like",
                    "\n     remove.packages(\"INLA\")")
                if (testing) {
                    cat("\n     install.packages(\"INLA\", repos=\"https://inla.r-inla-download.org/R/testing\")")
                } else {
                    cat("\n     install.packages(\"INLA\", repos=\"https://inla.r-inla-download.org/R/stable\")")
                }
                cat("\n *** and then restart R.", "\n")
            } else {
                suppressWarnings(update.packages(repos = repo, oldPkgs = "INLA", ask = ask))
                cat("\n *** If the INLA-package was updated, you need to restart R and load it again.\n\n")
            }
        } else {
            cat("\n *** You already have the latest version.\n\n")
        }
    } else {
        install.packages(pkgs = "INLA", lib = lib, repos = repo,
                         dependencies = TRUE)
        library("INLA")
    }

    return (invisible())
}
