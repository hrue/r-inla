## Export: inla.prune

##!\name{inla.prune}
##!\alias{inla.prune}
##!
##!\title{Prune the INLA-package}
##!
##!\description{Prune the INLA-package by removing binaries not supported by the running OS}
##!
##!\usage{
##!inla.prune()
##!}
##!\arguments{
##!}
##!\value{No value is returned.}
##!\author{Havard Rue \email{hrue@r-inla.org}}

`inla.prune` = function()
{
    dir.size <- function(d) {
        ## return the size of the directory (recursively) in Mb
        return (sum(file.info(list.files(d,
                                         all.files = TRUE,
                                         full.names = TRUE,
                                         recursive = TRUE))$size, na.rm = TRUE)/1024^2)
    }

    pkg <- installed.packages()
    if ("INLA" %in% pkg) {
        path <- pkg["INLA","LibPath"]
    } else {
        stop("Cannot find package 'INLA' using 'installed.packages()'")
    }
    bin.path <- paste0(path, "/INLA/bin")

    dd <- c()
    for (os in c("linux", "mac", "windows")) {
        for (arch in c("32bit", "64bit", "experimental")) {
            dd <- c(dd, paste0(bin.path, "/", os, "/", arch))
        }
    }
    native <- paste0(bin.path, "/", inla.os.type(), "/", inla.os.32or64bit(), "bit")
    dd <- setdiff(dd, native)

    size <- 0
    found <- FALSE
    for(d in dd) {
        info <- file.info(d)
        if (!is.na(info$isdir) && info$isdir) {
            found <- TRUE
            siz <- dir.size(d)
            if (siz < 0.1) siz <- 0
            cat("---> Found directory (size ", format(siz, dig = 1), "Mb): ", d, "\n", sep = "")

            ans <- askYesNo("Remove directory?", default = FALSE)
            if (is.na(ans)) {
                cat("-----> Cancel\n")
                return (invisible())
            }
            if (ans) {
                cat("-----> Remove...", "\n")
                size <- size + siz
                unlink(d, recursive = TRUE, force = TRUE)
            } 
        }
    }

    if (!found) {
        cat("---> Nothing to be pruned.\n")
    }
    if (size > 0) {
        cat("---> Removed ", format(size, dig = 1), "Mb in total.\n", sep = "")
    }

    return (invisible())
}
