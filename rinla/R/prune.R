`inla.dir.size` <- function(d) {
    ## return the disk usage of a directory, recursively, in Mb
    return(sum(file.info(list.files(d,
        all.files = TRUE,
        full.names = TRUE,
        recursive = TRUE
    ))$size, na.rm = TRUE) / 1024^2)
}



#' Prune the INLA-package
#' 
#' Prune the INLA-package by deleting binary files not supported by the running
#' OS
#' 
#' 
#' @aliases inla.prune prune
#' @param ask Logical. If TRUE, then ask for user confirmation before deleting.
#' If FALSE, then delete without user confirmation.
#' @return No value is returned.
#' @author Havard Rue \email{hrue@@r-inla.org}
#' @rdname prune
#' @export inla.prune
`inla.prune` <- function(ask = TRUE) {
    pkg <- installed.packages()
    if ("INLA" %in% pkg) {
        path <- pkg["INLA", "LibPath"]
    } else {
        stop("Cannot find package 'INLA' using 'installed.packages()'")
    }
    bin.path <- paste0(path, "/INLA/bin")

    arm <- paste0(bin.path, "/mac.arm64")
    dd <- arm
    for (os in c("linux", "mac", "windows")) {
        for (arch in c("32bit", "64bit", "experimental")) {
            dd <- c(dd, paste0(bin.path, "/", os, "/", arch))
        }
    }
    if (inla.os("mac.arm64")) {
        native <- arm
    } else {
        native <- paste0(bin.path, "/", inla.os.type(), "/", inla.os.32or64bit(), "bit")
    }
    dd <- setdiff(dd, native)

    size <- 0
    found <- FALSE
    for (d in dd) {
        info <- file.info(d)
        if (!is.na(info$isdir) && info$isdir) {
            found <- TRUE
            siz <- inla.dir.size(d)
            if (siz < 0.1) siz <- 0
            cat("---> Found directory (size ", format(siz, dig = 1), "Mb): ", d, "\n", sep = "")
            ans <- if (ask) askYesNo("Remove directory?", default = FALSE) else TRUE
            if (is.na(ans)) {
                cat("---> Cancel\n")
                return(invisible())
            }
            if (ans) {
                cat("---> Remove...", "\n")
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

    return(invisible())
}

`inla.prune.check` <- function() {
    ## return the ammount that could be removed in Mb
    pkg <- installed.packages()
    if ("INLA" %in% pkg) {
        path <- pkg["INLA", "LibPath"]
    } else {
        return(0)
    }
    bin.path <- paste0(path, "/INLA/bin")
    arm <- paste0(bin.path, "/mac.arm64")
    dd <- arm
    for (os in c("linux", "mac", "windows")) {
        for (arch in c("32bit", "64bit", "experimental")) {
            dd <- c(dd, paste0(bin.path, "/", os, "/", arch))
        }
    }
    if (inla.os("mac.arm64")) {
        native <- arm
    } else {
        native <- paste0(bin.path, "/", inla.os.type(), "/", inla.os.32or64bit(), "bit")
    }
    dd <- setdiff(dd, native)

    size <- 0
    for (d in dd) {
        info <- file.info(d)
        if (!is.na(info$isdir) && info$isdir) {
            size <- size + inla.dir.size(d)
        }
    }

    return(round(dig = 1, size))
}
