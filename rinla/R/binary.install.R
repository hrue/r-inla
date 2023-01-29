## Export: inla.binary.install

## ! \name{inla.binary.install}
## ! \alias{inla.binary.install}
## ! \alias{binary.install}
## ! \title{Install alternative binary builds}
## ! \description{
## !    Install alternative binary builds.
## ! }
## ! \usage{
## !    inla.binary.install(os = c("CentOS Linux-6", "CentOS Linux-7", "CentOS Linux-8",
## !                               "CentOS Stream-8", "Rocky Linux-9", 
## !                               "Fedora-33", "Fedora-34", "Fedora Linux-35",
## !                               "Fedora Linux-36", "Fedora Linux-37",
## !                               "Manjaro Linux",
## !                               "Ubuntu-16.04", "Ubuntu-18.04", "Ubuntu-20.04", "Ubuntu-22.04"), 
## !                        path = NULL, verbose = TRUE, md5.check = TRUE)
## ! }
## ! \arguments{
## !    \item{os}{If \code{os} is given, install binary build for this \code{os}.
## !              If \code{os} is not given, chose \code{os} interactively
## !              among available builds.}
## !    \item{path}{character. The install path. If \code{NULL} the path is derived from \code{INLA} package}
## !    \item{verbose}{Logical. Verbose output if \code{TRUE}}
## !    \item{md5.check}{Logical. If \code{TRUE}, stop if md5-checksum-file is not present
## !                     or md5-checksum fail. If \code{FALSE}, ignore md5-checksum check.}
## ! }
## ! \details{
## ! Install a new binary for \code{os} unless
## ! \code{missing(os)}, for which the \code{os} is chosen
## ! interactively among the available builds.
## ! }
## ! \value{
## ! Return \code{TRUE} if installation was sucessful and \code{FALSE} if not.
## ! }
## ! \author{Havard Rue \email{hrue@r-inla.org}}
## ! \examples{
## !   \dontrun{
## !     inla.binary.install()
## !     inla.binary.install(os = "CentOS Linux-7")
## !     inla.binary.install(os = "CentOS Linux-7",  path = "~/local/bin/inla.binary")
## !   }
## ! }

`inla.binary.install` <- function(os = c("CentOS Linux-6",
                                         "CentOS Linux-7", 
                                         "CentOS Linux-8", 
                                         "CentOS Stream-8", 
                                         "Rocky Linux-9", 
                                         "Fedora-33", 
                                         "Fedora-34", 
                                         "Fedora Linux-35", 
                                         "Fedora Linux-36", 
                                         "Fedora Linux-37", 
                                         "Manjaro Linux", 
                                         "Ubuntu-16.04", 
                                         "Ubuntu-18.04", 
                                         "Ubuntu-20.04", 
                                         "Ubuntu-22.04"),
                                  path = NULL, 
                                  verbose = TRUE,
                                  md5.check = TRUE) {
    show <- function(...) {
        if (verbose) {
            msg <- paste(unlist(list(...)), sep = "", collapse = "")
            cat("* ", msg, "\n", sep = "")
        }
    }

    map.filename <- function(fnm) {
        return(gsub(" ", "%20", fnm))
    }

    random.num <- gsub("\\.", "", as.character(abs(rnorm(1))))
    os <- if (missing(os)) NULL else match.arg(os)
    stopifnot(inla.os.type() == "linux")
    version <- paste("Version_", inla.version("version"), sep = "")
    show("Looking for ", version, " and os=", if (!is.null(os)) os else "'<choose interactively>'")

    address <- "https://inla.r-inla-download.org/Linux-builds"
    Files <- paste0(address, "/FILES")
    fp <- url(Files, open = "r")
    ff <- readLines(fp)
    close(fp)
    ff <- ff[grep(version, ff)]
    nf <- length(ff)

    if (is.null(os)) {
        cat("  Available alternatives:\n")
        for (i in seq_len(nf)) {
            cat("  \t", paste0("Alternative ", i), " is ", ff[i], "\n")
        }
        cat("  ", "Chose alternative [", 1, ":", nf, "]", sep = "", "\n\t")
        ans <- scan(file = "", what = integer(), n = 1, quiet = TRUE)
        if (length(ans) == 0) {
            return (invisible(FALSE))
        }
    } else {
        ans <- grep(os, ff)
        if (length(ans) == 0) {
            stop(paste0("Sorry, os=", os, " is not available for ", version))
        }
        if (length(ans) > 1) {
            stop("Internal error. Please report to <help@r-inla.org>")
        }
    }
    if (!(ans %in% seq_len(nf))) stop("Not a valid choice. Exit.")

    fnm <- paste0(address, "/", ff[ans])
    show("Install file [", fnm, "]")
    external.path <- FALSE
    if (is.null(path)) {
        pa <- searchpaths()
        pa <- pa[grep("/INLA$", pa)]
        if (length(pa) == 0) {
            stop("I cannot find '/INLA' in the searchpath(), so please retry after doing library('INLA')")
        }
        stopifnot(file.info(pa)$isdir)
        show("INLA is installed in [", pa, "]")
        pa <- paste0(pa, "/bin/linux")
    } else {
        external.path <- TRUE
        pa <- normalizePath(path)
        if (!file.info(pa)$isdir) {
            stopifnot(dir.create(pa, recursive = TRUE))
        }
    }

    show("Checking for write access...")
    test.fnm <- paste0(pa, "/test-file---", random.num, ".txt")
    test.result <- file.create(test.fnm, showWarnings = FALSE)
    if (test.result) {
        unlink(test.fnm, force = TRUE)
    } else {
        show(paste0("ERROR: No write access to [", pa, "]"))
    }

    show("Download file, please wait...")
    to.file <- paste0(pa, "/64bit-download-", random.num, ".tgz")
    ret <- download.file(map.filename(fnm), to.file, quiet = TRUE)
    
    if (md5.check) {
        fnm.md5 <- gsub("/64bit.tgz", "/md5sum.txt", fnm)
        md5.file <- paste0(pa, "/64bit-download-md5sum-", random.num, ".txt")
        ret.md5 <- try(download.file(map.filename(fnm.md5), md5.file, quiet = TRUE),  silent = TRUE)
        if (!inherits(ret.md5, "try-error")) {
            md5.checksum <- scan(file=md5.file, what=character(), n=1, quiet = TRUE)
            if (md5.checksum == tools::md5sum(to.file)) {
                show("md5-checksum [", md5.checksum, "] OK.")
            } else {
                stop(paste0("md5-checksum [", md5.checksum, "] FAILED. Stop."))
            }
        } else {
            stop("No md5-checksum found. Run with 'md5.check=FALSE' to force install.")
        }
    }

    if (ret == 0) {
        ##
    } else {
        unlink(to.file, force = TRUE)
        if (md5.check) {
            unlink(md5.file, force = TRUE)
        }
        stop("Error downloading file. Abort.")
    }

    my.restore <- function() {
        show("Error. Will try to restore old configuration.")
        show("If unsuccessful, then reinstall R-INLA.")
        unlink(from.dir, recursive = TRUE, force = TRUE)
        unlink(to.file, force = TRUE)
        file.rename(to.dir, from.dir)
    }
    
    show("Rename old 64bit directory...")
    from.dir <- paste0(pa, "/64bit")
    to.dir <- paste0(pa, "/64bit-", random.num)
    ret <- file.rename(from.dir, to.dir)
    if (ret == TRUE) {
        ##
    } else {
        my.restore()
        stop("Error renaming old 64bit directory. Abort.")
    }
    show("Unpack file...")
    ret <- untar(to.file, exdir = dirname(to.file), verbose = FALSE)
    if (ret == 0) {
        ##
    } else {
        my.restore()
        stop("Error unpacking file. Abort.")
    }

    show("Remove temporary file...")
    unlink(to.file, force = TRUE)
    if (md5.check) {
        unlink(md5.file, force = TRUE)
    }
    show("Remove old 64bit directory...")
    unlink(to.dir, recursive = TRUE, force = TRUE)
    show("Done!")

    if (external.path) {
        cat("* Examples of usage:\n")
        cat("* \tinla.setOption(inla.call = \"", paste0(from.dir, "/inla.mkl.run"), "\")\n")
        cat("* \tinla.setOption(inla.call = \"", paste0(from.dir, "/inla.run"), "\")\n")
        cat("* \tinla.setOption(fmesher.call = \"", paste0(from.dir, "/fmesher.run"), "\")\n")
    }

   return(invisible(TRUE))
}
