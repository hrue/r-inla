## Export: inla.binary.install

## ! \name{inla.binary.install}
## ! \alias{inla.binary.install}
## ! \alias{binary.install}
## ! \title{Install alternative binary builds}
## ! \description{
## !    Install alternative binary builds.
## ! }
## ! \usage{
## !    inla.binary.install(debug = TRUE)
## ! }
## ! \arguments{
## !    \item{debug}{Logical. Turn on debugging messages if \code{TRUE}}
## ! }
## ! \details{
## ! \code{inla.binary.install()} will offer a menu of alternative
## ! (Linux) binary builds to be installed. Currently offered,  are builds
## ! for Ubuntu, CentOS, Manjaro (Arch Linux) and Fedora
## ! }
## ! \value{
## ! No value returned.
## ! }
## ! \author{Havard Rue \email{hrue@r-inla.org}}
## ! \examples{
## !   \dontrun{
## !     inla.binary.install()
## !   }
## ! }

`inla.binary.install` <- function(debug = TRUE) {
    ## install alternative binary builds

    show <- function(...) {
        if (debug) {
            msg <- paste(unlist(list(...)), sep = "", collapse = "")
            cat("* ", msg, "\n", sep = "")
        }
    }

    map.filename <- function(fnm) {
        return(gsub(" ", "%20", fnm))
    }

    stopifnot(inla.os.type() == "linux")
    version <- paste("Version_", inla.version("version"), sep = "")
    show("Looking for ", version)

    address <- "https://inla.r-inla-download.org/Linux-builds"
    Files <- paste0(address, "/FILES")
    fp <- url(Files, open = "r")
    ff <- readLines(fp)
    close(fp)
    ff <- ff[grep(version, ff)]
    nf <- length(ff)
    cat("  Available alternatives:\n")
    for (i in seq_len(nf)) {
        cat("  \t", paste0("Alternative ", i), " is ", ff[i], "\n")
    }
    cat("  ", "Chose alternative [", 1, ":", nf, "]", sep = "", "\n\t")
    ans <- scan(file = "", what = integer(), n = 1, quiet = TRUE)
    if (length(ans) == 0) {
        return(invisible())
    }
    if (!(ans %in% seq_len(nf))) stop("Not a valid choice. Exit.")
    fnm <- paste0(address, "/", ff[ans])
    show("Install file [", fnm, "]")
    pa <- searchpaths()
    pa <- pa[grep("/INLA$", pa)]
    if (length(pa) == 0) {
        stop("I cannot find '/INLA' in the searchpath(), so please retry after doing library('INLA')")
    }
    stopifnot(file.info(pa)$isdir)
    show("INLA is installed in [", pa, "]")
    pa <- paste0(pa, "/bin/linux")

    show("Checking for write access...")
    test.fnm <- paste0(pa, "/test-file---", as.integer(runif(1) * .Machine$integer.max), ".txt")
    test.result <- file.create(test.fnm, showWarnings = FALSE)
    if (test.result) {
        unlink(test.fnm, force = TRUE)
    } else {
        show(paste0("ERROR: No write access to [", pa, "]"))
    }

    show("Download file, please wait...")
    to.file <- paste0(pa, "/64bit-download-", date(), ".tgz")
    ret <- download.file(map.filename(fnm), to.file, quiet = TRUE)
    if (ret == 0) {
        ##
    } else {
        unlink(to.file, force = TRUE)
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
    to.dir <- paste0(pa, "/64bit-", date())
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

    show("Remove old 64bit directory...")
    unlink(to.dir, recursive = TRUE, force = TRUE)
    show("Done!")

    return(invisible())
}
