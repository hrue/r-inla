## copied from the mgcv-code...

## this does not work...
##  invisible(utils::globalVariables(c("low", "high", "spde", "internal")))
##  invisible(utils::suppressForeignCheck(c("low", "high", "spde", "internal")))


inla.print.version <- function()
{
    info = library(help=INLA)$info[[1]]
    if (!is.null(info)) {
        version <- info[pmatch("Version", info)]
        built = strsplit(info[pmatch("Built", info)], "; ")[[1]][3]
        date <- info[pmatch("Date", info)]
        um <- strsplit(version," ")[[1]]
        version <- um[nchar(um)>0][2]
        um <- strsplit(date," ")[[1]]
        date <- um[nchar(um)>0]

        hello <- paste0("This is INLA_", version, 
                        " built ", built, ".", "\n",
                        "See www.r-inla.org/contact-us for how to get help.")
        if (!inla.os("windows")) {
            hello = paste0(hello,
                           "\n",
                           "To enable PARDISO sparse library; see inla.pardiso()")
        }
        packageStartupMessage(hello)
        if (!inla.os("windows") && inla.os.is.32bit()) {
            warning(paste0("INLA_",  version,
                           ": 32bit binaries are no longer supported!"))
        }
    }
}

.onLoad <- function(...)
{
    ## nothing for the moment
}

.onAttach <- function(...)
{ 
    inla.print.version()
}

.onUnload <- function(libpath) {
    ## library.dynam.unload("INLA", libpath)
}

