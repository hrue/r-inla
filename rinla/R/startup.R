## copied from the mgcv-code...

inla.print.version <- function()
{
    info = library(help=INLA)$info[[1]]
    version <- info[pmatch("Version", info)]
    date <- info[pmatch("Date", info)]
    um <- strsplit(version," ")[[1]]
    version <- um[nchar(um)>0][2]
    um <- strsplit(date," ")[[1]]
    date <- um[nchar(um)>0]

    hello <- paste("This is INLA ", version, ", dated ",
                   date[2], " (", date[9], date[11], ".", "\n",
                   "See www.r-inla.org/contact-us for how to get help.", 
                   sep="")
    packageStartupMessage(hello)
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

