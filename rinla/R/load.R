## Export: inla.load

##! \name{inla.load}
##! \alias{inla.load}
##! \title{Load or source a file}
##! \description{Load or source a file: (internal use)}
##! \usage{
##!     inla.load(filename, debug = TRUE)
##! }
##! \arguments{
##!   \item{filename}{The name of the file to be loaded, alternatively, sourced. }
##!   \item{debug}{Logical. Turn on/off debug information.}
##!  }
##! \value{
##!   None
##! }
##! \details{
##!   Try to \code{load} the file into the global environment,
##!   if that fail,  try to \code{source} the file into the
##!   global environment.
##! }    
##! \author{Havard Rue \email{hrue@r-inla.org}}

`inla.load` = function(filename, debug = TRUE)
{
    msg = function(...) {
        if (debug) {
            cat("inla.load: ", ...,  "\n", sep="")
        }
   }

    w = getOption("warn")
    options(warn = -1L)
    val = try(load(filename, envir = globalenv()), silent=TRUE)
    options(warn = w)

    if (inherits(val, "try-error")) {
        msg("source file [", filename, "] in the global environment")
        source(filename, echo = TRUE, local = FALSE)
    } else {
        msg("load file [", filename, "] in the global environment")
    }

    return (invisible())
}



    

