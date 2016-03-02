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
##! \description{
##!   Try to \code{load} the file,  if that fail,  try to \code{source} the file. 
##! }    
##! \author{Havard Rue \email{hrue@math.ntnu.no}}

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
        msg("source file [", filename, "]")
        source(filename, echo = TRUE)
    } else {
        msg("load file [", filename, "]")
    }

    return (invisible())
}



    

