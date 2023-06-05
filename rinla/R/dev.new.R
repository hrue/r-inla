#' Opens a new device
#' 
#' Open a new device using \code{\link{dev.new}} unless using RStudio
#' 
#' 
#' @param ... Optional arguments to \code{\link{dev.new}}
#' @return The value of \code{\link{dev.new}} if not running RStudio, otherwise
#' \code{NULL}
#' @author Havard Rue \email{hrue@@r-inla.org}
#' @rdname dev.new
#' @export inla.dev.new
`inla.dev.new` <- function(...) {
    ## If running in RStudio then don't open a new device,  otherwise,  do.
    dev <- getOption("device")
    if (is.character(dev) && inla.strncasecmp(dev, "RStudioGD")) {
        ret <- NULL
    } else {
        ret <- dev.new(...)
    }
    if (exists("inla.dev.new.hook") && is.function(inla.dev.new.hook)) {
        do.call("inla.dev.new.hook")
    }
    return(invisible(ret))
}
