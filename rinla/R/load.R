#' Load or source a file
#' 
#' Load or source a file: (internal use)
#' 
#' Try to `load` the file into the global environment, if that fail, try
#' to `source` the file into the global environment.
#' 
#' @param filename The name of the file to be loaded, alternatively, sourced.
#' @param debug Logical. Turn on/off debug information.
#' @return None
#' @author Havard Rue \email{hrue@@r-inla.org}
#' @rdname load
#' @export inla.load
#' @keywords internal
`inla.load` <- function(filename, debug = TRUE) {
    msg <- function(...) {
        if (debug) {
            cat("inla.load: ", ..., "\n", sep = "")
        }
    }

    w <- getOption("warn")
    options(warn = -1L)
    val <- try(load(filename, envir = globalenv()), silent = TRUE)
    options(warn = w)

    if (inherits(val, "try-error")) {
        msg("source file [", filename, "] in the global environment")
        source(filename, echo = TRUE, local = FALSE)
    } else {
        msg("load file [", filename, "] in the global environment")
    }

    return(invisible())
}
