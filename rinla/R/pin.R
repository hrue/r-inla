#' Control main thread pinning for INLA (experimental)
#'
#' Control main thread pinning for INLA (experimental)
#' 
#' @aliases inla.pin inla.unpin
#' @return No value is returned.
#' @author Havard Rue \email{hrue@@r-inla.org}
#' @examples
#' inla.pin()
#' inla.unpin()

#' @rdname inla.pin
#' @export
#' @details * `inla.pin` set OMP variables for pinning
`inla.pin` <- function() {
    Sys.setenv(OMP_PLACES='cores')
    Sys.setenv(OMP_PROC_BIND='spread,close')
    cat("pin: set OMP_PLACES and OMP_PROC_BIND\n")
    return (invisible())
}

#' @rdname inla.pin
#' @export
#' @details * `inla.unpin` unset OMP variables for pinning
`inla.unpin` <- function() {
    Sys.unsetenv('OMP_PLACES')
    Sys.unsetenv('OMP_PROC_BIND')
    cat("unpin: unset OMP_PLACES and OMP_PROC_BIND\n")
    return (invisible())
}
