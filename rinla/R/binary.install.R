#' Install the inla-binary
#' 
#' Install the inla-binary for the current system
#' 
#' @param ... Arguments passed onto `inla.stiles.install()`
#' @return This function is a wrapper to `inla.stiles.install` and return its value.
#' #' @author Havard Rue \email{hrue@@r-inla.org}
#' @examples
#' 
#'    \dontrun{
#'      inla.binary.install()
#' }
#'  
#' @rdname binary.install
#' @export inla.binary.install

`inla.binary.install` <- function(...) {
    return (inla.stiles.install(...))
}
