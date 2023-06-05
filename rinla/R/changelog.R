#' inla.changelog
#' 
#' List the recent changes in the inla-program and its R-interface
#' 
#' 
#' @aliases changelog inla.changelog inla.changes
#' @author Havard Rue \email{hrue@@r-inla.org}
#' @seealso \code{\link{inla}}
#' @rdname changelog
#' @export inla.changelog
`inla.changelog` <- function() {
    browseURL("https://github.com/hrue/r-inla/commits/devel")
    return(invisible())
}
