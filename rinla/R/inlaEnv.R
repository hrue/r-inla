## Define the environment for INLA used to store options and the
## model-list. Reuse the environment if it is there already.

if (!(exists(".inlaEnv") && is.environment(.inlaEnv))) {
    .inlaEnv <- new.env()
}



#' Return the internal environment used by INLA
#' 
#' A function which return the internal environment used by INLA
#' 
#' 
#' @aliases inla.get.inlaEnv get.inlaEnv
#' @returns This function returns the internal environment used by INLA to
#' keep internal variables.
#' @author Havard Rue \email{hrue@@r-inla.org}
#' @rdname inlaEnv
#' @export inla.get.inlaEnv
`inla.get.inlaEnv` <- function() {
    if (exists(".inlaEnv") && is.environment(.inlaEnv)) {
          return(.inlaEnv)
      }
    stop("Environment '.inlaEnv' does not exists and is required for INLA to work. Restart 'R'.")
}
