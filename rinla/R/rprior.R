#' @title Define prior in R
#' 
#' @description
#' A(n experimental) framework for defining a prior in R
#' 
#' @aliases rprior rprior.define inla.rprior.define
#' @param rprior An R-function returning the log-prior evaluated at its argument
#' @param ... Named list of variables that will be in the environment of `rprior`
#' @return An `inla.rprior`-object which can be used as a prior
#' @author Havard Rue \email{hrue@@r-inla.org}
#' @name rprior
#' @rdname rprior.define
#' @export inla.rprior.define
#' @examples
#'
#'  ## see example in inla.doc("rprior")
#'

`inla.rprior.define` <- function(rprior = NULL, ...) {
    stopifnot(!missing(rprior))
    stopifnot(is.function(rprior))
    args <- list(...)
    if (any(names(args) == "")) {
        stop("The '...' argument in 'inla.rprior.define()' needs _named_ arguments.")
    }

    env <- if (length(args) > 0) as.environment(args) else new.env()
    parent.env(env) <- .GlobalEnv
    environment(rprior) <- env

    fun.c <- inla.cmpfun(rprior)
    prior <- list(rprior = fun.c, code = inla.rprior.code())
    class(prior) <- "inla.rprior"
    return(prior)
}

`inla.rprior.code` <- function() {
    return ("c5c4fee74dc9299b6753b8605e303f59a1236bfa")
}

`inla.is.rprior` <- function(prior) {
    ## need this as often the class is ``lost''
    return (
        inherits(prior, "inla.rprior")
        ||
    (inherits(prior, "list")
        &&
        (
            (is.null(names(prior)) &&
             (length(prior) == 2L) &&
             is.function(prior[[1]]) &&
             is.character(prior[[2]]))
            ||
            (!is.null(prior$code) &&
             (prior$code == inla.rprior.code()))
        )
    )
    )
}
