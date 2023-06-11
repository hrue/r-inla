#' Joint-prior models
#' 
#' A framework for defining joint priors in R
#' 
#' 
#' @aliases jp jp.define inla.jp.define
#' @param jp The `jp`-function which returns the joint log-prior as a
#' function of argument `theta`. There is an optional second argument that
#' is a vector of `theta`-names. If second argument is not present,
#' argument `.theta.desc` will be added.
#' @param ... Named list of variables that defines the environment of `jp`
#' @returns This allows joint priors to be defined in `R`.
#' 
#' This function is for internal use only.
#' @author Havard Rue \email{hrue@@r-inla.org}
#' 
#' @name jp
#' @rdname jp
#' @export

`inla.jp.define` <- function(jp = NULL, ...) {
    stopifnot(!missing(jp))
    stopifnot(is.function(jp))
    args <- list(...)
    if (any(names(args) == "")) {
        stop("The '...' argument in 'inla.jp.define()' needs *named* arguments.")
    }

    a <- names(formals(jp))
    if (length(a) == 1) {
        ## add second argument
        fun <- inla.eval(paste0("function(", a[1], ", .theta.desc = NULL) NULL"))
        body(fun) <- body(jp)
    } else if (length(a) == 2) {
        fun <- jp
    } else {
        stop("Number of arguments in 'jp'-function must be 1 or 2.")
    }

    env <- if (length(args) > 0) as.environment(args) else new.env()
    parent.env(env) <- .GlobalEnv
    environment(fun) <- env

    fun.c <- inla.cmpfun(fun)
    rjp <- list(model = fun.c) ## maybe we need something additional later
    class(rjp) <- "inla.jp"

    return(rjp)
}
