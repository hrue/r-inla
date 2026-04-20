#' cloglike models
#' 
#' A framework for defining likelihood models in C
#' 
#' @aliases cloglike cloglike.define
#' @param model The name of the loglike function
#' @param shlib Name of the compiled object-file with `model`
#' @param debug Logical. Turn on/off debugging
#' @param ... Additional arguments, required by `inla.cloglike.define()` to be named arguments
#' @author Havard Rue \email{hrue@@r-inla.org}
#' 
#' @name cloglike
#' @rdname cloglike
#' @export
`inla.cloglike.define` <- function(model = NULL, shlib = NULL, debug = FALSE, ...)
{
    stopifnot(!missing(model))
    stopifnot(!missing(shlib))
    shlib <- normalizePath(shlib)
    if (!file.exists(shlib)) {
        shlib <- paste0(getwd(), "/", shlib)
        stopifnot(file.exists(shlib))
    }

    args <- list(...)
    if (any(names(args) == "")) {
        stop("The '...' argument in 'inla.cloglike.define()' needs *named* arguments.")
    }

    cmodel <- list(
        model = "cloglike",
        cloglike = list(
            model = model, 
            shlib = shlib, 
            debug = debug,
            data = inla.cloglike.convert.arguments(model = model, shlib = shlib, debug = as.integer(debug), ...)
        )
    )
    class(cmodel) <- "inla.cloglike"

    return(cmodel)
}

`inla.cloglike.convert.arguments` <-  function(model = NULL, shlib = NULL, debug = 0L, ...)
{
    a <- c(model = model, shlib = shlib, debug = debug, list(...))
    if (any(names(a) == "")) {
        stop("Unnamed arguments are not allowed")
    }

    res <- list(ints = list(), doubles = list(), characters = list(), matrices = list(), smatrices = list())
    for(idx in seq_along(a)) {
        xx <- a[idx][[1]]
        if (is.matrix(xx)) {
            ## store columnwise and not rowwise, so we need to revert the internal storage
            xx <- list(c(nrow(xx), ncol(xx), as.numeric(matrix(t(xx), nrow = ncol(xx), ncol = nrow(xx)))))
            names(xx) <- names(a)[idx]
            res$matrices <- c(res$matrices, xx)
        } else if (is(xx, "Matrix")) {
            if (!inherits(xx,"dgTMatrix")) {
                xx <- inla.as.sparse(xx)
            }
            xx <- list(c(nrow(xx), ncol(xx), length(xx@x), xx@i, xx@j, xx@x))
            names(xx) <- names(a)[idx]
            res$smatrices <- c(res$smatrices, xx)
        } else if (is.integer(xx)) {
            res$ints <- c(res$ints, a[idx])
        } else if (is.character(xx)) {
            res$characters <- c(res$characters, a[idx])
        } else if (is.double(xx)) {
            res$doubles <- c(res$doubles, a[idx])
        } else {
            stop(paste0("Unknown type: ", names(a)[idx], " ", idx))
        }
    }
    return (res)
}
