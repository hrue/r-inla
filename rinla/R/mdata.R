#' Create an mdata-object for INLA
#' 
#' This defines an mdata-object for matrix valued response-families
#' 
#' 
#' @aliases mdata inla.mdata is.inla.mdata as.inla.mdata print.inla.mdata
#' @param y The response vector/matrix
#' @param ... Additional vectors/matrics of same length as `y`
#' @param object Any `R`-object
#' @param x An mdata object
#' @return An object of class `inla.mdata`.  There is method for
#' `print`.
#' 
#' `is.inla.mdata` returns `TRUE` if `object` inherits from
#' class `inla.mdata`, otherwise `FALSE`.
#' 
#' `as.inla.mdata` returns an object of class `inla.mdata`
#' @note It is often required to set `Y=inla.mdata(...)` and then define
#' the formula as `Y~...`, especially when used with `inla.stack`.
#' @author Havard Rue
#' @seealso [inla()]
#' @rdname mdata
#' @export inla.mdata
`inla.mdata` <- function(y, ...) {
    names.ori <- as.list(match.call())[-1]
    y.obj <- as.list(as.data.frame(list(y)))
    x.obj <- as.list(as.data.frame(list(...)))
    if (length(list(...)) == 0) {
        ncols <- ncol(as.data.frame(list(y)))
    } else {
        ncols <- c(
            ncol(as.data.frame(list(y))),
            unlist(lapply(
                list(...),
                function(x) if (is.null(ncol(x))) 1 else ncol(x)
            ))
        )
    }
    names(y.obj) <- paste("Y", 1:length(y.obj), sep = "")
    if (length(x.obj) > 0) {
        names(x.obj) <- paste("X", 1:length(x.obj), sep = "")
        obj <- c(y.obj, x.obj)
    } else {
        obj <- y.obj
    }
    obj <- as.data.frame(obj)
    attr(obj, "inla.ncols") <- c(length(ncols), ncols)
    class(obj) <- c("inla.mdata", "data.frame")
    attr(obj, "names.ori") <- names.ori
    return(obj)
}

#' @rdname mdata
#' @method print inla.mdata
#' @export
`print.inla.mdata` <- function(x, ...) {
    cat("inla.cols = ", attr(x, "inla.ncols", exact = TRUE), "\n")
    print(as.data.frame(unclass(x)), ...)
}

#' @rdname mdata
#' @export
`as.inla.mdata` <- function(object) {
    if (is.inla.mdata(object)) {
        return(object)
    }
    object <- as.list(as.data.frame(object))
    if (length(object) == 1) {
        names(object) <- "Y1"
        ncols <- c(1, 1)
    } else {
        names(object) <- c("Y1", paste("X", 1:(length(object) - 1), sep = ""))
        ncols <- c(2, 1, length(object) - 1)
    }
    object <- as.data.frame(object)
    warning("Guess that ncol(response) == 1. Otherwise,  please modify 'names(object)'.")
    attr(object, "inla.ncols") <- ncols
    class(object) <- c("inla.mdata", "data.frame")
    return(object)
}

#' @rdname mdata
#' @export
`is.inla.mdata` <- function(object) {
    return(!is.null(attr(object, "inla.ncols", exact = TRUE)) ||
        inherits(object, "inla.mdata"))
}
