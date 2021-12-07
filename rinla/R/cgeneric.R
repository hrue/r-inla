## Export: inla.cgeneric.define

## !\name{cgeneric.define}
## !\alias{cgeneric}
## !\alias{cgeneric.define}
## !\alias{inla.cgeneric.define}
## !
## !\title{cgeneric models}
## !
## !\description{A framework for defining latent models in C}
## !
## !\usage{
## !inla.cgeneric.define(model = NULL, shlib = NULL, n = 0L, debug = FALSE, ...)
## !}
## !
## !\arguments{
## !  \item{model}{The name of the model function}
## !  \item{shlib}{Name of the compiled object-file with \code{model}}
## !  \item{n}{The size of the model}
## !  \item{debug}{Logical. Turn on/off debugging}
## !}
## !
## !\author{Havard Rue \email{hrue@r-inla.org}}


`inla.cgeneric.define` <- function(model = NULL, shlib = NULL, n = 0L, debug = FALSE, ...) {
    stopifnot(!missing(n))
    stopifnot(as.integer(n) > 0)
    stopifnot(!missing(model))
    stopifnot(!missing(shlib))
    shlib <- normalizePath(shlib)
    if (!file.exists(shlib)) {
        shlib <- paste0(getwd(), "/", shlib)
        stopifnot(file.exists(shlib))
    }

    args <- list(...)
    if (any(names(args) == "")) {
        stop("The '...' argument in 'inla.cgeneric.define()' needs *named* arguments.")
    }

    cmodel <- list(
        f = list(
            model = "cgeneric",
            n = as.integer(n), 
            cgeneric = list(
                model = model, 
                shlib = shlib, 
                n = as.integer(n), 
                debug = debug,
                data = inla.cgeneric.convert.arguments(...)
            )
        )
    )
    class(cmodel) <- "inla.cgeneric"
    class(cmodel$f$cgeneric) <- "inla.cgeneric"
    return(cmodel)
}

`inla.cgeneric.convert.arguments` <-  function(...) {

    a <- list(...)
    if (any(names(a) == "")) {
        stop("Unnamed arguments are not allowed")
    }

    res <- list(ints = list(), doubles = list(), characters = list(), matrices = list(), sparse.matrices = list())
    for(idx in seq_along(a)) {
        xx <- a[idx][[1]]
        if (is.matrix(xx)) {
            xx <- list(c(nrow(xx), ncol(xx), as.numeric(xx)))
            names(xx) <- names(a)[idx]
            res$matrices <- c(res$matrices, xx)
        } else if (is(xx, "Matrix")) {
            xx <- inla.as.sparse(xx)
            xx <- list(c(nrow(xx), ncol(xx), length(xx@x), xx@i, xx@j, xx@x))
            names(xx) <- names(a)[idx]
            res$sparse.matrices <- c(res$sparse.matrices, xx)
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
