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
## !inla.cgeneric.define(model = NULL, filename = NULL, n = 0L, debug = FALSE, ...)
## !}
## !
## !\arguments{
## !  \item{model}{The name of the model function}
## !  \item{filename}{Name of the compiled object-file with \code{model}}
## !  \item{n}{The size of the model}
## !  \item{debug}{Logical. Turn on/off debugging}
## !}
## !
## !\author{Havard Rue \email{hrue@r-inla.org}}


`inla.cgeneric.define` <- function(model = NULL, filename = NULL, n = 0L, debug = FALSE, ...) {
    stopifnot(!missing(n))
    stopifnot(as.integer(n) > 0)
    stopifnot(!missing(model))
    stopifnot(!missing(filename))
    filename <- normalizePath(filename)
    if (!file.exists(filename)) {
        filename <- paste0(getwd(), "/", filename)
        stopifnot(file.exists(filename))
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
                file = filename, 
                n = as.integer(n), 
                debug = debug
            )
        )
    )
    class(cmodel) <- "inla.cgeneric"
    class(cmodel$f$cgeneric) <- "inla.cgeneric"
    return(cmodel)
}
