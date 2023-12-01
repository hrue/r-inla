#' @title Describe and check the PARDISO support in R-INLA
#' 
#' @description
#' `inla.pardiso()` describes the `PARDISO` support in R-INLA, how to
#' get the license key and enable it in the `R-INLA` package.
#' `inla.pardiso.check()` check if the `PARDISO` support is working.
#' 
#' @name pardiso
#' @aliases pardiso inla.pardiso inla.pardiso.check
#' @author Havard Rue \email{hrue@@r-inla.org}
#' @rdname pardiso
#' @export inla.pardiso

`inla.pardiso` <- function() {
    browseURL("https://www.panua.ch")
}

#' @rdname pardiso
#' @export inla.pardiso.check
`inla.pardiso.check` <- function() {
    t.dir <- inla.tempdir()
    inla.set.sparselib.env(inla.dir = t.dir)
    if (inla.os("linux") || inla.os("mac") || inla.os("mac.arm64")) {
        ret <- system(paste(shQuote(inla.call.no.remote()), "-m pardiso"), intern = TRUE)
    } else if (inla.os("windows")) {
        ret <- system(paste(shQuote(inla.call.no.remote()), "-m pardiso"), intern = TRUE)
    } else {
        stop("\n\tNot supported architecture.")
    }
    unlink(t.dir, recursive = TRUE)

    return(ret)
}
