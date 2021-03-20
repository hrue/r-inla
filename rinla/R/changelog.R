## Export: inla.changelog

## !\name{inla.changelog}
## !\alias{changelog}
## !\alias{inla.changelog}
## !\alias{inla.changes}
## !
## !\title{inla.changelog}
## !
## !\description{List the recent changes in the inla-program and its R-interface}
## !
## !\usage{
## !inla.changelog()
## !}
## !\author{Havard Rue \email{hrue@r-inla.org} }
## !\seealso{\code{\link{inla}}}

`inla.changelog` <- function() {
    browseURL("https://github.com/hrue/r-inla/commits/devel")
    return(invisible())
}
