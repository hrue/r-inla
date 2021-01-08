## Export: inla.get.inlaEnv

## !\name{inla.get.inlaEnv}
## !\alias{inla.get.inlaEnv}
## !\alias{get.inlaEnv}
## !\title{Return the internal environment used by INLA}
## !\description{A function which return the internal environment used by INLA}
## !\usage{
## !inla.get.inlaEnv()
## !}
## !\arguments{
## !}
## !
## !\value{%%
## !  This function returns the internal environment used by INLA to
## !  keep internal variables.
## !}
## !\author{Havard Rue \email{hrue@r-inla.org}}


## Define the environment for INLA used to store options and the
## model-list. Reuse the environment if it is there already.

if (!(exists(".inlaEnv") && is.environment(.inlaEnv))) {
    .inlaEnv <- new.env()
}

`inla.get.inlaEnv` <- function() {
    if (exists(".inlaEnv") && is.environment(.inlaEnv)) {
          return(.inlaEnv)
      }
    stop("Environment '.inlaEnv' does not exists and is required for INLA to work. Restart 'R'.")
}
