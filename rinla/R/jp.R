## Export: inla.jp.define

##!\name{jp.define}
##!\alias{jp}
##!\alias{jp.define}
##!\alias{inla.jp.define}
##!\title{Joint-prior models}
##!\description{A framework for defining joint priors in R}
##!
##!\usage{
##!    inla.jp.define(jp = NULL, ...)
##!}
##!
##!\arguments{
##!  \item{jp}{The \code{jp}-function}
##!  \item{...}{Named list of variables that defines the environment of \code{jp}}
##!}
##!
##!\value{%%
##!  This allows joint priors to be defined in \code{R}.
##!
##!  This function is for internal use only.}
##!\author{Havard Rue \email{hrue@r-inla.org}}


`inla.jp.define` = function(jp = NULL, ...)
{
    stopifnot(!missing(jp))
    args = list(...)
    if (any(names(args) == "")) {
        stop("The '...' argument in 'inla.jp.define()' needs *named* arguments.")
    }
    env = if (length(args) > 0) as.environment(args) else new.env()
    parent.env(env) = .GlobalEnv
    environment(jp) = env
    
    rjp = list(model = jp) ## maybe we need something additional later
    class(rjp) = "inla.jp"

    return (rjp)
}

