## Export: inla.mdata is.inla.mdata as.inla.mdata
## Export: print!inla.mdata

##!\name{inla.mdata}
##!\alias{inla.mdata}
##!\alias{is.inla.mdata}
##!\alias{as.inla.mdata}
##!\alias{print.inla.mdata}
##!
##!\title{
##!Create an mdata-object for INLA
##!}
##!
##!\description{
##!This defines an mdata-object for matrix valued response-families
##!}
##!\usage{
##!inla.mdata(y, ...)
##!is.inla.mdata(object)
##!as.inla.mdata(object)
##!}
##!
##!\arguments{
##!  \item{y}{The response vector/matrix}
##!  \item{...}{Additional vectors/matrics of same length as \code{y}}
##!  \item{object}{Any \code{R}-object}
##!  \item{x}{An mdata object}
##!  }
##!
##!\value{
##!     An object of class \code{inla.mdata}.  There is method for \code{print}.
##!
##!     \code{is.inla.mdata} returns \code{TRUE} if \code{object}
##!     inherits from class \code{inla.mdata}, otherwise \code{FALSE}.
##!
##!     \code{as.inla.mdata} returns an object of class \code{inla.mdata}
##!}
##!
##!\author{
##! Havard Rue
##!}
##!
##!\seealso{
##!\code{\link{inla}}
##!}

`inla.mdata` = function(y, ...)
{
    y.obj = as.list(as.data.frame(list(y)))
    x.obj = as.list(as.data.frame(list(...)))
    names(y.obj) = paste("Y", 1:length(y.obj), sep="")
    if (length(x.obj) > 0) {
        names(x.obj) = paste("X", 1:length(x.obj), sep="")
    }
    obj = c(y.obj,  x.obj)
    class(obj) = "inla.mdata"
    return (obj)
}

`print.inla.mdata` = function(object, ...)
{
    print(as.data.frame(unclass(object)), ...)
}

`as.inla.mdata` = function(object)
{
    if (is.inla.mdata(object)) {
        return (object)
    }
    object = as.list(as.data.frame(object))
    if (length(object) == 1) {
        names(object) = "Y1"
    } else {
        names(object) = c("Y1",  paste("X", 1:(length(object)-1), sep=""))
    }
    warning("Guess that ncol(response) == 1. Otherwise,  please modify 'names(object)'.")
    class(object) = "inla.mdata"
    return (object)
}

`is.inla.mdata` = function(object)
{
    return (inherits(object, "inla.mdata"))
}
