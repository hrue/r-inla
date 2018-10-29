## Export: inla.link.log inla.link.invlog
## Export: inla.link.neglog inla.link.invneglog
## Export: inla.link.logit inla.link.invlogit inla.link.probit
## Export: inla.link.invprobit inla.link.cloglog inla.link.invcloglog 
## Export: inla.link.loglog inla.link.invloglog inla.link.tan inla.link.invtan
## Export: inla.link.identity inla.link.invidentity inla.link.invalid
## Export: inla.link.cauchit inla.link.invcauchit
## Export: inla.link.inverse inla.link.invinverse

##! \name{link}
##! \alias{link}
##! \alias{inla.link}
##! \alias{inla.link.log}
##! \alias{inla.link.invlog}
##! \alias{inla.link.neglog}
##! \alias{inla.link.invneglog}
##! \alias{inla.link.logit}
##! \alias{inla.link.invlogit}
##! \alias{inla.link.probit}
##! \alias{inla.link.invprobit}
##! \alias{inla.link.cloglog}
##! \alias{inla.link.invcloglog}
##! \alias{inla.link.loglog}
##! \alias{inla.link.invloglog}
##! \alias{inla.link.tan}
##! \alias{inla.link.invtan}
##! \alias{inla.link.identity}
##! \alias{inla.link.invidentity}
##! \alias{inla.link.invalid}
##! \alias{inla.link.invinvalid}
##! \alias{inla.link.cauchit}
##! \alias{inla.link.invcauchit}
##! \alias{inla.link.inverse}
##! \alias{inla.link.invinverse}
##! 
##! \title{Link functions in INLA}
##! 
##! \description{Define link-functions and its inverse}
##!
##! \usage{
##! inla.link.log(x, inverse=FALSE)
##! inla.link.invlog(x, inverse=FALSE)
##! inla.link.neglog(x, inverse=FALSE)
##! inla.link.invneglog(x, inverse=FALSE)
##! inla.link.logit(x, inverse=FALSE)
##! inla.link.invlogit(x, inverse=FALSE)
##! inla.link.probit(x, inverse=FALSE)
##! inla.link.invprobit(x, inverse=FALSE)
##! inla.link.cloglog(x, inverse=FALSE)
##! inla.link.invcloglog(x, inverse=FALSE)
##! inla.link.loglog(x, inverse=FALSE)
##! inla.link.invloglog(x, inverse=FALSE)
##! inla.link.tan(x, inverse=FALSE)
##! inla.link.invtan(x, inverse=FALSE)
##! inla.link.cauchit(x, inverse=FALSE)
##! inla.link.invcauchit(x, inverse=FALSE)
##! inla.link.identity(x, inverse=FALSE)
##! inla.link.invidentity(x, inverse=FALSE)
##! inla.link.inverse(x, inverse=FALSE)
##! inla.link.invinverse(x, inverse=FALSE)
##! inla.link.invalid(x, inverse=FALSE)
##! inla.link.invinvalid(x, inverse=FALSE)
##! }
##! 
##! \arguments{
##!     \item{x}{The argument. A numeric vector.}
##!     \item{inverse}{Logical. Use the link (\code{inverse=FALSE})
##!                    or its inverse (\code{inverse=TRUE})}
##!}
##! 
##! \value{Return the values of the link-function or its inverse.}
##! \note{The \code{inv}-functions are redundant,  as 
##!       \code{inla.link.invlog(x) = inla.link.log(x, inverse=TRUE)}
##!       and so on,  but they are simpler to use a arguments
##!       to other functions.}
##! \author{Havard Rue \email{hrue@r-inla.org}}

`inla.link.cauchit` = function(x, inverse = FALSE)
{
    if (!inverse) {
        return (tan(pi * (x - 0.5)))
    } else {
        return (1.0/pi * atan(x) + 0.5)
    }
}
`inla.link.invcauchit` = function(x, inverse = FALSE)
{
    return (inla.link.cauchit(x, inverse = !inverse))
}


`inla.link.log` = function(x, inverse = FALSE)
{
    if (!inverse) {
        return (log(x))
    } else {
        return (exp(x))
    }
}
`inla.link.invlog` = function(x, inverse = FALSE)
{
    return (inla.link.log(x, inverse = !inverse))
}

`inla.link.neglog` = function(x, inverse = FALSE)
{
    if (!inverse) {
        return (-log(x))
    } else {
        return (exp(-x))
    }
}
`inla.link.invneglog` = function(x, inverse = FALSE)
{
    return (inla.link.neglog(x, inverse = !inverse))
}

`inla.link.logit` = function(x, inverse = FALSE)
{
    if (!inverse) {
        return (log(x/(1.0-x)))
    } else {
        return (1.0/(1.0+exp(-x)))
    }
}
`inla.link.invlogit` = function(x, inverse = FALSE)
{
    return (inla.link.logit(x, inverse = !inverse))
}

`inla.link.probit` = function(x, inverse = FALSE)
{
    if (!inverse) {
        return (qnorm(x))
    } else {
        return (pnorm(x))
    }
}
`inla.link.invprobit` = function(x, inverse = FALSE)
{
    return (inla.link.probit(x, inverse = !inverse))
}

`inla.link.loglog` = function(x, inverse = FALSE)
{
    if (!inverse) {
        return (-log(-log(x)))
    } else {
        return (exp(-exp(-x)))
    }
}
`inla.link.invloglog` = function(x, inverse = FALSE)
{
    return (inla.link.loglog(x, inverse = !inverse))
}

`inla.link.cloglog` = function(x, inverse = FALSE)
{
    if (!inverse) {
        return (log(-log(1-x)))
    } else {
        return (1.0 - exp(-exp(x)))
    }
}
`inla.link.invcloglog` = function(x, inverse = FALSE)
{
    return (inla.link.cloglog(x, inverse = !inverse))
}

`inla.link.tan` = function(x, inverse = FALSE)
{
    if (!inverse) {
        return (tan(x/2.0))
    } else {
        return (2.0*atan(x))
    }
}
`inla.link.invtan` = function(x, inverse = FALSE)
{
    return (inla.link.tan(x, inverse = !inverse))
}

`inla.link.identity` = function(x, inverse = FALSE)
{
    return (x)
}
`inla.link.invidentity` = function(x, inverse = FALSE)
{
    return (x)
}

`inla.link.inverse` = function(x, inverse = FALSE)
{
    return (1/x)
}
`inla.link.invinverse` = function(x, inverse = FALSE)
{
    return (1/x)
}


## These are the invalid one
`inla.link.invalid` = function(x, inverse = FALSE)
{
    stop("The invalid link-function is used.")
}
`inla.link.invinvalid` = function(x, inverse = FALSE)
{
    stop("The invinvalid link-function is used.")
}
