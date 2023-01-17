## Export: inla.link.log inla.link.invlog
## Export: inla.link.neglog inla.link.invneglog
## Export: inla.link.logit inla.link.invlogit inla.link.probit
## Export: inla.link.invprobit inla.link.cloglog inla.link.invcloglog
## Export: inla.link.loglog inla.link.invloglog inla.link.tan inla.link.invtan
## Export: inla.link.identity inla.link.invidentity inla.link.invalid
## Export: inla.link.cauchit inla.link.invcauchit
## Export: inla.link.inverse inla.link.invinverse
## Export: inla.link.robit inla.link.invrobit
## Export: inla.link.sn inla.link.invsn

## ! \name{link}
## ! \alias{link}
## ! \alias{inla.link}
## ! \alias{inla.link.log}
## ! \alias{inla.link.invlog}
## ! \alias{inla.link.neglog}
## ! \alias{inla.link.invneglog}
## ! \alias{inla.link.logit}
## ! \alias{inla.link.invlogit}
## ! \alias{inla.link.probit}
## ! \alias{inla.link.invprobit}
## ! \alias{inla.link.cloglog}
## ! \alias{inla.link.invcloglog}
## ! \alias{inla.link.ccloglog}
## ! \alias{inla.link.invccloglog}
## ! \alias{inla.link.loglog}
## ! \alias{inla.link.invloglog}
## ! \alias{inla.link.tan}
## ! \alias{inla.link.invtan}
## ! \alias{inla.link.identity}
## ! \alias{inla.link.invidentity}
## ! \alias{inla.link.invalid}
## ! \alias{inla.link.invinvalid}
## ! \alias{inla.link.cauchit}
## ! \alias{inla.link.invcauchit}
## ! \alias{inla.link.inverse}
## ! \alias{inla.link.invinverse}
## ! \alias{inla.link.robit}
## ! \alias{inla.link.invrobit}
## ! \alias{inla.link.sn}
## ! \alias{inla.link.invsn}
## ! \alias{inla.link.invqpoisson}
## !
## ! \title{Link functions in INLA}
## !
## ! \description{Define link-functions and its inverse}
## !
## ! \usage{
## ! inla.link.log(x, inverse=FALSE)
## ! inla.link.invlog(x, inverse=FALSE)
## ! inla.link.neglog(x, inverse=FALSE)
## ! inla.link.invneglog(x, inverse=FALSE)
## ! inla.link.logit(x, inverse=FALSE)
## ! inla.link.invlogit(x, inverse=FALSE)
## ! inla.link.probit(x, inverse=FALSE)
## ! inla.link.invprobit(x, inverse=FALSE)
## ! inla.link.cloglog(x, inverse=FALSE)
## ! inla.link.invcloglog(x, inverse=FALSE)
## ! inla.link.ccloglog(x, inverse=FALSE)
## ! inla.link.invccloglog(x, inverse=FALSE)
## ! inla.link.loglog(x, inverse=FALSE)
## ! inla.link.invloglog(x, inverse=FALSE)
## ! inla.link.tan(x, inverse=FALSE)
## ! inla.link.invtan(x, inverse=FALSE)
## ! inla.link.cauchit(x, inverse=FALSE)
## ! inla.link.invcauchit(x, inverse=FALSE)
## ! inla.link.identity(x, inverse=FALSE)
## ! inla.link.invidentity(x, inverse=FALSE)
## ! inla.link.inverse(x, inverse=FALSE)
## ! inla.link.invinverse(x, inverse=FALSE)
## ! inla.link.robit(x, df=7, inverse=FALSE)
## ! inla.link.invrobit(x, df=7, inverse=FALSE)
## ! inla.link.sn(x, intercept=0.5, skew=0, a=0, inverse=FALSE)
## ! inla.link.invsn(x, intercept=0.5, skew=0, a=0, inverse=FALSE)
## ! inla.link.invalid(x, inverse=FALSE)
## ! inla.link.invinvalid(x, inverse=FALSE)
## ! inla.link.invqpoisson(x, inverse = FALSE, quantile = 0.5)
## ! }
## !
## ! \arguments{
## !     \item{x}{The argument. A numeric vector.}
## !     \item{df}{The degrees of freedom for the Student-t}
## !     \item{inverse}{Logical. Use the link (\code{inverse=FALSE})
## !                    or its inverse (\code{inverse=TRUE})}
## !     \item{intercept}{The quantile level for the intercept in the Skew-Normal link}
## !     \item{skew}{The skewness in the Skew-Normal.
## !                 Not both of \code{skew} and \code{a} can be given.}
## !     \item{a}{The \code{a}-paraeter in the Skew-Normal.
## !              Not both of \code{skew} and \code{a} can be given.}
## !     \item{quantile}{The quantile level for quantile links}
## !}
## !
## ! \value{Return the values of the link-function or its inverse.}
## ! \note{The \code{inv}-functions are redundant,  as
## !       \code{inla.link.invlog(x) = inla.link.log(x, inverse=TRUE)}
## !       and so on,  but they are simpler to use a arguments
## !       to other functions.}
## ! \author{Havard Rue \email{hrue@r-inla.org}}

`inla.link.cauchit` <- function(x, inverse = FALSE) {
    if (!inverse) {
        return(tan(pi * (x - 0.5)))
    } else {
        return(1.0 / pi * atan(x) + 0.5)
    }
}
`inla.link.invcauchit` <- function(x, inverse = FALSE) {
    return(inla.link.cauchit(x, inverse = !inverse))
}


`inla.link.log` <- function(x, inverse = FALSE) {
    if (!inverse) {
        return(log(x))
    } else {
        return(exp(x))
    }
}
`inla.link.invlog` <- function(x, inverse = FALSE) {
    return(inla.link.log(x, inverse = !inverse))
}

`inla.link.neglog` <- function(x, inverse = FALSE) {
    if (!inverse) {
        return(-log(x))
    } else {
        return(exp(-x))
    }
}
`inla.link.invneglog` <- function(x, inverse = FALSE) {
    return(inla.link.neglog(x, inverse = !inverse))
}

`inla.link.logit` <- function(x, inverse = FALSE) {
    if (!inverse) {
        return(log(x / (1.0 - x)))
    } else {
        return(1.0 / (1.0 + exp(-x)))
    }
}
`inla.link.invlogit` <- function(x, inverse = FALSE) {
    return(inla.link.logit(x, inverse = !inverse))
}

`inla.link.probit` <- function(x, inverse = FALSE) {
    if (!inverse) {
        return(qnorm(x))
    } else {
        return(pnorm(x))
    }
}
`inla.link.invprobit` <- function(x, inverse = FALSE) {
    return(inla.link.probit(x, inverse = !inverse))
}

`inla.link.robit` <- function(x, df = 7, inverse = FALSE) {
    s <- sqrt(df / (df - 2.0))
    if (!inverse) {
        return(qt(x, df = df) * s)
    } else {
        return(pt(x / s, df = df))
    }
}
`inla.link.invrobit` <- function(x, df = 7, inverse = FALSE) {
    return(inla.link.robit(x, df = df, inverse = !inverse))
}

`inla.link.loglog` <- function(x, inverse = FALSE) {
    if (!inverse) {
        return(-log(-log(x)))
    } else {
        return(exp(-exp(-x)))
    }
}
`inla.link.invloglog` <- function(x, inverse = FALSE) {
    return(inla.link.loglog(x, inverse = !inverse))
}

`inla.link.cloglog` <- function(x, inverse = FALSE) {
    if (!inverse) {
        return(log(-log(1 - x)))
    } else {
        return(1.0 - exp(-exp(x)))
    }
}
`inla.link.invcloglog` <- function(x, inverse = FALSE) {
    return(inla.link.cloglog(x, inverse = !inverse))
}

`inla.link.ccloglog` <- function(x, inverse = FALSE) {
    if (!inverse) {
        return(log(-log(x)))
    } else {
        return(exp(-exp(x)))
    }
}
`inla.link.invccloglog` <- function(x, inverse = FALSE) {
    return(inla.link.ccloglog(x, inverse = !inverse))
}

`inla.link.tan` <- function(x, inverse = FALSE) {
    if (!inverse) {
        return(tan(x / 2.0))
    } else {
        return(2.0 * atan(x))
    }
}
`inla.link.invtan` <- function(x, inverse = FALSE) {
    return(inla.link.tan(x, inverse = !inverse))
}

`inla.link.identity` <- function(x, inverse = FALSE) {
    return(x)
}
`inla.link.invidentity` <- function(x, inverse = FALSE) {
    return(inla.link.identity(x, inverse = !inverse))
}

`inla.link.inverse` <- function(x, inverse = FALSE) {
    return(1 / x)
}
`inla.link.invinverse` <- function(x, inverse = FALSE) {
    return(inla.link.inverse(x, inverse = !inverse))
}

`inla.link.invqpoisson` <- function(x, inverse = FALSE, quantile = 0.5) {
    stopifnot(inverse == FALSE)
    return (qgamma(1.0 - quantile, shape = exp(x)+1.0, rate=1))
}

`inla.link.sn` <- function(x, intercept = 0.5, skew = 0, a = NULL, inverse = FALSE) {
    inla.require("sn", stop.on.error = TRUE)

    if (is.null(a)) {
        cache <- inla.pc.sn.cache()
        a <- sign(skew) * cache$pos$a(abs(skew))
    } else {
        stopifnot(all(skew == 0))
    }

    delta <- a / sqrt(1 + a^2)
    omega <- 1 / sqrt(1 - 2 * delta^2 / pi)
    xi <- -omega * delta * sqrt(2 / pi)

    ## skewness
    ## print((4-pi)/2 * (delta*sqrt(2/pi))^3 / ((1-2*delta^2/pi)^(3/2)))

    stopifnot(length(intercept) == 1)
    ## can turn this feature off...
    if (intercept > 0 && intercept < 1.0) {
        intcept <- sn::qsn(intercept, xi = xi, omega = omega, alpha = a)
    } else {
        intcept <- 0
    }

    ## the sn's qsn and psn handle 'a' in a strange way
    if (length(a) == 1) {
        ## then it vectorize...
        if (!inverse) {
            ret <- sn::qsn(x, xi = xi, omega = omega, alpha = a) - intcept
        } else {
            ret <- sn::psn(x + intcept, xi = xi, omega = omega, alpha = a)
        }
    } else {
        ## then it does not vectorize... but does not give any error/warning!
        X <- cbind(x = x, xi = xi, omega = omega, alpha = a)
        ret <- numeric(nrow(X))
        if (!inverse) {
            for (i in 1:nrow(X)) {
                ret[i] <- sn::qsn(X[i, "x"],
                    xi = X[i, "xi"],
                    omega = X[i, "omega"], alpha = X[i, "alpha"]
                ) -
                    intcept
            }
        } else {
            for (i in 1:nrow(X)) {
                ret[i] <- sn::psn(X[i, "x"] + intcept,
                    xi = X[i, "xi"],
                    omega = X[i, "omega"], alpha = X[i, "alpha"]
                )
            }
        }
    }
    return(ret)
}
`inla.link.invsn` <- function(x, intercept = 0.5, skew = 0, a = NULL, inverse = FALSE) {
    return(inla.link.sn(x, intercept = intercept, skew = skew, a = a, inverse = !inverse))
}


## These are the invalid ones
`inla.link.invalid` <- function(x, inverse = FALSE) {
    stop("The invalid link-function is used.")
}
`inla.link.invinvalid` <- function(x, inverse = FALSE) {
    stop("The invinvalid link-function is used.")
}

`inla.link.sn.test` <- function() {
    n <- 4
    a <- runif(n, min = -10, max = 10)
    intercept <- runif(n)
    x <- rnorm(n)

    for (i in 1:n) {
        ans <- c(
            x = x[i], intercept = intercept[i], a = a[i],
            prob = inla.link.sn(x[i], a = a[i], intercept = intercept[i], inverse = TRUE),
            prob.inv = inla.link.sn(inla.link.sn(x[i], a = a[i], intercept = intercept[i], inverse = TRUE),
                a = a[i], intercept = intercept[i], inverse = FALSE
            ),
            "is 0?" = inla.link.sn(inla.link.sn(x[i], a = a[i], intercept = intercept[i], inverse = TRUE),
                a = a[i], intercept = intercept[i], inverse = FALSE
            ) - x[i]
        )
        print(round(dig = 3, ans))
    }

    cat("\n\n")
    n <- 4
    s <- runif(n, min = -0.98, max = 0.98)
    intercept <- runif(n)
    x <- rnorm(n)

    for (i in 1:n) {
        ans <- c(
            x = x[i], intercept = intercept[i], skew = s[i],
            prob = inla.link.sn(x[i], skew = s[i], intercept = intercept[i], inverse = TRUE),
            prob.inv = inla.link.sn(inla.link.sn(x[i], skew = s[i], intercept = intercept[i], inverse = TRUE),
                skew = s[i], intercept = intercept[i], inverse = FALSE
            ),
            "is 0?" = inla.link.sn(inla.link.sn(x[i], skew = s[i], intercept = intercept[i], inverse = TRUE),
                skew = s[i], intercept = intercept[i], inverse = FALSE
            ) - x[i]
        )
        print(round(dig = 3, ans))
    }
}
