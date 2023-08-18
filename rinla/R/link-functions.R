#' Link functions in INLA
#' 
#' Define link-functions and its inverse
#' 
#' 
#' @aliases link inla.link
#' @param x The argument. A numeric vector.
#' @param df The degrees of freedom for the Student-t
#' @param inverse Logical. Use the link (`inverse=FALSE`) or its inverse
#' (`inverse=TRUE`)
#' @param intercept The quantile level for the intercept in the Skew-Normal
#' link
#' @param skew The skewness in the Skew-Normal.  Only one of `skew` and
#' `a` can be given.
#' @param a The `a`-parameter in the Skew-Normal.  Only one of `skew`
#' and `a` can be given.
#' @param quantile The quantile level for quantile links
#' @return Return the values of the link-function or its inverse.
#' @note The `inv`-functions are redundant, as `inla.link.invlog(x) =
#' inla.link.log(x, inverse=TRUE)` and so on, but they are simpler to use a
#' arguments to other functions.
#' @author Havard Rue \email{hrue@@r-inla.org}
#' 
#' @name link
#' @rdname link-functions
NULL



#' @rdname link-functions
#' @export
`inla.link.cauchit` <- function(x, inverse = FALSE) {
    if (!inverse) {
        return(tan(pi * (x - 0.5)))
    } else {
        return(1.0 / pi * atan(x) + 0.5)
    }
}
#' @rdname link-functions
#' @export
`inla.link.invcauchit` <- function(x, inverse = FALSE) {
    return(inla.link.cauchit(x, inverse = !inverse))
}


#' @rdname link-functions
#' @export
`inla.link.log` <- function(x, inverse = FALSE) {
    if (!inverse) {
        return(log(x))
    } else {
        return(exp(x))
    }
}
#' @rdname link-functions
#' @export
`inla.link.invlog` <- function(x, inverse = FALSE) {
    return(inla.link.log(x, inverse = !inverse))
}

#' @rdname link-functions
#' @export
`inla.link.neglog` <- function(x, inverse = FALSE) {
    if (!inverse) {
        return(-log(x))
    } else {
        return(exp(-x))
    }
}
#' @rdname link-functions
#' @export
`inla.link.invneglog` <- function(x, inverse = FALSE) {
    return(inla.link.neglog(x, inverse = !inverse))
}

#' @rdname link-functions
#' @export
`inla.link.logit` <- function(x, inverse = FALSE) {
    if (!inverse) {
        return(log(x / (1.0 - x)))
    } else {
        return(1.0 / (1.0 + exp(-x)))
    }
}
#' @rdname link-functions
#' @export
`inla.link.invlogit` <- function(x, inverse = FALSE) {
    return(inla.link.logit(x, inverse = !inverse))
}

#' @rdname link-functions
#' @export
`inla.link.probit` <- function(x, inverse = FALSE) {
    if (!inverse) {
        return(qnorm(x))
    } else {
        return(pnorm(x))
    }
}
#' @rdname link-functions
#' @export
`inla.link.invprobit` <- function(x, inverse = FALSE) {
    return(inla.link.probit(x, inverse = !inverse))
}

#' @rdname link-functions
#' @export
`inla.link.robit` <- function(x, df = 7, inverse = FALSE) {
    s <- sqrt(df / (df - 2.0))
    if (!inverse) {
        return(qt(x, df = df) * s)
    } else {
        return(pt(x / s, df = df))
    }
}
#' @rdname link-functions
#' @export
`inla.link.invrobit` <- function(x, df = 7, inverse = FALSE) {
    return(inla.link.robit(x, df = df, inverse = !inverse))
}

#' @rdname link-functions
#' @export
`inla.link.loglog` <- function(x, inverse = FALSE) {
    if (!inverse) {
        return(-log(-log(x)))
    } else {
        return(exp(-exp(-x)))
    }
}
#' @rdname link-functions
#' @export
`inla.link.invloglog` <- function(x, inverse = FALSE) {
    return(inla.link.loglog(x, inverse = !inverse))
}

#' @rdname link-functions
#' @export
`inla.link.cloglog` <- function(x, inverse = FALSE) {
    if (!inverse) {
        return(log(-log(1 - x)))
    } else {
        return(1.0 - exp(-exp(x)))
    }
}
#' @rdname link-functions
#' @export
`inla.link.invcloglog` <- function(x, inverse = FALSE) {
    return(inla.link.cloglog(x, inverse = !inverse))
}

#' @rdname link-functions
#' @export
`inla.link.ccloglog` <- function(x, inverse = FALSE) {
    if (!inverse) {
        return(-log(-log(x)))
    } else {
        return(exp(-exp(-x)))
    }
}
#' @rdname link-functions
#' @export
`inla.link.invccloglog` <- function(x, inverse = FALSE) {
    return(inla.link.ccloglog(x, inverse = !inverse))
}

#' @rdname link-functions
#' @export
`inla.link.tan` <- function(x, inverse = FALSE) {
    if (!inverse) {
        return(tan(x / 2.0))
    } else {
        return(2.0 * atan(x))
    }
}
#' @rdname link-functions
#' @export
`inla.link.invtan` <- function(x, inverse = FALSE) {
    return(inla.link.tan(x, inverse = !inverse))
}

#' @rdname link-functions
#' @export
`inla.link.identity` <- function(x, inverse = FALSE) {
    return(x)
}
#' @rdname link-functions
#' @export
`inla.link.invidentity` <- function(x, inverse = FALSE) {
    return(inla.link.identity(x, inverse = !inverse))
}

#' @rdname link-functions
#' @export
`inla.link.inverse` <- function(x, inverse = FALSE) {
    return(1 / x)
}
#' @rdname link-functions
#' @export
`inla.link.invinverse` <- function(x, inverse = FALSE) {
    return(inla.link.inverse(x, inverse = !inverse))
}

#' @rdname link-functions
#' @export
`inla.link.invqpoisson` <- function(x, inverse = FALSE, quantile = 0.5) {
    stopifnot(inverse == FALSE)
    return (qgamma(1.0 - quantile, shape = exp(x)+1.0, rate=1))
}

#' @rdname link-functions
#' @export
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
#' @rdname link-functions
#' @export
`inla.link.invsn` <- function(x, intercept = 0.5, skew = 0, a = NULL, inverse = FALSE) {
    return(inla.link.sn(x, intercept = intercept, skew = skew, a = a, inverse = !inverse))
}


## These are the invalid ones
#' @rdname link-functions
#' @export
`inla.link.invalid` <- function(x, inverse = FALSE) {
    stop("The invalid link-function is used.")
}
#' @rdname link-functions
#' @export
`inla.link.invinvalid` <- function(x, inverse = FALSE) {
    stop("The invinvalid link-function is used.")
}


## Not exported
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
        print(round(digits = 3, ans))
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
        print(round(digits = 3, ans))
    }
}
