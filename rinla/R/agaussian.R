#' @title Aggregate Gaussian into an equivalent observation
#' 
#' @description
#' Aggregate Gaussians observed with the same mean and precision,
#' into an equivalent triplet, for use with `family="agaussian"`
#' 
#' @param y Repeated observations. If \code{y} is a matrix,  then each row represents
#'            repeated observations. if \code{y} is a list,  then each element of the list is
#'            a vector of repeated observations. If \code{y} is a vector,  then the whole vector
#'            represents repeated observations. The optional scaling \code{s},  must have 
#'            the same format as \code{y},  ie \code{matrix} or \code{vector}.
#'            \code{NA}'s in \code{y} (and \code{s}) are removed and not used or counted.
#'            If \code{s} is given, then the \code{NA}-pattern in \code{y} and \code{s} must be the same.
#' @param s Optional fixed scaling of the precisions. Must be in the same format as \code{y},  and
#'             have the same \code{NA}-pattern. See the documentation for details.
#'
#' @returns The output is a \code{inla.mdata}-object ready for use
#'   with \code{family="agaussian"}. See the example in the documentation.
#'
#' @author Havard Rue \email{hrue@@r-inla.org}
#' @export
#' @aliases inla.agaussian
#' @aliases aggregate.gaussian
#' @aliases agaussian
#' @rdname agaussian
#' 
#' @examples
#'  A = matrix(1:25,5,5)
#'  inla.agaussian(A)
#' 
#'  A[1,-1] = NA
#'  A[2,-(2:3)] = NA
#'  inla.agaussian(A)

`inla.agaussian` <- function(y, s = NULL) {

    agg <- function(i, yy, ss) {
        if (is.matrix(yy)) {
            y <- yy[i, ]
            nn <- length(y)
            s <- if (is.null(ss)) rep(1, nn) else ss[i, ]
        } else if (is.list(yy)) {
            y <- yy[[i]]
            nn <- length(y)
            s <- if (is.null(ss)) rep(1, nn) else ss[[i]]
        } else {
            y <- yy
            nn <- length(y)
            s <- if (is.null(ss)) rep(1, nn) else ss
        }
        
        idx <- !is.na(y)
        y <- y[idx]
        s <- s[idx]
        n <- length(y)
        m <- sum(s)
        y.bar <- sum(s*y)/m
        ldet.s <- 0.5 * sum(log(s))
        v <- sum(s * y^2) / m - y.bar^2
        
        return (list(v, ldet.s, m, n, y.bar))
    }
    
    if (is.matrix(y)) {
        return (inla.mdata(matrix(unlist(lapply(seq_len(nrow(y)), agg, yy = y, ss = s)),
                                  ncol = 5L, byrow = TRUE)))
    } else if (is.list(y)) {
        return (inla.mdata(matrix(unlist(lapply(seq_along(y), agg, yy = y, ss = s)),
                                  ncol = 5L, byrow = TRUE)))
    } else {
        return (inla.mdata(agg(yy = y, ss = s)))
    }
}
        
`inla.agaussian.test` <- function() {
    for(i in 1:10) {
        n <- sample(1:5, 1)
        s <- runif(n)
        y <- rnorm(n)
        mu <- rnorm(1)
        tau <- exp(rnorm(1))
        d <- inla.agaussian(y, s)
        v1 <- sum(dnorm(y, mean = mu, sd = 1/sqrt(s*tau), log = TRUE))
        v2 <- -d$Y4/2*log(2*pi) + 0.5 * d$Y4 * log(tau) + d$Y2 -
            0.5 * tau * d$Y3 * ( (d$Y5 - mu)^2 + d$Y1 )
        cat("Test ", i, "n=", n, "value=", v1, " error=", v2-v1, "\n")
    }
}
