#' Computes the mean and stdev for the spline from `scopy`
#' 
#' This function computes the mean and stdev for the spline function that is
#' implicite from an `scopy` model component
#' 
#' 
#' @aliases inla.scopy.summary scopy.summary
#' @param result An `inla`-object, ie the output from an `inla()`
#' call
#' @param name The name of the `scopy` model component see `?INLA::f`
#' and argument `scopy`
#' @param mean.value In case where the mean of the spline is fixed
#' and not estimated, you have to give it here
#' @param slope.value In case where the slope of the spline is fixed
#' and not estimated, you have to give it here
#' @param by The resolution of the results, in the scale where 
#' `diff(range(locations))` is 1
#' @param range The range of the locations, as `c(from, to)`
#' @param debug If `TRUE` then enable some debug output
#' @return A `data.frame` with locations, mean and stdev.  If `name`
#' is not found, NULL is returned.
#' @author Havard Rue \email{hrue@r-inla.org}
#' @examples
#' 
#'  ## see example in inla.doc("scopy")
#'  
#' @name scopy.summary
#' @export

`inla.scopy.summary` <- function(result, name, mean.value = NULL, slope.value = NULL,
                                 by = 0.01, range = c(0, 1), debug = FALSE)
{
    stopifnot(!missing(result) && inherits(result, "inla"))
    if (is.null(result$misc$configs)) {
        stop("you need an inla-object computed with option 'control.compute=list(config = TRUE)'.")
    }
    stopifnot(range[1] < range[2])

    cs <- result$misc$configs
    ld <- numeric(cs$nconfig)
    for (i in 1:cs$nconfig) {
        ld[i] <- cs$config[[i]]$log.posterior
    }
    p <- exp(ld - max(ld))
    p <- p / sum(p)

    k <- 2:length(inla.models()$latent$scopy$hyper)
    nms <- c(paste0("Beta", 0, " for ", name, " (scopy mean)"),
             paste0("Beta", 1, " for ", name, " (scopy slope)"),
             paste0("Beta", k, " for ", name, " (scopy theta)"))

    idx <- c()
    theta <- names(cs$config[[1]]$theta)
    for(nm in nms) {
        j <- which(nm == theta)
        stopifnot(length(j) <= 1)
        if (length(j) == 1) {
            idx <- c(idx, j)
        }
    }
    if (length(idx) == 0 && is.null(slope.value) && is.null(mean.value)) {
        return (NULL)
    }

    fixed.mean <- FALSE
    fixed.slope <- FALSE
    if ((length(theta) > 0) && (theta[idx[1]] %in% nms[1])) {
        ## then mean is there,  all ok
        stopifnot(missing(mean.value))
    } else {
        if (missing(mean.value))
            mean.value <- 0
        fixed.mean <- TRUE
    }

    if ((length(theta) > 0) &&
        ((theta[idx[1]] %in% nms[2]) ||
         (theta[idx[2]] %in% nms[2]))) {
        ## then slope is there,  all ok
        stopifnot(missing(slope.value))
    } else {
        if (missing(slope.value))
            slope.value <- 0
        fixed.slope <- TRUE
    }

    n <- length(idx) + as.numeric(fixed.slope)+ as.numeric(fixed.mean)
    prop <- inla.scopy.define(n)
    eps <- 1e-6
    stopifnot(diff(range) >= eps * n)
    by <- max(eps, min(1.0, by))
    by <- diff(range) * by / (1 - by)
    xx <- seq(range[1], range[2], by = by)
    xx.loc <- range[1] + diff(range) * seq(0, 1.0, len = n)
    xx.std <- seq(-0.5, 0.5, len = n)

    ex <- numeric(length(xx))
    exx <- numeric(length(xx))
    
    if (debug) {
        print(paste0("length of spline is ", n))
        print(paste0("fixed.mean ", fixed.mean))
        print(paste0("fixed.slope ", fixed.slope))
        print(paste0("cs$nconfig ", cs$nconfig))
    }

    if (cs$nconfig > 1) {
        for(i in 1:cs$nconfig) {
            th <- cs$config[[i]]$theta[idx]
            tth <- c()
            if (fixed.mean) {
                tth <- c(tth, mean.value)
            } else {
                tth <- c(tth, th[1])
                th <- th[-1]
            }
            if (fixed.slope) {
                tth <- c(tth, slope.value)
            } else {
                tth <- c(tth, th[1])
                th <- th[-1]
            }
            tth <- c(tth, th)
            spline.vals <- prop$W %*% tth
            fun <- splinefun(xx.loc, spline.vals, method = "natural")
            vals <- fun(xx)
            ex <- ex + p[i] * vals
            exx <- exx + p[i] * vals^2
        }
        m <- ex 
        s <- sqrt(pmax(0, exx - ex^2))
    } else {
        th <- cs$config[[1]]$theta[idx]
        tth <- c()
        if (fixed.mean) {
            tth <- c(tth, mean.value)
        } else {
            tth <- c(tth, th[1])
            th <- th[-1]
        }
        if (fixed.slope) {
            tth <- c(tth, slope.value)
        } else {
            tth <- c(tth, th[1])
            th <- th[-1]
        }
        tth <- c(tth, th)
        spline.vals <- prop$W %*% tth
        fun <- splinefun(xx.loc, spline.vals, method = "natural")
        m <- fun(xx)
        s <- rep(0, length(m))
    }

    return(data.frame(x = xx, mean = m, sd = s))
}

`inla.scopy.define` <- function(n = 5L) {
    stopifnot(n == 2 || n >= 5)
    if (n == 2) {
        VV = cbind(1, c(-0.5, 0.5))
    } else {
        Q <- inla.rw2(n, scale.model = TRUE)
        e <- eigen(Q)
        idx.remove <- (n-1):n
        V <- e$vectors[, -idx.remove, drop = FALSE]
        lambda <- e$values[ -idx.remove]
        ## this is the sqrt(covariance.matrix) which makes the remaining easy. put them in reverse
        ## order, which make more sense
        VV <- matrix(0, n, n)
        VV[, 1:2] <- cbind(rep(1, n), seq(-0.5, 0.5, len = n))
        m <- n - 2
        for (i in seq_len(m)) {
            j <- m - i + 3
            VV[, j] <- V[, i] * sqrt(1 / lambda[i])
        }
    }
    return (list(n = n, W = VV))
}
