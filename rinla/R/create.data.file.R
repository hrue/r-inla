## Nothing to Export

`inla.create.data.file` <- function(
                                    y.orig = NULL,
                                    MPredictor = NULL,
                                    mf = NULL,
                                    scale = NULL,
                                    weights = NULL,
                                    E = NULL,
                                    Ntrials = NULL,
                                    strata = NULL,
                                    lp.scale = NULL,
                                    event = NULL,
                                    family = NULL,
                                    data.dir = NULL,
                                    file = NULL,
                                    debug = FALSE)
{
    my.stop <- function(msg) {
        for (f in c(file, data.dir))
            if (file.exists(f))
                file.remove(f)
        stop(msg)
    }

    y.attr <- attr(y.orig, "inla.ncols", exact = TRUE)
    if (is.null(y.attr)) {
        y.attr <- 0
    }

    if (is.null(y.orig)) {
        y.orig <- c(mf[, 1L])
    } else if (is.inla.surv(y.orig)) {
        y.orig$.special <- NULL
        ## this only applies if we have no cure-model. make sure to remove it
        idx <- which(names(y.orig) == "cure")
        if (length(idx) == 1) {
            if (is.null(y.orig$cure)) {
                y.orig[[idx]] <- NULL
            }
        }
        y.orig <- as.data.frame(unclass(y.orig))
    } else if (is.inla.mdata(y.orig)) {
        y.orig <- as.data.frame(unclass(y.orig))
    } else {
        y.orig <- as.data.frame(y.orig)
    }
    n.data <- dim(y.orig)[1L]
    ind <- seq(0L, n.data - 1L)

    if (debug) {
        cat("inla.create.data.file: n.data = ", n.data, "\n")
    }

    if (!is.null(weights) && !is.function(weights)) {
        if (length(weights) == 1L) {
            weights <- rep(weights, n.data)
        }
        if (length(weights) != n.data) {
            my.stop(paste("Length of 'weights' has to be the same as the length of the response:", length(weights), n.data))
        }
    }

    if (!is.null(lp.scale)) {
        if (length(lp.scale) == 1L) {
            lp.scale <- rep(lp.scale, n.data)
        }
        if (length(lp.scale) != n.data) {
            my.stop(paste("Length of 'lp.scale' has to be the same as the length of the response:", length(lp.scale), n.data))
        }
        lp.scale <- as.numeric(lp.scale)
        lp.scale[is.na(lp.scale)] <- 0
        lp.scale[lp.scale < 0] <- 0
        lp.max <- max(lp.scale)
    }

    if (inla.one.of(family, c(
                                "gaussian",
                                "stdgaussian",
                                "normal",
                                "stdnormal",
                                "lognormal",
                                "t",
                                "sn",
                                "gev",
                                "logistic",
                                "circularnormal",
                                "wrappedcauchy",
                                "iidgamma",
                                "simplex",
                                "gamma",
                                "beta",
                                "tweedie",
                                "fmri"))) {
        if (is.null(scale)) {
            scale <- rep(1.0, n.data)
        }
        if (length(scale) == 1L) {
            scale <- rep(scale, n.data)
        }

        if (length(scale) != n.data) {
            my.stop(paste("Length of scale has to be the same as the length of the response:", length(scale), n.data))
        }

        response <- cbind(ind, scale, y.orig)
        null.dat <- is.na(response[, 3L])
        response <- response[!null.dat, ]
        
        if (any(is.na(response))) {
            my.stop(paste0("family:", family, ". NA's in argument 'scale', are not allowed"))
        }

    } else if (inla.one.of(family, c("tstrata"))) {
        if (is.null(scale)) {
            scale <- rep(1.0, n.data)
        }
        if (length(scale) == 1L) {
            scale <- rep(scale, n.data)
        }

        if (is.null(strata)) {
            strata <- rep(1L, n.data)
        }
        if (length(strata) == 1L) {
            strata <- rep(strata, n.data)
        }

        stopifnot(all(!is.na(strata)))
        stopifnot(all(as.integer(strata) == strata))
        stopifnot(all(strata > 0L))
        if (length(scale) != n.data || length(strata) != n.data) {
            my.stop(paste(
                "Length of scale and strata has to be the same as the length of the response:",
                length(scale), length(strata), n.data
            ))
        }

        ## strata goes from 0L to nstrata-1L,  therefore subtract 1L.
        response <- cbind(ind, scale, as.integer(strata) - 1L, y.orig)
        null.dat <- is.na(response[, 4L])
        response <- response[!null.dat, ]

        if (any(is.na(response))) {
            my.stop(paste0("family:", family, ". NA's in argument 'scale' or 'strata', are not allowed"))
        }

    } else if (inla.one.of(family, c(
                                       "poisson",
                                       "nzpoisson",
                                       "cenpoisson",
                                       "gammacount",
                                       "gpoisson",
                                       "xpoisson",
                                       "zeroinflatedcenpoisson0",
                                       "zeroinflatedcenpoisson1",
                                       "zeroinflatednbinomial0",
                                       "zeroinflatednbinomial1",
                                       "zeroinflatednbinomial2",
                                       "zeroinflatedpoisson0",
                                       "zeroinflatedpoisson1",
                                       "zeroinflatedpoisson2",
                                       "poisson.special1",
                                       "bell"))) {
        if (is.null(E)) {
            E <- rep(1.0, n.data)
        }
        if (length(E) == 1L) {
            E <- rep(E, n.data)
        }

        response <- cbind(ind, E, y.orig)

        if (length(E) != n.data) {
            my.stop(paste("Length of E has to be the same as the length of the response:", length(E), n.data))
        }

        null.dat <- is.na(response[, 3L])
        response <- response[!null.dat, ]

        if (any(is.na(response))) {
            my.stop(paste0("family:", family, ". NA's in argument 'E' are not allowed"))
        }

    } else if (inla.one.of(family, c("nbinomial"))) {
        if (is.null(E)) {
            E <- rep(1.0, n.data)
        }
        if (length(E) == 1L) {
            E <- rep(E, n.data)
        }

        if (is.null(scale)) {
            scale <- rep(1.0, n.data)
        }
        if (length(scale) == 1L) {
            scale <- rep(scale, n.data)
        }

        response <- cbind(ind, E, scale, y.orig)

        if (length(E) != n.data || length(scale) != n.data) {
            my.stop(paste("Length of E and scale has to be the same as the length of the response: E scale data ", length(E), length(scale), n.data))
        }
        null.dat <- is.na(response[, 4L])
        response <- response[!null.dat, ]

        if (any(is.na(response))) {
            my.stop(paste0("family:", family, ". NA's in argument 'E' or 'scale', are not allowed"))
        }

    } else if (inla.one.of(family, c("exponential", "weibull", "loglogistic", "gammajw", "gompertz"))) {
        response <- cbind(ind, y.orig)
        null.dat <- is.na(response[, 2L])
        response <- response[!null.dat, ]
    } else if (inla.one.of(family, c("zeroinflatednbinomial1strata2", "zeroinflatednbinomial1strata3"))) {
        if (is.null(E)) {
            E <- rep(1.0, n.data)
        }
        if (length(E) == 1L) {
            E <- rep(E, n.data)
        }

        if (is.null(strata)) {
            strata <- rep(1L, n.data)
        }
        if (length(strata) == 1L) {
            strata <- rep(strata, n.data)
        }

        stopifnot(all(!is.na(strata)))
        stopifnot(all(strata %in% 1:10))

        response <- cbind(ind, E, strata - 1L, y.orig)

        if (length(E) != n.data) {
            my.stop(paste("Length of E has to be the same as the length of the response:", length(E), n.data))
        }

        if (length(strata) != n.data) {
            my.stop(paste("Length of strata has to be the same as the length of the response:", length(strata), n.data))
        }

        null.dat <- is.na(response[, 4L])
        response <- response[!null.dat, ]

        if (any(is.na(response))) {
            my.stop(paste0("family:", family, ". NA's in argument 'scale' or 'strata', are not allowed"))
        }

    } else if (inla.one.of(family, c("xbinomial"))) {
        if (is.null(scale)) {
            scale <- rep(1.0, n.data)
        }
        if (length(scale) == 1L) {
            scale <- rep(scale, n.data)
        }

        if (is.null(Ntrials)) {
            Ntrials <- rep(1L, n.data)
        }
        if (length(Ntrials) == 1L) {
            Ntrials <- rep(Ntrials, n.data)
        }

        response <- cbind(ind, Ntrials, scale, y.orig)
        null.dat <- is.na(response[, 4L])
        response <- response[!null.dat, ]

        if (any(is.na(response))) {
            my.stop("NA's in argument 'Ntrials' or 'scale',  are not allowed")
        }

    } else if (inla.one.of(family, c("binomial",
                                     "binomialtest",
                                     "betabinomial",
                                     "nbinomial2",
                                     "zeroinflatedbinomial0",
                                     "zeroinflatedbinomial1",
                                     "zeroinflatedbinomial2",
                                     "zeroninflatedbinomial2",
                                     "zeroninflatedbinomial3",
                                     "zeroinflatedbetabinomial0",
                                     "zeroinflatedbetabinomial1",
                                     "zeroinflatedbetabinomial2"))) {
        if (is.null(Ntrials)) {
            Ntrials <- rep(1L, n.data)
        }
        if (length(Ntrials) == 1L) {
            Ntrials <- rep(Ntrials, n.data)
        }

        response <- cbind(ind, Ntrials, y.orig)
        null.dat <- is.na(response[, 3L])
        response <- response[!null.dat, ]

        if (any(is.na(response))) {
            my.stop("NA's in argument 'Ntrials' are not allowed")
        }

    } else if (inla.one.of(family, c("betabinomialna"))) {
        if (is.null(Ntrials)) {
            Ntrials <- rep(1L, n.data)
        }
        if (length(Ntrials) == 1L) {
            Ntrials <- rep(Ntrials, n.data)
        }

        if (is.null(scale)) {
            scale <- rep(1.0, n.data)
        }
        if (length(scale) == 1L) {
            scale <- rep(scale, n.data)
        }

        if (length(scale) != n.data) {
            my.stop(paste("Length of scale has to be the same as the length of the response:", length(scale), n.data))
        }

        response <- cbind(ind, Ntrials, scale, y.orig)
        null.dat <- is.na(response[, 4L])
        response <- response[!null.dat, ]

        if (any(is.na(response))) {
            my.stop(paste0("family:", family, ". NA's in argument 'Ntrials' or 'scale', are not allowed"))
        }

    } else if (inla.one.of(family, c("cbinomial"))) {
        if (!(is.matrix(Ntrials) && all(dim(Ntrials) == c(n.data, 2)))) {
            my.stop(paste("Argument 'Ntrials' for family='cbinomial' must be a", n.data, "x", 2, "-matrix; see the documentation."))
        }
        response <- cbind(ind, Ntrials, y.orig)
        null.dat <- is.na(response[, 4L])
        response <- response[!null.dat, ]

        if (any(is.na(response))) {
            my.stop(paste0("family:", family, ". NA's in argument 'Ntrials' are not allowed"))
        }

    } else if (inla.one.of(family, c("exponentialsurv", "weibullsurv", 
                                     "loglogisticsurv", "qloglogisticsurv", "lognormalsurv",
                                     "gammasurv", "gammajwsurv", "fmrisurv", "gompertzsurv"))) {
        if (!inla.model.properties(family, "likelihood")$survival) {
            my.stop("This should not happen.")
        }
        if (is.null(y.orig$time)) {
            my.stop("Responce does not contain variable `time'.")
        }
        len <- length(y.orig$time)

        if (is.null(y.orig$truncation)) {
            y.orig$truncation <- rep(0, len)
        }
        if (is.null(y.orig$lower)) {
            y.orig$lower <- rep(0, len)
        }
        if (is.null(y.orig$upper)) {
            y.orig$upper <- rep(Inf, len)
        }
        if (is.null(y.orig$event)) {
            y.orig$event <- rep(1, len)
        }

        idx.cure <- grep("^cure[.]?[1-999]*", names(y.orig))
        if (length(idx.cure) > 0) {
            for(i in idx.cure) {
                yy <- y.orig[, i]
                yy[is.na(yy)] <- 0
                y.orig[, i] <- yy
            }
        }

        idx <- !is.na(y.orig$time)
        response <- cbind(
            ind[idx],
            y.orig$event[idx],
            y.orig$truncation[idx],
            y.orig$lower[idx],
            y.orig$upper[idx],
            y.orig[idx, idx.cure], 
            y.orig$time[idx]
        )

        if (any(is.na(response))) {
            my.stop(paste0("family:", family, ". NA's in truncation/event/lower/upper/time is not allowed"))
        }

    } else if (inla.one.of(family, c("stochvol", "stochvolt", "stochvolnig", "stochvolsn", "loggammafrailty",
                                     "iidlogitbeta", "qkumar", "qloglogistic", "gp", "dgp", "pom",
                                     "logperiodogram"))) {
        response <- cbind(ind, y.orig)
        null.dat <- is.na(response[, 2L])
        response <- response[!null.dat, ]
    } else if (inla.one.of(family, c("nmix", "nmixnb"))) {
        if (inla.one.of(family, "nmix")) {
            mmax <- length(inla.model.properties(model = family, section = "likelihood")$hyper)
        } else if (inla.one.of(family, "nmixnb")) {
            ## remove the overdispersion parameter
            mmax <- length(inla.model.properties(model = family, section = "likelihood")$hyper) - 1
        } else {
            my.stop("This should not happen.")
        }

        response <- cbind(IDX = ind, y.orig)
        col.idx <- grep("^IDX$", names(response))
        col.x <- grep("^X[0-9]+", names(response))
        col.y <- grep("^Y[0-9]+", names(response))
        m.x <- length(col.x)
        m.y <- length(col.y)
        stopifnot(m.x >= 1 && m.x <= mmax)

        ## remove entries with NA's in all responses
        na.y <- apply(response[, col.y, drop = FALSE], 1, function(x) all(is.na(x)))
        response <- response[!na.y, , drop = FALSE]

        X <- response[, col.x, drop = FALSE]
        Y <- response[, col.y, drop = FALSE]
        idx <- response[, col.idx, drop = FALSE]
        yfake <- rep(-1, nrow(Y))

        ## replace NA's in the covariates with 0's
        X[is.na(X)] <- 0
        ## augment X til the maximum allowed,  padding with NA's
        X <- cbind(X, matrix(NA, nrow = nrow(response), ncol = mmax - m.x))

        ## sort each row so that the NA's are at the end. Sort numerically the non-NA's as well
        ## although it is not required
        Y <- matrix(c(apply(Y, 1, function(x) c(sort(x[!is.na(x)]), x[is.na(x)]))), ncol = ncol(Y), byrow = TRUE)
        response <- cbind(idx, X, Y, yfake)
    } else if (inla.one.of(family, c("agaussian"))) {
        response <- cbind(IDX = ind, y.orig)
        col.idx <- grep("^IDX$", names(response))
        col.y <- grep("^Y[0-9]+", names(response))
        m.y <- length(col.y)
        stopifnot(m.y == 5L)
        ## remove entries with NA's in all responses
        na.y <- apply(response[, col.y, drop = FALSE], 1, function(x) all(is.na(x)))
        response <- response[!na.y, , drop = FALSE]
        Y <- response[, col.y, drop = FALSE]
        idx <- response[, col.idx, drop = FALSE]
        response <- cbind(idx, Y)
    } else if (inla.one.of(family, c("cenpoisson2"))) {
        if (is.null(E)) {
            E <- rep(1, n.data)
        }
        if (length(E) == 1L) {
            E <- rep(E, n.data)
        }
        response <- cbind(ind, E, y.orig)
        stopifnot(ncol(response) == 5)
        null.dat <- is.na(response[, 3L])
        response <- response[!null.dat,, drop = FALSE]
        colnames(response) <- c("IDX", "E", "Y1", "Y2", "Y3")
        idx.inf <- (is.infinite(response$Y3) | (response$Y3 < 0))
        response[idx.inf, "Y3"] <- -1 ## code for infinite
        idx.inf <- (is.infinite(response$Y2) | (response$Y2 < 0))
        response[idx.inf, "Y2"] <- -1 ## code for infinite
        col.idx <- grep("^IDX$", names(response))
        col.y <- grep("^Y[0-9]+", names(response))
        m.y <- length(col.y)
        stopifnot(m.y == 3L)
        ## remove entries with NA's in all responses
        na.y <- apply(response[, col.y, drop = FALSE], 1, function(x) all(is.na(x)))
        response <- response[!na.y, , drop = FALSE]
        ## format: IDX, E, LOW, HIGH, Y
        response <- cbind(IDX = response$IDX, E = response$E, LOW = response$Y2, HIGH = response$Y3, Y = response$Y1)

        if (any(is.na(response))) {
            my.stop(paste0("family:", family, ". NA's in arguments 'E', 'LOW', 'HIGH', are not allowed"))
        }

    } else if (inla.one.of(family, c("cennbinomial2"))) {

        if (is.null(scale)) {
            scale <- rep(1.0, n.data)
        }
        if (length(scale) == 1L) {
            scale <- rep(scale, n.data)
        }

        if (is.null(E)) {
            E <- rep(1, n.data)
        }
        if (length(E) == 1L) {
            E <- rep(E, n.data)
        }

        response <- cbind(ind, E, scale, y.orig)
        stopifnot(ncol(response) == 6)
        null.dat <- is.na(response[, 4L])
        response <- response[!null.dat,, drop = FALSE ]
        colnames(response) <- c("IDX", "E", "S", "Y1", "Y2", "Y3")
        idx.inf <- (is.infinite(response$Y3) | (response$Y3 < 0))
        response[idx.inf, "Y3"] <- -1 ## code for infinite
        idx.inf <- (is.infinite(response$Y2) | (response$Y2 < 0))
        response[idx.inf, "Y2"] <- -1 ## code for infinite
        col.idx <- grep("^IDX$", names(response))
        col.y <- grep("^Y[0-9]+", names(response))
        m.y <- length(col.y)
        stopifnot(m.y == 3L)
        ## remove entries with NA's in all responses
        na.y <- apply(response[, col.y, drop = FALSE], 1, function(x) all(is.na(x)))
        response <- response[!na.y, , drop = FALSE]
        ## format: IDX, E, S, LOW, HIGH, Y
        response <- cbind(IDX = response$IDX, E = response$E, S = response$S, LOW = response$Y2, HIGH = response$Y3, Y = response$Y1)

        if (any(is.na(response))) {
            my.stop(paste0("family:", family, ". NA's in arguments 'E', 'S', 'LOW', 'HIGH', are not allowed"))
        }

    } else if (inla.one.of(family, c("gaussianjw"))) {

        response <- cbind(ind, y.orig)
        stopifnot(ncol(response) == 5)
        null.dat <- is.na(response[, 2L])
        response <- response[!null.dat,, drop = FALSE]
        null.dat <- is.na(response[, 3L])
        response <- response[!null.dat,, drop = FALSE]
        
        colnames(response) <- c("IDX", "Y1", "Y2", "Y3", "Y4")
        col.y <- 2:5
        na.y <- apply(response[, col.y, drop = FALSE], 1, function(x) all(is.na(x)))
        response <- response[!na.y, , drop = FALSE]
        ## format: IDX, N, DF, VAR, Y
        response <- cbind(IDX = response$IDX, N = response$Y3, DF = response$Y4, VAR = response$Y2, Y = response$Y1)

        if (any(is.na(response))) {
            my.stop(paste0("family:", family, ". NA's in arguments 'N', 'DF', 'VAR', are not allowed"))
        }

    } else if (inla.one.of(family, c("ggaussian"))) {

        response <- cbind(ind, y.orig)
        na.dat <- is.na(response[, 2L])
        response <- response[!na.dat,, drop = FALSE]
        ncovariates <- ncol(response) - 3
        if (!(ncovariates >= 0)) {
            my.stop(paste0("family=", family, ". Number of covariates in simple model is ",  ncovariates,
                           ". Maybe you forgot to add the 's' argument?"))
        }
        if (ncovariates > 0) {
            X <- response[, 4:(4 + ncovariates - 1), drop = FALSE]
            X[is.na(X)] <- 0
            cov.names <- paste0("X", 1:ncovariates)
            colnames(X) <- cov.names
        } else {
            cov.names <- NULL
            X <- NULL
        }

        colnames(response) <- c("IDX", "Y", "SCALE", cov.names)
        ## format: IDX, SCALE, X1, ...XN, Y
        response <- cbind(IDX = response$IDX, SCALE = response$SCALE, X, Y = response$Y)

        if (any(is.na(response))) {
            my.stop(paste0("family:", family, ". NA's in argument 's', are not allowed"))
        }

    } else if (inla.one.of(family, c("ggaussianS"))) {

        response <- cbind(ind, y.orig)
        na.dat <- is.na(response[, 2L])
        response <- response[!na.dat,, drop = FALSE]
        ncovariates <- ncol(response) - 3
        if (!(ncovariates >= 0)) {
            my.stop(paste0("family=", family, ". Number of covariates in simple model is ",  ncovariates,
                           ". Maybe you forgot to add the 'offset' argument?"))
        }
        if (ncovariates > 0) {
            X <- response[, 4:(4 + ncovariates - 1), drop = FALSE]
            X[is.na(X)] <- 0
            cov.names <- paste0("X", 1:ncovariates)
            colnames(X) <- cov.names
        } else {
            cov.names <- NULL
            X <- NULL
        }

        colnames(response) <- c("IDX", "Y", "OFFSET", cov.names)
        ## format: IDX, OFFSET, X1, ...XN, Y
        response <- cbind(IDX = response$IDX, OFFSET = response$OFFSET, X, Y = response$Y)

        if (any(is.na(response))) {
            my.stop(paste0("family:", family, ". NA's in argument 'OFFSET', are not allowed"))
        }

    } else if (inla.one.of(family, c("0poisson", "0poissonS", "0binomial", "0binomialS"))) {

        response <- cbind(ind, y.orig)
        na.dat <- is.na(response[, 2L])
        response <- response[!na.dat,, drop = FALSE]
        ncovariates <- ncol(response) - 3
        stopifnot(ncovariates >= 0)
        if (ncovariates > 0) {
            X <- response[, 4:(4 + ncovariates - 1), drop = FALSE]
            X[is.na(X)] <- 0
            cov.names <- paste0("X", 1:ncovariates)
            colnames(X) <- cov.names
        } else {
            cov.names <- NULL
            X <- NULL
        }

        colnames(response) <- c("IDX", "Y", "EorNtrials", cov.names)
        ## format: IDX, E/Ntrials, X1, ...XN, Y
        response <- cbind(IDX = response$IDX, EorNtrials = response$EorNtrials, X, Y = response$Y)

        if (any(is.na(response))) {
            my.stop(paste0("family:", family, ". NA's in arguments 'E/Ntrials', are not allowed"))
        }

    } else if (inla.one.of(family, c("bgev"))) {

        if (is.null(scale)) {
            scale <- rep(1.0, n.data)
        }
        if (length(scale) == 1L) {
            scale <- rep(scale, n.data)
        }

        if (length(scale) != n.data) {
            my.stop(paste("Length of scale has to be the same as the length of the response:", length(scale), n.data))
        }

        mmax <- length(inla.model.properties(model = family, section = "likelihood")$hyper) - 2L
        response <- cbind(IDX = ind, y.orig)
        col.idx <- grep("^IDX$", names(response))
        col.x <- grep("^X[0-9]+", names(response))
        col.y <- grep("^Y[0-9]+", names(response))
        m.x <- length(col.x)
        m.y <- length(col.y)

        ## remove entries with NA's in all responses
        na.y <- apply(response[, col.y, drop = FALSE], 1, function(x) all(is.na(x)))
        response <- response[!na.y, , drop = FALSE]
        scale <- scale[!na.y]

        X <- response[, col.x, drop = FALSE]
        Y <- response[, col.y, drop = FALSE]
        stopifnot(ncol(Y) == 1)
        idx <- response[, col.idx, drop = FALSE]
        ## replace NA's in the covariates with 0's
        X[is.na(X)] <- 0

        response <- cbind(idx, scale, X, Y)
        if (any(is.na(response))) {
            my.stop(paste0("family:", family, ". NA's in arguments 'scale', are not allowed"))
        }

        ## fix attr, so the order corresponds to (X, Y) and not (Y, X) as in the inla.mdata() input.
        ## y.attr[1] is number of attributes
        stopifnot(y.attr[1] > 0)
        if (y.attr[1] == 1) {
            y.attr <- c(3, y.attr[2], 0, 0)
        } else if (y.attr[1] == 2) {
            y.attr <- c(3, y.attr[2], y.attr[3], 0)
        } else if (y.attr[1] == 3) {
            y.attr <- c(3, y.attr[2], y.attr[3], y.attr[4])
        } else {
            my.stop("FIX THIS with y.attr")
        }
        y.attr <- c(y.attr[1], y.attr[-c(1, 2)], y.attr[2])
    } else {
        my.stop(paste("Family", family, ", not recognised in 'create.data.file.R'"))
    }

    file.data <- inla.tempfile(tmpdir = data.dir)
    inla.write.fmesher.file(as.matrix(response), filename = file.data, debug = debug)
    file.data <- gsub(data.dir, "$inladatadir", file.data, fixed = TRUE)

    file.weights <- inla.tempfile(tmpdir = data.dir)
    if (!is.null(weights) && !is.function(weights)) {
        inla.write.fmesher.file(as.matrix(weights), filename = file.weights, debug = debug)
    } else {
        file.create(file.weights)
    }
    file.weights <- gsub(data.dir, "$inladatadir", file.weights, fixed = TRUE)

    file.lp.scale <- inla.tempfile(tmpdir = data.dir)
    if (!is.null(lp.scale)) {
        inla.write.fmesher.file(as.matrix(lp.scale), filename = file.lp.scale, debug = debug)
    } else {
        file.create(file.lp.scale)
    }
    file.lp.scale <- gsub(data.dir, "$inladatadir", file.lp.scale, fixed = TRUE)

    file.attr <- inla.tempfile(tmpdir = data.dir)
    inla.write.fmesher.file(as.matrix(y.attr, ncol = 1), filename = file.attr, debug = debug)
    file.attr <- gsub(data.dir, "$inladatadir", file.attr, fixed = TRUE)

    return(list(file.data = file.data, file.weights = file.weights, file.attr = file.attr,
                file.lp.scale = file.lp.scale))
}
