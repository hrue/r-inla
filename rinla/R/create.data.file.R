## Nothing to Export

`inla.create.data.file` = function(
        y.orig = NULL,
        MPredictor = NULL,
        mf = NULL,
        scale = NULL,
        weights = NULL, 
        E = NULL,
        Ntrials = NULL,
        strata = NULL, 
        event = NULL,
        family = NULL,
        data.dir = NULL,
        file = NULL,
        debug = FALSE)
{
    if (is.null(y.orig)) {
        y.orig = c(mf[, 1L])
    } else if (is.inla.surv(y.orig)) {
        y.orig = as.data.frame(unclass(y.orig))
    } else if (is.inla.mdata(y.orig)) {
        y.orig = as.data.frame(unclass(y.orig))
    } else {
        y.orig = as.data.frame(y.orig)
    }
    n.data = dim(y.orig)[1L]
    ind=seq(0L, n.data-1L)

    if (debug) {
        cat("inla.create.data.file: n.data = ", n.data, "\n")
    }
    
    if (!is.null(weights) && !is.function(weights)) {
        if (length(weights) == 1L) {
            weights = rep(weights, n.data)
        }
        if (length(weights) != n.data) {
            file.remove(file)
            file.remove(data.dir)
            stop(paste("Length of 'weights' has to be the same as the length of the response:", length(weights), n.data))
        }
    }

    if (inla.one.of(family, c("gaussian",
                              "normal",
                              "lognormal",
                              "t",
                              "sn",
                              "sn2",
                              "skewnormal",
                              "gev",
                              "logistic",
                              "circularnormal",
                              "wrappedcauchy",
                              "iidgamma",
                              "simplex", 
                              "gamma",
                              "beta"))) {

        if (is.null(scale)) {
            scale = rep(1.0, n.data)
        }
        if (length(scale) == 1L) {
            scale = rep(scale, n.data)
        }

        if (length(scale) != n.data) {
            file.remove(file)
            file.remove(data.dir)
            stop(paste("Length of scale has to be the same as the length of the response:", length(scale), n.data))
        }

        response = cbind(ind, scale, y.orig)
        null.dat = is.na(response[, 3L])
        response = response[!null.dat,]

    } else if (inla.one.of(family, c("tstrata"))) {

        if (is.null(scale)) {
            scale = rep(1.0, n.data)
        }
        if (length(scale) == 1L) {
            scale = rep(scale, n.data)
        }

        if (is.null(strata)) {
            strata = rep(1L, n.data)
        }
        if (length(strata) == 1L) {
            strata = rep(strata, n.data)
        }

        stopifnot(all(!is.na(strata)))
        stopifnot(all(as.integer(strata) == strata))
        stopifnot(all(strata > 0L))
        if (length(scale) != n.data || length(strata) != n.data) {
            file.remove(file)
            file.remove(data.dir)
            stop(paste("Length of scale and strata has to be the same as the length of the response:", length(scale), length(strata), n.data))
        }

        ## strata goes from 0L to nstrata-1L,  therefore subtract 1L.
        response = cbind(ind, scale, as.integer(strata)-1L, y.orig)
        null.dat = is.na(response[, 4L])
        response = response[!null.dat,]

    } else if (inla.one.of(family, c("poisson",
                                     "cenpoisson", 
                                     "contpoisson", 
                                     "qcontpoisson", 
                                     "gpoisson", 
                                     "zeroinflatedpoisson0",
                                     "zeroinflatedpoisson1", 
                                     "zeroinflatedpoisson2", 
                                     "nbinomial",
                                     "zeroinflatednbinomial0",
                                     "zeroinflatednbinomial1",
                                     "zeroinflatednbinomial2"))) {
        
        if (is.null(E)) {
            E = rep(1.0, n.data)
        }
        if (length(E) == 1L) {
            E = rep(E, n.data)
        }

        response = cbind(ind, E, y.orig)

        if (length(E) != n.data) {
            file.remove(file)
            file.remove(data.dir)
            stop(paste("Length of E has to be the same as the length of the response:", length(E), n.data))
        }

        null.dat = is.na(response[, 3L])
        response = response[!null.dat,]

    } else if (inla.one.of(family, c("gammacount", "exponential", "weibull", "loglogistic"))) {

        response = cbind(ind, y.orig)
        null.dat = is.na(response[, 2L])
        response = response[!null.dat,]

    } else if (inla.one.of(family, c("zeroinflatednbinomial1strata2", "zeroinflatednbinomial1strata3"))) {
        if (is.null(E)) {
            E = rep(1.0, n.data)
        }
        if (length(E) == 1L) {
            E = rep(E, n.data)
        }

        if (is.null(strata)) {
            strata = rep(1L, n.data)
        }
        if (length(strata) == 1L) {
            strata = rep(strata, n.data)
        }

        stopifnot(all(!is.na(strata)))
        stopifnot(all(strata %in% 1:10))

        response = cbind(ind, E, strata-1L, y.orig)

        if (length(E) != n.data) {
            file.remove(file)
            file.remove(data.dir)
            stop(paste("Length of E has to be the same as the length of the response:", length(E), n.data))
        }

        if (length(strata) != n.data) {
            file.remove(file)
            file.remove(data.dir)
            stop(paste("Length of strata has to be the same as the length of the response:", length(strata), n.data))
        }

        null.dat = is.na(response[, 4L])
        response = response[!null.dat,]

    } else if (inla.one.of(family,
                           c("binomial",
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
            Ntrials = rep(1L, n.data)
        }
        if (length(Ntrials) == 1L) {
            Ntrials = rep(Ntrials, n.data)
        }

        response = cbind(ind, Ntrials, y.orig)
        null.dat = is.na(response[, 3L])
        response = response[!null.dat,]

    } else if (inla.one.of(family, c("cbinomial"))) {
        if (!(is.matrix(Ntrials) && all(dim(Ntrials) == c(n.data, 2)))) {
            stop(paste("Argument 'Ntrials' for family='cbinomial' must be a", n.data, "x", 2, "-matrix; see the documentation."))
        }
        response = cbind(ind, Ntrials, y.orig)
        null.dat = is.na(response[, 4L])
        response = response[!null.dat,]

    } else if (inla.one.of(family, c("exponentialsurv", "weibullsurv", "weibullcure",
                                     "loglogisticsurv",  "qloglogisticsurv", "lognormalsurv"))) {

        if (!inla.model.properties(family, "likelihood")$survival) {
            file.remove(file)
            file.remove(data.dir)
            stop("This should not happen.")
        }
        if (is.null(y.orig$time)) {
            file.remove(file)
            file.remove(data.dir)
            stop("Responce does not contain variable `time'.")
        }
        len = length(y.orig$time)

        if (is.null(y.orig$truncation)) {
            y.orig$truncation = rep(0, len)
        }
        if (is.null(y.orig$lower)) {
            y.orig$lower = rep(0, len)
        }
        if (is.null(y.orig$upper)) {
            y.orig$upper = rep(Inf, len)
        }
        if (is.null(y.orig$event)) {
            y.orig$event = rep(1, len)
        }

        idx = !is.na(y.orig$time)
        response = cbind(ind[idx],
                         y.orig$event[idx],
                         y.orig$truncation[idx],
                         y.orig$lower[idx],
                         y.orig$upper[idx],
                         y.orig$time[idx])
        if (any(is.na(response))) {
            file.remove(file)
            file.remove(data.dir)
            stop("NA in truncation/event/lower/upper/time is not allowed")
        }

    } else if (inla.one.of(family, c("stochvol", "stochvolt", "stochvolnig", "loggammafrailty",
                                     "iidlogitbeta", "qkumar", "qloglogistic", "gp", "pom"))) {
        response = cbind(ind, y.orig)
        null.dat = is.na(response[, 2L])
        response = response[!null.dat,]

    } else if (inla.one.of(family, c("nmix", "nmixnb"))) {

        ## yes. it must be the same for both model nmix and nmixnb
        mmax = length(inla.model.properties(model="nmix", section="likelihood")$hyper)

        response = cbind(IDX=ind, y.orig)
        col.idx = grep("^IDX$", names(response))
        col.x = grep("^X[0-9]+", names(response))
        col.y = grep("^Y[0-9]+", names(response))
        m.x = length(col.x)
        m.y = length(col.y)
        stopifnot(m.x >= 1 && m.x <= mmax)

        ## remove entries with NA's in all responses
        na.y = apply(response[, col.y, drop=FALSE], 1, function(x) all(is.na(x)))
        response = response[!na.y,, drop=FALSE]

        X = response[, col.x, drop=FALSE]
        Y = response[, col.y, drop=FALSE]
        idx = response[, col.idx, drop=FALSE]
        yfake = rep(-1, nrow(Y))

        ## replace NA's in the covariates with 0's
        X[is.na(X)] = 0
        ## augment X til the maximum allowed,  padding with NA's
        X = cbind(X, matrix(NA, nrow = nrow(response), ncol = mmax -m.x)) 

        ## sort each row so that the NA's are at the end. Sort numerically the non-NA's as well
        ## although it is not required
        Y = matrix(c(apply(Y, 1, function(x) c(sort(x[!is.na(x)]), x[is.na(x)]))), ncol = ncol(Y), byrow=TRUE)
        response = cbind(idx, X, Y, yfake)
        
    } else {

        file.remove(file)
        file.remove(data.dir)
        stop(paste("Family", family, ", not recognised in 'create.data.file.R'"))

    }

    file.data = inla.tempfile(tmpdir=data.dir)
    if (inla.getOption("internal.binary.mode")) {
        inla.write.fmesher.file(as.matrix(response), filename = file.data, debug=debug)
    } else {
        file.create(file.data)
        write(t(response), ncolumns=ncol(response), file=file.data, append=FALSE)
    }
    file.data = gsub(data.dir, "$inladatadir", file.data, fixed=TRUE)

    file.weights = inla.tempfile(tmpdir=data.dir)
    if (!is.null(weights) && !is.function(weights)) {
        if (inla.getOption("internal.binary.mode")) {
            inla.write.fmesher.file(as.matrix(weights), filename = file.weights, debug=debug)
        } else {
            file.create(file.weights)
            write(weights, ncolumns = 1, file=file.weights, append=FALSE)
        }
    } else {
        file.create(file.weights)
    }
    file.weights = gsub(data.dir, "$inladatadir", file.weights, fixed=TRUE)
    
    return(list(file.data = file.data, file.weights = file.weights))
}
