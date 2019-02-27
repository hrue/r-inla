## Export: merge!inla
## Export: inla.merge

##! \name{inla.merge}
##! \alias{inla.merge}
##! \alias{merge.inla}
##! \title{Merge a mixture of \code{inla}-objects}
##! \description{Merge a mixture of \code{inla}-objects}
##! \usage{
##!     \method{merge}{inla}(x, y, ..., prob = rep(1,  length(loo)), verbose = FALSE)
##!     inla.merge(loo, prob = rep(1,  length(loo)), verbose = FALSE)
##! }
##! \arguments{
##!   \item{x}{An \code{inla}-object to be merged}
##!   \item{y}{An \code{inla}-object to be merged}
##!   \item{...}{Additional \code{inla}-objects to be merged}
##!   \item{loo}{List of \code{inla}-objects to be merged}
##!   \item{prob}{The mixture of (possibly unnormalized) probabilities}
##!   \item{verbose}{Turn on verbose-output or not}
##!  }
##! \value{
##!   A merged \code{inla}-object.
##! }
##! \details{
##!    The function \code{merge.inla} implements method \code{merge} for
##!    \code{inla}-objects. \code{merge.inla} is a wrapper for the function
##!    \code{inla.merge}. The interface is slightly different, \code{merge.inla}
##!    is more tailored for interactive use, whereas \code{inla.merge} is better 
##!    in general code.
##!   
##!    \code{inla.merge} is intented for merging a mixture of \code{inla}-objects,
##!    each run with the same formula and settings, except for a set of
##!    hyperparameters that are fixed to different values. Using this function,
##!    we can then integrate over these hyperparameters using (unnormalized)
##!    integration weights
##!    \code{prob}. The main objects to be merged,  are the summary statistics
##!    and marginal densities (like for hyperparameters, fixed,  random,  etc).
##!    Not all entries in the object can be merged, and by default these
##!    are inheritated from the first object in the list,
##!    while some are just set to  \code{NULL}.
##!    Those objectes that are merged,
##!    will be listed if run with option \code{verbose=TRUE}.
##!  
##!    Note that merging hyperparameter in the user-scale is prone to
##!    discretization error in general, so it is more stable to convert
##!    the marginal of the hyperparameter from the merged internal scale
##!    to the user-scale. (This is not done by this function.)
##! }    
##! \author{Havard Rue \email{hrue@r-inla.org}}
##! \examples{
##! set.seed(123)
##! n = 100
##! y = rnorm(n)
##! y[1:10] = NA
##! x = rnorm(n)
##! z1 = runif(n)
##! z2 = runif(n)*n
##! idx = 1:n
##! idx2 = 1:n
##! lc1 = inla.make.lincomb(idx = c(1, 2, 3))
##! names(lc1) = "lc1"
##! lc2 = inla.make.lincomb(idx = c(0, 1, 2, 3))
##! names(lc2) = "lc2"
##! lc3 = inla.make.lincomb(idx = c(0, 0, 1, 2, 3))
##! names(lc3) = "lc3"
##! lc = c(lc1, lc2, lc3)
##! rr = list()
##! for (logprec in c(0, 1, 2))
##!     rr[[length(rr)+1]] = inla(y ~ 1 + x + f(idx, z1) + f(idx2, z2),
##!              lincomb = lc, 
##!              control.family = list(hyper = list(prec = list(initial = logprec))), 
##!              control.predictor = list(compute = TRUE, link = 1), 
##!              data = data.frame(y, x, idx, idx2, z1, z2))
##! r = inla.merge(rr, prob = seq_along(rr), verbose=TRUE)
##! summary(r)
##!}

`merge.inla` = function(x, y, ...,  prob = rep(1,  length(loo)), verbose = FALSE)
{
    return (inla.merge(loo = list(x, y, ...), prob = prob, verbose = verbose))
}

`inla.merge` = function(loo, prob = rep(1,  length(loo)), verbose = FALSE)
{
    warning("This function is experimental.", immediate. = TRUE)

    verboze = function(..., sep="") {
        if (verbose) {
            cat("inla.merge: ", ..., "\n", sep=sep)
        }
    }

    merge.marginals = function(lom, prob, nx = 128, eps.y = 1e-6) {
        ## compute the (min, median, max) for each marginal, and then the min, median, max of
        ## those again.
        xs = matrix(unlist(lapply(
            lom,
            function(m) {
            xx = m[, "x"]
            return (c(min(xx), median(xx), max(xx)))
        })), nrow = 3)
        x.min = min(xs[1, ])
        x.med = median(xs[2, ])
        x.max = max(xs[3, ])
        ## we guess the transformation... In any case, it should not to be to bad.
        if (x.min > 0 && x.max < 1) {
            ## probability
            m1 = function(x) 1.0/(1.0 + exp(-x))
            m1i = function(x) log(x/(1.0 - x))
        } else if (x.min > 0) {
            ## this is a positive marginal. Choose the scale that make is more symmetric
            if (abs((x.med - x.min)/(x.max - x.min) - 0.5) <
                abs((log(x.med) - log(x.min))/(log(x.max) - log(x.min)) - 0.5)) {
                ## linear
                m1 = m1i = function(x) x
            } else {
                ## log
                m1 = function(x) exp(x)
                m1i = function(x) log(x)
            }
        } else {
            ## use a linear mapping
            m1 = m1i = function(x) x
        }
        xx = m1(seq(m1i(x.min), m1i(x.max), len = nx))
        yy = rep(0, nx)
        for(i in seq_along(lom)) {
            yy = yy + prob[i] * inla.dmarginal(xx, lom[[i]])
        }
        ## remove points not needed
        idx.keep = (yy > eps.y * max(yy))
        marg = inla.smarginal(cbind(x = xx[idx.keep], y = yy[idx.keep]), factor = 2L, keep.type = TRUE)
        return (marg)
    }

    merge.summary = function(los, prob) {
        n = length(los)
        stopifnot(length(prob) == n)
        m = nrow(los[[1]])
        m1 = numeric(m)
        m2 = numeric(m)
        prob = prob / sum(prob)
        for(i in 1:n) {
            m1 = m1 + prob[i] * los[[i]][, "mean"]
            m2 = m2 + prob[i] * (los[[i]][, "sd"]^2 + los[[i]][, "mean"]^2)
        }
        df = data.frame(mean = m1,
                        sd = sqrt(pmax(.Machine$double.eps, m2-m1^2)))
        if (!is.null(los[[1]]$ID)) {
            df$ID = los[[1]]$ID
        }
        rownames(df) = rownames(los[[1]])
        return (df)
    }
    
    stopifnot(length(prob) == length(loo))
    stopifnot(all(prob > 0))
    prob = prob / sum(prob)
    m = length(loo)
    res = loo[[1]]
    
    verboze("Enter with prob = [", round(prob, dig = 3), "], and", m,  "models", sep=" ")

    ## list
    marginals = c(paste0("marginals.", c("hyperpar", "fixed", "lincomb", "lincomb.derived", "linear.predictor")),
                  "internal.marginals.hyperpar")                  
    ## list of list
    marginals2 = paste0("marginals.", c("random", "spde2.blc", "spde3.blc"))
    ## data.frame
    summaries = c(paste0("summary.", c("hyperpar", "fixed", "lincomb", "lincomb.derived", "linear.predictor")), 
                  "internal.summary.hyperpar")
    ## list of data.frame
    summaries2 = paste0("summary.", c("random", "spde2.blc", "spde3.blc"))
    ## items to remove
    remove = c("marginals.fitted.values", "summary.fitted.values", "dic", "cpo", "waic", "po",
               "neffp", "mode", ".args", "model.matrix")
    ## items to remove in $misc
    misc.remove = c("cov.intern", "cov.intern.eigenvalues", "cov.intern.eigenvectors",
                    "lincomb.derived.correlation.matrix", "lincomb.derived.covariance.matrix", 
                    "log.posterior.mode", "stdev.corr.negative", "stdev.corr.negative",
                    "stdev.corr.positive", "theta.mode")

    for(nm in marginals) {
        idx = which(names(res) == nm)
        if (length(idx) > 0) {
            verboze("Merge '$", nm, "'")
            nm = names(res[[idx]])
            res[[idx]] = inla.mclapply(
                seq_along(res[[idx]]),
                function(k) {
                margs = lapply(loo, function(x, kk) x[[idx]][[kk]], kk = k)
                return (merge.marginals(margs, prob))
            })
            names(res[[idx]]) = nm
        }
    }
    
    for (nm in marginals2) {
        idx = which(names(res) == nm)
        if (length(idx) > 0) {
            verboze("Merge '$", nm, "'")
            for(k in seq_along(res[[idx]])) {
                verboze("      '$", nm, "$", names(res[[idx]])[k], "'")
                nm = names(res[[idx]][[k]])
                res[[idx]][[k]] = inla.mclapply(seq_along(res[[idx]][[k]]),
                                                function(kk) {
                    margs = lapply(loo, function(x, k, kkk) x[[idx]][[k]][[kkk]], k = k, kkk = kk)
                    return (merge.marginals(margs, prob))
                })
                names(res[[idx]][[k]]) = nm
            }
        }
    }

    for (nm in summaries) {
        idx = which(names(res) == nm)
        if (length(idx) > 0) {
            verboze("Merge '$", nm, "'")
            if (!is.null(res[[idx]]) && nrow(res[[idx]]) > 0) {
                zum = lapply(loo, function(x) x[[idx]])
                res[[idx]] = merge.summary(zum, prob)
            }
        }
    }

    for (nm in summaries2) {
        idx = which(names(res) == nm)
        if (length(idx) > 0) {
            for(k in seq_along(res[[idx]])) {
                verboze(paste0("      '$", nm, "$", names(res[[idx]])[k], "'"))
                zum = lapply(loo, function(x) x[[idx]][[k]])
                res[[idx]][[k]] = merge.summary(zum, prob)
            }
        }
    }

    ## merge configs, if they are there
    if (!is.null(res$misc$configs)) {
        verboze("Merge '$misc$configs'")
        res$misc$configs$nconfig = m * res$misc$configs$nconfig
        res$misc$configs$config = as.list(seq_len(res$misc$configs$nconfig)) ## create the list to be filled
        res$misc$configs$max.log.posterior = NA ## to be computer later
        count = 1
        for(k in seq_along(loo)) {
            conf = loo[[k]]$misc$configs
            for(kk in seq_len(conf$nconfig)) {
                conf$config[[kk]]$log.posterior = conf$config[[kk]]$log.posterior +
                    conf$max.log.posterior + log(prob[k])
                conf$config[[kk]]$log.posterior.orig = conf$config[[kk]]$log.posterior.orig +
                    conf$max.log.posterior + log(prob[k])
                res$misc$configs$config[[count]] = conf$config[[kk]]
                count = count + 1
            }
        }
        ## rescale the log.posterior's, similar code as in 'inla.collect.misc()'
        res$misc$configs$max.log.posterior = max(sapply(res$misc$configs$config, function(x) x$log.posterior.orig))
        for(k in seq_len(res$misc$configs$nconfig)) {
            res$misc$configs$config[[k]]$log.posterior = res$misc$configs$config[[k]]$log.posterior -
                res$misc$configs$max.log.posterior
            res$misc$configs$config[[k]]$log.posterior.orig = res$misc$configs$config[[k]]$log.posterior.orig - 
                res$misc$configs$max.log.posterior
        }
    }

    verboze(paste0("Merge '$misc$nfunc"))
    res$misc$nfunc = sum(unlist(lapply(loo, function(x) x$misc$nfunc)))

    verboze(paste0("Merge '$cpu.used"))
    res$cpu.used = rowSums(sapply(loo, function(x) x$cpu.used))

    verboze(paste0("Merge '$logfile"))
    res$logfile = c()
    for (k in seq_along(loo)) {
        res$logfile = c(res$logfile,
                        "###",
                        paste0("### CONFIGURATION number ", k, ", prob = ", round(prob[k], dig=5)),
                        "###",
                        loo[[k]]$logfile)
    }

    if (!is.null(res$joint.hyper)) {
        verboze(paste0("Merge '$joint.hyper"))
        res$joint.hyper = data.frame()
        for (k in seq_along(loo)) {
            tmp = loo[[k]]$joint.hyper
            tmp[, ncol(tmp)] = tmp[, ncol(tmp)] + log(prob[k])
            res$joint.hyper = rbind(res$joint.hyper, tmp)
        }
    }

    for (nm in remove) {
        verboze(paste0("Remove '$", nm, "'"))
        idx = which(names(res) == nm)
        if (length(idx) >0) {
            res[[idx]] = NULL
        }
    }

    for (nm in misc.remove) {
        verboze(paste0("Remove '$misc$", nm, "'"))
        idx = which(names(res$misc) == nm)
        if (length(idx) >0) {
            res$misc[[idx]] = NULL
        }
    }

    return (res)
}
