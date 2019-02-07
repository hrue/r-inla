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
##!    \code{inla.merge}. The interface is slightly different.
##!   
##!    \code{inla.merge} is intented for merging a mixture of \code{inla}-objects,
##!    each run with the same formula and settings, except for a set of
##!    hyperparameters that are fixed to different values. Using this function,
##!    we can then integrate over these hyperparameters using (unnormalized)
##!    integration weights
##!    \code{prob}. Not all entries in the object can be merged, and some are
##!    inheritated from the first object in the list,  others set to \code{NULL}.
##!    Those objectes that are merged,
##!    will be listed with option \code{verbose=TRUE}.
##!  
##!    Note that merging hyperparameter in the user-scale can be prone to
##!    discretization error, so more stable results can be found converting
##!    the marginal of the hyperparameter in the internal scale to the user-scale.
##! }    
##! \author{Havard Rue \email{hrue@r-inla.org}}

`merge.inla` = function(x, y, ...,  prob = rep(1,  length(loo)), verbose = FALSE)
{
    return (inla.merge(loo = list(x, y, ...), prob = prob, verbose = verbose))
}

`inla.merge` = function(loo, prob = rep(1,  length(loo)), verbose = FALSE)
{
    verboze = function(...) {
        if (verbose) {
            cat("inla.merge: ", ..., "\n", sep="")
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
            ## this is positive marginal. Choose the scale that make is more symmetric
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
        marg = inla.smarginal(list(x = xx[idx.keep], y = yy[idx.keep]), factor = 2L)
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
    
    verboze("enter with ", m, " models")
    verboze("enter with prob = ", round(prob, dig = 3))

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
                     
    remove = c("marginals.fitted.values", "summary.fitted.values", "dic", "cpo", "waic", "po", "neffp")


    for(nm in marginals) {
        idx = which(names(res) == nm)
        if (length(idx) > 0) {
            verboze("Merge '$", nm, "'")
            dummy = INLA:::inla.mclapply(seq_along(res[[idx]]),
                                         function(k) {
                               margs = lapply(loo, function(x, kk) x[[idx]][[kk]], kk = k)
                               res[[idx]][[k]] = merge.marginals(margs, prob)
                               return (NULL)
                           })
        }
    }
    
    for (nm in marginals2) {
        idx = which(names(res) == nm)
        if (length(idx) > 0) {
            verboze("Merge '$", nm, "'")
            for(k in seq_along(res[[idx]])) {
                verboze("      '$", nm, "$", names(res[[idx]])[k], "'")
                dummy = INLA:::inla.mclapply(seq_along(res[[idx]][[k]]),
                                             function(kk) {
                                   margs = lapply(loo, function(x, k, kkk) x[[idx]][[k]][[kkk]], k = k, kkk = kk)
                                   res[[idx]][[k]][[kk]] = merge.marginals(margs, prob)
                                   return (NULL)
                               })
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

    for (nm in remove) {
        verboze(paste0("Remove '$", nm, "'"))
        idx = which(names(res) == nm)
        if (length(idx) >0) {
            res[[idx]] = NULL
        }
    }

    return (res)
}
