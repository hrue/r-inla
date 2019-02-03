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
    loo = list(x, y, ...)
    
    verboze = function(...) {
        if (verbose) {
            cat("inla.merge: ", ..., "\n")
        }
    }

    merge.marginals = function(lom, prob, nx = 64) {
        n = length(lom)
        eps = 0.001
        x.range = range(unlist(lapply(
            lom, function(m) inla.qmarginal(c(eps, 1-eps), m))))
        xx = seq(x.range[1], x.range[2], len = nx)
        yy = rep(0, nx)
        for(i in 1:n) {
            yy = yy + prob[i] * inla.dmarginal(xx, lom[[i]])
        }
        mm = list(x = xx, y = yy)
        marg = inla.smarginal(mm, factor = 2L)
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

    verboze("Merge '$internal.marginals.hyperpar'")
    for(k in seq_along(res$internal.marginals.hyperpar)) {
        margs = lapply(loo, function(x, k) x$internal.marginals.hyperpar[[k]], k = k)
        res$internal.marginals.hyperpar[[k]] = merge.marginals(margs, prob)
    }
    verboze("Merge '$internal.summary.hyperpar'")
    if (!is.null(res$internal.summary.hyperpar) && nrow(res$internal.summary.hyperpar) > 0) {
        zum = lapply(loo, function(x) x$internal.summary.hyperpar)
        res$internal.summary.hyperpar = merge.summary(zum, prob)
    }

    verboze("Merge '$marginals.hyperpar'")
    for(k in seq_along(res$marginals.hyperpar)) {
        margs = lapply(loo, function(x, k) x$marginals.hyperpar[[k]], k = k)
        res$marginals.hyperpar[[k]] = merge.marginals(margs, prob)
    }
    verboze("Merge '$summary.hyperpar'")
    if (!is.null(res$summary.hyperpar) && nrow(res$summary.hyperpar) > 0) {
        zum = lapply(loo, function(x) x$summary.hyperpar)
        res$summary.hyperpar = merge.summary(zum, prob)
    }

    verboze("Merge '$marginals.random'...")
    for(k in seq_along(res$marginals.random)) {
        verboze(paste0("      '$marginals.random$", names(res$marginals.random)[k], "'"))
        for(kk in seq_along(res$marginals.random[[k]])) {
            margs = lapply(loo, function(x, k, kk) x$marginals.random[[k]][[kk]], k = k, kk = kk)
            res$marginals.random[[k]][[kk]] = merge.marginals(margs, prob)
        }
    }

    verboze("Merge '$summary.random'...")
    for(k in seq_along(res$summary.random)) {
        verboze(paste0("      '$summary.random$", names(res$summary.random)[k], "'"))
        zum = lapply(loo, function(x) x$summary.random[[k]])
        res$summary.random[[k]] = merge.summary(zum, prob)
    }

    verboze("Merge '$marginals.fixed'...")
    for(k in seq_along(res$marginals.fixed)) {
        verboze(paste0("      '$marginals.fixed$", names(res$marginals.fixed)[k], "'"))
        margs = lapply(loo, function(x, k) x$marginals.fixed[[k]], k = k)
        res$marginals.fixed[[k]] = merge.marginals(margs, prob)
    }

    verboze("Merge '$summary.fixed'...")
    if (!is.null(res$summary.fixed) && nrow(res$summary.fixed) > 0) {
        zum = lapply(loo, function(x) x$summary.fixed)
        res$summary.fixed = merge.summary(zum, prob)
    }

    verboze("Merge '$marginals.lincomb'...")
    for(k in seq_along(res$marginals.lincomb)) {
        verboze(paste0("      '$marginals.lincomb$", names(res$marginals.lincomb)[k], "'"))
        margs = lapply(loo, function(x, k) x$marginals.lincomb[[k]], k = k)
        res$marginals.lincomb[[k]] = merge.marginals(margs, prob)
    }

    verboze("Merge '$summary.lincomb'...")
    if (!is.null(res$summary.lincomb) && nrow(res$summary.lincomb) > 0) {
        zum = lapply(loo, function(x) x$summary.lincomb)
        res$summary.lincomb = merge.summary(zum, prob)
    }

    verboze("Merge '$marginals.lincomb.derived'...")
    for(k in seq_along(res$marginals.lincomb.derived)) {
        verboze(paste0("      '$marginals.lincomb.derived$", names(res$marginals.lincomb.derived)[k], "'"))
        margs = lapply(loo, function(x, k) x$marginals.lincomb.derived[[k]], k = k)
        res$marginals.lincomb.derived[[k]] = merge.marginals(margs, prob)
    }

    verboze("Merge '$summary.lincomb.derived'...")
    if (!is.null(res$summary.lincomb.derived) && nrow(res$summary.lincomb.derived) > 0) {
        zum = lapply(loo, function(x) x$summary.lincomb.derived)
        res$summary.lincomb.derived = merge.summary(zum, prob)
    }


    verboze("Merge '$marginals.spde2.blc'...")
    for(k in seq_along(res$marginals.spde2.blc)) {
        verboze(paste0("      '$marginals.spde2.blc$", names(res$marginals.spde2.blc)[k], "'"))
        for(kk in seq_along(res$marginals.spde2.blc[[k]])) {
            margs = lapply(loo, function(x, k, kk) x$marginals.spde2.blc[[k]][[kk]], k = k, kk = kk)
            res$marginals.spde2.blc[[k]][[kk]] = merge.marginals(margs, prob)
        }
    }

    verboze("Merge '$summary.spde2.blc'...")
    for(k in seq_along(res$summary.spde2.blc)) {
        verboze(paste0("      '$summary.spde2.blc$", names(res$summary.spde2.blc)[k], "'"))
        zum = lapply(loo, function(x) x$summary.spde2.blc[[k]])
        res$summary.spde2.blc[[k]] = merge.summary(zum, prob)
    }

    verboze("Merge '$marginals.spde3.blc'...")
    for(k in seq_along(res$marginals.spde3.blc)) {
        verboze(paste0("      '$marginals.spde3.blc$", names(res$marginals.spde3.blc)[k], "'"))
        for(kk in seq_along(res$marginals.spde3.blc[[k]])) {
            margs = lapply(loo, function(x, k, kk) x$marginals.spde3.blc[[k]][[kk]], k = k, kk = kk)
            res$marginals.spde3.blc[[k]][[kk]] = merge.marginals(margs, prob)
        }
    }

    verboze("Merge '$summary.spde3.blc'...")
    for(k in seq_along(res$summary.spde3.blc)) {
        verboze(paste0("      '$summary.spde3.blc$", names(res$summary.spde3.blc)[k], "'"))
        zum = lapply(loo, function(x) x$summary.spde3.blc[[k]])
        res$summary.spde3.blc[[k]] = merge.summary(zum, prob)
    }

    verboze("Merge '$marginals.predictor'...")
    for(k in seq_along(res$marginals.predictor)) {
        margs = lapply(loo, function(x, k) x$marginals.predictor[[k]], k = k)
        res$marginals.random[[k]] = merge.marginals(margs, prob)
    }

    verboze("Merge '$summary.linear.predictor'")
    if (!is.null(res$summary.linear.predictor) && nrow(res$summary.linear.predictor) > 0) {
        zum = lapply(loo, function(x) x$summary.linear.predictor)
        res$summary.linear.predictor = merge.summary(zum, prob)
    }
    
    verboze("Merge '$mlik'...")
    if (!is.null(res$mlik)) {
        val = matrix(unlist(lapply(loo, function(x) x$mlik)), nrow = 2)
        for(i in 1:nrow(res$mlik)) {
            m = max(val[i, ])
            res$mlik[i, 1] = log(sum(prob * (exp(val[i, ] - m)))) + m
        }
    }

    for (nm in c("marginals.fitted.values", "summary.fitted.values",
                 "dic", "cpo", "waic", "po", "neffp")) {
        verboze(paste0("Remove '$", nm, "'..."))
        idx = which(names(res) == nm)
        if (length(idx) >0) {
            res[idx] = NULL
        }
    }

    return (res)
}
