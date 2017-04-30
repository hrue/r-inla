inla.cut = function(formula, data, split.by, correct=FALSE, debug=FALSE, ...)
{
    my.debug = function(...) if (debug) cat("*** inla.cut: ", ... , "\n")

    stopifnot(!missing(formula))
    stopifnot(!missing(data))
    stopifnot(!missing(split.by))

    ## use an INLA internal to get the variable names of the whole model, ie, the fixed effects
    ## + f(idx)'s
    intf = INLA:::inla.interpret.formula(formula)

    ## check that fixed/random effects exists, otherwise you can get an error e.g.
    ## for a model with only fixed effects:
    if (!is.null(intf$fixf)) {
        fixf = as.character(attr(terms.formula(intf$fixf), "variables"))[-c(1:2)]
    } else {
        fixf = NULL
    }
    if (!is.null(intf$randf)) {
        randf = as.character(attr(terms.formula(intf$randf), "variables"))[-c(1:2)]
    } else {
        randf = NULL
    }
    vars = c(fixf, randf)
    
    ## find the name of the response, store the response and remove it from the data.frame.
    ## merge the two dataframes and add the reponse back in
    resp = formula[[2]]
    stopifnot(is.name(resp))
    resp = as.character(resp)
    
    ## extract the variable which define how to split
    n = dim(data)[1]
    split.val = data[, split.by == names(data), drop=FALSE]
    split.uval = sort(unique(split.val))
    split.len = length(split.uval)
    stopifnot(split.len > 0)
    
    ## store the results, 
    result = rep(list(list()), split.len)

    ## store inla-results
    res = list()
    r = NULL
    p.linpred = numeric(split.len)
    for(split.idx in seq_len(split.len)) {
        my.debug(split.idx)
        idx = (split.val == split.uval[split.idx])
        idx.num = which(idx)
        
        ## First run REP version, everything except part i, "between" group:
        data.rep = data
        data.rep[idx.num, resp] = NA
        
        ## prepare the arguments for inla()
        args = list(...)
        args$data = data.rep
        args$formula = formula
        cont.compute = args$control.compute
        if (is.null(cont.compute)) cont.compute = list()
        cont.compute$config = TRUE
        cont.compute$return.marginals = TRUE
        args$control.compute = cont.compute

        ## make linear combinations for calculating linear predictor-based p-values:
        cont.inla = args$control.inla
        cont.inla$lincomb.derived.correlation.matrix = TRUE
        args$control.inla = cont.inla
        n.pred = length(idx.num)
        lc.rep = c()
        for (i in 1:n.pred) {
            pred.idx = rep(NA, n)
            pred.idx[idx.num[i]] = 1
            lci = inla.make.lincomb(Predictor = pred.idx)
            names(lci) = paste0("lc.pred",i)
            lc.rep = c(lc.rep, lci)
        }
        args$lincomb = lc.rep
        r.rep = do.call("inla", args = args)
        ## this case we do not do
        if (!is.null(r.rep$.args$control.predictor$A)) {
            stop("A-matrix in the linear predictor is not supported")
        }

        if (correct) {
            mu.rep =  r.rep$summary.linear.predictor$mean[idx.num]
            for (j in 1:length(mu.rep)) {
                mu.rep[j] = qnorm(inla.pmarginal(mu.rep[j],
                                                 r.rep$marginals.linear.predictor[[idx.num[j]]]))
            }
            sigma.rep = r.rep$misc$lincomb.derived.correlation.matrix
        } else {
            mu.rep = r.rep$summary.lincomb.derived$mean
            sigma.rep = r.rep$misc$lincomb.derived.covariance.matrix
        }
        
        ## Then run LIK version, group i, within-group:
        data.lik = data[idx.num,]
        if (!is.null(args$E)) args$E = args$E[idx.num]
        if (!is.null(args$offset)) args$offset = args$offset[idx.num]
        if (!is.null(args$scale)) args$scale = args$scale[idx.num]
        if (!is.null(args$weights)) args$weights = args$weights[idx.num]
        if (!is.null(args$Ntrials)) args$Ntrials = args$Ntrials[idx.num]
        if (!is.null(args$strata)) args$strata = args$strata[idx.num]
        if (!is.null(args$link.covariates)) args$link.covariates = args$link.covariates[idx.num, ]

        args$data = data.lik
        lc.lik = c()
        ## make new linear combinations:
        lc.names = paste0("lc.pred", inla.num(1:n.pred))
        for (i in 1:n.pred) {
            pred.idx = rep(NA, nrow(data.lik))
            pred.idx[i] = 1
            lci = inla.make.lincomb(Predictor = pred.idx)
            names(lci) = lc.names[i]
            lc.lik = c(lc.lik, lci)
        }

        args$lincomb = lc.lik
        ## set prior of hyperpar to posterior from REP version:
        args$control.update = list(result = r.rep)
        r.lik = do.call("inla", args = args)

        if (correct) {
            mu.lik = r.lik$summary.linear.predictor$mean
            for (j in 1:length(mu.lik)) {
                mu.lik[j] = qnorm(inla.pmarginal(mu.lik[j],
                                                 r.lik$marginals.linear.predictor[[j]]))
            }
            sigma.lik = r.lik$misc$lincomb.derived.correlation.matrix
        } else {
            mu.lik = r.lik$summary.lincomb.derived$mean
            sigma.lik = r.lik$misc$lincomb.derived.covariance.matrix
        }
        mu.diff = mu.rep - mu.lik
        sigma.diff = sigma.rep + sigma.lik
        respIdx = which(names(data.lik)==resp)
        curr.data = data.lik[,-respIdx]
        curr.data[is.na(curr.data)] = 0
        ## curr.data: change any factors to numeric:
        for (kk in 1:ncol(curr.data)) {
            if (is.factor(curr.data[,kk])) {
                curr.data[,kk] = as.numeric(as.character(curr.data[,kk]))
            }
        }
        
        ## Calculate p-values from linear predictor:
        lin.pred.mu = mu.diff
        lin.pred.sigma = sigma.diff
        lin.pred.df = sum(eigen(lin.pred.sigma)$values > 1e-3)
        lin.pred.Delta = t(lin.pred.mu) %*% inla.ginv(lin.pred.sigma, tol=1e-3) %*% lin.pred.mu
        p.linpred[split.idx] = 1 - pchisq(as.numeric(lin.pred.Delta), df=lin.pred.df) 
        my.debug(split.idx, "of", split.len,": p-value:", round(p.linpred[split.idx],4),
                 "df:", lin.pred.df)
    }

    return(p.linpred)
}
