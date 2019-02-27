## Export: inla.cut

##! \name{cut}
##! \alias{inla.cut}
##! \alias{cut}
##!
##! \title{Group-wise model criticism using node-splitting}
##!
##! \description{This function performs group-wise, cross-validatory
##! model assessment for an INLA model using so-called node-splitting
##! (Marshall and Spiegelhalter, 2007; Presanis et al, 2013).
##! The user inputs an object of class \code{inla} (i.e. a result
##! of a call to \code{inla()}) as well as  a variable name (\code{split.by}) specifying a grouping:
##! Data points that share the same value of \code{split.by} are in the same group.
##! The function then checks whether each group is an "outlier", or in conflict with
##! the remaining groups, using the methodology described in Ferkingstad et al (2017).
##! The result is a vector containing a p-value for each group, corresponding to a test
##! for each group \code{i}, where the null hypothesis is that group \code{i} is
##! consistent with the other groups except \code{i} (so a small p-value is evidence
##! that the group is an "outlier"). See Ferkingstad et al (2017) for further details.
##! }
##! 
##! \usage{
##!     inla.cut(result, split.by, debug=FALSE)
##! } 
##!
##! \arguments{
##!   \item{result}{An object of class \code{inla}, i.e. a result of a call to \code{inla()}}
##!   \item{split.by}{The name of the variable to group by. Data points that have
##!                   the same value of \code{split.by} are in the same group.}
##!   \item{debug}{Print debugging information if \code{TRUE}, default is \code{FALSE}}
##!  }
##!
##! \value{
##!  A numeric vector of p-values, corresponding to a test
##! for each group \code{i} where the null hypothesis is that group \code{i} is
##! consistent with the other groups except \code{i}. A small p-value for a group
##! indicates that the group is an "outlier" (in conflict with remaining groups).
##!
##!  This function is EXPERIMENTAL!!!
##! }
##!
##! \author{Egil Ferkingstad \email{egil.ferkingstad@gmail.com} and Havard Rue \email{hrue@r-inla.org}}
##!
##! \examples{
##! ## See http://www.r-inla.org/examples/case-studies/ferkingstad-2017 and Ferkingstad et al (2017).
##! }
##!
##! \references{
##! Ferkingstad, E., Held, L. and Rue, H. (2017). Fast and accurate Bayesian
##! model criticism and conflict diagnostics using R-INLA. arXiv preprint
##! arXiv:1708.03272, available at http://arxiv.org/abs/1708.03272.
##! Published in Stat, 6:331-344 (2017).
##! 
##! Marshall, E. C. and Spiegelhalter, D. J. (2007). Identifying outliers
##! in Bayesian hierarchical models: a simulation-based approach.
##! Bayesian Analysis, 2(2):409-444.
##! 
##! Presanis, A. M., Ohlssen, D., Spiegelhalter,
##! D. J., De Angelis, D., et al. (2013). Conflict diagnostics
##! in directed acyclic graphs, with applications in Bayesian evidence synthesis. 
##! Statistical Science, 28(3):376-397.
##! }

inla.cut = function(result, split.by, debug=FALSE)
{
    my.debug = function(...) if (debug) cat("*** inla.cut: ", ... , "\n")

    stopifnot(!missing(result))
    stopifnot(class(result) == "inla")
    stopifnot(!missing(split.by))
    formula <- result$.args$formula
    data <- result$.args$data
    if (!is.data.frame(data)) {
        stop(paste(match.call()[[1]], "is not yet implemented when 'data' is not a data.frame"))
    }

    ## find the name of the response, store the response and remove it from the data.frame.
    ## merge the two dataframes and add the reponse back in
    resp = formula[[2]]
    stopifnot(is.name(resp))
    resp = as.character(resp)
    
    ## extract the variable which define how to split
    n = dim(data)[1]
    split.val = data[, split.by == names(data)]
    split.uval = sort(unique(split.val))
    split.len = length(split.uval)
    stopifnot(split.len > 0)
    
    ## store inla-results
    res = list()
    r = NULL
    p.linpred = c(unlist(
        inla.mclapply(
            seq_len(split.len), 
            FUN = (function(split.idx) {
                idx = (split.val == split.uval[split.idx])
                idx.num = which(idx)
                
                ## First run REP version, everything except part i, "between" group:
                data.rep = data
                data.rep[idx.num, resp] = NA
                
                ## prepare the arguments for inla()
                args = result$.args
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
                cont.pred <- args$control.predictor
                cont.pred$link = 1
                args$control.predictor = cont.pred
                r.rep = do.call("inla", args = args)
                ## this case we do not do
                if (!is.null(r.rep$.args$control.predictor$A)) {
                    stop("A-matrix in the linear predictor is not supported")
                }
                mu.rep = r.rep$summary.lincomb.derived$mean
                sigma.rep = r.rep$misc$lincomb.derived.covariance.matrix
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
                lc.names = paste0("lc.pred", INLA:::inla.num(1:n.pred))
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
                mu.lik = r.lik$summary.lincomb.derived$mean
                sigma.lik = r.lik$misc$lincomb.derived.covariance.matrix
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
                min.eigen.value = 1E-3
                lin.pred.mu = mu.diff
                lin.pred.sigma = sigma.diff
                lin.pred.df = sum(eigen(lin.pred.sigma)$values > min.eigen.value)
                lin.pred.Delta = t(lin.pred.mu) %*% inla.ginv(lin.pred.sigma, tol=min.eigen.value) %*% lin.pred.mu
                pval = 1 - pchisq(as.numeric(lin.pred.Delta), df=lin.pred.df) 
                my.debug(split.idx, "of", split.len,": p-value:", round(pval,4), "df:", lin.pred.df)
                return (pval)
            }))))

    return(p.linpred)
}
