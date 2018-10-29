## Export: inla.cpo

##!\name{inla.cpo}
##!\alias{inla.cpo}
##!\alias{cpo.inla}
##!
##!\title{Improved estimates for the CPO/PIT-values}
##!\description{
##! Improve the estimates of the  CPO/PIT-values
##! be recomputing the model-fit by removing data-points. 
##!}
##!\usage{
##!inla.cpo(result,
##!         force = FALSE,
##!         verbose = TRUE,
##!         recompute.mode = TRUE)
##!}
##!
##!\arguments{
`inla.cpo` = function(
        ##!\item{result}{An object of class \code{inla}, ie a result
        ##!of a call to \code{inla()}}
        result,

        ##!\item{force}{If TRUE,  then recompute all CPO/PIT values and not just 
        ##!those with \code{result$cpo$failure > 0}.}
        force = FALSE,

        ##!\item{verbose}{Run in verbose mode?}
        verbose = TRUE,

        ##!\item{recompute.mode}{Should be mode (and the integration points) be recomputed when a data-point is removed or not?}
        recompute.mode = TRUE)
{
    ##!}
    ##!\value{
    ##!The object returned is the same as  \code{result} but the
    ##!new improved estimates of the CPO/PIT values replaced.}
    ##!\author{Havard Rue \email{hrue@r-inla.org}}
    ##!\examples{
    ##!n = 10
    ##!y = rnorm(n)
    ##!r = inla(y ~ 1, data = data.frame(y), control.compute = list(cpo=TRUE))
    ##!\donttest{
    ##!rr = inla.cpo(r, force=TRUE)
    ##!}
    ##!}
    ##!\seealso{\code{\link{inla}}}

    stopifnot(!missing(result))
    if (class(result) != "inla") {
        warning("The function 'inla.cpo()' has changed; see ?inla.cpo for futher details.\n  The old version is available as 'inla.cpo.OLD()' for a while.")
    }
    stopifnot(class(result) == "inla")

    ## if there is no cpo, then done
    if (is.null(result$cpo))
        return(result)

    ## loop over those with failure > 0
    if (!force) {
        idx.fail = which(result$cpo$failure > 0)
    } else {
        idx.fail = 1:length(result$cpo$failure)
    }

    result.new = result
    if (length(idx.fail) > 0L) {
        k=1L
        cpo.old = result$cpo$cpo[idx.fail]
        pit.old = result$cpo$pit[idx.fail]

        if (verbose) {
            cat("Compute new CPO/PIT values manually, for", length(idx.fail), "cases...\n")
        }

        res = inla.mclapply(
                        idx.fail,
                        function(idx, result) {
                            result$.args$control.expert = list(cpo.manual = TRUE, cpo.idx = idx)
                            result$.args$control.mode = list(
                                    result = NULL,
                                    fixed = result$.args$control.mode$fixed,
                                    theta = result$mode$theta,
                                    x = result$mode$x,
                                    restart = recompute.mode)
                            result$.args$control.compute$dic = result$.args$control.compute$waic = FALSE
                            rr = inla.rerun(result, plain = TRUE)$cpo
                            return (list(cpo = rr$cpo[idx], pit = rr$pit[idx], failure = rr$failure[idx]))
                        },  result = result, mc.cores = result$.args$num.threads)
        
        result.new$cpo$cpo[idx.fail] = unlist(lapply(res, function(xx) return (xx$cpo)))
        result.new$cpo$pit[idx.fail] = unlist(lapply(res, function(xx) return (xx$pit)))
        result.new$cpo$failure[idx.fail] = unlist(lapply(res, function(xx) return (xx$failure)))
        
        if (verbose) {
            dig = getOption("digits")
            options(digits = 6)
            print(cbind(index = idx.fail,
                        cpo.old = cpo.old, cpo.new = result.new$cpo$cpo[idx.fail],
                        pit.old = pit.old, pit.new = result.new$cpo$pit[idx.fail]))
            cat("\nThe retured result contain the new values.\n")
            options(digits = dig)
        }
    }

    return(result.new)
}
