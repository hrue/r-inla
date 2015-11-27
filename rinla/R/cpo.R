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
'inla.cpo' = function(
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
    ##!\author{Havard Rue \email{hrue@math.ntnu.no}}
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

'inla.cpo.OLD' = function(...,  force.cpo.manual = FALSE, be.verbose = TRUE, recompute.mode = TRUE)
{
    ## evaluate inla() using the same arguments, but compute the
    ## cpo-values manually for those with a 'failure'

    lib = "INLA::"
    arg.char = as.character(as.expression(match.call()))
    just.args = gsub(paste("^(", lib, ")?inla.cpo[(]", sep=""), "", arg.char)
    just.args = gsub("[)]$", "", just.args)

    ## remove the local arguments
    just.args = gsub(",[ ]*force\\.cpo\\.manual[ ]*=[^,]*", "", just.args)
    just.args = gsub(",[ ]*be\\.verbose[ ]*=[^,]*", "", just.args)
    just.args = gsub(",[ ]*recompute\\.mode[ ]*=[^,]*", "", just.args)

    ## call the standard inla()
    r = inla.eval(paste(lib, "inla(", just.args, ")", sep=""))

    ## if there is no cpo, then done
    if (is.null(r$cpo))
        return(r)

    ## loop over those with failure > 0
    if (!force.cpo.manual)
        idx.fail = which(r$cpo$failure > 0)
    else
        idx.fail = 1:length(r$cpo$failure)
    if (length(idx.fail) > 0) {
        cpo.old = r$cpo$cpo[idx.fail]
        pit.old = r$cpo$pit[idx.fail]
        k=1
        for(idx in idx.fail) {
            if (be.verbose) {
                cat("\r                                       ", inla.tictac(k))
                cat("\rCompute manually CPO/PIT for ", length(idx.fail), " nodes:\r", sep="")
                k = k+1
            }
            
            argument = paste(sep="", lib, "inla(", just.args, ",",
                    "control.expert = list(cpo.manual = TRUE, cpo.idx =", idx, "),",
                    "control.mode = list(result = r, restart=", inla.ifelse(recompute.mode, "TRUE", "FALSE"), "))")
            rr = inla.eval(argument)
            r$cpo$cpo[idx] = rr$cpo$cpo[idx]
            r$cpo$pit[idx] = rr$cpo$pit[idx]
            r$cpo$failure[idx] = rr$cpo$failure[idx]
        }

        if (be.verbose)
            cat("\n")

        cpo.new = r$cpo$cpo[idx.fail]
        pit.new = r$cpo$pit[idx.fail]

        if (be.verbose) {
            print(cbind(index = idx.fail,
                        cpo.old = cpo.old, cpo.new = cpo.new,
                        pit.old = pit.old, pit.new = pit.new))
            
            cat("\nThe retured object contain the new values\n")
        }
    }

    return(r)
}
