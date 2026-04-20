#' Improved estimates for the hyperparameters (classic-mode only)
#' 
#' Improve the estimates of the posterior marginals for the hyperparameters of
#' the model using the grid integration strategy. (classic-mode only)
#' 
#' 
#' @aliases inla.hyperpar hyperpar.inla
#' @return The object returned is the same as `object` but the estimates
#' of the hyperparameters are replaced by improved estimates.
#' @note This function might take a long time or if the number of
#' hyperparameters in the model is large. If it complains and says `I
#' cannot get enough memory`, try to increase the value of the argument
#' `dz` or decrease `diff.logdens`.
#' @author Havard Rue \email{hrue@@r-inla.org}
#' @seealso [inla()]
#' @references See the references in `inla`
#' @rdname hyperpar
#' @export inla.hyperpar
`inla.hyperpar` <-
    function(
             #' @param result An object of class `inla`, ie a result of a call to
             #' `inla()` in `classic` mode
             result,

             #' @param skip.configurations A boolean variable; skip configurations if the
             #' values at the main axis are to small. (Default TRUE)
             skip.configurations = TRUE,

             #' @param verbose Boolean indicating wheather the inla program should run in a
             #' verbose mode.
             verbose = FALSE,

             #' @param dz Step length in the standardized scale used in the construction of
             #' the grid, default 0.75.
             dz = 0.75,

             #' @param diff.logdens The difference of the log.density for the hyperpameters
             #' to stop numerical integration using int.strategy='grid'.  Default 15
             diff.logdens = 15,

             #' @param h The step-length for the gradient calculations for the
             #' hyperparameters. Default 0.01.
             h = NULL,

             #' @param restart A boolean defining whether the optimizer should start again
             #' to ind the mode or if it should use the mode contained in the `object`
             restart = FALSE,

             #' @param quantiles A vector of quantiles, to compute for each posterior
             #' marginal.
             quantiles = NULL,

             #' @param keep A boolean variable indicating the working files (ini file, data
             #' files and results files) should be kept
             keep = FALSE) 
{

    result.tmp <- result
    stopifnot(inherits(result.tmp, "inla"))

    if (result$.args$inla.mode != "classic") {
        warning("inla.hyperpar: Only for inla.mode='classic', skip...")
        return (result)
    }

    warning("inla.hyperpar: This function is deprecated, please do not use.")

    ## only arguments that we need to change,  is changed...
    result.tmp$.args$quantiles <- inla.ifelse(missing(quantiles) || is.null(quantiles), result.tmp$.args$quantiles, quantiles)
    result.tmp$.args$verbose <- verbose
    result.tmp$.args$only.hyperparam <- TRUE
    result.tmp$.args$keep <- keep

    result.tmp$.args$control.compute$hyperpar <- TRUE
    result.tmp$.args$control.compute$return.marginals <- TRUE
    result.tmp$.args$control.compute$dic <- FALSE
    result.tmp$.args$control.compute$mlik <- TRUE
    result.tmp$.args$control.compute$cpo <- FALSE
    result.tmp$.args$control.compute$po <- FALSE
    result.tmp$.args$control.compute$q <- FALSE
    result.tmp$.args$control.compute$graph <- FALSE

    result.tmp$.args$control.predictor$compute <- FALSE

    result.tmp$.args$control.inla$int.strategy <- "grid"
    result.tmp$.args$control.inla$strategy <- "gaussian"
    result.tmp$.args$control.inla$h <- inla.ifelse(is.null(h), result.tmp$.args$control.inla$h, h)
    result.tmp$.args$control.inla$dz <- dz
    result.tmp$.args$control.inla$diff.logdens <- diff.logdens
    result.tmp$.args$control.inla$print.joint.hyper <- TRUE
    result.tmp$.args$control.inla$force.diagonal <- TRUE
    result.tmp$.args$control.inla$adjust.weights <- FALSE
    result.tmp$.args$control.inla$use.directions <- FALSE
    result.tmp$.args$control.mode$result <- result.tmp
    result.tmp$.args$control.mode$restart <- restart

    ## cannot use this function with inla.call="submit", if so, then
    ## replace this into inla.call="remote"
    result.tmp$.args$inla.call <- sub("inla.submit", "inla.remote", result.tmp$.args$inla.call)

    result.tmp$.args$keep <- T
    ## call itself
    result.tmp <- inla.rerun(result.tmp, plain = TRUE)

    ## these are the entries that we want to replace
    replace.names <- c(
        "summary.hyperpar", "marginals.hyperpar", "internal.marginals.hyperpar",
        "internal.summary.hyperpar", "joint.hyper", "mlik", "version", "cpu.used"
    )

    for (nm in replace.names) {
        idx.result <- which(names(result) == nm)
        idx.result.tmp <- which(names(result.tmp) == nm)
        if (length(idx.result) > 0 && length(idx.result.tmp) > 0) {
            result[[idx.result]] <- result.tmp[[idx.result.tmp]]
        }
    }
    ## from $misc, we only want this one.
    result$misc$configs <- result.tmp$misc$configs

    return(result)
}
