#' Compute group.cv-values
#' 
#' From a fitted model, compute and add the `group.cv`-values
#' 
#' 
#' @return The object returned is list related to leave-group-out
#' cross-validation. See the vignette for details.
#' @author Havard Rue \email{hrue@@r-inla.org}
#' @seealso [control.compute()]
#' @rdname group.cv
#' @export inla.group.cv
`inla.group.cv` <-
    function(
             #' @param result An object of class `inla`, ie a result of a call to
             #' `inla()`.
             result,
             
             #' @param group.cv If given, the groups are taken from this argument.
             #' `group.cv` must be the output of previous call to
             #' `inla.group.cv()`.
             group.cv = NULL,
             
             #' @param num.level.sets Number of level.sets to use. The default value
             #' `-1` corresponds to leave-one-out cross-validation.
             num.level.sets = -1, 

             #' @param strategy One of `"posterior"` or `"prior"`. See the
             #' vignette for details.
             strategy = c("posterior", "prior"), 
             
             #' @param size.max The maximum size of a group. If the computed group-size is
             #' larger, it will be truncated to `size.max`.
             size.max = 32,
             
             #' @param groups An (optional) predefined list of groups.  See the vignette for
             #' details.
             groups = NULL, 

             #' @param selection An optional list of data-indices to use. If not given, then
             #' all data are used.
             selection = NULL, 

             #' @param group.selection An optional list of data-indices to use when building the
             #' groups. If given, each group beyond the observation itself, must be a subset of
             #' `group.selection`. If not given, then all data are used. 
             group.selection = NULL, 

             #' @param friends An optional list of lists of indices to use a friends
             friends = NULL, 

             #' @param verbose Run with `verbose` output of some of the internals in
             #' the calculations. This option will also enable `inla(...,
             #' verbose=TRUE)` if its not enabled already.
             verbose = FALSE, 

             #' @param epsilon Two correlations with a difference less than `epsilon`,
             #' will be classified as identical.
             epsilon = 0.005, 

             #' @param prior.diagonal When `strategy="prior"`, `prior.diagonal` is
             #' added to the diagonal of the prior precision matrix to avoid singularities
             prior.diagonal = 1e-4, 

             #' @param keep For `strategy="prior"`, then this gives a vector of the
             #' name of model-components TO USE when computing the groups. See the vignette
             #' for details. Not both of `keep` and `remove` can be defined.
             keep = NULL,

             #' @param remove For `strategy="prior"`, then this gives a vector of the
             #' name of model-components NOT TO USE when computing the groups. See the
             #' vignette for details. Not both of `keep` and `remove` can be
             #' defined.
             remove = NULL,
             
             #' @param remove.fixed For `strategy="prior"`, this is the default option
             #' which is in effect if both `keep` and `remove` are `NULL`. If
             #' `TRUE`, it will remove (or condition on) all fixed effects when
             #' computing the groups. See the vignette for details.
             remove.fixed = TRUE)
{
    stopifnot(!missing(result))
    stopifnot(inherits(result, "inla"))
    if (!is.null(group.cv)) {
        stopifnot(inherits(group.cv, "inla.group.cv"))
        get.groups <- function(cv) {
            lapply(seq_along(cv$groups),
                   function(ii, cv.arg) cv.arg$groups[[ii]]$idx,
                   cv.arg=cv)
        }
        return (inla.group.cv(result, groups = get.groups(group.cv), verbose = verbose))
    }

    cont.gcpo <- list(enable = TRUE,
                      num.level.sets = num.level.sets, 
                      size.max = size.max, 
                      strategy = match.arg(strategy), 
                      groups = groups, 
                      selection = selection, 
                      group.selection = group.selection, 
                      friends = friends, 
                      verbose = verbose, 
                      epsilon = epsilon, 
                      prior.diagonal = prior.diagonal, 
                      correct.hyperpar = FALSE, 
                      keep = keep, 
                      remove = remove, 
                      remove.fixed = remove.fixed)

    r <- result
    r$.args$control.compute$control.gcpo <- cont.gcpo
    r$.args$control.compute$return.marginals <- FALSE
    r$.args$control.compute$return.marginals.predictor <- FALSE
    r$.args$control.compute$dic <- FALSE
    r$.args$control.compute$cpo <- FALSE
    r$.args$control.compute$po <- FALSE
    r$.args$control.compute$waic <- FALSE
    r$.args$control.compute$config <- FALSE
    r$.args$control.compute$graph <- FALSE
    r$.args$control.compute$hyperpar <- FALSE
    r$.args$control.compute$q <- FALSE
    
    r$.args$control.fixed$correlation.matrix <- FALSE

    r$.args$control.inla$int.strategy <- "eb"
    r$.args$control.inla$use.directions <- r$misc$opt.directions

    r$.args$control.mode$result <- NULL
    r$.args$control.mode$restart <- FALSE
    r$.args$control.mode$theta <- r$mode$theta
    r$.args$control.mode$x <- r$mode$x
    r$.args$control.mode$fixed <- TRUE

    r$.args$inla.mode <- "compact"
    r$.args$verbose <- if (verbose) TRUE else r$.args$verbose
    r$.args$lincomb <- NULL
    r$.args$quantiles <- numeric(0)
    
    group.cv <- do.call("inla", args = r$.args)$gcpo
    group.cv$cv <- group.cv$gcpo
    group.cv$gcpo <- NULL
    r <- NULL

    class(group.cv) <- "inla.group.cv"
    return (group.cv)
}

