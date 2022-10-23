## Export: inla.group.cv

## !\name{inla.group.cv}
## !\alias{inla.group.cv}
## !
## !\title{Compute group.cv-values}
## !\description{
## ! From a fitted model, compute and add the \code{group.cv}-values
## !}
## !\usage{
## !inla.group.cv(result,
## !          num.level.sets = -1, 
## !          size.max = 32, 
## !          strategy = c("posterior", "prior"), 
## !          groups = NULL, 
## !          selection = NULL, 
## !          verbose = FALSE, 
## !          epsilon = 0.0025, 
## !          prior.diagonal = 1e-04, 
## !          keep = NULL, 
## !          remove = NULL, 
## !          remove.fixed = TRUE)
## !}
## !\arguments{
`inla.group.cv` <-
    function(
             ## !\item{result}{An object of class \code{inla}, ie a result
             ## !of a call to \code{inla()}.}
             result,
             
             ## !\item{num.level.sets}{Number of level.sets to use. The default value
             ## !\code{-1} corresponds to leave-one-out cross-validation.}
             num.level.sets = -1, 

             ## !\item{strategy}{One of \code{"posterior"} or \code{"prior"}.
             ## !See the  vignette for details.}
             strategy = c("posterior", "prior"), 
             
             ## !\item{size.max}{The maximum size of a group. If the computed
             ## !group-size is larger, it will be truncated to \code{size.max}.}
             size.max = 32,
             
             ## !\item{groups}{An (optional) predefined list of groups. 
             ## !See the  vignette for details.}
             groups = NULL, 

             ## !\item{selection}{An optional list of data-indices to use.
             ## !If not given, then all data are used.}
             selection = NULL, 

             ## !\item{verbose}{Run with \code{verbose} output of some of the internals
             ## !in the  calculations. This option will also
             ## !enable \code{inla(...,  verbose=TRUE)} if its not enabled already.}
             verbose = FALSE, 

             ## !\item{epsilon}{Two correlations with a difference less than
             ## !\code{epsilon}, will be classified as identical.}
             epsilon = 0.005, 

             ## !\item{prior.diagonal}{When \code{strategy="prior"}, \code{prior.diagonal}
             ## !is added to the diagonal of the prior precision matrix to avoid singularities}
             prior.diagonal = 1e-4, 

             ## !\item{keep}{For \code{strategy="prior"}, then this gives a vector of 
             ## !the name of model-components TO USE when computing the groups.
             ## !See the  vignette for details.
             ## !Not both of \code{keep} and \code{remove} can be defined.}
             keep = NULL,

             ## !\item{remove}{For \code{strategy="prior"}, then this gives a vector of 
             ## !the name of model-components NOT TO USE when computing the groups.
             ## !See the  vignette for details.
             ## !Not both of \code{keep} and \code{remove} can be defined.}
             remove = NULL,
             
             ## !\item{remove.fixed}{For \code{strategy="prior"}, this is the default
             ## !option which is in effect if both \code{keep} and \code{remove} are 
             ## !\code{NULL}. If \code{TRUE}, it will remove (or condition on) all fixed effects
             ## !when computing the groups.
             ## !See the  vignette for details.}
             remove.fixed = TRUE)
{
    ## !}
    ## !\value{The object returned is list related to leave-group-out cross-validation. See the vignette for details.}
    ## !\author{Havard Rue \email{hrue@r-inla.org}}
    ## !\seealso{\code{\link{control.compute}}}

    stopifnot(!missing(result))
    stopifnot(inherits(result, "inla"))

    cont.gcpo <- list(enable = TRUE,
                      num.level.sets = num.level.sets, 
                      size.max = size.max, 
                      strategy = match.arg(strategy), 
                      groups = groups, 
                      selection = selection, 
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

    r$.args$inla.mode <- "experimental"
    r$.args$verbose <- if (verbose) TRUE else r$.args$verbose
    r$.args$lincomb <- NULL
    r$.args$quantiles <- numeric(0)
    
    group.cv <- do.call("inla", args = r$.args)$gcpo
    group.cv$cv <- group.cv$gcpo
    group.cv$gcpo <- NULL
    r <- NULL

    return (group.cv)
}

