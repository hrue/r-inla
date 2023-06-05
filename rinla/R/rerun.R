#' Rerun an analysis
#' 
#' Rerun \code{\link{inla}} on an inla-object (output from \code{link{inla}})
#' 
#' 
#' @aliases inla.rerun rerun
#' @return This function will take the result in \code{object}, and rerun
#' \code{inla} again.  If \code{plain} is \code{FALSE}, start the optimization
#' from the mode in \code{object} so that we can obtain an improvement the mode
#' for the hyperparameters.  Otherwise, start from the same configuration as
#' for \code{object}.  The returned value is an \code{inla}-object.
#' @seealso \code{\link{inla}}
#' @examples
#' 
#' r = inla(y ~ 1,  data = data.frame(y=1:10))
#' r = inla.rerun(r)
#' 
#' @rdname rerun
#' @export inla.rerun
`inla.rerun` <- function(
                         #' @param object An \code{inla}-object, ie the output from an \code{inla}-call
                         object,
                         #' @param plain Logical. If \code{FALSE} (default), then make changes in
                         #' \code{object} to improve the performance
                         plain = FALSE)
{
    stopifnot(any(inherits(object, "inla")))

    ## need to do this, as if its true it will regenerate the linear combinations
    object$.args$control.fixed$correlation.matrix <- FALSE

    if (!plain) {
        object$.args$control.mode$result <- NULL
        object$.args$control.mode$restart <- TRUE
        object$.args$control.mode$theta <- object$mode$theta
        object$.args$control.mode$x <- object$mode$x
        ## to do a proper restart, we should keep the direction is any. if they are not used,
        ## they will be NULL, which is the same as FALSE.
        object$.args$control.inla$use.directions <- object$misc$opt.directions
    }

    new.obj <- do.call("inla", args = object$.args)
    ## revert, as variables are expanded...
    new.obj$call <- object$call

    return(new.obj)
}
